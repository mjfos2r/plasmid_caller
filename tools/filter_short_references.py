#!/usr/bin/env python3
"""Drop short and excluded WP references, keeping the caller's metadata in step with the FASTA.

Fragments below --min-bp are WGS scaffolds, not replicons. They produce spurious
low-coverage hits, so they must not reach the BLAST database. --exclude lists
contaminants, which return on every rebuild unless they are removed here.
Writes a NEW output directory; does not modify inputs. PF32 taxonomy rows pass through.
"""
import argparse
import csv
import pickle
from pathlib import Path


# Every genuine chromosome in the database falls well inside this range.
CHROMOSOME_BP = (500_000, 1_200_000)


class DataOnlyUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        raise pickle.UnpicklingError('Only plain-data parsing dictionaries are accepted')


def read_fasta(path):
    """Yield (id, header line, sequence lines) without altering the record text."""
    header, lines = None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if header is not None:
                    yield header[1:].split()[0], header, lines
                header, lines = line, []
            elif header is not None:
                lines.append(line)
    if header is not None:
        yield header[1:].split()[0], header, lines


def read_exclusions(path):
    with open(path, newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not {'reference_id', 'reason'}.issubset(reader.fieldnames or []):
            raise ValueError(f'{path}: require reference_id, reason columns')
        return {row['reference_id']: row['reason'] for row in reader}


def filter_fasta(source, output, min_bp, excluded=()):
    kept, dropped = {}, {}
    with open(output, 'w') as handle:
        for ident, header, lines in read_fasta(source):
            length = sum(len(line.strip()) for line in lines)
            if length < min_bp or ident in excluded:
                dropped[ident] = length
                continue
            if ident in kept:
                raise ValueError(f'{source}: repeated reference ID {ident}')
            kept[ident] = length
            handle.write(header)
            handle.writelines(lines)
    if not kept:
        raise ValueError(f'{source}: no references of at least {min_bp} bp')
    return kept, dropped


def filter_taxa(source, output, kept):
    with open(source, newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames or []
        if not {'database', 'reference_id', 'genospecies'}.issubset(fields):
            raise ValueError(f'{source}: require database, reference_id, genospecies columns')
        rows = list(reader)
    wp_ids = {row['reference_id'] for row in rows if row['database'] == 'wp'}
    if kept - wp_ids:
        raise ValueError(f'{source}: {len(kept - wp_ids)} kept references lack taxonomy, '
                         f'e.g. {sorted(kept - wp_ids)[:5]}')
    rows = [row for row in rows if row['database'] != 'wp' or row['reference_id'] in kept]
    with open(output, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def check_chromosomes(lookup, kept, lengths):
    low, high = CHROMOSOME_BP
    wrong = sorted(key for key in kept if lookup[key].get('name') == 'chromosome'
                   and not low <= lengths[key] <= high)
    if wrong:
        raise ValueError(f'{len(wrong)} references labelled chromosome are outside {low}-{high} bp, '
                         f'e.g. {wrong[:5]}')


def filter_lookup(source, output, kept, lengths):
    with open(source, 'rb') as handle:
        lookup = DataOnlyUnpickler(handle).load()
    if kept - set(lookup):
        raise ValueError(f'{source}: {len(kept - set(lookup))} kept references missing, '
                         f'e.g. {sorted(kept - set(lookup))[:5]}')
    check_chromosomes(lookup, kept, lengths)
    with open(output, 'wb') as handle:
        pickle.dump({key: value for key, value in lookup.items() if key in kept}, handle)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wp', required=True, type=Path, help='WP FASTA')
    parser.add_argument('--taxa', required=True, type=Path, help='reference_taxa.tsv')
    parser.add_argument('--lookup', required=True, type=Path, help='blast_parsing_dict.pkl')
    parser.add_argument('--min-bp', type=int, default=1000)
    parser.add_argument('--exclude', type=Path, help='TSV of reference_id, reason to remove')
    parser.add_argument('--out', required=True, type=Path, help='Staged database directory')
    args = parser.parse_args()

    (args.out / 'ref_seqs').mkdir(parents=True, exist_ok=True)
    excluded = read_exclusions(args.exclude) if args.exclude else {}
    lengths, dropped = filter_fasta(args.wp, args.out / 'ref_seqs' / 'wp.fna', args.min_bp, excluded)
    kept = set(lengths)
    taxa_rows = filter_taxa(args.taxa, args.out / 'reference_taxa.tsv', kept)
    filter_lookup(args.lookup, args.out / 'blast_parsing_dict.pkl', kept, lengths)
    print(f'Kept {len(kept)} WP references; dropped {len(dropped)} short or excluded; '
          f'{taxa_rows} taxonomy rows')
    for ident, length in sorted(dropped.items()):
        print(f'  dropped {ident}\t{length}')


if __name__ == '__main__':
    main()
