#!/usr/bin/env python3
"""Sync manually edited taxonomy into the caller's metadata. No inference/downloads.

Accepts the complete reference_taxa.tsv or an edited taxa_unresolved.tsv subset.
The installed reference_taxa.tsv supplies rows not present in the edited file.
Writes a NEW output directory; does not modify inputs or BLAST databases.
"""
import argparse
import csv
import pickle
import re
from pathlib import Path


class DataOnlyUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        raise pickle.UnpicklingError('Only plain-data parsing dictionaries are accepted')


def read_table(path):
    with open(path, newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames or []
        if not {'database', 'reference_id', 'genospecies'}.issubset(fields):
            raise ValueError(f'{path}: require database, reference_id, genospecies columns')
        result = {}
        for row in reader:
            if None in row:
                raise ValueError(f'{path}: row has more values than column headers')
            row = {key: (value or '').strip() for key, value in row.items()}
            key = row['database'], row['reference_id']
            if key[0] not in {'wp', 'pf32'} or not key[1]:
                raise ValueError(f'{path}: invalid reference key {key}')
            if key in result:
                raise ValueError(f'{path}: duplicate reference key {key}')
            result[key] = row
    return fields, result


def label(value):
    value = ' '.join(value.strip().lower().split())
    if value == 'unknown':
        return 'unclassified'
    match = re.fullmatch(r'(?:borrelia|borreliella) ([a-z][a-z-]+)(?: sensu stricto)?', value)
    if match:
        value = match[1]
    if not re.fullmatch(r'[a-z][a-z-]+', value) or value in {'sp', 'bacterium', 'uncultured'}:
        raise ValueError(f'Invalid genospecies label: {value!r}; use an epithet such as afzelii')
    return value


def normalize(row):
    value = row['genospecies']
    if ';' in value:
        candidates = {label(part) for part in value.split(';') if part.strip()}
        if not candidates or candidates & {'ambiguous', 'unclassified'}:
            raise ValueError(f'{row["reference_id"]}: explicit candidate list must contain genospecies names')
        name = next(iter(candidates)) if len(candidates) == 1 else 'ambiguous'
    else:
        name = label(value)
        if name == 'ambiguous':
            candidates = {label(part) for part in row.get('genospecies_candidates', '').split(';') if part.strip()}
            candidates -= {'ambiguous', 'unclassified'}
        elif name == 'unclassified':
            candidates = set()
        else:
            # A manually assigned single taxon overrides ALL stale candidates.
            candidates = {name}
    row['genospecies'] = name
    row['genospecies_candidates'] = ';'.join(sorted(candidates))
    return row


def sync(edited, db, out, reference=None):
    db, out = Path(db), Path(out)
    if out.exists():
        raise FileExistsError(f'Choose a new output directory: {out}')
    baseline = Path(reference) if reference else db / 'reference_taxa.tsv'
    fields, rows = read_table(baseline)
    edit_fields, edits = read_table(edited)
    unexpected = set(edits) - set(rows)
    if unexpected:
        raise ValueError(f'Edited IDs absent from the full reference table: {sorted(unexpected)[:10]}')
    for key, edit in edits.items():
        rows[key].update(edit)
    for row in rows.values():
        normalize(row)
    fields = list(dict.fromkeys(fields + edit_fields + ['genospecies_candidates']))
    with open(db / 'blast_parsing_dict.pkl', 'rb') as handle:
        lookup = DataOnlyUnpickler(handle).load()
    wp_ids = {ident for database, ident in rows if database == 'wp'}
    if wp_ids != set(lookup):
        missing = sorted(wp_ids - set(lookup))
        stale = sorted(set(lookup) - wp_ids)
        raise ValueError(f'WP table/dictionary mismatch: missing keys {missing[:5]}, extra keys {stale[:5]}. '
                         'Use a full reference table and dictionary from the same build.')
    for ident in wp_ids:
        row = rows['wp', ident]
        lookup[ident]['genospecies'] = row['genospecies']
        raw = row['genospecies_candidates']
        lookup[ident]['genospecies_candidates'] = raw.split(';') if raw else []
    unresolved = [row for row in rows.values() if row['genospecies'] in {'ambiguous', 'unclassified'}]
    out.mkdir(parents=True, exist_ok=False)
    for filename, data in [('reference_taxa.tsv', rows.values()), ('taxa_unresolved.tsv', unresolved)]:
        with open(out / filename, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(data)
    with open(out / 'blast_parsing_dict.pkl', 'wb') as handle:
        pickle.dump(lookup, handle, protocol=4)
    print(f'Applied {len(edits)} edited rows; wrote {len(rows)} reference rows and {len(lookup)} WP dictionary entries to {out}')
    print(f'{len(unresolved)} references remain ambiguous/unclassified. Input files were preserved.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--edited', required=True, help='Your edited full or partial taxon TSV')
    parser.add_argument('--db', default='src/plasmid_caller/db', help='Directory containing the matching parsing dictionary')
    parser.add_argument('--reference', help='Full baseline reference_taxa.tsv; defaults to --db/reference_taxa.tsv')
    parser.add_argument('--out', default='taxon_metadata_manual', help='New output directory')
    args = parser.parse_args()
    try:
        sync(args.edited, args.db, args.out, args.reference)
    except (OSError, ValueError, KeyError, pickle.UnpicklingError) as exc:
        parser.exit(1, f'ERROR: {exc}\n')


if __name__ == '__main__':
    main()
