#!/usr/bin/env python3
"""Build taxonomy and a matching WP dictionary from FINAL local reference FASTAs.
No downloads. Python standard library only. Existing files are not overwritten.
"""
import argparse
import csv
import gzip
import pickle
import re
from collections import defaultdict
from pathlib import Path

UNKNOWN = 'unclassified'
AMBIGUOUS = 'ambiguous'


class DataOnlyUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        raise pickle.UnpicklingError('Expected a plain-data parsing dictionary')


def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def fasta_metadata(path):
    header, length, seen = None, 0, set()
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if header is not None:
                    yield header.split()[0], header, length
                header, length = line[1:].strip(), 0
                if not header or header.split()[0] in seen:
                    raise ValueError(f'{path}: empty or duplicate FASTA identifier')
                seen.add(header.split()[0])
            elif line.strip():
                if header is None:
                    raise ValueError(f'{path}: sequence before FASTA header')
                length += len(''.join(line.split()))
        if header is not None:
            yield header.split()[0], header, length


def table(path):
    with open(path, newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_table(path, rows, fields):
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def genospecies(value):
    """Use the species epithet so genus spelling differences do not split bins."""
    value = (value or '').strip()
    match = re.search(r'\b(?:Borrelia|Borreliella)\s+([a-z][a-z-]+)\b', value)
    if match and match[1] not in {'sp', 'bacterium', 'uncultured'}:
        # Broad sensu-lato organism labels do not establish a genospecies.
        if re.search(r'\b(?:sensu lato|s\.?\s*l\.?)\b', value, re.I):
            return None
        return match[1]
    if re.fullmatch(r'[a-z][a-z-]+', value) and value not in {UNKNOWN, AMBIGUOUS, 'unknown', 'sp'}:
        return value
    return None


def strain_key(value):
    return re.sub(r'[^A-Z0-9]', '', (value or '').upper())


def resolve(candidates):
    known = set(candidates) - {UNKNOWN, '', None}
    if not known:
        return UNKNOWN, ''
    label = next(iter(known)) if len(known) == 1 and UNKNOWN not in candidates else AMBIGUOUS
    return label, ';'.join(sorted(known))


def build(args):
    out = Path(args.out)
    if out.exists():
        raise FileExistsError(f'Choose a new output directory: {out}')
    with open(args.lookup, 'rb') as handle:
        old_lookup = DataOnlyUnpickler(handle).load()
    strains = defaultdict(set)
    strain_names = defaultdict(set)

    def add_strain(name, label):
        if name and label:
            key = strain_key(name)
            strains[key].add(label)
            strain_names[key].add(name)

    if args.assemblies:
        for row in table(args.assemblies):
            label = genospecies(row.get('genospecies') or row.get('species'))
            add_strain(row.get('strain'), label)
            add_strain(row.get('ncbi_strain'), label)

    wp = list(fasta_metadata(args.wp))
    pf = list(fasta_metadata(args.pf32))
    if not wp or not pf:
        raise ValueError('Both reference FASTAs must contain records')
    for ident, header, _ in wp:
        if ident not in old_lookup:
            raise ValueError(f'{ident}: missing from --lookup; repair its strain/plasmid metadata first')
        add_strain(old_lookup[ident].get('strain'), genospecies(header))
    proteins = {}
    if args.genpept:
        with open_text(args.genpept) as handle:
            for record in re.split(r'^//\s*$', handle.read(), flags=re.M):
                accession = re.search(r'^VERSION\s+(\S+)', record, re.M)
                organism = re.search(r'^  ORGANISM\s+(.+)', record, re.M)
                strain = re.search(r'/(?:strain|isolate)="([^"]+)"', record)
                if accession and organism:
                    label = genospecies(organism[1])
                    proteins[accession[1]] = label
                    if strain:
                        add_strain(strain[1], label)

    def labels_for_protein(header):
        accessions = re.findall(r'(?<![A-Z0-9])(?:[A-Z]{3}\d+|[A-Z]{2}_\d+)\.\d+', header)
        direct = {proteins[acc] for acc in accessions if proteins.get(acc)}
        if direct:
            return direct
        # Older IDs often start with a locus tag, so match whole strain tokens
        # anywhere in the ID; do not mistake the first token for the strain.
        ident = header.split()[0]
        matched = []
        for key, names in strain_names.items():
            for name in names:
                pattern = r'(?:^|_)' + r'[-_]'.join(map(re.escape, re.split(r'[-_\s]+', name))) + r'(?:_|$)'
                if re.search(pattern, ident, re.I):
                    matched.append(key)
        if not matched:
            return {UNKNOWN}
        longest = max(map(len, matched))
        return set().union(*(strains[key] for key in matched if len(key) == longest))

    records = {}
    for ident, header, length in wp:
        info = old_lookup[ident]
        direct = genospecies(header)
        prior = genospecies(info.get('genospecies'))
        candidates = {direct} if direct else {prior} if prior else set(strains.get(strain_key(info.get('strain')), {UNKNOWN}))
        # Preserve ambiguity already recorded by an earlier metadata build.
        if info.get('genospecies') == AMBIGUOUS:
            raw = info.get('genospecies_candidates', [])
            candidates.update(raw.split(';') if isinstance(raw, str) else raw)
            candidates.add(UNKNOWN)
        records['wp', ident] = dict(candidates=candidates, strain=info.get('strain', ''), length=length)
    for ident, header, length in pf:
        records['pf32', ident] = dict(candidates=labels_for_protein(header), strain='', length=length)

    # Preserve taxon associations lost when identical sequences were collapsed.
    # Feed the formatter's headers.tsv here, including its duplicate alias rows.
    if args.audit:
        for row in table(args.audit):
            if row.get('status') == 'quarantined':
                continue
            database = {'wp.fna': 'wp', 'pf32.faa': 'pf32'}.get(row.get('file'))
            ident = row.get('retained_id') or row.get('new_id')
            target = records.get((database, ident))
            if target is None:
                continue
            labels = set(strains.get(strain_key(row.get('strain')), set()))
            if database == 'pf32':
                source_labels = labels_for_protein(row.get('old_id') or row.get('formatted_id') or ident)
                labels.update(source_labels - {UNKNOWN})
            if not labels and row.get('formatted_id') != ident:
                labels.add(UNKNOWN)
            target['candidates'].update(labels)
    else:
        print('NOTE: no --audit supplied; taxonomy describes retained references only, not removed aliases.')

    if args.overrides:
        for row in table(args.overrides):
            key = row['database'], row['reference_id']
            if key not in records:
                raise ValueError(f'Taxon override is not in final FASTAs: {key}')
            labels = {genospecies(value) for value in row['genospecies'].split(';')}
            records[key]['candidates'] = {value or UNKNOWN for value in labels}

    rows, lookup = [], {}
    for (database, ident), record in records.items():
        label, candidates = resolve(record['candidates'])
        rows.append(dict(database=database, reference_id=ident, genospecies=label,
                         genospecies_candidates=candidates, strain=record['strain']))
        if database == 'wp':
            lookup[ident] = dict(old_lookup[ident], ID=ident, length=record['length'],
                                 genospecies=label, genospecies_candidates=candidates.split(';') if candidates else [])
    out.mkdir(parents=True, exist_ok=False)
    fields = ['database', 'reference_id', 'genospecies', 'genospecies_candidates', 'strain']
    write_table(out / 'reference_taxa.tsv', rows, fields)
    unresolved = [row for row in rows if row['genospecies'] in {UNKNOWN, AMBIGUOUS}]
    write_table(out / 'taxa_unresolved.tsv', unresolved, fields)
    with open(out / 'blast_parsing_dict.pkl', 'wb') as handle:
        pickle.dump(lookup, handle, protocol=4)
    print(f'Wrote metadata for {len(wp)} WP and {len(pf)} PF32 references to {out}')
    print(f'{len(unresolved)} unknown/ambiguous references: {out / "taxa_unresolved.tsv"}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wp', required=True, help='Exact final WP FASTA used to build BLAST database')
    parser.add_argument('--pf32', required=True, help='Exact final protein FASTA used to build BLAST database')
    parser.add_argument('--lookup', required=True, help='Current whole-genome parsing dictionary')
    parser.add_argument('--assemblies', help='Assembly TSV containing species, strain and optional ncbi_strain')
    parser.add_argument('--audit', help='Formatter headers.tsv, including duplicate aliases')
    parser.add_argument('--genpept', help='Local pfam32_proteins.gp metadata')
    parser.add_argument('--overrides', help='TSV: database, reference_id, genospecies (semicolon separates ambiguous labels)')
    parser.add_argument('--out', default='taxon_metadata', help='New output directory')
    args = parser.parse_args()
    try:
        build(args)
    except (OSError, ValueError, KeyError, pickle.UnpicklingError) as exc:
        parser.exit(1, f'ERROR: {exc}\n')


if __name__ == '__main__':
    main()
