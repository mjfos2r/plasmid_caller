#!/usr/bin/env python3
"""Format local references for plasmid_caller 6.2.4. Standard library only.

No downloads. Inputs are the files produced by the earlier localization commands.
Outputs are merged references, deduplicated within each FASTA against the existing
DB and local inputs. Original files are never modified. DNA reverse complements
also count as duplicates; distinct alleles and circular rotations are not merged.
By default, all sides of conflicts go to review/ and are excluded from the clean
FASTAs and lookup. Review FASTAs contain one copy per sequence; entries.tsv keeps
all labels and source headers. Related supplied PF32 CDS are separated too.
--label-conflicts error stops instead; keep-first retains the first label
(existing DB first), but same-ID/different-sequence conflicts still stop.
"""
import argparse
import csv
import gzip
import hashlib
import json
import pickle
import re
from pathlib import Path

TYPE = re.compile(r'chromosome|(?:lp|cp)\d{1,2}(?:[-+#](?:lp|cp)?\d{1,2})*')
PROTEIN = re.compile(r'(?<![A-Z0-9])(?:[A-Z]{3}\d+|[A-Z]{2}_\d+)\.\d+')


def text_open(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def fasta(path):
    with text_open(path) as f:
        header, seq = None, []
        for line in f:
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq)
                header, seq = line[1:].strip(), []
            elif line.strip():
                if header is None:
                    raise ValueError(f'{path}: sequence before first header')
                seq.append(line.strip())
        if header is not None:
            yield header, ''.join(seq)


def write_fasta(handle, header, sequence):
    handle.write('>' + header + '\n')
    for i in range(0, len(sequence), 80):
        handle.write(sequence[i:i+80] + '\n')


def table(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_table(path, rows, fields):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def genbank_metadata(path):
    """Read source qualifiers only; never infer strain from a protein title."""
    if not path.exists():
        return {}
    result = {}
    with text_open(path) as f:
        content = f.read()
    for record in re.split(r'^//\s*$', content, flags=re.M):
        accession = re.search(r'^VERSION\s+(\S+)', record, re.M)
        if not accession:
            continue
        source = re.search(r'^     source\s+.*?(?=^     \S|^ORIGIN|^CONTIG|\Z)',
                           record, re.M | re.S)
        qualifiers = {}
        if source:
            for key, value in re.findall(r'/(\w+)="([^"]*)"', source[0]):
                qualifiers[key] = ' '.join(value.split())
        result[accession[1]] = qualifiers
    return result


def canonical_type(value):
    value = (value or '').replace('_', '-').replace('\u2013', '-').strip()
    return value if TYPE.fullmatch(value) else None


def named_type(value):
    matches = TYPE.findall((value or '').replace('_', '-').replace('\u2013', '-'))
    return matches[0] if len(set(matches)) == 1 else None


def scaffold_type(header):
    """Name the replicon in a WGS scaffold title, such as A91S_cp32-3lp25_to_PBi.

    NCBI files unplaced scaffolds under location "Chromosome", so that field says
    nothing about them. A title with no replicon token returns None and is skipped.
    """
    title = header.split(',')[0]
    if re.search(r'_chr(?:_|$)', title):
        return 'chromosome'
    return named_type(re.sub(r'(?<=\d)(?=lp|cp)', '-', title))


def pf32_header(row, metadata):
    accession = row['protein_accession']
    info = metadata.get(accession, {})
    original_strain = row.get('strain') or info.get('strain') or info.get('isolate')
    if not original_strain:
        raise ValueError(f'{accession}: source strain absent from GenPept and manifest')
    # The caller takes the first underscore-delimited token as the strain.
    strain = re.sub(r'[^A-Za-z0-9.-]+', '-', original_strain).strip('-')
    kind = canonical_type(row['pfam32_type'])
    if not strain or not kind:
        raise ValueError(f'{accession}: invalid strain/plasmid label')
    # The final ParA suffix is opaque to the caller. Accession makes it unique.
    # Avoid underscores in RefSeq accessions being mistaken for strain fields.
    suffix = accession.replace('_', '-')
    ident = f'{strain}_{kind}_ParA_{suffix}'
    nucleotide = row.get('nucleotide_accession') or 'NA'
    return ident, f'{ident} {nucleotide}-{accession}', original_strain, strain, kind


class DataOnlyUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        raise pickle.UnpicklingError('Only plain-data lookup dictionaries are accepted')


class References:
    """Keep the first sequence representative, but retain every alias in the audit."""
    def __init__(self, filename, nucleotide=False):
        self.filename = filename
        self.nucleotide = nucleotide
        self.records, self.by_id, self.by_sequence = [], {}, {}
        self.audit, self.conflicts, self.entries = [], [], []

    def add(self, header, sequence, label, strain, source, old_id=None,
            caller_strain=None, source_header=None):
        ident = header.split()[0]
        sequence = ''.join(sequence.split()).upper()
        # Ignore one terminal translation stop; internal stops remain significant.
        if not self.nucleotide and sequence.endswith('*'):
            sequence = sequence[:-1]
        allowed = set('ACGTRYSWKMBDHVN' if self.nucleotide else 'ABCDEFGHIJKLMNOPQRSTUVWXYZ*')
        if not sequence or set(sequence) - allowed:
            raise ValueError(f'{source}: empty or invalid sequence for {ident}')
        if not label:
            raise ValueError(f'{source}: missing plasmid label for {ident}')
        normalized = sequence
        if self.nucleotide:
            reverse = sequence.translate(str.maketrans('ACGTRYSWKMBDHVN',
                                                      'TGCAYRSWMKVHDBN'))[::-1]
            normalized = min(sequence, reverse)
        digest = hashlib.sha256(normalized.encode('ascii')).hexdigest()
        previous_id = self.by_id.get(ident)
        existing = self.by_sequence.get(digest)
        if previous_id and previous_id['digest'] != digest:
            self.conflicts.append(dict(file=self.filename, id=ident,
                retained_id=previous_id['id'], source=source,
                retained_source=previous_id['source'], plasmid_type=label,
                retained_type=previous_id['label'], reason='Same ID, different sequence'))
        if existing and existing['label'] != label:
            self.conflicts.append(dict(file=self.filename, id=ident,
                retained_id=existing['id'], source=source,
                retained_source=existing['source'], plasmid_type=label,
                retained_type=existing['label'], reason='Identical sequence, conflicting plasmid labels'))
        record = dict(id=ident, header=header, sequence=sequence, digest=digest,
                      label=label, strain=strain, source=source, old_id=old_id or ident,
                      source_header=source_header or header)
        self.entries.append(record)
        representative = existing or record
        self.by_id.setdefault(ident, record)
        if not existing:
            self.records.append(record)
            self.by_sequence[digest] = record
        self.audit.append(dict(file=self.filename, source=source,
            old_id=old_id or ident, formatted_id=ident, retained_id=representative['id'],
            status='duplicate' if existing else 'retained', strain=strain,
            caller_strain=caller_strain or strain, plasmid_type=label,
            sequence_sha256=digest))

    def write(self, out):
        with open(out / self.filename, 'w') as f:
            for record in self.records:
                write_fasta(f, record['header'], record['sequence'])


def conflict_members(refs, extra_ids=()):
    """Exclude whole connected groups: shared sequence OR shared identifier."""
    ids = set(extra_ids)
    for conflict in refs.conflicts:
        ids.update((conflict['id'], conflict['retained_id']))
    digests = set()
    while True:
        before = len(ids), len(digests)
        for record in refs.entries:
            if record['id'] in ids or record['digest'] in digests:
                ids.add(record['id'])
                digests.add(record['digest'])
        if before == (len(ids), len(digests)):
            return ids, digests


def separate_conflicts(sets, out):
    members = {refs.filename: conflict_members(refs) for refs in sets}
    # Link supplied protein/CDS records by their formatted IDs. Exclude linked
    # aliases recursively if either side has a conflict, including synonymous CDS.
    pf_sets = [refs for refs in sets if refs.filename in ('pf32.faa', 'pf32_cds.fna')]
    while True:
        before = {name: (frozenset(ids), frozenset(digests))
                  for name, (ids, digests) in members.items()}
        linked_ids = set().union(*(members[refs.filename][0] for refs in pf_sets))
        for refs in pf_sets:
            members[refs.filename] = conflict_members(refs, linked_ids)
        after = {name: (frozenset(ids), frozenset(digests))
                 for name, (ids, digests) in members.items()}
        if before == after:
            break
    review = out / 'review'
    review.mkdir()
    aliases = []
    for refs in sets:
        _, digests = members[refs.filename]
        representatives = {}
        seen_aliases = set()
        for record in refs.entries:
            if record['digest'] not in digests:
                continue
            digest = record['digest']
            representatives.setdefault(digest, record)
            alias = dict(file=refs.filename, review_id='review_' + digest,
                         original_id=record['old_id'], formatted_id=record['id'],
                         original_header=record['source_header'], source=record['source'],
                         plasmid_type=record['label'], strain=record['strain'],
                         orientation=('same' if record['sequence'] == representatives[digest]['sequence']
                                      else 'reverse-complement'))
            key = tuple(alias.values())
            if key not in seen_aliases:
                aliases.append(alias)
                seen_aliases.add(key)
        if representatives:
            with open(review / refs.filename, 'w') as f:
                for digest, record in representatives.items():
                    write_fasta(f, 'review_' + digest, record['sequence'])
        refs.records = [r for r in refs.records if r['digest'] not in digests]
        for row in refs.audit:
            if row['sequence_sha256'] in digests:
                row['status'] = 'quarantined'
                row['retained_id'] = ''
        print(f'{refs.filename}: {len(representatives)} unique sequences separated for review')
    write_table(review / 'entries.tsv', aliases,
                ['file', 'review_id', 'original_id', 'formatted_id', 'original_header',
                 'source', 'plasmid_type', 'strain', 'orientation'])
    (review / 'README.txt').write_text(
        'These records are excluded from the clean FASTAs and lookup.\n'
        'Each review FASTA stores one copy per normalized sequence (DNA reverse complements merged).\n'
        'entries.tsv maps every source/label to its review sequence and orientation.\n'
        'conflicts.tsv in the parent directory lists the directly detected conflicts.\n'
        'Protein/CDS aliases linked to a conflict are also excluded.\n'
        'Resolve labels or identifiers in the source inputs, then rerun into a new output directory.\n')


def existing_pf32_label(ident):
    # Mirror the existing caller's special-case label and type parsing.
    if 'HO14_NovelPFam32_DatabaseShortContig' in ident:
        return 'HO14_NovelPFam32'
    matches = TYPE.findall(ident.replace('_', '-'))
    return matches[0] if matches else None


def read_pf32_corrections(path):
    if not path:
        return {}
    corrections = {row['old_id']: row['new_id'] for row in table(path)}
    if len(set(corrections.values())) != len(corrections):
        raise ValueError(f'{path}: two entries share a new_id')
    return corrections


def correct_pf32_header(header, corrections, applied):
    """Rename a seed PF32 record. The caller reads the plasmid name from this ID."""
    ident, _, rest = header.partition(' ')
    if ident not in corrections:
        return header
    applied.add(ident)
    return f'{corrections[ident]} {rest}'.rstrip()


def load_existing(args, pf32, wp):
    db = Path(args.existing_db)
    with open(db / 'blast_parsing_dict.pkl', 'rb') as f:
        lookup = DataOnlyUnpickler(f).load()
    pf_path = db / args.existing_pf32
    wp_path = db / args.existing_wp
    corrections, applied = read_pf32_corrections(args.pf32_corrections), set()
    for header, sequence in fasta(pf_path):
        header = correct_pf32_header(header, corrections, applied)
        ident = header.split()[0]
        strain = '_'.join(ident.split('_')[:2]) if ident.startswith('NE_') else ident.split('_')[0]
        pf32.add(header, sequence, existing_pf32_label(ident), strain, str(pf_path))
    if set(corrections) - applied:
        raise ValueError(f'PF32 corrections matched nothing: {sorted(set(corrections) - applied)}')
    for header, sequence in fasta(wp_path):
        ident = header.split()[0]
        if ident not in lookup:
            raise ValueError(f'{wp_path}: {ident} absent from existing lookup dictionary')
        info = lookup[ident]
        wp.add(header, sequence, info['name'], info['strain'], str(wp_path))


def format_references(args):
    localized, out = Path(args.localized), Path(args.out)
    if out.exists():
        raise FileExistsError(f'Choose a new output directory: {out}')
    pf32 = References('pf32.faa')
    cds = References('pf32_cds.fna', nucleotide=True)
    wp = References('wp.fna', nucleotide=True)
    load_existing(args, pf32, wp)
    rows = table(args.pfam32)
    required = {'pfam32_type', 'source', 'protein_accession'}
    if not rows or not required.issubset(rows[0]):
        raise ValueError('Use the expanded pfam32_references.tsv manifest')
    metadata = genbank_metadata(localized / 'pfam32_proteins.gp')
    converted = {}
    for row in rows:
        if row['source'] != 'GenBank' or not row['protein_accession']:
            continue
        accession = row['protein_accession']
        value = pf32_header(row, metadata)
        if accession in converted and converted[accession] != value:
            raise ValueError(f'Conflicting manifest rows for {accession}')
        converted[accession] = value
    skipped = []
    for row in rows:
        if row['source'] != 'GenBank' or not row['protein_accession']:
            skipped.append({'source': 'PFam32 manifest', 'id': row['pfam32_type'],
                            'reason': 'No supported local GenBank protein accession'})

    for input_name, references in [('pfam32_proteins.faa', pf32),
                                    ('pfam32_cds.fna', cds)]:
        path = localized / input_name
        if not path.exists():
            if input_name.endswith('.faa'):
                raise FileNotFoundError(path)
            continue
        found = set()
        for old_header, sequence in fasta(path):
            candidates = set(PROTEIN.findall(old_header)) & set(converted)
            if len(candidates) != 1:
                raise ValueError(f'{path}: cannot map header unambiguously: {old_header}')
            accession = candidates.pop()
            found.add(accession)
            ident, new_header, original_strain, strain, kind = converted[accession]
            references.add(new_header, sequence, kind, original_strain, str(path),
                           old_id=old_header.split()[0], caller_strain=strain, source_header=old_header)
        missing = set(converted) - found
        if missing:
            raise ValueError(f'{path}: missing manifest accessions: {sorted(missing)}')

    # Optional S3 accession-specific names override the names deposited in NCBI.
    overrides = {}
    if args.renames:
        for row in table(args.renames):
            name = row['paper_plasmid_name'].strip()
            if not name:
                raise ValueError('Empty paper plasmid name in corrections table')
            for accession in re.findall(r'\b[A-Z]{1,4}\d+(?:\.\d+)?\b', row['nucleotide_accession']):
                if accession in overrides and overrides[accession] != name:
                    raise ValueError(f'Conflicting name corrections for {accession}')
                overrides[accession] = name
    # Explicit accession-specific corrections apply to both old and new WP records.
    # Rebuild the collection so corrected labels participate in conflict detection.
    if overrides:
        corrected = References('wp.fna', nucleotide=True)
        db = Path(args.existing_db)
        with open(db / 'blast_parsing_dict.pkl', 'rb') as f:
            old_lookup = DataOnlyUnpickler(f).load()
        for header, sequence in fasta(db / args.existing_wp):
            ident = header.split()[0]
            info = old_lookup[ident]
            label = overrides.get(ident) or overrides.get(ident.split('.')[0]) or info['name']
            corrected.add(header, sequence, label, info['strain'], str(db / args.existing_wp))
        wp = corrected

    data = localized / 'genomes' / 'ncbi_dataset' / 'data'
    for assembly in table(args.assemblies):
        accession = assembly['assembly_accession']
        if not accession:
            skipped.append(dict(source='assemblies.tsv', id=assembly['strain'], reason='No assembly accession'))
            continue
        directory = data / accession
        genomes = [p for p in directory.glob('*_genomic.fna') if p.name != 'cds_from_genomic.fna']
        if len(genomes) != 1:
            skipped.append(dict(source=accession, id=assembly['strain'], reason='Expected one local genomic FASTA'))
            continue
        reports = {}
        report_file = directory / 'sequence_report.jsonl'
        if report_file.exists():
            with open(report_file) as f:
                for line in f:
                    if line.strip():
                        rec = json.loads(line)
                        for key in ('genbankAccession', 'refseqAccession', 'genbank_accession', 'refseq_accession'):
                            if rec.get(key):
                                reports[rec[key]] = rec
        source_info = genbank_metadata(directory / 'genomic.gbff')
        for header, sequence in fasta(genomes[0]):
            record_id = header.split()[0]
            info = source_info.get(record_id, {})
            report = reports.get(record_id, {})
            location = report.get('assignedMoleculeLocationType', report.get('assigned_molecule_location_type', ''))
            name = overrides.get(record_id) or overrides.get(record_id.split('.')[0])
            # Unplaced scaffolds carry location and chrName "chromosome" whatever they are.
            unplaced = location.lower() == 'chromosome' and report.get('role') != 'assembled-molecule'
            if not name:
                name = scaffold_type(header) if unplaced else 'chromosome' if location.lower() == 'chromosome' else None
            if not name:
                chr_name = '' if unplaced else report.get('chrName', report.get('chr_name', ''))
                name = named_type(info.get('plasmid')) or named_type(chr_name)
            if not name:
                match = re.search(r'\b(plasmid\s+.*|chromosome\b)', header, re.I)
                name = named_type(match[0].lower()) if match else None
            if not name:
                skipped.append(dict(source=accession, id=record_id, reason='No unambiguous chromosome/plasmid type'))
                continue
            wp.add(header, sequence, name, assembly['strain'], str(genomes[0]))
    if not wp.records:
        raise ValueError('No whole-genome references resolved')

    sets = [pf32, wp] + ([cds] if cds.records else [])
    out.mkdir(parents=True, exist_ok=False)
    conflicts = [row for refs in sets for row in refs.conflicts]
    if args.label_conflicts == 'separate':
        separate_conflicts(sets, out)
    audit = [row for refs in sets for row in refs.audit]
    write_table(out / 'headers.tsv', audit,
                ['file', 'source', 'old_id', 'formatted_id', 'retained_id', 'status',
                 'strain', 'caller_strain', 'plasmid_type', 'sequence_sha256'])
    write_table(out / 'duplicates.tsv', [row for row in audit if row['status'] == 'duplicate'],
                ['file', 'source', 'old_id', 'formatted_id', 'retained_id', 'status',
                 'strain', 'caller_strain', 'plasmid_type', 'sequence_sha256'])
    write_table(out / 'unresolved.tsv', skipped, ['source', 'id', 'reason'])
    write_table(out / 'conflicts.tsv', conflicts,
                ['file', 'id', 'retained_id', 'source', 'retained_source',
                 'plasmid_type', 'retained_type', 'reason'])
    fatal = [c for c in conflicts if args.label_conflicts == 'error' or
             (args.label_conflicts == 'keep-first' and c['reason'] == 'Same ID, different sequence')]
    if fatal:
        raise ValueError(f'{len(fatal)} conflicts: review {out / "conflicts.tsv"}. '
                         'No FASTAs or lookup written. Resolve labels/IDs in inputs, '
                         'then rerun with a new --out directory.')
    if conflicts and args.label_conflicts == 'keep-first':
        print(f'WARNING: kept first labels for {len(conflicts)} conflicting entries; '
              f'review {out / "conflicts.tsv"}')
    for refs in sets:
        refs.write(out)
    lookup = {r['id']: dict(ID=r['id'], name=r['label'], length=len(r['sequence']),
                           strain=r['strain']) for r in wp.records}
    with open(out / 'blast_parsing_dict.pkl', 'wb') as f:
        pickle.dump(lookup, f, protocol=4)
    for refs in sets:
        removed = sum(row['status'] == 'duplicate' for row in refs.audit)
        print(f'{refs.filename}: {len(refs.records)} retained; {removed} duplicates removed')
    print(f'Merged FASTAs and matching lookup: {out}')
    print(f'{len(skipped)} unresolved entries: {out / "unresolved.tsv"}')
    print('pf32_cds.fna, when present, contains only locally supplied CDS references; '
          'its IDs need not match the retained protein representatives. See headers.tsv.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--existing-db', required=True,
                   help='Existing plasmid_caller db directory containing ref_seqs and blast_parsing_dict.pkl')
    p.add_argument('--existing-pf32', default='ref_seqs/pf32_v4.faa',
                   help='Existing protein FASTA path relative to --existing-db (or absolute)')
    p.add_argument('--existing-wp', default='ref_seqs/wp_v2.fna',
                   help='Existing whole-genome FASTA path relative to --existing-db (or absolute)')
    p.add_argument('--localized', default='localized', help='Directory from the localization commands')
    p.add_argument('--assemblies', default='assemblies.tsv')
    p.add_argument('--pfam32', default='pfam32_references.tsv', help='Expanded PFam32 TSV')
    p.add_argument('--pf32-corrections', help='TSV of old_id, new_id for seed PF32 records')
    p.add_argument('--renames', help='Optional accession-specific plasmid_name_corrections.tsv from S3')
    p.add_argument('--label-conflicts', choices=['separate', 'error', 'keep-first'], default='separate',
                   help='Default: separate all conflicting entries into review/ and write clean merged references. error stops; keep-first explicitly retains the first label.')
    p.add_argument('--out', default='formatted_refs', help='New output directory; must not already exist')
    args = p.parse_args()
    try:
        format_references(args)
    except (ValueError, OSError, KeyError, pickle.UnpicklingError) as exc:
        p.exit(1, f'ERROR: {exc}\n')


if __name__ == '__main__':
    main()
