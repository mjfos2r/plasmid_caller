"""Reference-based genospecies annotation and per-assembly composition."""
import csv
import math
import warnings
from pathlib import Path

import pandas as pd

UNCLASSIFIED = 'unclassified'
AMBIGUOUS = 'ambiguous'


def reference_tokens(*values):
    """Accept bare, lcl|, gb| and ref| IDs without dropping accession versions."""
    for value in values:
        if not value:
            continue
        token = str(value).split()[0]
        yield token
        if token.startswith('gnl|BL_ORD_ID|'):
            continue
        for part in token.split('|'):
            if part and part not in {'lcl', 'gb', 'ref', 'emb', 'dbj', 'sp', 'pdb', 'gi'}:
                yield part


def lookup_reference(hit_id, lookup, hit_def='', accession=''):
    for token in reference_tokens(hit_id, hit_def, accession):
        if token in lookup:
            return lookup[token]
    raise KeyError(f'Reference missing from blast_parsing_dict.pkl: {hit_id!r}. '
                   'Build the dictionary from the exact FASTA used for this BLAST DB.')


def load_reference_taxa(db_dir):
    path = Path(db_dir) / 'reference_taxa.tsv'
    if not path.exists():
        warnings.warn(f'{path} absent; missing reference taxonomy will be unclassified.')
        return {}
    with path.open(newline='') as handle:
        rows = csv.DictReader(handle, delimiter='\t')
        if not {'database', 'reference_id', 'genospecies'}.issubset(rows.fieldnames or []):
            raise ValueError(f'{path}: invalid taxonomy columns')
        result = {}
        for row in rows:
            key = row['database'], row['reference_id']
            if key in result and result[key] != row:
                raise ValueError(f'{path}: conflicting taxonomy rows for {key}')
            result[key] = row
        return result


def hit_taxonomy(database, alignment, taxa, lookup):
    for token in reference_tokens(alignment.hit_id, getattr(alignment, 'hit_def', ''),
                                  getattr(alignment, 'accession', '')):
        info = taxa.get((database, token))
        if info is not None:
            return info.get('genospecies') or UNCLASSIFIED, info.get('genospecies_candidates', '')
    if database == 'wp':
        info = lookup_reference(alignment.hit_id, lookup,
                                getattr(alignment, 'hit_def', ''), getattr(alignment, 'accession', ''))
        candidates = info.get('genospecies_candidates', '')
        if isinstance(candidates, (list, tuple)):
            candidates = ';'.join(candidates)
        return info.get('genospecies') or UNCLASSIFIED, candidates
    return UNCLASSIFIED, ''


def _candidate_set(row):
    label = row.get('genospecies')
    raw = row.get('genospecies_candidates')
    labels = set(str(raw).split(';')) if pd.notna(raw) else set()
    if pd.notna(label) and label not in {UNCLASSIFIED, AMBIGUOUS}:
        labels.add(str(label))
    return labels - {'', UNCLASSIFIED, AMBIGUOUS}


def annotate_best_hits(full, best, database, min_cov_pct=50, min_cov_bp=1000):
    """Preserve no-hit queries, and flag co-best WP taxon disagreements."""
    result = best.copy()
    if database == 'wp':
        cov = pd.to_numeric(full['query_coverage_percent'], errors='coerce')
        bp = pd.to_numeric(full['query_covered_length'], errors='coerce')
        pid = pd.to_numeric(full['overall_percent_identity'], errors='coerce')
        eligible = full.loc[(cov >= min_cov_pct) | (bp >= min_cov_bp)].copy()
        eligible['_taxon_score'] = (cov * pid).loc[eligible.index]
        calls = {}
        for contig_id, group in eligible.groupby('contig_id', sort=False):
            peak = group['_taxon_score'].max()
            tied = group.loc[group['_taxon_score'].apply(
                lambda score: pd.notna(score) and math.isclose(float(score), float(peak), rel_tol=1e-9, abs_tol=1e-9))]
            candidates = set()
            unresolved = False
            for _, row in tied.iterrows():
                candidates.update(_candidate_set(row))
                value = row.get('genospecies')
                unresolved |= pd.isna(value) or value in {UNCLASSIFIED, AMBIGUOUS}
            label = (next(iter(candidates)) if len(candidates) == 1 and not unresolved
                     else AMBIGUOUS if candidates else UNCLASSIFIED)
            calls[contig_id] = label, ';'.join(sorted(candidates))
        result['genospecies'] = result['contig_id'].map(lambda key: calls.get(key, (UNCLASSIFIED, ''))[0])
        result['genospecies_candidates'] = result['contig_id'].map(lambda key: calls.get(key, (UNCLASSIFIED, ''))[1])
    query_columns = ['assembly_id', 'contig_id', 'contig_len', 'query_length']
    queries = full[query_columns].drop_duplicates('contig_id')
    missing = queries.loc[~queries['contig_id'].isin(result['contig_id'])].copy()
    missing['genospecies'], missing['genospecies_candidates'] = UNCLASSIFIED, ''
    if not missing.empty:
        result = pd.concat([result, missing], ignore_index=True)
    result['genospecies'] = result['genospecies'].fillna(UNCLASSIFIED)
    result['genospecies_candidates'] = result['genospecies_candidates'].fillna('')
    return result


def add_genospecies_calls(summary, min_contig_bp=1000):
    """Assign each contig from WP evidence; PF32 taxonomy remains hit metadata."""
    result = summary.copy()
    if 'genospecies_wp' in result:
        result['genospecies'] = result['genospecies_wp'].fillna(UNCLASSIFIED)
        result['genospecies_candidates'] = result['genospecies_candidates_wp'].fillna('')
    else:
        result['genospecies'], result['genospecies_candidates'] = UNCLASSIFIED, ''
    short = pd.to_numeric(result['contig_len'], errors='coerce') < min_contig_bp
    result.loc[short, ['genospecies', 'genospecies_candidates']] = [UNCLASSIFIED, '']
    return result


def write_genospecies_composition(summary, output):
    """One contig once; full contig length, not summed alignments or read abundance."""
    frame = summary[['assembly_id', 'contig_id', 'contig_len', 'genospecies']].copy()
    if frame.duplicated(['assembly_id', 'contig_id']).any():
        raise ValueError('Duplicate assembly/contig rows would inflate taxon composition')
    frame['contig_len'] = pd.to_numeric(frame['contig_len'], errors='coerce')
    if frame['contig_len'].isna().any() or (frame['contig_len'] < 0).any():
        raise ValueError('Missing or negative query lengths; cannot calculate composition')
    frame['genospecies'] = frame['genospecies'].fillna(UNCLASSIFIED)
    composition = frame.groupby(['assembly_id', 'genospecies'], as_index=False).agg(
        contigs=('contig_id', 'size'), assigned_bp=('contig_len', 'sum'))
    totals = composition.groupby('assembly_id')['assigned_bp'].transform('sum')
    composition['total_query_bp'] = totals
    composition['percent_composition'] = (100 * composition['assigned_bp'] / totals.where(totals > 0)).fillna(0)

    composition = composition.sort_values(
        ["assembly_id", "percent_composition", "genospecies"],
        ascending=[True, False, True],
    )
    
    composition.to_csv(output, sep='\t', index=False, float_format='%.6f')
    print(f"Output genospecies composition to {output}")
    
    return composition
