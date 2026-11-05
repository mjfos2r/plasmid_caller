"""Replicon name handling: which reference names are fusions, and of what.

Fusion names vary across sources, so every spelling is resolved here:
cp32-1+3, cp32-5-1, cp32-3-lp25, cp32-3lp25, lp28-4_cp32-1, lp36+32-4, cp32-28-4.
"""
import re

import pandas as pd

from .scoring import _parse_intervals

FUSION_SPLIT = re.compile(r"\+|[-_]?(?=lp|cp)")
CP32_PAIR = re.compile(r"cp32-(\d+)-(\d+)")
TILE_MIN_IDENTITY = 90
TILE_MIN_GAIN_BP = 3_000
# The longest fused reference is 112 kb; anything far larger is a chromosome.
TILE_MAX_CONTIG_BP = 200_000


def replicon_parts(name):
    """Return the replicons a name is made of; an unfused name returns itself."""
    parts = []
    for piece in filter(None, FUSION_SPLIT.split(name)):
        pair = CP32_PAIR.fullmatch(piece)
        if piece in {"28-4", "cp32-28-4"}:
            parts.append("lp28-4")
        elif piece.isdigit():
            parts.append(f"cp32-{piece}")
        elif re.fullmatch(r"32-\d+", piece):
            parts.append(f"cp{piece}")
        elif pair:
            parts += [f"cp32-{pair[1]}", f"cp32-{pair[2]}"]
        else:
            parts.append(piece)
    return list(dict.fromkeys(parts))


def is_fusion(name):
    return replicon_parts(name) != [name]


def _merge(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


def _size(intervals):
    return sum(end - start + 1 for start, end in intervals)


def tile_contig(hits):
    """Greedily cover one contig with unfused plasmid references.

    Returns [(name, new_bp), ...] in pick order. Two or more entries mean no single
    plasmid explains the contig. Names are reliable for lp+cp pairs only: cp32s are
    similar enough that one can stand in for another.
    """
    names = hits["plasmid_name"].astype(str)
    usable = hits[~names.map(is_fusion) & (names != "chromosome") & hits["plasmid_id"].notna()
                  & (hits["overall_percent_identity"] >= TILE_MIN_IDENTITY)]
    remaining = {}
    for name, group in usable.groupby("plasmid_name"):
        best = group.sort_values("query_covered_length", ascending=False).iloc[0]
        remaining[name] = _merge(tuple(pair) for pair in _parse_intervals(best["query_intervals"]))
    covered, picks = [], []
    while remaining:
        gains = {name: _size(_merge(covered + spans)) - _size(covered)
                 for name, spans in remaining.items()}
        name = max(gains, key=gains.get)
        if gains[name] < TILE_MIN_GAIN_BP:
            break
        covered = _merge(covered + remaining.pop(name))
        picks.append((name, gains[name]))
    return picks


def fusion_candidates(hits):
    """Map contig_id to its tiling parts, for contigs where two or more plasmids are needed."""
    hits = hits[pd.to_numeric(hits["contig_len"], errors="coerce") < TILE_MAX_CONTIG_BP]
    tiles = {contig: tile_contig(group) for contig, group in hits.groupby("contig_id")}
    return {contig: picks for contig, picks in tiles.items() if len(picks) >= 2}
