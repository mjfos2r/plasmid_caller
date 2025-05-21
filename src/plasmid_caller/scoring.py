# scoring.py
from ctypes.macholib import dyld
from typing import Any
import pandas

#### scoring parameters #####
PF32_MIN_BP         = 200
WP_MIN_COV_PCT      = 50
WP_MIN_COV_BP       = 1_000
WP_OVERRIDE_COV_PCT = 95
WP_OVERRIDE_PID_PCT = 98
#### Chromosomal Params #####
MIN_CALL_BP         = 1_000
CHROMOSOME_MIN_BP   = 100_000
#############################
# TODO: Catch divergent cases.(where pf32's best hit has stats below the min thresholds)

def best_pf32_hit(df: pandas.DataFrame) -> pandas.DataFrame:
    """Make sure the covered length is over 200, then sort by percent identity and choose the highest."""
    df = df[df["query_covered_length"] >= PF32_MIN_BP]
    if df.empty:
        return df
    return (
        df.sort_values(["overall_percent_identity"], ascending=False).groupby("contig_id", as_index=False).first()
    )

def best_wp_hit(df: pandas.DataFrame) -> pandas.DataFrame:
    """score each wp hit by the product of the percent coverage and the percent identity. (as long as coverage is either over 10% or greater than 1000bp) """
    # todo add modular params to tweak these cutoffs # done
    df = df[
            (df["query_coverage_percent"] >= WP_MIN_COV_PCT) |
            (df["query_covered_length"] >= WP_MIN_COV_BP) ]
    if df.empty:
        return df
    df["_score"] = (df["query_coverage_percent"] *
                    df["overall_percent_identity"])
    best = (df.sort_values(["_score"], ascending=False).groupby("contig_id", as_index=False).first())
    return best


def _valid_pf32(row):
    return row["query_covered_length_pf32"] >= PF32_MIN_BP


def _valid_wp(row):
    return (
        row["query_coverage_percent_wp"] >= WP_MIN_COV_PCT
        or row["query_covered_length_wp"] >= WP_MIN_COV_BP
    )

def _parse_intervals(raw):
    """
    this accepts in whatever value is stored at row["query_intervals"] and
    evaluates it to the type we need for intervaltree using ast.literal_eval
    """
    if isinstance(raw, (list, tuple)):
        return raw
    if isinstance(raw, str):
        try:
            val = ast.literal_eval(raw)
            if isinstance(val, (list, tuple)):
                return val
        except Exception:
            pass
    return []

from intervaltree import Interval, IntervalTree
import pandas as pd
from .scoring import best_pf32_hit, PF32_MIN_BP


def summarise_multiple_pf32(full_pf32_df: pandas.DataFrame, min_cov_bp: int = PF32_MIN_BP):
    """
    For each contig:
      • filter for PF32 hits passing min_cov_bp
      • build an IntervalTree of its query_intervals
      • merge overlaps → distinct loci
      • for each locus, pick the best hit via best_pf32_hit()
      • collect that best-hit row + locus_idx/start/end

    Returns:
      - multi_best_df: DataFrame with one row per locus,
                       original columns + locus_idx, locus_start, locus_end
      - concat_series: per-contig concatenation of those locus' plasmid_names in descending PID order
    """
    records = []

    # 1 QC filter up front
    qual = full_pf32_df[
        (full_pf32_df["query_covered_length"] >= min_cov_bp) &
        full_pf32_df["query_intervals"].notna()
    ]

    # 2 process contig by contig
    for contig_id, contig_df in qual.groupby("contig_id"):
        contig_df = contig_df.reset_index(drop=True)
        tree = IntervalTree()
        # add each hit's intervals carrying its row-index
        for ix, val in enumerate(contig_df["query_intervals"]):
            ivals = _parse_intervals(val)
            for start, end in ivals:
                tree.add(Interval(int(start), int(end), {ix}))
        # merge overlaps for this contig only
        tree.merge_overlaps(data_reducer=lambda a, b: a.union(b))

        if len(tree) < 2:
            locus_idx = 1
            continue

        # 3 pick best hit per locus
        for locus_idx, iv in enumerate(sorted(tree), start=1):
            hit_idxs = iv.data
            locus_subset = contig_df.loc[list(hit_idxs)]
            best = best_pf32_hit(locus_subset)
            if best.empty:
                continue
            row = best.iloc[0].to_dict()
            # attach locus metadata
            row.update({
                "locus_idx":    locus_idx,
                "locus_start":  iv.begin,
                "locus_end":    iv.end
            })
            records.append(row)

    # assemble the per-locus DataFrame
    if records:
        multi_best_df = pandas.DataFrame.from_records(records)
    else:
        cols = list(full_pf32_df.columns) + ["locus_idx", "locus_start", "locus_end"]
        multi_best_df = pandas.DataFrame(columns=cols)

    def _concat_pf32(names):
        uniq = pandas.unique(names)
        if len(uniq) > 1:
            return ":::".join(uniq)
        return f"{uniq[0]}*"

    # 4 concatenated call series
    if not multi_best_df.empty:
        concat = (
            multi_best_df
            .sort_values(["contig_id", "overall_percent_identity"], ascending=[True, False])
            .groupby("contig_id")["plasmid_name"]
            .apply(_concat_pf32)
            .rename("concat_call")
        )
    else:
        concat = pandas.Series(dtype="object", name="concat_call")

    return multi_best_df, concat

def choose_final_call(row):
    """after determining the best pf32 and wp hits, choose the pf32 hit first (if present) and only use the wp hit if it's very convincing. (set params above)"""
    if pandas.isna(row.get("plasmid_name_pf32")) and pandas.isna(row.get("plasmid_name_wp")):
        return "unclassified"

    # if contig_len is less than MIN_CALL_BP, we don't want to call it since it's obviously a fragment.
    if row["contig_len"] < MIN_CALL_BP:
        return "unclassified"

    pf32_ok = not pandas.isna(row.get("plasmid_name_pf32")) and _valid_pf32(row)
    wp_ok = not pandas.isna(row.get("plasmid_name_wp")) and _valid_wp(row)

    # check for really long contigs, likely chromosomal fragments.
    # in prior versions, anything over 100,000 was automatically classified as chromosome.
    # let's flesh out this logic a bit more.
    if row["contig_len"] >= CHROMOSOME_MIN_BP:
        # if we have a chromosome call, return that.
        if (pf32_ok and row["plasmid_name_pf32"].lower() == "chromosome") or \
            (wp_ok and row["plasmid_name_wp"].lower() == "chromosome"):
                return "chromosome"
        else:
            return row["plasmid_name_wp"]
        ## allow a really good wp_hit to override. (this probably won't happen. but could in cases of tremendous misassembly)
        #if wp_ok and (row["query_coverage_percent_wp"] >= WP_OVERRIDE_COV_PCT and
        #              row["overall_percent_identity_wp"] >= WP_OVERRIDE_PID_PCT):
        #    return row["plasmid_name_wp"]

    if row.get("multiple_loci_pf32", False):
        return row["concat_call_pf32"]

    if pf32_ok and not wp_ok:
        return row["plasmid_name_pf32"]

    if wp_ok and not pf32_ok:
        return row["plasmid_name_wp"]

    # identity is primarily based on pf32, I don't feel like dealing with the nuance of this right now.
    if pf32_ok and wp_ok:
        return row["plasmid_name_pf32"]
    #if pf32_ok and wp_ok:
    #    if (row["query_coverage_percent_wp"] >= WP_OVERRIDE_COV_PCT and row["overall_percent_identity_wp"] >= WP_OVERRIDE_COV_PCT):
    #        return row["plasmid_name_wp"]
    #    return row["plasmid_name_pf32"]

    return "unclassified"

# unused function to determine the best match
def get_best_match(matches, key):
    """Highest percent identity takes the cake. specify which feature to compare.
    ex: get_best_match(matches, "percent_identity")"""
    best_match = None
    best_score = -1
    for match in matches:
        if match[key] > best_score:
            best_score = match[key]
            best_match = match
    return best_match

# take in a dataframe of results, choose the best one.
