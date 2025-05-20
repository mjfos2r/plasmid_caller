# scoring.py
import pandas

### scoring parameters ###
PF32_MIN_BP         = 200
WP_MIN_COV_PCT      = 10
WP_MIN_COV_BP       = 1000
WP_OVERRIDE_COV_PCT = 90
WP_OVERRIDE_PID_PCT = 95
##########################


def best_pf32_hit(df: pandas.DataFrame) -> pandas.Dataframe:
    """Make sure the covered length is over 200, then sort by percent identity and choose the highest."""
    df = df[df["query_covered_length"] >= PF32_MIN_BP]
    if df.empty:
        return df
    return (
        df.sort_values(["overall_percent_identity", "e_value"], ascending=[False, True]).groupby("contig_id", as_index=False).first()
    )

def best_wp_hit(df: pandas.DataFrame) -> pandas.Dataframe:
    """score each wp hit by the product of the percent coverage and the percent identity. (as long as coverage is either over 10% or greater than 1000bp) """
    # todo add modular params to tweak these cutoffs # done
    df = df[df
            (df["query_coverage_percent"] >= WP_MIN_COV_PCT) |
            (df["query_covered_length"] >= WP_MIN_COV_BP) ]
    if df.empty:
        return df
    df["_score"] = (df["query_coverage_percent"] *
                    df["overall_percent_identity"])
    best = (df.sort_values(["_score", "e_value"], ascending=[False, True]).groupby("contig_id", as_index=False).first())
    return best


def _valid_pf32(row):
    return row["query_covered_length_pf32"] >= PF32_MIN_BP


def _valid_wp(row):
    return (
        row["query_covered_percent_wp"] >= WP_MIN_COV_PCT
        or row["query_covered_length_wp"] >= WP_MIN_COV_BP
    )


def choose_final_call(row):
    """after determining the best pf32 and wp hits, choose the pf32 hit first (if present) and only use the wp hit if it's very convincing. (set params above)"""
    if pandas.isna(row.get("plasmid_name_pf32")) and pandas.isna(row.get("plasmid_name_wp")):
        return "unclassified"

    pf32_ok = not pandas.isna(row.get("plasmid_name_pf32")) and _valid_pf32(row)
    wp_ok = not pandas.isna(row.get("plasmid_name_wp")) and _valid_wp(row)

    if pf32_ok and not wp_ok:
        return row["plasmid_name_pf32"]

    if wp_ok and not pf32_ok:
        return row["plasmid_name_wp"]

    if pf32_ok and wp_ok:
        if (row["query_coverage_percent_wp"] >= WP_OVERRIDE_COV_PCT and row["overall_percent_identity_wp"] >= WP_OVERRIDE_COV_PCT):
            return row["plasmid_name_wp"]
        return row["plasmid_name_pf32"]

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
