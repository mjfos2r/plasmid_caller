import pandas as pd
import pytest

from plasmid_caller.replicons import fusion_candidates, is_fusion, replicon_parts, tile_contig


@pytest.mark.parametrize("name, parts", [
    ("cp32-1+3", ["cp32-1", "cp32-3"]),
    ("cp32-4+7+28-4", ["cp32-4", "cp32-7", "lp28-4"]),
    ("cp32-5-1", ["cp32-5", "cp32-1"]),
    ("cp32-9-4", ["cp32-9", "cp32-4"]),
    ("cp32-3-lp25", ["cp32-3", "lp25"]),
    ("cp32-3lp25", ["cp32-3", "lp25"]),
    ("lp28-4cp32-1", ["lp28-4", "cp32-1"]),
    ("lp28-4_cp32-1", ["lp28-4", "cp32-1"]),
    ("lp21-cp9", ["lp21", "cp9"]),
    ("lp36+32-4", ["lp36", "cp32-4"]),
    ("cp32-28-4", ["lp28-4"]),
])
def test_fusion_spellings(name, parts):
    assert replicon_parts(name) == parts and is_fusion(name)


@pytest.mark.parametrize("name", ["chromosome", "cp9", "cp9-1", "cp26", "cp32-1", "cp32-13", "lp28-4",
                                  "lp28-14", "lp29-9", "lp32-10", "lp54"])
def test_unfused_names(name):
    assert replicon_parts(name) == [name] and not is_fusion(name)


def hits(*rows):
    return pd.DataFrame([
        {"contig_id": "c1", "contig_len": 40000, "plasmid_id": "gb|A|", "plasmid_name": name,
         "overall_percent_identity": identity, "query_covered_length": end - start + 1,
         "query_intervals": str([(start, end)])}
        for name, identity, start, end in rows])


def test_two_plasmids_tile_a_contig():
    frame = hits(("lp21", 99, 1, 20000), ("cp9", 98, 20001, 30000), ("lp17", 99, 1, 19000))
    assert tile_contig(frame) == [("lp21", 20000), ("cp9", 10000)]
    assert list(fusion_candidates(frame)) == ["c1"]


def test_fused_chromosome_and_divergent_hits_do_not_tile():
    frame = hits(("lp21", 99, 1, 20000), ("lp21-cp9", 99, 1, 30000), ("chromosome", 99, 20001, 30000),
                 ("cp9", 80, 20001, 30000), ("lp17", 99, 18001, 21000 - 1))
    assert tile_contig(frame) == [("lp21", 20000)]
    assert fusion_candidates(frame) == {}
