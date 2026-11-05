import gzip

import pytest

from plasmid_caller.scoring import decide_call
from plasmid_caller import main as pc
from plasmid_caller.scoring import _parse_intervals


def make_db(root, *indexes):
    for index in indexes:
        path = root / index.split(".")[0] / index
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()


def test_get_db_type_resolves_exactly_two_and_ignores_nested_copies(tmp_path):
    make_db(tmp_path, "pf32.psq", "wp.nsq")
    (tmp_path / "wp" / "wp").mkdir()
    (tmp_path / "wp" / "wp" / "wp.nsq").touch()
    assert pc.get_db_type(tmp_path) == [
        (str(tmp_path.resolve() / "pf32" / "pf32"), "blastx"),
        (str(tmp_path.resolve() / "wp" / "wp"), "blastn"),
    ]


def test_get_db_type_exits_when_a_database_is_missing(tmp_path):
    make_db(tmp_path, "wp.nsq")
    with pytest.raises(SystemExit, match="pf32.psq"):
        pc.get_db_type(tmp_path)


def test_run_blast_raises_on_failure(tmp_path):
    params = pc.get_blast_command("blastn", "x", tmp_path / "absent.fna", tmp_path,
                                  tmp_path / "absent_db", 1, 5)
    with pytest.raises(RuntimeError, match="BLAST command failed"):
        pc.run_blast(params)
    assert params["program"] == "blastn"


@pytest.mark.parametrize("name, expected", [
    ("C41.fna", "C41"),
    ("C41.FASTA", "C41"),
    ("GCF_000008685.2_ASM868v2_genomic.fna.gz", "GCF_000008685.2_ASM868v2_genomic"),
    ("sample.fa.gz", "sample"),
])
def test_get_file_id(name, expected):
    assert pc.get_file_id(name) == expected


def test_get_contig_lengths_plain_and_gzipped(tmp_path):
    fasta = ">a desc\nACGT\nAC\n>b\nA\n"
    plain = tmp_path / "x.fna"
    plain.write_text(fasta)
    packed = tmp_path / "x.fna.gz"
    with gzip.open(packed, "wt") as handle:
        handle.write(fasta)
    assert pc.get_contig_lengths(plain) == {"a": 6, "b": 1}
    assert pc.get_contig_lengths(packed) == {"a": 6, "b": 1}


def test_parse_intervals_reads_tsv_strings():
    assert _parse_intervals("[(1, 246), (300, 400)]") == [(1, 246), (300, 400)]


def call_row(pf32, wp, length=24_000):
    return {"plasmid_name_pf32": pf32, "plasmid_name_wp": wp, "contig_len": length,
            "query_covered_length_pf32": 785, "query_coverage_percent_wp": 99,
            "query_covered_length_wp": length, "multiple_loci_pf32": False}


@pytest.mark.parametrize("pf32, wp, length, expected", [
    ("chromosome", "lp21-cp9", 24_000, "lp21-cp9"),
    ("chromosome", "chromosome", 24_000, "chromosome"),
    ("lp21", "lp21-cp9", 24_000, "lp21"),
    ("chromosome", "lp21", 500_000, "chromosome"),
])
def test_chromosomal_para_does_not_override_a_plasmid_hit(pf32, wp, length, expected):
    assert decide_call(call_row(pf32, wp, length))[0] == expected


@pytest.mark.parametrize("hit_id, strain, name", [
    ("HO14_NovelPFam32_DatabaseShortContig_ParA_X", "HO14", "unknown"),
    ("RS01295_JD1_JD1_cp32-5_ParA_X", "RS01295", "cp32-5"),
    ("NE_5261_chromosome_ParA_2", "NE_5261", "chromosome"),
])
def test_parse_hit_id(hit_id, strain, name):
    assert pc.parse_hit_id(hit_id) == (strain, name)
