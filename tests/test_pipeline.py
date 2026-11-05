import gzip

import pytest

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
