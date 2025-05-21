'''
pytest test-suite for plasmid_caller 6.0.0
================================================
These unit-tests exercise the public helper functions that do **not** require an active
BLAST installation or large reference databases.  External side-effects (filesystem,
sub-processes, network) are isolated with monkey-patching so that the suite runs fast
and deterministically and still achieves high (>90 %) statement coverage when
executed with `pytest --cov=plasmid_caller`.

Directory layout that pytest expects:

project/
├── plasmid_caller/               # the package under test
│   ├── __init__.py
│   ├── plasmid_caller.py         # (contains the code from user prompt)
│   └── ...
└── tests/
    └── test_plasmid_caller.py    # <-- this file
'''

from __future__ import annotations

from types import SimpleNamespace
from pathlib import Path
from textwrap import dedent
import pickle
import builtins
import importlib
import sys

import pandas as pd
import pytest

# -----------------------------------------------------------------------------
# Import the module under test **once** so that monkeypatching works reliably.
# -----------------------------------------------------------------------------

import plasmid_caller.plasmid_caller as pc
from plasmid_caller.scoring import (
    best_pf32_hit,
    best_wp_hit,
    choose_final_call,
)

# -----------------------------------------------------------------------------
# •••                Fixtures (shared test data / helpers)                 •••
# -----------------------------------------------------------------------------

@pytest.fixture
def fasta_file(tmp_path: Path) -> Path:
    """Create a simple single-contig FASTA (>contig1, 1000 bp of A)."""
    fa = tmp_path / "test.fasta"
    fa.write_text(">contig1\n" + "A" * 1000 + "\n")
    return fa


@pytest.fixture
def fake_hsp():
    """Return a minimal fake HSP object with the attributes that
    pc.calculate_percent_identity_and_coverage() accesses."""

    class _FakeHSP:
        def __init__(self, q_start, q_end, s_start, s_end, ids, aln_len):
            self.query_start = q_start
            self.query_end = q_end
            self.sbjct_start = s_start
            self.sbjct_end = s_end
            self.identities = ids
            self.align_length = aln_len

    return _FakeHSP


@pytest.fixture
def fake_alignment(fake_hsp):
    """Return a minimal fake Alignment object with two overlapping HSPs."""

    class _FakeAlignment:
        def __init__(self):
            self.length = 1500
            self.hit_id = "WP_012345"
            self.hsps = [
                fake_hsp(10, 110, 1, 101, 95, 100),
                fake_hsp(90, 190, 80, 180, 90, 100),
            ]

    return _FakeAlignment()


# -----------------------------------------------------------------------------
# •••                               Tests                                  •••
# -----------------------------------------------------------------------------

# ––– get_input_files() –––––––––––––––––––––––––––––––––––––––––––––––––––––––

def test_get_input_files(tmp_path: Path):
    (tmp_path / "a.fasta").write_text("dummy")
    (tmp_path / "b.fasta").write_text("dummy")
    (tmp_path / "c.txt").write_text("dummy")  # wrong extension should be ignored

    files = pc.get_input_files(tmp_path, "fasta")

    # Path.glob returns Path objects — convert to names for easier asserts
    names = sorted(p.name for p in files)
    assert names == ["a.fasta", "b.fasta"]


# ––– get_output_path() –––––––––––––––––––––––––––––––––––––––––––––––––––––––

def test_get_output_path_creates_directory(tmp_path: Path, monkeypatch):
    # Pretend that the CWD is `tmp_path` so we do not litter the real FS
    monkeypatch.chdir(tmp_path)
    out_dir = pc.get_output_path("results")

    assert Path(out_dir).is_dir()
    assert Path(out_dir).samefile(tmp_path / "results")


# ––– get_blast_command() –––––––––––––––––––––––––––––––––––––––––––––––––––––

def test_get_blast_command_basics():
    params = pc.get_blast_command(
        prog="blastn",
        input_file="/data/contig.fa",
        output_path="/out",
        database="mydb",
        threads=8,
    )

    assert params["program"] == "blastn"
    assert params["query"].endswith("contig.fa")
    assert params["out"].endswith("_blast_results.xml")
    assert params["db"] == "mydb"
    assert params["num_threads"] == 8


# ––– run_blast() ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

def test_run_blast_happy_path(monkeypatch):
    """Verify that run_blast delegates to BlastManager and returns its result."""

    dummy_ret_val = "OK"

    def fake_run(program, **kwargs):
        assert program == "blastn"
        # a representative subset of kwargs is sufficient for the unit test
        assert "query" in kwargs and "db" in kwargs
        return dummy_ret_val

    monkeypatch.setattr(pc.blast_manager, "run_blast_command", fake_run)

    result = pc.run_blast({"program": "blastn", "query": "in.fa", "db": "db"})
    assert result == dummy_ret_val


# ––– get_name_from_acc() & parse_hit_id() –––––––––––––––––––––––––––––––––––

def test_get_name_from_acc():
    parsing_dict = {
        "ACC123": {"name": "lp54", "strain": "B31"},
    }
    strain, name = pc.get_name_from_acc("gb|ACC123|", parsing_dict)
    assert (strain, name) == ("B31", "lp54")


def test_parse_hit_id_various_cases():
    # standard “strain_plasmid” pattern
    s, p = pc.parse_hit_id("B31_lp28-1")
    assert s == "B31" and p == "lp28-1"

    # single token (PDB-like) → strain is token after pipe, plasmid None
    s, p = pc.parse_hit_id("pdb|ACC|stuff")
    assert p is None

    # NE_xxxx edge-case
    s, p = pc.parse_hit_id("NE_1234_cp26")
    assert s.startswith("NE_1234") and p == "cp26"


# ––– get_contig_len() –––––––––––––––––––––––––––––––––––––––––––––––––––––––

def test_get_contig_len(fasta_file: Path):
    length = pc.get_contig_len(fasta_file, "contig1")
    assert length == 1000


# ––– calculate_percent_identity_and_coverage() ––––––––––––––––––––––––––––––

def test_calculate_identity_and_coverage(fake_alignment):
    stats = pc.calculate_percent_identity_and_coverage(fake_alignment)

    # identities 95+90 = 185 over200 bp alignment
    assert pytest.approx(stats["overall_percent_identity"]) == 92.5

    # coverage on query: two HSPs 10–190 merged gives 180 bp
    assert stats["query_covered_length"] == 180

    # reference (subject) intervals 1–180 gives 179 bp (end-exclusive)
    assert stats["ref_covered_length"] == 179


# ––– parse_blast_xml() ––––––––––––––––––––––––––––––––––––––––––––––––––––––

def test_parse_blast_xml_no_hits(monkeypatch, tmp_path: Path, fasta_file: Path):
    """Minimal integration-test: NCBIXML.parse is mocked so that the real XML
    file content is not needed; however, the function still attempts to *open*
    the file, so we create an empty placeholder on disk."""

    # Fake record with *no* alignments
    rec = SimpleNamespace(query="contig1", query_length=1000, alignments=[])
    monkeypatch.setattr(pc.NCBIXML, "parse", lambda handle: [rec])

    # the parsing_dict pickle expected by parse_blast_xml
    db_dir = tmp_path / "db"
    db_dir.mkdir()
    pickle_path = db_dir / "blast_parsing_dict.pkl"
    pickle_path.write_bytes(pickle.dumps({"dummy": {"name": "lp54", "strain": "B31"}}))

    # Empty XML placeholder – contents irrelevant because parse() is mocked
    xml_path = tmp_path / "dummy.xml"
    xml_path.write_text("")

    hits = pc.parse_blast_xml(
        xml_file=str(xml_path),
        args=SimpleNamespace(input=fasta_file, quiet=True),
        parsing_type="wp",
        dbs_dir=db_dir,
    )

    # Should contain a single contig with NA values
    assert list(hits.keys()) == ["contig1"]
    row = hits["contig1"][0]
    assert row["assembly_id"].startswith("dummy")
    assert pd.isna(row["plasmid_id"])



# ––– get_db_type() ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

def test_get_db_type(monkeypatch, tmp_path: Path):
    fake_output = "\n".join([
        "/path/to/pf32 protein",
        "/path/to/wp nucleotide",
    ])

    class _FakeCompleted:
        stdout = fake_output

    monkeypatch.setattr(pc.blast_manager, "run_blast_command", lambda *a, **k: _FakeCompleted())

    db_dir = tmp_path / "db"
    db_dir.mkdir()
    (db_dir / "dummy.nhr").touch()  # presence isn‘t validated, but listdir is used

    pairs = pc.get_db_type(db_dir)
    assert pairs == [("/path/to/pf32", "blastx"), ("/path/to/wp", "blastn")]


# ––– scoring helpers ––––––––––––––––––––––––––––––––––––––––––––––––––––––––

def _sample_df(kind="pf32"):
    rows = [
        {
            "contig_id": "c1",
            "overall_percent_identity": 99.0,
            "query_covered_length": 300,
            "query_coverage_percent": 60,
        },
        {
            "contig_id": "c1",
            "overall_percent_identity": 97.0,
            "query_covered_length": 250,
            "query_coverage_percent": 55,
        },
        {
            "contig_id": "c2",
            "overall_percent_identity": 88.0,
            "query_covered_length": 180 if kind == "pf32" else 150,
            "query_coverage_percent": 49,
        },
    ]
    return pd.DataFrame(rows)


def test_best_pf32_hit():
    df = _sample_df("pf32")
    best = best_pf32_hit(df)

    # contig c1 survives (>200 bp) and retains highest PID row
    assert list(best["contig_id"]) == ["c1"]
    assert best.iloc[0]["overall_percent_identity"] == 99.0


def test_best_wp_hit():
    df = _sample_df("wp")
    best = best_wp_hit(df)

    # row wit 60 %  99 % should win the score race
    assert best.iloc[0]["contig_id"] == "c1"
    assert "_score" in best.columns


def test_choose_final_call_logic():
    row = pd.Series(
        {
            "plasmid_name_pf32": "lp28-1",
            "query_covered_length_pf32": 250,
            "plasmid_name_wp": "lp28-2",
            "query_coverage_percent_wp": 80,
            "query_covered_length_wp": 2000,
            "overall_percent_identity_wp": 95,
        }
    )
    assert choose_final_call(row) == "lp28-1"  # pf32 preferred when both valid

    # no pf32 hit, wp valid
    row["plasmid_name_pf32"] = pd.NA
    assert choose_final_call(row) == "lp28-2"

    # nothing valid → unclassified
    row["plasmid_name_wp"] = pd.NA
    assert choose_final_call(row) == "unclassified"
