import gzip
import io
import sys
import time

import pytest

fastapi = pytest.importorskip("fastapi")
from fastapi.testclient import TestClient

from plasmid_caller.web.app import create_app
from plasmid_caller.web.jobs import JobStore, Settings, safe_filename, validate_fasta

FASTA = b">contig_1\nACGTACGTAC\n>contig_2\nACGT\n"
FAKE_CLI = """
import pathlib, sys
out = pathlib.Path(sys.argv[1]); out.mkdir(parents=True)
(out / "summary_best_hits.tsv").write_text("contig_id\\tcontig_len\\tfinal_call\\nc1\\t10\\tlp17\\n")
(out / "summary_all_calls.tsv").write_text("contig_id\\nc1\\n")
(out / "genospecies_composition.tsv").write_text("genospecies\\tpercent_composition\\nafzelii\\t100\\n")
(out / "x_renamed.fasta").write_text(">c\\nA\\n")
xml = out / "wp" / "xml_files"; xml.mkdir(parents=True); (xml / "a.xml").write_text("x")
"""


def fake_command(script):
    def command(self, input_path, output):
        return [sys.executable, "-c", script, str(output)]
    return command


def client_for(tmp_path, monkeypatch, script=FAKE_CLI, **settings):
    monkeypatch.setattr(JobStore, "_command", fake_command(script))
    app = create_app(Settings(runs_dir=tmp_path / "runs", **settings))
    return TestClient(app, base_url="http://localhost")


def upload(client, name="a.fna", content=FASTA):
    return client.post("/api/batches", files=[("files", (name, io.BytesIO(content)))])


def wait_for(client, job_id):
    for _ in range(200):
        job = client.get(f"/api/jobs/{job_id}").json()
        if job["status"] not in {"queued", "running"}:
            return job
        time.sleep(0.05)
    raise AssertionError("job did not finish")


@pytest.mark.parametrize("name, expected", [
    ("C41.fna", "C41.fna"),
    ("../../etc/passwd.fna", "passwd.fna"),
    ("..\\..\\evil.FASTA.GZ", "evil.fasta.gz"),
    ("my assembly (1).fa", "my_assembly__1.fa"),
    ("...fna", "assembly.fna"),
])
def test_safe_filename(name, expected):
    assert safe_filename(name) == expected


@pytest.mark.parametrize("name, content, message", [
    ("a.txt", FASTA, "Expected"),
    ("a.fna", b"", "No FASTA"),
    ("a.fna", b"not a fasta\n", "No FASTA"),
    ("a.fna", b">a\nACGT\n>a\nACGT\n", "Repeated"),
    ("a.fna", b">a\nMKVLAAGIVG\n", "Non-nucleotide"),
    ("a.fna.gz", b"not gzip", "could not be read"),
])
def test_validate_fasta_rejects(tmp_path, name, content, message):
    path = tmp_path / name
    path.write_bytes(content)
    with pytest.raises(ValueError, match=message):
        validate_fasta(path, 1_000_000)


def test_validate_fasta_stops_gzip_bombs(tmp_path):
    path = tmp_path / "bomb.fna.gz"
    with gzip.open(path, "wb") as handle:
        handle.write(b">a\n" + b"A" * 5_000_000)
    with pytest.raises(ValueError, match="exceeds"):
        validate_fasta(path, 1_000_000)


def test_upload_run_results_and_cleanup(tmp_path, monkeypatch):
    with client_for(tmp_path, monkeypatch) as client:
        batch = upload(client).json()
        job = wait_for(client, batch["jobs"][0]["id"])
        assert job["status"] == "done" and job["contigs"] == 2 and job["bases"] == 14
        results = client.get(f"/api/jobs/{job['id']}/results").json()
        assert results["contigs"][0]["final_call"] == "lp17"
        assert results["contigs"][0]["genospecies"] is None
        assert client.get(f"/api/jobs/{job['id']}/download/renamed.fasta").status_code == 200
        assert not list((tmp_path / "runs").rglob("xml_files"))


def test_download_only_serves_known_files(tmp_path, monkeypatch):
    with client_for(tmp_path, monkeypatch) as client:
        job = wait_for(client, upload(client).json()["jobs"][0]["id"])
        for name in ["job.json", "..%2Fjob.json", "log.txt"]:
            assert client.get(f"/api/jobs/{job['id']}/download/{name}").status_code == 404
        assert client.get("/api/jobs/../batch.json").status_code == 404


def test_failed_pipeline_is_reported_without_details(tmp_path, monkeypatch):
    with client_for(tmp_path, monkeypatch, script="import sys; sys.exit(1)") as client:
        job = wait_for(client, upload(client).json()["jobs"][0]["id"])
        assert job["status"] == "failed" and "/" not in job["error"].replace("log.txt", "")
        assert client.get(f"/api/jobs/{job['id']}/results").status_code == 409


def test_exit_zero_without_outputs_is_a_failure(tmp_path, monkeypatch):
    with client_for(tmp_path, monkeypatch, script="pass") as client:
        job = wait_for(client, upload(client).json()["jobs"][0]["id"])
        assert job["status"] == "failed"


def test_rejected_upload_keeps_rest_of_batch(tmp_path, monkeypatch):
    with client_for(tmp_path, monkeypatch) as client:
        response = client.post("/api/batches", files=[
            ("files", ("good.fna", io.BytesIO(FASTA))),
            ("files", ("bad.fna", io.BytesIO(b"nope"))),
        ])
        jobs = {job["name"]: job for job in response.json()["jobs"]}
        assert jobs["bad.fna"]["status"] == "rejected"
        assert wait_for(client, jobs["good.fna"]["id"])["status"] == "done"


def test_runs_survive_restart_and_unfinished_jobs_requeue(tmp_path, monkeypatch):
    with client_for(tmp_path, monkeypatch) as client:
        job = wait_for(client, upload(client).json()["jobs"][0]["id"])
    path = next((tmp_path / "runs").rglob("job.json"))
    path.write_text(path.read_text().replace('"done"', '"running"'))
    with client_for(tmp_path, monkeypatch) as client:
        assert wait_for(client, job["id"])["status"] == "done"
        assert len(client.get("/api/batches").json()) == 1


def test_limits_and_origin_checks(tmp_path, monkeypatch):
    with client_for(tmp_path, monkeypatch, max_upload_mb=1, max_queued=1) as client:
        assert upload(client, content=b"A" * 2_000_000).status_code == 413
        files = [("files", (f"{i}.fna", io.BytesIO(FASTA))) for i in range(2)]
        assert client.post("/api/batches", files=files).status_code == 429
        evil = {"origin": "http://evil.example"}
        assert client.post("/api/batches", files=files, headers=evil).status_code == 403
        assert client.get("/", headers={"host": "evil.example"}).status_code == 403


def test_token_required_when_configured(tmp_path, monkeypatch):
    with client_for(tmp_path, monkeypatch, token="secret") as client:
        assert client.get("/api/batches").status_code == 401
        assert client.get("/?token=wrong").status_code == 401
        assert client.get("/?token=secret", follow_redirects=False).status_code == 303
        assert client.get("/api/batches").status_code == 200
