import gzip
import io
import json
import shutil
import sys
import time
from types import SimpleNamespace

import pandas as pd
import pytest

fastapi = pytest.importorskip("fastapi")
from fastapi.testclient import TestClient

from plasmid_caller.web import auth
from plasmid_caller.web.app import create_app, extract_contig_fasta, guard_cell, identity_cells, replicon_totals
from plasmid_caller.web.jobs import JobStore, Settings, safe_filename, validate_fasta

jwt = pytest.importorskip("jwt")
from cryptography.hazmat.primitives.asymmetric import rsa

ACCESS_KEY = rsa.generate_private_key(public_exponent=65537, key_size=2048)
FASTA = b">contig_1\nACGTACGTAC\n>contig_2\nACGT\n"
FAKE_CLI = """
import pathlib, sys
out = pathlib.Path(sys.argv[1]); out.mkdir(parents=True)
(out / "summary_best_hits.tsv").write_text(
    "contig_id\\tcontig_len\\tfinal_call\\tfiltered_call\\n001\\t10\\tlp17\\tlp17\\n")
(out / "summary_all_calls.tsv").write_text("contig_id\\tcontig_len\\tfinal_call\\n001\\t10\\tunclassified\\n")
(out / "genospecies_composition.tsv").write_text("genospecies\\tpercent_composition\\nafzelii\\t100\\n")
(out / "x_renamed.fasta").write_text(">c\\nA\\n")
xml = out / "wp" / "xml_files"; xml.mkdir(parents=True); (xml / "a.xml").write_text("x")
pf32 = out / "pf32" / "tables"; pf32.mkdir(parents=True)
(pf32 / "a_multi_loci.tsv").write_text(
    "contig_id\\tplasmid_name\\tlocus_idx\\tlocus_start\\tlocus_end\\n001\\tlp17\\t1\\t10\\t90\\n")
tables = out / "wp" / "tables"; tables.mkdir()
(tables / "a_all.tsv").write_text(
    "contig_id\\tcontig_len\\tplasmid_id\\tplasmid_name\\tstrain\\tgenospecies\\toverall_percent_identity\\t"
    "query_coverage_percent\\tquery_covered_length\\tquery_intervals\\n"
    "001\\t10000\\tgb|A|\\tlp17\\tB31\\tburgdorferi\\t98\\t90\\t9000\\t[(1, 9000)]\\n")
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


def hit(contig, name, species, identity, coverage, length):
    return {"contig_id": contig, "plasmid_name": name, "genospecies": species, "strain": "s",
            "plasmid_id": "gb|A|", "overall_percent_identity": identity,
            "query_coverage_percent": coverage, "query_covered_length": length}


def cells(rows, by="plasmid_name"):
    return {(cell["contig_id"], cell["label"]): cell for cell in identity_cells(pd.DataFrame(rows), by)}


def test_identity_cell_keeps_the_hit_the_caller_prefers():
    result = cells([
        hit("c1", "lp28", "burgdorferi", 100, 2, 300),
        hit("c1", "lp28", "garinii", 98, 95, 19000),
        hit("c1", "lp17", "afzelii", 99, 10, 2000),
    ])
    lp28, lp17 = result["c1", "lp28"], result["c1", "lp17"]
    assert lp28["overall_percent_identity"] == 98 and lp28["other"] == "garinii"
    assert lp28["top"] and not lp17["top"] and lp17["passes"]


def test_identity_cell_below_gates_is_kept_but_never_top():
    result = cells([hit("c1", "lp28", "burgdorferi", 100, 5, 300), hit("c2", None, None, None, None, None)])
    assert list(result) == [("c1", "lp28")]
    assert not result["c1", "lp28"]["passes"] and not result["c1", "lp28"]["top"]


def test_identity_ties_are_all_top():
    result = cells([hit("c1", "lp28", "burgdorferi", 99, 90, 9000),
                    hit("c1", "lp28", "garinii", 99, 90, 9000)], "genospecies")
    assert result["c1", "burgdorferi"]["top"] and result["c1", "garinii"]["top"]


def test_identity_endpoint_matches_summary_contig_ids(tmp_path, monkeypatch):
    with client_for(tmp_path, monkeypatch) as client:
        job = wait_for(client, upload(client).json()["jobs"][0]["id"])
        contig = client.get(f"/api/jobs/{job['id']}/results").json()["contigs"][0]
        identity = client.get(f"/api/jobs/{job['id']}/identity").json()
        assert contig["contig_id"] == "001"
        assert identity["plasmid"][0]["contig_id"] == "001" and identity["plasmid"][0]["top"] is True
        assert identity["genospecies"][0]["label"] == "burgdorferi"
        assert identity["gates"]["covered_bp"] == 1000
        # Parsed into real [start, end] pairs, not left as the TSV's stringified tuple list --
        # labelContigsTable draws tracks straight from this endpoint.
        assert identity["plasmid"][0]["query_intervals"] == [[1, 9000]]


def test_fusions_share_one_cell_and_keep_their_name():
    result = cells([
        hit("c1", "lp21-cp9", "burgdorferi", 99.9, 95, 18000),
        hit("c1", "cp32-1+3", "garinii", 90, 20, 4000),
        hit("c1", "lp21", "burgdorferi", 99, 90, 17000),
    ])
    assert set(result) == {("c1", "fusion"), ("c1", "lp21")}
    fusion = result["c1", "fusion"]
    assert fusion["fusion"] == "lp21-cp9" and fusion["parts"] == ["lp21", "cp9"] and fusion["top"]


def test_replicon_totals_split_multi_locus_calls():
    totals = replicon_totals([
        {"filtered_call": "lp56:::cp32-10", "contig_len": 50000},
        {"filtered_call": "cp32-12*", "contig_len": 30000},
        {"filtered_call": "lp56", "contig_len": 2000},
        {"filtered_call": "unclassified", "contig_len": 999},
        {"filtered_call": None, "contig_len": 999},
    ])
    assert totals == {"lp56": {"contigs": 2, "bp": 52000}, "cp32-10": {"contigs": 1, "bp": 50000},
                      "cp32-12": {"contigs": 1, "bp": 30000}}


def test_batch_overview_lists_finished_jobs_only(tmp_path, monkeypatch):
    with client_for(tmp_path, monkeypatch) as client:
        response = client.post("/api/batches", files=[
            ("files", ("good.fna", io.BytesIO(FASTA))), ("files", ("bad.fna", io.BytesIO(b"nope")))]).json()
        good = next(job for job in response["jobs"] if job["name"] == "good.fna")
        wait_for(client, good["id"])
        overview = client.get(f"/api/batches/{response['id']}/overview").json()
        assert [job["name"] for job in overview] == ["good.fna"]
        assert overview[0]["replicons"] == {"lp17": {"contigs": 1, "bp": 10}}
        assert overview[0]["composition"][0]["genospecies"] == "afzelii"
        assert client.get("/api/batches/nope/overview").status_code == 404


def test_contig_endpoint(tmp_path, monkeypatch):
    with client_for(tmp_path, monkeypatch) as client:
        job = wait_for(client, upload(client).json()["jobs"][0]["id"])
        url = f"/api/jobs/{job['id']}/contig"
        body = client.get(url, params={"contig_id": "001"}).json()
        assert body["rule"] == "No PF32 or WP hit." and body["contig"]["filtered_call"] == "lp17"
        assert body["plasmids"][0]["query_intervals"] == [[1, 9000]] and body["tiling"] == [["lp17", 9000]]
        assert body["loci"][0]["plasmid_name"] == "lp17" and body["loci"][0]["locus_end"] == 90
        assert client.get(url, params={"contig_id": "nope"}).status_code == 404


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


@pytest.mark.parametrize("cell, expected", [
    ("=cmd()", "'=cmd()"), ("+1+1", "'+1+1"), ("@sum", "'@sum"), ("-x", "'-x"),
    ("-1.5", "-1.5"), ("contig_1", "contig_1"), ("cp32-1+3", "cp32-1+3"), ("", ""),
])
def test_guard_cell(cell, expected):
    assert guard_cell(cell) == expected


def test_tsv_downloads_are_guarded(tmp_path, monkeypatch):
    script = FAKE_CLI.replace("001", "=cmd()")
    with client_for(tmp_path, monkeypatch, script=script) as client:
        job = wait_for(client, upload(client).json()["jobs"][0]["id"])
        text = client.get(f"/api/jobs/{job['id']}/download/summary_best_hits.tsv").text
        assert "\n'=cmd()\t" in text
        assert client.get(f"/api/jobs/{job['id']}/results").json()["contigs"][0]["contig_id"] == "=cmd()"


def test_jobs_are_stamped_with_provenance(tmp_path, monkeypatch):
    monkeypatch.setenv("WEB_VERSION", "v0.5.0")
    with client_for(tmp_path, monkeypatch) as client:
        version = client.get("/api/version").json()
        assert version["web_version"] == "v0.5.0"
        assert len(version["database_sha256"]) == 64 and version["blast_version"] != "unknown"
        job = wait_for(client, upload(client).json()["jobs"][0]["id"])
        assert job["provenance"] == version
        stamp = client.get(f"/api/jobs/{job['id']}/download/provenance.json").json()
        assert stamp["database_sha256"] == version["database_sha256"]


class FakeJWKS:
    def __init__(self, *args, **kwargs):
        pass

    def get_signing_key_from_jwt(self, token):
        return SimpleNamespace(key=ACCESS_KEY.public_key())


def access_client(tmp_path, monkeypatch):
    monkeypatch.setenv("CF_TEAM_DOMAIN", "team.cloudflareaccess.com")
    monkeypatch.setenv("CF_ACCESS_AUD", "aud")
    monkeypatch.setattr(auth, "PyJWKClient", FakeJWKS)
    return client_for(tmp_path, monkeypatch, access=True, public_host="portal.example.org")


def access_headers():
    now = int(time.time())
    token = jwt.encode({"aud": ["aud"], "iss": "https://team.cloudflareaccess.com", "iat": now,
                        "exp": now + 600, "email": "pi@example.org"}, ACCESS_KEY, algorithm="RS256")
    return {"Cf-Access-Jwt-Assertion": token, "host": "portal-abc.a.run.app"}


@pytest.mark.parametrize("path", ["/", "/static/app.js", "/api/batches", "/api/version",
                                  "/api/jobs/x", "/nonexistent", "/docs"])
def test_every_route_needs_access_token(tmp_path, monkeypatch, path):
    with access_client(tmp_path, monkeypatch) as client:
        assert client.get(path, headers={"host": "portal-abc.a.run.app"}).status_code == 403


def test_access_serves_run_app_host_and_public_origin(tmp_path, monkeypatch):
    with access_client(tmp_path, monkeypatch) as client:
        headers = access_headers()
        assert client.get("/api/version", headers=headers).status_code == 200
        files = [("files", ("a.fna", io.BytesIO(FASTA)))]
        response = client.post("/api/batches", files=files,
                               headers={**headers, "origin": "https://portal.example.org"})
        assert response.status_code == 201
        files = [("files", ("a.fna", io.BytesIO(FASTA)))]
        response = client.post("/api/batches", files=files,
                               headers={**headers, "origin": "https://evil.example"})
        assert response.status_code == 403


def test_access_refuses_to_start_half_configured(tmp_path, monkeypatch):
    monkeypatch.delenv("CF_TEAM_DOMAIN", raising=False)
    monkeypatch.delenv("CF_ACCESS_AUD", raising=False)
    monkeypatch.delenv("ALLOW_INSECURE_DEV", raising=False)
    with pytest.raises(RuntimeError):
        create_app(Settings(runs_dir=tmp_path / "runs", access=True))


def test_token_required_when_configured(tmp_path, monkeypatch):
    with client_for(tmp_path, monkeypatch, token="secret") as client:
        assert client.get("/api/batches").status_code == 401
        assert client.get("/?token=wrong").status_code == 401
        assert client.get("/?token=secret", follow_redirects=False).status_code == 303
        assert client.get("/api/batches").status_code == 200


def test_extract_contig_fasta_matches_by_header_suffix(tmp_path):
    path = tmp_path / "x_renamed.fasta"
    path.write_text(">file_id_001__lp17\nACGTACGTAC\nACGT\n>file_id_002__cp26\nTTTT\n")
    assert extract_contig_fasta(path, "001") == ">001\nACGTACGTACACGT\n"
    assert extract_contig_fasta(path, "002") == ">002\nTTTT\n"
    assert extract_contig_fasta(path, "nope") is None


def test_contig_fasta_endpoint(tmp_path, monkeypatch):
    script = FAKE_CLI.replace(
        '(out / "x_renamed.fasta").write_text(">c\\nA\\n")',
        '(out / "x_renamed.fasta").write_text(">file_001__lp17\\nACGTACGTAC\\n")')
    with client_for(tmp_path, monkeypatch, script=script) as client:
        job = wait_for(client, upload(client).json()["jobs"][0]["id"])
        response = client.get(f"/api/jobs/{job['id']}/contig/001/fasta")
        assert response.status_code == 200 and response.text == ">001\nACGTACGTAC\n"
        assert client.get(f"/api/jobs/{job['id']}/contig/nope/fasta").status_code == 404


def test_reference_db_endpoint(tmp_path, monkeypatch):
    with client_for(tmp_path, monkeypatch) as client:
        db = client.get("/api/reference-db").json()
        assert db["wp"] and db["pf32"]
        assert {"accession", "name", "strain", "genospecies", "length"} <= db["wp"][0].keys()
        assert {"accession", "name", "strain"} <= db["pf32"][0].keys()
        assert any(row["name"] != "unknown" for row in db["pf32"])


def test_duplicate_upload_is_found_not_done(tmp_path, monkeypatch):
    with client_for(tmp_path, monkeypatch) as client:
        first = wait_for(client, upload(client, name="a.fna").json()["jobs"][0]["id"])
        assert first["status"] == "done" and first["genome_hash"]
        second = upload(client, name="b.fna").json()["jobs"][0]
        assert second["status"] == "found"  # no polling needed: it's already finished
        assert second["genome_hash"] == first["genome_hash"]
        first_results = client.get(f"/api/jobs/{first['id']}/results").json()
        second_results = client.get(f"/api/jobs/{second['id']}/results").json()
        assert first_results["contigs"] == second_results["contigs"]


def test_bundle_downloads_are_built_and_served(tmp_path, monkeypatch):
    with client_for(tmp_path, monkeypatch) as client:
        job = wait_for(client, upload(client).json()["jobs"][0]["id"])
        for name in ["bundle.tar.gz", "bundle.zip"]:
            response = client.get(f"/api/jobs/{job['id']}/download/{name}")
            assert response.status_code == 200 and response.content


def test_visibility_toggle_controls_the_public_surface(tmp_path, monkeypatch):
    with client_for(tmp_path, monkeypatch) as client:
        job = wait_for(client, upload(client).json()["jobs"][0]["id"])
        genome_hash = job["genome_hash"]
        assert client.get("/api/public/genomes").json() == []
        assert client.get(f"/api/public/genomes/{genome_hash}").status_code == 404
        toggled = client.post(f"/api/jobs/{job['id']}/visibility", json={"visible": True}).json()
        assert toggled == {"genome_hash": genome_hash, "visible": True}
        assert [g["hash"] for g in client.get("/api/public/genomes").json()] == [genome_hash]
        public = client.get(f"/api/public/genomes/{genome_hash}")
        assert public.status_code == 200 and public.json()["contigs"][0]["final_call"] == "lp17"
        assert public.headers["x-robots-tag"].startswith("noindex")
        assert client.get(f"/api/public/genomes/{genome_hash}/identity").status_code == 200
        assert client.get(f"/api/public/genomes/{genome_hash}/contig",
                          params={"contig_id": "001"}).status_code == 200
        assert client.get(f"/api/public/genomes/{genome_hash}/download/summary_best_hits.tsv").status_code == 200
        assert "Disallow: /" in client.get("/robots.txt").text
        assert client.post(f"/api/jobs/{job['id']}/visibility", json={"visible": False}).status_code == 200
        assert client.get(f"/api/public/genomes/{genome_hash}").status_code == 404


def test_batch_overview_includes_found_jobs_not_just_done(tmp_path, monkeypatch):
    with client_for(tmp_path, monkeypatch) as client:
        wait_for(client, upload(client, name="a.fna").json()["jobs"][0]["id"])
        batch = upload(client, name="b.fna").json()
        dedup_job = batch["jobs"][0]
        assert dedup_job["status"] == "found"
        overview = client.get(f"/api/batches/{batch['id']}/overview").json()
        assert [row["name"] for row in overview] == ["b.fna"]
        assert overview[0]["replicons"] == {"lp17": {"contigs": 1, "bp": 10}}


def test_batch_visibility_controls_the_public_batch_overview(tmp_path, monkeypatch):
    with client_for(tmp_path, monkeypatch) as client:
        batch = upload(client).json()
        wait_for(client, batch["jobs"][0]["id"])
        assert client.get("/api/public/batches").json() == []
        assert client.get(f"/api/public/batches/{batch['id']}").status_code == 404
        toggled = client.post(f"/api/batches/{batch['id']}/visibility", json={"visible": True}).json()
        assert toggled == {"batch_id": batch["id"], "visible": True}
        assert [b["id"] for b in client.get("/api/public/batches").json()] == [batch["id"]]
        public_batch = client.get(f"/api/public/batches/{batch['id']}")
        assert public_batch.status_code == 200 and public_batch.json()["visible"] is True
        overview = client.get(f"/api/public/batches/{batch['id']}/overview").json()
        assert overview[0]["replicons"] == {"lp17": {"contigs": 1, "bp": 10}}
        assert client.post(f"/api/batches/{batch['id']}/visibility", json={"visible": False}).status_code == 200
        assert client.get(f"/api/public/batches/{batch['id']}").status_code == 404
        assert client.post("/api/batches/nope/visibility", json={"visible": True}).status_code == 404


def test_publishing_a_batch_cascades_to_its_genomes_but_unpublishing_does_not(tmp_path, monkeypatch):
    with client_for(tmp_path, monkeypatch) as client:
        batch = upload(client).json()
        job = wait_for(client, batch["jobs"][0]["id"])
        genome_hash = job["genome_hash"]
        assert client.get("/api/public/genomes").json() == []
        client.post(f"/api/batches/{batch['id']}/visibility", json={"visible": True})
        assert [g["hash"] for g in client.get("/api/public/genomes").json()] == [genome_hash]
        client.post(f"/api/batches/{batch['id']}/visibility", json={"visible": False})
        assert [g["hash"] for g in client.get("/api/public/genomes").json()] == [genome_hash]


def test_public_responses_do_not_expose_publisher_or_job_ids(tmp_path, monkeypatch):
    with access_client(tmp_path, monkeypatch) as client:
        headers = access_headers()
        files = [("files", ("a.fna", io.BytesIO(FASTA)))]
        batch = client.post("/api/batches", files=files, headers=headers).json()
        job = batch["jobs"][0]
        for _ in range(200):
            state = client.get(f"/api/jobs/{job['id']}", headers=headers).json()
            if state["status"] not in {"queued", "running"}:
                break
            time.sleep(0.05)
        hash_ = state["genome_hash"]
        client.post(f"/api/jobs/{job['id']}/visibility", json={"visible": True}, headers=headers)
        client.post(f"/api/batches/{batch['id']}/visibility", json={"visible": True}, headers=headers)
        paths = ["/api/public/genomes", f"/api/public/genomes/{hash_}", "/api/public/batches",
                 f"/api/public/genomes/{hash_}/contig?contig_id=001",
                 f"/api/public/batches/{batch['id']}", f"/api/public/batches/{batch['id']}/overview",
                 f"/api/public/genomes/{hash_}/download/provenance.json"]
        for path in paths:
            response = client.get(path, headers={"host": "portal-abc.a.run.app"})
            assert response.status_code == 200, path
            assert "pi@example.org" not in response.text and job["id"] not in response.text, path
            assert "visible_by" not in response.text and "job_id" not in response.text, path


def test_jobs_from_before_dedup_keep_their_results(tmp_path, monkeypatch):
    with client_for(tmp_path, monkeypatch) as client:
        batch = upload(client).json()
        job = wait_for(client, batch["jobs"][0]["id"])
        output = next((tmp_path / "runs" / "genomes").glob("*/output"))
    legacy = tmp_path / "runs" / job["batch_id"] / job["id"]
    shutil.move(output, legacy / "output")
    saved = json.loads((legacy / "job.json").read_text())
    del saved["genome_hash"]
    (legacy / "job.json").write_text(json.dumps(saved))
    (tmp_path / "runs" / "genomes.sqlite").unlink(missing_ok=True)
    with client_for(tmp_path, monkeypatch) as client:
        assert client.get(f"/api/jobs/{job['id']}/results").status_code == 200
        assert client.get(f"/api/batches/{job['batch_id']}/overview").status_code == 200
        response = client.post(f"/api/jobs/{job['id']}/visibility", json={"visible": True})
        assert response.status_code == 409


def test_public_batch_hides_unpublished_and_failed_jobs(tmp_path, monkeypatch):
    with client_for(tmp_path, monkeypatch) as client:
        batch = client.post("/api/batches", files=[
            ("files", ("shown.fna", io.BytesIO(FASTA))),
            ("files", ("hidden.fna", io.BytesIO(b">other\nACGTTGCA\n"))),
            ("files", ("rejected.fna", io.BytesIO(b"nope"))),
        ]).json()
        jobs = {job["name"]: job for job in batch["jobs"]}
        shown = wait_for(client, jobs["shown.fna"]["id"])
        hidden = wait_for(client, jobs["hidden.fna"]["id"])
        client.post(f"/api/batches/{batch['id']}/visibility", json={"visible": True})
        client.post(f"/api/jobs/{hidden['id']}/visibility", json={"visible": False})
        published = [f"/api/public/batches/{batch['id']}", f"/api/public/batches/{batch['id']}/overview",
                     "/api/public/batches"]
        for path in published:
            text = client.get(path).text
            assert "shown.fna" in text, path
            assert "hidden.fna" not in text and hidden["genome_hash"] not in text, path
            assert "rejected.fna" not in text, path
        overview = client.get(published[1]).json()
        assert [row["name"] for row in overview] == ["shown.fna"]
        client.post(f"/api/jobs/{shown['id']}/visibility", json={"visible": False})
        assert client.get(published[0]).status_code == 404
        assert client.get("/api/public/batches").json() == []


def test_public_routes_bypass_token_auth_but_the_toggle_does_not(tmp_path, monkeypatch):
    monkeypatch.setattr(JobStore, "_command", fake_command(FAKE_CLI))
    app = create_app(Settings(runs_dir=tmp_path / "runs", token="secret"))
    with TestClient(app, base_url="http://localhost") as admin, \
         TestClient(app, base_url="http://localhost") as anon:
        admin.get("/?token=secret", follow_redirects=False)
        job = wait_for(admin, upload(admin).json()["jobs"][0]["id"])
        admin.post(f"/api/jobs/{job['id']}/visibility", json={"visible": True})
        admin.post(f"/api/batches/{job['batch_id']}/visibility", json={"visible": True})
        genome_hash = job["genome_hash"]
        assert anon.get("/api/public/genomes").status_code == 200
        assert anon.get(f"/api/public/genomes/{genome_hash}").status_code == 200
        assert anon.get("/api/public/batches").status_code == 200
        assert anon.get(f"/api/public/batches/{job['batch_id']}/overview").status_code == 200
        assert anon.get("/api/batches").status_code == 401
        assert anon.post(f"/api/jobs/{job['id']}/visibility", json={"visible": True}).status_code == 401
        assert anon.post(f"/api/batches/{job['batch_id']}/visibility", json={"visible": True}).status_code == 401
