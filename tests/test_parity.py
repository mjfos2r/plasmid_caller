"""Check a running portal against the local CLI, row for row.

Set PORTAL_URL and PARITY_FASTAS (space-separated paths). Add PORTAL_HEADERS as JSON for
a deployed portal. Skipped otherwise, because each assembly takes a real BLAST run.
"""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import pytest

httpx = pytest.importorskip("httpx")

URL = os.environ.get("PORTAL_URL", "").rstrip("/")
FASTAS = os.environ.get("PARITY_FASTAS", "").split()
HEADERS = json.loads(os.environ.get("PORTAL_HEADERS", "{}"))
pytestmark = pytest.mark.skipif(not URL or not FASTAS, reason="Set PORTAL_URL and PARITY_FASTAS")


def run_on_portal(path):
    with httpx.Client(base_url=URL, headers=HEADERS, timeout=120) as client:
        with open(path, "rb") as handle:
            batch = client.post("/api/batches", files=[("files", (Path(path).name, handle))])
        batch.raise_for_status()
        job_id = batch.json()["jobs"][0]["id"]
        for _ in range(1800):
            job = client.get(f"/api/jobs/{job_id}").json()
            if job["status"] not in {"queued", "running"}:
                break
            time.sleep(1)
        assert job["status"] == "done", job
        results = client.get(f"/api/jobs/{job_id}/results").json()
        download = client.get(f"/api/jobs/{job_id}/download/summary_best_hits.tsv").content
    print(f"{Path(path).name}: {job['contigs']} contigs, {job['bases']} bp, "
          f"{job['finished_at'] - job['started_at']:.0f} s")
    return results, download


@pytest.mark.parametrize("path", FASTAS)
def test_portal_matches_cli(path, tmp_path):
    output = tmp_path / "cli"
    subprocess.run([sys.executable, "-m", "plasmid_caller.main", "-i", path, "-o", str(output),
                    "-t", "8", "-q"], check=True, capture_output=True)
    expected = (output / "summary_best_hits.tsv").read_bytes()
    results, download = run_on_portal(path)
    assert download == expected
    table = pd.read_csv(output / "summary_best_hits.tsv", sep="\t", dtype={"contig_id": str})
    calls = {row["contig_id"]: row["filtered_call"] for row in results["contigs"]}
    assert calls == dict(zip(table["contig_id"], table["filtered_call"].where(table["filtered_call"].notna(), None)))
