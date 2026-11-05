"""Job store and queue. Each assembly runs the plasmid_caller CLI; state lives on disk."""
import gzip
import json
import logging
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

from ..__about__ import __version__
from .provenance import collect

log = logging.getLogger("plasmid_caller.web")

FASTA_NAME = re.compile(r"\.(fa|fna|fasta)(\.gz)?$", re.IGNORECASE)
ID_PATTERN = re.compile(r"^[0-9a-f]{32}$")
NUCLEOTIDES = set("ACGTUNRYSWKMBDHV-")
UNFINISHED = {"queued", "running"}
REQUIRED_OUTPUTS = ("summary_best_hits.tsv", "summary_all_calls.tsv", "genospecies_composition.tsv")


@dataclass(frozen=True)
class Settings:
    runs_dir: Path
    jobs: int = 4
    threads: int = 4
    max_upload_mb: int = 500
    max_assembly_mb: int = 200
    max_queued: int = 200
    job_timeout: int = 3600
    retention_days: int = 0
    database: Path | None = None
    token: str | None = None
    access: bool = False
    public_host: str | None = None


class QueueFull(Exception):
    pass


def safe_filename(name):
    """Keep a readable name on disk. Never trust the client's path."""
    name = Path(str(name).replace("\\", "/")).name
    match = FASTA_NAME.search(name)
    suffix = match.group(0).lower() if match else ""
    stem = name[: match.start()] if match else name
    stem = re.sub(r"[^A-Za-z0-9._-]", "_", stem).strip("._")[:100]
    return f"{stem or 'assembly'}{suffix}"


def _open_binary(path):
    return gzip.open(path, "rb") if path.suffix.lower() == ".gz" else open(path, "rb")


def validate_fasta(path, max_bytes):
    """Return (contigs, bases) or raise ValueError with a message safe to show the uploader."""
    path = Path(path)
    if not FASTA_NAME.search(path.name):
        raise ValueError("Expected .fa, .fna or .fasta, optionally gzipped")
    try:
        # Measure before parsing: a small gzip can expand far past the limit.
        total = 0
        with _open_binary(path) as handle:
            while chunk := handle.read(1 << 20):
                total += len(chunk)
                if total > max_bytes:
                    raise ValueError(f"Assembly exceeds {max_bytes // 1_000_000} MB uncompressed")
        seen, bases = set(), 0
        opener = gzip.open if path.suffix.lower() == ".gz" else open
        with opener(path, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id in seen:
                    raise ValueError(f"Repeated contig ID: {record.id[:80]}")
                if not record.seq:
                    raise ValueError(f"Empty contig: {record.id[:80]}")
                if set(str(record.seq).upper()) - NUCLEOTIDES:
                    raise ValueError(f"Non-nucleotide characters in contig: {record.id[:80]}")
                seen.add(record.id)
                bases += len(record.seq)
    except (OSError, EOFError, UnicodeDecodeError) as exc:
        raise ValueError("File could not be read as FASTA") from exc
    if not seen:
        raise ValueError("No FASTA records found")
    return len(seen), bases


class JobStore:
    def __init__(self, settings):
        self.settings = settings
        self.runs_dir = Path(settings.runs_dir).resolve()
        self.lock = threading.Lock()
        self.jobs = {}
        self.batches = {}
        self.running = {}
        self.stopping = False
        self.executor = None
        self.provenance = None

    # ---------- lifecycle ----------
    def start(self):
        self.runs_dir.mkdir(parents=True, exist_ok=True)
        self.provenance = collect(self.settings.database)
        self.executor = ThreadPoolExecutor(self.settings.jobs, thread_name_prefix="job")
        self._load()
        self.prune()
        for job in sorted(self.jobs.values(), key=lambda j: j["created_at"]):
            if job["status"] in UNFINISHED:
                self._update(job["id"], status="queued", started_at=None)
                self.executor.submit(self._run, job["id"])

    def stop(self):
        self.stopping = True
        if self.executor:
            self.executor.shutdown(wait=False, cancel_futures=True)
        with self.lock:
            procs = list(self.running.values())
        for proc in procs:
            self._kill(proc)

    def _load(self):
        for batch_file in self.runs_dir.glob("*/batch.json"):
            batch = json.loads(batch_file.read_text())
            self.batches[batch["id"]] = batch
            for job_file in batch_file.parent.glob("*/job.json"):
                job = json.loads(job_file.read_text())
                self.jobs[job["id"]] = job

    def prune(self):
        days = self.settings.retention_days
        if not days:
            return
        cutoff = time.time() - days * 86_400
        with self.lock:
            expired = [
                batch["id"] for batch in self.batches.values()
                if batch["created_at"] < cutoff
                and all(j["status"] not in UNFINISHED for j in self._batch_jobs(batch["id"]))
            ]
            for batch_id in expired:
                for job in self._batch_jobs(batch_id):
                    del self.jobs[job["id"]]
                del self.batches[batch_id]
        for batch_id in expired:
            shutil.rmtree(self.runs_dir / batch_id, ignore_errors=True)

    # ---------- reads ----------
    def _batch_jobs(self, batch_id):
        jobs = [job for job in self.jobs.values() if job["batch_id"] == batch_id]
        return sorted(jobs, key=lambda job: job["name"])

    def list_batches(self):
        with self.lock:
            batches = [self._batch_view(batch) for batch in self.batches.values()]
        return sorted(batches, key=lambda batch: batch["created_at"], reverse=True)

    def _batch_view(self, batch):
        return {**batch, "jobs": [dict(job) for job in self._batch_jobs(batch["id"])]}

    def get_batch(self, batch_id):
        with self.lock:
            batch = self.batches.get(batch_id)
            return self._batch_view(batch) if batch else None

    def get_job(self, job_id):
        with self.lock:
            job = self.jobs.get(job_id)
            return dict(job) if job else None

    def job_dir(self, job):
        return self.runs_dir / job["batch_id"] / job["id"]

    def output_dir(self, job):
        return self.job_dir(job) / "output"

    # ---------- writes ----------
    def create_batch(self, uploads):
        """uploads: iterable of (client filename, readable binary file)."""
        uploads = list(uploads)
        with self.lock:
            queued = sum(job["status"] in UNFINISHED for job in self.jobs.values())
        if queued + len(uploads) > self.settings.max_queued:
            raise QueueFull(f"Queue is full ({self.settings.max_queued} assemblies)")

        batch = {"id": uuid.uuid4().hex, "created_at": time.time(), "version": __version__}
        batch_dir = self.runs_dir / batch["id"]
        batch_dir.mkdir()
        _write_json(batch_dir / "batch.json", batch)
        jobs = [self._create_job(batch["id"], name, handle) for name, handle in uploads]
        with self.lock:
            self.batches[batch["id"]] = batch
            self.jobs.update({job["id"]: job for job in jobs})
            view = self._batch_view(batch)
        for job in jobs:
            if job["status"] == "queued":
                self.executor.submit(self._run, job["id"])
        return view

    def _create_job(self, batch_id, client_name, handle):
        filename = safe_filename(client_name)
        job = {
            "id": uuid.uuid4().hex, "batch_id": batch_id, "name": filename,
            "status": "queued", "error": None, "contigs": None, "bases": None,
            "created_at": time.time(), "started_at": None, "finished_at": None,
            "provenance": self.provenance,
        }
        input_dir = self.job_dir(job) / "input"
        input_dir.mkdir(parents=True)
        path = input_dir / filename
        with open(path, "wb") as out:
            shutil.copyfileobj(handle, out, 1 << 20)
        try:
            job["contigs"], job["bases"] = validate_fasta(
                path, self.settings.max_assembly_mb * 1_000_000
            )
        except ValueError as exc:
            job.update(status="rejected", error=str(exc), finished_at=time.time())
            path.unlink()
        _write_json(self.job_dir(job) / "job.json", job)
        return job

    def _update(self, job_id, **fields):
        with self.lock:
            job = self.jobs[job_id]
            job.update(fields)
            snapshot = dict(job)
        _write_json(self.job_dir(snapshot) / "job.json", snapshot)

    # ---------- execution ----------
    def _command(self, input_path, output):
        cmd = [
            sys.executable, "-m", "plasmid_caller.main",
            "-i", str(input_path), "-o", str(output),
            "-t", str(self.settings.threads), "-q",
        ]
        if self.settings.database:
            cmd += ["-db", str(self.settings.database)]
        return cmd

    def _kill(self, proc):
        # BLAST is a grandchild, so signal the whole process group.
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass

    def _run(self, job_id):
        job = self.get_job(job_id)
        if job is None or self.stopping:
            return
        self._update(job_id, status="running", started_at=time.time())
        job_dir, output = self.job_dir(job), self.output_dir(job)
        shutil.rmtree(output, ignore_errors=True)
        input_path = job_dir / "input" / job["name"]
        error = None
        try:
            with open(job_dir / "log.txt", "w") as log_file:
                proc = subprocess.Popen(
                    self._command(input_path, output), stdout=log_file,
                    stderr=subprocess.STDOUT, start_new_session=True,
                )
                with self.lock:
                    self.running[job_id] = proc
                try:
                    code = proc.wait(self.settings.job_timeout)
                except subprocess.TimeoutExpired:
                    self._kill(proc)
                    proc.wait()
                    code, error = None, "Timed out"
                finally:
                    with self.lock:
                        self.running.pop(job_id, None)
            if self.stopping:
                return  # Requeued on the next start.
            if code:
                error = "The pipeline failed. Details are in this job's log.txt on the server."
            elif code == 0 and not all((output / name).is_file() for name in REQUIRED_OUTPUTS):
                error = "The pipeline finished without writing its results."
        except Exception:
            log.exception("Job %s failed", job_id)
            error = "The server could not run this job."
        if error is None:
            self._clean(output)
        self._update(job_id, status="failed" if error else "done", error=error,
                     finished_at=time.time())
        self.prune()

    def _clean(self, output):
        """The viewer reads tables, so drop the large intermediates (WP XML reaches ~500 MB)."""
        for path in [*output.glob("*/xml_files"), output / "sanitized"]:
            shutil.rmtree(path, ignore_errors=True)


def _write_json(path, data):
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data))
    tmp.replace(path)
