"""Job store and queue. Each assembly runs the plasmid_caller CLI; state lives on disk."""
import gzip
import hashlib
import io
import json
import logging
import os
import re
import shutil
import signal
import subprocess
import sys
import tarfile
import threading
import time
import uuid
import zipfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

from ..__about__ import __version__
from . import reference_db
from .genome_index import GenomeIndex
from .provenance import collect
from .storage import build_storage

log = logging.getLogger("plasmid_caller.web")

FASTA_NAME = re.compile(r"\.(fa|fna|fasta)(\.gz)?$", re.IGNORECASE)
ID_PATTERN = re.compile(r"^[0-9a-f]{32}$")
HASH_PATTERN = re.compile(r"^[0-9a-f]{64}$")
NUCLEOTIDES = set("ACGTUNRYSWKMBDHV-")
UNFINISHED = {"queued", "running"}
REQUIRED_OUTPUTS = ("summary_best_hits.tsv", "summary_all_calls.tsv", "genospecies_composition.tsv")
BUNDLE_SOURCE_NAMES = ("summary_best_hits.tsv", "summary_all_calls.tsv",
                       "genospecies_composition.tsv", "summary_best_hits.json")


@dataclass(frozen=True)
class Settings:
    runs_dir: Path
    jobs: int = 4
    threads: int = 4
    max_upload_mb: int = 500
    max_assembly_mb: int = 250
    max_queued: int = 200
    job_timeout: int = 3600
    retention_days: int = 0
    database: Path | None = None
    bucket: str | None = None
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


def contig_digest(seq_upper):
    """sha256 of one contig's uppercased sequence. Pure, no file I/O."""
    return hashlib.sha256(seq_upper.encode()).hexdigest()


def genome_hash(contig_digests):
    """Combine per-contig digests into one hash, independent of header names and contig order.

    Sorting the digests is what makes order not matter; hashing the uppercased sequence
    (not raw bytes) is what makes whitespace/line-wrapping/gzip not matter, since SeqIO
    already normalizes those away before we see `record.seq`.

    NOT folded in: reverse-complement orientation. A whole-genome reverse-complement of
    the same organism hashes differently from the forward-strand assembly today -- a
    deliberate, conservative default.
    """
    return hashlib.sha256("".join(sorted(contig_digests)).encode()).hexdigest()


def validate_fasta(path, max_bytes):
    """Return (contigs, bases, genome_hash) or raise ValueError with a message safe to show the uploader."""
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
        seen, bases, digests = set(), 0, []
        opener = gzip.open if path.suffix.lower() == ".gz" else open
        with opener(path, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id in seen:
                    raise ValueError(f"Repeated contig ID: {record.id[:80]}")
                if not record.seq:
                    raise ValueError(f"Empty contig: {record.id[:80]}")
                seq_upper = str(record.seq).upper()
                if set(seq_upper) - NUCLEOTIDES:
                    raise ValueError(f"Non-nucleotide characters in contig: {record.id[:80]}")
                seen.add(record.id)
                bases += len(record.seq)
                digests.append(contig_digest(seq_upper))
    except (OSError, EOFError, UnicodeDecodeError) as exc:
        raise ValueError("File could not be read as FASTA") from exc
    if not seen:
        raise ValueError("No FASTA records found")
    return len(seen), bases, genome_hash(digests)


def _tar_add_bytes(tar, name, data):
    info = tarfile.TarInfo(name)
    info.size = len(data)
    tar.addfile(info, io.BytesIO(data))


def build_bundles(output, job):
    """Build bundle.tar.gz/.zip + checksums.txt once per genome, from the files the viewer already serves."""
    files = {name: output / name for name in BUNDLE_SOURCE_NAMES if (output / name).is_file()}
    renamed = list(output.glob("*_renamed.fasta"))
    if renamed:
        files["renamed.fasta"] = renamed[0]
    provenance_bytes = json.dumps({"job": job["id"], "finished_at": job["finished_at"],
                                   **(job.get("provenance") or {})}, indent=2).encode()

    checksums = [f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {name}"
                 for name, path in sorted(files.items())]
    checksums.append(f"{hashlib.sha256(provenance_bytes).hexdigest()}  provenance.json")
    checksums.append(f"{job['genome_hash']}  (genome content hash)")
    checksums_text = "\n".join(checksums) + "\n"

    with tarfile.open(output / "bundle.tar.gz", "w:gz") as tar:
        for name, path in files.items():
            tar.add(path, arcname=name)
        _tar_add_bytes(tar, "provenance.json", provenance_bytes)
        _tar_add_bytes(tar, "checksums.txt", checksums_text.encode())

    with zipfile.ZipFile(output / "bundle.zip", "w", zipfile.ZIP_DEFLATED) as zf:
        for name, path in files.items():
            zf.write(path, arcname=name)
        zf.writestr("provenance.json", provenance_bytes)
        zf.writestr("checksums.txt", checksums_text)


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
        self.reference_db = None
        self.storage = None
        self.index = None

    # ---------- lifecycle ----------
    def start(self):
        self.runs_dir.mkdir(parents=True, exist_ok=True)
        self.provenance = collect(self.settings.database)
        self.reference_db = reference_db.load(self.settings.database)
        self.storage = build_storage(self.settings.bucket)
        self.index = GenomeIndex(self.runs_dir / "index.sqlite3", self.storage)
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
            batch.setdefault("visible", False)  # pre-visibility batch.json, see plan §8
            self.batches[batch["id"]] = batch
            for job_file in batch_file.parent.glob("*/job.json"):
                job = json.loads(job_file.read_text())
                job.setdefault("genome_hash", None)  # pre-dedup job.json, see plan §8
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

    def _view(self, job):
        """A job dict plus whether its genome is currently public, for display only."""
        view = dict(job)
        row = self.index.lookup(job["genome_hash"]) if job.get("genome_hash") else None
        view["visible"] = row["visible"] if row else None
        return view

    def _batch_view(self, batch):
        return {**batch, "jobs": [self._view(job) for job in self._batch_jobs(batch["id"])]}

    def get_batch(self, batch_id):
        with self.lock:
            batch = self.batches.get(batch_id)
            return self._batch_view(batch) if batch else None

    def list_visible_batches(self):
        with self.lock:
            batches = [self._batch_view(batch) for batch in self.batches.values() if batch.get("visible")]
        return sorted(batches, key=lambda batch: batch["created_at"], reverse=True)

    def set_batch_visible(self, batch_id, visible, by=None):
        """Publishing a batch also publishes every genome already called in it -- otherwise its
        own table is just a page of names with no links. Un-publishing does NOT reverse this:
        a genome made visible this way may be shared (by content hash) with another batch that's
        still public, so hiding it back is left as a deliberate per-genome action."""
        with self.lock:
            batch = self.batches.get(batch_id)
            if batch is None:
                return False
            batch["visible"] = visible
            snapshot = dict(batch)
            jobs = self._batch_jobs(batch_id) if visible else []
        _write_json(self.runs_dir / batch_id / "batch.json", snapshot)
        for job in jobs:
            if job["status"] in {"done", "found"} and job.get("genome_hash"):
                self.index.set_visible(job["genome_hash"], True, by=by)
        return True

    def get_job(self, job_id):
        with self.lock:
            job = self.jobs.get(job_id)
            return self._view(job) if job else None

    def finished_jobs(self):
        """Every done/found job -- for cross-referencing by name, e.g. reference-strain lookup."""
        with self.lock:
            return [self._view(job) for job in self.jobs.values() if job["status"] in {"done", "found"}]

    def job_dir(self, job):
        return self.runs_dir / job["batch_id"] / job["id"]

    def _staging_dir(self, job):
        """Where _run writes a fresh pipeline run, before it's promoted to the canonical slot."""
        return self.job_dir(job) / "output"

    def output_dir_for_hash(self, genome_hash):
        """The durable, hash-addressed result directory, restored from the bucket on a local miss."""
        if genome_hash is None:
            return self.runs_dir / "genomes" / "_unknown" / "output"
        path = self.runs_dir / "genomes" / genome_hash / "output"
        if not path.is_dir():
            self.storage.download_genome(genome_hash, path)
        return path

    def output_dir(self, job):
        # Jobs from before dedup have no hash, and their results stayed in their own folder.
        if job.get("genome_hash") is None:
            return self.job_dir(job) / "output"
        return self.output_dir_for_hash(job["genome_hash"])

    # ---------- writes ----------
    def create_batch(self, uploads):
        """uploads: iterable of (client filename, readable binary file)."""
        uploads = list(uploads)
        with self.lock:
            queued = sum(job["status"] in UNFINISHED for job in self.jobs.values())
        if queued + len(uploads) > self.settings.max_queued:
            raise QueueFull(f"Queue is full ({self.settings.max_queued} assemblies)")

        batch = {"id": uuid.uuid4().hex, "created_at": time.time(), "version": __version__, "visible": False}
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
            "provenance": self.provenance, "genome_hash": None,
        }
        input_dir = self.job_dir(job) / "input"
        input_dir.mkdir(parents=True)
        path = input_dir / filename
        with open(path, "wb") as out:
            shutil.copyfileobj(handle, out, 1 << 20)
        try:
            job["contigs"], job["bases"], job["genome_hash"] = validate_fasta(
                path, self.settings.max_assembly_mb * 1_000_000
            )
        except ValueError as exc:
            job.update(status="rejected", error=str(exc), finished_at=time.time())
            path.unlink()
        else:
            existing = self.index.lookup(job["genome_hash"])
            if existing:
                # Already called: report "found", not "done" -- nothing ran for this upload.
                job.update(status="found", started_at=time.time(), finished_at=time.time(),
                           provenance=existing["provenance"])
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
        job_dir, staging = self.job_dir(job), self._staging_dir(job)
        shutil.rmtree(staging, ignore_errors=True)
        input_path = job_dir / "input" / job["name"]
        error = None
        try:
            with open(job_dir / "log.txt", "w") as log_file:
                proc = subprocess.Popen(
                    self._command(input_path, staging), stdout=log_file,
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
            elif code == 0 and not all((staging / name).is_file() for name in REQUIRED_OUTPUTS):
                error = "The pipeline finished without writing its results."
        except Exception:
            log.exception("Job %s failed", job_id)
            error = "The server could not run this job."
        finished_at = time.time()
        if error is None:
            self._clean(staging)
            build_bundles(staging, {**job, "finished_at": finished_at})
            self._promote(job, staging)
        self._update(job_id, status="failed" if error else "done", error=error, finished_at=finished_at)
        self.prune()

    def _promote(self, job, staging):
        """Claim the canonical slot for this genome hash, or discard a redundant duplicate run."""
        genome_hash = job["genome_hash"]
        canonical = self.output_dir_for_hash(genome_hash)
        claimed = self.index.try_insert(genome_hash, job_id=job["id"], name=job["name"],
                                         contigs=job["contigs"], bases=job["bases"],
                                         provenance=job["provenance"])
        if claimed:
            canonical.parent.mkdir(parents=True, exist_ok=True)
            shutil.rmtree(canonical, ignore_errors=True)
            shutil.move(str(staging), str(canonical))
            self.storage.upload_genome(genome_hash, canonical)
        else:
            shutil.rmtree(staging, ignore_errors=True)  # someone else's run won the race

    def _clean(self, output):
        """The viewer reads tables, so drop the large intermediates (WP XML reaches ~500 MB)."""
        for path in [*output.glob("*/xml_files"), output / "sanitized"]:
            shutil.rmtree(path, ignore_errors=True)


def _write_json(path, data):
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data))
    tmp.replace(path)
