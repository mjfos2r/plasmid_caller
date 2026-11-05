"""FastAPI app: uploads, job status, and read-only access to CLI output."""
import secrets
from contextlib import asynccontextmanager
from pathlib import Path
from urllib.parse import urlsplit

import pandas as pd
from fastapi import FastAPI, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, JSONResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles

from .jobs import ID_PATTERN, JobStore, QueueFull

STATIC = Path(__file__).parent / "static"
TOKEN_COOKIE = "plasmid_caller_token"
SECURITY_HEADERS = {
    "Content-Security-Policy": "default-src 'self'; frame-ancestors 'none'; base-uri 'none'",
    "X-Content-Type-Options": "nosniff",
    "Referrer-Policy": "no-referrer",
}
SUMMARY_COLUMNS = [
    "contig_id", "contig_len", "final_call", "filtered_call", "genospecies",
    "genospecies_candidates", "plasmid_name_pf32", "plasmid_name_wp", "strain_wp",
    "overall_percent_identity_wp", "query_coverage_percent_wp",
    "overall_percent_identity_pf32", "query_coverage_percent_pf32", "multiple_loci_pf32",
]
DOWNLOADS = {
    "summary_best_hits.tsv": "text/tab-separated-values",
    "summary_all_calls.tsv": "text/tab-separated-values",
    "genospecies_composition.tsv": "text/tab-separated-values",
    "summary_best_hits.json": "application/json",
}


def read_records(path, columns=None):
    frame = pd.read_csv(path, sep="\t", keep_default_na=True, na_values=[""])
    if columns:
        frame = frame.reindex(columns=columns)
    frame = frame.astype(object).where(frame.notna(), None)
    return frame.to_dict(orient="records")


def create_app(settings):
    store = JobStore(settings)

    @asynccontextmanager
    async def lifespan(app):
        store.start()
        yield
        store.stop()

    app = FastAPI(title="plasmid_caller", lifespan=lifespan, docs_url=None, redoc_url=None,
                  openapi_url=None)
    app.state.store = store
    max_body = settings.max_upload_mb * 1_000_000

    def authorised(request):
        if not settings.token:
            return True
        supplied = request.cookies.get(TOKEN_COOKIE) or ""
        return secrets.compare_digest(supplied.encode(), settings.token.encode())

    @app.middleware("http")
    async def guard(request: Request, call_next):
        host = request.headers.get("host", "")
        if not settings.token and urlsplit(f"//{host}").hostname not in {"localhost", "127.0.0.1", "::1"}:
            # DNS rebinding: an unauthenticated server only answers to loopback names.
            return JSONResponse({"detail": "Unknown host"}, 403)
        if request.url.path == "/" and settings.token and "token" in request.query_params:
            supplied = request.query_params["token"]
            if not secrets.compare_digest(supplied.encode(), settings.token.encode()):
                return JSONResponse({"detail": "Invalid token"}, 401)
            response = RedirectResponse("/", 303)
            response.set_cookie(TOKEN_COOKIE, supplied, httponly=True, samesite="strict",
                                secure=request.url.scheme == "https")
            return response
        if not authorised(request):
            return JSONResponse({"detail": "Open the link containing the access token"}, 401)
        if request.method == "POST":
            origin = request.headers.get("origin")
            if origin and urlsplit(origin).netloc != host:
                return JSONResponse({"detail": "Cross-origin request refused"}, 403)
            length = request.headers.get("content-length")
            if length is None:
                return JSONResponse({"detail": "Content-Length required"}, 411)
            if int(length) > max_body:
                return JSONResponse(
                    {"detail": f"Upload exceeds {settings.max_upload_mb} MB"}, 413)
        response = await call_next(request)
        response.headers.update(SECURITY_HEADERS)
        return response

    def job_or_404(job_id):
        job = store.get_job(job_id) if ID_PATTERN.match(job_id) else None
        if job is None:
            raise HTTPException(404, "Job not found")
        return job

    def finished_job(job_id):
        job = job_or_404(job_id)
        if job["status"] != "done":
            raise HTTPException(409, f"Job is {job['status']}")
        return job

    @app.get("/")
    def index():
        return FileResponse(STATIC / "index.html")

    @app.get("/api/batches")
    def list_batches():
        return store.list_batches()

    # Plain def: FastAPI runs it in a worker thread, so FASTA validation won't block the server.
    @app.post("/api/batches", status_code=201)
    def create_batch(files: list[UploadFile]):
        try:
            return store.create_batch((f.filename, f.file) for f in files)
        except QueueFull as exc:
            raise HTTPException(429, str(exc)) from exc

    @app.get("/api/batches/{batch_id}")
    def get_batch(batch_id: str):
        batch = store.get_batch(batch_id) if ID_PATTERN.match(batch_id) else None
        if batch is None:
            raise HTTPException(404, "Batch not found")
        return batch

    @app.get("/api/jobs/{job_id}")
    def get_job(job_id: str):
        return job_or_404(job_id)

    @app.get("/api/jobs/{job_id}/results")
    def get_results(job_id: str):
        job = finished_job(job_id)
        output = store.output_dir(job)
        return {
            "job": job,
            "contigs": read_records(output / "summary_best_hits.tsv", SUMMARY_COLUMNS),
            "composition": read_records(output / "genospecies_composition.tsv"),
            "downloads": [*DOWNLOADS, "renamed.fasta"],
        }

    @app.get("/api/jobs/{job_id}/download/{name}")
    def download(job_id: str, name: str):
        job = finished_job(job_id)
        output = store.output_dir(job)
        stem = Path(job["name"]).name.split(".")[0]
        if name == "renamed.fasta":
            matches = list(output.glob("*_renamed.fasta"))
            if not matches:
                raise HTTPException(404, "File not found")
            return FileResponse(matches[0], media_type="text/plain",
                                filename=f"{stem}_renamed.fasta")
        if name not in DOWNLOADS:
            raise HTTPException(404, "File not found")
        return FileResponse(output / name, media_type=DOWNLOADS[name],
                            filename=f"{stem}_{name}")

    app.mount("/static", StaticFiles(directory=STATIC), name="static")
    return app
