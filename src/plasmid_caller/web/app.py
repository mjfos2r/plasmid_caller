"""FastAPI app: uploads, job status, and read-only access to CLI output."""
import secrets
from contextlib import asynccontextmanager
from pathlib import Path
from urllib.parse import urlsplit

import pandas as pd
from fastapi import FastAPI, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, JSONResponse, PlainTextResponse, RedirectResponse, Response
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from ..replicons import fusion_candidates, is_fusion, replicon_parts, tile_contig
from ..scoring import WP_MIN_COV_BP, WP_MIN_COV_PCT, _parse_intervals, decide_call, wp_hit_passes
from . import reference_db
from .auth import install_access_auth
from .jobs import HASH_PATTERN, ID_PATTERN, JobStore, QueueFull

STATIC = Path(__file__).parent / "static"
TOKEN_COOKIE = "plasmid_caller_token"
# Unauthenticated, read-only, kept out of search indexes (see robots() and the noindex header below).
PUBLIC_PATHS = frozenset({"/public", "/robots.txt", "/static/public.js", "/static/common.js", "/static/style.css"})
PUBLIC_PREFIX = "/api/public/"
SECURITY_HEADERS = {
    "Content-Security-Policy": "default-src 'self'; frame-ancestors 'none'; base-uri 'none'",
    "X-Content-Type-Options": "nosniff",
    "Referrer-Policy": "no-referrer",
}
NOINDEX_HEADER = {"X-Robots-Tag": "noindex, nofollow"}
SUMMARY_COLUMNS = [
    "contig_id", "contig_len", "final_call", "filtered_call", "genospecies",
    "genospecies_candidates", "plasmid_name_pf32", "plasmid_name_wp", "strain_wp",
    "overall_percent_identity_wp", "query_coverage_percent_wp",
    "overall_percent_identity_pf32", "query_coverage_percent_pf32", "multiple_loci_pf32",
]
IDENTITY_COLUMNS = [
    "contig_id", "label", "overall_percent_identity", "query_coverage_percent",
    "query_covered_length", "strain", "plasmid_id", "passes", "top", "score", "fusion", "parts",
    "ref_length", "ref_covered_length", "query_intervals",
]
LOCUS_COLUMNS = ["plasmid_name", "locus_idx", "locus_start", "locus_end"]
PF32_COLUMNS = [
    "plasmid_name", "strain", "plasmid_id", "overall_percent_identity", "query_covered_length",
    "query_intervals",
]
FORMULA_PREFIXES = ("=", "+", "-", "@")
DOWNLOADS = {
    "summary_best_hits.tsv": "text/tab-separated-values",
    "summary_all_calls.tsv": "text/tab-separated-values",
    "genospecies_composition.tsv": "text/tab-separated-values",
    "summary_best_hits.json": "application/json",
    "bundle.tar.gz": "application/gzip",
    "bundle.zip": "application/zip",
}
# The bundles get their own buttons at the top of the results view, not a line in this list.
RESULT_DOWNLOADS = [name for name in [*DOWNLOADS, "renamed.fasta", "provenance.json"]
                    if name not in {"bundle.tar.gz", "bundle.zip"}]


def guard_cell(cell):
    """Stop a contig name like =cmd() running as a formula when a TSV opens in a spreadsheet."""
    if cell[:1] not in FORMULA_PREFIXES:
        return cell
    try:
        float(cell)
        return cell
    except ValueError:
        return f"'{cell}"


def guarded_tsv(path):
    lines = Path(path).read_text().splitlines()
    return "".join("\t".join(map(guard_cell, line.split("\t"))) + "\n" for line in lines)


def read_table(path):
    # Contig IDs stay text: "001" must not become 1.
    return pd.read_csv(path, sep="\t", keep_default_na=True, na_values=[""],
                       dtype={"contig_id": str})


def to_records(frame):
    return frame.astype(object).where(frame.notna(), None).to_dict(orient="records")


def with_intervals(records):
    return [{**row, "query_intervals": _parse_intervals(row["query_intervals"])} for row in records]


def extract_contig_fasta(renamed_fasta, contig_id):
    """Pull one contig's own record out of the renamed multi-FASTA, by its original contig_id.

    Headers there look like "{file_id}_{contig_id}__{call}" (rename_fasta_headers in main.py),
    so a record belongs to contig_id when its header, with the call suffix stripped, ends with
    it. Streamed line by line rather than loaded whole, since an assembly can run to ~200 MB.
    """
    suffix = f"_{contig_id}"
    sequence, found = [], False
    with open(renamed_fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                if found:
                    break
                header = line[1:].strip().rsplit("__", 1)[0]
                found = header == contig_id or header.endswith(suffix)
                continue
            if found:
                sequence.append(line.strip())
    if not found:
        return None
    bases = "".join(sequence)
    wrapped = "\n".join(bases[i:i + 70] for i in range(0, len(bases), 70))
    return f">{contig_id}\n{wrapped}\n"


def read_records(path, columns=None):
    frame = read_table(path)
    return to_records(frame.reindex(columns=columns) if columns else frame)


def replicon_totals(contigs):
    """Sum called bases per replicon. A multi-locus contig counts towards each of its names."""
    totals = {}
    for row in contigs:
        call = row.get("filtered_call")
        if not call or call == "unclassified":
            continue
        for name in {part.rstrip("*") for part in call.split(":::")}:
            entry = totals.setdefault(name, {"contigs": 0, "bp": 0})
            entry["contigs"] += 1
            entry["bp"] += int(row["contig_len"])
    return totals


def identity_cells(hits, by, fold_fusions=True):
    """One row per contig and `by` value: the hit the caller would prefer among them.

    `top` marks the highest-scoring hit on the contig, ties included, as in
    taxonomy.annotate_best_hits. Hits below the WP gates can never be top.
    With `fold_fusions`, fused plasmids share a single "fusion" label.
    """
    hits = hits.dropna(subset=[by]).copy()
    label = hits[by].astype(object)
    hits["fusion"], hits["parts"] = None, [[] for _ in range(len(hits))]
    if by == "plasmid_name":
        fused = hits["plasmid_name"].map(is_fusion)
        hits.loc[fused, "fusion"] = hits.loc[fused, "plasmid_name"]
        hits["parts"] = [replicon_parts(name) if flag else []
                         for name, flag in zip(hits["plasmid_name"], fused)]
        if fold_fusions:
            label = label.where(~fused, "fusion")
    hits["label"] = label
    hits["passes"] = wp_hit_passes(hits)
    hits["score"] = hits["query_coverage_percent"] * hits["overall_percent_identity"]
    best = (hits.sort_values(["passes", "score"], ascending=False, kind="stable")
            .drop_duplicates(["contig_id", "label"]))
    peak = best.loc[best["passes"]].groupby("contig_id")["score"].max()
    peak = best["contig_id"].map(peak)
    best["top"] = best["passes"] & ((best["score"] - peak).abs() <= 1e-9 * peak)
    best["other"] = best["genospecies" if by == "plasmid_name" else "plasmid_name"]
    return to_records(best.reindex(columns=IDENTITY_COLUMNS + ["other"]))


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
        loopback = urlsplit(f"//{host}").hostname in {"localhost", "127.0.0.1", "::1"}
        if not settings.token and not settings.access and not loopback:
            # DNS rebinding: an unauthenticated server only answers to loopback names.
            return JSONResponse({"detail": "Unknown host"}, 403)
        path = request.url.path
        is_public = path in PUBLIC_PATHS or path.startswith(PUBLIC_PREFIX)
        if path == "/" and settings.token and "token" in request.query_params:
            supplied = request.query_params["token"]
            if not secrets.compare_digest(supplied.encode(), settings.token.encode()):
                return JSONResponse({"detail": "Invalid token"}, 401)
            response = RedirectResponse("/", 303)
            response.set_cookie(TOKEN_COOKIE, supplied, httponly=True, samesite="strict",
                                secure=request.url.scheme == "https")
            return response
        if not is_public and not authorised(request):
            return JSONResponse({"detail": "Open the link containing the access token"}, 401)
        if request.method == "POST":
            origin = request.headers.get("origin")
            # Behind the Worker, Host is the run.app address but Origin is the public name.
            if origin and urlsplit(origin).netloc not in {host, settings.public_host}:
                return JSONResponse({"detail": "Cross-origin request refused"}, 403)
            length = request.headers.get("content-length")
            if length is None:
                return JSONResponse({"detail": "Content-Length required"}, 411)
            if int(length) > max_body:
                return JSONResponse(
                    {"detail": f"Upload exceeds {settings.max_upload_mb} MB"}, 413)
        response = await call_next(request)
        response.headers.update(SECURITY_HEADERS)
        if is_public:
            response.headers.update(NOINDEX_HEADER)
        return response

    def attachment(filename):
        return {"Content-Disposition": f'attachment; filename="{filename}"'}

    def job_or_404(job_id):
        job = store.get_job(job_id) if ID_PATTERN.match(job_id) else None
        if job is None:
            raise HTTPException(404, "Job not found")
        return job

    def finished_job(job_id):
        job = job_or_404(job_id)
        if job["status"] not in {"done", "found"}:
            raise HTTPException(409, f"Job is {job['status']}")
        return job

    def genome_or_404(genome_hash):
        row = store.index.lookup(genome_hash) if HASH_PATTERN.match(genome_hash) else None
        if row is None or not row["visible"]:
            raise HTTPException(404, "Not found")
        return {**row, **reference_db.genome_reference_info(row["name"])}

    @app.get("/api/version")
    def version():
        return store.provenance

    @app.get("/api/reference-db")
    def get_reference_db():
        db = store.reference_db
        links = reference_db.strain_job_links(store.finished_jobs())
        return {**db, "wp": [{**row, "job_id": links.get(row["strain"])} for row in db["wp"]]}

    @app.get("/robots.txt", response_class=PlainTextResponse)
    def robots():
        return "User-agent: *\nDisallow: /\n"

    @app.get("/")
    def index():
        return FileResponse(STATIC / "index.html")

    @app.get("/public")
    def public_index():
        return FileResponse(STATIC / "public.html")

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

    def batch_overview(jobs):
        overview = []
        for job in jobs:
            if job["status"] not in {"done", "found"}:
                continue
            output = store.output_dir(job)
            overview.append({
                "id": job["id"], "name": job["name"], "genome_hash": job["genome_hash"],
                "visible": job["visible"], "contigs": job["contigs"], "bases": job["bases"],
                "composition": read_records(output / "genospecies_composition.tsv"),
                "replicons": replicon_totals(
                    read_records(output / "summary_best_hits.tsv", SUMMARY_COLUMNS)),
            })
        return overview

    @app.get("/api/batches/{batch_id}/overview")
    def get_overview(batch_id: str):
        return batch_overview(get_batch(batch_id)["jobs"])

    @app.get("/api/jobs/{job_id}")
    def get_job(job_id: str):
        return job_or_404(job_id)

    class VisibilityIn(BaseModel):
        visible: bool

    @app.post("/api/jobs/{job_id}/visibility")
    def set_visibility(job_id: str, payload: VisibilityIn, request: Request):
        job = finished_job(job_id)
        claims = getattr(request.state, "access_claims", None)
        by = claims.get("email") if claims else None
        store.index.set_visible(job["genome_hash"], payload.visible, by=by)
        return {"genome_hash": job["genome_hash"], "visible": payload.visible}

    @app.post("/api/batches/{batch_id}/visibility")
    def set_batch_visibility(batch_id: str, payload: VisibilityIn):
        if not store.set_batch_visible(batch_id, payload.visible):
            raise HTTPException(404, "Batch not found")
        return {"batch_id": batch_id, "visible": payload.visible}

    def results_payload(output):
        return {
            "contigs": read_records(output / "summary_best_hits.tsv", SUMMARY_COLUMNS),
            "composition": read_records(output / "genospecies_composition.tsv"),
            "downloads": RESULT_DOWNLOADS,
        }

    def hit_table(output, database, suffix="all"):
        tables = list((output / database / "tables").glob(f"*_{suffix}.tsv"))
        return read_table(tables[0]) if tables else None

    def identity_payload(output):
        hits = hit_table(output, "wp")
        if hits is None or hits.empty:
            return {"gates": None, "plasmid": [], "genospecies": [], "candidates": {}}
        return {
            "gates": {"coverage_pct": WP_MIN_COV_PCT, "covered_bp": WP_MIN_COV_BP},
            "plasmid": with_intervals(identity_cells(hits, "plasmid_name")),
            "genospecies": with_intervals(identity_cells(hits, "genospecies")),
            "candidates": fusion_candidates(hits),
        }

    def contig_payload(output, contig_id):
        summary = read_table(output / "summary_all_calls.tsv")
        rows = summary[summary["contig_id"] == contig_id]
        calls = read_table(output / "summary_best_hits.tsv")
        calls = calls[calls["contig_id"] == contig_id]
        if rows.empty or calls.empty:
            raise HTTPException(404, "Contig not found")
        call, rule = decide_call(rows.iloc[0])
        wp, pf32 = hit_table(output, "wp"), hit_table(output, "pf32")
        wp = wp[(wp["contig_id"] == contig_id) & wp["plasmid_id"].notna()] if wp is not None else None
        pf32 = pf32[(pf32["contig_id"] == contig_id) & pf32["plasmid_id"].notna()] if pf32 is not None else None

        has_wp = wp is not None and not wp.empty
        loci = hit_table(output, "pf32", "multi_loci")
        loci = loci[loci["contig_id"] == contig_id] if loci is not None else pd.DataFrame()
        return {
            "contig": to_records(calls.reindex(columns=SUMMARY_COLUMNS))[0],
            "rule": rule if call == rows.iloc[0].get("final_call") else None,
            "plasmids": with_intervals(identity_cells(wp, "plasmid_name", fold_fusions=False)) if has_wp else [],
            "genospecies": with_intervals(identity_cells(wp, "genospecies")) if has_wp else [],
            "pf32": with_intervals(to_records(pf32.reindex(columns=PF32_COLUMNS))) if pf32 is not None else [],
            "loci": to_records(loci.reindex(columns=LOCUS_COLUMNS)),
            "tiling": tile_contig(wp) if has_wp else [],
        }

    def download_payload(output, name, stem):
        if name == "renamed.fasta":
            matches = list(output.glob("*_renamed.fasta"))
            if not matches:
                raise HTTPException(404, "File not found")
            return FileResponse(matches[0], media_type="text/plain",
                                filename=f"{stem}_renamed.fasta")
        if name not in DOWNLOADS:
            raise HTTPException(404, "File not found")
        if name.endswith(".tsv"):
            return Response(guarded_tsv(output / name), media_type=DOWNLOADS[name],
                            headers=attachment(f"{stem}_{name}"))
        return FileResponse(output / name, media_type=DOWNLOADS[name],
                            filename=f"{stem}_{name}")

    @app.get("/api/jobs/{job_id}/results")
    def get_results(job_id: str):
        job = finished_job(job_id)
        return {"job": job, **results_payload(store.output_dir(job))}

    @app.get("/api/jobs/{job_id}/identity")
    def get_identity(job_id: str):
        return identity_payload(store.output_dir(finished_job(job_id)))

    @app.get("/api/jobs/{job_id}/contig")
    def get_contig(job_id: str, contig_id: str):
        job = finished_job(job_id)
        return {"job": job, **contig_payload(store.output_dir(job), contig_id)}

    @app.get("/api/jobs/{job_id}/contig/{contig_id}/fasta", response_class=PlainTextResponse)
    def get_contig_fasta(job_id: str, contig_id: str):
        matches = list(store.output_dir(finished_job(job_id)).glob("*_renamed.fasta"))
        fasta = extract_contig_fasta(matches[0], contig_id) if matches else None
        if fasta is None:
            raise HTTPException(404, "Contig not found")
        return fasta

    @app.get("/api/jobs/{job_id}/download/{name}")
    def download(job_id: str, name: str):
        job = finished_job(job_id)
        stem = Path(job["name"]).name.split(".")[0]
        if name == "provenance.json":
            return JSONResponse({"job": job["id"], "finished_at": job["finished_at"],
                                 **(job.get("provenance") or {})},
                                headers=attachment(f"{stem}_{name}"))
        return download_payload(store.output_dir(job), name, stem)

    def visible_batch_or_404(batch_id):
        batch = store.get_batch(batch_id) if ID_PATTERN.match(batch_id) else None
        if batch is None or not batch["visible"]:
            raise HTTPException(404, "Not found")
        return batch

    @app.get("/api/public/batches")
    def list_public_batches():
        return store.list_visible_batches()

    @app.get("/api/public/batches/{batch_id}")
    def get_public_batch(batch_id: str):
        return visible_batch_or_404(batch_id)

    @app.get("/api/public/batches/{batch_id}/overview")
    def get_public_overview(batch_id: str):
        catalog = reference_db.assembly_catalog()
        overview = batch_overview(visible_batch_or_404(batch_id)["jobs"])
        return [{**row, **reference_db.genome_reference_info(row["name"], catalog)} for row in overview]

    @app.get("/api/public/genomes")
    def list_public_genomes():
        catalog = reference_db.assembly_catalog()
        return [{**row, **reference_db.genome_reference_info(row["name"], catalog)}
                for row in store.index.list_visible()]

    @app.get("/api/public/genomes/{genome_hash}")
    def get_public_genome(genome_hash: str):
        row = genome_or_404(genome_hash)
        return {"genome": row, **results_payload(store.output_dir_for_hash(genome_hash))}

    @app.get("/api/public/genomes/{genome_hash}/identity")
    def get_public_identity(genome_hash: str):
        genome_or_404(genome_hash)
        return identity_payload(store.output_dir_for_hash(genome_hash))

    @app.get("/api/public/genomes/{genome_hash}/contig")
    def get_public_contig(genome_hash: str, contig_id: str):
        row = genome_or_404(genome_hash)
        return {"genome": row, **contig_payload(store.output_dir_for_hash(genome_hash), contig_id)}

    @app.get("/api/public/genomes/{genome_hash}/contig/{contig_id}/fasta", response_class=PlainTextResponse)
    def get_public_contig_fasta(genome_hash: str, contig_id: str):
        genome_or_404(genome_hash)
        matches = list(store.output_dir_for_hash(genome_hash).glob("*_renamed.fasta"))
        fasta = extract_contig_fasta(matches[0], contig_id) if matches else None
        if fasta is None:
            raise HTTPException(404, "Contig not found")
        return fasta

    @app.get("/api/public/genomes/{genome_hash}/download/{name}")
    def download_public(genome_hash: str, name: str):
        row = genome_or_404(genome_hash)
        stem = Path(row["name"]).name.split(".")[0]
        if name == "provenance.json":
            return JSONResponse({"job": row["job_id"], "finished_at": row["created_at"],
                                 **(row.get("provenance") or {})},
                                headers=attachment(f"{stem}_{name}"))
        return download_payload(store.output_dir_for_hash(genome_hash), name, stem)

    app.mount("/static", StaticFiles(directory=STATIC), name="static")
    if settings.access:
        # Added last, so it runs first and covers every route except the public surface.
        install_access_auth(app, exempt_paths=PUBLIC_PATHS, exempt_prefixes=frozenset({PUBLIC_PREFIX}))
    return app
