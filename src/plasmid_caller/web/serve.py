"""`plasmid_caller serve` command."""
import argparse
import ipaddress
import os
import sys
from pathlib import Path

from .jobs import Settings


def _is_loopback(host):
    if host == "localhost":
        return True
    try:
        return ipaddress.ip_address(host).is_loopback
    except ValueError:
        return False


def serve(argv):
    defaults = Settings(runs_dir=Path("runs"))
    parser = argparse.ArgumentParser(
        prog="plasmid_caller serve",
        description="Web interface for submitting assemblies and exploring results.",
    )
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--runs-dir", type=Path, default=defaults.runs_dir,
                        help="Where uploads and results are kept (default: ./runs)")
    parser.add_argument("--jobs", type=int, default=defaults.jobs,
                        help="Assemblies to run at once")
    parser.add_argument("--threads", type=int, default=defaults.threads,
                        help="BLAST threads per assembly")
    parser.add_argument("-db", "--database", type=Path,
                        help="BLAST database directory (default: built-in dbs)")
    parser.add_argument("--max-upload-mb", type=int, default=defaults.max_upload_mb,
                        help="Largest upload request")
    parser.add_argument("--max-assembly-mb", type=int, default=defaults.max_assembly_mb,
                        help="Largest assembly once decompressed")
    parser.add_argument("--max-queued", type=int, default=defaults.max_queued,
                        help="Most assemblies waiting or running")
    parser.add_argument("--job-timeout", type=int, default=defaults.job_timeout,
                        help="Seconds before an assembly is stopped")
    parser.add_argument("--retention-days", type=int, default=defaults.retention_days,
                        help="Delete batches older than this (default: keep forever)")
    args = parser.parse_args(argv)

    token = os.environ.get("PLASMID_CALLER_TOKEN")
    if not token and not _is_loopback(args.host):
        parser.error("Set PLASMID_CALLER_TOKEN before binding to a non-loopback address.")

    try:
        import uvicorn
    except ImportError:
        sys.exit("ERROR: The web server needs extra packages: pip install 'plasmid_caller[web]'")
    from .app import create_app

    settings = Settings(
        runs_dir=args.runs_dir, jobs=args.jobs, threads=args.threads,
        max_upload_mb=args.max_upload_mb, max_assembly_mb=args.max_assembly_mb,
        max_queued=args.max_queued, job_timeout=args.job_timeout,
        retention_days=args.retention_days, database=args.database, token=token,
    )
    url = f"http://{args.host}:{args.port}/" + (f"?token={token}" if token else "")
    print(f"plasmid_caller web: {url}")
    uvicorn.run(create_app(settings), host=args.host, port=args.port, log_level="warning")
    return 0
