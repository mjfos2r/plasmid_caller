"""What produced a result: stamped on every job and served at /api/version."""
import hashlib
import logging
import os
import subprocess
from pathlib import Path

from ..__about__ import __version__
from ..blast_manager import BlastManager

log = logging.getLogger(__name__)
DB_ROOT = Path(__file__).resolve().parent.parent / "db"
DB_PATTERNS = ["wp/wp.*", "pf32/pf32.*", "blast_parsing_dict.pkl", "reference_taxa.tsv"]


def database_hash(root):
    """Hash every file the caller reads, so two equal hashes mean the same reference set."""
    digest = hashlib.sha256()
    files = sorted(path for pattern in DB_PATTERNS for path in root.glob(pattern))
    if not files:
        raise FileNotFoundError(root)
    for path in files:
        digest.update(path.relative_to(root).as_posix().encode())
        with open(path, "rb") as handle:
            while chunk := handle.read(1 << 20):
                digest.update(chunk)
    return digest.hexdigest()


def blast_version():
    blastn = BlastManager().blast_path / "blastn"
    output = subprocess.run([str(blastn), "-version"], capture_output=True, text=True, check=True)
    return output.stdout.splitlines()[0].split(":", 1)[1].strip()


def _attempt(function, *args):
    try:
        return function(*args)
    except Exception:
        log.exception("Could not determine %s", function.__name__)
        return "unknown"


def collect(database=None):
    return {
        "caller_version": __version__,
        "database_sha256": _attempt(database_hash, Path(database) if database else DB_ROOT),
        "blast_version": _attempt(blast_version),
        # Set at image build; a running container cannot read its own digest.
        "source_revision": os.environ.get("SOURCE_REVISION", "unknown"),
        "base_image": os.environ.get("BASE_IMAGE", "unknown"),
        "service_revision": os.environ.get("K_REVISION", "unknown"),
    }
