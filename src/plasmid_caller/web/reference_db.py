"""What's packed into the shipped reference database: every WP and PF32 reference."""
import csv
import logging
import pickle
import re
from pathlib import Path

from ..blast_manager import BlastManager
from ..main import parse_hit_id

log = logging.getLogger(__name__)
DB_ROOT = Path(__file__).resolve().parent.parent / "db"
# new_db/ is a build-time directory, not shipped with the installed package -- this is a
# dev-only convenience (linking a called reference strain back to its result) that degrades
# to no links wherever that manifest isn't present, e.g. in a deployed container.
ASSEMBLIES_MANIFEST = Path(__file__).resolve().parents[3] / "new_db" / "assemblies.tsv"
ASSEMBLY_ACCESSION = re.compile(r"^(GC[AF]_\d+\.\d+)")


def strain_accessions(manifest=ASSEMBLIES_MANIFEST):
    """strain -> assembly accession, from the database build manifest."""
    manifest = Path(manifest)
    if not manifest.exists():
        return {}
    with open(manifest, newline="") as handle:
        return {row["strain"]: row["assembly_accession"].strip()
                for row in csv.DictReader(handle, delimiter="\t") if row["assembly_accession"].strip()}


def assembly_catalog(manifest=ASSEMBLIES_MANIFEST):
    """assembly_accession -> {species, strain}, from the database build manifest."""
    manifest = Path(manifest)
    if not manifest.exists():
        return {}
    with open(manifest, newline="") as handle:
        return {row["assembly_accession"].strip(): {"species": row["species"], "strain": row["strain"]}
                for row in csv.DictReader(handle, delimiter="\t") if row["assembly_accession"].strip()}


def genome_reference_info(name, catalog=None):
    """species/strain/accession_id for a called genome's filename, if it matches a reference assembly."""
    catalog = assembly_catalog() if catalog is None else catalog
    match = ASSEMBLY_ACCESSION.match(name)
    accession = match.group(1) if match else None
    info = catalog.get(accession, {}) if accession else {}
    return {"species": info.get("species"), "strain": info.get("strain"), "accession_id": accession}


def strain_job_links(jobs, manifest=ASSEMBLIES_MANIFEST):
    """strain -> job_id, by joining the manifest's accessions against called jobs named after one."""
    accessions = strain_accessions(manifest)
    if not accessions:
        return {}
    by_accession = {}
    for job in jobs:
        match = ASSEMBLY_ACCESSION.match(job["name"])
        if match:
            by_accession.setdefault(match.group(1), job["id"])
    return {strain: by_accession[acc] for strain, acc in accessions.items() if acc in by_accession}


def _wp_references(root):
    with open(root / "blast_parsing_dict.pkl", "rb") as handle:
        parsing_dict = pickle.load(handle)
    rows = []
    for accession, info in parsing_dict.items():
        candidates = info.get("genospecies_candidates") or []
        rows.append({
            "accession": accession, "name": info.get("name"), "strain": info.get("strain"),
            "genospecies": info.get("genospecies"), "length": info.get("length"),
            "genospecies_candidates": ";".join(candidates),
        })
    return sorted(rows, key=lambda row: (row["genospecies"] or "", row["name"] or "", row["strain"] or ""))


def _pf32_references(root):
    db_path = root / "pf32" / "pf32"
    result = BlastManager().run_blast_command("blastdbcmd", db=str(db_path), entry="all", outfmt="%i")
    rows = []
    for line in result.stdout.splitlines():
        hit_id = line.strip()
        if not hit_id:
            continue
        strain, name = parse_hit_id(hit_id)
        rows.append({"accession": hit_id, "name": name, "strain": strain})
    return sorted(rows, key=lambda row: (row["name"] or "", row["strain"] or ""))


def _attempt(function, *args):
    try:
        return function(*args)
    except Exception:
        log.exception("Could not load %s", function.__name__)
        return []


def load(database=None):
    root = Path(database) if database else DB_ROOT
    return {"wp": _attempt(_wp_references, root), "pf32": _attempt(_pf32_references, root)}
