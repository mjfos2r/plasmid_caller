# plasmid_caller.py
# Version: 6.0.0
# May 20, 2025
# - Michael J. Foster
# - Ben Kotzen

import argparse
import glob
import importlib.util
import json
import os
import pickle
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from . import __about__ as about

__version__ = about.__version__

import pandas
from Bio import SeqIO
from Bio.Blast import NCBIXML
from intervaltree import IntervalTree

from .blast_manager import BlastManager
from .scoring import best_pf32_hit, best_wp_hit, choose_final_call, summarise_multiple_pf32

# init blast_manager
blast_manager = BlastManager()
os.environ["PATH"] = f"{blast_manager.blast_path}:os.environ['PATH']"

#todo: convert to logging, optimize visually
# convert paths to pathlib.Path()s
# write better tests for the scoring.
# also finish HMM integration
# also split code into submodules and tidy things up

def get_input_files(input_path, input_extension):
    """from an input path, return a list of filepaths to each input file"""
    files = []
    files.extend(Path(input_path).glob(f"*.{input_extension}"))
    return files

def _existing_file(path_str: str) -> Path:
    p = Path(path_str)
    if not p.is_file():
        raise argparse.ArgumentTypeError(f"Input FASTA not found: {p}")
    return p

def get_output_path(results_dir):
    """ensure the output path exists and create it if absent."""
    cwd = Path(os.getcwd())
    output_path = os.path.join(cwd, results_dir)
    if not os.path.exists(output_path):
        print(f"Creating output path: {output_path}")
        os.makedirs(output_path, exist_ok=True)
    return output_path


def _db_path(path_str: str = None) -> Path:
    """
    Find the database directory within the package installation.
    *OR*
    convert the passed path and return it.

    Returns:
        Path: Path to the database directory.
    """
    if path_str != "def_db":
        return Path(path_str)

    try:
        # print("Looking for bundled databases: [ 'wp', 'pf32' ]")
        spec = importlib.util.find_spec("plasmid_caller")

        if spec and spec.origin:
            package_dir = Path(spec.origin).parent
            db_path = package_dir / "db"
            # print(f"Looking in: {package_dir}")

            if db_path.exists():
                # print(f"Found DB directory: {db_path}")
                db_path_contents = os.listdir(db_path)
                # print(f"Contents: {db_path_contents}")
                return db_path

    except (ImportError, AttributeError):
        pass

    # if package detection failed, check if running from src
    print(
        "Could not detect package location. We're going to look relative to the location of this file."
    )
    current_dir = Path(__file__).parent.absolute()
    print(f"Current directory: {current_dir}")
    # is db in the current source directory structure?
    repo_db_path = current_dir.parent / "db"
    if repo_db_path.exists():
        repo_db_contents = os.listdir(repo_db_path)
        print("Found db directory relative to this file!")
        print(f"Path: {repo_db_path}")
        print(f"Contents: {repo_db_contents}")
        return repo_db_path
    print(
        "Could not find db relative to package location or relative to this file. Defaulting to the container path: '/db'"
    )
    return Path("/db")  # if all else fails return the container path. (BAD)


# this is when dbs exist and are of proper version.
def get_db_type(db_dir, quiet=True):
    """Get the database name and type for execution"""
    db_dir = db_dir.resolve() if isinstance(db_dir, Path) else Path(db_dir).resolve()
    print(f"DB_Dir: {db_dir}")
    db_dir_contents = os.listdir(db_dir)
    print(db_dir_contents)
    try:
        print(f"Getting database information from databases contained in: {db_dir}")
        print(f"Contents: {db_dir_contents}")
        result = blast_manager.run_blast_command(
            "blastdbcmd",
            list=db_dir,
            recursive=True,
        )
    except Exception as exc:
        if not quiet:
            print(f"[get_db_type] blastdbcmd failed: {exc}")
        return []

    db_prog_pairs: list[tuple[str, str]] = []

    for line in (
        result.stdout.splitlines()
    ):  # [:-1]:  # last line is empty so don't iterate to it!
        if not line.strip():
            continue

        parts = line.split()
        if len(parts) < 2:
            if not quiet:
                print(f"[get_db_type] Unrecognized line: {line!r}")
            continue

        name, db_type = parts[0], parts[1].lower()

        if db_type == "protein":
            program = "blastx"
        elif db_type == "nucleotide":
            program = "blastn"
        else:
            if not quiet:
                print(f"[get_db_type] Unknown DB type {db_type!r} for {name}")
            continue
        db_prog_pairs.append((name, program))

        if not quiet:
            print(f"[get_db_type] {name:<40} -> {program}")
    return db_prog_pairs


def get_blast_command(prog, input_file, output_path, database, threads):
    """Construct a params dict from the given input for use in run_blast()"""
    file_id = Path(input_file).stem
    output_file = f"{output_path}/{file_id}_blast_results.xml"
    return {
        "program": prog,
        "query": input_file,
        "db": database,
        "out": output_file,
        "evalue": "1e-100",
        "num_threads": threads,
        "outfmt": "5",
        "max_target_seqs": "5",
        "max_hsps": "10",
    }


def run_blast(params):
    """Run a blast command using blast_manager"""
    try:
        if isinstance(params, dict):
            program = params.pop(
                "program", "blastn"
            )  # blastn if we don't have one.. FIX
            return blast_manager.run_blast_command(program, **params)
        return (
            f"SUCCESS: {' '.join(command)}" if result.returncode == 0 else result.stderr
        )
    except subprocess.CalledProcessError as e:
        return f"FAILURE: '{command}'\nerror: {e.stderr}"
    except Exception as e:
        return f"An unexpected error occurred: {str(e)}"


def get_name_from_acc(hit_id, parsing_dict):
    """this is how strain and plasmid name are parsed using the previously assembled parsing dictionary"""
    acc_id = hit_id.split("|")[1]
    name = parsing_dict[acc_id]["name"]
    strain = parsing_dict[acc_id]["strain"]
    return strain, name


def parse_hit_id(hit_id):
    """this is how strain and plasmid name are parsed using the hit_id for pf32 hits"""
    gene_id_list = hit_id.strip().split("_")
    if len(gene_id_list) == 2:
        strain = gene_id_list[0]
    elif len(gene_id_list) == 1:  # this is the weird pdb case.
        strain = gene_id_list[0].split("|")[1]
    elif (
        gene_id_list[0] == "NE"
    ):  # this is the NE_1234 strains single case, hate this too.
        strain = "_".join(gene_id_list[0:2])
    else:
        strain = gene_id_list[0]
    pattern = r"chromosome|(?:lp|cp)\d{1,2}(?:[-+#](?:lp|cp)?\d{1,2})*"
    match = re.findall(pattern, hit_id.replace("_", "-"))
    plasmid_id = (
        match[0]
        if match
        else "HO14_NovelPFam32"
        if "HO14_NovelPFam32_DatabaseShortContig" in hit_id
        else None
    )  # bad hardcode. needs to be dropped.
    return strain, plasmid_id


def get_contig_len(input_file, contig_id):
    """Get contig length directly from the input FASTA file"""
    fasta_file = input_file
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id.split()[0] == contig_id.split()[0]:
            return len(record.seq)
        else:
            continue

def calculate_percent_identity_and_coverage(alignment):
    total_identities = 0
    total_alignment_length = 0
    query_intervals = IntervalTree()
    subject_intervals = IntervalTree()

    for hsp in alignment.hsps:
        total_identities += hsp.identities
        total_alignment_length += hsp.align_length

        # Add HSP interval to the query and subject intervals
        query_start, query_end = sorted([hsp.query_start, hsp.query_end])
        subject_start, subject_end = sorted([hsp.sbjct_start, hsp.sbjct_end])

        query_intervals.addi(query_start, query_end)
        subject_intervals.addi(subject_start, subject_end)

    # Merge overlapping intervals
    query_intervals.merge_overlaps()
    subject_intervals.merge_overlaps()

    # get total length now!
    ref_covered_length = sum(interval.length() for interval in subject_intervals)
    contig_covered_length = sum(interval.length() for interval in query_intervals)

    # Calculate overall percent identity
    overall_percent_identity = (
        (total_identities / total_alignment_length) * 100
        if total_alignment_length != 0
        else 0
    )

    return {
        "overall_percent_identity": overall_percent_identity,
        "query_covered_length": contig_covered_length,
        "ref_covered_length": ref_covered_length,
        "covered_intervals": [
            (interval.begin, interval.end) for interval in subject_intervals
        ],
        "query_intervals": [
            (interval.begin, interval.end) for interval in query_intervals
        ],
        "subject_hit_coords": [
            (hsp.sbjct_start, hsp.sbjct_end) for hsp in alignment.hsps
        ],
    }


def parse_blast_xml(xml_file, args, *, parsing_type: str, dbs_dir: Path) -> dict:
    """
    rewritten XML parsing function.
    """
    # ---------------- prep ----------------
    assembly_id = Path(xml_file).stem.replace("_blast_results", "")
    dbs_dir     = Path(dbs_dir)
    parsing_dict = pickle.load(open(dbs_dir / "blast_parsing_dict.pkl", "rb"))

    # master field list ( += contig_len which was missing in the old list)
    KEYS = [
        "assembly_id", "contig_id", "contig_len",
        "plasmid_id",  "plasmid_name", "strain",
        "query_length", "ref_length",
        "overall_percent_identity",
        "query_coverage_percent",      # renamed to match later code
        "query_covered_length", "ref_covered_length",
        "covered_intervals", "query_intervals", "subject_hit_coords",
    ]
    NA = pandas.NA   # convenient alias

    # ---------------- parse ----------------
    blast_hits: dict[str, list[dict]] = defaultdict(list)

    with open(xml_file) as handle:
        for record in NCBIXML.parse(handle):
            contig_id      = record.query
            contig_length  = record.query_length
            contig_len     = get_contig_len(args.input, contig_id)

            if not record.alignments:
                blast_hits[contig_id].append(
                    {k: NA for k in KEYS} | {"assembly_id": assembly_id,
                                            "contig_id": contig_id,
                                            "contig_len": contig_len,
                                            "query_length": contig_length}
                )
                continue

            for aln in record.alignments:
                plasmid_id = aln.hit_id
                ref_length = aln.length

                if parsing_type == "wp":
                    strain, plasmid_name = get_name_from_acc(plasmid_id, parsing_dict)
                elif parsing_type == "pf32":
                    strain, plasmid_name = parse_hit_id(plasmid_id)
                else:
                    strain, plasmid_name = "unknown", plasmid_id

                # metrics derived from HSPs
                aln_stats = calculate_percent_identity_and_coverage(aln)
                aln_stats["query_coverage_percent"] = (
                    aln_stats["query_covered_length"] / contig_length * 100
                    if contig_length else NA
                )

                row = {
                    "assembly_id": assembly_id,
                    "contig_id":   contig_id,
                    "contig_len":  contig_len,
                    "plasmid_id":  plasmid_id,
                    "plasmid_name": plasmid_name,
                    "strain":       strain,
                    "query_length": contig_length,
                    "ref_length":   ref_length,
                } | aln_stats

                # ensure *all* KEYS exist
                complete_row = {k: row.get(k, NA) for k in KEYS}
                blast_hits[contig_id].append(complete_row)

    return blast_hits

def get_results(results_dir, quiet=False):
    if not quiet:
        print(os.listdir(results_dir))
    return glob.glob(f"{results_dir}/**/*.xml", recursive=True)


def get_hits_table(hits):
    rows = []
    for _, hits in hits.items():
        for hit in hits:
            rows.append(hit)
    df = pandas.DataFrame(rows)
    return df


def sanitize_fa_headers(input_fa, output_dir):
    tdir = output_dir / "sanitized"
    tdir.mkdir(exist_ok=True)
    sanitized_fa = tdir / f"{Path(input_fa).stem}.fasta"
    with open(sanitized_fa, 'w') as handle:
        for record in SeqIO.parse(input_fa, 'fasta'):
            record.description = record.id # this is required for the tables to be readable.
            SeqIO.write(record, handle, 'fasta')
    return sanitized_fa


def parse_to_tsv(file_id, xml_file, full_table_path, args, parsing_type, db_path, best_table_path, tables_dir):
    hits = parse_blast_xml(xml_file, args, parsing_type=parsing_type, dbs_dir=db_path)
    full_hits_df = get_hits_table(hits)
    full_hits_df.to_csv(full_table_path, sep="\t", index=False)

    if parsing_type == "pf32":
        best_hits_df = best_pf32_hit(full_hits_df)
        multi_best_df, concat_series = summarise_multiple_pf32(full_hits_df)
        multi_table = tables_dir / f"{file_id}_multi_loci.tsv"
        if not multi_best_df.empty:
            multi_best_df.to_csv(multi_table, sep='\t', index=False)
            if not args.quiet:
                print(f"[pf32] wrote multi-PF32 loci table → {multi_table.name}")
        # attach concatenated name with all loci + loci flag to best_hits_df
        best_hits_df = (
            best_hits_df.set_index("contig_id")
                        .join(concat_series, how="left")
                        .assign(multiple_loci=lambda df: df["concat_call"].notna())
                        .reset_index()
        )

    elif parsing_type == "wp":
        best_hits_df = best_wp_hit(full_hits_df)

    best_hits_df.to_csv(best_table_path, sep="\t", index=False)

    if not args.quiet:
        print(
            f"[{parsing_type}] wrote: "
            f"{full_table_path.name} ({len(full_hits_df)} rows) | "
            f"{best_table_path.name} ({len(best_hits_df)} rows)"
        )

    return best_hits_df

def rename_fasta_headers(input_fa, output_fa, mapping):
    """
    Rewrite fasta headers of the input file using the newly determined classifications.
    header becomes:
        >mapping[old_id]
        NNNNNNHiNNNNNNNNNNNNNNN....
    """
    with open(output_fa, "w") as handle_out:
        for rec in SeqIO.parse(input_fa, "fasta"):
            old_id = rec.description
            rec.id = mapping.get(old_id, old_id)
            SeqIO.write(rec, handle_out, "fasta")


class FullVersion(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        print(f"{parser.prog} {__version__}")
        print(BlastManager().versions)
        parser.exit(0)


## Define main function logic.
def main(args=None):
    parser = argparse.ArgumentParser(
        prog="plasmid_caller",
        description="Plasmid/replicon classifier for a single FASTA file containing B.burgdorferi contigs.",
    )
    # Add the arguments
    parser.add_argument(
        "-v",
        "--version",
        nargs=0,
        action=FullVersion,
        help="Show program and BLAST versions, then exit.",
    )
    parser.add_argument(
        "-i",
        "--input",
        type=_existing_file,
        required=True,
        help="Input FASTA file path"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output directory"
    )
    parser.add_argument(
        "-t",
        "--threads",
        required=False,
        type=int,
        help="How many cores to use?",
        default=4
    )
    parser.add_argument(
        "-db",
        "--database",
        required=False,
        type=_db_path,
        help="Path to directory containing blast databases (default: built-in dbs)",
        default="def_db"
    )
    parser.add_argument(
        "--skip_blast",
        action="store_true",
        help="Skip running BLAST and only parse existing results",
        default=False,
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Run without terminal output for workflow integration",
        default=False,
    )

    # Parse the arguments
    args = parser.parse_args()

    # validate database location
    if not args.database.exists():
        parser.error(f"Database directory not found!")

        # now lets get our dbs to run against
    dbs_dir = args.database
    dbs = get_db_type(dbs_dir)
    print(f"Databases to run against: {dbs}")

    if not args.quiet:
        print(f"Input file: {args.input}")
        print(f"Output directory: {args.output}")
        print(f"Job threads: {args.threads}")
        print(f"Databases: {args.database}")
        print(f"Skip BLAST?: {args.skip_blast}")

    args.output.mkdir(parents=True, exist_ok=True)

    file_id = Path(args.input).stem

    combined_summary_parts = []

    for db in dbs:
        db_path = db[0]
        prog = db[1]
        db_name = Path(db_path).stem
        results_dir = args.output / db_name / "xml_files"
        tables_dir = args.output / db_name / "tables"
        results_dir.mkdir(parents=True, exist_ok=True)
        tables_dir.mkdir(parents=True, exist_ok=True)
        sanitized_fa= sanitize_fa_headers(args.input)
        if not args.skip_blast:
            blast_params = get_blast_command(prog, sanitized_fa, results_dir, db_path, args.threads)
            # Run BLAST without parallelization since there's only one input file
            run_blast(blast_params)

        xml_file = results_dir / f"{file_id}_blast_results.xml"
        full_table = tables_dir / f"{file_id}_all.tsv"
        best_table = tables_dir / f"{file_id}_best.tsv"
        best_df = parse_to_tsv(
            file_id=file_id,
            xml_file=xml_file,
            full_table_path=full_table,
            args=args,
            parsing_type=db_name,
            db_path=dbs_dir,
            best_table_path=best_table,
            tables_dir=tables_dir,
        )
        tag = f"_{db_name}"
        tagged_best_df = best_df.add_suffix(tag).rename(
                columns={f"contig_id{tag}": "contig_id", f"contig_len{tag}": "contig_len", f"assembly_id{tag}": "assembly_id"}
        )

        combined_summary_parts.append(tagged_best_df)

    summary_path = args.output / "summary_best_hits.tsv"

    frames = [df.set_index("contig_id") for df in combined_summary_parts]
    summary_df = pandas.concat(frames, axis=1, join="outer").reset_index()
    summary_df = summary_df.loc[:, ~summary_df.columns.duplicated(keep="first")]
    summary_df["final_call"] = summary_df.apply(choose_final_call, axis=1)
    summary_df.to_csv(summary_path, sep="\t", index=False)

    best_map = summary_df.set_index("contig_id")["final_call"].to_dict()
    json_path = args.output / "summary_best_hits.json"
    json_path.write_text(json.dumps(best_map, indent=4))

    renamed_fa = args.output / f"{file_id}_renamed.fasta"
    rename_fasta_headers(Path(args.input), renamed_fa, best_map)

    if not args.quiet:
        print(f"Wrote combined summary -> {summary_path}")
        print(f"Wrote dictionary of final calls -> {json_path}")
        print(f"Wrote renamed fasta file -> {renamed_fa}")

    return 0


if __name__ == "__main__":
    main()
