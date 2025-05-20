# plasmid_caller.py
# Version: 6.0.0
# March 20, 2025
# - Michael J. Foster
# - Ben Kotzen
# Modified to accept a single fasta file as input

import argparse
import glob
import importlib.util
import os
import pickle
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from plasmid_caller import __about__ as about

__version__ = about.__version__

import pandas
from Bio import SeqIO
from Bio.Blast import NCBIXML
from intervaltree import IntervalTree

from .blast_manager import BlastManager

blast_manager = BlastManager()
os.environ["PATH"] = f"{blast_manager.blast_path}:os.environ['PATH']"


def get_input_files(input_path, input_extension):
    """from an input path, return a list of filepaths to each input file"""
    files = []
    files.extend(Path(input_path).glob(f"*.{input_extension}"))
    return files


def get_output_path(results_dir):
    """ensure the output path exists and create it if absent."""
    cwd = Path(os.getcwd())
    output_path = os.path.join(cwd, results_dir)
    if not os.path.exists(output_path):
        print(f"Creating output path: {output_path}")
        os.makedirs(output_path, exist_ok=True)
    return output_path


def get_default_db_path():
    """
    Find the database directory within the package installation.

    Returns:
        Path: Path to the database directory.
    """
    try:
        print("Looking for bundled databases: [ 'wp', 'pf32' ]")
        spec = importlib.util.find_spec("plasmid_caller")

        if spec and spec.origin:
            package_dir = Path(spec.origin).parent
            db_path = package_dir / "db"
            print(f"Looking in: {package_dir}")

            if db_path.exists():
                print(f"Found DB directory: {db_path}")
                db_path_contents = os.listdir(db_path)
                print(f"Contents: {db_path_contents}")
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
    db_dir = Path(db_dir).expanduser().resolve()
    db_dir_contents = os.listdir(db_dir)
    try:
        print(f"Getting database information from databases contained in: {db_dir}")
        print(f"Contents: {db_dir_contents}")
        result = blast_manager.run_blast_command(
            "blastdbcmd",
            list=str(db_dir),
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


def get_results(results_dir, quiet=False):
    if not quiet:
        print(os.listdir(results_dir))
    return glob.glob(f"{results_dir}/**/*.xml", recursive=True)


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
    return 1  # TODO Raise exception. (y tho?)


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

    # Merge overlapping intervals to calculate the covered length
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


def parse_blast_xml(xml_file, args, **kwargs):
    """This is the new and improved blast results parsing function.
    - it takes kwargs for parsing_type which is either 'wp' or 'pf32'"""
    assembly_id = Path(xml_file).stem.replace("_blast_results", "")
    parsing_type = kwargs.get("parsing_type", "general")
    dbs_dir = kwargs.get("dbs_dir", "/app")
    dbs_dir = Path(dbs_dir) if isinstance(dbs_dir, str) else dbs_dir
    parsing_pkl = dbs_dir / "blast_parsing_dict.pkl"
    parsing_dict = pickle.load(open(parsing_pkl, "rb"))

    with open(xml_file, "r") as handle:
        records = NCBIXML.parse(handle)

        # set up dict for intervals
        blast_hits = defaultdict(dict)

        for record in records:
            keys = [
                "assembly_id",
                "contig_id",
                "plasmid_id",
                "plasmid_name",
                "strain",
                "query_length",
                "ref_length",
                "overall_percent_identity",
                "query_coverage_percentage",
                "query_covered_length",
                "ref_covered_length",
                "covered_intervals",
                "query_intervals",
                "subject_hit_coords",
            ]

            contig_id = record.query
            contig_length = record.query_length

            if contig_id not in blast_hits:
                blast_hits[contig_id] = []
            hits = blast_hits[contig_id]

            if len(record.alignments) == 0:
                hits = dict(zip(keys, "NaN" * len(keys)))
            for alignment in record.alignments:
                plasmid_id = alignment.hit_id
                ref_length = alignment.length

                if parsing_type == "wp":
                    strain, plasmid_name = get_name_from_acc(plasmid_id, parsing_dict)
                elif parsing_type == "pf32":
                    strain, plasmid_name = parse_hit_id(plasmid_id)
                else:
                    strain = "unknown"
                    plasmid_name = alignment.hit_id

                results = {
                    "assembly_id": assembly_id,
                    "contig_id": contig_id,
                    "contig_len": get_contig_len(args, contig_id),
                    "plasmid_id": plasmid_id,
                    "plasmid_name": plasmid_name,
                    "strain": strain,
                    "query_length": contig_length,
                    "ref_length": ref_length,
                }
                alignment_results = calculate_percent_identity_and_coverage(alignment)
                alignment_results["query_coverage_percent"] = (
                    alignment_results["query_covered_length"] / contig_length * 100
                )
                results.update(alignment_results)
                hits.append(results)
    return blast_hits


def get_hits_table(hits):
    rows = []
    for _, hits in hits.items():
        for hit in hits:
            rows.append(hit)
    df = pandas.DataFrame(rows)
    return df


def parse_to_tsv(file, output_path, args, parsing_type, db_path):
    hits = parse_blast_xml(file, args, parsing_type=parsing_type, dbs_dir=db_path)
    hits_df = get_hits_table(hits)
    hits_df.to_csv(str(output_path), sep="\t", index=False)
    message = f"{parsing_type} table for: {Path(file).stem.split('_')[0]} written to {output_path}"
    if not args.quiet:
        print(message)
    return message


class FullVersion(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        print(f"{parser.prog} {__version__}")
        print("###BLAST VERSIONS###")
        print(blast_manager.versions)
        print("###---###")


## Define main function logic.
def main(args=None):
    # get this setup before argparse needs it.
    default_db_path = get_default_db_path()
    # Create the parser
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
        required=True,
        type=str,
        help="Input FASTA file path"
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=str,
        help="The directory for outputs"
    )
    parser.add_argument(
        "-t",
        "--threads",
        required=False,
        type=int,
        help="How many cores to use?"
    )
    parser.add_argument(
        "-db",
        "--database",
        required=False,
        type=str,
        help="Path to directory containing blast databases",
        default=default_db_path,
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

    if len(sys.argv) == 1:
        parser.print_help()
        parser.error("Please specify relevant arguments!")
    else:
        # Parse the arguments
        args = parser.parse_args()
        # validate database location
        db_path = Path(args.db)
        if not db_path.exists():
            parser.error(f"Database directory not found: {db_path}")

        if not args.quiet:
            print(args)

    # Verify input file exists
    if not os.path.exists(args.input):
        parser.error("Input FASTA file does not exist!")

    # Get output paths
    output_path = Path(args.output)

    if not args.quiet:
        print(f"Input file: {args.input}")
        print(f"Output directory: {args.output}")
        print(f"Job threads: {args.threads}")
        print(f"Databases: {args.db}")
        print(f"Skip BLAST?: {args.skip_blast}")

    # now lets get our dbs to run against
    dbs_dir = args.db
    dbs = get_db_type(dbs_dir)
    print(f"Databases to run against: {dbs}")

    if not output_path.exists():
        if not args.quiet:
            print("Creating output directory!\n")
        output_path.mkdir(exist_ok=True)

    for db in dbs:
        db_path = db[0]
        prog = db[1]
        db_name = Path(db_path).stem
        results_dir = output_path / db_name / "xml_files"
        tables_dir = output_path / db_name / "tables"
        results_dir.mkdir(exist_ok=True)
        tables_dir.mkdir(exist_ok=True)

        if not args.skip_blast:
            blast_params = get_blast_command(
                prog, args.input, results_dir, db_path, args.threads
            )
            # Run BLAST without parallelization since there's only one input file
            run_blast(blast_params)

        # Parse results
        results = get_results(results_dir, args.quiet)
        for result_file in results:
            results_table = f"{Path(result_file).stem}_table.tsv"
            output_table = tables_dir / results_table
            parse_to_tsv(result_file, output_table, args, db_name, db_path=dbs_dir)

    if not args.quiet:
        print("Finished!")
    return 0


if __name__ == "__main__":
    main()
