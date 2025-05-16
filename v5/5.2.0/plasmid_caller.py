# plasmid_caller.py
# Version: 5.2.0
# March 20, 2025
# - Michael J. Foster
# - Ben Kotzen
# Modified to accept a single fasta file as input

import os
import re
import subprocess
import glob
import pandas
import pickle
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from tqdm import tqdm
from pathlib import Path
from intervaltree import Interval, IntervalTree
from concurrent.futures import ProcessPoolExecutor, as_completed

def get_output_path(results_dir):
    cwd = Path(os.getcwd())
    output_path = os.path.join(cwd, results_dir)
    if not os.path.exists(output_path):
        print(f"Creating output path: {output_path}")
        os.makedirs(output_path, exist_ok=True)
    return output_path

def get_blast_command(prog, input_file, output_path, database, threads):
    file_id = Path(input_file).stem
    output_file = f'{output_path}/{file_id}_blast_results.xml'
    command = [prog, '-query', input_file, '-db', database,
               '-out', output_file, '-evalue', '1e-100', '-num_threads', threads,
               '-outfmt', '5', '-max_target_seqs', '5', '-max_hsps', '10']
    return command

def make_blast_db(dbtype, input_file, db_output):
    logfile = f'{Path(input_file).stem}_blastdb_creation.log'
    command = ['makeblastdb','-in', input_file, '-dbtype', dbtype,
               '-parse_seqids', '-logfile', logfile, '-out', db_output]
    return command

def get_db_type(db_dir, quiet=True):
    """ Get the individual database for blast and the database type """
    command = ['blastdbcmd', '-list', f'{db_dir}', '-recursive',]
    results = subprocess.run(command, capture_output=True, text=True)
    dbs = []
    progs = []
    output = results.stdout
    for line in output.split('\n')[:-1]: # last line is empty so don't iterate to it!
        db_info = line.split(' ')
        if not quiet:
            print(db_info)
        db_name = db_info[0]
        db_type = db_info[1].lower()
        dbs.append(db_name)
        if db_type == 'protein':
            progs.append('blastx')
        elif db_type == 'nucleotide':
            progs.append('blastn')
    return list(zip(dbs, progs))

def get_contig_len(args, contig_id):
    """Get contig length directly from the input FASTA file"""
    fasta_file = args.input
    for record in SeqIO.parse(fasta_file, 'fasta'):
        if record.id.split()[0] == contig_id.split()[0]:
            return len(record.seq)
    return 1 # TODO Raise exception.

################################################################################
#:::::::::::::::::::::::::::::: RESULTS PARSING :::::::::::::::::::::::::::::::#
################################################################################

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
    overall_percent_identity = (total_identities / total_alignment_length) * 100 if total_alignment_length != 0 else 0

    return {
        "overall_percent_identity": overall_percent_identity,
        "query_covered_length": contig_covered_length,
        "ref_covered_length": ref_covered_length,
        "covered_intervals": [(interval.begin, interval.end) for interval in subject_intervals],
        "query_intervals": [(interval.begin, interval.end) for interval in query_intervals],
        "subject_hit_coords": [(hsp.sbjct_start, hsp.sbjct_end) for hsp in alignment.hsps]
    }

def get_name_from_acc(hit_id, parsing_dict):
    """ this is how strain and plasmid name are parsed using the previously assembled parsing dictionary"""
    acc_id = hit_id.split('|')[1]
    name = parsing_dict[acc_id]['name']
    strain = parsing_dict[acc_id]['strain']
    return strain, name

def parse_hit_id(hit_id):
    """ this is how strain and plasmid name are parsed using the hit_id for pf32 hits """
    gene_id_list = hit_id.strip().split('_')
    if len(gene_id_list) == 2:
        strain = gene_id_list[0]
    elif len(gene_id_list) == 1: # this is the weird pdb case.
        strain = gene_id_list[0].split('|')[1]
    elif gene_id_list[0] == 'NE': # this is the NE_1234 strains single case, hate this too.
        strain = '_'.join(gene_id_list[0:2])
    else:
        strain = gene_id_list[0]
    pattern = r"chromosome|(?:lp|cp)\d{1,2}(?:[-+#](?:lp|cp)?\d{1,2})*"
    match = re.findall(pattern, hit_id.replace('_', '-'))
    plasmid_id = match[0] if match else 'HO14_NovelPFam32' if 'HO14_NovelPFam32_DatabaseShortContig' in hit_id else None
    return strain, plasmid_id

def parse_blast_xml(xml_file, args, **kwargs):
    """ This is the new and improved blast results parsing function.
    - it takes kwargs for parsing_type which is either 'wp' or 'pf32'"""
    assembly_id = Path(xml_file).stem.replace('_blast_results', '')
    parsing_type = kwargs.get('parsing_type', 'general')

    with open(xml_file, 'r') as handle:
        records = NCBIXML.parse(handle)
        parsing_pkl = '/app/blast_parsing_dict.pkl'
        parsing_dict = pickle.load(open(parsing_pkl, 'rb'))

        # set up dict for intervals
        blast_hits = defaultdict(dict)

        for record in records:
            keys = [
                "assembly_id", "contig_id", "plasmid_id", "plasmid_name", "strain", "query_length",
                "ref_length", "overall_percent_identity", "query_coverage_percentage",
                "query_covered_length", "ref_covered_length", "covered_intervals", "query_intervals", "subject_hit_coords",
            ]

            contig_id = record.query
            contig_length = record.query_length

            if contig_id not in blast_hits:
                blast_hits[contig_id] = []
            hits = blast_hits[contig_id]

            if len(record.alignments) == 0:
                hits = dict(zip(keys, 'NaN'*len(keys)))
            for alignment in record.alignments:
                plasmid_id = alignment.hit_id
                ref_length = alignment.length

                if parsing_type == 'wp':
                    strain, plasmid_name = get_name_from_acc(plasmid_id, parsing_dict)
                elif parsing_type == 'pf32':
                    strain, plasmid_name = parse_hit_id(plasmid_id)
                else:
                    strain = 'unknown'
                    plasmid_name = alignment.hit_id

                results = {
                    "assembly_id": assembly_id,
                    "contig_id": contig_id,
                    "contig_len": get_contig_len(args, assembly_id, contig_id),
                    "plasmid_id": plasmid_id,
                    "plasmid_name": plasmid_name,
                    "strain": strain,
                    "query_length": contig_length,
                    "ref_length": ref_length,
                }
                alignment_results = calculate_percent_identity_and_coverage(alignment)
                alignment_results['query_coverage_percent'] = alignment_results['query_covered_length'] / contig_length * 100
                results.update(alignment_results)
                hits.append(results)
    return blast_hits

# Function to determine the best match
def get_best_match(matches, key):
    """Highest percent identity takes the cake. specify which feature to compare.
    ex: get_best_match(matches, "percent_identity")"""
    best_match = None
    best_score = -1
    for match in matches:
        if match[key] > best_score:
            best_score = match[key]
            best_match = match
    return best_match

################################################################################
#:::::::::::::::::::::::::::::TABLE OUTPUT STUFF:::::::::::::::::::::::::::::::#
################################################################################

def get_hits_table(hits):
    rows = []
    for _, hits in hits.items():
        for hit in hits:
            rows.append(hit)
    df = pandas.DataFrame(rows)
    return df

def parse_to_tsv(file, output_path, args, parsing_type):
    hits = parse_blast_xml(file, args, parsing_type=parsing_type)
    hits_df = get_hits_table(hits)
    hits_df.to_csv(str(output_path), sep='\t', index=False)
    message = f"{parsing_type} table for: {Path(file).stem.split('_')[0]} written to {output_path}"
    if not args.quiet:
        print(message)
    return message

################################################################################
#:::::::::::::::::::::::::::::: MULTIPROCESSING :::::::::::::::::::::::::::::::#
################################################################################
def run_command(command):
    try:
        command = [str(arg) for arg in command]
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        return f"SUCCESS: {' '.join(command)}" if result.returncode == 0 else result.stderr
    except subprocess.CalledProcessError as e:
        return f"FAILURE: '{command}'\nerror: {e.stderr}"
    except Exception as e:
        return f"An unexpected error occurred: {str(e)}"

def run_blast(prog, database, input_file, results_dir, threads, quiet=False):
    if not quiet:
        print(f'Running {prog} for {input_file}!')
    command = get_blast_command(prog, input_file, results_dir, database, threads)
    return run_command(command)

def get_results(results_dir, quiet=False):
    if not quiet:
        print(os.listdir(results_dir))
    return glob.glob(f'{results_dir}/**/*.xml', recursive=True)

## Define main function logic.
def main(args):
    # Verify input file exists
    if not os.path.exists(args.input):
        parser.error("Input FASTA file does not exist!")

    # Get output paths
    output_path = args.output

    if not args.quiet:
        print(f"Input file: {args.input}")
        print(f"Output directory: {args.output}")
        print(f"Job threads: {args.threads}")
        print(f"Databases: {args.db}")
        print(f"Skip BLAST?: {args.skip_blast}")

    # now lets get our dbs to run against
    dbs = get_db_type(args.db)
    print(f"Databases to run against: {dbs}")
    # ok now let's get our inputs.
    inputs = get_input_files(input_path, 'fasta')
    print(f"Input files: {len(inputs)}")
    print('\n')

    if not os.path.exists(output_path):
        os.makedirs(output_path)
        if not args.quiet:
            print('Creating output directory!\n')

    for db in dbs:
        db_path = db[0]
        prog = db[1]
        db_name = Path(db_path).stem
        results_dir = os.path.join(output_path, db_name, 'xml_files')
        tables_dir = os.path.join(output_path, db_name, 'tables')
        os.makedirs(results_dir, exist_ok=True)
        os.makedirs(tables_dir, exist_ok=True)

        if not args.skip_blast:
            # Run BLAST without parallelization since there's only one input file
            run_blast(prog, db_path, args.input, results_dir, args.threads, args.quiet)

        # Parse results
        results = get_results(results_dir, args.quiet)
        for result_file in results:
            results_table = f'{Path(result_file).stem}_table.tsv'
            output_table = str(os.path.join(tables_dir, results_table))
            parse_to_tsv(result_file, output_table, args, db_name)

    if not args.quiet:
        print("Finished!")
    return 0

if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(prog = "PlasmidCaller_v6.0.0",
                                    description = "Plasmid/replicon classifier for a single FASTA file")
    # Add the arguments
    parser.add_argument('--input', required=True, type=str, help='Input FASTA file path')
    parser.add_argument('--output', required=True, type=str, help='The directory for outputs')
    parser.add_argument('--threads', required=False, type=int, help='How many cores to use?')
    parser.add_argument('--db', required=False, type=str, help='Path to directory containing blast databases', default='/db')
    parser.add_argument('--skip_blast', action='store_true', help='Skip running BLAST and only parse existing results', default=False)
    parser.add_argument('--quiet', action='store_true', help='Run without terminal output for workflow integration', default=False)

    # Parse the arguments
    args = parser.parse_args()
    if not args.quiet:
        print(args)
    main(args)