# plasmid_caller.py
# Version: 5.1.2
# December 5, 2024
# - Michael J. Foster
# - Ben Kotzen

import os
import csv
import sys
import subprocess
import glob
import pandas
import pickle
import argparse
from math import floor
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from tqdm import tqdm
from pathlib import Path
from intervaltree import Interval, IntervalTree
from concurrent.futures import ProcessPoolExecutor, as_completed

def get_input_files(in_dir, ext):
    print(f'looking in {in_dir} for input files with extension: {ext}')
    files = glob.glob(f'{in_dir}/*.{ext}')
    if len(files) == 0:
        files = glob.glob(f'{in_dir}/**/*.{ext}', recursive=True)
    print(f'Input files found: {len(files)}')
    return files

def get_output_path(results_dir):
    cwd = Path(os.getcwd())
    output_path = os.path.join(cwd, results_dir)
    if not os.path.exists(output_path):
        print(f"Creating output path: {output_path}")
        os.makedirs(output_path, exist_ok=True) # too lazy to refactor this function rn i just need my calls 
    return output_path

def get_blast_command(prog, input_file, output_path, database, job_threads):
    file_id = Path(input_file).stem
    output_file = f'{output_path}/{file_id}_blast_results.xml'
    command = [prog, '-query', input_file, '-db', database, 
               '-out', output_file, '-evalue', '1e-100', '-num_threads', job_threads, 
               '-outfmt', '5', '-max_target_seqs', '5', '-max_hsps', '10'] # '-task', f'"{prog}"',
    return command

def make_blast_db(dbtype, input_file, db_output):
    logfile = f'{Path(input_file).stem}_blastdb_creation.log'
    command = ['makeblastdb','-in', input_file, '-dbtype', dbtype,
               '-parse_seqids', '-logfile', logfile, '-out', db_output]
    return command

def get_db_type(db_dir):
    """ Get the individual database for blast and the database type """
    command = ['blastdbcmd', '-list', f'{db_dir}', '-recursive',]
    results = subprocess.run(command, capture_output=True, text=True)
    dbs = []
    progs = []
    output = results.stdout
    for line in output.split('\n')[:-1]: # last line is empty so don't iterate to it!
        db_info = line.split(' ')
        print(db_info)
        db_name = db_info[0]
        db_type = db_info[1].lower()
        dbs.append(db_name)
        if db_type == 'protein':
            progs.append('blastx')
        elif db_type == 'nucleotide':
            progs.append('blastn')
    return list(zip(dbs, progs))

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
        plasmid_id = gene_id_list[-1]
    elif len(gene_id_list) == 1: # this is the weird pdb case.
        strain = gene_id_list[0].split('|')[1]
        plasmid_id = gene_id_list[0].split('|')[-1]
    elif gene_id_list[0] == 'NE': # this is the NE_1234 strains single case, hate this too.
        strain = '_'.join(gene_id_list[0:2])
        plasmid_id = gene_id_list[2]
    elif gene_id_list.startswith('RS'):
        strain = gene_id_list[1]
        plasmid_id = gene_id_list[-3]        
    else:
        strain = gene_id_list[0]
        plasmid_id = gene_id_list[1]
        
    return strain, plasmid_id

def parse_blast_xml(xml_file, **kwargs):
    """ This is the new and improved blast results parsing function.
    - it takes kwargs for parsing_type which is either 'wp' or 'pf32'"""
    assembly_id = Path(xml_file).stem.split('_')[0]
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
            hsp_counter = 0

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
                    strain = 'unknown' # this needs work
                    plasmid_name = alignment.hit_id

                results = {
                    "assembly_id": assembly_id,
                    "contig_id": contig_id,
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
    for contig, hits in hits.items():
        for hit in hits:
            rows.append(hit)
    df = pandas.DataFrame(rows)
    return df

def join_and_write_tables(tables):
    main_df = pandas.concat(tables, ignore_index = True)
    main_df.to_csv(output_table_path, sep='\t')

def parse_to_tsv(file, output_path, parsing_type):
    hits = parse_blast_xml(file, parsing_type=parsing_type)
    hits_df = get_hits_table(hits)
    hits_df.to_csv(str(output_path), sep='\t', index=False)
    return f"{parsing_type} table for: {Path(file).stem.split('_')[0]} written to {output_path}"
    
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

def progress_bar(futures):
    with tqdm(total=len(futures)) as pbar:
        for future in as_completed(futures):
            try:
                result = future.result()
                tqdm.write(result)
            except Exception as e:
                tqdm.write(f"Error: {e}")
            pbar.update(1)

def parallel_blast(num_workers, job_threads, prog, database, input_files, results_dir):
    print('Starting parallel_blast!')
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        #output_path = get_output_path(results_dir) BAD?
        for file in input_files:
            command = get_blast_command(prog, file, results_dir, database, job_threads)
            print(f'Queueing {prog} job for {file}!')
            futures.append(executor.submit(run_command, command))
        progress_bar(futures)

def parallel_parse(cpus, results, db_name, output_dir):
    print('Starting parallel_parse!')
    with ProcessPoolExecutor(max_workers=cpus) as executor:
        futures = []
        #output_path = get_output_path(output_dir)
        for file in results:
            results_table  = f'{Path(file).stem}_table.tsv'
            output_table = str(os.path.join(output_dir, results_table))
            print(f'Queueing parse job for {file}!')
            futures.append(executor.submit(parse_to_tsv, file, output_table, db_name))
        progress_bar(futures)

def get_results(results_dir):
    print(os.listdir(results_dir))
    return glob.glob(f'{results_dir}/**/*.xml', recursive=True)


## Define main function logic.
def main(args):
    # First things first, let's fix the input to path correctly within the container. 
    # User should not have to worry about anything other than binding the working directory and specifying locations.
    input_path = os.path.join('/data', args.input)
    output_path = os.path.join('/data', args.output)
    
    if not os.path.exists(input_path):
        parser.error("Must specify a directory containing fasta files!")
    
    # determine how many parallel workers we're gonna spin up
    num_workers = floor(args.cpus/args.job_threads)
    
    print(f"Input directory: {args.input}")
    print(f"Output directory: {args.output}")
    print(f"CPUs: {args.cpus}")
    print(f"Job threads: {args.job_threads}")
    print(f"Databases: {args.db}")
    print(f"Run BLAST?: {args.blast}")
    print(f"Number of parallel workers: {num_workers}")
    
    # now lets get our dbs to run against
    dbs = get_db_type(args.db)
    print(f"Databases to run against: {dbs}")
    # ok now let's get our inputs.
    inputs = get_input_files(input_path, 'fna')
    print(f"Input files: {len(inputs)}")
    print('\n')
    
    if not os.path.exists(output_path):
        os.mkdir(output_path)
        print('creating output directory!\n')
    
    for db in dbs:
        db_path = db[0]
        prog = db[1]
        db_name = Path(db_path).stem
        results_dir = os.path.join(output_path, db_name,'xml_files')
        tables_dir = os.path.join(output_path, db_name,'tables')
        os.makedirs(results_dir, exist_ok=True)
        os.makedirs(tables_dir, exist_ok=True)
        if args.blast:
            parallel_blast(num_workers, args.job_threads, prog, db_path, inputs, results_dir)
        results = get_results(results_dir)
        parallel_parse(args.cpus, results, db_name, tables_dir)
    print("Finished!")
    return 0

if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(prog = "PlasmidCaller_v5.1.2",
                                     description = "Michael's handy dandy plasmid/replicon classifier")
    # Add the arguments
    parser.add_argument('--input',         required=True,  type=str, help='A directory containing the fasta files for the assemblies to parse')
    #parser.add_argument('annotations_dir',required=True,  type=str, help='The directory containing the annotations (genbanks)') # Not yet!
    parser.add_argument('--output',        required=True,  type=str, help='The directory for outputs')
    parser.add_argument('--cpus',          required=True,  type=int, help='How many cores we rippin')
    parser.add_argument('--job_threads',   required=False,  type=int, help='How many cores per job?', default=2)
    parser.add_argument('--db',            required=False, type=str, help='path to directory containing *all* blast databases to use in analysis', default='/db',)
    parser.add_argument('--blast', action='store_true', required=False, help='Use this to run blast against all provided databases', default=True)
    # Parse the arguments
    args = parser.parse_args()
    print(args)
    main(args)