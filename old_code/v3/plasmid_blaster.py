import os
import glob
import pickle
import subprocess
from axe import axe
from Bio import SeqIO
from pathlib import Path
from collections import defaultdict

def parse_blast_results(blast_results : str) -> dict:
    file_size = os.path.getsize(blast_results)
    if file_size == 0:
        return None
    results = {}
    with open(blast_results, "r") as f:
        blast_records = NCBIXML.parse(f)
        for record in blast_records:
            results[record.query] = []
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    results[record.query].append((alignment.title, hsp.expect, hsp.score, hsp.query_start, hsp.query_end, hsp.sbjct_start, hsp.sbjct_end, hsp.identities, hsp.positives, hsp.gaps, hsp.frame, hsp.query, hsp.sbjct))
    return results

mjf = axe.mjf(logfile_name='plasmid_identification_v3_blastlog', logfile_level='DEBUG', console_level='DEBUG')

testing_set_file = 'dbs/Minimal-PF32-testor-set-March2021.faa'

################################################################################
#:::::::::::::::::::::::::::::::::define paths:::::::::::::::::::::::::::::::::#
################################################################################
# get working directory
init_cwd = os.getcwd()
# Set working directory
working_dir = '/Users/mf019/new_plasmid_id'
# path out the databases
databases_dir = f'{working_dir}/dbs'
whole_plasmid_db = f'{databases_dir}/plasmid_db/plasmid_db'
pf32_db = f'{databases_dir}/pf32_db/pf32'
# lets specify where parsing tables are to be found!
parsing_tables_dir = f'{working_dir}/parsing_tables'
# set assemblies directory # DO I WANT TO ARGV THIS?
assemblies_dir = '/Users/mf019/bioinformatics/longread_GWAS/assemblies/paired_assemblies/paired_only'
shortread_dir  = f'{assemblies_dir}/shortread'
longread_dir = f'{assemblies_dir}/longread'
shortread_contigs  = f'{shortread_dir}/contigs'
longread_contigs = f'{longread_dir}/contigs'
shortread_annotations  = f'{shortread_dir}/annotation'
longread_annotations = f'{longread_dir}/annotation'
output_dir = f'{working_dir}/output'
blast_results_dir = f'{output_dir}/blast_results'
pf32_blast_results = f'{blast_results_dir}/pf32'
plasmid_blast_results = f'{blast_results_dir}/whole_plasmid'
################################################################################
#::::::::::::::::::::::::::::::log paths to file:::::::::::::::::::::::::::::::#
################################################################################
mjf.logblock(msg='analysis inputs')
mjf.loginfo(f'Working directory: {working_dir}')
mjf.loginfo(f'Databases directory: {databases_dir}')
mjf.loginfo(f'Plasmid database: {whole_plasmid_db}')
mjf.loginfo(f'PF32 database: {pf32_db}')
mjf.loginfo(f'Parsing tables directory: {parsing_tables_dir}')
mjf.loginfo(f'Assemblies directory: {assemblies_dir}')
mjf.loginfo(f'Shortread directory: {shortread_dir}')
mjf.loginfo(f'Longread directory: {longread_dir}')
mjf.loginfo(f'Shortread contigs directory: {shortread_contigs}')
mjf.loginfo(f'Shortread annotations directory: {shortread_annotations}')
mjf.loginfo(f'Longread contigs directory: {longread_contigs}')
mjf.loginfo(f'Longread annotations directory: {longread_annotations}')
mjf.loginfo(f'Output directory: {output_dir}')
mjf.loginfo(f'Whole contig blast results directory: {blast_results_dir}')
mjf.loginfo(f'PF32 Blast results: {pf32_blast_results}')
mjf.logblock(msg='Find Assemblies in Input')
################################################################################
#::::::::::::::::::::::::::::::get genbank files:::::::::::::::::::::::::::::::#
################################################################################
# get annotations for each method
sr_gbs = glob.glob(f'{shortread_annotations}/*.gbff')
mjf.loginfo(f'Found {len(sr_gbs)} shortread gbff files')
lr_gbs = glob.glob(f'{longread_annotations}/*.gbff')
mjf.loginfo(f'Found {len(lr_gbs)} longread gbff files')
mjf.logblock(msg='Dictionary Assembly Begin')
assembly_dict = defaultdict(dict)
for assembly in lr_gbs: # iterate through each longread assembly
    assembly_name = Path(assembly).stem # get the name of the assembly
    mjf.logdebug(f'Processing assembly: {assembly_name}')
    assembly_dict[assembly_name] = defaultdict(dict)
    assembly_dict[assembly_name]['longread'] = defaultdict(dict)
    with open(assembly, 'r') as f: # get longread contigs and annotations
        num_contigs = 0
        num_genes = 0
        for contig in SeqIO.parse(handle=f, format='genbank'):
            num_contigs += 1
            contig_id = contig.id
            assembly_dict[assembly_name]['longread'][contig_id]['seq'] = contig.seq
            assembly_dict[assembly_name]['longread'][contig_id]['annotations'] = contig.annotations
            assembly_dict[assembly_name]['longread'][contig_id]['features'] = contig.features
            if len(contig.features) > 0:
                for feature in contig.features:
                    if feature.type == 'CDS':
                        num_genes += 1
        mjf.loginfo(f'longread:{assembly_name} has {num_contigs} contigs')
        mjf.loginfo(f'longread: {assembly_name} has {num_genes} annotated genes')
    # now we need to get the shortread contigs and annotations if they exist.
    shortread_name = assembly_name[:-1:]
    shortread_file = f'{shortread_annotations}/{shortread_name}.gbff'
    if shortread_file in sr_gbs:
        mjf.logdebug(f'Found shortread annotations for {shortread_name}')
        assembly_dict[assembly_name]['shortread'] = defaultdict(dict)
        with open(shortread_file, 'r') as f:
            num_genes = 0
            num_contigs = 0
            for contig in SeqIO.parse(handle=f, format='genbank'):
                num_contigs += 1
                contig_id = contig.id
                assembly_dict[assembly_name]['shortread'][contig_id]['seq'] = contig.seq
                assembly_dict[assembly_name]['shortread'][contig_id]['annotations'] = contig.annotations
                assembly_dict[assembly_name]['shortread'][contig_id]['features'] = contig.features
                if len(contig.features) > 0:
                    for feature in contig.features:
                        if feature.type == 'CDS':
                            num_genes += 1
                # now we need to blast this contig against the pf32 and whole plasmid databases
                # we will need to save these results in a dictionary
        mjf.loginfo(f'shortread:{shortread_name} has {num_contigs} contigs')
        mjf.loginfo(f'shortread: {shortread_name} has {num_genes} annotated genes')
    else:
        assembly_dict[assembly_name]['shortread'] = 'NOT AVAILABLE'
        mjf.logerror(f'No shortread annotations found for {assembly_name}')
mjf.logwarn('dictionary construction complete! pickling!')
pickle.dump(assembly_dict, open(f'{output_dir}/assembly_dict.pkl', 'wb'))
mjf.logfileop(f'dictionary pickled to: {output_dir}/assembly_dict.pkl')
################################################################################
#::::::::::::::::::::::::::::Begin Blast Execution:::::::::::::::::::::::::::::#
################################################################################
mjf.logblock(msg=str('Start Blastin'))
for isolate in assembly_dict:
    input_files = [] # empty list to store input files
    lr_data = assembly_dict[isolate]['longread']
    sr_data = assembly_dict[isolate]['shortread']
    lr_id = isolate # since we're operating on lr first
    input_files.append(f'{longread_contigs}/{lr_id}.fasta') # add lr contigs to list
    if sr_data == 'NOT AVAILABLE':
        mjf.logcritical(f'No shortread annotations found for {isolate}')
        sr_id = None
    else:
        sr_id = isolate[:-1:]
        input_files.append(f'{shortread_contigs}/{sr_id}.fasta')
    mjf.loginfo(f'setting up blast commands for {isolate}')
    for file in input_files:
        name = Path(file).stem
        wp_output_file = f'{plasmid_blast_results}/{name}_whole_plasmid.xml'
        wp_blast_cmd = f'blastn -query {file} -task "blastn" '
        wp_blast_cmd += f'-db {whole_plasmid_db} '
        wp_blast_cmd += f' -out {wp_output_file} '
        wp_blast_cmd += "-evalue 1e-100 -num_threads 8 -outfmt 5 -max_target_seqs 5 -max_hsps 5"
        mjf.logwarn(f'Running blast command:\n{wp_blast_cmd}')
        subprocess.run(wp_blast_cmd, shell=True) # run it
        mjf.logfileop(f'whole plasmid blast results written to: {wp_output_file}')
        mjf.loginfo(f'setting up blast command for pf32 blast')
        pf_output_file = f'{pf32_blast_results}/{name}_pf32.xml'
        pf_blast_cmd = f'blastx -query {file} -task "blastx" '
        pf_blast_cmd += f'-db {pf32_db} '
        pf_blast_cmd += f' -out {pf_output_file} '
        pf_blast_cmd += "-evalue 1e-100 -num_threads 8 -outfmt 5 -max_target_seqs 5 -max_hsps 5"
        subprocess.run(pf_blast_cmd, shell=True, check=False) # run it, but don't check for errors to avoid stopping the script
        mjf.logfileop(f'pf32 blast results written to: {pf_output_file}')
mjf.logblock(msg='saddle up, we are done here!')
mjf.choplog()