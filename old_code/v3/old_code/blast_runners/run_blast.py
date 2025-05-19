# # Blast DB creation and search execution 
# Okay so the following requires setting up a jupyter server in the x86 conda environment
# 
# to do this, open a new terminal and type the following:
# ```$bash
# >condax86
# >conda activate blast
# >cd /Users/mf019/bioinformatics/blastDB/BbPlasmid
# >jupyter notebook
# ```
# then grab the server URL and put that into the python kernel option in VSCode
# http://localhost:8888/tree?token=b783c19720292ce981bed187f7550a10bfebf745d7e14cb8
# 
# click the python3 kernel button in the top right
# -> select another kernel -> existing jupyter server 
# -> enter URL
# 

import subprocess
import os 
import Bio
import glob
import re
import tqdm
from time import sleep
from pathlib import Path

RUN_DB_BOOL = False
working_directory = '/Users/mf019/bioinformatics/blastDB/BbPlasmid'
fasta_for_db = 'db/db_attempt_2/all_plasmids.fasta'
blast_db = 'db/db_attempt_2/db/plasmid_db'
assembly_dir = 'assemblies/'
canu_blast_out = 'plasmid_seqs/fresh_attempt/output/blast_results/canu/'
flye_blast_out = 'plasmid_seqs/fresh_attempt/output/blast_results/flye/'

print("CWD: ",os.getcwd())
print("moving to working directory!")
os.chdir(working_directory)
print("CWD: ",os.getcwd())

# lets glob the assemblies from their respective directories!
canu_assemblies = glob.glob(assembly_dir+'canu/*.fasta')
flye_assemblies = glob.glob(assembly_dir+'flye/*.fasta')

canu_pbar = tqdm.tqdm(range(len(canu_assemblies)), desc="Canu Assemblies Found: ")
for file in canu_assemblies:
    sleep(0.1)
    canu_pbar.update(1)
canu_pbar.close()
    
flye_pbar = tqdm.tqdm(range(len(flye_assemblies)), desc="Flye assemblies found: ")
for file in flye_assemblies:
    sleep(0.1)
    flye_pbar.update(1)
flye_pbar.close()

if RUN_DB_BOOL is True:
    make_db_cmd = "makeblastdb -in "+fasta_for_db+" -out "+blast_db+" -dbtype 'nucl' -hash_index"
    print(make_db_cmd)
    subprocess.run(make_db_cmd, shell=True)

else:
    print("Blast DB already exists! Not creating it again!")
    print("Change RUN_DB_BOOL to True if you need to rebuild the DB.")

#loop over the files in the canu dir and then set up the blast command for each file

# Run Blast on all Canu Assemblies!
num_files = len(canu_assemblies)
canu_pbar = tqdm.tqdm(total=num_files, desc="Canu Assemblies Blasted: ")
for input_file in canu_assemblies:
    path = Path(input_file)
    suffix = path.suffixes[0]
    assembly_id = path.name.split(suffix)[0]
    output_file = canu_blast_out+assembly_id+"_canu_blast"
    blast_cmd = "blastn -query "+input_file+" -task 'blastn' "
    blast_cmd += "-db "+blast_db
    blast_cmd += " -out "+output_file+".xml"
    blast_cmd += " -evalue 1e-100 -num_threads 8 -outfmt 5 -max_target_seqs 5 -max_hsps 1"
    subprocess.run(blast_cmd, shell=True)
    canu_pbar.update(1)
canu_pbar.close()

# run Blast on all flye assemblies
num_files = len(flye_assemblies)
flye_pbar = tqdm.tqdm(total=num_files, desc="Flye Assemblies Blasted: ")
for input_file in flye_assemblies:
    path = Path(input_file)
    suffix = path.suffixes[0]
    assembly_id = path.name.split(suffix)[0]
    output_file = flye_blast_out+assembly_id+"_flye_blast"
    blast_cmd = "blastn -query "+input_file+" -task 'blastn' "
    blast_cmd += "-db "+blast_db
    blast_cmd += " -out "+output_file+".xml"
    blast_cmd += " -evalue 1e-100 -num_threads 8 -outfmt 5 -max_target_seqs 5 -max_hsps 1"
    subprocess.run(blast_cmd, shell=True)
    flye_pbar.update(1)
flye_pbar.close()





