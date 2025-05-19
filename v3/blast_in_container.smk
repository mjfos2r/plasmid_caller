# Snakefile
import glob
from pathlib import Path

Samples = [Path(asm).stem for asm in glob.glob('assemblies/*.fasta')]

rule all:
    input:
        expand('output/blast_results/whole_plasmid/{asm}_whole_plasmid.xml', asm=Samples),
        expand('output/blast_results/pf32/{asm}_pf32.xml', asm=Samples)

rule blast:
    input:
        'assemblies/{asm}.fasta'
    output:
        wp_out = 'output/blast_results/whole_plasmid/{asm}_whole_plasmid.xml',
        pf32_out = 'output/blast_results/pf32/{asm}_pf32.xml'
    params:
        pf32_db = 'dbs/pf32_db/pfam32db',
        wp_db = 'dbs/plasmid_db/wpdb'
    threads: 60
    shell:
        'blastx -query {input} -task "blastx" -db={params.pf32_db}  -out {output.pf32_out} -evalue 1e-100 -num_threads {threads} -outfmt 5 -max_target_seqs 5 -max_hsps 5 &&'
        'echo "blastn finished!" &&'
        'blastn -query {input} -task "blastn" -db={params.wp_db}  -out {output.wp_out} -evalue 1e-100 -num_threads {threads} -outfmt 5 -max_target_seqs 5 -max_hsps 5 &&'
        'echo "blastx finished!"'