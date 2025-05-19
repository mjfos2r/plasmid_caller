if RUN_BLAST == True:
    for input_file in Assemblies:
        path = Path(input_file)
        suffix = path.suffixes[0]
        assembly_id = path.name.split(suffix)[0]
        blast_cmd = "blastn -query "+input_file+" -task 'blastn' "
        blast_cmd += "-db "+blast_db
        blast_cmd += " -out "+output_file+".xml"
        blast_cmd += " -evalue 1e-100 -num_threads 8 -outfmt 5 -max_target_seqs 5 -max_hsps 1"
        subprocess.run(blast_cmd, shell=True)