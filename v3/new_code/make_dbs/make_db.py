if RUN_DB_BOOL is True:
    make_db_cmd = "makeblastdb -in "+fasta_for_db+" -out "+blast_db+" -dbtype 'nucl' -hash_index"
    print(make_db_cmd)
    subprocess.run(make_db_cmd, shell=True)

else:
    print("Blast DB already exists! Not creating it again!")
    print("Change RUN_DB_BOOL to True if you need to rebuild the DB.")