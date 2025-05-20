# unused functions for now.

def create_blast_db(dbtype, input_file, db_output):
    """Create a blast db of specified type from input_file, and saved to db_output"""
    params = {
        "in": input_file,
        "dbtype": dbtype,
        "parse_seqids": True,
        "logfile": f"{Path(input_file).stem}_blastdb_creation.log",
        "out": db_output,
    }
    try:
        result = blast_manager.run_blast_command("makeblastdb", **params)
        return f"SUCCESS: makeblastdb for {input_file}" if result.returncode == 0 else result.stderr
    except Exception as e:
        return f"ERROR: {e}"

### This does not appear to be used.
def run_blast__old(prog, database, input_file, results_dir, threads, quiet=False):
    if not quiet:
        print(f"Running {prog} for {input_file}!")
    # get params for this command
    params = get_blast_command(prog, input_file, results_dir, database, threads)

    try:
        result = blast_manager.run_blast_command(prog, **params)
        return f"Success: {prog} for {input_file}" if result.returncode == 0 else result.stderr
    except Exception as e:
        return f"ERROR: {e}"
###
