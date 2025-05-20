# unused functions for now.

def get_output_path(results_dir):
    cwd = Path(os.getcwd())
    output_path = os.path.join(cwd, results_dir)
    if not os.path.exists(output_path):
        print(f"Creating output path: {output_path}")
        os.makedirs(output_path, exist_ok=True)
    return output_path

def get_input_files(input_path, input_extension):
    files = []
    files.extend(Path(input_path).glob(f"*.{input_extension}"))
    return files

