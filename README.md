# Plasmid Caller

A bioinformatics tool for identifying and analyzing plasmids in B. burgdorferi genome assemblies.

## Installation

```bash
git clone https://github.com/mjfos2r/plasmid_caller.git
uv venv && source activate .venv/bin/activate
uv pip install plasmid_caller
```

## Usage

```bash
# always work under this venv
source activate .venv/bin/activate
```

```bash
# first usage to make sure bash is groovy
plasmid_caller --version
# to use from then on:
plasmid_caller --input <input_file> --output <output_dir> --threads <threads>
```