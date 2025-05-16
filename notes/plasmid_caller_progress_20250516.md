# Plasmid Caller

Michael J. Foster

2025-May-16:

I've been reviewing the code: re-familiarizing myself with current functionality and identifying which gaps still exist.
I'm currently re-writing the most recent version to greatly simplify execution and processing logic. There is quite a bit of dead weight from prior versions and I would like to have things run in a more modular fashion to accommodate future methods (HMM, synteny/gene presence, etc.) Inclusion of additional classification methods are not immediately planned but the current rigidity needs fixing regardless.

## Current Progress

- Organized repository
- Reviewing and optimizing the core processing logic.
    >It was originally written for bulk processing of many genomes (using multiprocessing/multithreading) and given the switch to single genome processing, much of it is dead weight. I am extending current functionality to perform single input processing by default and multiple genome processing when specified.
    >The switch to single-input processing was done to simplify parallel execution in pipelines but multi-input execution is still useful, especially when not in a pipeline.
- Reviewing and optimizing output generation (matrix, evidence table, etc.)
- Extending functionality to accomodate single/multi input
    >the shift from multi-input processing to single-input processing broke *many* things and while they are functional, the fixes are far from flexible.
- Implementing hatchling as the primary build backend
    >(To allow for `pip install plasmid_caller` using `pyproject.toml`)
- Implementing PyTest as the testing framework
- BLAST functionality is being modularized.
- implementing self-awareness for docker execution

## Additions to Functionality

### Immediate

- [ ] automate container builds / tag management ( already have the makefiles written, just need to test )
- [ ] creation of `pyproject.toml`
- [ ] implementation of build backend
- [ ] implementation of testing framework
- [ ] unit testing for all accepted inputs.
    >FASTA, genbank, ...
- [ ] input validation/input sanitization
    >(especially important if we want to deploy this as a public web service)
- [ ] automate core logic testing and generate test dataset. (B31 reference for positive tests, error datasets to check error handling is in progress)
- [ ] checks to determine if execution is happening within a docker container or native installation and modify paths accordingly.
    >this is to simplify things for the end user, mainly path composition to output results to a mounted volume. (else they're ephemeral and die with the container)

### Upcoming

- [ ] write `.wdl` for use in terra.
- [ ] documentation/wiki creation
- [ ] return input files with renamed headers ( very high priority )
- [ ] generate renaming_table for use in Bakta annotation
- [ ] pickle raw hits, pickle temp files, emit these after execution for debugging (when debug flag is set)
- [ ] dependency checking on installation
    >(python deps should be handled by pip/hatchling but I need to accomodate this for BLAST installation)
    >The docker container comes with everything pre-installed and ready to use but this still needs to be handled to simplify that.

### Planned

- [ ] Deploy as public facing web app for anyone to use. (this is a primary motivation for input sanitization and rigorous testing)
- [ ] read classification + binning for individual replicon assembly.
    >>(fastq processing is not functional, core logic is too rigid and blasting every read is far from optimal. This will likely require an HMM approach and/or DMND, I have several ideas but this is low priority)
- [ ] gene presence & order as another source of evidence. (I have some draft attempts at implementing this but this is also low-priority)

## Current Dependencies

[Python]
- Python 3.12

[Modules]
- biopython
- pandas
- numpy
- tqdm
- intervaltree

[Build Backend]
- Hatch/Hatchling

[Testing Backend]
- PyTest

[External Tools]
- ncbi-blast-2.16.0+

[Base Docker image]
- python:3.12-slim
