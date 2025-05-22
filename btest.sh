#!/bin/bash

source .venv/bin/activate

plasmid_caller -i tests/test_exec/test_B100609.fasta -o tests/test_exec/output -t $(nproc)
exit 0
