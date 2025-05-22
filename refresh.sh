#!/bin/bash

source .venv/bin/activate
uv pip uninstall plasmid_caller
git pull
uv pip install -e .
exit 0
