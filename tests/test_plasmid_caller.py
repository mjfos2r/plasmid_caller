# tests/test_plasmid_caller.py

import io
import json
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import pandas
import pytest

import plasmid_caller as pc

@pytest.fixture
def tmp_fasta(tmp_path):
    #TBD
    pass