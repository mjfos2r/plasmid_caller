# plasmid_caller.py
# Version: 6.0.0
# 2025-May-16
# - Michael J. Foster
# - Ben Kotzen
# Rewritten to incorporate best practices and simplify execution logic.
# To be added:
## operate on fastq reads.
## generate a renaming table that can be fed into Bakta.
## rename the input fasta file/genbank using the calls determined.
## output an evidence table.

import os
import re
import subprocess
import glob
import pandas
import pickle
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from tqdm import tqdm
from pathlib import Path
from intervaltree import Interval, IntervalTree
from concurrent.futures import ProcessPoolExecutor, as_completed