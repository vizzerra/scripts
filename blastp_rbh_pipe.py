# Show plots as part of the notebook and make tools available
%matplotlib inline
import matplotlib.pyplot as plt

# Standard library packages
import os

# Import Numpy, Pandas and Seaborn
import numpy as np
import pandas as pd
import seaborn as sns

# Import Biopython tools for running local BLASTX
from Bio.Blast.Applications import NcbiblastpCommandline

# Colour scale transformation
from matplotlib.colors import LogNorm

# Define input and output directories
datadir = "PATH/TO/DATA/DIR"
outdir = "PATH/TO/OUT/DIR"

# Define input file paths
s1 = os.path.join(datadir, 'FILE')
s2 = os.path.join(datadir, 'FILE')

# Define output BLAST results
fwd_out = os.path.join(outdir, 'FWD_OUTFILE')
rev_out = os.path.join(outdir, 'REV_OUTFILE')

# Create BLAST command-lines for forward and reverse BLAST searches
# Change outfmt to suit your needs
fwd_blastp = NcbiblastpCommandline(query=s1, subject=s2, out=fwd_out,
                                   outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue",
                                   max_target_seqs=1)
rev_blastp = NcbiblastpCommandline(query=s2, subject=s1, out=rev_out,
                                   outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue",
                                   max_target_seqs=1)

# Inspect command-lines
print("FORWARD: %s" % fwd_blastp)
print("REVERSE: %s" % rev_blastp)

# THIS CELL RUNS LARGE LOCAL BLAST SEARCHES
# IT IS SKIPPED BY DEFAULT

# Run BLAST searches
# !! Uncomment to run local BLAST searches !!
#fwd_stdout, fwd_stderr = fwd_blastp()
rev_stdout, rev_stderr = rev_blastp()

# Check STDOUT, STDERR
#print("FWD STDOUT: %s" % fwd_stdout)
#print("FWD STDERR: %s" % fwd_stderr)
#print("REV STDOUT: %s" % rev_stdout)
#print("REV STDERR: %s" % rev_stderr)

# PRECALCULATED BLAST RESULTS
# COMMENT OUT THESE TWO LINES IF YOU WANT TO USE THE RESULTS FROM THE CELL ABOVE
fwd_out = os.path.join('path','to', 'file')
rev_out = os.path.join('path','to', 'file')

# Load the BLAST results into Pandas dataframes
fwd_results = pd.read_csv(fwd_out, sep="\t", header=None)
rev_results = pd.read_csv(rev_out, sep="\t", header=None)

# Add headers to forward and reverse results dataframes
headers = ["query", "subject", "identity", "coverage",
           "qlength", "slength", "alength",
           "bitscore", "E-value"]
fwd_results.columns = headers
rev_results.columns = headers