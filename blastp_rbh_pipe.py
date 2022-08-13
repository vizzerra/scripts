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

# Color scale transformation
from matplotlib.colors import LogNorm

# Define input and output directories
# You may need to change this depedning on OS
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
# You can comment out below if you have previous BLAST results, see line 59
fwd_stdout, fwd_stderr = fwd_blastp()
rev_stdout, rev_stderr = rev_blastp()

# Check STDOUT, STDERR
print("FWD STDOUT: %s" % fwd_stdout)
print("FWD STDERR: %s" % fwd_stderr)
print("REV STDOUT: %s" % rev_stdout)
print("REV STDERR: %s" % rev_stderr)

# FOR PRECALCULATED BLAST RESULTS
# COMMENT OUT THESE TWO LINES IF YOU WANT TO USE THE RESULTS FROM THE CELL ABOVE
# fwd_out = os.path.join('path','to', 'file')
# rev_out = os.path.join('path','to', 'file')

# Load the BLAST results into Pandas dataframes
fwd_results = pd.read_csv(fwd_out, sep="\t", header=None)
rev_results = pd.read_csv(rev_out, sep="\t", header=None)

# Add headers to forward and reverse results dataframes
# You can change if your blast.tab are different
headers = ["query", "subject", "identity", "coverage",
           "qlength", "slength", "alength",
           "bitscore", "E-value"]
fwd_results.columns = headers
rev_results.columns = headers

# FIND RECIPROCAL BEST HITS (RBH)

# Merge foreward and reverse resutls
rbbh = pd.merge(fwd_results, rev_results[['query', 'subject']],
                left_on='subject', right_on='query',
                how='outer')

# Remove non-RBH
rbbh = rbbh.loc[rbbh.query_x == rbbh.subject_y]

# Group duplicate RBH rows, taking the maximum value in each column
rbbh = rbbh.groupby(['query_x', 'subject_x']).max()

print("Done.") 
