#! /bin/bash
# metagenomeQuast.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Calculate contig statistics and characteristics following assembly.
# Usage: bash metagenomeQuast.sh CONTIGS1 CONTIGS2 CONTIGS3 ... CONTIGSN

##################
# Set Script Env #
##################

# Variables defined by user
CONTIGS=${@:?ERROR: Need to define CONTIGS.}

# Other variables
OUTDIR=data/quast/



###############################################
# Calculating Assembly Statistics Using Quast #
###############################################

echo PROGRESS: Running Quast to calculate assembly statistics.

# Running Quast on all of the assembled contigs.fasta files at once
quast.py --space-efficient -L $(echo "${CONTIGS[@]}") -o "${OUTDIR}"

# Use transposed_report.tsv for analyzing in R if desired
