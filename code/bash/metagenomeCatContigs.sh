#! /bin/bash
# metagenomeCatContigs.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Combine all of the assembled contigs into a single contig library.
# Usage: bash metagenomeCatContigs.sh CONTIGS1 CONTIGS2 CONTIGS3 ... CONTIGSN

##################
# Set Script Env #
##################

# Variables defined by user
CONTIGS=${@:?ERROR: Need to define CONTIGS.}

# Other variables
OUTDIR=data/metaspades/library/



###################################
# Concatenating Assembled Contigs #
###################################

echo PROGRESS: Concatenating all of the assembled contigs into a single fasta.

# Making the output dir if it doesn't already exist
mkdir -p "${OUTDIR}"

# Combining all of the assembled contigs into a single output fasta
cat $(echo "${CONTIGS[@]}") > "${OUTDIR}"/contig_library.fa
