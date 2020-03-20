#! /bin/bash
# metagenomeDedupe.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Use BBMap dedupe.sh to remove exact replicates and containments from a contig library.
# Usage: bash metagenomeDedupe.sh FASTA

##################
# Set Script Env #
##################

# Variables defined by user
FASTA=${1:?ERROR: Need to define FASTA.}

# Other variables
OUTDIR=data/dedupe/



#############################################################
# Removing Duplicates and Containments with BBMap dedupe.sh #
#############################################################

echo PROGRESS: Running BBMap dedupe.sh to remove sequence duplicates.

# Create dir if it doesn't already exist
mkdir -p "${OUTDIR}"/

# Pulling FASTA file name for naming deduped output
FILENAME=$(echo "${FASTA}" | sed 's/.*\/\(.*\)\.fa.*/\1/')

# Deduplicating the sequences
dedupe.sh in="${FASTA}" out="${OUTDIR}"/"${FILENAME}"_deduped.fa
