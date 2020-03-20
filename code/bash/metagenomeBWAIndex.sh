#! /bin/bash
# metagenomeBWAIndex.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Use BWA index to create indices of contigs/bins before mapping reads to them.
# Usage: bash metagenomeBWAIndex.sh FASTA

##################
# Set Script Env #
##################

# Variables defined by user
FASTA=${1:?ERROR: Need to define FASTA.}

# Other variables
THREADS=$(nproc) # Automatically determines the number of cores based on local resources
MEM=$(echo "${THREADS}" | awk '{print int($1*4*0.98)}') # Determines total available mem by assuming 4000mb per thread/processor and 2% RAM usage by system processes
BLOCKSIZE=$(echo "${MEM}" | awk '{print int($1/8*10^9)}') # Total memory usage = 8 * $BLOCKSIZE (number of bases that can be processed at once) so defining based on that
OUTDIR=data/bwa/



####################################################
# Creating BWA Alignment Indices for Mapping Reads #
####################################################

echo PROGRESS: Creating BWA alignment indices.

# Pulling FASTA file name for naming indexed output
FILENAME=$(echo "${FASTA}" | sed 's/.*\/\(.*\)\.fa.*/\1/')

if echo "${FILENAME}" | grep -q "contig"; then

	# Create dir if it doesn't already exist
	mkdir -p "${OUTDIR}"/ "${OUTDIR}"/contigs/ "${OUTDIR}"/contigs/index/

	# Creating BWA index using bwtsw algorithm
	bwa index -a bwtsw -b "${BLOCKSIZE}" -p "${OUTDIR}"/contigs/index/"${FILENAME}" "${FASTA}"

elif echo "${FILENAME}" | grep -q "bin"; then

	# Create dir if it doesn't already exist
	mkdir -p "${OUTDIR}"/ "${OUTDIR}"/bins/ "${OUTDIR}"/bins/index/

	# Creating BWA index using bwtsw algorithm
	bwa index -a bwtsw -b "${BLOCKSIZE}" -p "${OUTDIR}"/bins/index/"${FILENAME}" "${FASTA}"

fi
