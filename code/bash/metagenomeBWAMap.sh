#! /bin/bash
# metagenomeBWAMap.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Use BWA-MEM to map reads to contigs/bins.
# Usage: bash metagenomeBWAMap.sh FASTQTRIMMEDPAIRED1 FASTQTRIMMEDPAIRED2 FASTAINDEX

##################
# Set Script Env #
##################

# Variables defined by user
FASTQTRIMMEDPAIRED1=${1:?ERROR: Need to define FASTQTRIMMEDPAIRED1.} 
FASTQTRIMMEDPAIRED2=${2:?ERROR: Need to define FASTQTRIMMEDPAIRED2.}
FASTAINDEX=${3:?ERROR: Need to define FASTAINDEX.} # Supply only one of the index files for pulling the name

# Other variables
THREADS=$(nproc) # Automatically determines the number of cores based on local resources
OUTDIR=data/bwa/



##############################################
# Mapping Reads to Contigs/Bins with BWA MEM #
##############################################

echo PROGRESS: Aligning reads from "${SAMPLE}".

# Pulling sample name for use in naming output files
SAMPLE=$(echo "${FASTQTRIMMEDPAIRED1}" | sed 's/.*\/\(.*\)_R.*\.fq\.gz/\1/')

# Pulling index name to determine the type of sequences (contigs vs bins)
INDEXNAME=$(echo "${FASTAINDEX}" | sed 's/\(.*\)\.[a-z]\{2,3\}/\1/')

# Mapping trimmed reads to contigs
if echo "${INDEXNAME}" | grep -q "contig"; then

	# Create dir if it doesn't already exist
	mkdir -p "${OUTDIR}"/ "${OUTDIR}"/contigs/ "${OUTDIR}"/contigs/bam/

	# Mapping reads
	bwa mem -t "${THREADS}" -M "${INDEXNAME}" \
		"${FASTQTRIMMEDPAIRED1}" \
		"${FASTQTRIMMEDPAIRED2}" |
	samtools sort -o "${OUTDIR}"/contigs/bam/"${SAMPLE}"_contig_aligned.bam # Sorts by coordinate and formats the output as a .bam file for compactness

# Mapping trimmed reads to bins
elif echo "${INDEXNAME}" | grep -q "bin"; then

	# Create dir if it doesn't already exist
	mkdir -p "${OUTDIR}"/ "${OUTDIR}"/bins/ "${OUTDIR}"/bins/bam/

	# Mapping reads
	bwa mem -t "${THREADS}" -M "${INDEXNAME}" \
		"${FASTQTRIMMEDPAIRED1}" \
		"${FASTQTRIMMEDPAIRED2}" |
	samtools sort -o "${OUTDIR}"/bins/bam/"${SAMPLE}"_bin_aligned.bam # Sorts by coordinate and formats the output as a .bam file for compactness

fi
