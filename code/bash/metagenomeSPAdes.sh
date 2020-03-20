#! /bin/bash
# metagenomeSPAdes.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Aligns reads and forms contigs using metaSPAdes implementation of SPAdes.
# Usage: bash metagenomeSPAdes.sh FASTQDECONPAIRED1 FASTQDECONPAIRED2 FASTQDECONUNPAIRED

##################
# Set Script Env #
##################

# Variables defined by user
FASTQDECONPAIRED1=${1:?ERROR: Need to define FASTQDECONPAIRED1.} 
FASTQDECONPAIRED2=${2:?ERROR: Need to define FASTQDECONPAIRED2.} 
FASTQDECONUNPAIRED=${3:?ERROR: Need to define FASTQDECONUNPAIRED2.} 

# Other variables
THREADS=$(nproc) # Automatically determines the number of cores based on local resources
MEM=$(echo "${THREADS}" | awk '{print int($1*4*0.98)}') # Limits max mem by assuming 4000mb per thread/processor and 2% RAM usage by system processes
OUTDIR=data/metaspades/



#############################################
# Assembling Trimmed Reads Using MetaSPAdes #
#############################################

echo PROGRESS: Assembling trimmed reads using metaSPAdes.

# Pulling sample name for use in naming output dir
SAMPLE=$(echo "${FASTQDECONPAIRED1}" | sed 's:.*\/\(.*\)_R.*\.fq\.gz:\1:')

# Making output directory if it doesn't already exist
mkdir -p "${OUTDIR}" "${OUTDIR}"/"${SAMPLE}"/

# Running metaSPAdes
metaspades.py \
	-t "${THREADS}" \
	-m "${MEM}" \
	-k 21,33,55,77 \
	--pe1-1 "${FASTQDECONPAIRED1}" \
	--pe1-2 "${FASTQDECONPAIRED2}" \
	--pe1-s "${FASTQDECONUNPAIRED}" \
	-o "${OUTDIR}"/"${SAMPLE}"/
