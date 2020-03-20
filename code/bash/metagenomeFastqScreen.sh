#! /bin/bash
# metagenomeFastqScreen.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Do basic check for contamination of reads using Fastq Screen.
# Usage: bash metagenomeFastqScreen.sh FASTQSCREENCONFIG FASTQ1 FASTQ2 FASTQ3 ... FASTQN

##################
# Set Script Env #
##################

# Variables defined by user
FASTQSCREENCONFIG=${1:?ERROR: Need to define FASTQSCREENCONFIG.} # Location of the fastq_screen.conf file with database locations, etc.
FASTQ=${@:2}


# Other variables
OUTDIR=data/fastq_screen/



##############################################
# Contaminant Read Mapping with Fastq Screen #
##############################################

# Running Fastq Screen on raw reads (pulls from path)
if echo "${FASTQ[0]}" | grep -q "raw"; then

	echo PROGRESS: Running Fastq Screen on raw reads.
	
	# Making output dir if it doesn't exist yet
	mkdir -p "${OUTDIR}"/raw/

	# Running fastq screen on the provided sequence files
	fastq_screen --subset 100000 --force --outdir "${OUTDIR}"/raw/ --conf "${FASTQSCREENCONFIG}" --aligner bowtie2 $(echo "${FASTQ[@]}")

# Running Fastq Screen on decontaminated reads (pulls from path)
elif echo "${FASTQ[0]}" | grep -q "bowtie2"; then 

	echo PROGRESS: Running Fastq Screen on trimmed and decontaminated reads.

	# Making output dir if it doesn't exist yet
	mkdir -p "${OUTDIR}"/after_decon/

	# Running fastq screen on the provided sequence files
	fastq_screen --subset 100000 --force --outdir "${OUTDIR}"/after_decon/ --conf "${FASTQSCREENCONFIG}" --aligner bowtie2 $(echo "${FASTQ[@]}")

fi
