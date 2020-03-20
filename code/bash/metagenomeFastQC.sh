#! /bin/bash
# metagenomeFastQC.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Run FastQC to check read characteristics such as quality, adapter content, repeats, etc.
# Usage: bash metagenomeFastQC.sh FASTQ1 FASTQ2 FASTQ3 ... FASTQN

##################
# Set Script Env #
##################

# Variables defined by user
FASTQ=${@:?ERROR: Need to define FASTQ}

# Other variables
OUTDIR=data/fastqc/



#####################################
# Checking Read Quality with FastQC #
#####################################

# Setting subdirectory based on name of reads
# Raw reads
if echo "${FASTQ[0]}" | grep -q "raw"; then

	echo PROGRESS: Running FastQC on raw reads to check initial read quality.
	
	# Setting subdirectory
	export OUTSUBDIR="${OUTDIR}"/raw/

# Trimmed reads
elif echo "${FASTQ[0]}" | grep -q "trimmomatic"; then 

	echo PROGRESS: Running FastQC on trimmed reads to check read quality.

	# Setting subdirectory
	export OUTSUBDIR="${OUTDIR}"/after_trimmomatic/

# Decontaminated reads
elif echo "${FASTQ[0]}" | grep -q "bowtie2"; then 

	echo PROGRESS: Running FastQC on trimmed reads after decontamination to check read quality.

	# Setting subdirectory
	export OUTSUBDIR="${OUTDIR}"/after_decon/

fi

# Making output dir if it doesn't exist yet
mkdir -p "${OUTSUBDIR}"

# Running FastQC on reads
fastqc -o "${OUTSUBDIR}" $(echo "${FASTQ[@]}")
