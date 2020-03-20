#! /bin/bash
# metagenomeMetaBAT.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Cluster contigs into metagenomic bins using read alignment statistics and MetaBAT2.
# Usage: bash metagenomeMetaBAT.sh CONTIGLIBRARY SEED CONTIGALIGNEDBAM1 CONTIGALIGNEDBAM2 CONTIGALIGNEDBAM3 ... CONTIGALIGNEDBAMN

##################
# Set Script Env #
##################

# Variables defined by user
CONTIGLIBRARY=${1:?ERROR: Need to define CONTIGLIBRARY.}
SEED=${2:?ERROR: Need to define SEED.} # Setting seed for reproducibility
CONTIGALIGNEDBAM=${@:3} # The list of bam files aligned to the deduped contig library

# Error message triggered when only one file was supplied which means either the library or bam files were left out
if [ "$#" -lt 3 ]; then
	echo ERROR: Need to define CONTIGALIGNEDBAM.
	exit 1
fi

# Other variables
OUTDIR=data/metabat/



##########################################
# De Novo Binning Contigs Using MetaBAT2 #
##########################################

echo PROGRESS: Generating MetaBAT2 depth file.

# Creating output dir if it doesn't exist yet
mkdir -p "${OUTDIR}"/ "${OUTDIR}"/bins/

# Creating depth file summarizing sample base coverage mean and variance across contig library
jgi_summarize_bam_contig_depths --outputDepth "${OUTDIR}"/bins/depth.txt $(echo "${CONTIGALIGNEDBAM[@]}")



echo PROGRESS: Binning contigs using MetaBAT2.

# Binning contigs
metabat2 -i "${CONTIGLIBRARY}" \
	-a "${OUTDIR}"/bins/depth.txt \
	-o "${OUTDIR}"/bins/bin \
	-m 1500 \
	-s 1 \
	--seed "${SEED}" \
	-v
