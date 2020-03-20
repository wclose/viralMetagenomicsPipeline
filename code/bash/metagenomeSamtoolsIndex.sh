#! /bin/bash
# metagenomeSamtoolsIndex.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Uses Samtools to create indices of bam files for faster calculation of read recruitment stats.
# Usage: bash metagenomeSamtoolsIndex.sh BAM

##################
# Set Script Env #
##################

# Variables defined by user
BAM=${1:?ERROR: Need to define BAM.}

# Don't need to define OUTDIR because only adds new suffix to input BAM


#################################################
# Indexing Aligned Bam Files for Downstream Use #
#################################################

echo PROGRESS: Indexing bam files.

samtools index -b "${BAM}" "${BAM}".bai
