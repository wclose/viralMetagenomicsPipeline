#! /bin/bash
# metagenomeSamtoolsIdxstats.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Uses Samtools to calculate read recruitment statistics after mapping reads to contigs/bins.
# Usage: bash metagenomeSamtoolsIdxstats.sh BAM

##################
# Set Script Env #
##################

# Variables defined by user
BAM=${1:?ERROR: Need to define BAM.}

# Other variables
OUTDIR=data/idxstats/



###########################################
# Calculating Read Recruitment Statistics #
###########################################

echo PROGRESS: Calculating read recruitment stats.

# Creating dir if it doesn't already exist
mkdir -p "${OUTDIR}"

# Pulling filename for naming output
FILENAME=$(echo "${BAM}" | sed 's/.*\/\(.*\)\.bam/\1/')

# Calculating recruitment stats
samtools idxstats "${BAM}" > "${OUTDIR}"/"${FILENAME}".idxstats.txt
