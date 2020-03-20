#! /bin/bash
# metagenomeCheckM.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Use CheckM to assess bacterial/archaeal contamination of metagenomic bins.
# Usage: bash metagenomeCheckM.sh BIN

##################
# Set Script Env #
##################

# Variables defined by user
BIN=${1:?ERROR: Need to define BIN.} # Need to supply the name of one of the bin files so the directory can be pulled from it

# Other variables
THREADS=$(nproc) # Automatically determines the number of cores based on local resources
OUTDIR=data/checkm/



##################################################
# Detecting Bacterial Contamination Using CheckM #
##################################################

# Pulling bin path from supplied bin
BINDIR=$(echo "${BIN}" | sed 's/\(.*\/\)bin.*/\1/')

echo PROGRESS: Detecting bacterial bin contamination.

# Creating output dir if it doesn't exist yet
mkdir -p "${OUTDIR}"/

# Assigning taxonomic classification to bins
checkm lineage_wf -t "${THREADS}" -x fa --force_overwrite -f "${OUTDIR}"/checkm_metrics.txt "${BINDIR}" "${OUTDIR}"
