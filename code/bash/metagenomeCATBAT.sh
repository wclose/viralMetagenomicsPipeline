#! /bin/bash
# metagenomeCATBAT.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Classify viral metagenomic bins by aligning predicted ORFs to the RefSeq viral database.
# Usage: bash metagenomeCATBAT.sh BIN DMND NODES

##################
# Set Script Env #
##################

# Variables defined by user
BIN=${1:?ERROR: Need to define BIN.} # Need to supply the name of one of the bin files so the directory can be pulled from it
DMND=${2:?ERROR: Need to define DMND.} # Location of RefSeq viral diamond database to be used for alignments
NODES=${3:?ERROR: Need to define NODES} # The NCBI taxonomy database nodes file detailing taxid classifications

# Other variables
OUTDIR=data/catbat/
THREADS=$(nproc) # Automatically determines the number of cores based on local resources



###########################################################
# Assigning Taxonomic Classifications to Metagenomic Bins #
###########################################################

echo PROGRESS: Taxonomically classifying vlp metagenomic bins using the RefSeq viral database and CAT/BAT.

# Creating output directory
mkdir -p "${OUTDIR}"/

# Finding paths to relevant directories needed as inputs for CAT/BAT
BINDIR=$(echo "${BIN}" | sed 's/\(.*\/\)bin.*/\1/')
DATABASEDIR=$(echo "${DMND}" | sed 's/\(.*\/\).*/\1/')
TAXONOMYDIR=$(echo "${NODES}" | sed 's/\(.*\/\).*/\1/')

# Classifying bins using viral RefSeq reference database
CAT bins -r 10 -f 0.5 --index_chunks 1 --block_size 6 --nproc "${THREADS}" --sensitive -b "${BINDIR}" -s .fa -d "${DATABASEDIR}" -t "${TAXONOMYDIR}" -o "${OUTDIR}"/metagenome --force --no_log

# Converting taxids into taxonomic classification names
CAT add_names --only_official -i "${OUTDIR}"/metagenome.bin2classification.txt -t "${TAXONOMYDIR}" -o "${OUTDIR}"/metagenome.taxonomy.txt
