#! /bin/bash
# metagenomeVirSorterDependencies.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Preparing databases for predicting viral contigs with VirSorter.
# Usage: bash metagenomeVirSorterDependencies.sh

##################
# Set Script Env #
##################

# Other variables
ENVSDATADIR=envs/share/virsorter/



############################
# Predicting Viral Contigs #
############################

echo PROGRESS: Downloading VIRSorter database files.

# Creating output dir if it doesn't exist yet
mkdir -p "${ENVSDATADIR}"/

# Downloading data files
wget -N -P "${ENVSDATADIR}"/ https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz

# Decompressing data files
tar -xvzf "${ENVSDATADIR}"/virsorter-data-v2.tar.gz -C "${ENVSDATADIR}"/
