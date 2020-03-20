#! /bin/bash
# metagenomeVirFinderDependencies.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Preparing the eukaryotic/prokaryotic model for predicting viral contigs with VirFinder.
# Usage: bash metagenomeVirFinderDependencies.sh

##################
# Set Script Env #
##################

# Other variables
ENVSDATADIR=envs/share/virfinder/



############################
# Predicting Viral Contigs #
############################

echo PROGRESS: Downloading VirFinder eukaryotic/prokaryotic model.

# Creating output dir if it doesn't exist yet
mkdir -p "${ENVSDATADIR}"/tmp/

# Downloading data files
wget -N -P "${ENVSDATADIR}"/tmp/ https://github.com/jessieren/VirFinder/archive/v1.1.tar.gz

# Decompressing data files
tar -xvzf "${ENVSDATADIR}"/tmp/v1.1.tar.gz -C "${ENVSDATADIR}"/tmp/

# Moving model file to top level
mv "${ENVSDATADIR}"/tmp/VirFinder-1.1/EPV/VF.modEPV_k8.rda "${ENVSDATADIR}"/

# Cleaning up unused files
rm -r "${ENVSDATADIR}"/tmp/
