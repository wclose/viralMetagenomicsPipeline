#! /bin/bash
# metagenomeCheckMDependencies.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Downloading and preparing databases for assessing bin contamination with CheckM.
# Usage: bash metagenomeCheckMDependencies.sh

##################
# Set Script Env #
##################

# Other variables
ENVSDATADIR=envs/share/checkM/
ENVSDATADIRFULL=$(realpath "${ENVSDATADIR}"/)



##########################################################
# Preparing CheckM Data Files and Setting Root Directory #
##########################################################

echo PROGRESS: Preparing CheckM data files.

# Creating output dir if it doesn't exist yet
mkdir -p "${ENVSDATADIR}"/

# Download the data to the data folder
wget -P "${ENVSDATADIR}"/ https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz

# Unpacking it to the desired location
tar -xzf "${ENVSDATADIR}"/checkm_data_2015_01_16.tar.gz -C "${ENVSDATADIR}"/

# Setting the data directory within CheckM
echo "${ENVSDATADIRFULL}"/ | checkm data setRoot "${ENVSDATADIRFULL}"/
