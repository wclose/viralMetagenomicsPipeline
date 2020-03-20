#! /bin/bash
# metagenomeCATBATDependencies.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Download and prepare databases for running CAT/BAT to classify metagenome bins.
# Usage: bash metagenomeCATBATDependencies.sh

##################
# Set Script Env #
##################

# Other variables
THREADS=$(nproc) # Automatically determines the number of cores based on local resources
ENVSDATADIR=envs/share/catbat/



##############################################################
# Preparing Databases Required for Classification by CAT/BAT #
##############################################################

echo PROGRESS: Preparing NCBI taxonomy database files.

# Creating output directory
mkdir -p "${ENVSDATADIR}"/taxonomy/

# Downloading NCBI tax database
wget -N -P "${ENVSDATADIR}"/taxonomy/ ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

# Extracting tax database
tar -xvzf "${ENVSDATADIR}"/taxonomy/taxdump.tar.gz -C "${ENVSDATADIR}"/taxonomy/

# Downloading protein accession to tax id file
wget -N -P "${ENVSDATADIR}"/taxonomy/ ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz



echo PROGRESS: Preparing RefSeq viral database files.

# Creating output directory
mkdir -p "${ENVSDATADIR}"/database/

# Downloading RefSeq viral database adding time stamp so I know when the database was created 
wget -r -np -nd -N -P "${ENVSDATADIR}"/database/ ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/ && date > "${ENVSDATADIR}"/database/timestamp.txt

# Combining all of the viral protein fasta files together into a single database
zcat $(find "${ENVSDATADIR}"/database/ -regex ".*protein.faa.gz" | grep -v "nonredundant") > "${ENVSDATADIR}"/database/refseq_viral.faa && gzip "${ENVSDATADIR}"/database/refseq_viral.faa

# Renaming RefSeq viral file to trick CAT/BAT into using it for downstream steps (CAT/BAT assumes you are using entire BLAST non-redundant database)
cp "${ENVSDATADIR}"/database/refseq_viral.faa.gz "${ENVSDATADIR}"/database/nr.gz



echo PROGRESS: Creating DIAMOND database from RefSeq viral database.

# Creating diamond database from custom faa library
diamond makedb --in "${ENVSDATADIR}"/database/refseq_viral.faa.gz -d "${ENVSDATADIR}"/database/refseq_viral -p "${THREADS}" --quiet



echo PROGRESS: Preparing remaining CAT/BAT database files.

# Using CAT/BAT to create the remaining database files needed for bin classification
CAT prepare --existing -d "${ENVSDATADIR}"/database/ -t "${ENVSDATADIR}"/taxonomy/ --no_log
