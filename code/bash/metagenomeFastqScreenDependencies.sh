#! /bin/bash
# metagenomeFastQScreenDependencies.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Download reference files needed for Fastq Screen.
# Usage: bash metagenomeFastQScreenDependencies.sh MOUSEBOWTIE2INDEX HUMANBOWTIE2INDEX

##################
# Set Script Env #
##################

# Variables defined by user
MOUSEBOWTIE2INDEX=${1:?ERROR: Need to define MOUSEBOWTIE2INDEX.} # Only supply name/path of one index file, will use to find base name
HUMANBOWTIE2INDEX=${2:?ERROR: Need to define HUMANBOWTIE2INDEX.} # Only supply name/path of one index file, will use to find base name

# Other variables
ENVSDATADIR=envs/share/fastq_screen/



####################################################
# Creating BWA Alignment Indices for Mapping Reads #
####################################################

echo PROGRESS: Setting index locations in FastQ Screen configuration file.

# Creating output directories
mkdir -p "${ENVSDATADIR}"/

# Copying example configuration file for use as working configuration file
cp $(realpath $(which fastq_screen) | sed 's:fastq_screen$::')/fastq_screen.conf.example "${ENVSDATADIR}"/fastq_screen.conf

# Pulling bowtie2 index paths/base names
MOUSEINDEXPATH=$(realpath $(echo "${MOUSEBOWTIE2INDEX}"  | sed 's/\(.*\)\..*\..*/\1/'))
HUMANINDEXPATH=$(realpath $(echo "${HUMANBOWTIE2INDEX}"  | sed 's/\(.*\)\..*\..*/\1/'))

# Modifying configuration file with database information
echo >> "${ENVSDATADIR}"/fastq_screen.conf
echo '## Mouse - sequence available from ftp://ftp.ensembl.org/pub/release-95/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz' >> "${ENVSDATADIR}"/fastq_screen.conf
echo -e "DATABASE\tMouse\t"${MOUSEINDEXPATH}"" >> "${ENVSDATADIR}"/fastq_screen.conf
echo >> "${ENVSDATADIR}"/fastq_screen.conf
echo '## Human - sequences available from ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz' >> "${ENVSDATADIR}"/fastq_screen.conf
echo -e "DATABASE\tHuman\t"${HUMANINDEXPATH}"" >> "${ENVSDATADIR}"/fastq_screen.conf

