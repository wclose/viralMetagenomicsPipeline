#! /bin/bash
# metagenomeBowtie2Dependencies.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Download and prepare Bowtie2 indices for removing contaminants from reads and/or mapping reads with Fastq Screen.
# Usage: bash metagenomeBowtie2Dependencies.sh

##################
# Set Script Env #
##################

# Other variables
THREADS=$(nproc) # Automatically determines the number of cores based on local resources
ENVSDATADIR=envs/share/bowtie2/



################################################################
# Downloading Contaminant Reference Files and Building Indices #
################################################################

echo PROGRESS: Downloading human and mouse reference genomes.

# Creating output directories
mkdir -p "${ENVSDATADIR}"/reference/

# Downloading Mus Musculus C57BL6 reference genome
wget -P "${ENVSDATADIR}"/reference/ ftp://ftp.ensembl.org/pub/release-95/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

# Downloading Homo Sapiens reference genome
wget -P "${ENVSDATADIR}"/reference/ ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz



echo PROGRESS: Creating Bowtie2 alignment indices from human and mouse reference genomes.

# Creating output directories
mkdir -p "${ENVSDATADIR}"/index/

# Making Mus Musculus C57BL6 index for FastQ Screen
bowtie2-build --threads $THREADS \
	"${ENVSDATADIR}"/reference/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz \
	"${ENVSDATADIR}"/index/m_musculus_index

# Making Homo Sapiens index for FastQ Screen
bowtie2-build --threads $THREADS \
	"${ENVSDATADIR}"/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
	"${ENVSDATADIR}"/index/h_sapiens_index

# Combining Mus Musculus C57BL6 and Homo Sapiens genomes into single reference for removing host contamination
bowtie2-build --threads $THREADS --large-index \
	"${ENVSDATADIR}"/reference/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz,"${ENVSDATADIR}"/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
	"${ENVSDATADIR}"/index/m_musculus_h_sapiens_combined_index
