#! /bin/bash
# metagenomeCurateContigs.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Remove non-viral reads, as determined by VirFinder and VirSorter, from combined contig library.
# Usage: bash metagenomeCurateContigs.sh FASTA VIRSORTERPRED VIRFINDERPRED

##################
# Set Script Env #
##################

# Variables defined by user
FASTA=${1:?ERROR: Need to define FASTA.} # Deduped contig library fasta file
VIRSORTERPRED=${2:?ERROR: Need to define VIRSORTERPRED.} # List of predicted viral contigs from VirSorter
VIRFINDERPRED=${3:?ERROR: Need to define VIRFINDERPRED.} # List of predicted viral contigs from VirFinder


# Other variables
OUTDIR=data/curate/



##################################################
# Removing Non-Viral Contigs from Contig Library #
##################################################

# Making output directory if it doesn't exist
mkdir -p "${OUTDIR}"/

echo PROGRESS: Combining lists of predicted viral contigs from VirSorter and VirFinder.

# Combining lists of putative "viral" contig names and removing duplicates
sort -u "${VIRSORTERPRED}" "${VIRFINDERPRED}" > "${OUTDIR}"/predicted_viral_contigs.txt



echo PROGRESS: Removing non-viral contigs from library based on VirSorter and VirFinder output.

# Using VirFinder/VirSorter output to remove non-viral contigs from library.
awk 'NR==FNR{viral[$1]; next} NF>1{f=($2 in viral)} f' FS=' ' "${OUTDIR}"/predicted_viral_contigs.txt  FS='>' "${FASTA}" > "${OUTDIR}"/viral_contig_library.fa
