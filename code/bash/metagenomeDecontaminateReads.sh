#! /bin/bash
# metagenomeBowtie2.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Use Bowtie2 to remove host, etc. contaminants from read libraries.
# Usage: bash metagenomeBowtie2.sh FASTQTRIMMEDPAIRED1 FASTQTRIMMEDPAIRED2 FASTQTRIMMEDUNPAIRED HOSTINDEX

##################
# Set Script Env #
##################

# Variables defined by user
FASTQTRIMMEDPAIRED1=${1:?ERROR: Need to define FASTQTRIMMEDPAIRED1.} 
FASTQTRIMMEDPAIRED2=${2:?ERROR: Need to define FASTQTRIMMEDPAIRED2.}
FASTQTRIMMEDUNPAIRED=${3:?ERROR: Need to define FASTQTRIMMEDUNPAIRED.}
HOSTINDEX=${4:?ERROR: Need to define HOSTINDEX.} # Supply only one of the index files for pulling the name

# Other variables
THREADS=$(nproc) # Automatically determines the number of cores based on local resources
OUTDIR=data/bowtie2/



###############################################################
# Using Bowtie2 to Remove Reads Mapping to Known Contaminants #
###############################################################

mkdir -p  "${OUTDIR}"

# Pulling name information for naming output files
INDEXBASENAME=$(echo "${HOSTINDEX}" | sed 's/\(.*\)\..*\..*/\1/')
FASTQBASENAME=$(echo "${FASTQTRIMMEDPAIRED1}" | sed 's/.*\/\(.*\)_R.*/\1/')



echo PROGRESS: Removing mouse and human contamination from paired reads.

# Mapping reads to contaminant reference genomes and converting to sam file
bowtie2 -q --very-sensitive-local -p "${THREADS}" \
	-x "${INDEXBASENAME}" \
	-1 "${FASTQTRIMMEDPAIRED1}" \
	-2 "${FASTQTRIMMEDPAIRED2}" \
	-S "${OUTDIR}"/"${FASTQBASENAME}"_paired.sam

# Pulling out unmapped reads with non-secondary alignments, sorting, and converting to bam
# After completing successfully, deletes input file to save space
samtools view -b -f 12 -F 256 "${OUTDIR}"/"${FASTQBASENAME}"_paired.sam | \
	samtools sort -n -o "${OUTDIR}"/"${FASTQBASENAME}"_paired_unmapped_sorted.bam \
	&& rm "${OUTDIR}"/"${FASTQBASENAME}"_paired.sam

# Converting paired read bam file into separate fastq output files
# After completing successfully, deletes input file to save space
bedtools bamtofastq -i "${OUTDIR}"/"${FASTQBASENAME}"_paired_unmapped_sorted.bam \
	-fq "${OUTDIR}"/"${FASTQBASENAME}"_R1_paired_decon.fq \
	-fq2 "${OUTDIR}"/"${FASTQBASENAME}"_R2_paired_decon.fq \
	&& rm "${OUTDIR}"/"${FASTQBASENAME}"_paired_unmapped_sorted.bam

# Compressing output files
gzip "${OUTDIR}"/"${FASTQBASENAME}"_R1_paired_decon.fq
gzip "${OUTDIR}"/"${FASTQBASENAME}"_R2_paired_decon.fq



echo PROGRESS: Removing mouse and human contamination from unpaired reads.

# Mapping reads to contaminant reference genomes and converting to sam file
bowtie2 -q --very-sensitive-local -p "${THREADS}" \
	-x "${INDEXBASENAME}" \
	-U "${FASTQTRIMMEDUNPAIRED}" \
	-S "${OUTDIR}"/"${FASTQBASENAME}"_unpaired.sam

# Pulling out unmapped reads with non-secondary alignments, sorting, and converting to bam
# After completing successfully, deletes input file to save space
samtools view -b -f 4 -F 256 "${OUTDIR}"/"${FASTQBASENAME}"_unpaired.sam | \
	samtools sort -n -o "${OUTDIR}"/"${FASTQBASENAME}"_unpaired_unmapped_sorted.bam \
	&& rm "${OUTDIR}"/"${FASTQBASENAME}"_unpaired.sam

# Converting unpaired read bam file to single fastq output file
# After completing successfully, deletes input file to save space
bedtools bamtofastq -i "${OUTDIR}"/"${FASTQBASENAME}"_unpaired_unmapped_sorted.bam \
	-fq "${OUTDIR}"/"${FASTQBASENAME}"_unpaired_decon.fq \
	&& rm "${OUTDIR}"/"${FASTQBASENAME}"_unpaired_unmapped_sorted.bam

# Compressing output file
gzip "${OUTDIR}"/"${FASTQBASENAME}"_unpaired_decon.fq
