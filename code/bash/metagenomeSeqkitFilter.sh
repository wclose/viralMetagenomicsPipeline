#! /bin/bash
# metagenomeSeqkitFilter.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Uses Seqkit to remove contigs less than 1.5 kb in length which speeds up downstream processes. Doesn't affect
# results as only 1.5 kb contigs were used in binnning, etc.
# Usage: bash metagenomeSeqkitFilter.sh FASTA MINLENGTH

##################
# Set Script Env #
##################

# Variables defined by user
FASTA=${1:?ERROR: Need to define FASTA.} # Deduped contig library fasta file
MINLENGTH=${2:?ERROR: Need to define MINLENGTH.} # Minimum length of contigs to keep

# Other variables
OUTDIR=data/seqkit/



########################################
# Removing Contigs Below Length Cutoff #
########################################

# Making output directory if it doesn't exist
mkdir -p "${OUTDIR}"/

echo PROGRESS: Filtering out sequences below "${MINLENGTH}" bp.

# Filtering sequences
seqkit seq -m "${MINLENGTH}" "${FASTA}" > "${OUTDIR}"/contig_library_deduped_filtered.fa

# Storing the number of sequences in the input and output for debugging
INPUTSEQNO=$(grep ">" "${FASTA}" | wc -l)
OUTPUTSEQNO=$(grep ">" "${OUTDIR}"/contig_library_deduped_filtered.fa | wc -l)

echo Number of sequences input: "${INPUTSEQNO}"
echo Number of sequences surviving: "${OUTPUTSEQNO}"
