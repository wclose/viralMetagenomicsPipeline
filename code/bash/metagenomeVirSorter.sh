#! /bin/bash
# metagenomeVirSorter.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Identify contigs identified as high-confidence by VirSorter.
# Usage: bash metagenomeVirSorter.sh FASTA VIRSORTERREADME

##################
# Set Script Env #
##################

# Variables defined by user
FASTA=${1:?ERROR: Need to define FASTA.} # Need to supply contig fasta
VIRSORTERREADME=${2:?ERROR: Need to define VIRSORTERREADME.} # Need specify location of README from VIRSorter database

# Other variables
THREADS=$(nproc) # Automatically determines the number of cores based on local resources
OUTDIR=data/virsorter/



############################
# Predicting Viral Contigs #
############################

echo PROGRESS: Predicting viral contigs with VIRSorter.

# Creating output dir if it doesn't exist yet
mkdir -p "${OUTDIR}"/

# Pulling name of VIRSorter database directory from database README
VIRSORTERDATA=$(echo "${VIRSORTERREADME}" | sed 's:\(.*/\).*:\1:')

# Predicting viral contigs using virome decontamination mode
wrapper_phage_contigs_sorter_iPlant.pl --virome -f "${FASTA}" --db 1 --wdir "${OUTDIR}"/ --ncpu "${THREADS}" --data-dir "${VIRSORTERDATA}"

# Creating list of high-confidence phage contigs (Categories 1 and 2) and reformatting names back to the original sequence names
awk '/\#\# 1/{flag=1;next}/\#\# 3/{flag=0}flag' "${OUTDIR}"/VIRSorter_global-phage-signal.csv | \
	grep -v "##" | \
	awk '{print $1}' FS=',' | \
	sed 's:VIRSorter_\(.*\):\1:;s:\(.*\)_:\1\.:;s:\(.*\)-circular:\1:' > "${OUTDIR}"/virsorter_predicted_contigs.txt

# Appending list of high-confidence prophage contigs (Categories 4 and 5) and reformatting names back to the original sequence names
awk '/\#\# 4/{flag=1;next}/\#\# 6/{flag=0}flag' "${OUTDIR}"/VIRSorter_global-phage-signal.csv | \
	grep -v "##" | \
	awk '{print $1}' FS=',' | \
	sed 's:VIRSorter_\(.*\):\1:;s:\(.*\)_:\1\.:;s:\(.*\)-circular:\1:' >> "${OUTDIR}"/virsorter_predicted_contigs.txt
