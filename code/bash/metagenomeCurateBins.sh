#! /bin/bash
# metagenomeCurateBins.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Remove viral bins contaminated by bacterial/archaea as determined by CheckM.
# Usage: bash metagenomeCurateBins.sh CHECKMMETRIC BINFASTA

##################
# Set Script Env #
##################

# Variables defined by user
CHECKMMETRIC=${1:?ERROR: Need to define CHECKMMETRIC.} # Location of CheckM report detailing which bins are free from contamination
BINFASTA=${2:?ERROR: Need to define BINFASTA.} # A single bin fasta file to be used for pulling the file path from

# Other variables
OUTDIR=data/metabat/



############################################################################
# Renaming and Concatenating Bin Sequences Without Bacterial Contamination #
############################################################################

echo PROGRESS: Adding bin names to FASTA headers and concatenating sequences from uncontaminated bins.

# Making output directories if they don't already exist
mkdir -p "${OUTDIR}"/ "${OUTDIR}"/library/

# Overwriting the final output file if it already exists or creating an empty file if it doesn't
echo -n > "${OUTDIR}"/library/bin_library_decon.fa

# Parses CheckM taxonomy metrics and filters out any bins with non-viral taxonomic assignments then adds '.fa' extension
DECONBINS=$(awk -F: '/^[^-]/' "${CHECKMMETRIC}" | awk '/root/ {print ""$1".fa"}' | sort)

# Pullting path to bin fasta files for finding later (avoids issues with argument list being too long)
BINPATH=$(echo "${BINFASTA}" | sed 's:\(.*/\)bin\..*fa:\1:')

# Labeling decontaminated bins with bin number and concatenating them together for downstream use
for BIN in $(find "${BINPATH}" -name "*fa"); do

	# Pulling bin name from file path
	BINNAME=$(echo "${BIN}" | sed 's/.*\/\(.*\.fa\)/\1/')

	# If the supplied bin fasta is in the decontaminated bin list, add the bin name to sequence headers and concatenate together
	if [[ "${DECONBINS}" =~ "${BINNAME}" ]]; then
		sed 's/>\(.*\)/>'"${BINNAME}"'_\1/' "${BIN}" >> "${OUTDIR}"/library/bin_library_decon.fa
	fi

done

# Creating list of contig names in the bin library
awk '/>/ { print substr($1,2) }' "${OUTDIR}"/library/bin_library_decon.fa > "${OUTDIR}"/library/bin_library_decon_names.txt

# Reporting number of bins/contigs as input vs output
echo Input Bins: $(find "${BINPATH}" -name "*fa" | wc -l)
echo Input Contigs: $(find "${BINPATH}" -name "*fa" -exec grep ">" {} + | wc -l)
echo Surviving Bins: $(echo "${DECONBINS}" | wc -l)
echo Surviving Contigs: $(grep ">" "${OUTDIR}"/library/bin_library_decon.fa | wc -l)
