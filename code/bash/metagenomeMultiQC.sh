#! /bin/bash
# metagenomeMultiQC.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Generate quality control reports for the various steps of read processing and assembly using MultiQC.
# Usage: bash metagenomeMultiQC.sh DATAFILE1 DATAFILE2 DATAFILE3 ... DATAFILEN

##################
# Set Script Env #
##################

# Variables defined by user
DATAFILES=${@:?ERROR: Need to define DATAFILES.}

# Other variables
OUTDIR=data/multiqc/



##########################################
# Combining Read QC Reports with MultiQC #
##########################################

echo PROGRESS: Running MultiQC and aggregating QC output.

# Create dir if it doesn't already exist, all reports will be deposited here
mkdir -p "${OUTDIR}"

# Naming output file based on data source
# Separate reports are created for each stage to keep file sizes manageable while forcing plots to be interactive
# Raw reads
if echo "${DATAFILES[0]}" | grep -q "raw"; then

	echo PROGRESS: Combining QC reports from raw reads.

	# Naming output file
	export OUTPUT="raw_multiqc_report.html"

# Creating separate report for trimmed reads
elif echo "${DATAFILES[0]}" | grep -q "after_trimmomatic"; then 

	echo PROGRESS: Combining QC reports from trimmed reads.

	# Naming output file
	export OUTPUT="trimmed_multiqc_report.html"

# Creating separate report for decontaminated reads
elif echo "${DATAFILES[0]}" | grep -q "after_decon"; then 

	echo PROGRESS: Combining QC reports from trimmed reads after decontamination.

	# Naming output file
	export OUTPUT="decon_multiqc_report.html"

# Creating separate report for assembled contigs
elif echo "${DATAFILES[0]}" | grep -q "quast"; then 

	echo PROGRESS: Combining QC reports from assembled contigs.

	# Naming output file
	export OUTPUT="assembly_multiqc_report.html"

fi

# Forcing plots to be interactive, removing multiqc data dir, and overwriting existing reports
multiqc $(echo "${DATAFILES[@]}") --interactive --no-data-dir  --force --outdir "${OUTDIR}"/ --filename "${OUTPUT}"
