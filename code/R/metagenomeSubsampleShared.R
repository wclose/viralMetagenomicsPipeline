#!/usr/bin/env Rscript
# metagenomeSubsampleShared.R
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Subsampling shared file to ensure equal sampling depth across samples.
# Usage: Rscript metagenomeSubsampleShared.R SHAREDFILE IDXSTATSFILE SUBTHRESH

# Setting environment -----------------------------------------------------

# Parsing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Variables defined by user
sharedFile <- args[1] # Raw shared file
idxstatsFile <- args[2] # Combined idxstats file
subthresh <- suppressWarnings(as.integer(args[3])) # Number of reads to subsample to

# Checking if inputs are set
if (!file.exists(sharedFile)) {
  stop(paste(sharedFile, "does not exist."))
} else if (!file.exists(idxstatsFile)) {
  stop(paste(idxstatsFile, "does not exist."))
} else if (is.na(subthresh)) {
  stop("Input subthresh must be defined and an integer")
}

# Other variables
outDir <- "data/shared/"



# Loading dependencies ----------------------------------------------------

library(tidyverse)

# Loading common functions
source("code/R/metagenomeFunctions.R")



# Analysis ----------------------------------------------------------------

message("PROGRESS: Subsampling shared file and normalizing output based on total bin length.")

# Creating output directory if it doesn't already exist
dir.create(outDir, recursive = TRUE, showWarnings=FALSE)

# Reading in data files
shared <- read_tsv(sharedFile, col_types = cols()) %>% 
  mutate(sample = as.character(sample))
idxstats <- read_tsv(idxstatsFile, col_types = cols()) %>% 
  mutate(sample = as.character(sample))

# Removing samples without enough reads for subsampling
filtered_shared <- filter_shared_by_reads(shared, subthresh)

# Setting seed to make analysis reproducible
set.seed(20170415)

# Subsampling the shared file and normalizing based on read abundances
sub_shared <- subsample_shared(filtered_shared, min(filtered_shared$numReads)) %>% 
  normalize_shared(idxstats)

# Creating name for output file
out_file <- str_replace(str_extract(sharedFile, "[[:alpha:].]*metagenome.raw.shared"), "(^.*).raw.shared", "\\1.subsample.norm.shared")

message("PROGRESS: Writing out subsampled shared file.")

# Write out the subsampled shared file
write_tsv(sub_shared, path = paste0(outDir, out_file))
