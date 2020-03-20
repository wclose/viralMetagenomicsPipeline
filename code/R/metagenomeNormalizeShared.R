#!/usr/bin/env Rscript
# metagenomeNormalizeShared.R
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Normalizing shared files based on total bin length.
# Usage: Rscript metagenomeNormalizeShared.R SHAREDFILE IDXSTATSFILE

# Setting environment -----------------------------------------------------

# Parsing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Variables defined by user
sharedFile <- args[1] # Raw shared file
idxstatsFile <- args[2] # Combined idxstats file

# Checking if inputs are set
if (!file.exists(sharedFile)) {
  stop(paste(sharedFile, "does not exist."))
} else if (!file.exists(idxstatsFile)) {
  stop(paste(idxstatsFile, "does not exist."))
}

# Other variables
outDir <- "data/shared/"



# Loading dependencies ----------------------------------------------------

library(tidyverse)

# Loading common functions
source("code/R/metagenomeFunctions.R")



# Analysis ----------------------------------------------------------------

message("PROGRESS: Normalizing metagenomic shared file based on total bin lengths.")

# Creating output directory if it doesn't already exist
dir.create(outDir, recursive = TRUE, showWarnings=FALSE)

# Reading in data files
shared <- read_tsv(sharedFile, col_types = cols()) %>% 
  mutate(sample = as.character(sample))
idxstats <- read_tsv(idxstatsFile, col_types = cols()) %>% 
  mutate(sample = as.character(sample))

# Tidying shared file in preparation for normalizing
tidy_shared <- shared %>% 
  gather(key = "bin", value = "count", matches("^Bin\\d+"))

# Normalizing shared file
norm_shared <- normalize_shared(tidy_shared, idxstats)

# Creating name for output file
out_file <- str_replace(str_extract(sharedFile, "[[:alpha:].]*metagenome.raw.shared"), "(^.*).raw.shared", "\\1.norm.shared")

message("PROGRESS: Writing out normalized shared file.")

# Write out the normalized shared file
write_tsv(norm_shared, path = paste0(outDir, out_file))
