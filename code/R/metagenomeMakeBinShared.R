#!/usr/bin/env Rscript
# metagenomeMakeBinShared.R
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Use idxstats output from Samtools to calculate read recruitment and generate a mothur-style shared file.
# Usage: Rscript metagenomeMakeBinShared.R COMBINEDIDXSTATSFILE

# Setting environment -----------------------------------------------------

# Parsing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Variables defined by user
combinedIdxstatsFile <- args[1]

# Checking if inputs are set
if (!file.exists(combinedIdxstatsFile)) {
  stop(paste(combinedIdxstatsFile, "does not exist."))
} 


# Other variables
outDir <- "data/shared/"



# Loading dependencies ----------------------------------------------------

library(tidyverse)



# Functions ---------------------------------------------------------------

# Function for creating shared file by reformatting idxstats df
make_shared <- function(combined_idxstats_df) {
  
  # Progress message
  message("PROGRESS: Creating metagenomic shared file.")
  
  # Making shared file
  shared <- combined_idxstats_df %>% 
    select(sample, bin, reads_mapped) %>% # Read counts will only be based on reads that mapped
    group_by(sample) %>% 
    mutate(numReads = sum(reads_mapped)) %>% # Calculating the total number of reads for each sample
    ungroup() %>% 
    mutate(numBins = length(unique(bin))) %>% # Calculating the total number of bins for each sample in the shared file
    spread(key = bin, value = reads_mapped) # Spreading data out into shared file format
  
  return(shared)
  
}



# Analysis ----------------------------------------------------------------

message("PROGRESS: Using read recruitment stats to generate master metagenomic shared file.")

# Creating output directory if it doesn't already exist
dir.create(outDir, recursive = TRUE, showWarnings=FALSE)

# Reading in idxstats file
idx <- read_tsv(combinedIdxstatsFile, col_types = cols())

# Creating shared file and replacing sequencing core sample names with actual sample names
shared <- make_shared(idx)

message("PROGRESS: Writing out new master shared file.")

# Write out the shared file
write_tsv(shared, path = paste0(outDir, "metagenome.raw.shared"))
