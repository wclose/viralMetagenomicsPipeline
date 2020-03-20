#!/usr/bin/env Rscript
# metagenomeCombineIdxstats.R
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Combine output from Samtools idxstats into a tidy format so it can be used by other scripts.
# Usage: Rscript metagenomeCombineIdxStats.R IDXFILE1 IDXFILE2 IDXFILE3 ... IDXFILEN

# Setting environment -----------------------------------------------------

# Parsing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Variables defined by user
idxFiles <- args

# Checking if inputs are set
if (!all(file.exists(idxFiles))) {
  stop(paste("The following files do not exist:", paste(idxFiles, collapse = ", ")))
}


# Other variables
outDir <- "data/idxstats/"



# Loading dependencies ----------------------------------------------------

library(tidyverse)



# Functions ---------------------------------------------------------------

# Function for reading and processing read recruitment output from samtools idxstats
tidy_idxstats <- function(idx_path_chr){
  
  idx_cols <- c("seq_name", "seq_length", "mapped_reads", "unmapped_reads") # The idxstats output doesn't contain col names so need to specify
  sample_name <- str_replace(idx_path_chr, ".*/(\\d+)_.*", "\\1") # Pulling the sample number from the file path
  
  idx_data <- read_tsv(idx_path_chr, col_names = idx_cols, col_types = cols()) %>% # Reading in the data and assigning col names
    mutate(sample = sample_name, # Creating col for sample names
           bin = paste0("Bin",str_replace(seq_name, "^bin\\.(\\d+)\\.fa_.*", "\\1")), # Creating col for bin number by extracting number from seq_name
           seq_name = str_replace(seq_name, "^bin\\.\\d+\\.fa_", "")) %>% # Removing bin info from seq_name (would be redundant otherwise)
    filter(seq_length != "0") %>% # Removing rows for unmapped reads
    select(sample, bin, everything()) # Making the output human friendly
  
  return(idx_data)
  
}


# Function for formatting and joining all of the idxstats output
combine_idxstats <- function(idx_path_chr) {
  
  # Progress message
  message("PROGRESS: Combining idxstats files.")
  
  # Combining idxstats files
  combined_idxstats <- map_df(idx_path_chr, tidy_idxstats) %>% 
    group_by(sample, bin) %>% 
    summarize(bin_length = sum(seq_length), # Calculating total bin length in bp
              reads_mapped = sum(mapped_reads), # Calculating total number of reads (NOTE: NOT read pairs) mapped to each bin
              reads_unmapped = sum(unmapped_reads)) %>% # Calculating total number of reads that didn't map but whose mate did map
    ungroup()
  
  return(combined_idxstats)
  
}



# Analysis ----------------------------------------------------------------

message("PROGRESS: Combining and tidying bin read recruitment stats.")

# Creating output directory if it doesn't already exist
dir.create(outDir, recursive = TRUE, showWarnings=FALSE)

# Combining idxstats files
combined_idxstats <- combine_idxstats(idxFiles)

message("PROGRESS: Writing out combined read recruitment stats.")

# Write out idxstats file
write_tsv(combined_idxstats, path = paste0(outDir, "combined_idxstats.tsv"))
