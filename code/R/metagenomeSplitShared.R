#!/usr/bin/env Rscript
# metagenomeRemoveSamples.R
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Remove specific samples from shared files similar to mothur `remove.groups()`.
# Usage: Rscript metagenomeRemoveSamples.R SAMPLE1 SAMPLE2 SAMPLE3 ... SAMPLEN

# Setting environment -----------------------------------------------------

# Parsing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Variables defined by user
sharedFile <- args[1]
sampleNames <- as.character(args[-1]) # Names of samples to remove from shared files

# Checking if inputs are set
if (!file.exists(sharedFile)) {
  stop(paste(sharedFile, "does not exist."))
} else if (all(is.na(sampleNames))) {
  stop("Need to define sampleNames.")
}


# Other variables
outDir <- "data/shared/"



# Loading dependencies ----------------------------------------------------

library(tidyverse)



# Functions ---------------------------------------------------------------

# Function for creating shared file of only samples
get_sample_shared <- function(shared_df, sample_names_chr) {
  
  # Progress message
  message("PROGRESS: Removing samples from shared file.")
  
  # Removing samples from shared file
  filtered_shared <- shared_df %>% 
    filter(!(sample %in% sample_names_chr))
  
  return(filtered_shared)
  
}


# Function for creating shared file of only controls
get_control_shared <- function(shared_df, sample_names_chr) {
  
  # Progress message
  message("PROGRESS: Removing samples from shared file.")
  
  # Removing samples from shared file
  filtered_shared <- shared_df %>% 
    filter(sample %in% sample_names_chr)
  
  return(filtered_shared)
  
}



# Analysis ----------------------------------------------------------------

message("PROGRESS: Splitting master shared file into separate shared files for samples versus controls.")

# Creating output directory if it doesn't already exist
dir.create(outDir, recursive = TRUE, showWarnings=FALSE)

# Reading in inputs
shared <- read_tsv(sharedFile, col_types = cols()) %>% 
  mutate(sample = as.character(sample))

# Removing controls from shared file
sample_shared <- get_sample_shared(shared, sampleNames)

# Removing samples from shared file
control_shared <- get_control_shared(shared, sampleNames)

message("PROGRESS: Writing out filtered shared files.")

# Parsing apart shared file name and creating name for output shared file
sample_output_file <- str_replace(str_extract(sharedFile, "metagenome[[:alpha:].]*shared$"), "(metagenome.*shared)", "sample.\\1")
control_output_file <- str_replace(str_extract(sharedFile, "metagenome[[:alpha:].]*shared$"), "(metagenome.*shared)", "control.\\1")

# Write out the shared files
write_tsv(sample_shared, path = paste0(outDir, sample_output_file))
write_tsv(control_shared, path = paste0(outDir, control_output_file))
