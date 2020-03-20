#!/usr/bin/env Rscript
# metagenomeBeta.R
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Calculate beta diversity by iteratively subsampling normalized shared files.
# Usage: Rscript metagenomeBeta.R SHAREDFILE IDXSTATS SUBTHRESH BETAMETRIC

# Setting environment -----------------------------------------------------

# Parsing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Variables defined by user
sharedFile <- args[1] # Raw metagenomic shared file
idxstatsFile <- args[2] # Combined idxstats file
subthresh <- suppressWarnings(as.integer(args[3])) # Minimum number of reads to subsample to
betaMetric <- args[4] # Beta diversity metric to use. See vegan:vegdist() for options

# Checking if inputs are set
if (!file.exists(sharedFile)) {
  stop(paste(sharedFile, "does not exist."))
} else if (!file.exists(idxstatsFile)) {
  stop(paste(idxstatsFile, "does not exist."))
} else if (is.na(subthresh)) {
  stop("Input subthresh must be defined and an integer")
} else if (is.na(betaMetric)) {
  stop("Need to define betaMetric")
}


# Other variables
outDir <- "data/diversity/"



# Loading dependencies ----------------------------------------------------

library(tidyverse)
library(furrr)
library(vegan)

# Loading common functions
source("code/R/metagenomeFunctions.R")



# Functions ---------------------------------------------------------------

# Function for calculating beta diversity using desired distance metric
# n_iter_int serves as a placeholder when iterating with map downstream
subsample_beta_diversity <- function(shared_df, combined_idxstats_df, n_reads_int, beta_metric_chr, n_iter_int=1) {
  
  # Subsample shared df
  sub_shared <- subsample_shared(shared_df, n_reads_int)
  
  # Normalize subsampled shared df
  norm_shared <- normalize_shared(sub_shared, combined_idxstats_df)

  # Transform shared df into matrix
  shared_mat <- norm_shared  %>%
    select(-numBins) %>% 
    column_to_rownames(var = "sample")
  
  # Calculating beta diversity distances
  beta_dist <- vegdist(shared_mat, method = beta_metric_chr)
  
  return(beta_dist)
  
}


# Function for calculating beta diversity by iteratively subsampling, computing beta diversity, and calculating the average over 1000 iterations
calc_beta_diversity <- function(shared_df, combined_idxstats_df, n_reads_int, beta_metric_chr) {
  
  # Iteratively calculating beta diversity while subsampling
  beta_list <- future_map(1:1000, subsample_beta_diversity, # Repeating subsampling, normalizing, and beta diversity calculation 1000 times
                         shared_df = filtered_shared, combined_idxstats_df = combined_idxstats_df, # Feeding in data files (it's a complex function)
                         n_reads_int = n_reads_int, beta_metric_chr = beta_metric_chr, # More inputs
                         .progress = TRUE, .options = future_options(seed = TRUE)) # Adding progress bar and making seed reproducible
  
  # Averaging distance matrices from all iterations
  beta_dist <- reduce(beta_list, `+`)/length(beta_list)
  
  # Formatting output for saving
  beta_mat <- as.matrix(beta_dist)
  
  return(beta_mat)
  
}



# Analysis ----------------------------------------------------------------

message("PROGRESS: Calculating beta diversity in parallel.")

# Creating output directory if it doesn't already exist
dir.create(outDir, recursive = TRUE, showWarnings=FALSE)

# Reading in data files
shared <- read_tsv(sharedFile, col_types = cols()) %>% 
  mutate(sample = as.character(sample))
idxstats <- read_tsv(idxstatsFile, col_types = cols()) %>% 
  mutate(sample = as.character(sample))

# Removing samples without enough reads for subsampling
filtered_shared <- filter_shared_by_reads(shared, subthresh)

# Setting resource plan for parallelizing calculations
plan(multiprocess)

# Setting seed to make analysis reproducible
set.seed(20170415)

# Calculating beta diversity
beta_div <- calc_beta_diversity(filtered_shared, idxstats, min(filtered_shared$numReads), betaMetric)

# Creating name for output file
out_file <- paste0(str_replace(str_extract(sharedFile, "[[:alpha:].]+shared"), "(^.+).raw.shared", "\\1.beta."), betaMetric, ".tsv")

message("PROGRESS: Writing beta diversity results to file.")

# Write out the beta diversity matrix
write.table(beta_div, file = paste0(outDir, out_file), sep = "\t", row.names = FALSE)
