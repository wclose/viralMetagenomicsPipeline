#!/usr/bin/env Rscript
# metagenomeAlpha.R
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Calculate alpha diversity by iteratively subsampling normalized shared files.
# Usage: Rscript metagenomeAlpha.R SHAREDFILE IDXSTATS SUBTHRESH ALPHAMETRIC1 ALPHAMETRIC2 ALPHAMETRIC3 ... ALPHAMETRICN

# Setting environment -----------------------------------------------------

# Parsing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Variables defined by user
sharedFile <- args[1] # Raw metagenomic shared file
idxstatsFile <- args[2] # Combined idxstats file
subthresh <- suppressWarnings(as.integer(args[3])) # Minimum number of reads to subsample to
alphaMetrics <- args[4:length(args)] # Alpha diversity metrics to use. Possible options: "shannon", "simpson", "invsimpson", "sobs"

# Checking if inputs are set
if (!file.exists(sharedFile)) {
  stop(paste(sharedFile, "does not exist."))
} else if (!file.exists(idxstatsFile)) {
  stop(paste(idxstatsFile, "does not exist."))
} else if (is.na(subthresh)) {
  stop("Input subthresh must be defined and an integer")
} else if (is.na(alphaMetrics)) {
  stop("Need to define alphaMetrics")
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

# Function for calculating alpha diversity using a single alpha diversity metric
calc_alpha_diversity_metric <- function(shared_df, alpha_metric_chr) {
  
  # Transform shared df into matrix
  shared_mat <- shared_df  %>%
    select(-numBins) %>% 
    column_to_rownames(var = "sample")
  
  # Setting which functions to use for which alpha metrics
  if (alpha_metric_chr %in% c("shannon", "simpson", "invsimpson")) {
    
    # Calculating shannon, simpson, and/or inverse simpson indices
    alpha_vec <- diversity(shared_mat, index = alpha_metric_chr)
  
  } else if (alpha_metric_chr == "sobs") {
    
    # Calculating number of observed otus
    alpha_vec <- rowSums(shared_mat != 0)
    
  }

  # Transforming into df
  alpha_df <- tibble(sample = names(alpha_vec),
                     !!alpha_metric_chr := alpha_vec)
  
  return(alpha_df)
  
}


# Function for calculating alpha diversity using all supplied alpha diversity metrics
# n_iter_int serves as a placeholder when iterating with map downstream
subsample_alpha_diversity <- function(shared_df, combined_idxstats_df, n_reads_int, alpha_metric_vec, n_iter_int=1) {
  
  # Subsample shared df
  sub_shared <- subsample_shared(shared_df, min(shared_df$numReads)) %>% 
    normalize_shared(combined_idxstats_df)

  # Calc alpha diversity
  alpha_df <- map(alpha_metric_vec, calc_alpha_diversity_metric, shared_df = sub_shared) %>% # Calculate values of each diversity metric
    reduce(full_join, by = "sample") %>% # Join all diversity metric dfs into a single df
    mutate(n_iter = n_iter_int) %>% # Add a column representing the iteration number (used when iteratively sampling)
    gather(key = "metric", value = "value", -sample, -n_iter) # Making it look nice
  
  return(alpha_df)
  
}


# Function for calculating alpha diversity by iteratively subsampling, computing alpha diversity, and calculating the average over 1000 iterations
calc_alpha_diversity <- function(shared_df, combined_idxstats_df, n_reads_int, alpha_metric_vec) {

  # Iteratively calculating alpha diversity while subsampling
  alpha_df <- future_map_dfr(1:1000, subsample_alpha_diversity, # Repeating subsampling, normalizing, and alpha diversity calculation 1000 times
                             shared_df = shared_df, combined_idxstats_df = combined_idxstats_df, # Feeding in data files (it's a complex function)
                             n_reads_int = n_reads_int, alpha_metric_vec = alpha_metric_vec, # More inputs
                             .progress = TRUE, .options = future_options(seed = TRUE)) %>% # Adding progress bar and making seed reproducible
    group_by(sample, metric) %>% # Setting up to calculate average of metrics across all iterations
    summarize(ave = mean(value)) %>% # Calculating final averages
    ungroup() %>%
    mutate(ave = trimws(format(round(ave, 6), nsmall = 6))) %>%
    spread(key = metric, value = ave) # Making it look nice

  return(alpha_df)

}



# Analysis ----------------------------------------------------------------

message("PROGRESS: Calculating alpha diversity in parallel.")

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

# Calculating alpha diversity
alpha_div <- calc_alpha_diversity(filtered_shared, idxstats, min(filtered_shared$numReads), alphaMetrics)

# Creating name for output file
out_file <- str_replace(str_extract(sharedFile, "[[:alpha:].]+shared"), "(^.+).raw.shared", "\\1.alpha.tsv")

message("PROGRESS: Writing alpha diversity results to file.")

# Write out the alpha diversity df
write_tsv(alpha_div, path = paste0(outDir, out_file))
