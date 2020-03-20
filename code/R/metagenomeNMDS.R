#!/usr/bin/env Rscript
# metagenomeNMDS.R
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Calculate axes for principal coordinate analysis ordination of beta diversity.
# Usage: Rscript metagenomeNMDS.R BETAFILE

# Setting environment -----------------------------------------------------

# Parsing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Variables defined by user
betaFile <- args[1] # Beta diversity matrix

# Checking if inputs are set
if (!file.exists(betaFile)) {
  stop(paste(betaFile, "does not exist."))
} 


# Other variables
outDir <- "data/diversity/"



# Loading dependencies ----------------------------------------------------

library(tidyverse)
library(vegan)
library(ecodist)



# Functions ---------------------------------------------------------------

# Function for tidying configuration dfs
tidy_conf <- function(conf_list, sample_vec, iter_int) {
  
  # Tidying conf tables
  conf_df <- as_tibble(conf_list[[as.integer(iter_int)]]) %>% # Selecting conf table
    add_column(sample = sample_vec, .before = "V1") %>% # Adding col for sample names
    add_column(iter = as.integer(iter_int), .before = "V1") %>% # Adding col for iteration number
    rename_all(str_replace, "V", "axis") # Reformatting col names
  
  return(conf_df)
  
}



# Analysis ----------------------------------------------------------------

message("PROGRESS: Calculating NMDS ordination of samples.")

# Creating output directory if it doesn't already exist
dir.create(outDir, recursive = TRUE, showWarnings=FALSE)

# Reading in the full distance matrix and formatting as a dist object
beta <- read_tsv(betaFile, col_types = cols()) %>% 
  as.dist()

# Pulling sample names from dist object
sample_names <- colnames(as.matrix(beta))

# Setting random seed for reproducibility
set.seed(20170415)

# Computing the nmds for the given matrix using mothur defaults
nmds_list <- nmds(beta, mindim = 2, maxdim = 2, nits = 10,
                  epsilon = 1e-12, maxit = 500)

# Pulling together nmds stats into a single df
nmds_stats_df <- tibble(iter = 1:nmds_list$nits,
                        stress = nmds_list$stress,
                        r2 = nmds_list$r2)

# Tidying all of the conf df and combining together
nmds_ordination_df <- map_dfr(1:nmds_list$nits, tidy_conf,
                              conf_list = nmds_list$conf, sample_vec = sample_names)

# Formatting output file names
stats_output_file <- str_replace(str_extract(betaFile, "[[:alpha:].]+beta.\\w+.tsv"), "(^.+).tsv", "\\1.nmds.stats.tsv")
ordination_output_file <- str_replace(str_extract(betaFile, "[[:alpha:].]+beta.\\w+.tsv"), "(^.+).tsv", "\\1.nmds.axes.tsv")

message("PROGRESS: Writing NMDS results to file.")

# Write out the nmds stats and axes dfs
write_tsv(nmds_stats_df, path = paste0(outDir, stats_output_file))
write_tsv(nmds_ordination_df, path = paste0(outDir, ordination_output_file))
