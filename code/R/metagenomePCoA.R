#!/usr/bin/env Rscript
# metagenomePCoA.R
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Calculate axes for principal coordinate analysis ordination of beta diversity.
# Usage: Rscript metagenomePCoA.R BETAFILE

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



# Analysis ----------------------------------------------------------------

message("PROGRESS: Calculating PCoA ordination of samples.")

# Creating output directory if it doesn't already exist
dir.create(outDir, recursive = TRUE, showWarnings=FALSE)

# Reading in the full distance matrix and formatting as a dist object
beta <- read_tsv(betaFile, col_types = cols()) %>% 
  as.dist()

# Pulling sample names from dist object
sample_names <- colnames(as.matrix(beta))

# Computing the pcoa for the given distance matrix
pcoa_list <- pco(beta)

# Subsetting eigenvalues to only those above 0
pcoa_eigen_subset <- subset(pcoa_list$values, pcoa_list$values > 0)

# Calculating the amount of variation for each axes using eigenvalues
pcoa_stats_df <- enframe(pcoa_eigen_subset/sum(pcoa_eigen_subset)*100,
                         name = "axis", value = "loading")

# Combining axes coordinates with variation explained by each axes into a df for plotting
pcoa_ordination_df <- as_tibble(pcoa_list$vectors[,1:length(pcoa_eigen_subset)]) %>%
  add_column(sample = sample_names, .before = "V1") %>% 
  rename_all(str_replace, "V", "axis")

# Formatting output file names
stats_output_file <- str_replace(str_extract(betaFile, "[[:alpha:].]+beta.\\w+.tsv"), "(^.+).tsv", "\\1.pcoa.stats.tsv")
ordination_output_file <- str_replace(str_extract(betaFile, "[[:alpha:].]+beta.\\w+.tsv"), "(^.+).tsv", "\\1.pcoa.axes.tsv")

message("PROGRESS: Writing PCoA results to file.")

# Write out the nmds stats and axes dfs
write_tsv(pcoa_stats_df, path = paste0(outDir, stats_output_file))
write_tsv(pcoa_ordination_df, path = paste0(outDir, ordination_output_file))
