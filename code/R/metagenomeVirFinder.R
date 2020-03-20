#!/usr/bin/env Rscript
# metagenomeVirFinder.R
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Use VirFinder to predict viral contigs at given false-discovery and false-positive rates.
# Usage: Rscript metagenomeVirFinder.R CONTIGFASTA MODELFILE FPRTHRESHOLD FDRTHRESHOLD

# Setting environment -----------------------------------------------------

# Parsing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Variables defined by user
contigFasta <- args[1] # Location of contig library file
modelFile <- args[2] # Location of eukaryotic/prokaryotic model file
fprThreshold <- suppressWarnings(as.numeric(args[3])) # False positive rate (fp/(fp+tn)) threshold
fdrThreshold <- suppressWarnings(as.numeric(args[4])) # False discovery rate (fp/(fp+tp)) threshold

# Checking if inputs are set
if (!file.exists(contigFasta)) {
  stop(paste(contigFasta, "does not exist."))
} else if (!file.exists(modelFile)) {
  stop(paste(modelFile, "does not exist."))
} else if (is.na(fprThreshold)) {
  stop("Input fprThreshold must be defined and numeric.")
} else if (is.na(fdrThreshold)) {
  stop("Input fdrThreshold must be defined and numeric.")
}


# Other variables
outDir <- "data/virfinder/"




# Loading dependencies ----------------------------------------------------

library(VirFinder)
library(Biostrings)
library(furrr)
library(tidyverse)



# Functions ---------------------------------------------------------------

# Functions here are for parallelizing VirFinder alogrithm using futures,
# results will be the same as running VF.pred.user().
# Based in part on code from R. Eric Collins (Github: rec3141) found here:
# https://github.com/rec3141/VirFinder/blob/master/linux/VirFinder/R/parVF.pred.R

# Function for predicting viral contigs one sequence at a time
predict_VF_single <- function (string_name, string_seq) {
  
  # Parsing sequence into character vector
  seqFa <- strsplit(x=as.character(string_seq),split="",fixed=T)[[1]]
  
  # Counting k-mer frequences (model features)
  featureOut <- countSeqFeatureCpp(seqFa, modEPV)[[1]]
  
  # Counting total sequence length
  seqLength <- length(seqFa)
  
  # Determining model parameters based on sequence length
  if (seqLength < 1 * 1000) {
    lasso.mod <- attr(modEPV, "lasso.mod_0.5k")
    nullDis <- attr(modEPV, "nullDis_0.5k")
  } else if (seqLength < 3 * 1000) {
    lasso.mod <- attr(modEPV, "lasso.mod_1k")
    nullDis <- attr(modEPV, "nullDis_1k")
  } else {
    lasso.mod <- attr(modEPV, "lasso.mod_3k")
    nullDis <- attr(modEPV, "nullDis_3k")
  }
  
  # Running model on data
  lasso.pred <- predict(lasso.mod, t(as.matrix(featureOut)), type = "response", s = "lambda.min")
  
  # Determining significance of prediction
  pvalue <- mean(nullDis > as.numeric(lasso.pred))
  
  # Printing status message
  print(paste(string_name,
              "len", seqLength,
              "score", round(lasso.pred, 4), 
              "pvalue", round(pvalue, 4)))
  
  # Formatting results as df
  predResult <- data.frame(name = as.character(string_name),
                           length = seqLength,
                           score = unname(lasso.pred), # Need to remove default name
                           pvalue = pvalue)
  
  # Converting names from factor to character for downstream processing
  predResult$name <- as.character(predResult$name)
  
  return(predResult)
  
}


# Function for predicting viral contigs from a dataframe of contig names and sequences
predict_VF <- function(string_df) {
  
  # Parallelizing prediction of each sequence in fasta file
  predResult <- map2_dfr(.x = string_df$name, .y = string_df$seq, .f = predict_VF_single)
  
  return(predResult)
  
}


# Function for predicting viral contigs in parallel using VirFinder
predict_VF_parallel <- function(seqFaIn_chr, cores=1) {
  
  # Status message
  message("PROGRESS: Reading contig sequences from fasta file.")
  
  # Reading the fasta file
  string_list <- readDNAStringSet(seqFaIn_chr)
  
  # Creating dataframe listing one group per core and number of contigs per group for creating group_no col later
  group_count <- tibble(group_no = 1:cores, count = floor(length(string_list)/cores))
  
  # If the number of contigs is not wholly divisible by the number of cores, adjust the group counts by 1 until it is
  if ((length(string_list) %% cores) != 0) {
    
    group_count <- group_count %>%
      mutate(count = case_when(group_no %in% seq(1:(length(string_list) %% cores)) ~ count + 1,
                               TRUE ~ count))
    
  }
  
  # Reformatting the sequence data to be in a df, divided into one group per core, and nested for mapping over
  string_df <- tibble(name = names(string_list), seq = as.character(string_list)) %>%
    mutate(group_no = rep(group_count$group_no, group_count$count)) %>% # Creating a grouping index for one group per core
    group_by(group_no) %>% # Grouping by new index
    group_nest() # Nesting groups for mapping over
  
  # Status message
  message("PROGRESS: Predicting viral contigs with VirFinder.")

  # Predicting viral contigs using parallelization
  predResult <- future_map_dfr(string_df$data, predict_VF, # Mapping over nested dfs
                               .progress = TRUE, .options = future_options(seed = TRUE)) %>% # Adding a progress bar and making seed reproducible
    as_tibble() # Converting output to tibble for final processing steps
  
  return(predResult)
  
}



# Analysis ----------------------------------------------------------------

# Creating output directory if it doesn't already exist
dir.create(outDir, recursive = TRUE, showWarnings=FALSE)

# Loading the prokaryotic/eukaryotic prediction model
load(modelFile)

# Determining the number of available cores for parallel processing
cores <- availableCores()

# Defining the type of parallelization to be used
plan(multiprocess)

# Predicting contigs using the model and calculating qvalues (fdr)
contig_pred <- predict_VF_parallel(contigFasta, cores) %>% # Predicting contigs
  mutate(qvalue = VF.qvalue(pvalue)) # Calculating qvalue (fdr)

# Filtering results based on fpr/fdr thresholds
filtered_contig_pred <- contig_pred %>% 
  filter(pvalue <= fprThreshold & qvalue <= fdrThreshold) # Filtering out contigs with too high fpr and/or fdr

message("PROGRESS: Writing VirFinder results to file.")

# Writing raw results to output tsv file
results_file <- paste0(outDir, "virfinder_results.tsv")
write_tsv(contig_pred, results_file)

# Writing results after quality filtration to output tsv file
filtered_results_file <- paste0(outDir, "virfinder_fpr", fprThreshold, "_fdr", fdrThreshold, "_results.tsv")
write_tsv(filtered_contig_pred, filtered_results_file)

# Writing list of predicted viral contigs to output txt file
name_file <- paste0(outDir, "virfinder_predicted_contigs.txt")
write_tsv(enframe(filtered_contig_pred$name, name = NULL), name_file, col_names = FALSE)
