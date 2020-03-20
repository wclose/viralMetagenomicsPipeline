#!/usr/bin/env Rscript
# metagenomeFunctions.R
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Repository for metagenome-related functions that are shared by more than one script.

# Loading dependencies ----------------------------------------------------

library(tidyverse)



# Functions ---------------------------------------------------------------

# Function for removing and reporting samples without enough reads before subsampling
filter_shared_by_reads <- function(shared_df, n_reads_int) {
  
  # Finding smallest number of sample reads above specified subsampling threshold
  subsample_reads <- shared_df %>% 
    filter(numReads >= n_reads_int) %>% 
    pull(numReads) %>% 
    min()
  
  # Pulling sample names WITHOUT enough reads for subsampling
  remove_samples <- shared_df %>% 
    filter(numReads < subsample_reads) %>% # Finding samples that don't have enough reads for subsample
    pull(sample)
  
  # Adding warning message if smallest number of reads is larger than the subsampling threshold
  # so the user knows the number of reads will be changed to maximize information.
  if (subsample_reads != n_reads_int) {
    
    message(paste("PROGRESS: Set subsamping threshold to", n_reads_int, "reads but able to increase to",
                  subsample_reads, "reads without losing additional samples.",
                  "Subsampling increased to", subsample_reads, "reads instead."))
    
  }
  
  # Adding warning message if samples are removed for having too low of counts
  if (length(remove_samples) > 0) {
    
    # Calculating percent of samples removed by subsampling and formatting to show two decimal points
    percent_removed <- paste0("(", format(round(length(remove_samples)/nrow(shared_df) * 100, 2), nsmall = 2), "%)")
    
    message(paste("PROGRESS: Removed", length(remove_samples), percent_removed, "samples for not having enough reads. These samples were removed:",
                  paste0(remove_samples, collapse = ", ")))
    
  }
  
  # Removing samples without enough reads from downstream analysis
  filtered_shared <- shared_df %>% 
    filter(!(sample %in% remove_samples))
  
  return(filtered_shared)
  
}


# Function for subsampling a specific row (sample) of shared file
subsample_shared_row <- function(shared_df, sample_chr, n_reads_int) {
  
  # Pulling out list of bins and sample counts for creating sampling vector
  bin_names <- str_subset(names(shared_df), "^Bin\\d+$") # Pulling bin names from col names
  sample_counts <- as.numeric(shared_df[shared_df$sample == sample_chr, bin_names]) # Pulling bin counts for one row at a time
  
  # Creating sampling vector
  sampling_vector <- rep(bin_names, sample_counts) # Repeats each bin name for as many read counts for that bin
  
  # Subsampling sampling vector to n_reads
  subsampled_row <- sample(x = sampling_vector, size = n_reads_int) %>% # Subsampling to desired read depth
    table() %>% # Counting instances of each bin
    enframe(name = "bin", value = sample_chr)
  
  return(subsampled_row)
  
}


# Function for subsampling entire shared file to desired read depth. Uses parallelization to speed up.
subsample_shared <- function(shared_df, n_reads_int) {
  
  # Subsampling each row and joining into final subsampled shared file
  subsampled_shared <- map(as.character(shared_df$sample), subsample_shared_row, shared_df = shared_df, n_reads_int = n_reads_int) %>% # Subsampling all samples
    reduce(full_join, by = "bin") %>% # Combines all of the subsampled sample dfs keepin all bins with counts from all samples
    right_join(tibble(bin = str_subset(names(shared_df), "^Bin\\d+$")), by = "bin") %>% # Joining to list of all bins to fill in missing bins
    gather(key = "sample", value = "count", -bin) %>% # Transforming to long df
    mutate(count = replace_na(count, 0), # Changing NAs from join to 0
           count = as.integer(count)) # Adding column of total bin counts

  return(subsampled_shared)
  
}


# Function for normalizing read counts by adjusting for bin size in bp
normalize_shared <- function(shared_df, combined_idxstats_df){
  
  # Normalizing shared file
  norm_shared <- shared_df %>% 
    left_join(select(combined_idxstats_df, sample, bin, bin_length), by = c("sample", "bin")) %>% # Joining bin length information to reformatted shared df
    mutate(norm_count = count/bin_length*1000) %>% # Calculating normalized read counts per kb of bin
    select(sample, bin, norm_count) %>% # Dropping unused cols
    mutate(numBins = length(unique(bin))) %>% # Adding col for number of bins
    spread(key = bin, value = norm_count) # Spreading df back into shared format
  
  return(norm_shared)
  
}
