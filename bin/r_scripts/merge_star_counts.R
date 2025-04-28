#!/usr/bin/env Rscript

# merge_star_counts.R
# -------------------
# Scan a directory for STAR ReadsPerGene.out.tab files,
# build a count matrix, and write it out as TSV.

# Usage:
#   Rscript merge_star_counts.R [input_dir] [output_file]
# Defaults:
#   input_dir    = "."
#   output_file  = "gene_counts_matrix.tsv"

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir   <- ifelse(length(args) >= 1, args[1], ".")
output_file <- ifelse(length(args) >= 2, args[2], "gene_counts_matrix.tsv")

# Find all STAR count files
count_files <- list.files(
  path       = input_dir,
  pattern    = "ReadsPerGene\\.out\\.tab$",
  full.names = TRUE
)
if (length(count_files)==0) {
  stop("No ReadsPerGene.out.tab files found in: ", input_dir)
}

# Function to read one STAR count file
read_star_counts <- function(filepath) {
  # Derive sample name from filename
  sample <- basename(filepath) %>%
    sub("ReadsPerGene\\.out\\.tab$", "", .)
  
  # Read without header: cols are (1) geneID, (2) unstranded, (3) forward, (4) reverse
  df <- read_tsv(filepath, col_names = FALSE, comment = "#",
                 show_col_types = FALSE)
  
  tibble(
    gene_id = df[[1]],
    !!sample := df[[2]]    # change to df[[3]] or df[[4]] if stranded
  )
}

# Read all and join by gene_id
count_list <- map(count_files, read_star_counts)
counts_df  <- reduce(count_list, full_join, by = "gene_id")

# Optionally: Replace NAs with zero
counts_df[is.na(counts_df)] <- 0

# Write out
write_tsv(counts_df, output_file)
message("Wrote merged count matrix to: ", output_file)