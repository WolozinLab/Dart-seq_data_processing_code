#!/usr/bin/env Rscript
#
# Merge and Filter Bullseye Sites
# This script merges multiple Bullseye output BED files and filters sites by DRACH motif
#

# Install required packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GenomicRanges", "BSgenome", "Biostrings", "gUtils"))

# Load required libraries
library(GenomicRanges)
library(BSgenome)
library(Biostrings)
library(gUtils)
library(dplyr)
library(purrr)
library(stringr)

# Function to merge and filter Bullseye sites
merge_filter_sites <- function(bed_dir, output_dir, genome_path) {
  # List BED files in directory
  fn <- list.files(bed_dir, pattern = "*.bed")
  
  if (length(fn) == 0) {
    stop("No BED files found in the specified directory")
  }
  
  # Read and process each BED file
  list_beds <- lapply(fn, function(i) {
    message("Processing file: ", i)
    
    bed_df <- read.table(file.path(bed_dir, i), 
                         sep = "\t", 
                         header = FALSE)
    
    # Filter out chromosomes starting with G or K
    bed_df <- bed_df[!grepl("^G|^K", bed_df$V1), ]
    
    # Set column names
    colnames(bed_df) <- c("chr", "start", "end", "gene", "score", "strand", 
                          "control_ratio", "control_total", "dart_ratio", 
                          "dart_total", "conversion")
    
    # Extract sample name from filename
    sample_name <- str_split_fixed(i, "_", 4)[, 3]
    
    # Add sample prefix to specific columns
    colnames(bed_df)[c(4, 5, 7, 8, 9, 10, 11)] <- paste0(
      sample_name, "_", colnames(bed_df)[c(4, 5, 7, 8, 9, 10, 11)]
    )
    
    # Create unique key for merging
    bed_df$key <- paste0(gsub("chr", "", bed_df$chr), "_", bed_df$start, "_",
                         bed_df$end, "_", bed_df$strand)
    
    return(bed_df)
  })
  
  # Merge all dataframes using full_join
  message("Merging all BED files...")
  merged_df <- reduce(list_beds, full_join, by = c("chr", "start", "end", "strand", "key"))
  
  # Save unfiltered merged results
  unfiltered_output <- file.path(output_dir, "Sum_of_all_no_DRACH_filter.csv")
  message("Writing unfiltered results to: ", unfiltered_output)
  write.table(merged_df, unfiltered_output, sep = "\t", quote = FALSE, 
              col.names = TRUE, row.names = FALSE)
  
  # Create GRanges object for sequence extraction
  gr <- GRanges(seqnames = merged_df$chr,
                ranges = IRanges(start = merged_df$start, end = merged_df$end),
                strand = merged_df$strand)
  
  # Adjust coordinates for BED format
  start(gr) <- start(gr) + 1
  
  # Export bed file for DRACH filtering
  bed_output <- file.path(output_dir, "unfilter_raw_bullseye.bed")
  message("Exporting BED file for DRACH filtering: ", bed_output)
  export(gr, bed_output)
  
  # Run DRACH filtering script
  message("Running DRACH filtering...")
  system(paste("bash", file.path(dirname(output_dir), "DRACH_filter.sh"), 
               bed_output, output_dir))
  
  # Read RAC filtered results
  rac_file <- file.path(output_dir, "RAC.bed")
  message("Reading RAC filtered results from: ", rac_file)
  RAC <- read.delim(rac_file, header = FALSE)
  
  # Filter by first nucleotide (A or G)
  RAC <- RAC[substr(RAC$V7, 1, 1) %in% c("A", "G"), ]
  
  # Create keys for filtering
  RAC$key <- paste0(gsub("chr", "", RAC$V1), "_", RAC$V2, "_",
                    RAC$V3, "_", RAC$V6)
  
  # Filter merged dataframe by RAC keys
  merged_df <- merged_df[merged_df$key %in% RAC$key, ]
  
  # Save DRACH filtered results
  filtered_output <- file.path(output_dir, "Sum_of_all_filter_by_DRACH.csv")
  message("Writing DRACH filtered results to: ", filtered_output)
  write.table(merged_df, filtered_output, sep = "\t", quote = FALSE, 
              col.names = TRUE, row.names = FALSE)
  
  # Extract editing ratio and gene information
  editing_ratio_df <- merged_df %>% 
    select(starts_with("chr"), 
           starts_with("start"), 
           starts_with("end"), 
           starts_with("strand"),
           grep("score", colnames(.), value = TRUE))
  
  count_df <- merged_df %>% 
    select(starts_with("chr"), 
           starts_with("start"), 
           starts_with("end"), 
           starts_with("strand"), 
           grep("gene", colnames(.), value = TRUE))
  
  # Save extracted data
  write.table(count_df, file.path(output_dir, "matrix_Col1.csv"), 
              quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
  
  write.table(editing_ratio_df, file.path(output_dir, "matrix_Col6.csv"), 
              quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
  
  message("Processing completed successfully!")
}

# Main execution
if (!interactive()) {
  # Get command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 3) {
    cat("Usage: Rscript merge_filter_sites.R <bed_directory> <output_directory> <genome_path>\n")
    cat("Example: Rscript merge_filter_sites.R /path/to/beds /path/to/output /path/to/genome.fa\n")
    quit(status = 1)
  }
  
  bed_dir <- args[1]
  output_dir <- args[2]
  genome_path <- args[3]
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  merge_filter_sites(bed_dir, output_dir, genome_path)
} 