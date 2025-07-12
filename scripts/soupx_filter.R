#!/usr/bin/env Rscript

# Load required packages
suppressPackageStartupMessages({
  library(argparse)
  library(SoupX)
  library(DropletUtils)
  library(Seurat)
})

# Argument parsing
parser <- ArgumentParser(description = "Run SoupX ambient RNA filtering per sample")
parser$add_argument("--base_path", required = TRUE, help = "Base path to sample directories")
parser$add_argument("--samples", required = TRUE, nargs = "+", help = "List of sample IDs to process")
args <- parser$parse_args()

base_path <- args$base_path
samples <- args$samples

if (interactive()) {
  # Test values for interactive use
  base_path <- "/home/azaid/vallabh_lab/work"
  samples <- c("SRR16502301", "SRR16502302", "SRR16502308", "SRR16502309")
}

# Function to process each sample
process_sample <- function(s_tmp) {
  # Paths for input and output
  filtered_path <- file.path(base_path, s_tmp, "filtered_matrix")
  raw_path <- file.path(base_path, s_tmp, "starsolo/Solo.out/Gene/raw")
  output_path <- file.path(base_path, s_tmp, "soupx_matrix")
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  
  # Load filtered and raw data
  toc <- Read10X(data.dir = filtered_path)
  tod <- Read10X(data.dir = raw_path)
  
  # Create SoupChannel object
  sc <- SoupChannel(tod, toc, calcSoupProfile = TRUE)
  
  # Manually set contamination fraction (default 0.1)
  sc <- setContaminationFraction(sc, 0.1)
  
  # Adjust counts to remove ambient RNA
  out <- adjustCounts(sc, roundToInt = FALSE)

  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = out, project = s_tmp, assay = "RNA")

  # Save Seurat object as RDS
  saveRDS(seurat_obj, file = file.path(output_path, paste0(s_tmp, "_soupx_seurat.rds")))
}

# Process each sample
for (s_tmp in samples) {
  cat("Processing sample:", s_tmp, "\n")
  process_sample(s_tmp)
  cat("Sample processed:", s_tmp, "\n")
}
