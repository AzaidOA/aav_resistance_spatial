#!/usr/bin/env Rscript

# ---- Load libraries ----
suppressPackageStartupMessages({
  library(argparse)
  library(DropletUtils)
  library(Matrix)
})

# ---- Argument parsing with argparse ----
parser <- ArgumentParser(
  description = 'Filter scRNA-seq droplets using DropletUtils defaultDrops with customizable paths and samples.'
)
parser$add_argument('-i', '--input_base', required=TRUE, help='Input base directory (where your STARsolo results are).')
parser$add_argument('-o', '--output_base', required=TRUE, help='Output base directory (where filtered matrices will be saved).')
parser$add_argument('-s', '--samples', required=TRUE, help='Comma-separated list of sample IDs (e.g. SRR12345,SRR67890).')
parser$add_argument('--lower_prop', type='double', default=0.5, help='Proportion of lowest-UMI barcodes used as empty (default: 0.1).')
parser$add_argument('--seed', type='integer', default=100, help='Random seed for reproducibility (default: 100).')

args <- parser$parse_args()

input_base <- args$input_base
output_base <- args$output_base
samples <- strsplit(args$samples, ",")[[1]]
lower_prop <- args$lower_prop
seed <- args$seed

if (interactive()) {
  # --- TEST VALUES FOR INTERACTIVE USE ---
  input_base <- "/home/azaid/vallabh_lab/work"
  output_base <- "/export/storage/users/azaid/vallabh_lab/work"
  samples <- c("SRR16502301", "SRR16502302", "SRR16502308", "SRR16502309")
  lower_prop <- 0.05
  seed <- 100
}
# ---- Function to process a sample ----
filter_droplets <- function(s_tmp, input_base, output_base, lower_prop, seed) {
  input_dir <- file.path(input_base, s_tmp, "starsolo", "Solo.out", "Gene", "raw")
  output_dir <- file.path(output_base, s_tmp, "filtered_matrix")
  
  sce <- read10xCounts(input_dir)
  
  if (sum(counts(sce)) == 0) {
    stop(paste("Error: No counts in matrix for sample", s_tmp))
  }
  
  set.seed(seed)
  called <- defaultDrops(counts(sce), lower.prop=0.05)
  
  sce <- sce[, called]
  
  write10xCounts(
    path = paste0("/export/storage/users/azaid/vallabh_lab/work/", s_tmp, "/filtered_matrix"),
    x = counts(sce),
    gene.id = rowData(sce)$ID,
    gene.symbol = rowData(sce)$Symbol,
    barcodes = colData(sce)$Barcode, # Usar códigos de barras correctos
    version = "3", # Formato v3
    overwrite = TRUE
  )
}

# ---- Main loop ----
for (s_tmp in samples) {
  cat("Processing sample:", s_tmp, "\n")
  filter_droplets(s_tmp, input_base, output_base, lower_prop, seed)
  cat("Filtering completed for sample:", s_tmp, "\n")
}
