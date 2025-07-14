suppressPackageStartupMessages({
  library(scds)
  library(Seurat)
  library(optparse)
})

option_list  <-  list(
  make_option(c("-b", "--base_path"), type="character", default=NULL, 
              help="Directorio base donde se guardarán los datos de salida. 
                    Ej: /export/storage/users/azaid/vallabh_lab/work", metavar="character"),
  make_option(c("-f", "--fig_dir"), type="character", default=NULL, 
              help="Directorio para guardar las figuras generadas. 
                    Ej: /home/azaid/vallabh_lab/work/figures", metavar="character"),
  make_option(c("-s", "--samples"), type="character", default=NULL, 
              help="IDs de las muestras a procesar, separadas por comas. 
                    Ej: SRR16502301,SRR16502302,SRR16502308", metavar="character")
)

opt_parser  <-  OptionParser(option_list=option_list)
opt  <-  parse_args(opt_parser)

base_path <- opt$base_path
fig_dir <- opt$fig_dir
samples <- unlist(strsplit(opt$samples, ","))

if (interactive()) {
  base_path <- "/export/storage/users/azaid/vallabh_lab/work"
  fig_dir <- '/home/azaid/vallabh_lab/work'
  samples <- c("SRR16502301", "SRR16502302")
}

for (s_tmp in samples) {
    # Load the preprocessed Seurat object
    seurat_obj <- readRDS(file.path(base_path, s_tmp, "soupx_matrix", paste0(s_tmp, "_soupx_seurat.rds")))
    
    # Convert to SingleCellExperiment for scds
    sce <- suppressWarnings({as.SingleCellExperiment(seurat_obj)})
    
    # Calculate doublet hybrid scores (scds)
    sce <- cxds_bcds_hybrid(sce)
    
    # Save the score as metadata in the Seurat object
    seurat_obj$scds_score <- sce$hybrid_score
    
    # Visualize the score distribution (optional but recommended)
    fig_path <- file.path(fig_dir, s_tmp, "soupx_matrix")
    dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)
    fig_path <- file.path(fig_path, "scds_score_hist.png")
    png(fig_path)
    hist(seurat_obj$scds_score, breaks = 50, main = paste("scds_score", s_tmp))
    abline(v = 1.3, col = "red")
    dev.off()
    
    # Filter cells directly using Seurat subset
    seurat_obj <- subset(seurat_obj, subset = scds_score < 1.3)
    
    # Save the filtered Seurat object
    dir.create(file.path(base_path, s_tmp, "scds_matrix"), recursive = TRUE, showWarnings = FALSE)
    saveRDS(seurat_obj, file.path(base_path, s_tmp, "scds_matrix", paste0(s_tmp, "_seurat.rds")))
    
    # Final summary of the filtering process
    cat(paste(
      "Sample:", s_tmp,
      "\nOriginal cells:", ncol(sce),
      "\nPost-filtering cells:", ncol(seurat_obj),
      "\nDeleted cells:", ncol(sce) - ncol(seurat_obj), "\n"
    ))
}