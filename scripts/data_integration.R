# AAV scRNA-seq and snRNA-seq Integration and Annotation
# Description :
# 1. Exploratory Integration: Correct batch effects between the AAV (scRNA-seq)
#    and brain atlas (snRNA-seq) datasets for joint visualization.
# 2. PC Determination: Use data-driven diagnostics to select the optimal
#    number of Principal Components (PCs) for downstream analysis.
# 3. Label Transfer: Use the AAV dataset as a reference to predict the
#    transduction status ('AAV+' or 'AAV-') onto cells in the brain atlas.

### ------------------------ Load libraries ------------------------ ###
library(Seurat)
library(data.table)
library(ggplot2)
library(stringr)

### --------------------------- Load paths --------------------------- ###
base_path <- "/export/storage/users/azaid/vallabh_lab/work"
atlas_path <- "/export/storage/users/azaid/spatial"
output_path <- "/home/azaid/vallabh_lab/outs/integrated_analysis"
srt_out_path <- "/export/storage/users/azaid/vallabh_lab/work/integrated/seurat_obj"
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# Increase memory limit for large objects
options(future.globals.maxSize = 1000 * 1024^2)

### --------------- Data Loading and Preparation --------------- ###

# Load pre-processed Seurat objects
aav_obj <- readRDS(file.path(srt_out_path, "aav_annotated_filtered_seurat.rds"))
snrna_obj <- readRDS(file.path(atlas_path, "processed", "mouseBrain.snRNAseq.308Clusters.seurat.20241021.rds"))

# Load additional metadata to filter by brain region
cell_zones <- fread('/home/azaid/vallabh_lab/data/reference/cell_zone.csv')

# Ensure gene consistency between both objects
shared_genes <- intersect(rownames(aav_obj), rownames(snrna_obj))
aav_obj <- aav_obj[shared_genes, ]
snrna_obj <- snrna_obj[shared_genes, ]
gc()

# Filter the snRNA-seq atlas to keep only relevant regions (Cortex & Hippocampus)
snrna_obj$combined <- paste(snrna_obj$Cell_class, snrna_obj$Cell_subclass, sep="_")
cell_zones$combined <- paste(cell_zones$Cell_class, cell_zones$Cell_subclass, sep="_")
snrna_obj$Cell_zone <- cell_zones$Cell_zone[match(snrna_obj$combined, cell_zones$combined)]
snrna_obj$combined <- NULL

snrna_obj_filtered <- subset(snrna_obj, subset = Cell_zone %in% c("Cortex", "Hippocampus", "Cortex and hippocampus"))
rm(snrna_obj); gc()

### --------------------- Exploratory Integration --------------------- ###

# ---> Merge datasets, determine the optimal number of PCs, and create a joint UMAP to visualize the alignment of cells.

# Merge objects
snrna_obj_filtered$experiment <- "snRNAseq_Atlas"
aav_obj$experiment <- "AAV_scRNAseq"
combined_obj <- merge(snrna_obj_filtered, aav_obj, add.cell.ids = c("snRNAseq_Atlas", "AAV_scRNAseq"))

DefaultAssay(combined_obj) <- "RNA"
combined_obj <- NormalizeData(combined_obj, verbose = FALSE)
combined_obj <- FindVariableFeatures(combined_obj, verbose = FALSE)
combined_obj <- ScaleData(combined_obj, verbose = FALSE)

# Run PCA with number of components to evaluate
combined_obj <- RunPCA(combined_obj, npcs = 100, verbose = FALSE)

# PC Selection
pc_output_path <- file.path(output_path, "pc_diagnostics")
dir.create(pc_output_path, recursive = TRUE, showWarnings = FALSE)

elbow_plot <- ElbowPlot(combined_obj, ndims = 100)
ggsave(file.path(pc_output_path, "pc_elbow_plot.pdf"), plot = elbow_plot, width = 8, height = 6)

pdf(file.path(pc_output_path, "pc_heatmaps.pdf"), width = 12, height = 16)
DimHeatmap(combined_obj, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(combined_obj, dims = 40:55, cells = 500, balanced = TRUE)
DimHeatmap(combined_obj, dims = 85:100, cells = 500, balanced = TRUE)
dev.off()

# ---> Manual step
# Inspect the 'pc_elbow_plot.pdf' and 'pc_heatmaps.pdf'
# Based on inspection, choose the number of PCs to use for downstream analysis
# A value of 50 is a reasonable starting point
n_pcs_to_use <- 50

# Run Integration and Visualization
combined_obj <- IntegrateLayers(
  object = combined_obj,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  dims = 1:n_pcs_to_use,
  verbose = FALSE
)

combined_obj <- FindNeighbors(combined_obj, reduction = "integrated.cca", dims = 1:n_pcs_to_use)
combined_obj <- RunUMAP(combined_obj, reduction = "integrated.cca", dims = 1:n_pcs_to_use)

# Save Visualizations
p1 <- DimPlot(combined_obj, reduction = "umap", group.by = "experiment", pt.size = 0.5) +
  ggtitle("Integrated UMAP by Experiment")
ggsave(file.path(output_path, "exploratory_umap_by_experiment.pdf"), plot = p1, width = 8, height = 6)

p2 <- DimPlot(combined_obj, reduction = "umap", group.by = "AAV_status", pt.size = 0.5) +
  ggtitle("Integrated UMAP by AAV Status")
ggsave(file.path(output_path, "exploratory_umap_by_aav_status.pdf"), plot = p2, width = 8, height = 6)

# Save the integrated object
saveRDS(combined_obj, file = file.path(srt_out_path, "full_integrated_seurat.rds"))

# Clean up memory
rm(combined_obj); gc()

### ------------------------ Transfer Annotation ------------------------ ###
# ---> Predict 'AAV_status' onto the snRNA-seq atlas.

# Set default assay and process each object independently
DefaultAssay(aav_obj) <- "RNA"
DefaultAssay(snrna_obj_filtered) <- "RNA"

# Process each object
aav_obj <- NormalizeData(aav_obj, verbose = FALSE)
aav_obj <- FindVariableFeatures(aav_obj, verbose = FALSE)
aav_obj <- ScaleData(aav_obj, verbose = FALSE)
aav_obj <- RunPCA(aav_obj, npcs = n_pcs_to_use, verbose = FALSE)

snrna_obj_filtered <- NormalizeData(snrna_obj_filtered, verbose = FALSE)
snrna_obj_filtered <- FindVariableFeatures(snrna_obj_filtered, verbose = FALSE)
snrna_obj_filtered <- ScaleData(snrna_obj_filtered, verbose = FALSE)
snrna_obj_filtered <- RunPCA(snrna_obj_filtered, npcs = n_pcs_to_use, verbose = FALSE)
snrna_obj_filtered <- RunUMAP(snrna_obj_filtered, dims = 1:n_pcs_to_use)

# Find transfer anchors and transfer labels
transfer_anchors <- FindTransferAnchors(
  reference = aav_obj,
  query = snrna_obj_filtered,
  dims = 1:n_pcs_to_use,
  reference.reduction = "pca",
  verbose = FALSE
)

predictions <- TransferData(
  anchorset = transfer_anchors,
  refdata = aav_obj$AAV_status,
  dims = 1:n_pcs_to_use,
  verbose = FALSE
)

snrna_obj_filtered <- AddMetaData(snrna_obj_filtered, metadata = predictions)

# Save the final annotated object
saveRDS(snrna_obj_filtered, file = file.path(srt_out_path, "snrna_annotated_aav_status_seurat.rds"))

cat("Analysis complete. Both the exploratory integrated object and the final annotated atlas have been saved.\n")