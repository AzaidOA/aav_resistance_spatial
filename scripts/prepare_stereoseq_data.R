# -----------------------------------------------------------------------------
# Script to Load and Assemble ALL Stereo-seq Data (Both Mice)
# Strategy: Load each slice as an individual Seurat object and store them in a list.
# -----------------------------------------------------------------------------

# ---- 1. Environment Setup ----
library(Seurat)
library(data.table)
library(Matrix)
library(stringr)
library(ggplot2)

# ---- 2. Define Input and Output Paths ----
# Main path containing the 'mouse1' and 'mouse2' subdirectories
base_data_path <- "/export/storage/users/azaid/spatial"
# Path to the master annotation file
annotations_file <- file.path(base_data_path, "processed", "stereoseq.celltypeTransfer.2mice.all.tsv.gz")
# Output path for the final object
output_path <- "/export/storage/users/azaid/vallabh_lab/work/spatial"
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# ---- 3. Define a Function to Process a Single Slice ----
# Encapsulating the logic in a function makes the code clean and reusable.
process_stereo_slice <- function(expression_file_path) {
  
  # Extract identifiers from the filename for tracking
  slice_name <- str_extract(basename(expression_file_path), "T\\d+")
  mouse_name <- str_extract(basename(expression_file_path), "mouse\\d+")
  cat(paste0("Processing: ", mouse_name, " - ", slice_name, "\n"))
  
  # Load expression data using the fast fread function
  spatial_data <- fread(expression_file_path)
  
  # Aggregate UMI counts per gene for each segmented cell ('cell_label')
  # Filter out 'cell_label' 0, which represents background noise.
  cell_counts <- spatial_data[cell_label!= 0,.(total_umi = sum(umi_count)), by =.(gene, cell_label)]
  
  # Create the sparse count matrix (genes x cells format)
  counts_matrix <- dcast(cell_counts, gene ~ cell_label, value.var = "total_umi", fun.aggregate = sum, fill = 0)
  gene_names <- counts_matrix$gene
  counts_matrix <- as.matrix(counts_matrix[, -1])
  rownames(counts_matrix) <- gene_names
  counts_sparse <- as(counts_matrix, "sparseMatrix")
  
  # Create the Seurat object
  seurat_obj <- CreateSeuratObject(counts = counts_sparse, project = paste0(mouse_name, "_", slice_name))
  
  # Calculate cell centroids and add the spatial information
  cell_centroids <- spatial_data[cell_label!= 0,.(x = mean(x), y = mean(y)), by =.(cell_label)]
  coords_matrix <- as.matrix(cell_centroids[,.(x, y)])
  rownames(coords_matrix) <- cell_centroids$cell_label
  coords_matrix <- coords_matrix[colnames(seurat_obj), ]
  
  # Paso 1: Crear un objeto 'Centroids' a partir de la matriz de coordenadas.
  # Este es el nuevo formato que Seurat v5 espera.
  centroids_obj <- CreateCentroids(coords = coords_matrix)

  # Paso 2: Crear el objeto FOV usando el objeto 'Centroids' que acabamos de hacer.
  # Nota que ya no necesitamos el argumento 'type = "centroids"' porque el tipo de objeto ya lo especifica.
  fov <- CreateFOV(
    coords = centroids_obj,
    assay = "RNA"
  )


  seurat_obj[["slice"]] <- fov # Store the FOV in an image slot named 'slice'
  
  return(seurat_obj)
}

# ---- 4. Find and Process All Expression Files ----

cat("Finding all Stereo-seq expression files...\n")
# Use a regular expression to find all 'total_gene...' files in the subdirectories
all_expression_files <- list.files(
  path = base_data_path,
  pattern = "^total_gene_.*\\.txt\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)

cat(paste0("Found ", length(all_expression_files), " slice files. Processing now...\n"))

# Apply the processing function to every file and store the results in a list.
# This is the most time-consuming step.
stereo_list <- lapply(all_expression_files, process_stereo_slice)
gc()

# Assign informative names to each element in the list for easy access later
names(stereo_list) <- str_extract(basename(all_expression_files), "T\\d+_mouse_f001_2D_mouse\\d+")

# ---- 5. Add Cell Type Annotations to All Objects in the List ----

cat("Loading the master annotation file...\n")
# Load the large annotation file only once for efficiency
master_annotations <- fread(annotations_file)

cat("Adding annotations to each Seurat object in the list...\n")
for (slice_id in names(stereo_list)) {
  
  # Extract the mouse name and slice ID from the list element's name
  current_mouse <- str_extract(slice_id, "mouse\\d+")
  current_t_id <- str_extract(slice_id, "T\\d+")
  
  cat(paste0("Annotating: ", current_mouse, " - ", current_t_id, "\n"))
  
  # Filter the master annotation table to get annotations for the current slice
  slice_annotations <- master_annotations[mouse == current_mouse & section_id == current_t_id]
  
  # Prepare data for merging
  slice_annotations[, cell_id := as.character(cell_id)]
  setkey(slice_annotations, cell_id)
  
  # Get current metadata and merge with the new annotations
  current_metadata <- stereo_list[[slice_id]]@meta.data
  current_metadata$cell_id <- rownames(current_metadata)
  
  new_metadata <- merge(current_metadata, slice_annotations, by = "cell_id", all.x = TRUE)
  rownames(new_metadata) <- new_metadata$cell_id
  
  # Add the enriched metadata back to the object in the list
  stereo_list[[slice_id]] <- AddMetaData(stereo_list[[slice_id]], metadata = new_metadata)

  # Contar cuántas células tienen anotaciones NA antes de filtrar
  cells_before_qc <- ncol(stereo_list[[slice_id]])
  cat(paste0("Cells before QC filtering: ", cells_before_qc, "\n"))

  # Filtrar el objeto Seurat para mantener solo las células que tienen una anotación (es decir, no son NA)
  # Usamos la columna 'cell_class' como referencia. Si es NA, la célula fue filtrada por los autores.
  stereo_list[[slice_id]] <- subset(stereo_list[[slice_id]], subset =!is.na(cell_class))

  cells_after_qc <- ncol(stereo_list[[slice_id]])
  cat(paste0("Cells after QC filtering: ", cells_after_qc, " (", cells_before_qc - cells_after_qc, " removed)\n"))
}

# ---- 6. Verification and Final Save ----

# This section uses a robust ggplot2-based method to generate a verification plot
# for a single example slice, confirming the success of our data processing.

cat("Generating verification plot for an example slice...\n")

# Select an example slice to plot (e.g., the first one in the list).
example_slice_name <- names(stereo_list)[1]
seurat_obj_example <- stereo_list[[example_slice_name]]

# --- Robust Plotting with ggplot2 ---

# 1. Extract the necessary data: spatial coordinates and cell metadata.
#    We access the FOV (Field of View) object, which is named 'slice' in your objects.
coords_df <- GetTissueCoordinates(seurat_obj_example[['slice']])
metadata_df <- seurat_obj_example@meta.data

# 2. Combine coordinates and metadata into a single data frame for plotting.
coords_df$cell_id <- rownames(coords_df)
metadata_df$cell_id <- rownames(metadata_df)
plot_data <- merge(metadata_df, coords_df, by = "cell_id")

# 3. Create the plot using ggplot2's fundamental functions.
#    This approach is more stable and bypasses the internal issues of SpatialDimPlot.
p <- ggplot(plot_data, aes(x = x, y = y, color = cell_class)) +
  geom_point(size = 0.5) + # Use small points for high-density data
  coord_fixed() +         # Maintain the correct aspect ratio of the tissue
  ggtitle(paste("Cell Type Distribution in Slice:", example_slice_name)) +
  guides(color = guide_legend(override.aes = list(size = 3), title = "Cell Class")) +
  theme_void() +          # Use a clean theme without axes for spatial plots
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "right"
  )

# 4. Save the generated plot to the specified output file.
output_plot_path <- file.path(output_path, paste0("plot_", example_slice_name, ".pdf"))
ggsave(filename = output_plot_path, plot = p, width = 12, height = 10, dpi = 300)
cat(paste("Verification plot saved to:", output_plot_path, "\n"))


# ---- 7. Save the Final Processed Object ----

# Save the complete list of processed Seurat objects to an.rds file.
output_rds_path <- file.path(output_path, "stereo_seq_full_dataset_list.rds")
cat(paste("Saving the final list of", length(stereo_list), "processed Seurat objects...\n"))
saveRDS(stereo_list, file = output_rds_path)
cat(paste("Final object saved successfully to:", output_rds_path, "\n"))
cat("Script finished.\n")