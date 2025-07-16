### --- Load libraries --- ###
library(Seurat)
library(ggplot2)
library(data.table)
library(harmony)
library(patchwork) 

### ----- Set up paths and parameters ----- ###
base_path <- "/export/storage/users/azaid/vallabh_lab/work"
samples <- c("SRR16502301", "SRR16502302")
fig_dir <- '/home/azaid/vallabh_lab/work'
processed_dir <- '/export/storage/users/azaid/aav_scrna/processed'
aav_vector <- c('AAV9', 'AAV9 2YF')
barcode_files <- list(
    B1_non_enriched = file.path(processed_dir, "GSM5641048_cell_by_aav_mouse_B1_brain_results.csv.gz"),
    B2_non_enriched = file.path(processed_dir, "GSM5641049_cell_by_aav_mouse_B2_brain_results.csv.gz"),
    B1_enriched = file.path(processed_dir, "GSM5641055_GFP_targeted_cell_by_aav_mouse_B1_brain_results.csv.gz"),
    B2_enriched = file.path(processed_dir, "GSM5641056_GFP_targeted_cell_by_aav_mouse_B2_brain_results.csv.gz")
)
key_files <- list(
    B1 = file.path(processed_dir, "GSM5641056_aav_barcode_key_mouse.txt.gz"),
    B2 = file.path(processed_dir, "GSM5641057_aav_barcode_key_mouse.txt.gz")
)
### -------- Quality control and filtering -------- ###
for (s_tmp in samples) {
    # Load Seurat object (doublets removed)
    seurat_obj <- readRDS(file.path(base_path, s_tmp, "scds_matrix", paste0(s_tmp, "_seurat.rds")))
    
    # Calculate mitochondrial content percentage
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
    
    # Create violin plot for each quality metric
    p1 <- VlnPlot(seurat_obj, features = "nCount_RNA", pt.size = 0.1) +
        geom_hline(yintercept = 500, linetype = "dashed", color = "red") +
        geom_hline(yintercept = 40000, linetype = "dashed", color = "red") +
        ggtitle(paste0(s_tmp, " - nCount_RNA"))
        NoLegend()

    p2 <- VlnPlot(seurat_obj, features = "nFeature_RNA", pt.size = 0.1) +
        geom_hline(yintercept = 500, linetype = "dashed", color = "red") +
        geom_hline(yintercept = 6500, linetype = "dashed", color = "red") +
        ggtitle(paste0(s_tmp, " - nFeature_RNA")) +
        NoLegend()

    p3 <- VlnPlot(seurat_obj, features = "percent.mt", pt.size = 0.1) +
        geom_hline(yintercept = 15, linetype = "dashed", color = "red") +
        ggtitle(paste0(s_tmp, " - percent.mt")) +
        NoLegend()

    p <- p1 + p2 + p3 + plot_layout(ncol = 3)

    # Ensure plot directory exists and save QC plot
    plot_dir <- file.path(fig_dir, s_tmp, "plots")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(plot_dir, paste0(s_tmp, "_qc_metrics.png")), p, width = 12, height = 5)
    
    # Filter cells based on QC thresholds (based on the violin plots)
    seurat_obj <- subset(seurat_obj, 
        subset = nCount_RNA > 500 & nCount_RNA < 40000 & 
        nFeature_RNA > 500 & nFeature_RNA < 6500 
        & percent.mt < 15)
    
    # Generate and save violin plots for filtered data
    p <- VlnPlot(seurat_obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)
    ggsave(file.path(plot_dir, paste0(s_tmp, "_qc_metrics_filtered.png")), p, width = 12, height = 5)
    
    # Ensure Seurat object directory exists and save filtered object
    srt_dir <- file.path(base_path, s_tmp, "seurat_obj")
    dir.create(srt_dir, recursive = TRUE, showWarnings = FALSE)
    
    saveRDS(seurat_obj, file.path(srt_dir, paste0(s_tmp, "_filtered_seurat.rds")))
}

### -------- Integration, Normalization, and Clustering -------- ###

# Load filtered Seurat objects
seurat_list <- lapply(samples, function(s_tmp) {
    readRDS(file.path(base_path, s_tmp, "seurat_obj", paste0(s_tmp, "_filtered_seurat.rds")))
})
names(seurat_list) <- samples

# Merge objects
combined_obj <- merge(seurat_list[[1]], seurat_list[2:length(seurat_list)], add.cell.ids = samples)
combined_obj <- JoinLayers(combined_obj, assay = "RNA")

# Normalize with SCTransform and perform PCA
combined_obj <- SCTransform(combined_obj, assay = "RNA", verbose = FALSE)
combined_obj <- RunPCA(combined_obj, npcs = 100)

# Run Harmony for batch effect correction
combined_obj <- RunHarmony(combined_obj, group.by.vars = "orig.ident", dims.use = 1:50)

### -------- AAV Barcode Annotation -------- ###
# Load and combine AAV barcode data
barcode_data <- lapply(names(barcode_files), function(s_tmp) {
    dt <- fread(barcode_files[[s_tmp]])
    # Asignar prefijo según muestra y enriquecimiento
    dt$sample <- ifelse(s_tmp %like% "B1", "SRR16502301", "SRR16502302")
    dt
    })
barcode_data <- rbindlist(barcode_data)
setnames(barcode_data, old = "V1", new = "barcode")
setnames(barcode_data, tolower(colnames(barcode_data)))

# Reformat barcodes to match combined object
barcode_data$barcode <- paste0(barcode_data$sample, "_", barcode_data$barcode)

# Load AAV variant keys
variant_keys <- lapply(names(key_files), function(s_tmp) {
    dt <- fread(key_files[[s_tmp]])
    dt$sample <- ifelse(s_tmp %like% "B1", "SRR16502301", "SRR16502302")
    dt
})
variant_keys <- rbindlist(variant_keys)
setnames(variant_keys, tolower(colnames(variant_keys)))
variant_keys$barcode <- as.character(variant_keys$barcode)

# Filter barcodes present in combined object
barcode_data <- barcode_data[barcode %in% colnames(combined_obj)]
# Create annotation columns
barcode_cols <- grep("^barcode[0-9]+$", colnames(barcode_data), value = TRUE)
barcode_data[,AAV_status:="AAV+"]
# Assign predominant AAV variant per cell
barcode_data[, AAV_variant := character()]
for (i in seq_along(barcode_data$barcode)) {
    sample_id <- barcode_data[i, sample]
    max_count <- max(barcode_data[i,..barcode_cols])
    barcode_name <- barcode_data[i,..barcode_cols] == max_count
    barcode_name <- barcode_cols[barcode_name]
    variant <- variant_keys[barcode %in% sub('^barcode','',barcode_name) & sample == sample_id, variant]
    if(any(!is.na(match(variant,aav_vector)))){
        variant <- aav_vector[as.integer(na.omit(match(variant,aav_vector)))][1]
    } else {
        variant <- variant[1]
    }
    barcode_data[i, AAV_variant := variant]
    
}

# Add annotations to Seurat object metadata
meta_data <- barcode_data[, .(barcode, AAV_status, AAV_variant)]
combined_obj$AAV_status <- meta_data[match(colnames(combined_obj), barcode), AAV_status]
combined_obj$AAV_status[is.na(combined_obj$AAV_status)] <- "AAV-"
combined_obj$AAV_variant <- meta_data[match(colnames(combined_obj), barcode), AAV_variant]

# Save annotated Seurat object
srt_dir <- file.path(base_path, "integrated", "seurat_obj")
dir.create(srt_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(combined_obj, file.path(srt_dir, "aav_annotated_filtered_seurat.rds"))