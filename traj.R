#!/usr/bin/env Rscript
# Monocle3 Trajectory Analysis from Raw 10x Matrix Files
# Based on pipeline.py preprocessing workflow

library(monocle3)
library(Seurat)
library(dplyr)
library(Matrix)

# Set working directory
# setwd("/home/bjx143_pkuhpc/WL-HW1")

# Create output directory
dir.create("figures", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

# %%
# 1. Load datasets from raw 10x matrix files
print("Loading raw 10x matrix files...")

data_path <- "data"
subdata_dirs <- c(
  "GSM6304413_EB_D0",
  "GSM6304414_iMK_D3", 
  "GSM6304415_iMK_D5",
  "GSM6304416_iMK_D7"
)
titles <- c("EB_D0", "iMK_D3", "iMK_D5", "iMK_D7")

# Load each dataset
adatas <- list()
for (i in seq_along(subdata_dirs)) {
  subdir <- subdata_dirs[i]
  print(paste("Loading", subdir, "..."))
  
  # Read 10x data using Seurat's Read10X function
  data <- Read10X(data.dir = file.path(data_path, subdir))
  
  # Create Seurat object
  adata <- CreateSeuratObject(counts = data, project = subdir, min.cells = 0, min.features = 0)
  
  # Add sample and timepoint information
  adata$sample <- subdir
  adata$timepoint <- strsplit(subdir, "_")[[1]][3]  # D0, D3, D5, D7
  
  adatas[[i]] <- adata
  print(paste("Loaded", subdir, ":", ncol(adata), "cells,", nrow(adata), "genes"))
}

# %%
# 2. Merge datasets BEFORE QC (same as pipeline.py)
print("Merging datasets...")

# Merge all Seurat objects
adata_raw_combined <- merge(adatas[[1]], y = adatas[-1], add.cell.ids = subdata_dirs, project = "combined")

# Calculate mitochondrial percentage
adata_raw_combined <- PercentageFeatureSet(adata_raw_combined, pattern = "^MT-", col.name = "percent.mt")

# %%
# 3. Quality Control (same as pipeline.py)
print("Performing quality control...")

# QC metrics
adata_raw_combined$total_counts <- colSums(adata_raw_combined@assays$RNA@counts)
adata_raw_combined$n_genes_by_counts <- colSums(adata_raw_combined@assays$RNA@counts > 0)

# Filter cells
adata_combined <- subset(adata_raw_combined, 
                         subset = total_counts >= 2500 & 
                           n_genes_by_counts >= 200 & 
                           percent.mt < 15)

# Filter genes
adata_combined <- subset(adata_combined, subset = nFeature_RNA >= 3)

print(paste("After QC:", ncol(adata_combined), "cells,", nrow(adata_combined), "genes"))

# %%
# 4. Normalization and preprocessing (same as pipeline.py)
print("Normalizing and preprocessing...")

# Normalize data
adata_combined <- NormalizeData(adata_combined, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
adata_combined <- FindVariableFeatures(adata_combined, selection.method = "vst", nfeatures = 2000)

# Scale data
adata_combined <- ScaleData(adata_combined, features = rownames(adata_combined))

# %%
# 5. Batch effect correction with Harmony
print("Applying Harmony batch correction...")

# Run PCA
adata_combined <- RunPCA(adata_combined, features = VariableFeatures(object = adata_combined))

# Integrate with Harmony
adata_combined <- IntegrateData(adata_combined, assay = "RNA", normalization.method = "LogNormalize", 
                                dims = 1:40, k.anchor = 5, reduction = "cca")

# Run UMAP on integrated data
adata_combined <- RunUMAP(adata_combined, reduction = "pca", dims = 1:40)

# %%
# 6. Prepare data for Monocle3
print("Preparing data for Monocle3...")

# Extract expression matrix, cell metadata, and gene metadata
expression_matrix <- as.matrix(GetAssayData(adata_combined, assay = "RNA", slot = "data"))
cell_metadata <- adata_combined@meta.data
gene_metadata <- data.frame(
  gene_short_name = rownames(adata_combined),
  biotype = "protein_coding",
  row.names = rownames(adata_combined)
)

# %%
# 7. Create Monocle3 CellDataSet
print("Creating Monocle3 CellDataSet...")

# Create CellDataSet
cds <- new_cell_data_set(
  expression_data = expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

# Preprocess the data
cds <- preprocess_cds(cds, num_dim = 50)

# Reduce dimensions
cds <- reduce_dimension(cds, reduction_method = "UMAP")

# Cluster cells
cds <- cluster_cells(cds, resolution = 1e-3)

# Learn trajectory graph
cds <- learn_graph(cds)

# %%
# 8. Order cells in pseudotime (using EB_D0 as root)
print("Ordering cells in pseudotime...")

# Get EB_D0 cells as root
eb_d0_cells <- rownames(cell_metadata)[cell_metadata$sample == "GSM6304413_EB_D0"]

if (length(eb_d0_cells) > 0) {
  # Order cells using EB_D0 as root
  cds <- order_cells(cds, root_cells = eb_d0_cells[1:5])  # Use first 5 cells as root
  print("Pseudotime calculation completed!")
} else {
  print("Warning: No EB_D0 cells found for root. Pseudotime may not be accurate.")
}

# %%
# 9. Visualize results
print("Creating visualizations...")

# Plot trajectory colored by sample
plot_cells(cds, 
           color_cells_by = "sample",
           show_trajectory_graph = TRUE,
           label_cell_groups = FALSE)
ggsave("figures/monocle3_trajectory_by_sample.png", width = 10, height = 8)

# Plot trajectory colored by pseudotime
plot_cells(cds,
           color_cells_by = "pseudotime",
           show_trajectory_graph = TRUE,
           label_cell_groups = FALSE)
ggsave("figures/monocle3_trajectory_pseudotime.png", width = 10, height = 8)

# Plot trajectory colored by clusters
plot_cells(cds,
           color_cells_by = "cluster",
           show_trajectory_graph = TRUE,
           label_cell_groups = FALSE)
ggsave("figures/monocle3_trajectory_clusters.png", width = 10, height = 8)

# Plot gene expression along trajectory
# Define marker genes
eryth_genes <- c("GYPA", "KLF1", "EPB42")
prog_genes <- c("ANXA1", "CD44", "LMO4")
mega_genes <- c("ITGA3", "ITGA6", "GP1BA", "GP9", "F2R", "CD53")
marker_genes <- c(eryth_genes, prog_genes, mega_genes)

# Plot gene expression trends
plot_genes_in_pseudotime(cds, 
                         marker_genes[1:6],  # Plot first 6 genes
                         color_by = "sample",
                         nrow = 2)
ggsave("figures/monocle3_gene_trends.png", width = 12, height = 8)

# %%
# 10. Save results
print("Saving results...")

# Save CDS object
saveRDS(cds, "data/monocle3_cds.rds")

# Save cell metadata with pseudotime
cell_metadata_with_pseudotime <- data.frame(
  cell_id = rownames(colData(cds)),
  pseudotime = pseudotime(cds),
  cell_type = colData(cds)$cell_type,
  sample = colData(cds)$sample,
  timepoint = colData(cds)$timepoint,
  cluster = clusters(cds)
)

write.csv(cell_metadata_with_pseudotime, "data/monocle3_pseudotime_results.csv", row.names = FALSE)

# Save gene metadata
gene_metadata <- data.frame(
  gene_id = rownames(fData(cds)),
  gene_short_name = fData(cds)$gene_short_name,
  num_cells_expressed = fData(cds)$num_cells_expressed,
  row.names = NULL
)

write.csv(gene_metadata, "data/monocle3_gene_metadata.csv", row.names = FALSE)

print("Monocle3 analysis completed successfully!")
print("Results saved in data/ and figures/ directories")

