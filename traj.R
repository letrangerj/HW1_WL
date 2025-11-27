suppressMessages(suppressWarnings(suppressPackageStartupMessages({
  library(monocle3)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
})))

# %%
# 加载python得到的数据
expression_matrix <- as.matrix(readMM("data/expression_matrix.mtx"))
cell_metadata <- read.csv("data/cell_metadata.csv", row.names = 1)
gene_metadata <- read.csv("data/gene_metadata.csv", row.names = 1)
umap_coords <- as.matrix(read.csv("data/umap_coords.csv", row.names = 1))

expression_matrix <- t(expression_matrix) #需要转置否则细胞和基因对不上
rownames(expression_matrix) <- rownames(gene_metadata) #由于readMM没有行列名，需要补上
colnames(expression_matrix) <- rownames(cell_metadata)
gene_metadata$gene_short_name <- rownames(gene_metadata) #确保gene_short_name列存在

# %%
# 创建 Monocle3 CellDataSet
cds <- new_cell_data_set(
  expression_data = expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

# 预处理数据
cds <- preprocess_cds(cds, num_dim = 50)
reducedDims(cds)[["UMAP"]] <- umap_coords
# 聚类细胞
cds <- cluster_cells(cds, resolution = 1e-3)
# 轨迹图
cds <- learn_graph(cds)

# %%
# 使用EB_D0作为根节点建立伪时间
eb_d0_cells <- rownames(cell_metadata)[cell_metadata$sample == "GSM6304413_EB_D0"]

cds <- order_cells(cds, root_cells = eb_d0_cells[1:5])  # Use first 5 cells as root
print("Pseudotime calculation completed!")


# %%
# 可视化
plot_cells(cds, 
           color_cells_by = "sample",
           show_trajectory_graph = TRUE,
           label_cell_groups = FALSE)
ggsave("figures/monocle3_by_sample.png", width = 10, height = 8)

plot_cells(cds,
           color_cells_by = "pseudotime",
           show_trajectory_graph = TRUE,
           label_cell_groups = FALSE)
ggsave("figures/monocle3_pseudotime.png", width = 10, height = 8)

plot_cells(cds,
           color_cells_by = "cluster",
           show_trajectory_graph = TRUE,
           label_cell_groups = FALSE)
ggsave("figures/monocle3_clusters.png", width = 10, height = 8)

plot_cells(cds,
           color_cells_by = "cell_type",
           show_trajectory_graph = TRUE,
           label_cell_groups = FALSE)
ggsave("figures/monocle3_cell_type.png", width = 10, height = 8)


# marker随轨迹变化
eryth_genes <- c("GYPA", "KLF1", "EPB42")
prog_genes <- c("ANXA1", "CD44", "LMO4")
mega_genes <- c("ITGA3", "ITGA6", "GP1BA", "GP9", "F2R", "CD53")
marker_genes <- c(eryth_genes, prog_genes, mega_genes)

cds_subset <- cds[marker_genes, ]
# 可视化
plot_genes_in_pseudotime(cds_subset,
                         label_by_short_name = TRUE
                         )
ggsave("figures/monocle3_gene_trends.png", width = 12, height = 8)

# Figure 6E中的基因
sub_genes2 <- c("HBA2", "GP9", "CDCA2")
plot_genes_in_pseudotime(cds[sub_genes2, ],
                         label_by_short_name = TRUE)
ggsave("figures/monocle3_subgene2_trends.png", width=8, height = 6)


# %%
# 保留monocle3对象（前面的分析十分耗时，防止需要其他可视化）
save_monocle_objects(cds, "data/monocle3_cds")