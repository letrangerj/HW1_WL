# %%
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore", message="dendrogram data not found")
warnings.filterwarnings("ignore", message="Groups are not reordered")
# %%
# Load datasets
result_file = "analysis_results.h5ad"
data_path = "data"
subdata_dirs = [
    "GSM6304413_EB_D0",
    "GSM6304414_iMK_D3",
    "GSM6304415_iMK_D5",
    "GSM6304416_iMK_D7",
]
titles = ["EB_D0", "iMK_D3", "iMK_D5", "iMK_D7"]

adatas = []
for subdir in subdata_dirs:
    adata = sc.read_10x_mtx(
        f"{data_path}/{subdir}/", var_names="gene_symbols", cache=True
    )
    adata.var_names_make_unique()
    adatas.append(adata)

print(f"Loaded {len(adatas)} datasets.")

# %%
# Merge datasets BEFORE QC
adata_raw_combined = adatas[0].concatenate(
    adatas[1:],
    batch_key="sample",  
    batch_categories=subdata_dirs,  
    index_unique="-",
)

# 计算线粒体比例并做 QC metrics
adata_raw_combined.var["mt"] = adata_raw_combined.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata_raw_combined, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

# 绘制 QC 前整体指标（不分组）
metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
fig, axes = plt.subplots(1, 3, figsize=(14, 4))
plt.suptitle("Before QC (Merged Overall)")
for i, metric in enumerate(metrics):
    sc.pl.violin(
        adata_raw_combined,
        [metric],
        jitter=0.4,
        ax=axes[i],
        show=False,
    )
    axes[i].set_ylabel(metric.replace("_", " ").title())
    axes[i].set_title(metric.replace("_", " ").title())
plt.tight_layout()
plt.savefig("figures/1c_1.png")

# 过滤（merged data）
adata_combined = adata_raw_combined[
    (adata_raw_combined.obs.total_counts >= 2500)
    & (adata_raw_combined.obs.pct_counts_mt < 15)
].copy()
sc.pp.filter_cells(adata_combined, min_genes=200)
sc.pp.filter_genes(adata_combined, min_cells=3)

# 绘制 QC 后整体指标
fig, axes = plt.subplots(1, 3, figsize=(14, 4))
plt.suptitle("After QC (Merged Overall)")
for i, metric in enumerate(metrics):
    sc.pl.violin(
        adata_combined,
        [metric],
        jitter=0.4,
        ax=axes[i],
        show=False,
    )
    axes[i].set_ylabel(metric.replace("_", " ").title())
    axes[i].set_title(metric.replace("_", " ").title())
plt.tight_layout()
plt.savefig("figures/1c_2.png")

# 对过滤后的 merged 数据进行标准化
sc.pp.normalize_total(adata_combined, target_sum=1e4)
sc.pp.log1p(adata_combined)

print("QC and normalization completed on merged data.")

# %%
# 2. Batch Effect (Harmony)
import scanpy.external as sce

print("Starting batch effect correction with Harmony...")

## High variable genes: 计算高变基因但不进行过滤，保留所有基因用于后续标记基因分析（否则可能找不到marker）
sc.pp.highly_variable_genes(
    adata_combined,
    n_top_genes=2000,  # Match paper
    batch_key="sample",
    flavor="seurat",
)

### 转化为稀疏矩阵以节省内存
import scipy.sparse as sp

adata_combined.X = sp.csr_matrix(adata_combined.X)

# 只使用高变基因进行降维，但保留所有基因在数据中
adata_combined.layers["scaled"] = adata_combined.X.copy()
sc.pp.scale(adata_combined, max_value=10, layer="scaled", zero_center=True)
sc.tl.pca(adata_combined, svd_solver="arpack", use_highly_variable=True)
adata_before = adata_combined.copy()

## Harmony 修正
sce.pp.harmony_integrate(adata_combined, key="sample")

## 校正后流程 - 使用高变基因进行PCA
sc.pp.neighbors(adata_combined, n_neighbors=10, n_pcs=40, use_rep="X_pca_harmony")
sc.tl.umap(adata_combined)

## 校正前流程（无需重复计算PCA）
sc.pp.neighbors(adata_before, n_neighbors=10, n_pcs=40, use_rep="X_pca")
sc.tl.umap(adata_before)

## UMAP/t-SNE
sc.pl.umap(
    adata_before,
    color="sample",
    title="Umap Before Batch Correlation",
    save="3b_1_UMAP_Before.png",
    show=True,
)
sc.pl.umap(
    adata_combined,
    color="sample",
    title="Umap After Batch Correlation",
    save="3b_2_UMAP_After.png",
    show=True,
)

sc.tl.tsne(adata_before, use_rep="X_pca")
sc.pl.tsne(adata_before, color="sample", save="3b_3_tSNE_Before.png", show=True)
del adata_before  # 释放内存

sc.tl.tsne(adata_combined, use_rep="X_pca_harmony")
sc.pl.tsne(adata_combined, color="sample", save="3b_3_tSNE_After.png", show=True)
adata_combined.write("./data/adata_processed.h5ad")
print("Batch effect correction completed and data saved.")

# 释放内存
del adata_combined
import gc

gc.collect()
print("Memory released.")

# %%
# 3. Cell Clustering and Annotation
# Using the processed data after batch correction and leiden clustering
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

## Load processed data
adata = sc.read("./data/adata_processed.h5ad")

eryth_genes = ["GYPA", "KLF1", "EPB42"]
prog_genes = ["ANXA1", "CD44", "LMO4"]
mega_genes = ["ITGA3", "ITGA6", "GP1BA", "GP9", "F2R", "CD53"]
marker_genes = {
    "Erythroid (EB)": eryth_genes,
    "Precursor (iPEM)": prog_genes,
    "Megakaryocyte (iMK)": mega_genes,
}

sc.tl.leiden(adata, resolution=0.35)
sc.tl.dendrogram(adata, groupby="leiden")
sc.pl.dotplot(adata, marker_genes, groupby="leiden", save="3c_1_dotplot.png")

# %%
# Cell Type Annotation based on Marker Genes
cluster_annotation = {
    "1": "C1_EB",
    "2": "C1_EB",
    "6": "C1_EB",
    "3": "C2_EB-like",
    "5": "C2_EB-like",
    "0": "C3_iPEM",
    "7": "C3_iPEM",
    "9": "C3_iPEM",
    "10": "C3_iPEM",
    "11": "C3_iPEM",  # CD53+ but MK markers are too weak to be C5
    "4": "C4_iMK-1",  # High MK markers, Low CD53
    "8": "C5_iMK-2",  # High MK markers, High CD53
}
adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_annotation).astype("category")
sc.pl.dotplot(
    adata, marker_genes, groupby="cell_type", standard_scale="var", dendrogram=True
)
sc.pl.tsne(
    adata,
    color="cell_type",
    legend_loc="on data",
    title="Cell Type Annotation in t-SNE",
    save="3b_1_tSNE_Cell_Type.png"
)
sc.pl.umap(
    adata,
    color="cell_type",
    legend_loc="on data",
    title="Cell Type Annotation in UMAP",
    save="3b_2_UMAP_Cell_Type.png"
)
# %%
# 3c. Plot heatmap and violin plot for marker genes within annotated cell types

# 改名
cluster_name_map = {
    "C1_EB": "C1",
    "C2_EB-like": "C2",
    "C3_iPEM": "C3",
    "C4_iMK-1": "C4",
    "C5_iMK-2": "C5",
}
sample_name_map = {
    "GSM6304413_EB_D0": "D0",
    "GSM6304414_iMK_D3": "D3",
    "GSM6304415_iMK_D5": "D5",
    "GSM6304416_iMK_D7": "D7",
}
adata.obs["cell_type"] = adata.obs["cell_type"].cat.rename_categories(cluster_name_map)

sc.pl.dotplot(
    adata,
    marker_genes,
    groupby="cell_type",
    standard_scale="var",
    cmap="viridis",
    swap_axes=True,
    dot_max=1.0,
    title="Heatmap of Marker Genes by Cell Type",
    show=True,
    save="3c_1_Dotplot_Marker_Genes.png",
)
sc.pl.stacked_violin(
    adata,
    marker_genes,
    groupby="cell_type",
    swap_axes=True,
    dendrogram=False,
    standard_scale=None,
    cmap="viridis",
    row_palette="tab20",
    show=True,
    title="Violin Plot of Marker Genes by Cell Type",
    save="3c_2_StackedViolin_Marker_Genes.png",
)

# Visualize sample proportions for each cell type cluster
# 计算频率矩阵
freq_table = pd.crosstab(adata.obs["cell_type"], adata.obs["sample"])
freq_table_norm = freq_table.div(freq_table.sum(axis=1), axis=0)
freq_table_norm.rename(index=cluster_name_map, columns=sample_name_map, inplace=True)
freq_table_norm.plot(
    kind="bar",
    stacked=True,
    colormap="Set3",
    figsize=(5, 4),
    width=0.8,
    edgecolor="none",
    title="Sample Proportions in Each Cell Type Cluster",
    ylabel="Proportion",
    xlabel="Cell Type Cluster",
)
plt.savefig("figures/3c_3_Sample_Proportions_by_Cell_Type.png")
adata.write("./data/adata_annotated.h5ad")
# %%
# 4. GO Enrichment Analysis
## 4a. 2 cluster's GO enrichment analysis
from gseapy import enrichr

adata = sc.read("./data/adata_annotated.h5ad")

sc.tl.rank_genes_groups(
    adata,
    groupby="cell_type",
    method="wilcoxon",  # or 'logreg' or 't-test'
    pts=True,  # Optional: computes fraction of cells expressing gene
)

clusters = ["C1", "C3"]
enrichr_results = {}
for cluster in clusters:
    ## 得到富集的基因的dataframe
    df_top = sc.get.rank_genes_groups_df(adata, group=cluster)

    filtered_df = df_top[
        (df_top["pvals_adj"] < 0.05) & (df_top["logfoldchanges"] > 0.5)
    ]
    gene_lst = filtered_df["names"].tolist()

    enr = enrichr(
        gene_list=gene_lst,
        gene_sets="GO_Biological_Process_2023",
        organism="Human",
        cutoff=0.05,
        outdir=f"enrichr_results/{cluster}",
    )
    print(f"{cluster}'s GO results: ", enr.results.head())
    enrichr_results[cluster] = enr.results

# %%
# 4b. GO可视化
import seaborn as sns

for cluster in clusters:
    go_df = enrichr_results[cluster].copy()
    top_terms = go_df.sort_values("Adjusted P-value").head(10).copy()

    # Calculate -log10 P-value
    top_terms["log_pval"] = -np.log10(top_terms["Adjusted P-value"])
    # Clean the Term names (Remove " (GO:XXXXXXX)")
    top_terms["Term_Clean"] = top_terms["Term"].str.split(" \(GO").str[0]

    plt.figure(figsize=(8, 6))
    barplot = sns.barplot(
        data=top_terms, x="log_pval", y="Term_Clean", dodge=False, legend=False
    )
    for i, (value, p_val) in enumerate(
        zip(top_terms["log_pval"], top_terms["Adjusted P-value"])
    ):
        label = f"p={p_val:.1e}"
        plt.text(value + 0.2, i, label, va="center", color="black", fontsize=9)
    plt.title(
        f"Top 10 Enriched GO Terms: {cluster}", fontsize=14, fontweight="bold", pad=15
    )
    plt.xlabel("Significance (-log10 Adj. P-value)", fontsize=12)
    plt.ylabel("")
    # plt.tight_layout()

    plt.savefig(f"figures/GO_Analysis_{cluster}.png", dpi=300, bbox_inches='tight')
    plt.show()

# %%
# Advanced Topic，使用monocle3进行轨迹推断分析，这里准备所需要的文件
import scanpy as sc
import pandas as pd
from scipy import io
import scipy.sparse as sp

adata = sc.read("data/adata_annotated.h5ad")

X = adata.X
io.mmwrite("data/expression_matrix.mtx", X)

adata.obs.to_csv("data/cell_metadata.csv", index=True)   # 行名为 cell IDs
adata.var.to_csv("data/gene_metadata.csv", index=True)   # 行名为 gene names

pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names).to_csv("data/umap_coords.csv")
