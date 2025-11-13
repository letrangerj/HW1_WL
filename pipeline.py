# %%
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

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
# filtering
## 计算线粒体比例
for i, ad in enumerate(adatas):
    # sc.pp.filter_cells(ad, min_genes=200)
    # sc.pp.filter_genes(ad, min_cells=3)
    ad.var["mt"] = ad.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        ad, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

## 绘制QC前指标
fig, axes = plt.subplots(3, 4, figsize=(12, 8))
plt.suptitle("Before QC")
metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
for row_idx, metric in enumerate(metrics):
    for col_idx, ad in enumerate(adatas):
        sc.pl.violin(ad, [metric], jitter=0.4, ax=axes[row_idx, col_idx], show=False)
        # Set titles only for top row and leftmost column for clarity
        if row_idx == 0:
            axes[row_idx, col_idx].set_title(titles[col_idx])
        if col_idx == 0:
            axes[row_idx, col_idx].set_ylabel(metric.replace("_", " ").title())
plt.tight_layout()
plt.savefig("figures/1c_1.png")

## 绘制QC后指标
adata_new = []
for i, ad in enumerate(adatas):
    sc.pp.filter_cells(ad, min_genes=200)
    sc.pp.filter_genes(ad, min_cells=3)
    ad = ad[ad.obs.n_genes_by_counts > 2500, :]
    ad = ad[ad.obs.pct_counts_mt < 15, :]

fig, axes = plt.subplots(3, 4, figsize=(12, 8))
plt.suptitle("After QC")
metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
for row_idx, metric in enumerate(metrics):
    for col_idx, ad in enumerate(adatas):
        sc.pl.violin(ad, [metric], jitter=0.4, ax=axes[row_idx, col_idx], show=False)
        # Set titles only for top row and leftmost column for clarity
        if row_idx == 0:
            axes[row_idx, col_idx].set_title(titles[col_idx])
        if col_idx == 0:
            axes[row_idx, col_idx].set_ylabel(metric.replace("_", " ").title())
plt.tight_layout()
plt.savefig("figures/1c_2.png")

for i, ad in enumerate(adatas):
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    adata_new.append(ad)

# %%
# 2. Batch Effect (Harmony)
import scanpy.external as sce

## 合并数据
adata_combined = adata_new[0].concatenate(
    adata_new[1:],
    batch_key="sample",  # 标记原始来源
    batch_categories=subdata_dirs,  # 名称
    index_unique="-",
)

## High variable genes: 根据scanpy文档，应当在合并数据之后，去除batch之前选取
## 合并后计算能更好地识别真实的细胞异质性
sc.pp.highly_variable_genes(
    adata_combined, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key="sample"
)
adata_combined = adata_combined[:, adata_combined.var.highly_variable]
# 在合并后但在应用 Harmony 之前先 scale + pca
sc.pp.scale(adata_combined, max_value=10)
sc.tl.pca(adata_combined, svd_solver="arpack")
adata_before = adata_combined.copy()

## Harmony 修正
sce.pp.harmony_integrate(adata_combined, key="sample")

## 校正后流程
sc.pp.scale(adata_combined, max_value=10)
sc.tl.pca(adata_combined, svd_solver="arpack")
sc.pp.neighbors(adata_combined, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_combined)

## 校正前流程（无需重复计算PCA）
sc.pp.neighbors(adata_before, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_before)

## UMAP/t-SNE
sc.pl.umap(
    adata_before,
    color="sample",
    title="Umap Before Batch Correlation",
    save="3b_1_UMAP_Before.png",
    show=False,
)
sc.pl.umap(
    adata_combined,
    color="sample",
    title="Umap After Batch Correlation",
    save="3b_2_UMAP_After.png",
    show=False,
)

sc.tl.tsne(adata_before, use_rep="X_pca")
sc.pl.tsne(adata_before, color="sample", save="3b_3_tSNE_Before.png", show=False)

sc.tl.tsne(adata_combined, use_rep="X_pca_harmony")
sc.pl.tsne(adata_combined, color="sample", save="3b_3_tSNE_After.png", show=False)

# %%
# 3. Cell Clustering and Annotation

## Marker Genes from Figure 5e
# 巨核细胞标记基因（来自原论文Figure 5e-g）
megakaryocyte_markers = [
    "PF4",
    "PPBP",
    "ITGA2B",
    "ITGB3",
    "GP9",
    "TUBB1",
    "VWF",
    "NRGN",
]
# 红细胞标记基因
erythroid_markers = ["HBB", "HBA1", "GYPA", "KLF1", "EPOR"]
# 前体细胞标记基因
progenitor_markers = ["CD34", "MPL", "THPO", "KIT"]
# 关键转录因子（原论文强调）
key_transcription_factors = ["FLI1", "MEIS1", "GATA1", "RUNX1", "NF-E2"]

# 可视化dotplot
sc.pl.dotplot(adata_combined,
              megakaryocyte_markers + erythroid_markers + progenitor_markers,
              groupby="leiden",
              save="3_1_rank_genes_groups_dotplot.png")


