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
    # 创建副本进行过滤
    ad_filtered = ad.copy()
    sc.pp.filter_cells(ad_filtered, min_genes=200)
    sc.pp.filter_genes(ad_filtered, min_cells=3)
    ad_filtered = ad_filtered[ad_filtered.obs.n_genes_by_counts > 2500, :]
    ad_filtered = ad_filtered[ad_filtered.obs.pct_counts_mt < 15, :]
    adata_new.append(ad_filtered)

fig, axes = plt.subplots(3, 4, figsize=(12, 8))
plt.suptitle("After QC")
metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
for row_idx, metric in enumerate(metrics):
    for col_idx, ad in enumerate(adata_new):  # 使用过滤后的数据
        sc.pl.violin(ad, [metric], jitter=0.4, ax=axes[row_idx, col_idx], show=False)
        # Set titles only for top row and leftmost column for clarity
        if row_idx == 0:
            axes[row_idx, col_idx].set_title(titles[col_idx])
        if col_idx == 0:
            axes[row_idx, col_idx].set_ylabel(metric.replace("_", " ").title())
plt.tight_layout()
plt.savefig("figures/1c_2.png")

# 对过滤后的数据进行标准化
for i, ad in enumerate(adata_new):
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)

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

## High variable genes: 计算高变基因但不进行过滤，保留所有基因用于后续标记基因分析（否则可能找不到marker）
sc.pp.highly_variable_genes(
    adata_combined, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key="sample"
)
# 只使用高变基因进行降维，但保留所有基因在数据中
adata_combined.layers["scaled"] = adata_combined.X.copy()
sc.pp.scale(adata_combined, max_value=10, layer="scaled")
sc.tl.pca(adata_combined, svd_solver="arpack", use_highly_variable=True)
adata_before = adata_combined.copy()

## Harmony 修正
sce.pp.harmony_integrate(adata_combined, key="sample")

## 校正后流程 - 使用高变基因进行PCA
sc.tl.pca(adata_combined, svd_solver="arpack", use_highly_variable=True)
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
megakaryocyte_markers = ["ITGA3", "ITGA6", "GP1BA", "GP9", "F2R", "SELP"]
erythroid_markers = ["GYPA", "KLF1", "EPB42"]
progenitor_markers = ["ANXA1", "CD44", "LMO4"]
cd53_marker = ["CD53"]


## Cluster
sc.tl.leiden(adata_combined, resolution=0.5)
sc.tl.rank_genes_groups(adata_combined, 'leiden', method='wilcoxon')

sc.pl.dotplot(adata_combined,
              megakaryocyte_markers + erythroid_markers +
              progenitor_markers + cd53_marker,
              groupby='leiden',
              save="3a_markers_dotplot.png")

def annotation_clusters(adata):
    cluster_annotation = {}
    for cluster_id in adata.obs['leiden'].cat.categories:
        cluster_data = adata[adata.obs['leiden'] == cluster_id]
        mega_expr = cluster_data[:, megakaryocyte_markers].X.mean()
        eryth_expr = cluster_data[:, erythroid_markers].X.mean()
        prog_expr = cluster_data[:, progenitor_markers].X.mean()
        cd53_expr = cluster_data[:, cd53_marker].X.mean()
        
        # 根据表达模式注释
        if eryth_expr > 0.5 and mega_expr < 0.2:
            cluster_annotation[cluster_id] = "C1: EB"
        elif eryth_expr > 0.2 and eryth_expr < 0.5 and prog_expr > 0.3:
            cluster_annotation[cluster_id] = "C2: EB-like"
        elif prog_expr > 0.5 and eryth_expr < 0.2 and mega_expr > 0.2:
            cluster_annotation[cluster_id] = "C3: iPEM"
        elif mega_expr > 0.5 and eryth_expr < 0.2:
            if cd53_expr > 0.3:
                cluster_annotation[cluster_id] = "C5: iMK-2"
            else:
                cluster_annotation[cluster_id] = "C4: iMK-1"
        else:
            cluster_annotation[cluster_id] = f"Unknown_{cluster_id}"
    
    return cluster_annotation

cluster_names_dict = annotation_clusters(adata_combined)
adata_combined.obs['cell_type'] = adata_combined.obs['leiden'].map(cluster_names_dict)

sc.pl.umap(adata_combined, color=['leiden', 'cell_type'], 
           save="3b_annotated_clusters.png", show=False)

