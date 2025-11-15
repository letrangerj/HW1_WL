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
    ad_filtered = ad[(ad.obs.n_genes_by_counts > 2500) & (ad.obs.pct_counts_mt < 15)].copy()
    sc.pp.filter_cells(ad_filtered, min_genes=200)
    sc.pp.filter_genes(ad_filtered, min_cells=3)
    adata_new.append(ad_filtered)
    del ad  # 释放内存

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
    
print("QC and normalization completed.")

# %%
# 2. Batch Effect (Harmony)
import scanpy.external as sce

print("Starting batch effect correction with Harmony...")
## 合并数据
adata_combined = adata_new[0].concatenate(
    adata_new[1:],
    batch_key="sample",  # 标记原始来源
    batch_categories=subdata_dirs,  # 名称
    index_unique="-",
)

## High variable genes: 计算高变基因但不进行过滤，保留所有基因用于后续标记基因分析（否则可能找不到marker）
sc.pp.highly_variable_genes(
    adata_combined, min_mean=0.0125, 
    max_mean=3, min_disp=0.5, 
    batch_key="sample"
)

### 转化为稀疏矩阵以节省内存
import scipy.sparse as sp
if not sp.issparse(adata_combined.X):
    print("Converting to sparse matrix format...")
    adata_combined.X = sp.csr_matrix(adata_combined.X)

# 只使用高变基因进行降维，但保留所有基因在数据中
adata_combined.layers["scaled"] = adata_combined.X.copy()
sc.pp.scale(adata_combined, max_value=10, layer="scaled", zero_center=False)
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
del adata_before  # 释放内存

sc.tl.tsne(adata_combined, use_rep="X_pca_harmony")
sc.pl.tsne(adata_combined, color="sample", save="3b_3_tSNE_After.png", show=False)

adata_combined.write('./data/adata_processed.h5ad')
print("Batch effect correction completed and data saved.")

# 释放内存
del adata_combined
import gc
gc.collect()
print("Memory released.")

# %%
# 3. Cell Clustering and Annotation

## Load processed data
adata = sc.read('./data/adata_processed.h5ad')

## Marker Genes from Figure 5e
megakaryocyte_markers = ["ITGA3", "ITGA6", "GP1BA", "GP9", "F2R"]
erythroid_markers = ["GYPA", "KLF1", "EPB42"]
progenitor_markers = ["ANXA1", "CD44", "LMO4"]
cd53_marker = ["CD53"]


## Cluster
sc.tl.leiden(adata, resolution=0.5)
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

sc.pl.dotplot(adata,
              megakaryocyte_markers + erythroid_markers +
              progenitor_markers + cd53_marker,
              groupby='leiden',
              save="3a_markers_dotplot.png")

## 用scanpy内置的scoring功能annotate
def annotation_with_scoring(adata):
    # 计算每组基因的score
    sc.tl.score_genes(adata, megakaryocyte_markers, score_name='mega_score')
    sc.tl.score_genes(adata, erythroid_markers, score_name='eryth_score')
    sc.tl.score_genes(adata, progenitor_markers, score_name='prog_score')
    
    # 基于score注释
    def assign_celltype(row):
        scores = {
            'mega': row['mega_score'],
            'eryth': row['eryth_score'],
            'prog': row['prog_score']
        }
        
        max_score_type = max(scores, key=scores.get)
        
        if max_score_type == 'eryth':
            if scores['prog'] > scores['eryth'] * 0.5:
                return "C2: EB-like"
            return "C1: EB"
        elif max_score_type == 'prog':
            return "C3: iPEM"
        elif max_score_type == 'mega':
            if row['CD53'] > adata.obs['CD53'].median():
                return "C5: iMK-2"
            return "C4: iMK-1"
        return "Unknown"
    
    # 获取CD53表达
    cd53_expr = adata[:, 'CD53'].X
    if hasattr(cd53_expr, 'toarray'):
        cd53_expr = cd53_expr.toarray().flatten()
    adata.obs['CD53'] = cd53_expr
    
    adata.obs['cell_type_v2'] = adata.obs.apply(assign_celltype, axis=1)
    
    return adata

adata = annotation_with_scoring(adata)

sc.pl.umap(adata, 
           color=['leiden', 'cell_type'], 
           save="3b_annotated_clusters.png")


# %%
