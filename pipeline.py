# %%
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import warnings

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
    batch_key="sample",  # 标记原始来源
    batch_categories=subdata_dirs,  # 名称
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
## 合并数据
# 数据已在 QC 前合并并过滤，直接复用 adata_combined

## High variable genes: 计算高变基因但不进行过滤，保留所有基因用于后续标记基因分析（否则可能找不到marker）
sc.pp.highly_variable_genes(
    adata_combined, 
    n_top_genes=2000,  # Match paper
    batch_key="sample",
    flavor='seurat'
)

### 转化为稀疏矩阵以节省内存
import scipy.sparse as sp
if not sp.issparse(adata_combined.X):
    print("Converting to sparse matrix format...")
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
adata_combined.write('./data/adata_processed.h5ad')
print("Batch effect correction completed and data saved.")

# 释放内存
del adata_combined
import gc
gc.collect()
print("Memory released.")

# %%
# 3. Cell Clustering and Annotation

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

## Load processed data
adata = sc.read('./data/adata_processed.h5ad')

eryth_genes = ["GYPA", "KLF1", "EPB42"]
prog_genes = ["ANXA1", "CD44", "LMO4"]
mega_genes = ["ITGA3", "ITGA6", "GP1BA", "GP9", "F2R"]
cd53_genes = ["CD53"]

# Leiden clustering, resolution is the same as paper
sc.tl.leiden(adata, resolution=1.5e-5, key_added='leiden')

print(f"Number of clusters: {adata.obs['leiden'].nunique()}")
print(f"Cluster sizes:\n{adata.obs['leiden'].value_counts().sort_index()}")

# CALCULATE GENE SCORES
sc.tl.score_genes(adata, eryth_genes, score_name='erythroid_score')
sc.tl.score_genes(adata, prog_genes, score_name='progenitor_score')
sc.tl.score_genes(adata, mega_genes, score_name='megakaryocyte_score')

# Get CD53 expression
if 'CD53' in adata.var_names:
    cd53_expr = adata[:, 'CD53'].X
    if hasattr(cd53_expr, 'toarray'):
        cd53_expr = cd53_expr.toarray().flatten()
    else:
        cd53_expr = cd53_expr.flatten()
    adata.obs['CD53_expression'] = cd53_expr
else:
    adata.obs['CD53_expression'] = 0
    print("Warning: CD53 not found")

# ANNOTATE CELL TYPES
# Normalize scores to 0-1 range for comparison
for score in ['erythroid_score', 'progenitor_score', 'megakaryocyte_score']:
    score_min = adata.obs[score].min()
    score_max = adata.obs[score].max()
    adata.obs[f'{score}_norm'] = (adata.obs[score] - score_min) / (score_max - score_min)

def assign_celltype_minimal(row):
    """
    Minimal annotation based only on gene expression scores.
    Logic:
    - High erythroid, low mega → C1 (EB)
    - Medium erythroid + progenitor → C2 (EB-like)
    - High progenitor OR balanced eryth+mega → C3 (iPEM)
    - High mega + low CD53 → C4 (iMK-1)
    - High mega + high CD53 → C5 (iMK-2)
    """
    
    eryth = row['erythroid_score_norm']
    prog = row['progenitor_score_norm']
    mega = row['megakaryocyte_score_norm']
    cd53 = row['CD53_expression']
    
    # Thresholds
    HIGH = 0.6
    MED = 0.35
    LOW = 0.25
    
    # C1: EB - Dominant erythroid
    if eryth > HIGH and mega < LOW:
        return "C1_EB"
    
    # C2: EB-like - Erythroid with emerging progenitor features
    if eryth > MED and prog > MED and mega < MED:
        return "C2_EB-like"
    
    # C3: iPEM - High progenitor OR balanced intermediate state
    if prog > HIGH:
        return "C3_iPEM"
    if eryth > MED and mega > MED and abs(eryth - mega) < 0.3:
        return "C3_iPEM"
    
    # C4/C5: iMK - Dominant megakaryocyte
    if mega > HIGH:
        cd53_median = adata.obs['CD53_expression'].median()
        if cd53 > cd53_median:
            return "C5_iMK-2"
        else:
            return "C4_iMK-1"

    return f"unknown"

adata.obs['cell_type'] = adata.obs.apply(assign_celltype_minimal, axis=1)

# 1. UMAP with annotations
fig, axes = plt.subplots(2, 2, figsize=(15, 10))

sc.pl.umap(adata, color='cell_type', ax=axes[0, 0], show=False, title='Cell Type')
sc.pl.umap(adata, color='erythroid_score', ax=axes[0, 1], show=False, 
           title='Erythroid Score', cmap='Reds')
sc.pl.umap(adata, color='progenitor_score', ax=axes[1, 0], show=False,
           title='Progenitor Score', cmap='Greens')
sc.pl.umap(adata, color='megakaryocyte_score', ax=axes[1, 1], show=False,
           title='Megakaryocyte Score', cmap='Blues')

plt.tight_layout()
plt.show()

# 2. Dotplot of marker genes
all_markers = eryth_genes + prog_genes + mega_genes + cd53_genes
sc.pl.dotplot(adata, all_markers, groupby='cell_type',
              figsize=(8, 4))

# 3. Stacked violin plots
sc.pl.stacked_violin(adata, all_markers, groupby='cell_type',
                     figsize=(10, 4))
# 4. Score distributions by cell type
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

cell_types = adata.obs['cell_type'].cat.categories
scores = ['erythroid_score', 'progenitor_score', 'megakaryocyte_score']
titles = ['Erythroid Score', 'Progenitor Score', 'Megakaryocyte Score']

for idx, (score, title) in enumerate(zip(scores, titles)):
    data = [adata.obs[adata.obs['cell_type'] == ct][score].values 
            for ct in cell_types]
    axes[idx].violinplot(data, positions=range(len(cell_types)), 
                         showmeans=True, showmedians=True)
    axes[idx].set_xticks(range(len(cell_types)))
    axes[idx].set_xticklabels(cell_types, rotation=45, ha='right')
    axes[idx].set_ylabel('Score')
    axes[idx].set_title(title)
    axes[idx].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

print("Analysis completed.")