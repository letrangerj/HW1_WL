# 单细胞RNA测序数据分析作业指导文档

## 背景简介
本作业分析红细胞(EB)向诱导性巨核细胞(iMK)重编程过程中的单细胞RNA测序数据，包含四个时间点：EB (d0)、iMKd3、iMKd5、iMKd7。

数据来源：GEO数据库 GSE207654

---

## 任务0：环境搭建 (10分)

### 关键步骤
使用Conda/Miniconda创建并配置Python环境

### 必需的包
```bash
# 创建环境
conda create -n scrna python=3.9
conda activate scrna

# 安装scanpy及相关包
pip install scanpy
pip install leidenalg  # 用于聚类
pip install python-igraph  # Leiden算法依赖
pip install harmonypy  # 批次效应校正（可选）
pip install scanorama  # 批次效应校正（可选）
pip install gseapy  # GO富集分析
```

### 参考文档
- `scRNA_seq_analyze_totorial_2025.html` - 环境配置部分
- Scanpy官方教程：https://scanpy-tutorials.readthedocs.io/

---

## 任务1：数据预处理与质量控制 (10分)

### 1a. 基本预处理流程

#### 关键函数与来源

**数据读取** (`scRNA_seq_analyze_totorial_2025.html`)
```python
import scanpy as sc
import numpy as np

# 读取10x数据
adata = sc.read_10x_mtx('path/to/filtered_feature_bc_matrix/')
```

**过滤低表达基因** 
```python
# 过滤在少于3个细胞中表达的基因
sc.pp.filter_genes(adata, min_cells=3)
```
- 函数来源：`scanpy.pp.filter_genes`
- 文档参考：`scRNA_seq_analyze_totorial_2025.html` QC部分

**计算QC指标**
```python
# 计算线粒体基因百分比
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], 
                           percent_top=None, 
                           log1p=False, 
                           inplace=True)
```
- 函数来源：`scanpy.pp.calculate_qc_metrics`
- 文档参考：`scRNA_seq_analyze_totorial_2025.html`

**细胞过滤**
```python
# 过滤低质量细胞
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.n_genes_by_counts > 200, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
```
- 阈值需根据实际数据调整
- 文档参考：`scRNA_seq_analyze_totorial_2025.html` filtering部分

**归一化**
```python
# 标准化每个细胞的总counts
sc.pp.normalize_total(adata, target_sum=1e4)
# log转换
sc.pp.log1p(adata)
```
- 函数来源：`scanpy.pp.normalize_total`, `scanpy.pp.log1p`
- 文档参考：`scRNA_seq_analyze_totorial_2025.html` normalization部分

### 1b. 线粒体基因百分比作为QC标准的原因

**关键概念**：
- 细胞死亡或破损时，细胞质mRNA会泄漏流失
- 线粒体基因的mRNA相对稳定，保留在细胞中
- 因此死亡/低质量细胞会表现出较高的线粒体基因表达比例
- 通常设定阈值为5-20%

**参考**：`scRNA_seq_analyze_totorial_2025.html` 质量控制章节

### 1c. QC指标可视化

**关键图表**：
```python
# QC前的可视化
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
             
# 或使用散点图
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
```

**指标含义**：
- `n_genes_by_counts` / `nFeature_RNA`: 每个细胞检测到的基因数量
- `total_counts` / `nCount_RNA`: 每个细胞的总UMI/reads数
- `pct_counts_mt` / `percent.mt`: 线粒体基因表达百分比

**文档参考**：`scRNA_seq_analyze_totorial_2025.html` visualization部分

---

## 任务2：批次效应去除 (30分)

### 2a. 批次效应校正算法选择

#### 推荐完整工作流程

对于多批次单细胞数据，推荐的完整分析流程如下：

**1. 数据预处理与合并**
```python
# 标准化所有数据集
for ad in adatas:
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)

# 合并数据
adata_combined = adata_new[0].concatenate(
    adata_new[1:],
    batch_key="sample",
    batch_categories=subdata_dirs,
    index_unique="-",
)

# 在合并数据上计算高可变基因
sc.pp.highly_variable_genes(
    adata_combined,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5,
    batch_key='sample'
)
adata_combined = adata_combined[:, adata_combined.var.highly_variable]

# 创建备份用于比较
adata_before = adata_combined.copy()
```

**2. 批次效应校正（Harmony）**
```python
import scanpy.external as sce

# 标准化和PCA
sc.pp.scale(adata_combined, max_value=10)
sc.tl.pca(adata_combined, svd_solver='arpack')

# 使用Harmony校正批次效应
sce.pp.harmony_integrate(adata_combined, key='sample')  # 'sample'是批次标签列名

# 计算邻居和降维（使用Harmony处理后的嵌入）
sc.pp.neighbors(adata_combined, n_neighbors=10, n_pcs=40, use_rep='X_pca_harmony')
sc.tl.umap(adata_combined)
sc.tl.tsne(adata_combined, use_rep='X_pca_harmony')

# 同样处理批处理前的数据
sc.pp.scale(adata_before, max_value=10)
sc.tl.pca(adata_before, svd_solver='arpack')
sc.pp.neighbors(adata_before, n_neighbors=10, n_pcs=40, use_rep='X_pca')
sc.tl.umap(adata_before)
sc.tl.tsne(adata_before, use_rep='X_pca')
```

**关键点说明：**
- Harmony在PCA空间中运行，需要先计算PCA
- 批处理后的降维使用`X_pca_harmony`嵌入
- 批处理前的降维使用`X_pca`嵌入
- 两者可视化比较可以评估批次效应移除效果

- 优点：速度快，效果好，保持生物学差异
- 适用：多个样本/批次整合
- 文档参考：Harmony原始论文及Scanpy文档

**Scanorama**
```python
import scanorama

# 按批次分割数据
adatas = [adata[adata.obs.sample == s] for s in adata.obs.sample.unique()]

# 整合
scanorama.integrate_scanpy(adatas, dimred=50)
```
- 优点：可处理大规模数据，跨数据集整合能力强
- 文档参考：`scRNA_seq_analyze_totorial_2025.html` batch correction部分

**BBKNN** (Graph-based)
```python
import bbknn

sc.tl.pca(adata)
bbknn.bbknn(adata, batch_key='sample')
```
- 优点：基于图的方法，可直接用于聚类
- 参考：BBKNN原始论文

**算法选择理由**：
- 样本数量：4个时间点 → Harmony或Scanorama均适用
- 生物学意义：需要保留重编程的时间梯度 → Harmony更合适
- 计算效率：样本量中等 → Harmony速度更快

### 2b. 批次效应前后的可视化

**降维与可视化**
```python
# 批次效应校正前
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.pl.umap(adata, color='sample', title='Before batch correction')

# 批次效应校正后（使用Harmony例子）
sce.pp.harmony_integrate(adata, key='sample')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40, use_rep='X_pca_harmony')
sc.tl.umap(adata)
sc.pl.umap(adata, color='sample', title='After batch correction')

# 也可以用t-SNE
sc.tl.tsne(adata, use_rep='X_pca_harmony')
sc.pl.tsne(adata, color='sample')
```

**参考原文图表**：Figure 5a, 5b - 展示不同样本的UMAP分布

**评估标准**：
- 校正前：不同批次（样本）分离明显
- 校正后：相同细胞类型混合，但保留生物学差异

**文档参考**：
- `scRNA_seq_analyze_totorial_2025.html` dimensionality reduction章节
- Scanpy tutorials - PBMC3k教程

---

## 任务3：细胞聚类与注释 (30分)

### 3a. 细胞聚类

**关键步骤**：

**1. 寻找高变基因（推荐方法）**

对于多批次数据，推荐流程是：
1. 先合并所有数据集
2. 在合并数据上统一计算高可变基因
3. 然后进行批次效应校正

```python
# 合并数据
adata_combined = adata_new[0].concatenate(
    adata_new[1:],
    batch_key="sample",  # 批次标识
    batch_categories=subdata_dirs,
    index_unique="-",
)

# 在合并数据上计算高可变基因，考虑批次效应
sc.pp.highly_variable_genes(
    adata_combined,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5,
    batch_key='sample'  # 关键参数：计算时考虑批次结构
)

# 仅保留高可变基因
adata_combined = adata_combined[:, adata_combined.var.highly_variable]
```

**为什么这种方法更好？**
- 避免了各批次选择不同的高可变基因
- 能够识别跨批次的真实生物学变异
- 统一的基因集确保后续分析的一致性
- `batch_key='sample'`参数确保批次结构不影响基因选择

- 函数来源：`scanpy.pp.highly_variable_genes`
- 参考：scanpy官方文档和单细胞分析最佳实践

**2. 数据缩放**
```python
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
```

**3. 降维**
```python
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
```

**4. 聚类（Leiden算法）**
```python
sc.tl.leiden(adata, resolution=0.5)
```
- 函数来源：`scanpy.tl.leiden`
- resolution参数控制聚类数量，需要调整
- 参考：`scRNA_seq_analyze_totorial_2025.html` clustering章节

**寻找marker基因**
```python
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
```

**根据原论文Figure 5e-g的标记基因进行注释**

基于原论文，巨核细胞分化过程中的关键标记基因包括：

```python
# 巨核细胞标记基因（来自原论文Figure 5e-g）
megakaryocyte_markers = ['PF4', 'PPBP', 'ITGA2B', 'ITGB3', 'GP9', 'TUBB1', 'VWF', 'NRGN']
# 红细胞标记基因
erythroid_markers = ['HBB', 'HBA1', 'GYPA', 'KLF1', 'EPOR']
# 前体细胞标记基因
progenitor_markers = ['CD34', 'MPL', 'THPO', 'KIT']
# 关键转录因子（原论文强调）
key_transcription_factors = ['FLI1', 'MEIS1', 'GATA1', 'RUNX1', 'NF-E2']

# 可视化标记基因表达
sc.pl.dotplot(adata, 
              megakaryocyte_markers + erythroid_markers + progenitor_markers,
              groupby='leiden', 
              save='_megakaryocyte_markers.png')

# 查看关键转录因子
sc.pl.violin(adata, key_transcription_factors, groupby='leiden', 
             save='_key_tfs.png')
```

### 3b. UMAP/t-SNE可视化

```python
# 按聚类结果着色
sc.pl.umap(adata, color='leiden', legend_loc='on data')

# 按样本来源着色（验证批处理效果）
sc.pl.umap(adata, color='sample', title='UMAP by sample')

# 按关键转录因子着色
sc.pl.umap(adata, color=['FLI1', 'MEIS1'], save='_tfs_on_umap.png')
```

# 按样本着色
sc.pl.umap(adata, color='sample')

# 多面板展示
sc.pl.umap(adata, color=['leiden', 'sample'], ncols=2)
```

**细胞类型注释**：
根据原文Figure 5e-g的marker基因进行注释：
- HSC/MPP：CD34, KIT, CD38
- MEP：ITGA2B, GP1BA, FLI1
- Erythroid：HBA1, HBB, GYPA
- Megakaryocyte：PF4, PPBP, ITGA2B

```python
# 手动注释示例
cluster2annotation = {
    '0': 'HSC/MPP',
    '1': 'MEP',
    '2': 'Erythroid',
    '3': 'Megakaryocyte',
    # ... 根据marker基因判断
}
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster2annotation)
sc.pl.umap(adata, color='cell_type')
```

**参考**：
- 原文Figure 5b
- `scRNA_seq_analyze_totorial_2025.html` cell type annotation

### 3c. Marker基因热图与小提琴图

**热图（Heatmap）**
```python
# 选择每个cluster的top marker基因
sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, groupby='cell_type', 
                                cmap='viridis', show_gene_labels=True)

# 或者指定特定marker基因
marker_genes = ['CD34', 'KIT', 'ITGA2B', 'GP1BA', 'HBA1', 'PF4']
sc.pl.heatmap(adata, marker_genes, groupby='cell_type', 
              cmap='RdYlBu_r', standard_scale='var')
```
- 参考：原文Figure 5d, 5h

**小提琴图（Violin plot）**
```python
sc.pl.violin(adata, marker_genes, groupby='cell_type', rotation=90)

# 或者点状小提琴图
sc.pl.stacked_violin(adata, marker_genes, groupby='cell_type')
```

**文档参考**：
- `scRNA_seq_analyze_totorial_2025.html` visualization章节
- Scanpy plotting tutorial

### 3d. 细胞类型在样本中的比例

```python
import pandas as pd
import matplotlib.pyplot as plt

# 计算每个样本中各细胞类型的比例
cell_counts = pd.crosstab(adata.obs['cell_type'], adata.obs['sample'], 
                         normalize='columns') * 100

# 堆叠柱状图（参考原论文Figure 5c）
cell_counts.T.plot(kind='bar', stacked=True, figsize=(10, 6), 
                  color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'])  # 使用不同颜色区分细胞类型
plt.ylabel('Percentage (%)')
plt.xlabel('Sample')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title('Cell type composition across samples\n(Accompanying EB→iMK reprogramming process)', fontsize=12)
plt.tight_layout()
plt.savefig('figures/3d_cell_composition.pdf', dpi=300)

# 创建额外的分析表格（与原论文类似）
cell_counts_abs = pd.crosstab(adata.obs['cell_type'], adata.obs['sample'])
cell_counts_abs.to_csv('figures/3d_cell_counts_absolute.csv')
```

**生物学意义解读**（参考原论文Figure 5c）：
- EB(d0)样本中主要包含红细胞谱系细胞
- 随着重编程进程(iMKd3→iMKd5→iMKd7)，巨核细胞比例逐渐增加
- 这种比例变化反映了红系到巨核系的转分化过程
- 可以量化重编程效率：后期样本中巨核细胞比例

**参考**：
- 原文Figure 5c
- `scRNA_seq_analyze_totorial_2025.html` composition analysis

---

## 任务4：GO富集分析 (10分)

### 4a. GO富集分析方法

**使用GSEApy进行富集分析**

参考原论文Figure 5j，选择巨核细胞相关cluster进行GO富集分析：

```python
import gseapy as gp

# 选择巨核细胞cluster的marker基因（根据聚类注释）
megakaryocyte_clusters = ['1', '2']  # 根据实际聚类结果调整
all_markers = []
for cluster in megakaryocyte_clusters:
    markers = adata.uns['rank_genes_groups']['names'][cluster][:50]
    all_markers.extend(markers)

# 去除重复基因
all_markers = list(set(all_markers))

# GO富集分析（关注造血和巨核细胞分化相关通路）
enr = gp.enrichr(gene_list=all_markers,
                 gene_sets=['GO_Biological_Process_2023',
                           'KEGG_2021_Human',  # 包含造血相关通路
                           'Reactome_2022'],   # 包含血小板形成通路
                 organism='Human',
                 outdir='./enrichment_results',
                 cutoff=0.05)

# 筛选与巨核细胞分化相关的显著富集项
# 如："platelet activation", "megakaryocyte differentiation", "hemopoiesis"
```

**生物通路筛选**（基于原论文发现）：
```python
# 筛选原论文中提到的关键通路
key_terms = ['platelet', 'megakaryocyte', 'hemopoiesis', 'blood coagulation', 
           'erythroid', 'cell differentiation', 'FLI1']

# 自定义筛选与原论文一致的GO terms
filtered_results = enr.results[enr.results['Term'].str.contains('|'.join(key_terms), case=False)]

# 保存筛选后的结果
filtered_results.to_csv('figures/4a_filtered_go_results.csv')
```

**或使用Scanpy内置工具**
```python
# 需要先准备gene set数据库
from gseapy import Biomart

# 获取基因注释
gene_sets = sc.queries.biomart_annotations('hsapiens', 
                                           host='www.ensembl.org')
                                           
# 富集分析
sc.tl.filter_rank_genes_groups(adata)
# 导出基因列表用于在线工具分析（如DAVID, Enrichr, Metascape）
```

**可视化**
```python
# 使用gseapy绘图
from gseapy.plot import barplot, dotplot

# 条形图
barplot(enr.res2d, title='GO Enrichment', 
        cutoff=0.05, figsize=(10, 8))

# 点图
dotplot(enr.res2d, title='GO Enrichment',
        cmap='viridis_r', size=10, figsize=(10, 8))
```

**参考**：
- 原文Figure 5j
- GSEApy文档
- `scRNA_seq_analyze_totorial_2025.html` functional enrichment部分

### 4b. 结果解读要点
**结果解读要点**

参考原论文Figure 5j的解读：

- 选择显著的GO terms（adjusted p-value < 0.05）
- 重点解释与巨核细胞分化相关的生物学过程：
  - "Platelet activation"（血小板激活）
  - "Megakaryocyte differentiation"（巨核细胞分化）
  - "Hemopoiesis"（造血过程）
- 对比EB和iMK阶段的功能特征差异，反映重编程过程
- 注意关键转录因子（FLI1, MEIS1）调控的通路富集情况
- 讨论这些富集结果如何支持"红系到巨核系"的细胞命运转换

**可视化建议**（参考原论文Figure 5j）：
```python
from gseapy.plot import barplot, dotplot

# 条形图（最显著的前10个terms）
barplot(filtered_results.head(10), title='Top GO Terms in Megakaryocyte Clusters',
       figsize=(10, 8), ofname='figures/4b_go_barplot.png')

# 点图（展示富集强度和基因比例）
dotplot(filtered_results.head(15), title='GO Enrichment for Megakaryocyte Markers',
       figsize=(12, 10), ofname='figures/4b_go_dotplot.png',
       cmap='viridis_r')
```

---

## 任务5：进阶分析 (10分，任选一项)

### 选项a：Doublet去除

**Scrublet算法**
```python
import scrublet as scr

# 对每个样本分别运行
scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.06)
doublet_scores, predicted_doublets = scrub.scrub_doublets()

# 添加到adata
adata.obs['doublet_score'] = doublet_scores
adata.obs['predicted_doublet'] = predicted_doublets

# 可视化
sc.pl.umap(adata, color=['doublet_score', 'predicted_doublet'])

# 过滤doublets
adata_filtered = adata[~adata.obs['predicted_doublet']].copy()
```

**DoubletFinder（需要R接口）**
```python
# 或使用Python的doubletdetection包
import doubletdetection

clf = doubletdetection.BoostClassifier()
doublets = clf.fit(adata.X).predict()
adata.obs['doublet'] = doublets
```

**评估影响**：
- 比较去除前后的聚类数量和质量
- 检查是否有"中间态"cluster消失
- 对比marker基因表达的清晰度

**参考**：Scrublet文档，DoubletFinder论文

### 选项b：轨迹推断（Trajectory Inference）

**为什么需要轨迹推断**（根据原论文）：

本研究的核心是探索"红系细胞到巨核细胞的转分化过程"，这是一个**连续的细胞状态转换**：

- **生物学背景**：红细胞(EB)和巨核细胞起源于共同的祖细胞(MEP)，正常情况下沿不同路径分化
- **重编程过程**：小分子 cocktail 驱动红系细胞逆转分化路径，转向巨核细胞命运
- **轨迹分析价值**：
  1. 捕获从EB → iMKd3 → iMKd5 → iMKd7的连续状态转换
  2. 识别关键分支点（红系 vs 巨核系分化决策点）
  3. 发现驱动转换的关键转录因子（FLI1, MEIS1等）
  4. 量化重编程效率和程度
  5. 可视化动态基因表达变化，揭示分子机制

**使用PAGA进行轨迹推断**
```python
# PAGA (Partition-based graph abstraction)
# 注意：使用经过Harmony批处理校正后的数据
sc.tl.paga(adata_combined, groups='cell_type')
sc.pl.paga(adata_combined, plot=False)  # 移除未连接的边
sc.tl.umap(adata_combined, init_pos='paga')  # 使用PAGA初始化UMAP

# 可视化轨迹（参考原论文图）
sc.pl.paga(adata_combined, color=['cell_type', 'sample'], 
          threshold=0.03, fontsize=12, save='_trajectory_paga.png')

# 在UMAP上可视化轨迹
sc.pl.umap(adata_combined, color=['cell_type', 'dpt_pseudotime', 'sample'],
           save='_trajectory_umap.png')
```

**使用Monocle3（通过Python接口或R）**
```python
# 导出数据用于Monocle3 (R)
adata.write_h5ad('for_monocle.h5ad')
```

在R中：
```r
library(monocle3)
# 读取数据并分析
# 参考原文方法部分的Monocle3流程
```

**Diffusion pseudotime**
```python
# 使用扩散伪时间
adata.uns['iroot'] = np.flatnonzero(adata.obs['cell_type'] == 'EB')[0]
sc.tl.diffmap(adata)
sc.tl.dpt(adata)

# 可视化
sc.pl.umap(adata, color=['cell_type', 'dpt_pseudotime'])

# 沿伪时间的基因表达变化
sc.pl.umap(adata, color=['FLI1', 'MEIS1', 'dpt_pseudotime'])
```

**主要发现**（根据原论文）：
- **重编程轨迹**：EB → iMKd3 → iMKd5 → iMKd7 的连续转换路径
- **分化分支**：在某个时间点出现红系和巨核系两个潜在分化方向
- **关键转录因子动态**：
  - FLI1和MEIS1表达逐渐增加，驱动巨核细胞分化
  - 红系特异因子（如GATA1, KLF1）表达逐渐下降
- **细胞命运转换**：成功重编程的细胞表现出完全的巨核细胞基因表达特征
- **效率评估**：不同阶段中成功转换的细胞比例
- **分子机制**：小分子 cocktail 通过表观遗传重塑（ATAC数据支持）促进基因表达重编程

**参考**：
- 原文trajectory analysis部分
- PAGA论文
- Monocle3文档
- `scRNA_seq_analyze_totorial_2025.html` pseudotime analysis

### 选项c：自定义分析

**可能的创新方向**（结合原论文）：

1. **跨时间点差异表达分析**
```python
# 比较不同时间点间的基因表达变化
sc.tl.rank_genes_groups(adata_combined, 'sample', 
                        groups=['GSM6304416_iMK_D7'], 
                        reference='GSM6304413_EB_D0')
# 识别重编程过程中的关键上调/下调基因
```

2. **基因调控网络推断**
```python
# 使用SCENIC推断调控网络
# 识别驱动重编程的关键转录因子及其靶基因
# 重点关注FLI1和MEIS1调控的网络变化
```

3. **RNA velocity分析**
```python
import scvelo as scv
# 分析转录动态，预测细胞状态转换方向
# 验证从EB到iMK的转换轨迹
```

4. **细胞周期评分**
```python
# 评估重编程过程中的细胞周期变化
sc.tl.score_genes_cell_cycle(adata_combined, s_genes=S_genes, g2m_genes=G2M_genes)
```

5. **亚群精细分类**
```python
# 对巨核细胞亚群进行更精细的分类
# 区分早期和晚期巨核细胞
# 识别不完全重编程的中间态细胞
```

**评估标准**：
- 分析是否揭示原论文未充分探索的生物学见解
- 是否有助于理解红系到巨核系转换的分子机制
- 是否提供了验证或扩展原论文发现的新证据

---

## 代码组织建议

### 推荐的Jupyter Notebook结构

```python
# 1. 导入库和设置
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

# 2. 数据读取
# 读取四个样本并合并

# 3. 质量控制
# QC前可视化
# 过滤
# QC后可视化

# 4. 归一化和标准化

# 5. 批次效应校正

# 6. 降维和聚类

# 7. 细胞类型注释

# 8. 可视化

# 9. GO富集分析

# 10. 进阶分析（选一项）

# 11. 保存结果
adata.write('final_result.h5ad')
```

---

## 关键参数调整建议

### QC阈值
- `n_genes`: 200-2500（根据数据分布调整）
- `pct_counts_mt`: <5% 或 <10%（取决于组织类型）

### 聚类参数
- `n_neighbors`: 10-30（影响UMAP和聚类粒度）
- `resolution`: 0.4-1.2（控制cluster数量）

### 批次校正
- 选择适合数据规模的算法
- 验证是否过度校正（丢失生物学差异）

---

## 常用函数快速参考

### 数据预处理
- `sc.pp.filter_cells()` - 过滤细胞
- `sc.pp.filter_genes()` - 过滤基因
- `sc.pp.normalize_total()` - 归一化
- `sc.pp.log1p()` - log转换
- `sc.pp.highly_variable_genes()` - 高变基因
- `sc.pp.scale()` - 缩放

### 降维与聚类
- `sc.tl.pca()` - PCA降维
- `sc.pp.neighbors()` - 构建邻域图
- `sc.tl.umap()` - UMAP降维
- `sc.tl.tsne()` - t-SNE降维
- `sc.tl.leiden()` - Leiden聚类

### Marker基因分析
- `sc.tl.rank_genes_groups()` - 寻找marker基因
- `sc.tl.filter_rank_genes_groups()` - 过滤marker基因

### 可视化
- `sc.pl.umap()` - UMAP图
- `sc.pl.violin()` - 小提琴图
- `sc.pl.heatmap()` - 热图
- `sc.pl.rank_genes_groups_heatmap()` - Marker基因热图
- `sc.pl.dotplot()` - 点图

---

## 文档来源总结

1. **scRNA_seq_analyze_totorial_2025.html** - 主要教程文档
   - 完整的scanpy工作流程
   - QC、归一化、聚类的详细说明
   
2. **HW1.html** - 作业要求
   - 任务详细描述
   - 评分标准
   - 提交要求

3. **Scanpy官方文档** (https://scanpy.readthedocs.io/)
   - 所有函数的详细API文档
   - 完整的PBMC3k教程
   - 最佳实践指南

4. **原始论文**
   - 细胞类型定义
   - Marker基因
   - 生物学解释

---

## 注意事项

1. **阈值设定**：所有过滤阈值需根据实际数据分布调整，不要盲目使用教程中的值

2. **批次效应**：必须先验证是否存在批次效应，再选择合适的校正方法

3. **计算资源**：大规模数据处理需要足够的内存，建议至少16GB RAM

4. **重复性**：设置随机种子确保结果可重复
   ```python
   np.random.seed(0)
   sc.settings.seed = 0
   ```

5. **保存中间结果**：在关键步骤后保存adata对象
   ```python
   adata.write('checkpoint_qc.h5ad')
   ```

6. **参数优化**：聚类、UMAP等参数需要多次尝试才能获得最佳结果

7. **生物学验证**：结果要与原文对比，确保生物学意义合理

---

## 报告撰写提示

### 结构建议
1. 引言 - 简述研究背景和分析目的
2. 方法 - 详细描述每个分析步骤和参数选择
3. 结果 - 展示关键图表和发现
4. 讨论 - 解释结果的生物学意义，与原文对比
5. 结论 - 总结主要发现

### 图表要求
- 高分辨率（300 dpi以上）
- 清晰的标签和图例
- 合适的配色方案
- 每个图配有详细说明

### 篇幅控制
- 报告不超过10页A4纸
- 重点展示关键结果
- 代码放在附录或单独文件

---

## 推荐学习资源

1. **Scanpy教程**
   - PBMC 3k tutorial
   - PBMC 68k tutorial
   - 整合分析tutorial

2. **单细胞最佳实践**
   - https://www.sc-best-practices.org/

3. **相关综述**
   - Current best practices in single-cell RNA-seq analysis: a tutorial
   - Integrative single-cell analysis

4. **在线课程**
   - Analysis of single cell RNA-seq data (Hemberg Lab)
   - Orchestrating Single-Cell Analysis with Bioconductor
