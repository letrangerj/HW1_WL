# HW1

## 1. Preprocessing and Quality Control
### b. 线粒体比例的作用

在细胞死亡或出现破损时，随着细胞膜的破碎，细胞质中的mRNA会被释放出细胞而流失，但线粒体中的mRNA通常能够保留在线粒体结构内。因此低质量的细胞（死亡或破损细胞）会出现线粒体基因含量高的特征。

### c. QC前后的nFeature_RNA, nCount_RNA与percent.mt

指标的含义（前面是代码中使用的名称）：
- `n_genes_by_counts` / `nFeature_RNA`: 每个细胞检测到的基因数量
- `total_counts` / `nCount_RNA`: 每个细胞的总UMI/reads数
- `pct_counts_mt` / `percent.mt`: 线粒体基因表达百分比
在原文章中，methodology部分有：
“cells with a very small
library size (<2,500) and a very high (>0.15) mitochondrial genome transcript ratio were removed”

## 2. Batch Effect Removal
### a. 使用Harmony去除batch effect
文章中使用的方法为Mutual Nearest Neighbor (MNN)，这是由于文章的预处理流程依赖于Monocle3。这个方法会在局域中寻找最近邻，因此在处理不同batch间相似性/重合较高时会较有优势，并且会给出一个“修正后的表达矩阵”。

我希望使用课上教的Harmony算法，它与scanpy的流程是十分契合且匹配，并且速度上具有优势。Harmony在处理连续轨迹时能很好的保存这部分信息，这对于后续进行轨迹分析有优势。

## 3. Cell Clustering and Annotation
### a. 根据figure 5e-g:
得到marker genes如下：
```python
eryth_genes = ["GYPA", "KLF1", "EPB42"]
prog_genes = ["ANXA1", "CD44", "LMO4"]
mega_genes = ["ITGA3", "ITGA6", "GP1BA", "GP9", "F2R", "CD53"]
```
并且由文章可知：
```
  # 伪代码逻辑
if eryth_genes == "高":
    cluster = "C1"
elif eryth_genes == "中" and prog_genes == "上调":
    cluster = "C2"
elif prog_genes == "高":
    cluster = "C3"
elif mega_genes == "高" and eryth_genes == "低":
    cluster = "C4 or C5"  # 需要CD53进一步区分, C5的CD53高
```
leiden之后绘制不同的dotplot，将它们分配至不同的cell type完成了cell type annotation，如下