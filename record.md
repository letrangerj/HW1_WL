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


## 3. Cell Clustering and Annotation
### a. 根据figure 5e-g:
得到marker genes如下：
```python
megakaryocyte_markers = ["ITGA3", "ITGA6", "GP1BA", "GP9", "F2R", "SELP"]
erythroid_markers = ["GYPA", "KLF1", "EPB42"]
progenitor_markers = ["ANXA1", "CD44", "LMO4"]
```
并且有
```
  # 伪代码逻辑
if erythroid_markers == "高" and megakaryocyte_markers == "低":
    cluster = "C1"
    
elif erythroid_markers == "中" and progenitor_markers == "上调":
    cluster = "C2"
    
elif progenitor_markers == "高" and erythroid_markers == "降低" and megakaryocyte_markers == "开始上调":
    cluster = "C3"
    
elif megakaryocyte_markers == "高" and erythroid_markers == "低":
    cluster = "C4 or C5"  # 需要CD53进一步区分, C5的CD53高
```
