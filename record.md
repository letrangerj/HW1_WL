# HW1

## 1. Preprocessing and Quality Control
### b. 线粒体比例的作用

在细胞死亡或出现破损时，随着细胞膜的破碎，细胞质中的mRNA会被释放出细胞而流失，但线粒体中的mRNA通常能够保留在线粒体结构内。因此低质量的细胞（死亡或破损细胞）会出现线粒体基因含量高的特征。

### c. QC前后的nFeature_RNA, nCount_RNA与percent.mt

指标的含义（前面是代码中使用的名称）：
- `n_genes_by_counts` / `nFeature_RNA`: 每个细胞检测到的基因数量
- `total_counts` / `nCount_RNA`: 每个细胞的总UMI/reads数
- `pct_counts_mt` / `percent.mt`: 线粒体基因表达百分比
