# Cell Type Marker Analysis from Paper PMID: 35931032

## Paper Information
**Title:** Direct chemical reprogramming of human cord blood erythroblasts to induced megakaryocytes that produce platelets

**Journal:** Cell Stem Cell, 2022

**DOI:** 10.1016/j.stem.2022.07.004

**Authors:** Qin et al.

**Data Source:** GEO GSE207654

---

## Cell Types Identified in Figure 5B (scRNA-seq Analysis)

Based on the paper's Figure 5B and associated text, **6 clusters (C1-C6)** were identified through t-SNE and Leiden clustering of 24,967 individual cells from four time points (Day 0, 3, 5, 7):

---

## Cell Type Marker Table

| Cluster | Cell Type | Key Markers (High Expression) | Additional Markers | Biological Characteristics | Validation in HW1_指导文档.md |
|---------|-----------|------------------------------|-------------------|---------------------------|-------------------------------|
| **C1** | **EB (Erythroblasts)** | HBA1, GYPA, KLF1, EPB42, AHSP, HBB | TFRC, CD235a, HBA2 | Mainly day 0 cells; Active proliferation (G2M/S phase); Erythroid lineage | ✅ KLF1, GYPA, EPB42 mentioned |
| **C2** | **EB-like cells** | HBA1, GYPA, KLF1, EPB42, AHSP | Similar to C1 but gradually decreasing | Closely related to C1; Gradual decrease of erythroid genes; Active proliferation (G2M/S phase) | ✅ Markers valid |
| **C3** | **iPEMs (induced Precursors for Erythrocytes and Megakaryocytes)** | HPGDS, LMO4, KIT, CD44, ANXA1, IFITM3 | MYB (erythroid TF), FLI1 (MK TF) | Intermediate/bipotent state; Express both erythroid and MK markers; Transition cluster | ✅ LMO4, ANXA1, CD44 mentioned |
| **C4** | **iMK-1 (Thrombopoiesis-biased MKs)** | RGS18, F2R, PPBP, PF4, GP1BA, GP9, SELP, CD9, TUBB1, THBS1 | ITGA3, ITGA6, CD41a, CD42b, CD53- | Day 7 cells; Platelet aggregation/formation functions; Thrombopoiesis-biased; Majority of CD41a+CD42b+ cells | ✅ ITGA3, ITGA6, GP1BA, GP9, F2R, SELP valid |
| **C5** | **iMK-2 (Immune MKs)** | RGS18, F2R, PPBP, PF4, GP1BA, GP9, SELP, CD9, TUBB1, THBS1 | CD41a, CD42b, CD53+ | Day 7 cells; Immune/inflammatory functions; CD53+ distinguishes from iMK-1 | ✅ Same megakaryocytic markers |
| **C6** | **Lymphocytes** | (Lymphocyte-specific genes) | - | Removed from downstream analysis | ❌ Not relevant for analysis |

---

## Key Markers by Cell State/Lineage

### **1. Erythroid Lineage Markers** (C1, C2 - EBs and EB-like cells)

**Core markers:**
- HBA1, HBA2 (Hemoglobin alpha)
- HBB (Hemoglobin beta)
- GYPA (Glycophorin A)
- KLF1 (Kruppel-like factor 1)
- EPB42 (Erythrocyte membrane protein band 4.2)
- AHSP (Alpha hemoglobin stabilizing protein)

**Surface marker:**
- CD235a (erythroid-specific, also known as GYPA)

**Transcription factors:**
- KLF1 (erythroid master regulator)
- MYB (myeloblastosis oncogene)

**Additional:**
- TFRC (Transferrin receptor)

---

### **2. Bipotent Precursor Markers** (C3 - iPEMs)

**Precursor-associated markers:**
- HPGDS (Hematopoietic prostaglandin D synthase)
- LMO4 (LIM domain only 4)
- KIT (c-Kit, stem cell factor receptor)
- CD44 (Cell surface glycoprotein)
- ANXA1 (Annexin A1)
- IFITM3 (Interferon-induced transmembrane protein 3)

**Dual transcription factor expression:**
- MYB (erythroid lineage TF)
- FLI1 (megakaryocyte lineage TF)

**Biological significance:** This cluster represents the critical transition state where cells express markers from both erythroid and megakaryocyte lineages, indicating bipotent precursor characteristics.

---

### **3. Megakaryocyte Markers** (C4, C5 - iMK-1 and iMK-2)

**Core MK genes:**
- GP1BA (Glycoprotein Ib platelet subunit alpha, also CD42b)
- GP9 (Glycoprotein IX platelet)
- SELP (P-selectin, CD62P)
- CD9 (Tetraspanin)
- TUBB1 (Tubulin beta-1 chain)
- THBS1 (Thrombospondin 1)
- RGS18 (Regulator of G protein signaling 18)
- F2R (Protease-activated receptor 1)
- PPBP (Pro-platelet basic protein)
- PF4 (Platelet factor 4)

**Surface markers (for flow cytometry):**
- CD41a (ITGA2B, integrin alpha-IIb)
- CD42b (GP1BA)
- CD53 (iMK-2 specific, immune MK marker)

**Integrins:**
- ITGA3 (Integrin alpha-3)
- ITGA6 (Integrin alpha-6)
- ITGA2B (CD41a)
- ITGB3 (Integrin beta-3)

**Transcription factors:**
- FLI1 (Friend leukemia integration 1)
- MEIS1 (Meis homeobox 1)
- GATA2 (GATA binding protein 2)
- FOG1 (Friend of GATA 1)
- GABPA (GA binding protein transcription factor subunit alpha)
- RUNX1 (Runt-related transcription factor 1)

**Platelet/secretion markers:**
- PF4 (Platelet factor 4)
- THBS1 (Thrombospondin 1)
- vWF (von Willebrand factor)
- PLEK (Pleckstrin)

---

### **4. iMK Subtype Distinguishing Markers**

**iMK-1 (Thrombopoiesis-biased MKs):**
- Phenotype: CD41a+ CD42b+ CD53-
- Function: Platelet aggregation and formation
- GO enrichment: PLT aggregation, PLT formation

**iMK-2 (Immune MKs):**
- Phenotype: CD41a+ CD42b+ CD53+
- Function: Immune and inflammatory responses
- CD53 is the key distinguishing surface marker

---

## Validation Against HW1_指导文档.md

The guidance document (lines 235-270) mentions some markers but doesn't provide complete lists. Based on the paper analysis:

### ✅ **VALID markers mentioned in guidance:**

**Erythroid:**
- KLF1, EPB42, GYPA ✓

**Megakaryocyte:**
- ITGA3, ITGA6, GP1BA, GP9, F2R, SELP ✓

**Precursor:**
- LMO4, ANXA1, CD44 ✓

### 📝 **Additional IMPORTANT markers from paper NOT in guidance:**

**Erythroid (critical markers missing):**
- HBA1, HBA2, AHSP (very important for erythroid identification)
- HBB (hemoglobin beta)

**Precursor (missing):**
- HPGDS, KIT, IFITM3

**Megakaryocyte (critical markers missing):**
- RGS18, PPBP, PF4, CD9, TUBB1, THBS1 (essential for MK identification)
- CD53 (critical for distinguishing iMK subtypes)

**Transcription factors (not mentioned but critical):**
- FLI1, MEIS1 (key TFs driving megakaryopoiesis)
- MYB (erythroid TF)

---

## How to Use These Markers for t-SNE/UMAP Plotting with Cell Type Annotations

### **Method 1: Visualizing Individual Marker Genes**

```python
import scanpy as sc
import matplotlib.pyplot as plt

# After preprocessing, normalization, and dimensionality reduction

# Visualize erythroid markers on UMAP/t-SNE
erythroid_markers = ['HBA1', 'GYPA', 'KLF1', 'EPB42', 'AHSP']
sc.pl.umap(adata, color=erythroid_markers, cmap='Reds', ncols=3)
sc.pl.tsne(adata, color=erythroid_markers, cmap='Reds', ncols=3)

# Visualize precursor markers
precursor_markers = ['LMO4', 'KIT', 'CD44', 'HPGDS', 'ANXA1']
sc.pl.umap(adata, color=precursor_markers, cmap='Greens', ncols=3)

# Visualize megakaryocyte markers
mk_markers = ['GP1BA', 'GP9', 'PF4', 'SELP', 'THBS1', 'RGS18']
sc.pl.umap(adata, color=mk_markers, cmap='Blues', ncols=3)

# Visualize key transcription factors
tfs = ['KLF1', 'MYB', 'FLI1', 'MEIS1']
sc.pl.umap(adata, color=tfs, cmap='Purples', ncols=2)
```

---

### **Method 2: Gene Signature Scoring for Cell Type Annotation**

```python
# Define marker gene sets for each cell type
marker_genes = {
    'EB': ['HBA1', 'HBA2', 'GYPA', 'KLF1', 'EPB42', 'AHSP', 'HBB'],
    'EB-like': ['HBA1', 'GYPA', 'KLF1', 'EPB42'],
    'iPEM': ['HPGDS', 'LMO4', 'KIT', 'CD44', 'ANXA1', 'IFITM3'],
    'iMK': ['GP1BA', 'GP9', 'SELP', 'PPBP', 'PF4', 'THBS1', 'RGS18', 'F2R', 'CD9', 'TUBB1']
}

# Calculate gene signature scores for each cell type
for cell_type, genes in marker_genes.items():
    # Filter genes present in dataset
    genes_use = [g for g in genes if g in adata.var_names]
    sc.tl.score_genes(adata, genes_use, score_name=f'{cell_type}_score')

# Visualize scores on UMAP
sc.pl.umap(adata, color=['EB_score', 'EB-like_score', 'iPEM_score', 'iMK_score'], 
           cmap='viridis', ncols=2, vmin=0)

# Assign cell types based on highest score
import numpy as np
score_cols = ['EB_score', 'EB-like_score', 'iPEM_score', 'iMK_score']
scores = adata.obs[score_cols].values
adata.obs['predicted_celltype'] = [score_cols[i].replace('_score', '') 
                                     for i in np.argmax(scores, axis=1)]

# Plot predicted cell types
sc.pl.umap(adata, color='predicted_celltype', legend_loc='on data')
sc.pl.tsne(adata, color='predicted_celltype', legend_loc='on data')
```

---

### **Method 3: Cluster-based Annotation Using Marker Expression**

```python
# Perform clustering
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata, resolution=0.5, key_added='leiden')

# Find marker genes for each cluster
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# Check marker expression in each cluster
sc.pl.dotplot(adata, marker_genes, groupby='leiden', dendrogram=True)

# Manual annotation based on marker expression
# Check top markers for each cluster and compare with known markers
cluster_annotations = {
    '0': 'EB',
    '1': 'EB-like',
    '2': 'iPEM',
    '3': 'iMK-1',
    '4': 'iMK-2'
}

# Apply annotations (adjust based on your actual cluster results)
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_annotations)

# Visualize annotated cell types
sc.pl.umap(adata, color='cell_type', legend_loc='on data', 
           palette={'EB': '#d62728', 'EB-like': '#ff7f0e', 
                   'iPEM': '#2ca02c', 'iMK-1': '#1f77b4', 'iMK-2': '#9467bd'})

sc.pl.tsne(adata, color='cell_type', legend_loc='on data', 
           palette={'EB': '#d62728', 'EB-like': '#ff7f0e', 
                   'iPEM': '#2ca02c', 'iMK-1': '#1f77b4', 'iMK-2': '#9467bd'})
```

---

### **Method 4: Comprehensive Marker Visualization**

```python
# Create dotplot showing all key markers across cell types
all_markers = {
    'Erythroid': ['HBA1', 'GYPA', 'KLF1', 'EPB42', 'AHSP', 'HBB'],
    'Precursor': ['HPGDS', 'LMO4', 'KIT', 'CD44', 'ANXA1'],
    'Megakaryocyte': ['GP1BA', 'GP9', 'SELP', 'PPBP', 'PF4', 'THBS1', 
                     'RGS18', 'F2R', 'CD9', 'TUBB1'],
    'TFs': ['KLF1', 'MYB', 'FLI1', 'MEIS1']
}

# Flatten marker dictionary
marker_list = []
for category, genes in all_markers.items():
    marker_list.extend(genes)

# Remove duplicates while preserving order
seen = set()
marker_list_unique = [x for x in marker_list if not (x in seen or seen.add(x))]

# Dotplot showing marker expression across cell types
sc.pl.dotplot(adata, marker_list_unique, groupby='cell_type', 
              dendrogram=True, standard_scale='var')

# Heatmap of top markers
sc.pl.heatmap(adata, marker_list_unique, groupby='cell_type', 
              swap_axes=True, cmap='RdBu_r', vmin=-2, vmax=2, 
              dendrogram=True, figsize=(8, 10))

# Stacked violin plot
sc.pl.stacked_violin(adata, marker_list_unique, groupby='cell_type', 
                     swap_axes=False, dendrogram=True)

# Matrix plot
sc.pl.matrixplot(adata, marker_list_unique, groupby='cell_type', 
                 dendrogram=True, cmap='viridis', standard_scale='var')
```

---

### **Method 5: Integration with Time Point Information**

```python
# Combine cell type annotation with time point information
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Plot colored by cell type
sc.pl.umap(adata, color='cell_type', ax=axes[0], show=False, 
           title='Cell Type Annotation')

# Plot colored by time point
sc.pl.umap(adata, color='sample', ax=axes[1], show=False, 
           title='Time Point')
plt.tight_layout()
plt.show()

# Show distribution of cell types across time points
import pandas as pd

# Create crosstab
cell_type_time = pd.crosstab(adata.obs['cell_type'], 
                             adata.obs['sample'], 
                             normalize='columns') * 100

# Plot stacked bar chart
ax = cell_type_time.T.plot(kind='bar', stacked=True, 
                           figsize=(10, 6), 
                           colormap='Set3')
plt.ylabel('Percentage of cells (%)')
plt.xlabel('Time point')
plt.title('Cell type composition across reprogramming time course')
plt.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# Split UMAP by time point
sc.pl.umap(adata, color='cell_type', split_show=True, 
           groups=['GSM6304413_EB_D0', 'GSM6304414_iMK_D3', 
                  'GSM6304415_iMK_D5', 'GSM6304416_iMK_D7'],
           legend_loc='on data', ncols=2)
```

---

### **Method 6: Distinguishing iMK-1 and iMK-2 Subtypes**

```python
# For cells annotated as iMK, further classify into iMK-1 and iMK-2
# Based on CD53 expression (CD53- = iMK-1, CD53+ = iMK-2)

# Check if CD53 is in the dataset
if 'CD53' in adata.var_names:
    # Get iMK cells
    imk_mask = adata.obs['cell_type'].isin(['iMK-1', 'iMK-2'])
    
    # Calculate CD53 expression threshold (e.g., median)
    cd53_threshold = adata[imk_mask, 'CD53'].X.toarray().flatten().median()
    
    # Classify based on CD53 expression
    adata.obs.loc[imk_mask, 'iMK_subtype'] = ['iMK-2' if x > cd53_threshold 
                                                else 'iMK-1' 
                                                for x in adata[imk_mask, 'CD53'].X.toarray().flatten()]
    
    # Visualize
    sc.pl.umap(adata, color=['iMK_subtype', 'CD53'], 
               palette={'iMK-1': '#1f77b4', 'iMK-2': '#9467bd'})
    
    # Compare marker expression between iMK-1 and iMK-2
    imk_markers = ['GP1BA', 'GP9', 'PF4', 'THBS1', 'CD53']
    sc.pl.violin(adata[imk_mask], imk_markers, groupby='iMK_subtype')
```

---

## Practical Workflow for Cell Type Annotation

### **Step-by-step approach:**

1. **Quality control and preprocessing** (as per HW1 guidance)

2. **Dimensionality reduction**
   ```python
   sc.pp.highly_variable_genes(adata, n_top_genes=2000)
   sc.pp.scale(adata, max_value=10)
   sc.tl.pca(adata, svd_solver='arpack')
   sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
   sc.tl.umap(adata)
   sc.tl.tsne(adata)
   ```

3. **Clustering**
   ```python
   sc.tl.leiden(adata, resolution=0.5)
   ```

4. **Visualize marker expression**
   ```python
   # Check all key markers
   sc.pl.umap(adata, color=['leiden', 'HBA1', 'GYPA', 'LMO4', 
                           'GP1BA', 'PF4', 'FLI1'])
   ```

5. **Find cluster-specific markers**
   ```python
   sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
   sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, groupby='leiden')
   ```

6. **Annotate clusters** based on marker expression patterns

7. **Validate annotations**
   ```python
   # Compare with time point
   sc.pl.umap(adata, color=['cell_type', 'sample'])
   
   # Check marker expression in annotated types
   sc.pl.dotplot(adata, marker_list_unique, groupby='cell_type')
   ```

---

## Expected Cell Type Distribution Across Time Points

Based on Figure 5C from the paper:

| Cell Type | Day 0 (EB) | Day 3 (iMK-D3) | Day 5 (iMK-D5) | Day 7 (iMK-D7) |
|-----------|------------|----------------|----------------|----------------|
| **C1 (EB)** | High (~90%) | Low | Very low | Minimal |
| **C2 (EB-like)** | Low | Medium | Low | Minimal |
| **C3 (iPEM)** | Minimal | Medium | High | Medium |
| **C4 (iMK-1)** | Minimal | Low | Medium | High (~60-70%) |
| **C5 (iMK-2)** | Minimal | Low | Low | Medium (~20-30%) |

This temporal pattern reflects the reprogramming trajectory: EB → EB-like → iPEM → iMK

---

## Key Transcription Factor Changes During Reprogramming

The paper emphasizes critical transcription factor transitions:

**Erythroid TFs (downregulated during reprogramming):**
- KLF1, KLF3, KLF4, KLF5, KLF6
- MYB

**Megakaryocyte TFs (upregulated during reprogramming):**
- FLI1 (most critical)
- MEIS1 (most critical)
- GATA2
- FOG1
- GABPA
- RUNX1

**Monitoring these TFs is crucial for validating successful reprogramming.**

---

## Summary and Recommendations

### **Cell Types in the Study:**

1. **C1 (EB):** Initial erythroblasts - express high levels of HBA1, GYPA, KLF1
2. **C2 (EB-like):** Transitioning erythroblasts - erythroid markers gradually decrease
3. **C3 (iPEM):** Bipotent precursors - express both MYB and FLI1
4. **C4 (iMK-1):** Thrombopoiesis-biased megakaryocytes - CD53-, high platelet markers
5. **C5 (iMK-2):** Immune megakaryocytes - CD53+, immune-related genes

### **Markers Validation:**

The markers mentioned in **HW1_指导文档.md are valid but incomplete**. To ensure robust cell type annotation, you should use:

**Minimum marker set for annotation:**
- **Erythroid:** HBA1, GYPA, KLF1, EPB42
- **Precursor:** LMO4, KIT, CD44
- **Megakaryocyte:** GP1BA, GP9, PF4, SELP, THBS1

**Recommended comprehensive marker set:**
- **Erythroid:** HBA1, HBA2, GYPA, KLF1, EPB42, AHSP, HBB
- **Precursor:** HPGDS, LMO4, KIT, CD44, ANXA1
- **Megakaryocyte:** GP1BA, GP9, SELP, PPBP, PF4, THBS1, RGS18, F2R
- **TFs:** MYB, KLF1, FLI1, MEIS1

### **For t-SNE/UMAP Visualization:**

1. Use gene signature scores for initial annotation
2. Perform clustering and find cluster-specific markers
3. Manually annotate clusters based on marker expression
4. Visualize with dotplots, heatmaps, and violin plots
5. Integrate time point information to validate reprogramming trajectory

---

## References

1. **Original Paper:** Qin et al. (2022). Direct chemical reprogramming of human cord blood erythroblasts to induced megakaryocytes that produce platelets. *Cell Stem Cell*, 29(8), 1229-1245.e7. PMID: 35931032

2. **GEO Data:** GSE207654

3. **Key Methods:** t-SNE, Leiden clustering, scRNA-seq analysis with Scanpy

---

*Document created: 2025-11-11*
*For: Tumor Computation HW1 - Single-cell RNA-seq analysis*
