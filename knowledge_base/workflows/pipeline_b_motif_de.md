# Pipeline B: Motif → Differential Expression (Single-FOV)

## When to use
Use this pipeline after running motif enrichment (Pipeline A) to ask what genes are differentially expressed between cells defined by their spatial context. This pipeline extends Pipeline A by requiring `return_cellID=True`.

Typical user requests:
- "What genes are upregulated in T cells that are near macrophages?"
- "Compare T cells with vs. without this motif in their neighborhood"
- "What's different about macrophages that are close to tumor cells?"
- "Compare anchor cells between motif A and motif B"

**Critical**: always run `motif_enrichment_knn/dist` with `return_cellID=True` before any DE sub-pipeline.

---

## Prerequisite: run motif enrichment with return_cellID=True

```python
import anndata
import numpy as np
import pandas as pd
from SpatialQuery import spatial_query

adata = anndata.read_h5ad("path/to/data.h5ad")
sp = spatial_query(
    adata=adata,
    dataset='sample1',
    spatial_key='X_spatial',
    label_key='cell_type',
    feature_name='gene',
    build_gene_index=False,
    if_lognorm=True,
)

ct = 'T cell'  # anchor cell type

# IMPORTANT: return_cellID=True is required for DE
result = sp.motif_enrichment_knn(
    ct=ct,
    motifs=None,           # or specify e.g. ['Macrophage', 'B cell']
    k=30,
    min_support=0.5,
    max_dist=20,
    return_cellID=True,    # adds center_id and neighbor_id columns
)
# result columns include: motifs, if_significant, center_id, neighbor_id
# center_id: list of anchor cell indices that have this motif nearby
# neighbor_id: list of neighbor cell indices participating in this motif

sig = result[result['if_significant']]
print(f"Significant motifs: {len(sig)}")
```

---

## Sub-pipeline B1: Motif+ anchor vs. Motif− anchor

**Biological question**: What changes in the anchor cell type due to the presence of a specific motif in its neighborhood?

```python
# Select one significant motif row
row = sig.iloc[0]
motif_name = str(row['motifs'])  # for display

# Group 1: anchor cells that HAVE this motif nearby
ind_group1 = row['center_id']   # list of cell indices

# Group 2: all anchor cells MINUS those with the motif
# IMPORTANT: use np.setdiff1d, NOT set subtraction (shapes differ)
all_anchor_ids = np.where(np.array(adata.obs['cell_type']) == ct)[0]
ind_group2 = np.setdiff1d(all_anchor_ids, ind_group1)

print(f"Group 1 (motif+): {len(ind_group1)} cells")
print(f"Group 2 (motif−): {len(ind_group2)} cells")

de = sp.de_genes(
    ind_group1=ind_group1,
    ind_group2=ind_group2,
    method='wilcoxon',   # options: 'wilcoxon', 't-test', 'fisher'
    min_fraction=0.05,
)
# de columns: gene, proportion_1, proportion_2, p_value, adj-pval, de_in
sig_de = de[de['adj-pval'] < 0.05].sort_values('adj-pval')
print(f"Significant DE genes: {len(sig_de)}")
print(sig_de.head(20).to_string())
```

---

## Sub-pipeline B2: Motif-proximal vs. motif-distal neighbor cells

**Biological question**: What changes in a neighbor cell type (type A) due to spatial proximity to the anchor cell?

```python
# Identify the neighbor cell type to analyze
neighbor_ct = 'Macrophage'  # a non-anchor cell type present in the motif

row = sig.iloc[0]
neighbor_id = row['neighbor_id']  # all neighbor cells in this motif

# Group 1: cells of neighbor_ct that ARE in the motif's neighbor_id
all_neighbor_ct_ids = np.where(np.array(adata.obs['cell_type']) == neighbor_ct)[0]
ind_group1 = np.intersect1d(all_neighbor_ct_ids, neighbor_id)

# Group 2: cells of neighbor_ct that are NOT in neighbor_id
ind_group2 = np.setdiff1d(all_neighbor_ct_ids, neighbor_id)

print(f"Motif-proximal {neighbor_ct}: {len(ind_group1)} cells")
print(f"Motif-distal {neighbor_ct}: {len(ind_group2)} cells")

de = sp.de_genes(
    ind_group1=ind_group1,
    ind_group2=ind_group2,
    method='wilcoxon',
    min_fraction=0.05,
)
sig_de = de[de['adj-pval'] < 0.05].sort_values('adj-pval')
print(f"Significant DE genes: {len(sig_de)}")
print(sig_de.head(20).to_string())
```

*If the motif contains multiple non-anchor cell types, loop over each and run B2 separately for each.*

---

## Sub-pipeline B3: Motif-a+ anchor vs. Motif-b+ anchor

**Biological question**: What distinguishes anchor cells in different microenvironmental contexts (surrounded by motif A vs. surrounded by motif B)?

```python
# Select two significant motifs to compare
row_a = sig.iloc[0]  # motif A
row_b = sig.iloc[1]  # motif B

center_a = set(row_a['center_id'])
center_b = set(row_b['center_id'])

# Compute overlap
shared = center_a & center_b
print(f"Motif A cells: {len(center_a)}, Motif B cells: {len(center_b)}, Shared: {len(shared)}")
print(f"Overlap fraction (A): {len(shared)/len(center_a):.2%}")
print(f"Overlap fraction (B): {len(shared)/len(center_b):.2%}")

# Unique cells per motif
unique_a = np.array(list(center_a - shared))
unique_b = np.array(list(center_b - shared))

# Warn if either group is too small
MIN_CELLS = 10
if len(unique_a) < MIN_CELLS or len(unique_b) < MIN_CELLS:
    print(f"WARNING: unique cell counts are small (A={len(unique_a)}, B={len(unique_b)}). "
          f"Consider lowering min_support or relaxing parameters.")

# Only proceed if both groups are non-empty
if len(unique_a) > 0 and len(unique_b) > 0:
    de = sp.de_genes(
        ind_group1=unique_a,
        ind_group2=unique_b,
        method='wilcoxon',
        min_fraction=0.05,
    )
    sig_de = de[de['adj-pval'] < 0.05].sort_values('adj-pval')
    print(f"Significant DE genes: {len(sig_de)}")
    print(sig_de.head(20).to_string())
```

---

## Data flow summary

```
motif_enrichment_knn/dist(ct, return_cellID=True) → result
  └─ result[result['if_significant']] → sig

B1: sig.iloc[i]['center_id'] vs. all_anchor − center_id
      └─ de_genes(ind_group1, ind_group2) → DE table

B2: intersect(all_neighbor_ct, neighbor_id) vs. setdiff(all_neighbor_ct, neighbor_id)
      └─ de_genes(ind_group1, ind_group2) → DE table

B3: setdiff(center_id_a, shared) vs. setdiff(center_id_b, shared)
      └─ de_genes(ind_group1, ind_group2) → DE table
```

---

## Common failure modes

| Problem | Cause | Fix |
|---------|-------|-----|
| `ind_group2` is empty | All anchor cells have the motif | Check `n_center` vs `n_center_motif`; lower `min_support` |
| Too few cells in B3 unique groups | Motifs largely overlap | Report overlap fraction; suggest using B1 instead or picking different motifs |
| `de_genes` returns empty | No expressed genes pass `min_fraction` | Lower `min_fraction` to 0.01; check that adata.X is log-normalized |
| `neighbor_id` does not contain expected cell type | Motif query may return mixed types | Filter `neighbor_id` to the specific cell type using `np.intersect1d` |
| Wrong group sizes | Using `-` on arrays instead of `np.setdiff1d` | Always use `np.setdiff1d(all_ids, group_ids)` |
