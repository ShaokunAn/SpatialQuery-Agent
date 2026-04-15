# Pipeline C: Motif → Gene Co-variation (Single-FOV)

## When to use
Use this pipeline after running motif enrichment (Pipeline A) to discover gene pairs that show coordinated expression between spatially proximal cells. This reveals cell-cell communication signatures and niche-specific gene regulatory programs.

Typical user requests:
- "What gene pairs are co-regulated between T cells and macrophages in this motif?"
- "Find ligand-receptor-like co-variation between anchor and neighbor cells"
- "Which genes in fibroblasts co-vary with neighboring tumor cells?"

Both C1 and C2 handle the neighbor/non-neighbor partitioning internally — you only need to provide `ct` and `motif`.

---

## Prerequisite: motif enrichment (return_cellID not required)

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

# Run enrichment to identify which motifs to analyze
result = sp.motif_enrichment_knn(ct=ct, motifs=None, k=30, min_support=0.5, return_cellID=False)
sig = result[result['if_significant']]
print(sig[['motifs', 'adj-pval']].to_string())
```

---

## Sub-pipeline C1: Pooled co-variation

**When to use**: motif contains only one non-anchor cell type, OR user explicitly asks for pooled analysis across all neighbor types.

**Biological question**: Which gene pairs between anchor cells and all motif neighbors show significantly different co-expression compared to non-neighbor pairs?

```python
# Select a significant motif
row = sig.iloc[0]
motif = list(row['motifs'])   # e.g. ['Macrophage', 'B cell']

results_df, ids = sp.compute_gene_gene_correlation(
    ct=ct,
    motif=motif,
    k=30,           # use k= OR max_dist=, not both
    min_size=0,
    min_nonzero=10,
    alpha=None,     # default 0.05
)
# results_df columns: gene_center, gene_motif, corr_neighbor, corr_non_neighbor,
#   p_value_test1, delta_corr_test1, corr_center_no_motif, p_value_test2,
#   delta_corr_test2, combined_score, q_value_test1, q_value_test2,
#   reject_test1_fdr, reject_test2_fdr, abs_combined_score, if_significant

sig_pairs = results_df[results_df['if_significant']]
print(f"Significant co-varying gene pairs: {len(sig_pairs)}")
print(sig_pairs[['gene_center', 'gene_motif', 'combined_score']].head(20).to_string())

# Visualize gene pair heatmap (biclustering)
cluster_df = sp.plot_gene_pair_heatmap(sig_pairs, figsize=(7, 5))

# Visualize spatial distribution of center/motif cells
sp.plot_all_center_motif(ct=ct, ids=ids, figsize=(6, 6))

# Extract cell IDs from ids if needed downstream
# ids['center_neighbor_motif_pair'] is ndarray shape (n_pairs, 2)
# pairs = ids['center_neighbor_motif_pair']
# center_ids = np.unique(pairs[:, 0])
# neighbor_ids = np.unique(pairs[:, 1])
```

---

## Sub-pipeline C2: Per-type co-variation (recommended default)

**When to use**: motif contains 2+ non-anchor cell types. Default recommendation — reveals type-specific interaction signatures that pooled analysis would average out.

**Biological question**: For each individual neighbor cell type in the motif, which gene pairs between anchor cells and that specific neighbor type show spatially-coordinated expression?

```python
row = sig.iloc[0]
motif = list(row['motifs'])   # e.g. ['Macrophage', 'B cell', 'NK cell']

# Note: motif must contain at least 2 non-anchor cell types for per-type analysis
# If motif has only 1 non-anchor type, C1 and C2 are equivalent

results_df = sp.compute_gene_gene_correlation_by_type(
    ct=ct,
    motif=motif,
    k=30,
    min_size=0,
    min_nonzero=10,
    alpha=None,
)
# Same columns as C1 plus: 'cell_type' (which neighbor type)

# Results per neighbor type
for neighbor_type in results_df['cell_type'].unique():
    sub = results_df[
        (results_df['cell_type'] == neighbor_type) &
        (results_df['if_significant'])
    ]
    print(f"\n{neighbor_type}: {len(sub)} significant gene pairs")
    print(sub[['gene_center', 'gene_motif', 'combined_score']].head(10).to_string())

# Visualize per type
sig_pairs = results_df[results_df['if_significant']]
cluster_df = sp.plot_gene_pair_heatmap(sig_pairs, figsize=(7, 5))
```

---

## Choosing C1 vs. C2

| Condition | Recommended |
|-----------|-------------|
| Motif has exactly 1 non-anchor cell type | C1 (C1 and C2 are equivalent) |
| Motif has 2+ non-anchor cell types | C2 (default) |
| User asks for "overall" or "pooled" | C1 |
| User asks "per cell type" or "type-specific" | C2 |

---

## Comparing co-variation between two conditions

If co-variation results from two conditions (e.g., treatment vs. control) are available, use `test_score_difference` to find gene pairs that differ:

```python
# result_A and result_B are DataFrames from compute_gene_gene_correlation
diff = spatial_query.test_score_difference(
    result_A, result_B,
    score_col='combined_score',
    significance_col='if_significant',
    gene_center_col='gene_center',
    gene_motif_col='gene_motif',
    percentile_threshold=95.0,
    background='Significant',
)
print(diff.head(20).to_string())
```

---

## Data flow summary

```
motif_enrichment_knn/dist(ct, return_cellID=False) → sig
  └─ motif = list(sig.iloc[i]['motifs'])

C1 (pooled):
  compute_gene_gene_correlation(ct, motif, k=30) → (results_df, ids)
    ├─ plot_gene_pair_heatmap(sig_pairs)
    └─ plot_all_center_motif(ct, ids)

C2 (per-type, default when motif has 2+ types):
  compute_gene_gene_correlation_by_type(ct, motif, k=30) → results_df
    └─ plot_gene_pair_heatmap(sig_pairs)

Cross-condition comparison:
  test_score_difference(result_A, result_B) → diff
```

---

## Common failure modes

| Problem | Cause | Fix |
|---------|-------|-----|
| Zero significant gene pairs | Sparse expression or weak spatial signal | Lower `alpha`; increase `k` or `max_dist`; lower `min_nonzero` |
| `compute_gene_gene_correlation_by_type` error | Motif has only 1 non-anchor type | Use C1 (`compute_gene_gene_correlation`) instead |
| Very few cells in `ids` | Motif is rare | Lower `min_support` in motif enrichment step |
| `plot_gene_pair_heatmap` fails | Empty `sig_pairs` DataFrame | Check `results_df['if_significant'].sum()` first |
| Cannot use both `k` and `max_dist` | API restriction | Pass only one: `k=30` OR `max_dist=20`, not both |
