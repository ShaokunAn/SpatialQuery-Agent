# Multi-FOV Pipelines

## When to use
Use `spatial_query_multi` when the analysis involves multiple tissue sections, samples, or experimental conditions. All single-FOV pipelines (A, B, C) have multi-FOV equivalents. Key advantages:
- Aggregate motif enrichment across all FOVs for higher statistical power
- Compare motif frequencies between two groups of samples (differential analysis)
- Run DE comparisons across conditions using a unified cell index format

Typical user requests:
- "Compare neighborhoods between healthy and disease samples"
- "Run motif enrichment across all my samples"
- "What motifs are more enriched in condition A vs condition B?"

---

## Initialization and pre-processing

```python
import anndata
import numpy as np
import pandas as pd
from SpatialQuery import spatial_query_multi

# Load multiple AnnData objects
adata1 = anndata.read_h5ad("sample1.h5ad")
adata2 = anndata.read_h5ad("sample2.h5ad")
# ... add more as needed

# Inspect keys (should be consistent across samples)
print(adata1.obsm.keys())
print(adata1.obs.columns.tolist())
print(adata1.obs['cell_type'].unique().tolist())

# Initialize multi-FOV object
sq_multi = spatial_query_multi(
    adatas=[adata1, adata2],
    datasets=['healthy1', 'healthy2'],   # underscores replaced with hyphens internally
    spatial_key='X_spatial',
    label_key='cell_type',
    feature_name='gene',
    build_gene_index=False,
    if_lognorm=True,
)
# Internally, datasets become: 'healthy1_0', 'healthy1_1', etc. if duplicated names
# Access FOV names via sq_multi.datasets
print(sq_multi.datasets)
```

---

## Pre-processing visualization

Visualize cell type composition across samples before running analysis.

```python
# Proportion of each cell type across all datasets (stacked bar or similar)
sq_multi.plot_cell_type_distribution(
    dataset=None,          # None = all datasets
    data_type='proportion', # 'proportion' or 'number'
    colormap='tab20c',
)

# Per-FOV cell type distribution for a specific dataset group
sq_multi.plot_cell_type_distribution_fov(
    dataset='healthy1',    # the dataset name prefix
    data_type='proportion',
    colormap='tab20c',
)
```

---

## Multi-FOV Pipeline A: Motif Discovery

Aggregate motif enrichment across all FOVs or a specified subset.

```python
ct = 'T cell'   # anchor cell type

# Aggregate across all FOVs (recommended starting point)
result = sq_multi.motif_enrichment_knn(
    ct=ct,
    motifs=None,       # None = auto-discover
    dataset=None,      # None = use all datasets
    k=30,
    min_support=0.5,
    max_dist=20,
)
# OR distance-based:
result = sq_multi.motif_enrichment_dist(
    ct=ct,
    motifs=None,
    dataset=None,
    max_dist=20,
    min_size=0,
    min_support=0.5,
)

sig = result[result['if_significant']]
print(f"Significant motifs: {len(sig)}")

# Visualize
sq_multi.plot_motif_enrichment_heatmap(sig, figsize=(7, 5),
                                        title=f'Motif enrichment around {ct}')
```

To restrict to specific datasets:
```python
result_cond_a = sq_multi.motif_enrichment_knn(ct=ct, dataset='healthy', k=30)
result_cond_b = sq_multi.motif_enrichment_knn(ct=ct, dataset='disease', k=30)
```

---

## Multi-FOV Pipeline B: Motif → DE

`sq_multi.de_genes` accepts cell groups as `Dict[str, List[int]]` keyed by dataset name (using the internal `name_N` format from `sq_multi.datasets`).

### B1: Motif+ anchor vs. Motif− anchor (within-dataset)

Same logic as single-FOV B1, repeated per dataset and packed into dictionaries.

```python
result = sq_multi.motif_enrichment_knn(
    ct=ct, motifs=None, dataset=None, k=30, min_support=0.5, return_cellID=True
)
sig = result[result['if_significant']]
row = sig.iloc[0]

# center_id and neighbor_id may be dicts keyed by dataset name, or lists
# Check the structure:
print(type(row['center_id']))

# Build ind_group dicts (structure depends on motif_enrichment return format)
# For each dataset, separate motif+ and motif- anchor cells
ind_group1 = {}  # motif+ anchor per dataset
ind_group2 = {}  # motif- anchor per dataset

for ds in sq_multi.datasets:
    # Get all anchor cells in this dataset
    sp_fov = sq_multi.spatial_queries[ds]  # individual spatial_query object
    all_anchor = np.where(np.array(sp_fov.labels) == ct)[0]
    motif_pos = np.array(row['center_id'].get(ds, []))
    motif_neg = np.setdiff1d(all_anchor, motif_pos)
    if len(motif_pos) > 0:
        ind_group1[ds] = motif_pos.tolist()
        ind_group2[ds] = motif_neg.tolist()

de = sq_multi.de_genes(
    ind_group1=ind_group1,
    ind_group2=ind_group2,
    method='wilcoxon',
    min_fraction=0.05,
)
sig_de = de[de['adj-pval'] < 0.05].sort_values('adj-pval')
print(sig_de.head(20).to_string())
```

### B2: Motif+ anchor in condition A vs. condition B (cross-dataset)

Compare motif-positive anchor cells between two conditions.

```python
# Assume: condition A = ['healthy1_0', 'healthy2_0'], condition B = ['disease1_0', 'disease2_0']
cond_a_datasets = [ds for ds in sq_multi.datasets if 'healthy' in ds]
cond_b_datasets = [ds for ds in sq_multi.datasets if 'disease' in ds]

# Use aggregated enrichment per condition
result_a = sq_multi.motif_enrichment_knn(ct=ct, dataset='healthy', k=30, return_cellID=True)
result_b = sq_multi.motif_enrichment_knn(ct=ct, dataset='disease', k=30, return_cellID=True)

sig_a = result_a[result_a['if_significant']]
sig_b = result_b[result_b['if_significant']]

# For a shared motif, gather motif+ anchor cells from each condition
target_motif = list(sig_a.iloc[0]['motifs'])

row_a = sig_a[sig_a['motifs'].apply(lambda x: set(x) == set(target_motif))].iloc[0]
row_b = sig_b[sig_b['motifs'].apply(lambda x: set(x) == set(target_motif))].iloc[0]

ind_group1 = {ds: ids for ds, ids in row_a['center_id'].items()}  # condition A
ind_group2 = {ds: ids for ds, ids in row_b['center_id'].items()}  # condition B

de = sq_multi.de_genes(ind_group1=ind_group1, ind_group2=ind_group2, method='wilcoxon')
sig_de = de[de['adj-pval'] < 0.05].sort_values('adj-pval')
print(sig_de.head(20).to_string())
```

### B3: Motif-proximal vs. motif-distal neighbor cells (within-dataset)

Single-FOV B2 logic, applied per dataset.

```python
neighbor_ct = 'Macrophage'
row = sig.iloc[0]

ind_group1 = {}  # proximal neighbor cells per dataset
ind_group2 = {}  # distal neighbor cells per dataset

for ds in sq_multi.datasets:
    sp_fov = sq_multi.spatial_queries[ds]
    all_neighbor_ct = np.where(np.array(sp_fov.labels) == neighbor_ct)[0]
    proximal = np.array(row['neighbor_id'].get(ds, []))
    proximal_of_type = np.intersect1d(all_neighbor_ct, proximal)
    distal_of_type = np.setdiff1d(all_neighbor_ct, proximal)
    if len(proximal_of_type) > 0:
        ind_group1[ds] = proximal_of_type.tolist()
        ind_group2[ds] = distal_of_type.tolist()

de = sq_multi.de_genes(ind_group1=ind_group1, ind_group2=ind_group2, method='wilcoxon')
sig_de = de[de['adj-pval'] < 0.05].sort_values('adj-pval')
print(sig_de.head(20).to_string())
```

### B4: Motif neighbor cells in condition A vs. condition B (cross-dataset)

Compare motif-participating neighbor cells across conditions.

```python
ind_group1 = {}   # neighbor cells from condition A
ind_group2 = {}   # neighbor cells from condition B

for ds in cond_a_datasets:
    neighbor_ids = np.array(row_a['neighbor_id'].get(ds, []))
    sp_fov = sq_multi.spatial_queries[ds]
    of_type = np.intersect1d(
        np.where(np.array(sp_fov.labels) == neighbor_ct)[0],
        neighbor_ids
    )
    if len(of_type) > 0:
        ind_group1[ds] = of_type.tolist()

for ds in cond_b_datasets:
    neighbor_ids = np.array(row_b['neighbor_id'].get(ds, []))
    sp_fov = sq_multi.spatial_queries[ds]
    of_type = np.intersect1d(
        np.where(np.array(sp_fov.labels) == neighbor_ct)[0],
        neighbor_ids
    )
    if len(of_type) > 0:
        ind_group2[ds] = of_type.tolist()

de = sq_multi.de_genes(ind_group1=ind_group1, ind_group2=ind_group2, method='wilcoxon')
sig_de = de[de['adj-pval'] < 0.05].sort_values('adj-pval')
print(sig_de.head(20).to_string())
```

**Note**: B1–B4 are representative scenarios, not exhaustive. `sq_multi.de_genes` accepts any two cell groups as dicts — users can construct arbitrary comparisons.

---

## Multi-FOV Pipeline C: Differential Motif Analysis

Identify which motifs are enriched in one group of samples vs. another.

```python
from SpatialQuery.plotting import plot_differential_pattern_heatmap

# datasets parameter: list of two lists of dataset names
# Group 1: condition A; Group 2: condition B
datasets = [cond_a_datasets, cond_b_datasets]

# Or pass dataset name strings directly — check the actual API signature
diff = sq_multi.differential_analysis_knn(
    ct=ct,
    datasets=datasets,
    motifs=None,        # None = auto-discover; or specify list
    k=30,
    min_support=0.5,
    max_dist=20,
)
# diff: Dict[str, DataFrame] keyed by dataset name
# Each DataFrame columns: itemsets, support_{ds0}_mean, support_{ds1}_mean, adj-pval

# OR distance-based:
diff = sq_multi.differential_analysis_dist(
    ct=ct,
    datasets=datasets,
    motifs=None,
    max_dist=20,
    min_support=0.5,
    min_size=0,
)

# Visualize
plot_differential_pattern_heatmap(diff, ct=ct, figsize=None, cmap='YlOrRd')
```

---

## Multi-FOV Pipeline D: Gene Co-variation

Structure mirrors single-FOV Pipeline C (D1 = pooled, D2 = per-type), with two key differences:
- The second return value of `compute_gene_gene_correlation` is `fov_info` (FOV-level statistics), not `ids`
- `compute_gene_gene_correlation_by_type` returns only a `DataFrame` (no second element)
- Column names for FDR differ slightly (see table below)

### Column name reference

| Column | Pooled (D1) | Per-type (D2) |
|--------|-------------|---------------|
| FDR test1 | `adj_p_value_test1` | `q_value_test1` + `reject_test1_fdr` |
| FDR test2 | `adj_p_value_test2` | `q_value_test2` + `reject_test2_fdr` |
| Significance flag | `if_significant` | `if_significant` (= `reject_test1_fdr & reject_test2_fdr`) |
| Combined score | `combined_score` | `combined_score` |

---

### Sub-pipeline D1: Pooled co-variation

**When to use**: motif has one non-anchor cell type, or user asks for pooled analysis.

```python
from SpatialQuery import spatial_query_multi

# Aggregate across all FOVs (dataset=None)
results_df, fov_info = sq_multi.compute_gene_gene_correlation(
    ct=ct,
    motif=['Macrophage'],     # can be single type or list
    dataset=None,             # None = all datasets; or 'condition_A' / list of dataset names
    k=30,                     # use k= OR max_dist=, not both
    min_size=0,
    min_nonzero=10,
    alpha=None,
)
# results_df columns: gene_center, gene_motif, corr_neighbor, corr_non_neighbor,
#   corr_center_no_motif, p_value_test1, p_value_test2, delta_corr_test1,
#   delta_corr_test2, combined_score, adj_p_value_test1, adj_p_value_test2, if_significant
#
# fov_info dict keys: fov_statistics, total_pairs_neighbor, total_pairs_non_neighbor,
#   total_pairs_no_motif, n_fovs_analyzed

print(f"FOVs analyzed: {fov_info['n_fovs_analyzed']}")
print(f"Total neighbor pairs: {fov_info['total_pairs_neighbor']}")

sig_pairs = results_df[results_df['if_significant']]
print(f"Significant co-varying gene pairs: {len(sig_pairs)}")
print(sig_pairs[['gene_center', 'gene_motif', 'combined_score']].head(20).to_string())

# Visualize
sp.plot_gene_pair_heatmap(sig_pairs, figsize=(7, 5))
```

---

### Sub-pipeline D2: Per-type co-variation (recommended default)

**When to use**: motif has 2+ non-anchor cell types. Reveals type-specific interaction signatures.

```python
results_df = sq_multi.compute_gene_gene_correlation_by_type(
    ct=ct,
    motif=['Macrophage', 'B cell', 'NK cell'],
    dataset=None,
    k=30,
    min_nonzero=10,
    alpha=None,
)
# results_df columns: cell_type, gene_center, gene_motif, corr_neighbor, corr_non_neighbor,
#   corr_center_no_motif, p_value_test1, p_value_test2, delta_corr_test1, delta_corr_test2,
#   combined_score, abs_combined_score, q_value_test1, q_value_test2,
#   reject_test1_fdr, reject_test2_fdr,
#   if_significant (= reject_test1_fdr & reject_test2_fdr)

# Filter significant pairs per neighbor type
for neighbor_type in results_df['cell_type'].unique():
    sub = results_df[
        (results_df['cell_type'] == neighbor_type) &
        (results_df['if_significant'])
    ]
    print(f"\n{neighbor_type}: {len(sub)} significant gene pairs")
    print(sub[['gene_center', 'gene_motif', 'combined_score']].head(10).to_string())

# Visualize all significant pairs together
sig_pairs = results_df[results_df['if_significant']]
sp.plot_gene_pair_heatmap(sig_pairs, figsize=(7, 5))
```

---

### D3: Cross-condition comparison with `test_score_difference`

Identify gene pairs with significantly different co-variation scores between two conditions. Results are ranked by score difference; group-enriched pairs are extracted via `outlier_direction`.

```python
# Step 1: Compute co-variation separately per condition
results_cond_a, fov_info_a = sq_multi.compute_gene_gene_correlation(
    ct=ct,
    motif=['Macrophage'],
    dataset='condition_A',   # restrict to condition A FOVs
    k=30,
    min_nonzero=10,
)
results_cond_b, fov_info_b = sq_multi.compute_gene_gene_correlation(
    ct=ct,
    motif=['Macrophage'],
    dataset='condition_B',
    k=30,
    min_nonzero=10,
)

# Step 2: Compare scores between conditions
diff = spatial_query_multi.test_score_difference(
    result_A=results_cond_a,
    result_B=results_cond_b,
    score_col='combined_score',
    significance_col='if_significant',
    gene_center_col='gene_center',
    gene_motif_col='gene_motif',
    percentile_threshold=95.0,
    background='Significant',  # or 'Overlapping' for more pairs
)
# diff columns: gene_center, gene_motif, score_A, score_B, score_diff,
#   percentile, is_outlier, significant_in_A, significant_in_B, outlier_direction
# outlier_direction values: 'higher_in_A', 'lower_in_A', 'not_outlier'

# Group-enriched pairs: significant in that condition AND score is an outlier in that direction
enriched_in_A = diff[
    (diff['outlier_direction'] == 'higher_in_A') &
    (diff['significant_in_A'] == True)
]
enriched_in_B = diff[
    (diff['outlier_direction'] == 'lower_in_A') &   # lower_in_A = higher_in_B
    (diff['significant_in_B'] == True)
]

print(f"Pairs enriched in condition A: {len(enriched_in_A)}")
print(enriched_in_A[['gene_center', 'gene_motif', 'score_A', 'score_B', 'score_diff']].head(10).to_string())

print(f"\nPairs enriched in condition B: {len(enriched_in_B)}")
print(enriched_in_B[['gene_center', 'gene_motif', 'score_A', 'score_B', 'score_diff']].head(10).to_string())
```

**For per-type results**, subset by `cell_type` before calling `test_score_difference`:

```python
results_a_by_type = sq_multi.compute_gene_gene_correlation_by_type(
    ct=ct, motif=['Macrophage', 'B cell'], dataset='condition_A', k=30)
results_b_by_type = sq_multi.compute_gene_gene_correlation_by_type(
    ct=ct, motif=['Macrophage', 'B cell'], dataset='condition_B', k=30)

# Compare per neighbor type
for neighbor_type in results_a_by_type['cell_type'].unique():
    sub_a = results_a_by_type[results_a_by_type['cell_type'] == neighbor_type]
    sub_b = results_b_by_type[results_b_by_type['cell_type'] == neighbor_type]
    if len(sub_a) == 0 or len(sub_b) == 0:
        continue

    diff = spatial_query_multi.test_score_difference(
        result_A=sub_a,
        result_B=sub_b,
        percentile_threshold=95.0,
        background='Significant',
    )
    enriched_A = diff[(diff['outlier_direction'] == 'higher_in_A') & (diff['significant_in_A'])]
    enriched_B = diff[(diff['outlier_direction'] == 'lower_in_A') & (diff['significant_in_B'])]
    print(f"\n[{neighbor_type}] A-enriched: {len(enriched_A)}, B-enriched: {len(enriched_B)}")
```

---

### Choosing D1 vs D2

| Condition | Recommended |
|-----------|-------------|
| Motif has exactly 1 non-anchor cell type | D1 |
| Motif has 2+ non-anchor cell types | D2 (default) |
| User asks for "overall" or "pooled" | D1 |
| User asks "per cell type" or "type-specific" | D2 |

---

## Common failure modes

| Problem | Cause | Fix |
|---------|-------|-----|
| Dataset name not found | Underscore replaced with hyphen internally | Use `sq_multi.datasets` to see actual names |
| Empty `ind_group1` or `ind_group2` dict | No matching cells in any dataset | Check that cell type exists in all datasets; lower `min_support` |
| `differential_analysis` returns empty | Few samples per condition | Ensure ≥2 samples per condition |
| `plot_differential_pattern_heatmap` import error | Not exposed at top level | Import from `SpatialQuery.plotting` |
