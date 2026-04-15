# Autopilot: End-to-End Single-FOV Analysis

## When to use this tutorial
Use this tutorial when the user provides a single AnnData file and requests a complete or "full" spatial analysis around an anchor cell type. Each step is self-contained with complete code. Execute steps in order.

Typical user requests:
- "Run a full spatial analysis on T cells"
- "Do a complete analysis around macrophages"
- "Run everything on this data for fibroblasts"
- "Analyze spatial neighborhoods of [cell type] end to end"

For more flexible DE comparisons (proximal vs. distal neighbors, cross-motif anchor comparison), see `pipeline_b_motif_de.md`.

---

## Overview of steps

```
Step 0: Load data and inspect keys
Step 1: Initialize spatial_query object
Step 2: Motif enrichment — radius-based (with retry)
Step 3: Visualize significant motifs
Step 4: Differential expression — motif+ anchor vs. motif− anchor
Step 5: Gene co-variation — per-type (compute_gene_gene_correlation_by_type)
Step 6: Print summary
```

---

## Step 0: Load data and inspect keys

Always inspect the AnnData keys before initializing. Do NOT assume key names — always read them from the data.

```python
import anndata
import numpy as np
import pandas as pd
from SpatialQuery import spatial_query

# Load the AnnData file (user provides the path)
adata = anndata.read_h5ad("path/to/data.h5ad")

# --- Inspect spatial coordinate key ---
# Look for keys like 'X_spatial', 'spatial', 'X_umap'
print("obsm keys (look for spatial coordinates):")
print(list(adata.obsm.keys()))

# --- Inspect cell label key ---
# Look for columns like 'cell_type', 'predicted_label', 'celltype', 'leiden'
print("\nobs columns (look for cell type label):")
print(adata.obs.columns.tolist())

# --- Inspect gene name key ---
# Look for columns like 'gene', 'gene_name', 'feature_name'
# If adata.var.columns is empty, gene names are in adata.var.index
print("\nvar columns (look for gene name column):")
print(adata.var.columns.tolist())
print("var index (gene names may be stored here if var.columns is empty):")
print(adata.var.index[:5].tolist())

# --- List available cell types ---
# IMPORTANT: the anchor cell type (ct) must be EXACTLY one of these strings
label_key = 'cell_type'  # replace with actual key from adata.obs.columns
print(f"\nAvailable cell types in '{label_key}':")
cell_types = adata.obs[label_key].unique().tolist()
print(cell_types)

# --- Basic data summary ---
print(f"\nData shape: {adata.shape[0]} cells x {adata.shape[1]} genes")
```

**Parameter inference rules**:
- `spatial_key`: use the key from `adata.obsm.keys()` that contains spatial coordinates (prefer `'X_spatial'`, then `'spatial'`)
- `label_key`: use the column from `adata.obs.columns` that contains cell type labels (prefer `'cell_type'`, then `'predicted_label'`, then `'celltype'`)
- `feature_name`: use the column from `adata.var.columns` that contains gene names (prefer `'gene'`, then `'gene_name'`, then `'feature_name'`). If `adata.var.columns` is empty, gene names are in the index — pass `feature_name=None`.
- `ct` (anchor cell type): must be EXACTLY one string from `adata.obs[label_key].unique()`. Never invent or abbreviate.

---

## Step 1: Initialize spatial_query object

Replace the parameter values with what you found in Step 0.

```python
# Replace these with actual values found in Step 0
spatial_key = 'X_spatial'      # key in adata.obsm
label_key = 'cell_type'        # key in adata.obs
feature_name = 'gene'          # key in adata.var (or None if genes are in adata.var.index)
ct = 'T cell'                  # anchor cell type — must match exactly

# Initialize
sp = spatial_query(
    adata=adata,
    dataset='sample1',               # arbitrary name for this sample
    spatial_key=spatial_key,
    label_key=label_key,
    feature_name=feature_name,
    build_gene_index=False,          # False = use raw expression matrix (recommended)
    if_lognorm=True,                 # True = data needs log-normalization; False = already done
    if_normalize_spatial_coord=True, # True = auto-normalize so 1 unit ≈ 1 cell diameter
)

# Verify initialization
print(f"Initialized spatial_query with {len(sp.labels)} cells")
print(f"Unique cell types: {list(set(sp.labels))}")
print(f"Number of genes: {len(sp.genes)}")

# Confirm anchor cell type exists
if ct not in set(sp.labels):
    raise ValueError(
        f"Anchor cell type '{ct}' not found. "
        f"Available types: {list(set(sp.labels))}"
    )
print(f"Anchor cell type '{ct}': {sum(l == ct for l in sp.labels)} cells")
```

---

## Step 2: Motif enrichment (radius-based, with retry)

Default method is **radius-based** (`motif_enrichment_dist` with `max_dist=20`). This is preferred when spatial coordinates have physical meaning or after normalization (1 unit ≈ 1 cell diameter).

`return_cellID=True` is required in autopilot — Steps 4 and 5 need the cell index columns (`center_id`, `neighbor_id`).

```python
# --- Configuration ---
max_dist = 20      # neighborhood radius in normalized units (≈ cell diameters)
min_support = 0.5  # start here; will retry at 0.3 and 0.1 if needed

# --- Run motif enrichment ---
# IMPORTANT: return_cellID=True adds 'center_id' and 'neighbor_id' columns to result
print(f"Running motif enrichment around '{ct}' (max_dist={max_dist}, min_support={min_support})...")

result = sp.motif_enrichment_dist(
    ct=ct,
    motifs=None,          # None = auto-discover all enriched motifs
    max_dist=max_dist,
    min_size=0,
    min_support=min_support,
    return_cellID=True,   # adds 'center_id' and 'neighbor_id' to each row
)
# result columns: center, motifs, n_center_motif, n_center, n_motif,
#                 expectation, p-values, adj-pval, if_significant,
#                 center_id (list of anchor cell indices WITH this motif),
#                 neighbor_id (list of neighbor cell indices IN this motif)

sig = result[result['if_significant']]
print(f"  Found {len(sig)} significant motifs at min_support={min_support}")

# --- Retry with lower min_support if no significant motifs ---
if len(sig) == 0:
    min_support = 0.3
    print(f"  Retrying with min_support={min_support}...")
    result = sp.motif_enrichment_dist(
        ct=ct, motifs=None, max_dist=max_dist, min_size=0,
        min_support=min_support, return_cellID=True,
    )
    sig = result[result['if_significant']]
    print(f"  Found {len(sig)} significant motifs at min_support={min_support}")

if len(sig) == 0:
    min_support = 0.1
    print(f"  Retrying with min_support={min_support}...")
    result = sp.motif_enrichment_dist(
        ct=ct, motifs=None, max_dist=max_dist, min_size=0,
        min_support=min_support, return_cellID=True,
    )
    sig = result[result['if_significant']]
    print(f"  Found {len(sig)} significant motifs at min_support={min_support}")

if len(sig) == 0:
    raise RuntimeError(
        f"No significant motifs found around '{ct}' at any min_support level.\n"
        f"Suggestions:\n"
        f"  1. Check that '{ct}' has enough cells: "
        f"{sum(l == ct for l in sp.labels)} cells found\n"
        f"  2. Try increasing max_dist (e.g., max_dist=30 or max_dist=50)\n"
        f"  3. Check that if_normalize_spatial_coord=True was set during init"
    )

# --- Display significant motifs ---
print(f"\n=== Significant motifs around '{ct}' ===")
for i, (_, row) in enumerate(sig.iterrows()):
    motif_list = sorted(list(row['motifs']))  # motifs column is a frozenset
    n_center = row['n_center_motif']
    padj = row['adj-pval']
    print(f"  Motif {i+1}: {motif_list}  |  n_center={n_center}  |  adj-pval={padj:.4f}")
```

---

## Step 3: Visualize significant motifs

```python
# Spatial plot: show where each motif occurs in tissue
# IMPORTANT: row['motifs'] is a frozenset — must convert to list before passing to plot functions
for i, (_, row) in enumerate(sig.iterrows()):
    motif_list = sorted(list(row['motifs']))
    print(f"Plotting spatial distribution of motif {i+1}: {motif_list}")
    sp.plot_motif_celltype(
        ct=ct,
        motif=motif_list,
        max_dist=max_dist,
        figsize=(5, 5),
    )

# Enrichment heatmap: all significant motifs summarized in one view
sp.plot_motif_enrichment_heatmap(
    sig,
    figsize=(7, 5),
    title=f'Motif enrichment around {ct}',
)
```

---

## Step 4: Differential expression — motif+ anchor vs. motif− anchor

For each significant motif, compare anchor cells that HAVE this motif in their neighborhood vs. anchor cells that do NOT.

**Biological question**: What genes are differentially expressed in the anchor cell type depending on whether it has this specific spatial context?

- `center_id`: list of anchor cell indices that have this motif nearby (motif+ group)
- `ind_group2`: all anchor cells MINUS the motif+ cells (motif− group)
- IMPORTANT: use `np.setdiff1d` to compute the difference — do NOT use set subtraction on arrays

```python
# All anchor cell indices in the whole dataset
# np.where returns a tuple — use [0] to get the index array
all_anchor_ids = np.where(np.array(adata.obs[label_key]) == ct)[0]
print(f"Total '{ct}' cells in dataset: {len(all_anchor_ids)}")

de_results = {}   # key → significant DE DataFrame

for i, (_, row) in enumerate(sig.iterrows()):
    motif_list = sorted(list(row['motifs']))
    motif_str = str(motif_list)

    # center_id is stored as a list in the DataFrame cell
    center_ids = np.array(row['center_id'])   # anchor cells WITH this motif nearby

    # motif− group: all anchor cells MINUS those with the motif
    # Use np.setdiff1d — NOT subtraction (array shapes differ)
    non_center_ids = np.setdiff1d(all_anchor_ids, center_ids)

    print(f"\n--- DE Motif {i+1}: {motif_list} ---")
    print(f"  motif+ anchors (group 1): {len(center_ids)} cells")
    print(f"  motif- anchors (group 2): {len(non_center_ids)} cells")

    if len(center_ids) < 5 or len(non_center_ids) < 5:
        print(f"  Skipped: too few cells in one group (need >= 5 in each)")
        continue

    de = sp.de_genes(
        ind_group1=center_ids,
        ind_group2=non_center_ids,
        method='wilcoxon',     # options: 'wilcoxon', 't-test', 'fisher'
        min_fraction=0.05,     # only test genes expressed in >= 5% of cells
    )
    # de columns: gene, proportion_1, proportion_2, p_value, adj-pval, de_in
    # de_in: 'group1' means higher in motif+ anchors; 'group2' means higher in motif- anchors

    de_sig = de[de['adj-pval'] < 0.05].sort_values('adj-pval')
    de_results[motif_str] = de_sig

    print(f"  Significant DE genes: {len(de_sig)}")
    if len(de_sig) > 0:
        print(de_sig[['gene', 'proportion_1', 'proportion_2', 'adj-pval', 'de_in']].head(10).to_string())
    else:
        print(f"  No significant DE genes. Try method='t-test' or lower min_fraction=0.01.")
```

---

## Step 5: Gene co-variation — per-type (compute_gene_gene_correlation_by_type)

For each significant motif, compute gene pairs that show spatially coordinated expression between anchor cells and their neighbors. We always use `compute_gene_gene_correlation_by_type` (per-type), which reports results separately for each neighbor cell type in the motif.

**Why per-type**: when the motif contains multiple cell types, pooled analysis averages out type-specific signals. Per-type reveals which neighbor cell type drives co-variation with the anchor.

`compute_gene_gene_correlation_by_type` returns only a **DataFrame** (no second element).

```python
corr_results = {}   # key → significant pairs DataFrame

for i, (_, row) in enumerate(sig.iterrows()):
    motif_list = sorted(list(row['motifs']))
    motif_str = str(motif_list)

    # Identify non-anchor cell types in the motif
    non_anchor_types = [m for m in motif_list if m != ct]

    print(f"\n--- Covariation Motif {i+1}: {motif_list} ---")
    print(f"  Non-anchor neighbor types: {non_anchor_types}")

    if len(non_anchor_types) == 0:
        print(f"  Skipped: motif contains no non-anchor cell types")
        continue

    try:
        results_df = sp.compute_gene_gene_correlation_by_type(
            ct=ct,
            motif=motif_list,
            max_dist=max_dist,   # use max_dist= to match radius-based enrichment
                                 # do NOT pass k= at the same time
            min_size=0,
            min_nonzero=10,      # require >= 10 non-zero cells per gene
            alpha=None,          # default significance threshold (0.05)
        )
        # results_df columns: cell_type, gene_center, gene_motif,
        #   corr_neighbor, corr_non_neighbor, corr_center_no_motif,
        #   p_value_test1, p_value_test2, delta_corr_test1, delta_corr_test2,
        #   combined_score, abs_combined_score,
        #   q_value_test1, q_value_test2,
        #   reject_test1_fdr, reject_test2_fdr,
        #   if_significant
        #
        # Each row = one gene pair (gene_center from anchor, gene_motif from neighbor)
        # cell_type = which neighbor cell type this pair involves

        sig_pairs = results_df[results_df['if_significant']]
        corr_results[motif_str] = sig_pairs

        print(f"  Significant co-varying gene pairs (all neighbor types): {len(sig_pairs)}")

        # Breakdown per neighbor cell type
        for neighbor_type in results_df['cell_type'].unique():
            sub = sig_pairs[sig_pairs['cell_type'] == neighbor_type]
            print(f"    [{neighbor_type}]: {len(sub)} significant pairs")
            if len(sub) > 0:
                print(sub[['gene_center', 'gene_motif', 'combined_score']].head(5).to_string())

        if len(sig_pairs) > 0:
            # Heatmap with biclustering of significant gene pairs
            cluster_df = sp.plot_gene_pair_heatmap(sig_pairs, figsize=(7, 5))
        else:
            print(f"  No significant pairs. Try: lower min_nonzero=5, or lower alpha=0.1.")

    except Exception as e:
        print(f"  Failed with error: {e}")
```

---

## Step 6: Summary

```python
print("\n" + "="*55)
print("AUTOPILOT ANALYSIS SUMMARY (Single-FOV)")
print("="*55)
print(f"Anchor cell type    : {ct}")
print(f"Dataset             : {adata.shape[0]} cells x {adata.shape[1]} genes")
print(f"Neighborhood radius : max_dist={max_dist}")
print(f"Significant motifs  : {len(sig)}")
for i, (_, row) in enumerate(sig.iterrows()):
    print(f"  Motif {i+1}: {sorted(list(row['motifs']))}  adj-pval={row['adj-pval']:.4f}")

print(f"\nDE results (motif+ vs motif- anchor):")
if de_results:
    for motif_key, df in de_results.items():
        print(f"  {motif_key}: {len(df)} significant genes")
else:
    print("  None (groups too small or no significant genes)")

print(f"\nCo-variation results (per-type):")
if corr_results:
    for motif_key, df in corr_results.items():
        print(f"  {motif_key}: {len(df)} significant gene pairs")
else:
    print("  None")
print("="*55)
```

---

## Common failure modes and fixes

| Problem | Likely cause | Fix |
|---------|-------------|-----|
| `ValueError: ct not found` | Cell type name typo or wrong `label_key` | Print `adata.obs[label_key].unique()` and use the exact string |
| Zero significant motifs after all retries | Sparse data, wrong `spatial_key`, or unnormalized coordinates | Check `adata.obsm.keys()`; ensure `if_normalize_spatial_coord=True`; try `max_dist=30` |
| DE skipped (too few cells) | Anchor type has few cells or motif is very sparse | Lower `min_support`; check anchor cell count |
| Zero DE genes | Weak signal or `min_fraction` too high | Lower `min_fraction=0.01`; try `method='t-test'` |
| Co-variation error "use k OR max_dist" | Both parameters passed simultaneously | Pass only one: `max_dist=20` OR `k=30`, never both |
| `if_significant` column missing in results_df | Version difference | Use `results_df['reject_test1_fdr']` as fallback |
| `adj-pval` column not found | Wrong column name | Use `adj-pval` with a hyphen (not underscore) |
