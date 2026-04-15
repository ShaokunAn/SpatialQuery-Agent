# Autopilot: End-to-End Multi-FOV Analysis

## When to use this tutorial
Use this tutorial when the user provides multiple AnnData files (multiple samples, conditions, or tissue sections) and requests a complete or "full" spatial analysis. Execute steps in order; cross-condition steps (Steps 5, 6b, 7, 8) only apply if two distinct condition groups exist.

Typical user requests:
- "Run a full analysis on my healthy and disease samples"
- "Compare spatial neighborhoods across conditions"
- "Do a complete multi-sample spatial analysis on T cells"
- "Run everything on these samples end to end"

For other DE comparisons (within-dataset proximal vs. distal, per-sample comparisons), see `pipeline_multi_fov.md`.

---

## Overview of steps

```
Step 0: Load multiple AnnData files and inspect keys
Step 1: Initialize spatial_query_multi object
Step 2: Visualize cell type distribution across samples
Step 3: Motif enrichment — radius-based, aggregated across all FOVs (with retry)
Step 4: Visualize significant motifs
Step 5: DE — motif+ anchor in condition A vs. condition B  [cross-condition only]
Step 6: Gene co-variation — per-type (compute_gene_gene_correlation_by_type)
Step 7: Cross-condition co-variation comparison with test_score_difference  [cross-condition only]
Step 8: Differential motif analysis between conditions  [cross-condition only]
Step 9: Print summary
```

---

## Step 0: Load multiple AnnData files and inspect keys

```python
import anndata
import numpy as np
import pandas as pd
from SpatialQuery import spatial_query_multi

# --- Load AnnData files ---
# Replace with actual file paths provided by the user
adata_list = [
    anndata.read_h5ad("path/to/sample1.h5ad"),
    anndata.read_h5ad("path/to/sample2.h5ad"),
    anndata.read_h5ad("path/to/sample3.h5ad"),
    # add more as needed
]

# Dataset names — use meaningful names with condition as a prefix
# IMPORTANT: underscores in names are replaced with hyphens internally
# Use consistent prefixes so conditions can be identified later
# e.g., 'healthy1', 'healthy2', 'disease1', 'disease2'
dataset_names = ['healthy1', 'healthy2', 'disease1']   # must match len(adata_list)

# --- Inspect keys using the first adata as reference ---
# Keys should be consistent across all samples
adata_ref = adata_list[0]

print("obsm keys (spatial coordinates):")
print(list(adata_ref.obsm.keys()))

print("\nobs columns (cell type labels):")
print(adata_ref.obs.columns.tolist())

print("\nvar columns (gene names):")
print(adata_ref.var.columns.tolist())
print("var index (first 5):", adata_ref.var.index[:5].tolist())

# Determine label_key and list all cell types
label_key = 'cell_type'   # replace with actual key from obs columns
print(f"\nAvailable cell types in '{label_key}':")
print(adata_ref.obs[label_key].unique().tolist())

# Verify all adatas have the same label_key
for i, ad in enumerate(adata_list):
    if label_key not in ad.obs.columns:
        raise KeyError(
            f"adata[{i}] ('{dataset_names[i]}') is missing label_key '{label_key}'"
        )

# Print per-sample cell counts
for name, ad in zip(dataset_names, adata_list):
    print(f"\n{name}: {ad.shape[0]} cells, {ad.shape[1]} genes")
    print(ad.obs[label_key].value_counts().to_string())
```

**Parameter inference rules**:
- `spatial_key`: key from `adata.obsm.keys()` (prefer `'X_spatial'`, then `'spatial'`)
- `label_key`: column from `adata.obs.columns` (prefer `'cell_type'`, then `'predicted_label'`)
- `feature_name`: column from `adata.var.columns` (prefer `'gene'`, then `'gene_name'`). Use `None` if gene names are in `adata.var.index`.
- `ct`: must be EXACTLY one string from `adata.obs[label_key].unique()`. Never invent.

---

## Step 1: Initialize spatial_query_multi object

```python
# Replace with actual values from Step 0
spatial_key = 'X_spatial'
label_key = 'cell_type'
feature_name = 'gene'    # or None if gene names are in adata.var.index
ct = 'T cell'            # anchor cell type — must match exactly

sq_multi = spatial_query_multi(
    adatas=adata_list,
    datasets=dataset_names,   # underscores → hyphens internally
    spatial_key=spatial_key,
    label_key=label_key,
    feature_name=feature_name,
    build_gene_index=False,   # False = use raw expression (recommended)
    if_lognorm=True,
)

# --- IMPORTANT: always use sq_multi.datasets for actual internal names ---
# After init, names may be modified:
#   underscores → hyphens: 'healthy_1' → 'healthy-1'
#   duplicate names get suffix: 'healthy' x2 → 'healthy_0', 'healthy_1'
# These internal names are required as dictionary keys for de_genes()
print("Internal dataset names (use these as dict keys):")
print(sq_multi.datasets)
all_datasets = sq_multi.datasets  # store for reuse

# Verify anchor cell type exists
anchor_found = any(ct in ad.obs[label_key].values for ad in adata_list)
if not anchor_found:
    raise ValueError(
        f"Anchor cell type '{ct}' not found in any dataset. "
        f"Available: {adata_ref.obs[label_key].unique().tolist()}"
    )
print(f"Anchor '{ct}' confirmed.")
```

---

## Step 2: Visualize cell type distribution across samples

```python
# Stacked bar chart: cell type proportion across all samples
sq_multi.plot_cell_type_distribution(
    dataset=None,             # None = all datasets
    data_type='proportion',   # 'proportion' (0–1) or 'number' (raw counts)
    colormap='tab20c',
)
```

---

## Step 3: Motif enrichment — radius-based, aggregated across all FOVs

Default method is **radius-based** (`motif_enrichment_dist` with `max_dist=20`).

When `return_cellID=True`, this returns a **3-tuple**:
- `result`: DataFrame with enrichment statistics
- `motif_cell_ids`: `{motif_str: {dataset_name: [cell_indices]}}` — neighbor cells per FOV
- `center_cell_ids`: `{motif_str: {dataset_name: [cell_indices]}}` — anchor cells with motif per FOV

The `dataset_name` keys inside these dicts are the **internal names** from `sq_multi.datasets`.

```python
max_dist = 20
min_support = 0.5

print(f"Running motif enrichment around '{ct}' across all FOVs "
      f"(max_dist={max_dist}, min_support={min_support})...")

result, motif_cell_ids, center_cell_ids = sq_multi.motif_enrichment_dist(
    ct=ct,
    motifs=None,          # None = auto-discover
    dataset=None,         # None = all datasets
    max_dist=max_dist,
    min_size=0,
    min_support=min_support,
    return_cellID=True,   # required — returns 3-tuple
)
# result columns: center, motifs, n_center_motif, n_center, n_motif,
#                 expectation, p-values, adj-pval, if_significant

sig = result[result['if_significant']]
print(f"  Found {len(sig)} significant motifs at min_support={min_support}")

# --- Retry with lower min_support if needed ---
if len(sig) == 0:
    min_support = 0.3
    print(f"  Retrying with min_support={min_support}...")
    result, motif_cell_ids, center_cell_ids = sq_multi.motif_enrichment_dist(
        ct=ct, motifs=None, dataset=None, max_dist=max_dist,
        min_size=0, min_support=min_support, return_cellID=True,
    )
    sig = result[result['if_significant']]
    print(f"  Found {len(sig)} significant motifs")

if len(sig) == 0:
    min_support = 0.1
    print(f"  Retrying with min_support={min_support}...")
    result, motif_cell_ids, center_cell_ids = sq_multi.motif_enrichment_dist(
        ct=ct, motifs=None, dataset=None, max_dist=max_dist,
        min_size=0, min_support=min_support, return_cellID=True,
    )
    sig = result[result['if_significant']]
    print(f"  Found {len(sig)} significant motifs")

if len(sig) == 0:
    raise RuntimeError(
        f"No significant motifs found around '{ct}' at any min_support level.\n"
        f"Suggestions:\n"
        f"  1. Verify '{ct}' exists across datasets\n"
        f"  2. Increase max_dist (e.g., max_dist=30)\n"
        f"  3. Check cell type counts per sample"
    )

# --- Inspect key structure of center_cell_ids ---
# center_cell_ids keys are string representations of the motif frozensets
# Use these keys to look up cell indices per dataset
print("\nMotif keys in center_cell_ids:")
for k in list(center_cell_ids.keys())[:3]:
    datasets_with_motif = list(center_cell_ids[k].keys())
    print(f"  Key: '{k}'  →  datasets: {datasets_with_motif}")

# --- Display significant motifs ---
print(f"\n=== Significant motifs around '{ct}' ===")
for i, (_, row) in enumerate(sig.iterrows()):
    motif_list = sorted(list(row['motifs']))
    print(f"  Motif {i+1}: {motif_list}  |  n_center={row['n_center_motif']}  |  adj-pval={row['adj-pval']:.4f}")
```

---

## Step 4: Visualize significant motifs

```python
# Enrichment heatmap: all significant motifs in one view
sq_multi.plot_motif_enrichment_heatmap(
    sig,
    figsize=(7, 5),
    title=f'Motif enrichment around {ct} (aggregated across FOVs)',
)
```

---

## Step 5: DE — motif+ anchor in condition A vs. condition B

**Skip this step if all samples are from the same condition.**

Compare anchor cells that have the motif in condition A vs. condition B. This reveals what changes in the anchor cell across conditions when the motif is present.

`sq_multi.de_genes` takes cell groups as `Dict[str, List[int]]` keyed by **internal dataset names** (from `sq_multi.datasets`).

```python
# --- Define condition groups using internal dataset names ---
# Adjust the filter strings to match your actual dataset name prefixes
cond_a_datasets = [ds for ds in all_datasets if 'healthy' in ds]
cond_b_datasets = [ds for ds in all_datasets if 'disease' in ds]

print(f"Condition A datasets: {cond_a_datasets}")
print(f"Condition B datasets: {cond_b_datasets}")

if len(cond_a_datasets) == 0 or len(cond_b_datasets) == 0:
    print("Skipping Step 5: could not identify two condition groups from dataset names.")
    print("To run manually, set cond_a_datasets and cond_b_datasets as lists of dataset names.")
else:
    de_results = {}

    for i, (_, row) in enumerate(sig.iterrows()):
        motif_list = sorted(list(row['motifs']))
        # center_cell_ids key: find the matching key for this motif
        motif_key = str(row['motifs'])   # default: string of frozenset
        if motif_key not in center_cell_ids:
            # frozenset string representations can vary — search by set equality
            for k in center_cell_ids.keys():
                try:
                    if set(eval(k)) == set(motif_list):
                        motif_key = k
                        break
                except Exception:
                    continue
            else:
                print(f"  Motif {i+1} key not found in center_cell_ids. Skipping.")
                continue

        print(f"\n--- DE Motif {i+1}: {motif_list} ---")

        # Group 1: motif+ anchor cells from condition A datasets
        ind_group1 = {}
        for ds in cond_a_datasets:
            cells = center_cell_ids[motif_key].get(ds, [])
            if len(cells) >= 3:
                ind_group1[ds] = list(cells)

        # Group 2: motif+ anchor cells from condition B datasets
        ind_group2 = {}
        for ds in cond_b_datasets:
            cells = center_cell_ids[motif_key].get(ds, [])
            if len(cells) >= 3:
                ind_group2[ds] = list(cells)

        total_a = sum(len(v) for v in ind_group1.values())
        total_b = sum(len(v) for v in ind_group2.values())
        print(f"  Cond A motif+ anchors: {total_a} cells across {len(ind_group1)} datasets")
        print(f"  Cond B motif+ anchors: {total_b} cells across {len(ind_group2)} datasets")

        if len(ind_group1) == 0 or len(ind_group2) == 0:
            print(f"  Skipped: motif not found in one of the conditions")
            continue

        de = sq_multi.de_genes(
            ind_group1=ind_group1,
            ind_group2=ind_group2,
            method='wilcoxon',    # options: 'wilcoxon', 't-test', 'fisher'
            min_fraction=0.05,
        )
        # de columns: gene, proportion_1, proportion_2, abs, difference, p_value, adj-pval, de_in
        # de_in: 'group1' = higher in condition A; 'group2' = higher in condition B

        de_sig = de[de['adj-pval'] < 0.05].sort_values('adj-pval')
        de_results[str(motif_list)] = de_sig

        print(f"  Significant DE genes: {len(de_sig)}")
        if len(de_sig) > 0:
            print(de_sig[['gene', 'proportion_1', 'proportion_2', 'adj-pval', 'de_in']].head(10).to_string())
```

---

## Step 6: Gene co-variation — per-type (compute_gene_gene_correlation_by_type)

Run per-type co-variation for each significant motif, aggregated across all FOVs.

**Key differences from single-FOV**:
- `compute_gene_gene_correlation_by_type` returns only a **DataFrame** (no second element)
- Significance column: use `if_significant` (= `reject_test1_fdr & reject_test2_fdr`), same as single-FOV
- FDR p-values are in `q_value_test1` and `q_value_test2`

```python
corr_results = {}

for i, (_, row) in enumerate(sig.iterrows()):
    motif_list = sorted(list(row['motifs']))
    non_anchor_types = [m for m in motif_list if m != ct]

    print(f"\n--- Covariation Motif {i+1}: {motif_list} ---")
    print(f"  Non-anchor types: {non_anchor_types}")

    if len(non_anchor_types) == 0:
        print(f"  Skipped: no non-anchor cell types in motif")
        continue

    try:
        results_df = sq_multi.compute_gene_gene_correlation_by_type(
            ct=ct,
            motif=motif_list,
            dataset=None,      # None = aggregate across all FOVs
            max_dist=max_dist, # use max_dist= to match radius-based enrichment
                               # do NOT pass k= at the same time
            min_size=0,
            min_nonzero=10,
            alpha=None,
        )
        # results_df columns: cell_type, gene_center, gene_motif,
        #   corr_neighbor, corr_non_neighbor, corr_center_no_motif,
        #   p_value_test1, p_value_test2, delta_corr_test1, delta_corr_test2,
        #   combined_score, abs_combined_score,
        #   q_value_test1, q_value_test2,
        #   reject_test1_fdr, reject_test2_fdr,
        #   if_significant (= reject_test1_fdr & reject_test2_fdr)

        sig_pairs = results_df[results_df['if_significant']]
        corr_results[str(motif_list)] = sig_pairs

        print(f"  Significant co-varying gene pairs (all neighbor types): {len(sig_pairs)}")

        # Breakdown per neighbor type
        for neighbor_type in results_df['cell_type'].unique():
            sub = sig_pairs[sig_pairs['cell_type'] == neighbor_type]
            print(f"    [{neighbor_type}]: {len(sub)} significant pairs")
            if len(sub) > 0:
                print(sub[['gene_center', 'gene_motif', 'combined_score']].head(5).to_string())

        if len(sig_pairs) > 0:
            sp.plot_gene_pair_heatmap(sig_pairs, figsize=(7, 5))
        else:
            print(f"  No significant pairs. Try: lower min_nonzero=5 or alpha=0.1.")

    except Exception as e:
        print(f"  Failed: {e}")
```

---

## Step 7: Cross-condition co-variation comparison with test_score_difference

**Skip this step if only one condition.**

For each significant motif, compute co-variation separately for each condition, then compare to find gene pairs enriched in one condition vs. the other.

`test_score_difference` is a **static method** — call it on the class directly.

Return value columns:
- `gene_center`, `gene_motif`: the gene pair
- `score_A`, `score_B`: `combined_score` in each condition
- `score_diff`: `score_A − score_B`
- `percentile`: percentile rank of |score_diff|
- `is_outlier`: True if in top/bottom `percentile_threshold`%
- `significant_in_A`, `significant_in_B`: whether significant in that condition
- `outlier_direction`:
  - `'higher_in_A'` = co-variation score is higher in condition A
  - `'lower_in_A'` = co-variation score is lower in condition A (= higher in B)
  - `'not_outlier'` = no significant difference

To get **condition A-enriched pairs**: `outlier_direction == 'higher_in_A'` AND `significant_in_A == True`
To get **condition B-enriched pairs**: `outlier_direction == 'lower_in_A'` AND `significant_in_B == True`

```python
if len(cond_a_datasets) > 0 and len(cond_b_datasets) > 0:

    for i, (_, row) in enumerate(sig.iterrows()):
        motif_list = sorted(list(row['motifs']))
        non_anchor_types = [m for m in motif_list if m != ct]

        if len(non_anchor_types) == 0:
            continue

        print(f"\n--- Cross-condition covariation Motif {i+1}: {motif_list} ---")

        try:
            # Compute co-variation per condition separately
            res_a = sq_multi.compute_gene_gene_correlation_by_type(
                ct=ct,
                motif=motif_list,
                dataset=cond_a_datasets,   # restrict to condition A FOVs
                max_dist=max_dist,
                min_nonzero=10,
            )
            res_b = sq_multi.compute_gene_gene_correlation_by_type(
                ct=ct,
                motif=motif_list,
                dataset=cond_b_datasets,   # restrict to condition B FOVs
                max_dist=max_dist,
                min_nonzero=10,
            )

            # Compare per neighbor cell type
            for neighbor_type in res_a['cell_type'].unique():
                sub_a = res_a[res_a['cell_type'] == neighbor_type].copy()
                sub_b = res_b[res_b['cell_type'] == neighbor_type].copy()

                if len(sub_a) == 0 or len(sub_b) == 0:
                    print(f"  [{neighbor_type}] Skipped: no results in one condition")
                    continue

                diff = spatial_query_multi.test_score_difference(
                    result_A=sub_a,
                    result_B=sub_b,
                    score_col='combined_score',
                    significance_col='if_significant',
                    gene_center_col='gene_center',
                    gene_motif_col='gene_motif',
                    percentile_threshold=95.0,
                    background='Significant',  # use 'Overlapping' for more pairs
                )

                # Condition A-enriched: high score in A, significant in A
                enriched_A = diff[
                    (diff['outlier_direction'] == 'higher_in_A') &
                    (diff['significant_in_A'] == True)
                ]
                # Condition B-enriched: high score in B (= low in A), significant in B
                enriched_B = diff[
                    (diff['outlier_direction'] == 'lower_in_A') &
                    (diff['significant_in_B'] == True)
                ]

                print(f"\n  [{neighbor_type}]")
                print(f"    Cond A-enriched pairs: {len(enriched_A)}")
                if len(enriched_A) > 0:
                    print(enriched_A[['gene_center', 'gene_motif', 'score_A', 'score_B', 'score_diff']].head(5).to_string())
                print(f"    Cond B-enriched pairs: {len(enriched_B)}")
                if len(enriched_B) > 0:
                    print(enriched_B[['gene_center', 'gene_motif', 'score_A', 'score_B', 'score_diff']].head(5).to_string())

        except Exception as e:
            print(f"  Failed: {e}")
```

---

## Step 8: Differential motif analysis between conditions

**Skip this step if only one condition.**

Identifies which motifs are significantly more enriched in one condition vs. the other.

```python
from SpatialQuery.plotting import plot_differential_pattern_heatmap

if len(cond_a_datasets) > 0 and len(cond_b_datasets) > 0:
    print(f"\n--- Differential motif analysis: condition A vs. condition B ---")

    try:
        diff_motif = sq_multi.differential_analysis_dist(
            ct=ct,
            datasets=[cond_a_datasets, cond_b_datasets],  # list of two lists
            motifs=None,          # None = auto-discover; or pass a list of motifs
            max_dist=max_dist,
            min_support=min_support,
            min_size=0,
        )
        # diff_motif: Dict[str, DataFrame] keyed by dataset name
        # Each DataFrame columns: itemsets, support_{ds0}_mean, support_{ds1}_mean, adj-pval

        # Visualize
        plot_differential_pattern_heatmap(
            diff_motif, ct=ct, figsize=None, cmap='YlOrRd'
        )

        # Print differentially enriched motifs
        for ds_key, df in diff_motif.items():
            diff_sig = df[df['adj-pval'] < 0.05] if 'adj-pval' in df.columns else df
            print(f"  [{ds_key}]: {len(diff_sig)} differentially enriched motifs")
            if len(diff_sig) > 0:
                print(diff_sig.head(10).to_string())

    except Exception as e:
        print(f"  Differential analysis failed: {e}")
```

---

## Step 9: Summary

```python
print("\n" + "="*60)
print("AUTOPILOT ANALYSIS SUMMARY (Multi-FOV)")
print("="*60)
print(f"Anchor cell type     : {ct}")
print(f"Total FOVs           : {len(all_datasets)} — {all_datasets}")
print(f"Total cells          : {sum(ad.shape[0] for ad in adata_list)}")
print(f"Neighborhood radius  : max_dist={max_dist}")
print(f"Significant motifs   : {len(sig)}")
for i, (_, row) in enumerate(sig.iterrows()):
    print(f"  Motif {i+1}: {sorted(list(row['motifs']))}  adj-pval={row['adj-pval']:.4f}")

if len(cond_a_datasets) > 0 and len(cond_b_datasets) > 0:
    print(f"\nConditions:")
    print(f"  Condition A: {cond_a_datasets}")
    print(f"  Condition B: {cond_b_datasets}")

    print(f"\nDE results (motif+ anchor: condA vs condB):")
    if de_results:
        for motif_key, df in de_results.items():
            print(f"  {motif_key}: {len(df)} significant genes")
    else:
        print("  None")

print(f"\nCo-variation results (per-type):")
if corr_results:
    for motif_key, df in corr_results.items():
        print(f"  {motif_key}: {len(df)} significant gene pairs")
else:
    print("  None")
print("="*60)
```

---

## Common failure modes and fixes

| Problem | Likely cause | Fix |
|---------|-------------|-----|
| `motif_enrichment_dist(return_cellID=True)` returns only 1 value | Missing `return_cellID=True` | Ensure the call includes `return_cellID=True`; it returns a 3-tuple |
| `KeyError` in `de_genes` | Dict key not in `sq_multi.datasets` | Print `sq_multi.datasets` and use exact names |
| `center_cell_ids` key not found | frozenset string format varies | Loop over keys and compare using `set(eval(k)) == set(motif_list)` |
| Conditions not identified automatically | Dataset names don't contain condition prefix | Manually set `cond_a_datasets` and `cond_b_datasets` |
| `compute_gene_gene_correlation_by_type` fails | Motif has 0 non-anchor types | Check `non_anchor_types` is non-empty before calling |
| `test_score_difference` finds no outliers | Too few shared gene pairs | Switch `background='Overlapping'`; lower `percentile_threshold=90.0` |
| `if_significant` column missing | Very old version | Use `reject_test1_fdr & reject_test2_fdr` as fallback |
| `plot_differential_pattern_heatmap` ImportError | Not at top-level export | Import from `SpatialQuery.plotting` |
| Co-variation error "use k OR max_dist" | Both passed simultaneously | Pass only `max_dist=max_dist`; never pass `k=` as well |
