# Pipeline A: Motif Discovery (Single-FOV)

## When to use
Use this pipeline when the user wants to identify which cell types tend to co-occur in the spatial neighborhood of a given anchor cell type. This is the most common entry point for spatial analysis. Output is a table of statistically enriched neighborhood motifs and spatial visualizations.

Typical user requests:
- "What cell types are enriched around T cells?"
- "Find spatial neighborhoods of macrophages"
- "Which cell types co-occur near fibroblasts?"

---

## Method call sequence

### Step 0 (optional): Discover frequent patterns automatically
If the user does not specify a motif, run `find_fp_knn` or `find_fp_dist` first to obtain candidate motifs.

```python
import anndata
import numpy as np
import pandas as pd
from SpatialQuery import spatial_query

# Load data
adata = anndata.read_h5ad("path/to/data.h5ad")

# Inspect keys
print(adata.obsm.keys())          # find spatial_key
print(adata.obs.columns.tolist()) # find label_key
print(adata.var.columns.tolist()) # find feature_name
print(adata.obs['cell_type'].unique().tolist())  # available cell types

# Initialize
sp = spatial_query(
    adata=adata,
    dataset='sample1',
    spatial_key='X_spatial',     # replace with actual key
    label_key='cell_type',       # replace with actual key
    feature_name='gene',         # replace with actual key
    build_gene_index=False,
    if_lognorm=True,
    if_normalize_spatial_coord=True,
)
```

**KNN-based (default):**
```python
ct = 'T cell'  # anchor cell type — must be in adata.obs[label_key].unique()

# Option A: auto-discover motifs via frequent patterns
fp = sp.find_fp_knn(ct=ct, k=30, min_support=0.5, max_dist=20)
# fp columns: itemsets, support
# itemsets is a frozenset; convert to list for display
fp['itemsets_list'] = fp['itemsets'].apply(list)
print(fp)
```

**Distance-based (default recommendation):**
```python
fp = sp.find_fp_dist(ct=ct, max_dist=20, min_size=0, min_support=0.5)
```

---

### Step 1: Motif enrichment
Run enrichment to get significance statistics. Use `return_cellID=False` unless this is followed immediately by DE or covariation analysis.

**KNN-based:**
```python
result = sp.motif_enrichment_knn(
    ct=ct,
    motifs=None,       # None = auto-discover; or pass list e.g. ['B cell', 'NK cell']
    k=30,
    min_support=0.5,
    max_dist=20,
    return_cellID=False,  # False for pure motif discovery
)
# result columns: center, motifs, n_center_motif, n_center, n_motif,
#                 expectation, p-values, adj-pval, if_significant
```

**Distance-based (recommended default):**
```python
result = sp.motif_enrichment_dist(
    ct=ct,
    motifs=None,
    max_dist=20,
    min_size=0,
    min_support=0.5,
    return_cellID=False,
)
```

Filter to significant motifs:
```python
sig = result[result['if_significant']]
print(f"Found {len(sig)} significant motifs")
print(sig[['motifs', 'n_center_motif', 'adj-pval']].to_string())
```

---

### Step 2: Spatial visualization per motif
Visualize where each significant motif occurs in tissue space.

```python
for _, row in sig.iterrows():
    motif = list(row['motifs'])  # motifs is a frozenset
    sp.plot_motif_celltype(ct=ct, motif=motif, max_dist=20, figsize=(5, 5))
```

---

### Step 3: Enrichment heatmap
Summarize all significant motifs in a single heatmap.

```python
sp.plot_motif_enrichment_heatmap(sig, figsize=(7, 5), title=f'Motif enrichment around {ct}')
```

---

## KNN vs. distance: when to use each

| Method | Use when |
|--------|----------|
| `motif_enrichment_dist` | Coordinates are meaningful (e.g., physical µm); you want a fixed-radius neighborhood |
| `motif_enrichment_knn` | You want a consistent number of neighbors regardless of local cell density |

Default to distance-based. Switch to KNN if the user asks for "k nearest neighbors" or if density is highly variable.

---

## Data flow summary

```
adata
  └─ sp = spatial_query(adata, ...)
       ├─ [optional] find_fp_knn/dist(ct) → fp (candidate motifs)
       ├─ motif_enrichment_knn/dist(ct, motifs=None) → result
       │    └─ result[result['if_significant']] → sig
       ├─ plot_motif_celltype(ct, motif) [per motif in sig]
       └─ plot_motif_enrichment_heatmap(sig)
```

---

## Common failure modes

| Problem | Cause | Fix |
|---------|-------|-----|
| Zero significant motifs | `min_support` too high or `max_dist` too small | Lower `min_support` to 0.3 then 0.1; increase `max_dist` to 30 or 50 |
| Cell type not found error | `ct` not in `adata.obs[label_key].unique()` | Always verify ct against `adata.obs[label_key].unique()` before calling |
| Empty `fp` from find_fp | `min_support` too stringent for sparse tissue | Lower `min_support`; check cell count for this ct |
| `motifs` frozenset display issue | frozenset is not serializable | Use `.apply(list)` or `str()` for display |
