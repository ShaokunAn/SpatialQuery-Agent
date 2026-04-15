# motif_enrichment_dist (spatial_query)

## Description
Perform motif enrichment analysis within a radius-based neighborhood of a center cell type in a single FOV. This method can unbiasedly discover frequent patterns around anchor cells and quantify their enrichment significance using the hypergeometric test, with FDR correction for multiple testing. Users can also specify custom motifs to test their significance in the neighborhood of a center cell type.

## Function Signature
```python
sp.motif_enrichment_dist(
    ct: str,
    motifs: Union[str, List[str]] = None,
    max_dist: float = 20,
    min_size: int = 0,
    min_support: float = 0.5,
    return_cellID: bool = False,
) -> pd.DataFrame
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| ct | str | Required | Cell type of the center cells |
| motifs | Union[str, List[str]] | None | Motifs to test. If None, automatically discovers frequent patterns first |
| max_dist | float | 20 | Maximum radius distance for considering a cell as a neighbor |
| min_size | int | 0 | Minimum neighborhood size for each center cell to be included |
| min_support | float | 0.5 | Minimum frequency threshold for pattern discovery (only used when motifs=None) |
| return_cellID | bool | False | If True, include cell indices in the output for downstream analysis |

## Return Value
`pd.DataFrame` with columns:
- `center`: Center cell type name
- `motifs`: The motif (list of cell type names)
- `n_center_motif`: Number of center cells with the motif in their neighborhood
- `n_center`: Total number of center cells
- `n_motif`: Total number of cells with the motif in their neighborhood
- `expectation`: Expected count under the null (hypergeometric mean)
- `p-values`: Raw p-values from hypergeometric test
- `adj-pval`: FDR-corrected p-values (when multiple motifs tested)
- `if_significant`: Boolean indicating statistical significance

When `return_cellID=True`, additional columns:
- `neighbor_id`: Array of unique cell indices of motif cells in the neighborhood
- `center_id`: Array of center cell indices with the motif

## Usage Example
```python
from SpatialQuery import spatial_query

# Initialize spatial_query object
spatial_key = 'X_spatial'
label_key = 'cell_type'
feature_name = 'gene'

sp = spatial_query(
    adata=adata,
    spatial_key=spatial_key,
    label_key=label_key,
    feature_name=feature_name
)

# Define parameters
ct = 'your_anchor_cell_type'
max_dist = 30
min_support = 0.3

# Auto-discover motifs and test enrichment
result = sp.motif_enrichment_dist(
    ct=ct,
    max_dist=max_dist,
    min_support=min_support
)

# Test specific motifs
motifs = ['cell_type1', 'cell_type2']

result = sp.motif_enrichment_dist(
    ct=ct,
    motifs=motifs,
    max_dist=max_dist
)

# Get cell IDs for downstream DE analysis
motifs_single = ['cell_type1']

result = sp.motif_enrichment_dist(
    ct=ct,
    motifs=motifs_single,
    max_dist=max_dist,
    return_cellID=True
)
```

## Notes
- When `motifs=None`, the method first calls `find_fp_dist()` to discover frequent patterns above the frequency threshold (`min_support`) and returns only maximal patterns, then tests each for enrichment. Due to the frequency threshold and maximal pattern constraints, this automatic/unbiased discovery does not list all significant motifs but serves as a reference. Users can identify cell types appearing in the discovered motifs and customize their own motifs of interest to further test for significance.
- Uses KD-tree `query_ball_point` for efficient radius-based neighbor search.
- P-values are corrected using FDR (positive correlation method) when multiple motifs are tested.
- **Only set `return_cellID=True` when planning downstream DE (`de_genes()`) or gene co-variation (`compute_gene_gene_correlation()`) analysis. For pure motif enrichment queries or visualization, leave it `False`.**
- When computing non-center IDs (anchor cells WITHOUT the motif) for DE analysis, use `np.setdiff1d`:
  ```python
  center_ids = sig['center_id']
  all_anchor_ids = np.where(np.array(adata.obs[label_key]) == ct)[0]
  non_center_ids = np.setdiff1d(all_anchor_ids, center_ids)
  de = sp.de_genes(ind_group1=center_ids, ind_group2=non_center_ids)
  ```
  **Do NOT use subtraction (`all_anchor_ids - center_ids`) — this fails when shapes differ.**

## Troubleshooting
- **Empty result**: The specified cell type may not exist in the data. Check with `sp.get_labels().unique()`.
- **All motifs non-significant**: Try adjusting `max_dist` to change the neighborhood radius.

## Related Methods
- `motif_enrichment_knn()`: KNN-based alternative for motif enrichment
- `find_fp_dist()`: Discover frequent patterns before enrichment
- `de_genes()`: Perform DE using cell IDs from enrichment
- `plot_motif_celltype()`: Visualize motif distribution around a cell type
- `plot_motif_enrichment_heatmap()`: Visualize enrichment results as heatmap
- `compute_gene_gene_correlation()`, `compute_gene_gene_correlation_by_type()`: Compute gene-gene co-varying patterns between anchor cell and motif cells
