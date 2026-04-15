# compute_gene_gene_correlation (spatial_query)

## Description
Compute gene-gene cross correlation between anchor and neighboring motif cells in a single FOV. This analysis identifies co-varying gene pairs between center cells and their neighboring motif cells to reveal spatial interaction-dependent gene regulatory patterns.

The method computes three types of correlations:
1. **Correlation 1 (neighbor)**: Between center cells and their neighboring motif cells (excluding center type from neighbors)
2. **Correlation 2 (non-neighbor)**: Between center cells with motif and distant motif cells (background)
3. **Correlation 3 (no-motif)**: Between center cells without the motif and their neighbors (control)

Uses shifted correlation (subtracting cell type means) to enable comparison across different niches/motifs.

## Function Signature
```python
sp.compute_gene_gene_correlation(
    ct: str,
    motif: Union[str, List[str]],
    genes: Optional[Union[str, List[str]]] = None,
    max_dist: Optional[float] = None,
    k: Optional[int] = None,
    min_size: int = 0,
    min_nonzero: int = 10,
    alpha: Optional[float] = None,
) -> Tuple[pd.DataFrame, dict]
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| ct | str | Required | Cell type of the center cells |
| motif | Union[str, List[str]] | Required | Motif cell type(s) to analyze. Include all cell types for neighbor finding |
| genes | Optional[Union[str, List[str]]] | None | Genes to analyze. If None, uses all genes |
| max_dist | Optional[float] | None | Maximum distance for radius-based neighbor search. Use either max_dist or k |
| k | Optional[int] | None | Number of nearest neighbors. Use either max_dist or k |
| min_size | int | 0 | Minimum neighborhood size (only used with max_dist) |
| min_nonzero | int | 10 | Minimum non-zero expression values required for a gene to be included |
| alpha | Optional[float] | None | Significance threshold for multi-testing correction |

## Return Value
Tuple of `(results_df, cell_groups)`:

**results_df** (`pd.DataFrame`):
- `gene_center`: Gene in center cells
- `gene_motif`: Gene in motif cells of the specific type
- `corr_neighbor`: Correlation with neighboring cells of this type
- `corr_non_neighbor`: Correlation with distant cells of this type
- `corr_center_no_motif`: Correlation with neighbors when no motif is present
- `p_value_test1`: P-value for test1 (neighbor vs non-neighbor)
- `p_value_test2`: P-value for test2 (neighbor vs no_motif)
- `q_value_test1`: FDR-corrected q-value for test1
- `q_value_test2`: FDR-corrected q-value for test2
- `delta_corr_test1`: Correlation difference (neighbor - non_neighbor)
- `delta_corr_test2`: Correlation difference (neighbor - no_motif)
- `reject_test1_fdr`: Whether test1 passes FDR threshold
- `reject_test2_fdr`: Whether test2 passes FDR threshold
- `combined_score`: Combined significance score
- `abs_combined_score`: Absolute value of combined score
- `if_significant`: Whether both tests pass FDR threshold

**cell_groups** (`dict`):
- `center_neighbor_motif_pair`: Array of shape (n_pairs, 2) with [center_idx, neighbor_idx] pairs. Each row is [center_cell_idx, neighbor_cell_idx].
- `non-neighbor_motif_cells`: Array of distant motif cell indices. Correlation 2 uses all combinations of center cells (from corr1) × these cells.
- `non_motif_center_neighbor_pair`: Array of shape (n_pairs, 2) for centers without motif. Each row is [center_cell_idx, neighbor_cell_idx]. 

Individual cell IDs can be extracted from pairs using np.unique() like:
- center_cells = np.unique(center_neighbor_motif_pair[:, 0])
- neighbor_cells = np.unique(center_neighbor_motif_pair[:, 1])

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
    feature_name=feature_name,
    build_gene_index=False,  # use raw data by default with build_gene_index=False, else use scfind index
    if_lognorm=True  # Set true if build_gene_index=False and adata.X is raw count data
)

# Define parameters
ct = 'your_anchor_cell_type'
motif = ['cell_type1', 'cell_type2']
genes = ['gene1', 'gene2', 'gene3', 'gene4']
max_dist = 5
min_nonzero = 10
min_size = 0
alpha = 0.05

# Analyze gene co-variation using radius-based neighborhood
results_df, cell_groups = sp.compute_gene_gene_correlation(
    ct=ct,
    motif=motif,
    genes=genes,  # Optional. Use all genes if not specified
    max_dist=max_dist,
    min_size=min_size,
    min_nonzero=min_nonzero,
    alpha=alpha
)

# If using kNN neighborhood
k = 10

results_df, cell_groups = sp.compute_gene_gene_correlation(
    ct=ct,
    motif=motif,
    genes=genes,  # Optional. Use all genes if not specified
    k=k,
    min_size=min_size,
    min_nonzero=min_nonzero,
    alpha=alpha
)

# Filter significant results
sig_results = results_df[results_df['if_significant'] == True]
```

## Notes
- Must specify either `max_dist` or `k`, but not both. `max_dist` is used for radius-based neighbor search, while `k` is used for kNN.
- When `build_gene_index=True`, uses binary expression data (scfind index) and `alpha=0.1` is recommended since binary data results tend to be conservative. When `build_gene_index=False`, uses continuous expression values and `alpha=0.05` is recommended.
- The shifted correlation approach subtracts cell type means before computing correlations, enabling fair comparison across different spatial contexts.
- Only inter-cell-type interactions are considered (center type cells are removed from neighbor groups).

## Troubleshooting
- **No results**: Ensure there are enough cells of both center and motif types. Check `min_size` threshold and try various neighborhood size.
- **All non-significant**: Consider relaxing `alpha` or ensuring sufficient cell numbers.

## Related Methods
- `get_anchor_motif_cell_ids()`: Get cell IDs only (same `cell_groups` dict) without computing correlations — use this when you only need cell IDs for downstream analysis (DE, plotting, etc.)
- `compute_gene_gene_correlation_by_type()`: Separate analysis for each cell type in the motif
- `test_score_difference()`: Compare correlation results between two conditions
- `plot_gene_pair_heatmap()`: Visualize gene pair correlation results
- `plot_all_center_motif()`: Visualize the spatial distribution of center/anchor-motif cell groups
