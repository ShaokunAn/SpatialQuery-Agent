# compute_gene_gene_correlation_by_type (spatial_query)

## Description
Compute gene-gene cross correlation separately for each cell type in the motif within a single FOV. It computes a shifted Pearson correlation using cell type-specific global means as reference baselines in each field of view (FOV), then compare the observed correlation against two controls—a non-spatial baseline pairing anchor cells with distal cells of the same motif cell types, and a non-motif background from anchor cells lacking the specific motif. Fisher's Z-transformation is applied to test whether correlations differ significantly between contexts. Unlike `compute_gene_gene_correlation()` which pools all motif cell types together, this method breaks down the analysis by individual non-center cell types in the motif.

For each non-center cell type in the motif:
- **Correlation 1**: Center cells with motif vs neighboring motif cells of THIS specific type
- **Correlation 2**: Center cells with motif vs distant motif cells of THIS specific type
- **Correlation 3**: Center cells without motif vs their neighbors (same for all types)

Only analyzes motifs with cell types besides the center type. If there is only one cell type except the center type in the motif, this computation is the same as `compute_gene_gene_correlation()`.

## Function Signature
```python
sp.compute_gene_gene_correlation_by_type(
    ct: str,
    motif: Union[str, List[str]],
    genes: Optional[Union[str, List[str]]] = None,
    max_dist: Optional[float] = None,
    k: Optional[int] = None,
    min_size: int = 0,
    min_nonzero: int = 10,
    alpha: Optional[float] = None,
) -> pd.DataFrame
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| ct | str | Required | Cell type of the center cells |
| motif | Union[str, List[str]] | Required | Motif cell type(s) to analyze (should contain >= 2 non-center types) |
| genes | Optional[Union[str, List[str]]] | None | Genes to analyze. If None, uses all genes |
| max_dist | Optional[float] | None | Maximum distance for radius-based neighbor search. Use either max_dist or k |
| k | Optional[int] | None | Number of nearest neighbors. Use either max_dist or k |
| min_size | int | 0 | Minimum neighborhood size (only used with max_dist) |
| min_nonzero | int | 10 | Minimum non-zero expression values required for a gene |
| alpha | Optional[float] | None | Significance threshold for multi-testing correction |

## Return Value
`pd.DataFrame` with columns:
- `cell_type`: The non-center cell type in motif
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

## Usage Example
```python
from SpatialQuery import spatial_query

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

# Analyze 
ct = 'your_anchor_cell_type'
motif = ['cell_type1', 'cell_type2', 'cell_type3']
genes = ['gene1', 'gene2', 'gene3'] 
max_dist = 5
# k = 10 # Optional. Use kNN neighborhood instead of radius-based neighborhood
min_size = 0
min_nonzero = 10
alpha = 0.05

results = sp.compute_gene_gene_correlation_by_type(
    ct=ct
    motif=motif,
    genes=genes,  # Optional. Use all features if not specified
    max_dist=max_dist,
    min_size=min_size,
    min_nonzero=min_nonzero,
    alpha=alpha
)

# If using kNN neighborhood
k = 10
results = sp.compute_gene_gene_correlation_by_type(
    ct=ct
    motif=motif,
    genes=genes,  # Optional. Use all features if not specified
    k=k,
    min_size=min_size,
    min_nonzero=min_nonzero,
    alpha=alpha
)

# Filter for significant results of a specific cell type
ct1_results = results[
    (results['cell_type'] == 'cell_type1') & (results['if_significant'] == True)
]
```

## Notes
- This method is designed for motifs with multiple non-center cell types where you want to understand type-specific coavariation.
- Must specify either `max_dist` or `k`, but not both. `max_dist` is used for radius-based neighbor search, while `k` is used for kNN.
- When `build_gene_index=True`, uses binary expression data (scfind index) and `alpha=0.1` is recommended since binary data results tend to be conservative. When `build_gene_index=False`, uses continuous expression values and `alpha=0.05` is recommended.
- Only inter-cell-type interactions are considered (center type cells are removed from neighbor groups).

## Related Methods
- `compute_gene_gene_correlation()`: Pooled analysis across all motif types
- `test_score_difference()`: Ranking covarying gene pairs between two covariation results. 
- `plot_gene_pair_heatmap()`: Visualize gene pair results
- `plot_all_center_motif()`: Visualize the spatial distribution of center/anchor-motif cell groups
