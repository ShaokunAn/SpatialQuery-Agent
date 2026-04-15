# compute_gene_gene_correlation_by_type (spatial_query_multi)

## Description
Compute gene-gene cross correlation separately for each cell type in the motif across multiple FOVs. 
Unlike `compute_gene_gene_correlation()` which pools all motif cell types together, this method breaks down the analysis by individual non-center cell types in the motif, aggregating statistics across FOVs.

For each non-center cell type in the motif:
- **Correlation 1**: Center cells with motif vs neighboring motif cells of THIS type
- **Correlation 2**: Center cells with motif vs distant motif cells of THIS type
- **Correlation 3**: Center cells without motif vs their neighbors

## Function Signature
```python
sp_multi.compute_gene_gene_correlation_by_type(
    ct: str,
    motif: Union[str, List[str]],
    dataset: Union[str, List[str]] = None,
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
| ct | str | Required | Center cell type |
| motif | Union[str, List[str]] | Required | Motif cell type(s) (should contain >= 2 non-center types) |
| dataset | Union[str, List[str]] | None | Datasets to include. If None, uses all |
| genes | Optional[Union[str, List[str]]] | None | Genes to analyze. If None, uses intersection across FOVs |
| max_dist | Optional[float] | None | Maximum radius for neighbors. Use either max_dist or k |
| k | Optional[int] | None | Number of nearest neighbors. Use either max_dist or k |
| min_size | int | 0 | Minimum neighborhood size (only with max_dist) |
| min_nonzero | int | 10 | Minimum non-zero values for gene inclusion |
| alpha | Optional[float] | None | Significance threshold |

## Return Value
`pd.DataFrame` with columns:
- `cell_type`: Non-center cell type being analyzed
- `gene_center`, `gene_motif`: Gene pair names
- `corr_neighbor`: Correlation with neighboring cells of this type
- `corr_non_neighbor`: Correlation with distant cells of this type
- `corr_center_no_motif`: Correlation with neighbors when no motif present
- `p_value_test1`, `p_value_test2`: P-values for statistical tests
- `q_value_test1`, `q_value_test2`: FDR-corrected q-values
- `delta_corr_test1`, `delta_corr_test2`: Correlation differences
- `reject_test1_fdr`, `reject_test2_fdr`: FDR significance flags (bool)
- `combined_score`, `abs_combined_score`: Combined significance scores
- `if_significant`: Whether both tests pass FDR threshold (`reject_test1_fdr & reject_test2_fdr`)

## Notes
- This method is designed for motifs with multiple non-center cell types where you want to understand type-specific coavariation.
- FOV-specific centering ensures that batch effects between FOVs are minimized.
- Must specify either `max_dist` or `k`, but not both. `max_dist` is used for radius-based neighbor search, while `k` is used for kNN.
- When `build_gene_index=True`, uses binary expression data (scfind index) and `alpha=0.1` is recommended since binary data results tend to be conservative. When `build_gene_index=False`, uses continuous expression values and `alpha=0.05` is recommended.
- Only inter-cell-type interactions are considered (center type cells are removed from neighbor groups).

## Usage Example
```python
# Define parameters
ct = 'your_anchor_cell_type'
motif = ['cell_type1', 'cell_type2', 'cell_type3']
dataset = 'condition_A'
max_dist = 10

# Analyze per-type contributions in a multi-type motif
results = sp_multi.compute_gene_gene_correlation_by_type(
    ct=ct,
    motif=motif,
    dataset=dataset,
    max_dist=max_dist
)

# Filter significant results
sig_results = results[results['if_significant']]

# Filter for a specific cell type
filter_cell_type = 'cell_type1'
ct1_results = results[(results['cell_type'] == filter_cell_type) & results['if_significant']]
```

## Related Methods
- `compute_gene_gene_correlation()`: Pooled analysis across all motif types
- `test_score_difference()`: Compare results between conditions
