# compute_gene_gene_correlation (spatial_query_multi)

## Description
Compute gene-gene co-varying patterns between motif and center cells across multiple FOVs. Aggregates center-neighbor pairs from all specified FOVs, using FOV-specific cell type means for centering (not global means) to ensure fair cross-FOV comparisons.

This function calculates cross correlation between gene expression in:
1. **Correlation 1 (neighbor)**: Motif cells that are neighbors of center cell type (paired across FOVs)
2. **Correlation 2 (non-neighbor)**: Motif cells that are NOT neighbors of center cell type (all-to-all across FOVs)
3. **Correlation 3 (no-motif)**: Neighboring cells of center cell type without nearby motif (paired across FOVs)

## Function Signature
```python
sp_multi.compute_gene_gene_correlation(
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
| motif | Union[str, List[str]] | Required | Motif cell type(s) to analyze |
| dataset | Union[str, List[str]] | None | Datasets to include. If None, uses all |
| genes | Optional[Union[str, List[str]]] | None | Genes to analyze. If None, uses intersection across FOVs |
| max_dist | Optional[float] | None | Maximum radius for neighbors. Use either max_dist or k |
| k | Optional[int] | None | Number of nearest neighbors. Use either max_dist or k |
| min_size | int | 0 | Minimum neighborhood size (only with max_dist) |
| min_nonzero | int | 10 | Minimum non-zero values for gene inclusion |
| alpha | Optional[float] | None | Significance threshold |

## Return Value
**results_df** (`pd.DataFrame`):
- `gene_center`, `gene_motif`: Gene pair names
- `corr_neighbor`: Correlation in neighbor group (across FOVs)
- `corr_non_neighbor`: Correlation in non-neighbor group
- `corr_center_no_motif`: Correlation for centers without motif
- `p_value_test1`: P-value for neighbor vs non-neighbor
- `p_value_test2`: P-value for with motif vs without motif
- `delta_corr_test1`, `delta_corr_test2`: Correlation differences
- `combined_score`: Combined significance score
- `q_value_test1`, `q_value_test2`: FDR-corrected q-values
- `reject_test1_fdr`, `reject_test2_fdr`: Whether each test passes FDR threshold
- `abs_combined_score`: Absolute value of combined score
- `if_significant`: Whether both tests pass FDR threshold

Note: Unlike the single-FOV version, the multi-FOV `compute_gene_gene_correlation` returns only the DataFrame, not a tuple of (DataFrame, dict).

## Usage Example
```python
# Define parameters
ct = 'your_anchor_cell_type'
motif = ['cell_type1', 'cell_type2']
dataset = 'condition_A'
max_dist = 10

# Analyze gene co-variation across all condition_A FOVs
results_df = sp_multi.compute_gene_gene_correlation(
    ct=ct,
    motif=motif,
    dataset=dataset,
    max_dist=max_dist
)

# Using KNN across all FOVs
motif_multi = ['cell_type1', 'cell_type2']
k = 30

results_df = sp_multi.compute_gene_gene_correlation(
    ct=ct,
    motif=motif_multi,
    k=k
)

# With specific genes
genes = ['gene1', 'gene2', 'gene3']

results_df = sp_multi.compute_gene_gene_correlation(
    ct=ct,
    motif=motif,
    genes=genes,
    max_dist=max_dist
)
```

## Notes
- FOV-specific centering ensures that batch effects between FOVs are minimized.
- Must specify either `max_dist` or `k`, but not both. `max_dist` is used for radius-based neighbor search, while `k` is used for kNN.
- When `build_gene_index=True`, uses binary expression data (scfind index) and `alpha=0.1` is recommended since binary data results tend to be conservative. When `build_gene_index=False`, uses continuous expression values and `alpha=0.05` is recommended.
- Only inter-cell-type interactions are considered (center type cells are removed from neighbor groups).

## Related Methods
- `compute_gene_gene_correlation_by_type()`: Per-cell-type analysis in motif 
- `test_score_difference()`: Compare results between groups
