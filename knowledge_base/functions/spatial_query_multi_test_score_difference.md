# test_score_difference (spatial_query_multi)

## Description
Test whether gene-pairs have significantly different correlation scores between two conditions. This is the same static method available on both `spatial_query` and `spatial_query_multi` classes. See the shared documentation at `spatial_query_test_score_difference.md` for full details.

## Function Signature
```python
spatial_query_multi.test_score_difference(
    result_A: pd.DataFrame,
    result_B: pd.DataFrame,
    score_col: str = 'combined_score',
    significance_col: str = 'if_significant',
    gene_center_col: str = 'gene_center',
    gene_motif_col: str = 'gene_motif',
    percentile_threshold: float = 95.0,
    background: Literal['Overlapping', 'Significant'] = 'Significant',
) -> pd.DataFrame
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| result_A | pd.DataFrame | Required | Results from condition A |
| result_B | pd.DataFrame | Required | Results from condition B |
| score_col | str | 'combined_score' | Column with correlation scores |
| significance_col | str | 'if_significant' | Column indicating significance |
| gene_center_col | str | 'gene_center' | Column for center gene names |
| gene_motif_col | str | 'gene_motif' | Column for motif gene names |
| percentile_threshold | float | 95.0 | Percentile for outlier identification |
| background | str | 'Significant' | Background set: 'Overlapping' or 'Significant' |

## Return Value
`pd.DataFrame` with columns: `gene_center`, `gene_motif`, `score_A`, `score_B`, `score_diff`, `percentile`, `is_outlier`, `significant_in_A`, `significant_in_B`, `outlier_direction`

## Usage Example
```python
from SpatialQuery import spatial_query_multi

# Define parameters
ct = 'your_anchor_cell_type'
motif = ['cell_type1']
max_dist = 30

# Compute correlations for two conditions
results_condition_A = sp_multi.compute_gene_gene_correlation(
    ct=ct,
    motif=motif,
    dataset='condition_A',
    max_dist=max_dist
)
results_condition_B = sp_multi.compute_gene_gene_correlation(
    ct=ct,
    motif=motif,
    dataset='condition_B',
    max_dist=max_dist
)

# Compare
percentile_threshold = 95.0

diff = spatial_query_multi.test_score_difference(
    result_A=results_condition_A,
    result_B=results_condition_B,
    percentile_threshold=percentile_threshold
)

higher_in_A = diff[diff['outlier_direction'] == 'higher_in_A']
```

## Related Methods
- `compute_gene_gene_correlation()`: Generate input data
- `compute_gene_gene_correlation_by_type()`: Generate type-specific input data
