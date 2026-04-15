# test_score_difference (spatial_query / spatial_query_multi)

## Description
Identify gene pairs with large score differences between two covariation pattern results.
This function compares covariation scores between two groups (e.g., disease vs control,
treatment vs baseline, covarying gene pairs by distinct motif types) and identifies gene 
pairs with the largest score differences using percentile-based ranking. Gene pairs in the 
top percentile_threshold% and bottom (100 - percentile_threshold)% of score differences are 
flagged as emperical significant ones.

Available as a static method on both `spatial_query` and `spatial_query_multi` classes.

## Function Signature
```python
spatial_query.test_score_difference(
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
| result_A | pd.DataFrame | Required | Results from `compute_gene_gene_correlation` or `compute_gene_gene_correlation_by_type` for condition A |
| result_B | pd.DataFrame | Required | Results from the same method for condition B |
| score_col | str | 'combined_score' | Column name containing correlation scores to compare |
| significance_col | str | 'if_significant' | Column name indicating significance status |
| gene_center_col | str | 'gene_center' | Column name for center gene names |
| gene_motif_col | str | 'gene_motif' | Column name for motif gene names |
| percentile_threshold | float | 95.0 | Percentile threshold for outlier identification (e.g., 95 = top/bottom 5%) |
| background | str | 'Significant' | Background set for comparison: 'Overlapping' (all overlapping gene pairs) or 'Significant' (only significant pairs) |

## Return Value
`pd.DataFrame` with columns:
- `gene_center`: Center gene name
- `gene_motif`: Motif gene name
- `score_A`: Score in condition A
- `score_B`: Score in condition B
- `score_diff`: Difference (score_A - score_B)
- `percentile`: Percentile rank of score_diff in the distribution
- `is_outlier`: Whether this pair is an outlier
- `significant_in_A`: Whether the pair is significant in condition A
- `significant_in_B`: Whether the pair is significant in condition B
- `outlier_direction`: 'higher_in_A' (>95th percentile), 'lower_in_A' (<5th percentile), or 'not_outlier'

## Usage Example
```python
from SpatialQuery import spatial_query

# Define two motifs to explore covariation patterns under various spatial contexts
ct = 'your_anchor_cell_type'
motif1 = ['cell_type1', 'cell_type2']
motif2 = ['cell_type1', 'cell_type4']
max_dist = 10

# Compute correlations between anchor cells and cell_type1 under various spatial contexts (assuming sp is initialized)
results1 = sp.compute_gene_gene_correlation_by_type(
    ct=ct,
    motif=motif1,
    max_dist=max_dist
)
results2 = sp.compute_gene_gene_correlation_by_type(
    ct=ct,
    motif=motif2,
    max_dist=max_dist
)

results1_ct1 = results1[results1['cell_type'] == 'cell_type1']
results2_ct1 = results2[results2['cell_type'] == 'cell_type1']

# Compare the two conditions
percentile_threshold = 95.0
background = 'Overlapping'

diff_results = spatial_query.test_score_difference(
    result_A=results1_ct1,
    result_B=results2_ct1,
    percentile_threshold=percentile_threshold,
    background=background
)

# Get gene pairs that are 1) significant between anchor cells and neighboring motif cells, and 2) significantly higher in group A than in group B
higher_in_A = diff_results[(diff_results['outlier_direction'] == 'higher_in_A') & (diff_results['significant_in_A'] == True)]
# Or gene pairs that are 1) significant between anchor cells and neighboring motif cells, and 2) significantly lower in group A than in group B
lower_in_A = diff_results[(diff_results['outlier_direction'] == 'lower_in_A') & (diff_results['significant_in_A'] == True)]
```

## Notes
- This is a static method, so it can be called on the class directly without an instance.
- The `background` parameter controls which gene pairs form the reference distribution:
  - `'Significant'`: Only uses gene pairs significant in at least one condition as background. Tends to output fewer gene pairs.
  - `'Overlapping'`: Uses all gene pairs present in both results. Tends to output more gene pairs.
- The gene pairs are ranked by the absolute value of the score difference between the two conditions.

## Related Methods
- `compute_gene_gene_correlation()`: Generate input data for this method
- `compute_gene_gene_correlation_by_type()`: Generate type-specific input data
