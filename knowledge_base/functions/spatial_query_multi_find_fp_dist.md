# find_fp_dist (spatial_query_multi)

## Description
Find frequent patterns within the radius-based neighborhood of a specified cell type across multiple FOVs. Aggregates neighbor transactions from all specified FOVs before applying FP-Growth.

## Function Signature
```python
sp_multi.find_fp_dist(
    ct: str,
    dataset: Union[str, List[str]] = None,
    max_dist: float = 20,
    min_size: int = 0,
    min_support: float = 0.5,
) -> pd.DataFrame
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| ct | str | Required | Cell type name to use as center cells |
| dataset | Union[str, List[str]] | None | Dataset name(s) to include. If None, uses all datasets |
| max_dist | float | 20 | Maximum radius distance for neighbors |
| min_size | int | 0 | Minimum neighborhood size for each center cell |
| min_support | float | 0.5 | Minimum frequency threshold for frequent patterns |

## Return Value
`pd.DataFrame` with columns:
- `support`: Frequency of the pattern across all FOVs
- `itemsets`: Cell type combinations (as sorted list)

## Usage Example
```python
# Define parameters
ct = 'your_anchor_cell_type'
max_dist = 10
min_support = 0.3

# Find patterns across all FOVs
fp = sp_multi.find_fp_dist(
    ct=ct,
    max_dist=max_dist,
    min_support=min_support
)

# Find patterns across all FOVs in specific datasets
dataset = 'condition_A'

fp = sp_multi.find_fp_dist(
    ct=ct,
    dataset=dataset,
    max_dist=max_dist
)
```

## Notes
- Transactions from all FOVs are pooled before running FP-Growth.
- Only maximal frequent patterns are returned.

## Related Methods
- `find_fp_knn()`: KNN-based alternative
- `find_fp_dist_fov()`: Per-FOV pattern discovery
- `motif_enrichment_dist()`: Additional enrichment analysis is performed for discovered motifs
