# find_fp_knn (spatial_query_multi)

## Description
Find frequent patterns within the k-nearest neighbors of a specified cell type across multiple FOVs. Aggregates neighbor transactions from all specified FOVs before applying FP-Growth, giving a global view of co-localization patterns.

## Function Signature
```python
sp_multi.find_fp_knn(
    ct: str,
    dataset: Union[str, List[str]] = None,
    k: int = 30,
    min_support: float = 0.5,
    max_dist: float = 20,
) -> pd.DataFrame
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| ct | str | Required | Cell type name to use as center cells |
| dataset | Union[str, List[str]] | None | Dataset name(s) to include. If None, uses all datasets |
| k | int | 30 | Number of nearest neighbors |
| min_support | float | 0.5 | Minimum frequency threshold for frequent patterns |
| max_dist | float | 20 | Maximum distance for valid neighbors |

## Return Value
`pd.DataFrame` with columns:
- `support`: Frequency of the pattern across all FOVs
- `itemsets`: Cell type combinations (as sorted list)

## Usage Example
```python
# Define parameters
ct = 'your_anchor_cell_type'
k = 30
min_support = 0.3

# Find patterns across all FOVs
fp = sp_multi.find_fp_knn(
    ct=ct,
    k=k,
    min_support=min_support
)

# Find patterns only in a specific condition
dataset = 'condition_A'

fp = sp_multi.find_fp_knn(
    ct=ct,
    dataset=dataset,
    k=k,
    min_support=min_support
)

# Find patterns in multiple specific datasets
datasets = ['condition_A', 'condition_B']

fp = sp_multi.find_fp_knn(
    ct=ct,
    dataset=datasets,
    min_support=min_support
)
```

## Notes
- Transactions are collected from all FOVs of the specified datasets, then FP-Growth is applied on the combined set.
- Only maximal frequent patterns are returned.

## Related Methods
- `find_fp_dist()`: Distance-based alternative
- `find_fp_knn_fov()`: Per-FOV pattern discovery
- `motif_enrichment_knn()`: Additional enrichment analysis is performed for discovered motifs
