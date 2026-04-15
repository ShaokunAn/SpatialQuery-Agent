# find_fp_knn_fov (spatial_query_multi)

## Description
Find frequent patterns within the KNN neighborhood of a specific cell type in a single specified FOV. Unlike `find_fp_knn()` which aggregates across FOVs, this method analyzes one FOV at a time. Useful for per-FOV pattern discovery in differential analysis workflows.

## Function Signature
```python
sp_multi.find_fp_knn_fov(
    ct: str,
    dataset_i: str,
    k: int = 30,
    min_support: float = 0.5,
    max_dist: float = 20,
) -> pd.DataFrame
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| ct | str | Required | Cell type name |
| dataset_i | str | Required | Modified dataset name in `dataset_fovindex` format (e.g., 'disease_0') |
| k | int | 30 | Number of nearest neighbors |
| min_support | float | 0.5 | Minimum frequency threshold |
| max_dist | float | 20 | Maximum distance for valid neighbors |

## Return Value
`pd.DataFrame` with columns:
- `support`: Frequency of the pattern
- `itemsets`: Cell type combinations

## Usage Example
```python
# Get available dataset names
print(sp_multi.datasets)  # e.g., ['condition_A_0', 'condition_A_1', 'condition_B_0', 'condition_B_1']

# Define parameters
ct = 'your_anchor_cell_type'
dataset_i = 'condition_A_0'
k = 30
min_support = 0.3

# Find patterns in a specific FOV
fp = sp_multi.find_fp_knn_fov(
    ct=ct,
    dataset_i=dataset_i,
    k=k,
    min_support=min_support
)
```

## Notes
- The `dataset_i` must be the modified dataset name (with FOV index), not the original name.
- Returns empty DataFrame if the cell type does not exist in the specified FOV.
- This analysis is equivalent to creating a `spatial_query` object for the specified FOV and calling `sp.find_fp_knn()`.

## Related Methods
- `find_fp_dist_fov()`: Distance-based per-FOV alternative
- `find_fp_knn()`: Aggregate pattern discovery across FOVs
