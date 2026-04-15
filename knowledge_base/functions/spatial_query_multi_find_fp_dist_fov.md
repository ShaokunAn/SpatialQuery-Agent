# find_fp_dist_fov (spatial_query_multi)

## Description
Find frequent patterns within the radius-based neighborhood of a specific/anchor cell type in a single specified FOV. Analyzes one FOV at a time for per-FOV pattern discovery.

## Function Signature
```python
sp_multi.find_fp_dist_fov(
    ct: str,
    dataset_i: str,
    max_dist: float = 20,
    min_size: int = 0,
    min_support: float = 0.5,
) -> pd.DataFrame
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| ct | str | Required | Cell type name |
| dataset_i | str | Required | Modified dataset name in `dataset_fovindex` format (e.g., 'disease_0') |
| max_dist | float | 20 | Maximum radius for neighbors |
| min_size | int | 0 | Minimum neighborhood size |
| min_support | float | 0.5 | Minimum frequency threshold |

## Return Value
`pd.DataFrame` with columns:
- `support`: Frequency of the pattern
- `itemsets`: Cell type combinations

## Usage Example
```python
# Define parameters
ct = 'your_anchor_cell_type'
dataset_i = 'condition_A_0'
max_dist = 5
min_support = 0.3

# Find patterns in a specific FOV
fp = sp_multi.find_fp_dist_fov(
    ct=ct,
    dataset_i=dataset_i,
    max_dist=max_dist,
    min_support=min_support
)
```

## Notes
- The `dataset_i` must be the modified dataset name (with FOV index).
- Returns empty DataFrame if the cell type does not exist in the specified FOV.
- This analysis is equivalent to creating a `spatial_query` object for the specified FOV and calling `sp.find_fp_dist()`.

## Related Methods
- `find_fp_knn_fov()`: KNN-based per-FOV alternative
- `find_fp_dist()`: Aggregate pattern discovery across FOVs
