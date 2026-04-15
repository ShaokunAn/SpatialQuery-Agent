# find_patterns_grid (spatial_query)

## Description
Create a regular grid across the tissue and find frequent cell type patterns in the radius-based neighborhood of each grid point. This is a location-agnostic approach that identifies spatial motifs across the entire FOV without focusing on a specific center cell type.

## Function Signature
```python
sp.find_patterns_grid(
    max_dist: float = 20,
    min_size: int = 0,
    min_support: float = 0.5,
    if_display: bool = True,
    figsize: tuple = (10, 5),
    return_cellID: bool = False,
    return_grid: bool = False,
)
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| max_dist | float | 20 | Grid spacing and radius for neighborhood search |
| min_size | int | 0 | Minimum number of cells in a grid point's neighborhood to be included |
| min_support | float | 0.5 | Minimum frequency threshold to consider a pattern as frequent |
| if_display | bool | True | Whether to display a scatter plot of cells involved in frequent patterns |
| figsize | tuple | (10, 5) | Figure size for the display plot |
| return_cellID | bool | False | If True, include cell indices for each frequent pattern |
| return_grid | bool | False | If True, return grid points along with the result |

## Return Value
- If `return_grid=True`: `pd.DataFrame` sorted by support (descending), with columns:
  - `support`: Frequency of the pattern
  - `itemsets`: Cell type combinations (as frozenset)
  - `neighbor_id`: (only when `return_cellID=True`) Set of cell indices in each pattern

- If `return_grid=True`: Tuple of `(pd.DataFrame, np.ndarray)` where the second element is the grid point coordinates.

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
    feature_name=feature_name
)

# Define parameters
max_dist = 10
min_support = 0.3

# Find global patterns with visualization
fp = sp.find_patterns_grid(
    max_dist=max_dist,
    min_support=min_support,
    if_display=True
)

# Get patterns with cell IDs and grid points
fp, grid = sp.find_patterns_grid(
    max_dist=max_dist,
    min_support=min_support,
    if_display=False,
    return_cellID=True,
    return_grid=True
)
```

## Notes
- The grid spacing equals `max_dist`, creating a regular lattice across the tissue.
- This method is useful for exploratory analysis when you do not have a specific center cell type in mind.
- Only maximal frequent patterns are returned.
- The display shows the spatial distribution of cells belonging to the discovered frequent patterns, color-coded by cell type.

## Troubleshooting
- **No patterns found**: Lower `min_support` or increase `max_dist`.
- **Too many patterns**: Increase `min_support` or decrease `max_dist`.
- **Slow performance**: Increase `max_dist` to reduce the number of grid points.

## Related Methods
- `find_patterns_rand()`: Random-point-based alternative
- `plot_motif_grid()`: Visualize a specific motif around grid points
