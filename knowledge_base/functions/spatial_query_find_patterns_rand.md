# find_patterns_rand (spatial_query)

## Description
Randomly generate reference points across the tissue and find frequent cell type patterns in the radius-based neighborhood of each random point. This provides a stochastic, location-agnostic approach to identify spatial motifs across the entire FOV without focusing on a specific center cell type.

## Function Signature
```python
sp.find_patterns_rand(
    max_dist: float = 20,
    n_points: int = 1000,
    min_support: float = 0.5,
    min_size: int = 0,
    if_display: bool = True,
    figsize: tuple = (10, 5),
    return_cellID: bool = False,
    seed: int = 2023,
) -> pd.DataFrame
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| max_dist | float | 20 | Radius for neighborhood search around each random point |
| n_points | int | 1000 | Number of random points to generate |
| min_support | float | 0.5 | Minimum frequency threshold to consider a pattern as frequent |
| min_size | int | 0 | Minimum number of cells in a point's neighborhood to be included |
| if_display | bool | True | Whether to display a scatter plot of cells in frequent patterns |
| figsize | tuple | (10, 5) | Figure size for the display plot |
| return_cellID | bool | False | If True, include cell indices for each frequent pattern |
| seed | int | 2023 | Random seed for reproducibility |

## Return Value
`pd.DataFrame` sorted by support (descending), with columns:
- `support`: Frequency of the pattern
- `itemsets`: Cell type combinations (as frozenset)
- `neighbor_id`: (only when `return_cellID=True`) Set of cell indices in each pattern

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
max_dist = 5
n_points = 2000
min_support = 0.3
seed = 42

# Find patterns using random points
fp = sp.find_patterns_rand(
    max_dist=max_dist,
    n_points=n_points,
    min_support=min_support,
    seed=seed
)

# Without visualization, with cell IDs
n_points_alt = 1000

fp = sp.find_patterns_rand(
    max_dist=max_dist,
    n_points=n_points_alt,
    min_support=min_support,
    if_display=False,
    return_cellID=True
)
```

## Notes
- Random points are uniformly distributed within the bounding box of the tissue.
- Use the `seed` parameter to ensure reproducible results.
- Only maximal frequent patterns are returned.
- Compared to `find_patterns_grid()`, this method avoids grid artifacts and may better capture patterns at irregular tissue boundaries.

## Troubleshooting
- **No patterns found**: Lower `min_support`, increase `n_points`, or increase `max_dist`.
- **Results vary between runs**: Set the same `seed` for reproducibility.
- **Slow performance**: Reduce `n_points`.

## Related Methods
- `find_patterns_grid()`: Grid-based alternative
- `plot_motif_rand()`: Visualize a specific motif around random points
