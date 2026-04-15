# plot_motif_rand (spatial_query)

## Description
Display randomly sampled reference points that have a specified motif in their radius-based neighborhood, along with the motif cells in those neighborhoods.

## Function Signature
```python
sp.plot_motif_rand(
    motif: Union[str, List[str]],
    max_dist: float = 100,
    n_points: int = 1000,
    figsize: tuple = (10, 5),
    seed: int = 2023,
    save_path: Optional[str] = None,
)
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| motif | Union[str, List[str]] | Required | Cell type name(s) to display |
| max_dist | float | 100 | Radius for neighborhood search around random points |
| n_points | int | 1000 | Number of random points to generate |
| figsize | tuple | (10, 5) | Figure size |
| seed | int | 2023 | Random seed for reproducibility |
| save_path | Optional[str] | None | Path to save the figure. If None, figure is not saved |

## Return Value
A matplotlib figure.

## Usage Example
```python
# Define parameters
motif = ['cell_type1', 'cell_type2']
max_dist = 10
n_points = 2000
seed = 42

# Plot motif distribution with random points
sp.plot_motif_rand(
    motif=motif,
    max_dist=max_dist,
    n_points=n_points,
    seed=seed
)
```

## Related Methods
- `find_patterns_rand()`: Discover patterns before visualization
- `plot_motif_grid()`: Grid-based visualization alternative
