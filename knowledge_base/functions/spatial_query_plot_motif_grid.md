# plot_motif_grid (spatial_query)

## Description
Display the spatial distribution of a specific motif (cell type combination) around grid points. Shows cells belonging to the specified motif types overlaid on a grid that matches the spacing used in `find_patterns_grid()`.

## Function Signature
```python
sp.plot_motif_grid(
    motif: Union[str, List[str]],
    figsize: tuple = (10, 5),
    max_dist: float = 20,
    save_path: Optional[str] = None,
)
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| motif | Union[str, List[str]] | Required | Cell type name(s) to display |
| figsize | tuple | (10, 5) | Figure size |
| max_dist | float | 20 | Grid spacing distance (should match `find_patterns_grid()` parameter) |
| save_path | Optional[str] | None | Path to save the figure. If None, figure is not saved |

## Return Value
A matplotlib figure.

## Usage Example
```python
# Define parameters
motif = ['cell_type1', 'cell_type2']
max_dist = 10
figsize = (12, 6)

# Plot motif distribution on grid
sp.plot_motif_grid(
    motif=motif,
    max_dist=max_dist,
    figsize=figsize
)
```

## Related Methods
- `find_patterns_grid()`: Discover patterns before visualization
- `plot_motif_rand()`: Random-point-based visualization alternative
