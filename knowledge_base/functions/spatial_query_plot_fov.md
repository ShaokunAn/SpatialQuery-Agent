# plot_fov (spatial_query)

## Description
Plot the spatial distribution of cell types in a single FOV. Displays a scatter plot with cells color-coded by their cell type annotation.

## Function Signature
```python
sp.plot_fov(
    min_cells_label: int = 50,
    title: str = 'Spatial distribution of cell types',
    figsize: tuple = (10, 5),
    save_path: Optional[str] = None,
)
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| min_cells_label | int | 50 | Minimum number of cells for a cell type to be displayed in the legend |
| title | str | 'Spatial distribution of cell types' | Figure title |
| figsize | tuple | (10, 5) | Figure size |
| save_path | Optional[str] | None | Path to save the figure. If None, figure is not saved |

## Return Value
A matplotlib figure showing the spatial scatter plot of cell types.

## Usage Example
```python
# Define parameters
min_cells_label = 10
title = 'Spatial distribution of cell types'
figsize = (12, 6)
save_path = 'fov_plot.png'

# Plot FOV
sp.plot_fov(
    min_cells_label=min_cells_label,
    title=title,
    figsize=figsize,
    save_path=save_path
)
```

## Related Methods
- `plot_motif_grid()`, `plot_motif_rand()`, `plot_motif_celltype()`: Motif-specific visualizations
