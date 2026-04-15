# plot_all_center_motif (spatial_query)

## Description
Plot the spatial distribution of center cells, their neighboring motif cells, distant motif cells, and anchor cells without motif. This visualization helps understand the spatial context of the three correlation groups used in `compute_gene_gene_correlation()`.

## Function Signature
```python
sp.plot_all_center_motif(
    ct: str,
    ids: dict,
    figsize: tuple = (6, 6),
    save_path: Optional[str] = None,
)
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| ct | str | Required | Center cell type |
| ids | dict | Required | Cell group dictionary returned by `compute_gene_gene_correlation()` (the second return value) |
| figsize | tuple | (6, 6) | Figure size |
| save_path | Optional[str] | None | Path to save the figure. If None, figure is not saved |

## Return Value
A matplotlib figure.

## Usage Example
```python
# Define parameters
ct = 'your_anchor_cell_type'
motif = ['cell_type1', 'cell_type2']
max_dist = 10

# First compute correlations
results_df, cell_groups = sp.compute_gene_gene_correlation(
    ct=ct,
    motif=motif,
    max_dist=max_dist # or kNN neighborhood with k
)

# Visualize the cell groups
figsize = (8, 8)
save_path = 'all_center_motif.png'

sp.plot_all_center_motif(
    ct=ct,
    ids=cell_groups,
    figsize=figsize,
    save_path=save_path  # Optional. If None, the figure will not be saved.
)
```

## Notes
The `ids` dictionary should contain:
- `center_neighbor_motif_pair`: Center-neighbor pairs for correlation 1
- `non-neighbor_motif_cells`: Distant motif cells for correlation 2
- `non_motif_center_neighbor_pair`: Center-neighbor pairs without motif for correlation 3

## Related Methods
- `compute_gene_gene_correlation()`: Generate the `ids` input for this function
