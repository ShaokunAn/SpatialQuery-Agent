# plot_motif_celltype (spatial_query)

## Description
Display the spatial distribution of a specific motif in the radius-based neighborhood of a given center cell type. This function could be used to visualize significant motifs by `motif_enrichment_dist()` or `motif_enrichment_knn()`, or customized motifs by users.

## Function Signature
```python
sp.plot_motif_celltype(
    ct: str,
    motif: Union[str, List[str]],
    max_dist: float = 20,
    figsize: tuple = (5, 5),
    save_path: Optional[str] = None,
)
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| ct | str | Required | Center cell type |
| motif | Union[str, List[str]] | Required | Motif cell type(s) to display |
| max_dist | float | 20 | Neighborhood radius (should match `motif_enrichment_dist()` parameter) |
| figsize | tuple | (5, 5) | Figure size |
| save_path | Optional[str] | None | Path to save the figure. If None, figure is not saved |

## Return Value
A matplotlib figure.

## Usage Example
```python
# Define parameters
ct = 'your_anchor_cell_type'
motif = ['cell_type1', 'cell_type2']
max_dist = 10

# Visualize motif distribution around center cell type
sp.plot_motif_celltype(
    ct=ct,
    motif=motif,
    max_dist=max_dist
)
```

## Related Methods
- `motif_enrichment_dist()`: Run enrichment analysis first to get enriched motifs
