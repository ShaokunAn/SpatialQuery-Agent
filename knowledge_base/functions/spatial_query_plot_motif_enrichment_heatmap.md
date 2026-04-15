# plot_motif_enrichment_heatmap (spatial_query)

## Description
Plot a heatmap showing the distribution of cell types in enriched motifs. Provides a visual summary of motif enrichment results from `motif_enrichment_dist()` or `motif_enrichment_knn()`.

## Function Signature
```python
sp.plot_motif_enrichment_heatmap(
    enrich_df: pd.DataFrame,
    figsize: tuple = (7, 5),
    save_path: Optional[str] = None,
    title: Optional[str] = None,
    cmap: str = 'GnBu',
)
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| enrich_df | pd.DataFrame | Required | Output DataFrame from `motif_enrichment_dist()` or `motif_enrichment_knn()` |
| figsize | tuple | (7, 5) | Figure size |
| save_path | Optional[str] | None | Path to save the figure. If None, figure is not saved |
| title | Optional[str] | None | Figure title. If None, uses a default title based on center cell type |
| cmap | str | 'GnBu' | Colormap for the heatmap |

## Return Value
A matplotlib figure showing the heatmap.

## Usage Example
```python
# Define parameters
ct = 'your_anchor_cell_type'
max_dist = 10
min_support = 0.3

# Run enrichment analysis
enrich_result = sp.motif_enrichment_dist(
    ct=ct,
    max_dist=max_dist,
    min_support=min_support
)

# Plot heatmap
title = f'Enriched Motifs surrounding {ct}'
save_path = 'enrichment_heatmap.png'

sp.plot_motif_enrichment_heatmap(
    enrich_df=enrich_result,
    title=title,
    save_path=save_path
)
```

## Related Methods
- `motif_enrichment_knn()` / `motif_enrichment_dist()`: Generate input data
