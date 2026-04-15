# plot_gene_pair_heatmap (spatial_query)

## Description
Plot a heatmap showing co-varying gene pairs from gene-gene correlation analysis. Visualizes the results from `compute_gene_gene_correlation()` or `compute_gene_gene_correlation_by_type()`.

## Function Signature
```python
sp.plot_gene_pair_heatmap(
    gene_pair_df: pd.DataFrame,
    figsize: tuple = (7, 5),
    save_path: Optional[str] = None,
    cmap: str = 'GnBu',
)
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| gene_pair_df | pd.DataFrame | Required | Output DataFrame from `compute_gene_gene_correlation()` or `compute_gene_gene_correlation_by_type()` |
| figsize | tuple | (7, 5) | Figure size |
| save_path | Optional[str] | None | Path to save the figure. If None, figure is not saved |
| cmap | str | 'GnBu' | Colormap for the heatmap |

## Return Value
A matplotlib figure showing the gene pair heatmap.

## Usage Example
```python
# Define parameters
ct = 'your_anchor_cell_type'
motif = ['cell_type1', 'cell_type2']
max_dist = 10

# Compute gene-gene correlation between anchor cells and each type of motif cells
results_df = sp.compute_gene_gene_correlation_by_type(
    ct=ct,
    motif=motif,
    max_dist=max_dist
)

# Plot heatmap of significant gene pairs
sig_results_ct1 = results_df[(results_df['if_significant'] == True) & (results_df['cell_type'] == 'cell_type1')]

save_path = 'gene_pair_heatmap.png'

sp.plot_gene_pair_heatmap(
    gene_pair_df=sig_results_ct1,
    save_path=save_path
)

# Compute gene-gene correlation between anchor cells and pooled motif cells
results_df = sp.compute_gene_gene_correlation(
    ct=ct,
    motif=motif,
    max_dist=max_dist
)

sig_results = results_df[results_df['if_significant'] == True]

title = 'Gene Co-variation Heatmap'
save_path = 'gene_pair_heatmap.png'

sp.plot_gene_pair_heatmap(
    gene_pair_df=sig_results,
    title=title,
    save_path=save_path
)
```

## Related Methods
- `compute_gene_gene_correlation()`: Generate input data
- `compute_gene_gene_correlation_by_type()`: Generate type-specific input data
