# plot_cell_type_distribution_fov (spatial_query_multi)

## Description
Visualize the distribution of cell types across individual FOVs within a single dataset using a stacked bar plot. Useful for assessing FOV-to-FOV variability within a condition.

## Function Signature
```python
sp_multi.plot_cell_type_distribution_fov(
    dataset: str,
    data_type: str = 'number',
    colormap: str = 'tab20c',
    save_path: str = None,
)
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| dataset | str | Required | Dataset name (original name, e.g., 'disease') |
| data_type | str | 'number' | Plot by 'number' (absolute counts) or 'proportion' |
| colormap | str | 'tab20c' | Matplotlib colormap name |
| save_path | str | None | Path to save the figure |

## Return Value
A stacked bar plot figure.

## Usage Example
```python
# Define parameters
dataset = 'condition_A'
data_type = 'number'

# Plot cell count distribution across FOVs in the dataset
sp_multi.plot_cell_type_distribution_fov(
    dataset=dataset,
    data_type=data_type
)

# Plot proportions
data_type_prop = 'proportion'
save_path = 'condition_A_fov_dist.png'

sp_multi.plot_cell_type_distribution_fov(
    dataset=dataset,
    data_type=data_type_prop,
    save_path=save_path
)
```

## Related Methods
- `plot_cell_type_distribution()`: Cross-dataset distribution
