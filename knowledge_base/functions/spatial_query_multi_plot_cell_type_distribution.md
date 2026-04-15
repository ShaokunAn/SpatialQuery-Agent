# plot_cell_type_distribution (spatial_query_multi)

## Description
Visualize the distribution of cell types across datasets using a stacked bar plot. Shows either absolute cell counts or proportions for each dataset/condition.

## Function Signature
```python
sp_multi.plot_cell_type_distribution(
    dataset: Optional[Union[str, List[str]]] = None,
    data_type: Literal['number', 'proportion'] = 'proportion',
    colormap: str = 'tab20c',
    save_path: Optional[str] = None,
)
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| dataset | Optional[Union[str, List[str]]] | None | Dataset(s) to include. If None, uses all |
| data_type | str | 'proportion' | Plot by 'number' (absolute counts) or 'proportion' |
| colormap | str | 'tab20c' | Matplotlib colormap name |
| save_path | Optional[str] | None | Path to save the figure |

## Return Value
A stacked bar plot figure.

## Usage Example
```python
# Define parameters
data_type = 'proportion'

# Plot proportions across all datasets
sp_multi.plot_cell_type_distribution(data_type=data_type)

# Plot cell counts for specific datasets
datasets = ['condition_A', 'condition_B']
data_type_number = 'number'
save_path = 'cell_distribution.png'

sp_multi.plot_cell_type_distribution(
    dataset=datasets,
    data_type=data_type_number,
    save_path=save_path
)
```

## Related Methods
- `plot_cell_type_distribution_fov()`: Per-FOV distribution within a single dataset
