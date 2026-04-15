# get_anchor_motif_cell_ids (spatial_query)

## Description
Get cell grouping information (anchor/center cells, neighbor motif cells, non-neighbor motif cells) without computing gene-gene correlations. This is a lightweight alternative to `compute_gene_gene_correlation` when you only need the cell IDs for downstream analysis (e.g., DE analysis, custom visualization, or passing to external tools).

Returns the same `cell_groups` dict that `compute_gene_gene_correlation` returns as its second output, but skips the expensive correlation computation.

## Function Signature
```python
sp.get_anchor_motif_cell_ids(
    ct: str,
    motif: Union[str, List[str]],
    max_dist: Optional[float] = None,
    k: Optional[int] = None,
    min_size: int = 0,
) -> dict
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| ct | str | Required | Cell type of the center/anchor cells |
| motif | Union[str, List[str]] | Required | Motif cell type(s) to find in neighborhoods |
| max_dist | Optional[float] | None | Maximum distance for radius-based neighbor search. Use either max_dist or k |
| k | Optional[int] | None | Number of nearest neighbors. Use either max_dist or k |
| min_size | int | 0 | Minimum neighborhood size (only used with max_dist) |

## Return Value
**cell_groups** (`dict`) with the same structure as the second return value of `compute_gene_gene_correlation`:
- `center_neighbor_motif_pair`: Array of shape (n_pairs, 2) with [center_idx, neighbor_idx] pairs. Each row is [center_cell_idx, neighbor_cell_idx].
- `non-neighbor_motif_cells`: Array of distant motif cell indices.
- `non_motif_center_neighbor_pair`: Array of shape (n_pairs, 2) for centers without motif. Each row is [center_cell_idx, neighbor_cell_idx].

Individual cell IDs can be extracted from pairs using np.unique() like:
- center_cells = np.unique(center_neighbor_motif_pair[:, 0])
- neighbor_cells = np.unique(center_neighbor_motif_pair[:, 1])

## Usage Example
```python
from SpatialQuery import spatial_query
import numpy as np

sp = spatial_query(
    adata=adata,
    spatial_key='X_spatial',
    label_key='cell_type',
    feature_name='gene',
    build_gene_index=False,
    if_lognorm=True,
)

# Get cell IDs for anchor-motif spatial relationships
ids = sp.get_anchor_motif_cell_ids(
    ct='TypeA',
    motif=['TypeB', 'TypeC'],
    max_dist=20,
)

# Visualize spatial distribution of center/anchor-motif cell groups
sp.plot_all_center_motif(ct='TypeA', ids=ids)

# Extract individual cell indices if needed
pairs = ids['center_neighbor_motif_pair']
center_ids = np.unique(pairs[:, 0])
neighbor_ids = np.unique(pairs[:, 1])
non_neighbor_ids = ids['non-neighbor_motif_cells']
```

## Notes
- Must specify either `max_dist` or `k`, but not both.
- Use this instead of `compute_gene_gene_correlation` when you only need cell IDs (e.g., for DE analysis, plotting, or custom downstream analysis). It avoids the expensive gene-gene correlation computation.
- The returned dict is fully compatible with plotting functions like `plot_all_center_motif` and `plot_gene_pair_spatial`.

## Related Methods
- `compute_gene_gene_correlation()`: Full correlation analysis (also returns cell IDs as second output)
- `plot_all_center_motif()`: Visualize center/anchor-motif cell groups using the returned ids
- `de_genes()`: Differential expression using extracted cell indices
