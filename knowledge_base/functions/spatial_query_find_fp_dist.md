# find_fp_dist (spatial_query)

## Description
Find frequent patterns (cell type co-localization motifs) within a radius-based neighborhood of a specified/anchor cell type in a single FOV. Uses FP-Growth algorithm on binarized neighbor transactions to identify cell types that frequently co-occur within a given distance.

## Function Signature
```python
sp.find_fp_dist(
    ct: str,
    max_dist: float = 20,
    min_size: int = 0,
    min_support: float = 0.5,
) -> pd.DataFrame
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| ct | str | Required | Cell type name to use as center cells |
| max_dist | float | 20 | Maximum radius distance for considering a cell as a neighbor |
| min_size | int | 0 | Minimum neighborhood size (number of neighbors) for each center cell to be included |
| min_support | float | 0.5 | Minimum frequency threshold to consider a pattern as frequent (range: 0 to 1) |

## Return Value
`pd.DataFrame` with columns:
- `support`: Frequency of the pattern across all neighborhoods of the center cell type
- `itemsets`: Cell type combinations forming the frequent pattern (as frozenset)

## Usage Example
```python
from SpatialQuery import spatial_query

# Custom parameters based on your data
spatial_key = 'X_spatial'
label_key = 'cell_type'
feature_name = 'gene'

sp = spatial_query(
    adata=adata, 
    spatial_key=spatial_key, 
    label_key=label_key, 
    feature_name=feature_name
    )

# Find frequent patterns surrounding a specific cell type
# Custom neighborhood size (max_dist) and frequency threshold (min_support)
ct = 'your_anchor_cell_type'
max_dist = 5
min_support = 0.5
min_size = 0

fp = sp.find_fp_dist(
    ct=ct,
    max_dist=max_dist,
    min_support=min_support,
    min_size=min_size
)
print(fp)
```

## Notes
- Only maximal frequent patterns are returned (subsets of larger patterns are removed).
- The `min_size` parameter filters out center cells with very few neighbors, which can reduce noise.
- Raises `ValueError` if the specified cell type is not found in the data.
- `motif_enrichment_dist()` also supports this function and additionally quantify the enrichmnt significance of frequent patterns.

## Troubleshooting
- **No patterns found**: Try lowering `min_support` or increasing `max_dist`.
- **Too many patterns**: Increase `min_support` or decrease `max_dist`.

## Related Methods
- `find_fp_knn()`: KNN-based alternative for finding frequent patterns/motifs
- `motif_enrichment_dist()`: Statistical enrichment analysis of discovered motifs
