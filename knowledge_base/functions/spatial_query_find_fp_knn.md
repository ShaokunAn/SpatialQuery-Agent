# find_fp_knn (spatial_query)

## Description
Find frequent patterns (cell type co-localization motifs) within the k-nearest neighbors (KNN) of a specified cell type in a single FOV. Uses FP-Growth algorithm on binarized neighbor transactions to identify cell types that frequently co-occur in the neighborhood.

## Function Signature
```python
sp.find_fp_knn(
    ct: str,
    k: int = 30,
    min_support: float = 0.5,
    max_dist: float = 20,
) -> pd.DataFrame
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| ct | str | Required | Cell type name to use as center cells |
| k | int | 30 | Number of nearest neighbors to consider |
| min_support | float | 0.5 | Minimum frequency threshold to consider a pattern as frequent (range: 0 to 1) |
| max_dist | float | 20 | Maximum distance for considering a cell as a valid neighbor |

## Return Value
`pd.DataFrame` with columns:
- `support`: Frequency of the pattern across all neighborhoods of the center cell type
- `itemsets`: Cell type combinations forming the frequent pattern (as frozenset)

## Usage Example
```python
from SpatialQuery import spatial_query

# Initialize spatial_query object
spatial_key = 'X_spatial'
label_key = 'cell_type'
feature_name = 'gene'

sp = spatial_query(
    adata=adata,
    spatial_key=spatial_key,
    label_key=label_key,
    feature_name=feature_name
)

# Define parameters for finding frequent patterns
ct = 'your_anchor_cell_type'
k = 10
min_support = 0.5
max_dist = 20

# Find frequent patterns in the KNN neighborhood
fp = sp.find_fp_knn(
    ct=ct,
    k=k,
    min_support=min_support,
    max_dist=max_dist
)
print(fp)
```

## Notes
- Only maximal frequent patterns are returned (subsets of larger patterns are removed).
- The `max_dist` parameter filters out neighbors that are too far away even within the KNN.
- Raises `ValueError` if the specified cell type is not found in the data.
- `motif_enrichment_knn()` also supports this function and additionally quantify the enrichmnt significance of frequent patterns.

## Troubleshooting
- **No patterns found**: Try lowering `min_support` or increasing `k` / `max_dist`.
- **Too many patterns**: Increase `min_support` or decrease `k` / `max_dist`.

## Related Methods
- `find_fp_dist()`: Distance-based alternative for finding frequent patterns
- `motif_enrichment_knn()`: Statistical enrichment analysis of discovered motifs
