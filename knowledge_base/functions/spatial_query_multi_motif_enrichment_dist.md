# motif_enrichment_dist (spatial_query_multi)

## Description
Perform motif enrichment analysis within radius-based neighborhoods across multiple FOVs. Aggregates statistics from all FOVs of specified datasets and performs a single hypergeometric test per motif, with FDR correction.

## Function Signature
```python
sp_multi.motif_enrichment_dist(
    ct: str,
    motifs: Union[str, List[str], List[List[str]]] = None,
    dataset: Union[str, List[str]] = None,
    max_dist: float = 20,
    min_size: int = 0,
    min_support: float = 0.5,
    return_cellID: bool = False,
) -> Union[pd.DataFrame, Tuple[pd.DataFrame, dict, dict]]
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| ct | str | Required | Center cell type |
| motifs | Union[str, List[str], List[List[str]]] | None | Motifs to test. Accepts flexible formats (see Notes). If None, auto-discovers |
| dataset | Union[str, List[str]] | None | Dataset name(s) to include. If None, uses all |
| max_dist | float | 20 | Maximum radius for neighbors |
| min_size | int | 0 | Minimum neighborhood size |
| min_support | float | 0.5 | Minimum frequency for pattern discovery (when motifs=None) |
| return_cellID | bool | False | If True, return cell indices per FOV |

## Return Value
If `return_cellID=False`:
- `pd.DataFrame` with columns: `center`, `motifs`, `n_center_motif`, `n_center`, `n_motif`, `expectation`, `p-values`, `adj-pval`, `if_significant`

If `return_cellID=True`:
- Tuple of `(pd.DataFrame, motif_cell_ids, center_cell_ids)` where:
  - `motif_cell_ids`: `{motif_str: {dataset_i: [cell_indices]}}` per FOV
  - `center_cell_ids`: `{motif_str: {dataset_i: [cell_indices]}}` per FOV

## Usage Example
```python
# Define parameters
ct = 'your_anchor_cell_type'
max_dist = 10
min_support = 0.3

# Auto-discover motifs
result = sp_multi.motif_enrichment_dist(
    ct=ct,
    max_dist=max_dist,
    min_support=min_support
)

# Test a single motif (single cell type)
motifs_single = 'cell_type1'

result = sp_multi.motif_enrichment_dist(
    ct=ct,
    motifs=motifs_single,
    max_dist=max_dist
)

# Test a single motif (multiple cell types together)
motifs_multi = ['cell_type1', 'cell_type2']

result = sp_multi.motif_enrichment_dist(
    ct=ct,
    motifs=motifs_multi,
    max_dist=max_dist
)

# Test multiple motifs
motifs_list = [['cell_type1'], ['cell_type2', 'cell_type3']]

result = sp_multi.motif_enrichment_dist(
    ct=ct,
    motifs=motifs_list,
    max_dist=max_dist
)

# Get cell IDs for downstream analysis
motifs_for_ids = ['cell_type1']

result, motif_ids, center_ids = sp_multi.motif_enrichment_dist(
    ct=ct,
    motifs=motifs_for_ids,
    return_cellID=True
)
```

## Notes
- **Flexible motif input formats**:
  - Single cell type string: `'cell_type1'` 
  - Single motif list: `['cell_type1', 'cell_type2']`
  - Multiple motifs: `[['cell_type1'], ['cell_type2', 'cell_type3']]`
- Statistics are aggregated across all specified FOVs.

## Related Methods
- `motif_enrichment_knn()`: KNN-based alternative
- `de_genes()`: Use cell IDs for DE analysis
- `find_fp_dist()`: Discover frequent patterns
