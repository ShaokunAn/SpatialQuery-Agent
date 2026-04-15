# motif_enrichment_knn (spatial_query_multi)

## Description
Perform motif enrichment analysis using k-nearest neighbors across multiple FOVs. Aggregates statistics from all specified FOVs and performs a single hypergeometric test per motif, with FDR correction for multiple testing.

## Function Signature
```python
sp_multi.motif_enrichment_knn(
    ct: str,
    motifs: Union[str, List[str]] = None,
    dataset: Union[str, List[str]] = None,
    k: int = 30,
    min_support: float = 0.5,
    max_dist: float = 20,
    return_cellID: bool = False,
) -> Union[pd.DataFrame, Tuple[pd.DataFrame, dict, dict]]
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| ct | str | Required | Center cell type |
| motifs | Union[str, List[str]] | None | Motifs to test. If None, auto-discovers frequent patterns first |
| dataset | Union[str, List[str]] | None | Dataset name(s) to include. If None, uses all datasets |
| k | int | 30 | Number of nearest neighbors |
| min_support | float | 0.5 | Minimum frequency for pattern discovery (only when motifs=None) |
| max_dist | float | 20 | Maximum distance for valid neighbors |
| return_cellID | bool | False | If True, return cell indices per FOV |

## Return Value
If `return_cellID=False`:
- `pd.DataFrame` with columns: `center`, `motifs`, `n_center_motif`, `n_center`, `n_motif`, `expectation`, `p-values`, `adj-pval`, `if_significant`

If `return_cellID=True`:
- Tuple of `(pd.DataFrame, motif_cell_ids, center_cell_ids)` where:
  - `motif_cell_ids`: `{motif_str: {dataset_i: [cell_indices]}}` - motif cell IDs per FOV
  - `center_cell_ids`: `{motif_str: {dataset_i: [cell_indices]}}` - center cell IDs per FOV

## Usage Example
```python
# Define parameters
ct = 'your_anchor_cell_type'
k = 30
min_support = 0.3

# Auto-discover and test motifs across all FOVs
result = sp_multi.motif_enrichment_knn(
    ct=ct,
    k=k,
    min_support=min_support
)

# Test specific motifs in a specific condition only
motifs = ['cell_type1', 'cell_type2']
dataset = 'condition_A'

result = sp_multi.motif_enrichment_knn(
    ct=ct,
    motifs=motifs,
    dataset=dataset,
    k=k
)

# Get cell IDs for downstream DE analysis
motifs_single = ['cell_type1']

result, motif_ids, center_ids = sp_multi.motif_enrichment_knn(
    ct=ct,
    motifs=motifs_single,
    return_cellID=True
)
```

## Notes
- Statistics are aggregated across FOVs: total cells, center cells, and motif occurrences are summed.
- Cell IDs are organized by FOV (using modified dataset names like 'disease_0', 'disease_1').
- The cell IDs from `return_cellID=True` can be directly used with `de_genes()` for multi-FOV DE analysis.
- FDR correction uses positive correlation method.

## Related Methods
- `motif_enrichment_dist()`: Radius-based alternative
- `de_genes()`: Use cell IDs from this method for DE analysis
- `find_fp_knn()`: Discover frequent patterns first
