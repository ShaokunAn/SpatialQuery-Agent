# differential_analysis_dist (spatial_query_multi)

## Description
Perform differential analysis of spatial motif patterns between two datasets using radius-based neighborhood.
This function identifies motif patterns that are differentially enriched in the radius-based
neighborhood of a center cell type between two conditions (e.g., disease vs control). It supports
two modes:
1. Unbiased discovery mode (motifs=None): Automatically discovers frequent patterns of each FOV in both 
    datasets, then tests for differential enrichment.
2. Hypothesis-driven mode (motifs specified): Tests user-specified motifs for differential
    enrichment, allowing validation of known or hypothesized spatial patterns.

## Function Signature
```python
sp_multi.differential_analysis_dist(
    ct: str,
    datasets: List[str],
    motifs: Optional[Union[str, List[str], List[List[str]]]] = None,
    max_dist: float = 20,
    min_support: float = 0.5,
    min_size: int = 0,
)
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| ct | str | Required | Center cell type of interest |
| datasets | List[str] | Required | Exactly 2 dataset names to compare |
| motifs | Optional[Union[str, List[str], List[List[str]]]] | None | User-specified motifs to test. If None, auto-discovers first |
| max_dist | float | 20 | Maximum radius for neighbors |
| min_support | float | 0.5 | Minimum frequency for pattern discovery (when motifs=None) |
| min_size | int | 0 | Minimum neighborhood size |

## Return Value
`dict` with keys as dataset names and values as DataFrames:
- Each DataFrame has columns:
  - `itemsets`: The motif pattern (as tuple)
  - `adj-pval`: FDR-corrected p-value
- Only significant patterns (adj-pval < 0.05) for each dataset are included.

## Usage Example
```python
# Define parameters
ct = 'your_anchor_cell_type'
datasets = ['condition_A', 'condition_B']
max_dist = 5
min_support = 0.3

# Auto-discover and compare
diff_result = sp_multi.differential_analysis_dist(
    ct=ct,
    datasets=datasets,
    max_dist=max_dist,
    min_support=min_support
)

# Test specific motifs
motifs = [['cell_type1', 'cell_type2'], ['cell_type3']]

diff_result = sp_multi.differential_analysis_dist(
    ct=ct,
    datasets=datasets,
    motifs=motifs,
    max_dist=max_dist
)

print(diff_result['condition_A'])  # condition_A-specific patterns
print(diff_result['condition_B'])  # condition_B-specific patterns
```

## Notes
- Requires exactly 2 datasets for comparison.
- **Motif input formats** (when specified):
  - Single cell type: `'CellTypeA'`
  - Single motif: `['CellTypeA', 'CellTypeB']`
  - Multiple motifs: `[['CellTypeA'], ['CellTypeB', 'CellTypeC']]`
- When `motifs=None`, the method discovers maximal frequent patterns above the frequency threshold (`min_support`) in each FOV of both datasets, collects the union of all discovered motifs as candidates, then tests each for differential enrichment. Due to the frequency threshold and maximal pattern constraints, this automatic discovery does not return all differential patterns but serves as a reference. Users can identify cell types appearing in the differential patterns and define their own patterns of interest to further validate whether they are differentially enriched.
- Statistical testing uses hypergeometric test with FDR correction within each condition.

## Related Methods
- `differential_analysis_knn()`: KNN-based alternative
- `motif_enrichment_dist()`: Enrichment without differential comparison
