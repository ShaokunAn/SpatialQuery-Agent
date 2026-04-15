# differential_analysis_knn (spatial_query_multi)

## Description
Perform differential analysis of spatial motif patterns between two datasets using KNN neighborhood.
This function identifies motif patterns that are differentially enriched in the KNN neighborhood
of a center cell type between two conditions (e.g., disease vs control). It supports two modes:
1. Unbiased discovery mode (motifs=None): Automatically discovers frequent patterns of each FOV in both 
  datasets, then tests for differential enrichment.
2. Hypothesis-driven mode (motifs specified): Tests user-specified motifs for differential
  enrichment, allowing validation of known or hypothesized spatial patterns.

## Function Signature
```python
sp_multi.differential_analysis_knn(
    ct: str,
    datasets: List[str],
    motifs: Optional[Union[str, List[str], List[List[str]]]] = None,
    k: int = 30,
    min_support: float = 0.5,
    max_dist: float = 20,
)
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| ct | str | Required | Center cell type of interest |
| datasets | List[str] | Required | Exactly 2 dataset names to compare |
| motifs | Optional[Union[str, List[str], List[List[str]]]] | None | User-specified motifs to test (see Notes for formats). If None, auto-discovers patterns first |
| k | int | 30 | Number of nearest neighbors |
| min_support | float | 0.5 | Minimum frequency for pattern discovery (only used when motifs=None) |
| max_dist | float | 20 | Maximum distance for valid neighbors |

## Return Value
`dict` with keys as dataset names and values as DataFrames containing significantly enriched patterns:
- Each DataFrame has columns:
  - `itemsets`: The motif pattern (as tuple)
  - `adj-pval`: FDR-corrected p-value
- Only significant patterns (adj-pval < 0.05) for each dataset are included.

## Usage Example
```python
# Define parameters
ct = 'your_anchor_cell_type'
datasets = ['condition_A', 'condition_B']
k = 30
min_support = 0.3

# Auto-discover and compare patterns between conditions
diff_result = sp_multi.differential_analysis_knn(
    ct=ct,
    datasets=datasets,
    k=k,
    min_support=min_support
)

# Check condition_A-specific patterns
print(diff_result['condition_A'])

# Check condition_B-specific patterns
print(diff_result['condition_B'])

# Test specific motifs
motifs = [['cell_type1', 'cell_type2'], ['cell_type3']]

diff_result = sp_multi.differential_analysis_knn(
    ct=ct,
    datasets=datasets,
    motifs=motifs,
    k=k
)
```

## Notes
- Requires exactly 2 datasets for comparison.
- **Motif input formats** (when specified):
  - Single cell type: `'CellTypeA'`
  - Single motif: `['CellTypeA', 'CellTypeB']`
  - Multiple motifs: `[['CellTypeA'], ['CellTypeB', 'CellTypeC']]`
- When `motifs=None`, the method discovers maximal frequent patterns above the frequency threshold (`min_support`) in each FOV of both datasets, collects the union of all discovered motifs as candidates, then tests each for differential enrichment. Due to the frequency threshold and maximal pattern constraints, this automatic discovery does not return all differential patterns but serves as a reference. Users can identify cell types appearing in the differential patterns and define their own patterns of interest to further validate whether they are differentially enriched.
- Statistical testing uses hypergeometric test with FDR correction within each condition.

## Troubleshooting
- **No differential patterns**: Lower `min_support` to discover more patterns, or adjust `k`/`max_dist`.
- **ValueError about datasets**: Ensure the dataset names match those used during initialization.

## Related Methods
- `differential_analysis_dist()`: Radius-based alternative
- `motif_enrichment_knn()`: Enrichment analysis without differential comparison
