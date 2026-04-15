# de_genes (spatial_query_multi)

## Description
Perform differential expression analysis between two groups of cells across multiple FOVs. The cell groups are specified as dictionaries mapping modified dataset names to lists of cell indices, enabling flexible cross-FOV and cross-condition comparisons.

## Function Signature
```python
sp_multi.de_genes(
    ind_group1: Dict[str, List[int]],
    ind_group2: Dict[str, List[int]],
    genes: Optional[Union[str, List[str]]] = None,
    min_fraction: float = 0.05,
    method: Literal['fisher', 't-test', 'wilcoxon'] = 'fisher',
    alpha: Optional[float] = None,
) -> pd.DataFrame
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| ind_group1 | Dict[str, List[int]] | Required | Dictionary mapping modified dataset names to cell indices for group 1 |
| ind_group2 | Dict[str, List[int]] | Required | Dictionary mapping modified dataset names to cell indices for group 2 |
| genes | Optional[Union[str, List[str]]] | None | Genes to test. If None, tests all overlapping genes |
| min_fraction | float | 0.05 | Minimum fraction of cells expressing a gene for it to be tested |
| method | str | 'fisher' | Statistical method: 'fisher', 't-test', or 'wilcoxon' |
| alpha | Optional[float] | None | Significance threshold. Defaults to 0.1 for Fisher's, 0.05 otherwise |

## Return Value
`pd.DataFrame` with columns:
- `gene`: Gene name
- `proportion_1`: Expression proportion in group 1
- `proportion_2`: Expression proportion in group 2
- `abs`: Absolute difference in proportions
- `difference`: Proportion difference (group1 - group2)
- `p_value`: Raw p-value
- `adj-pval`: FDR-corrected p-value
- `de_in`: Which group shows higher expression ('group1' or 'group2')

## Usage Example
```python
# Define two groups of cells by their indices per FOV
# Keys must be modified dataset names (e.g., 'condition_A_0', not 'condition_A')
ind_group1 = {'condition_A_0': [0, 1, 2, 3], 'condition_A_1': [0, 1, 2]}
ind_group2 = {'condition_A_0': [10, 11, 12, 13], 'condition_A_1': [10, 11, 12]}

# Perform DE analysis with Wilcoxon test
method = 'wilcoxon'

de_result = sp_multi.de_genes(
    ind_group1=ind_group1,
    ind_group2=ind_group2,
    method=method
)

# Using Fisher's exact test with scfind index
method_fisher = 'fisher'

de_result = sp_multi.de_genes(
    ind_group1=ind_group1,
    ind_group2=ind_group2,
    method=method_fisher
)
```

## Notes
- **Recommended**: Use `motif_enrichment_knn()` or `motif_enrichment_dist()` with `return_cellID=True` to obtain cell indices for spatially-defined groups. See the tutorial for a complete example of using `de_genes()` with motif enrichment results.
- **Dictionary keys must be modified dataset names** (e.g., 'condition_A_0', not 'condition_A'). These are the names stored in `sp_multi.datasets`.
- When `build_gene_index=True`: Only Fisher's exact test is supported. Uses scfind index for efficient binary expression queries and `alpha = 0.1` is recommended.
- When `build_gene_index=False`: All methods are available. Cells are collected from each FOV's adata, concatenated, and tested together. Same analysis as in the pipeline of Scanpy. `alpha = 0.05` is recommended.

## Troubleshooting
- **ValueError about datasets**: Ensure keys in ind_group1/ind_group2 match `sp_multi.datasets`.
- **No DE genes**: Try increasing `alpha`, decreasing `min_fraction`, or ensuring sufficient cells in both groups.

## Related Methods
- `motif_enrichment_knn()` / `motif_enrichment_dist()`: Get cell IDs (use `return_cellID=True`)
- `compute_gene_gene_correlation()` / `compute_gene_gene_correlation_by_type()`: Gene co-variation analysis between anchor cells and motif cells.
