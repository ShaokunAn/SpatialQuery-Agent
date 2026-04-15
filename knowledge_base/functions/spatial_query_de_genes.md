# de_genes (spatial_query)

## Description
Identify differentially expressed genes between two groups of cells within a single FOV. Supports multiple statistical methods depending on how the data was indexed during initialization. This is typically used after motif enrichment analysis to compare gene expression between cells in different spatial contexts.

## Function Signature
```python
sp.de_genes(
    ind_group1: List[int],
    ind_group2: List[int],
    genes: Optional[Union[str, List[str]]] = None,
    min_fraction: float = 0.05,
    method: Literal['fisher', 't-test', 'wilcoxon'] = 'fisher',
    alpha: Optional[float] = None,
) -> pd.DataFrame
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| ind_group1 | List[int] | Required | List of cell indices for group 1 |
| ind_group2 | List[int] | Required | List of cell indices for group 2 |
| genes | Optional[Union[str, List[str]]] | None | Gene names to test. If None, tests all genes |
| min_fraction | float | 0.05 | Minimum fraction of cells expressing a gene for it to be tested |
| method | str | 'fisher' | Statistical test method: 'fisher', 't-test', or 'wilcoxon' |
| alpha | Optional[float] | None | Significance threshold for adjusted p-values. Defaults to 0.1 for Fisher's test, 0.05 otherwise |

## Return Value
`pd.DataFrame` containing differentially expressed genes with columns:
- `gene`: Gene name
- `proportion_1`: Proportion of cells expressing the gene in group 1
- `proportion_2`: Proportion of cells expressing the gene in group 2
- `p_value`: Raw p-value
- `adj-pval`: FDR-corrected p-value (BH method)
- `de_in`: Which group shows higher expression ('group1' or 'group2')

Additional columns when `build_gene_index=True`:
- `a`, `b`, `c`, `d`: Contingency table values

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

# Define two groups of cells by their indices
ind_group1 = [0, 1, 2, 3, 4, 5]
ind_group2 = [10, 11, 12, 13, 14, 15]

# Perform DE analysis with Wilcoxon test
method = 'wilcoxon'
alpha = 0.05

de_result = sp.de_genes(
    ind_group1=ind_group1,
    ind_group2=ind_group2,
    method=method,
    alpha=alpha
)

# Test specific genes with Fisher's exact test
genes = ['gene1', 'gene2', 'gene3']
method = 'fisher'

de_result = sp.de_genes(
    ind_group1=ind_group1,
    ind_group2=ind_group2,
    genes=genes,
    method=method,
    alpha=alpha
)
```

## Notes
- **Recommended**: Use `motif_enrichment_knn()` or `motif_enrichment_dist()` with `return_cellID=True` to obtain cell indices for spatially-defined groups. See the tutorial for a complete example of using `de_genes()` to perform motif-associated DE analysis.
- When `build_gene_index=True` during initialization, only Fisher's exact test is supported. Other methods will trigger a warning and default to Fisher's. `alpha = 0.1` is recommended.
- When `build_gene_index=False`, all three methods (fisher, t-test, wilcoxon) are available using the stored adata.X matrix. `alpha = 0.05` is recommended.
- Results are filtered by the `alpha` significance threshold before returning.
- The `min_fraction` filter removes genes expressed in too few cells to ensure reliable statistical testing.

## Troubleshooting
- **No DE genes found**: Try increasing `alpha`, decreasing `min_fraction`, or ensuring group sizes are sufficient.
- **Method warning**: If you see a warning about method being ignored, it means `build_gene_index=True` was set during initialization.

## Related Methods
- `motif_enrichment_knn()` / `motif_enrichment_dist()`: Get cell IDs for DE analysis (use `return_cellID=True`)
