# spatial_query_multi (Class Initialization)

## Description
Initialize a `spatial_query_multi` object for spatial analysis across multiple fields of view (FOVs), supporting FOVs from multiple conditions/groups. This class creates individual `spatial_query` objects for each FOV and provides methods that aggregate analysis results across FOVs for specified groups.

## Class
`spatial_query_multi` (from `SpatialQuery`)

## Constructor Signature
```python
spatial_query_multi(
    adatas: List[AnnData],
    datasets: List[str],
    spatial_key: str,
    label_key: str,
    leaf_size: int = 10,
    build_gene_index: bool = False,
    feature_name: str = None,
    if_lognorm: bool = True,
)
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| adatas | List[AnnData] | Required | List of AnnData objects, one per FOV |
| datasets | List[str] | Required | List of dataset names corresponding to each AnnData. FOVs from the same group should share the same dataset name |
| spatial_key | str | Required | Key in `adata.obsm` for spatial coordinates |
| label_key | str | Required | Key in `adata.obs` for cell type labels |
| leaf_size | int | 10 | Leaf size for KDTree construction |
| build_gene_index | bool | False | If True, build scfind index to save memory; if False, use adata.X directly |
| feature_name | str | None | Key in `adata.var` for feature/gene names |
| if_lognorm | bool | True | Whether to log-normalize expression data (only when build_gene_index=False) |

## Return Value
A `spatial_query_multi` object with key attributes:
- `spatial_queries`: List of `spatial_query` objects, one per FOV
- `datasets`: List of modified dataset names (format: `datasetname_fovindex`)
- `spatial_key`, `label_key`, `build_gene_index`: Stored configuration

## Usage Example
```python
from SpatialQuery import spatial_query_multi

# Multiple FOVs from multiple conditions
adatas = [adata_group_A_fov1, adata_group_A_fov2, adata_group_B_fov1, adata_group_B_fov2]
datasets = ['group_A', 'group_A', 'group_B', 'group_B']
spatial_key = 'X_spatial'
label_key = 'cell_type'
feature_name = 'gene'
build_gene_index = False # will use adata.X directly or build scfind index
if_lognorm = True # will perform log-normalization for adata.X. set to False if data is already log-normalized

sp_multi = spatial_query_multi(
    adatas=adatas,
    datasets=datasets,
    spatial_key=spatial_key,
    label_key=label_key,
    feature_name=feature_name,
    build_gene_index=build_gene_index,
    if_lognorm=if_lognorm
)

# Multiple FOVs from a single group
adatas = [adata_fov1, adata_fov2, adata_fov3]
datasets = ['sample', 'sample', 'sample']
build_gene_index = False

sp_multi = spatial_query_multi(
    adatas=adatas,
    datasets=datasets,
    spatial_key=spatial_key,
    label_key=label_key,
    feature_name=feature_name,
    build_gene_index=build_gene_index,
    if_lognorm=if_lognorm
)
```

## Notes
- Dataset names are automatically modified to include FOV indices (e.g., 'disease_0', 'disease_1', 'control_0', 'control_1').
- Underscores in dataset names are replaced with hyphens to maintain internal naming conventions.
- Each FOV gets its own KD-tree and optional gene index.

## Related Methods
After initialization, you can use:
- `find_fp_knn()` / `find_fp_dist()`: Find frequent patterns across FOVs of specified group
- `motif_enrichment_knn()` / `motif_enrichment_dist()`: Enrichment analysis across FOVs of specified group
- `differential_analysis_knn()` / `differential_analysis_dist()`: Compare patterns between groups. Two groups are required for differential analysis.
- `de_genes()`: Differential expression between given populations
- `compute_gene_gene_correlation()` / `compute_gene_gene_correlation_by_type()`: Gene co-variation analysis
- `test_score_difference()`: Compare scores between covarying gene-pairs
- `plot_cell_type_distribution()`: Visualize the distribution of cell types across datasets using a stacked bar plot
- `plot_cell_type_distribution_fov()`: Visualize the distribution of cell types across FOVs in the dataset using a stacked bar plot
