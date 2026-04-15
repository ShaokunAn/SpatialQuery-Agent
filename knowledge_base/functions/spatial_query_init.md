# spatial_query (Class Initialization)

## Description
Initialize a `spatial_query` object for spatial analysis of a single field of view (FOV). During initialization, the following preprocessing steps are performed:

1. **Spatial coordinate normalization**: Spatial coordinates are normalized so that the mean nearest neighbor distance equals 1. This means 1 unit of distance is approximately the length of one cell. The purpose is to unify the distance scale across different spatial technologies (e.g., MERFISH, Visium, CODEX) so that distance-based parameters are comparable regardless of the original resolution or coordinate units.
2. **KD-tree construction**: A KD-tree is built on the normalized 2D spatial coordinates for fast neighbor search (both KNN and radius-based queries).
3. **Gene expression storage** (two options):
   - `build_gene_index=True`: Uses scfind to build a compressed binary index of gene expression data, which is memory-efficient for large datasets. Only supports Fisher's exact test for downstream differential expression analysis.
   - `build_gene_index=False` (default): Stores the original AnnData expression matrix directly. If `if_lognorm=True`, the data is log-normalized following the standard scanpy pipeline (`scanpy.pp.normalize_total` + `scanpy.pp.log1p`). This mode supports Fisher's exact test, t-test, and Wilcoxon test for downstream analysis.

## Class
`spatial_query` (from `SpatialQuery`)

## Constructor Signature
```python
spatial_query(
    adata: AnnData,
    dataset: str = 'ST',
    spatial_key: str = 'X_spatial',
    label_key: str = 'predicted_label',
    leaf_size: int = 10,
    build_gene_index: bool = False,
    feature_name: str = None,
    if_lognorm: bool = True,
    if_normalize_spatial_coord: bool = True,
)
```

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| adata | AnnData | Required | AnnData object containing spatial coordinates and cell type annotations |
| dataset | str | 'ST' | Dataset name identifier for this single FOV |
| spatial_key | str | 'X_spatial' | Key in `adata.obsm` for spatial coordinates |
| label_key | str | 'predicted_label' | Key in `adata.obs` for cell type labels |
| leaf_size | int | 10 | Leaf size for KDTree construction |
| build_gene_index | bool | False | If True, build scfind index for gene expression queries; if False, use adata.X directly |
| feature_name | str | None | Key in `adata.var` for feature/gene names. Required for gene expression analysis |
| if_lognorm | bool | True | Whether to log-normalize expression data (only when build_gene_index=False). Set to False if data is already log-normalized |
| if_normalize_spatial_coord | bool | True | Whether to normalize spatial coordinates so that 1 unit ≈ 1 cell diameter (mean nearest-neighbor distance = 1). Set to False to preserve original coordinate units. **Affects how distance parameters are interpreted**: with normalization on, `max_dist=20` means ~20 cell diameters; with normalization off, `max_dist` is in the original coordinate units of the data |

## Return Value
A `spatial_query` object with the following key attributes:
- `kd_tree`: KDTree built from spatial coordinates
- `labels`: Cell type labels (categorical)
- `spatial_pos`: Normalized spatial positions
- `genes`: List of gene names
- `adata`: AnnData object (when build_gene_index=False)
- `index`: SCFind index (when build_gene_index=True)

## Usage Example
```python
from SpatialQuery import spatial_query

# Define parameters
dataset = 'ST'
spatial_key = 'X_spatial'
label_key = 'cell_type'
leaf_size = 10
feature_name = 'gene'
if_lognorm = True  # will perform log-normalization. set to False if data is already log-normalized

# Basic initialization with raw data
sp = spatial_query(
    adata=adata,
    dataset=dataset,
    spatial_key=spatial_key,
    label_key=label_key,
    leaf_size=leaf_size,
    feature_name=feature_name,
    if_lognorm=if_lognorm,
    if_normalize_spatial_coord=True,   # default: normalize so 1 unit ≈ 1 cell diameter
)

# With scfind gene index to handle large data
build_gene_index = True

sp = spatial_query(
    adata=adata,
    spatial_key=spatial_key,
    label_key=label_key,
    leaf_size=leaf_size,
    build_gene_index=build_gene_index,
    feature_name=feature_name
)
```

## Notes
- **`if_normalize_spatial_coord`**: When `True` (default), spatial coordinates are normalized at init time so that 1 unit ≈ 1 cell diameter. This makes `max_dist` interpretable as "number of cell diameters" regardless of the original technology's coordinate units (e.g., microns for MERFISH, pixel units for Visium). Set to `False` only if you need to work in original coordinate units — in that case `max_dist` values must be specified in those units and will not be cross-technology comparable.
- **Spatial normalization**: The normalization rescales coordinates so the mean nearest-neighbor distance equals 1.
- Features with NA values in `adata.var[feature_name]` are filtered out.
- Duplicated features are removed (keeping the first occurrence).
- When `build_gene_index=False` and `if_lognorm=True`, log-normalization is performed using scanpy's standard pipeline (`normalize_total` + `log1p`). Set `if_lognorm=False` if data is already log-normalized to avoid double normalization.
- When `build_gene_index=True`, scfind compresses expression data into a binary index, significantly reducing memory usage. This is recommended for very large datasets but limits DE analysis to Fisher's exact test only.

## Related Methods
After initialization, you can use the following methods:
- `find_fp_knn()` / `find_fp_dist()`: Find frequent patterns around anchor cells
- `motif_enrichment_knn()` / `motif_enrichment_dist()`: Motif enrichment analysis
- `find_patterns_grid()` / `find_patterns_rand()`: Grid/random-based pattern finding frequent pattern across the FOV
- `de_genes()`: Differential expression analysis
- `compute_gene_gene_correlation()` / `compute_gene_gene_correlation_by_type()`: Gene-gene co-variation analysis
- `plot_*()`: Plotting functions
