import anndata
from typing import Optional

_SPATIAL_KEY_CANDIDATES  = ["X_spatial", "spatial", "X_umap"]
_LABEL_KEY_CANDIDATES    = ["cell_type", "predicted_label", "celltype", "cell_label"]
_FEATURE_NAME_CANDIDATES = ["gene", "gene_name", "feature_name", "gene_ids"]


def load_adata(path: str) -> anndata.AnnData:
    return anndata.read_h5ad(path)


def infer_spatial_key(adata) -> Optional[str]:
    return next((k for k in _SPATIAL_KEY_CANDIDATES if k in adata.obsm), None)


def infer_label_key(adata) -> Optional[str]:
    return next((k for k in _LABEL_KEY_CANDIDATES if k in adata.obs.columns), None)


def infer_feature_name(adata) -> Optional[str]:
    return next((k for k in _FEATURE_NAME_CANDIDATES if k in adata.var.columns), None)


def inspect_adata(adata) -> dict:
    spatial_key  = infer_spatial_key(adata)
    label_key    = infer_label_key(adata)
    feature_name = infer_feature_name(adata)
    cell_types   = sorted(adata.obs[label_key].dropna().unique().tolist()) if label_key else None
    return {
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "obsm_keys": list(adata.obsm.keys()),
        "obs_columns": list(adata.obs.columns),
        "var_columns": list(adata.var.columns),
        "inferred_spatial_key": spatial_key,
        "inferred_label_key": label_key,
        "inferred_feature_name": feature_name,
        "cell_types": cell_types,
    }


def format_params_for_confirmation(summary: dict) -> str:
    """Return a human-readable explanation of inferred parameters with available options."""
    lines = [
        "**Inferred parameters** (please verify):\n",
        f"- `spatial_key = {summary['inferred_spatial_key']!r}`",
        f"  → key in `adata.obsm` for spatial coordinates",
        f"  → available: `{summary['obsm_keys']}`\n",
        f"- `label_key = {summary['inferred_label_key']!r}`",
        f"  → column in `adata.obs` for cell type labels",
        f"  → available: `{summary['obs_columns']}`\n",
        f"- `feature_name = {summary['inferred_feature_name']!r}`",
        f"  → column in `adata.var` for gene names (set `None` if using index)",
        f"  → available: `{summary['var_columns']}`",
    ]
    return "\n".join(lines)


def update_summary_params(summary: dict, spatial_key: str, label_key: str,
                          feature_name: Optional[str], adata) -> dict:
    """Return a new summary dict with updated parameters and refreshed cell_types."""
    updated = dict(summary)
    updated["inferred_spatial_key"] = spatial_key
    updated["inferred_label_key"] = label_key
    updated["inferred_feature_name"] = feature_name if feature_name and feature_name.lower() != "none" else None
    updated["cell_types"] = sorted(adata.obs[label_key].dropna().unique().tolist())
    return updated


def format_summary_for_llm(summary: dict, filename: str = "") -> str:
    lines = []
    if filename:
        lines.append(f"Dataset: {filename}")
    lines.append(f"Shape: {summary['n_cells']} cells x {summary['n_genes']} genes")
    lines.append(f"adata.obsm keys: {summary['obsm_keys']}")
    lines.append(f"adata.obs columns: {summary['obs_columns']}")
    lines.append(f"adata.var columns: {summary['var_columns']}")
    lines.append("Inferred parameters:")
    lines.append(f"  spatial_key  = {summary['inferred_spatial_key']!r}")
    lines.append(f"  label_key    = {summary['inferred_label_key']!r}")
    lines.append(f"  feature_name = {summary['inferred_feature_name']!r}")
    if summary["cell_types"] is not None:
        lines.append(f"  Cell types ({len(summary['cell_types'])}): {', '.join(summary['cell_types'])}")
    return "\n".join(lines)
