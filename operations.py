# operations.py
"""Operation registry — declarative definitions for all supported analysis operations."""
from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class Operation:
    name: str
    description: str
    rag_query: Optional[str]
    required_ctx: List[str]
    produces: List[str]
    next_suggestions: List[str]
    tier: int  # 1 = template assembly, 3 = RAG free gen


# ---------------------------------------------------------------------------
# Tier 1 structured templates — assembled by code_gen.assemble_code()
# Each is a dict of named code blocks. code_gen selects and fills them.
# ---------------------------------------------------------------------------

MOTIF_TEMPLATE = {
    'postprocess': (
        "significant_motifs = motif_result[motif_result['if_significant']]\n"
        "if len(significant_motifs) > 0:\n"
        "    sp.plot_motif_enrichment_heatmap(significant_motifs)\n"
    ),
}

DE_TEMPLATE = {
    'body': (
        "de_results = {{}}\n"
        "for idx, row in significant_motifs.{motif_slice}.iterrows():\n"
        "    motif_key = str(sorted(list(row['motifs'])))\n"
        "    center_ids = np.array(row['center_id'])\n"
        "    non_center_ids = np.setdiff1d(all_anchor_ids, center_ids)\n"
        "    if len(center_ids) < 5 or len(non_center_ids) < 5:\n"
        "        continue\n"
        "    de_df = sp.de_genes(\n"
        "        ind_group1=center_ids, ind_group2=non_center_ids,\n"
        "        method='{de_method}', min_fraction=0.05,\n"
        "    )\n"
        "    de_results[motif_key] = de_df\n"
    ),
}

CORR_TEMPLATE = {
    'body': (
        "corr_results = {{}}\n"
        "corr_ids = {{}}\n"
        "for idx, row in significant_motifs.{motif_slice}.iterrows():\n"
        "    motif_key = str(sorted(list(row['motifs'])))\n"
        "    motif_list = sorted(list(row['motifs']))\n"
        "    corr_df = sp.compute_gene_gene_correlation_by_type(\n"
        "        ct=anchor_ct, motif=motif_list,\n"
        "        {neighborhood_arg}, min_size=0, min_nonzero=10, alpha=None,\n"
        "    )\n"
        "    corr_results[motif_key] = corr_df\n"
        "    # Also capture cell IDs so plot_gene_pair_spatial can use them later\n"
        "    _ids = sp.get_anchor_motif_cell_ids(\n"
        "        ct=anchor_ct, motif=motif_list, {neighborhood_arg})\n"
        "    corr_ids[motif_key] = _ids\n"
        "    sig_pairs = corr_df[corr_df['if_significant']]\n"
        "    if len(sig_pairs) > 0:\n"
        "        try:\n"
        "            sp.plot_gene_pair_heatmap(sig_pairs)\n"
        "        except ValueError:\n"
        "            print(f'Skipping heatmap for {{motif_key}}: too few significant pairs ({{len(sig_pairs)}})')\n"
    ),
}

FIND_PATTERNS_TEMPLATE = {
    'postprocess': (
        "if len(fp_result) > 0:\n"
        "    sp.plot_fp_heatmap(fp_result)\n"
    ),
}

PLOT_FOV_TEMPLATE = {
    'body': "sp.plot_fov(min_cells_label=50, figsize=(10,5))\n",
}

PLOT_MOTIF_TEMPLATE = {
    'from_enrichment': (
        "for idx, row in significant_motifs.{motif_slice}.iterrows():\n"
        "    motif_list = sorted(list(row['motifs']))\n"
        "    sp.plot_motif_celltype(\n"
        "        ct=anchor_ct, motif=motif_list, max_dist={max_dist}, figsize=(5,5))\n"
    ),
    'specific': (
        "sp.plot_motif_celltype(\n"
        "    ct='{anchor_ct}', motif={motif}, max_dist={max_dist}, figsize=(5,5))\n"
    ),
}

DE_CUSTOM_TEMPLATE = {
    'by_cell_type': (
        "group1_ids = np.where(np.array(sp.labels) == '{ct1}')[0]\n"
        "group2_ids = np.where(np.array(sp.labels) == '{ct2}')[0]\n"
        "_de_df = sp.de_genes(\n"
        "    ind_group1=group1_ids, ind_group2=group2_ids,\n"
        "    method='{de_method}', min_fraction=0.05)\n"
        "de_results = {{'{ct1} vs {ct2}': _de_df}}\n"
    ),
    'from_ids': (
        "pairs = ids['center_neighbor_motif_pair']\n"
        "center_ids = np.unique(pairs[:, 0])\n"
        "neighbor_ids = np.unique(pairs[:, 1])\n"
        "_de_df = sp.de_genes(\n"
        "    ind_group1=center_ids, ind_group2=neighbor_ids,\n"
        "    method='{de_method}', min_fraction=0.05)\n"
        "de_results = {{'center vs neighbor': _de_df}}\n"
    ),
}

SWEEP_TEMPLATE = {
    'preamble': (
        "anchor_ct = '{anchor_ct}'\n"
        "all_anchor_ids = np.where(np.array(sp.labels) == anchor_ct)[0]\n"
        "sweep_param = '{sweep_param}'\n"
        "sweep_values = {sweep_values}\n"
        "sweep_results = {{}}\n"
    ),
    'loop_knn': (
        "neighborhood_method = 'knn'\n"
        "for _sv in sweep_values:\n"
        "    _k = int(_sv) if sweep_param == 'k' else {k_default}\n"
        "    _md = float(_sv) if sweep_param == 'max_dist' else {max_dist_default}\n"
        "    _ms = float(_sv) if sweep_param == 'min_support' else {min_support_default}\n"
        "    _res = sp.motif_enrichment_knn(\n"
        "        ct=anchor_ct, motifs={motifs},\n"
        "        k=_k, min_support=_ms, max_dist=_md,\n"
        "        return_cellID=True,\n"
        "    )\n"
        "    sweep_results[_sv] = _res\n"
    ),
    'loop_dist': (
        "neighborhood_method = 'dist'\n"
        "for _sv in sweep_values:\n"
        "    _md = float(_sv) if sweep_param == 'max_dist' else {max_dist_default}\n"
        "    _ms = float(_sv) if sweep_param == 'min_support' else {min_support_default}\n"
        "    _res = sp.motif_enrichment_dist(\n"
        "        ct=anchor_ct, motifs={motifs},\n"
        "        max_dist=_md, min_size=0, min_support=_ms,\n"
        "        return_cellID=True,\n"
        "    )\n"
        "    sweep_results[_sv] = _res\n"
    ),
    'postprocess': (
        "# Build trajectory: count significant motifs per sweep value\n"
        "sweep_summary = pd.DataFrame([\n"
        "    {'sweep_value': sv,\n"
        "     'n_total': len(df),\n"
        "     'n_significant': int(df['if_significant'].sum()) if 'if_significant' in df.columns else 0}\n"
        "    for sv, df in sweep_results.items()\n"
        "])\n"
        "_fig, _ax = plt.subplots(figsize=(6, 4))\n"
        "_ax.plot(sweep_summary['sweep_value'], sweep_summary['n_significant'], 'o-', label='significant')\n"
        "_ax.plot(sweep_summary['sweep_value'], sweep_summary['n_total'], 's--', alpha=0.5, label='total')\n"
        "_ax.set_xlabel(sweep_param)\n"
        "_ax.set_ylabel('number of motifs')\n"
        "_ax.set_title(f'Sweep of {sweep_param} around {anchor_ct}')\n"
        "_ax.legend()\n"
        "_ax.grid(True, alpha=0.3)\n"
    ),
}

PLOT_GENE_PAIR_TEMPLATE = {
    # Uses corr_results (dict of DataFrames) + corr_ids (dict of ids dicts)
    # that must already be in exec_ctx from a prior corr run.
    'body': (
        "import ast\n"
        "# Pick a motif key from corr_results (first one if user didn't specify)\n"
        "_motif_key = next(iter(corr_results.keys()))\n"
        "corr_df = corr_results[_motif_key]\n"
        "_ids = corr_ids[_motif_key]\n"
        "# Select gene pairs to plot\n"
        "{pair_selection}"
        "_motif_list = ast.literal_eval(_motif_key)\n"
        "sp.plot_gene_pair_spatial(\n"
        "    gene_pairs=_pairs_to_plot, gene_pair_df=corr_df,\n"
        "    ids=_ids, ct=anchor_ct, motif=_motif_list,\n"
        "    figsize=(20, 5),\n"
        ")\n"
    ),
}

FIND_PATTERNS_GRID_TEMPLATE = {
    'body': (
        "fp_grid_result = sp.find_patterns_grid(\n"
        "    max_dist={max_dist}, min_size=0, min_support={min_support},\n"
        "    if_display=True, figsize=(10, 5),\n"
        "    return_cellID=False, return_grid=False,\n"
        ")\n"
    ),
}

# Tier 1 templates registry
# filter_result and niche_freq have no static template — their bodies depend
# on exec_ctx at generation time, so they are assembled entirely by code_gen.
TEMPLATES: Dict[str, dict] = {
    'motif_enrichment': MOTIF_TEMPLATE,
    'de': DE_TEMPLATE,
    'corr': CORR_TEMPLATE,
    'find_patterns': FIND_PATTERNS_TEMPLATE,
    'plot_fov': PLOT_FOV_TEMPLATE,
    'plot_motif': PLOT_MOTIF_TEMPLATE,
    'de_custom': DE_CUSTOM_TEMPLATE,
    'sweep': SWEEP_TEMPLATE,
    'plot_gene_pair_spatial': PLOT_GENE_PAIR_TEMPLATE,
    'find_patterns_grid': FIND_PATTERNS_GRID_TEMPLATE,
}

# ---------------------------------------------------------------------------
# Operation registry
# ---------------------------------------------------------------------------

OPERATIONS: Dict[str, Operation] = {}


def _register(op: Operation):
    OPERATIONS[op.name] = op


_register(Operation(
    name='motif_enrichment',
    description='Run motif enrichment to find enriched cell type neighborhoods around an anchor cell type',
    rag_query=None,
    required_ctx=[],
    produces=['motif_result', 'significant_motifs', 'anchor_ct', 'all_anchor_ids',
              'neighborhood_method'],
    next_suggestions=['de', 'corr', 'plot_motif'],
    tier=1,
))

_register(Operation(
    name='de',
    description='Run differential expression between motif-positive and motif-negative anchor cells',
    rag_query='differential expression de_genes motif positive negative anchor cells',
    required_ctx=['significant_motifs', 'all_anchor_ids'],
    produces=['de_results'],
    next_suggestions=['corr'],
    tier=1,
))

_register(Operation(
    name='corr',
    description='Run gene co-variation analysis between anchor and neighboring cells',
    rag_query='gene co-variation compute_gene_gene_correlation_by_type anchor neighboring cells per type',
    required_ctx=['significant_motifs', 'anchor_ct'],
    produces=['corr_results', 'corr_ids'],
    next_suggestions=['plot_gene_pair_spatial'],
    tier=1,
))

_register(Operation(
    name='find_patterns',
    description='Find all frequent co-occurrence patterns without significance testing',
    rag_query='find_fp_knn find_fp_dist frequent patterns co-occurrence',
    required_ctx=[],
    produces=['fp_result', 'anchor_ct', 'neighborhood_method', 'neighborhood_max_dist',
              'neighborhood_k'],
    next_suggestions=['motif_enrichment'],
    tier=1,
))

_register(Operation(
    name='plot_fov',
    description='Plot the tissue field of view showing all cell types',
    rag_query='plot_fov visualize tissue field of view',
    required_ctx=[],
    produces=[],
    next_suggestions=[],
    tier=1,
))

_register(Operation(
    name='plot_motif',
    description='Plot the spatial distribution of a specific motif in tissue',
    rag_query='plot_motif_celltype spatial distribution motif',
    required_ctx=[],
    produces=[],
    next_suggestions=[],
    tier=1,
))

_register(Operation(
    name='cell_ids',
    description='Get cell IDs for anchor-motif groups without computing correlations',
    rag_query='get_anchor_motif_cell_ids cell indices groups',
    required_ctx=[],
    produces=['ids'],
    next_suggestions=['de', 'plot_motif'],
    tier=1,
))

_register(Operation(
    name='de_custom',
    description='Run differential expression between two user-specified cell groups',
    rag_query='de_genes differential expression custom groups',
    required_ctx=[],
    produces=['de_results'],
    next_suggestions=[],
    tier=1,
))

_register(Operation(
    name='sweep',
    description='Sweep a motif enrichment parameter (max_dist / k / min_support) across multiple values and plot the trajectory',
    rag_query='motif enrichment parameter sweep trajectory max_dist k min_support sensitivity',
    required_ctx=[],
    produces=['sweep_results', 'sweep_summary', 'anchor_ct', 'all_anchor_ids',
              'sweep_param', 'sweep_values', 'neighborhood_method'],
    next_suggestions=['motif_enrichment', 'de', 'corr'],
    tier=1,
))

_register(Operation(
    name='plot_gene_pair_spatial',
    description='Plot the spatial distribution of specific gene pairs from a covariation result',
    rag_query='plot_gene_pair_spatial spatial visualization gene pair covariation',
    required_ctx=['corr_results', 'corr_ids', 'anchor_ct'],
    produces=[],
    next_suggestions=[],
    tier=1,
))

_register(Operation(
    name='find_patterns_grid',
    description='Scan the whole tissue for frequent cell-type patterns (no anchor required)',
    rag_query='find_patterns_grid whole tissue frequent patterns scan grid',
    required_ctx=[],
    produces=['fp_grid_result'],
    next_suggestions=[],
    tier=1,
))

_register(Operation(
    name='filter_result',
    description='Filter the most recent result (motifs / corr / de) by a criterion',
    rag_query='filter result motifs containing cell type top N gene pairs DE genes threshold',
    required_ctx=[],
    produces=['filtered_result'],
    next_suggestions=[],
    tier=1,
))

_register(Operation(
    name='niche_freq',
    description='Compute and plot niche pattern composition frequency around motif-positive anchor cells',
    rag_query='retrieve_niche_pattern_freq niche pattern frequency composition heatmap',
    required_ctx=[],
    produces=['niche_freqs'],
    next_suggestions=[],
    tier=1,
))

_register(Operation(
    name='pipeline',
    description='Run the full analysis pipeline: motif enrichment → DE → gene co-variation',
    rag_query=None,
    required_ctx=[],
    produces=[],
    next_suggestions=[],
    tier=1,
))

_register(Operation(
    name='question',
    description='Answer a general question about the SpatialQuery API',
    rag_query=None,
    required_ctx=[],
    produces=[],
    next_suggestions=[],
    tier=3,
))
