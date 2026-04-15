# code_gen.py
"""Tiered code generation for SpatialQuery operations."""
import re
from typing import Optional

from llama_index.llms.ollama import Ollama

from intent import IntentParams
from operations import Operation, TEMPLATES


# ---------------------------------------------------------------------------
# Tier 1: Template assembly (no LLM needed for the code body)
# ---------------------------------------------------------------------------

# Default parameter values when user doesn't specify
_DEFAULTS = {
    'method': 'dist',
    'k': 30,
    'max_dist': 20,
    'min_support': 0.5,
    'motifs': 'None',
    'n_motifs': 5,
    'de_method': 'wilcoxon',
    'sweep_param': 'max_dist',
    'sweep_values_max_dist':    [5.0, 10.0, 15.0, 20.0],
    'sweep_values_min_support': [0.3, 0.5, 0.7],
    'sweep_values_k':           [10.0, 20.0, 30.0, 50.0],
}


# ---------------------------------------------------------------------------
# Shared rendering helpers (used by _assemble_* functions)
# ---------------------------------------------------------------------------

def _render_motif_list(cts) -> str:
    """Render a list of cell type names as a Python list literal, or 'None'."""
    if not cts:
        return "None"
    return "[" + ", ".join(f"'{ct}'" for ct in cts) + "]"


def _render_anchor_preamble(anchor_ct: str, with_anchor_ids: bool = True) -> str:
    """Render the anchor cell-type declaration, optionally with all_anchor_ids.

    Args:
        anchor_ct: cell type name (goes into a string literal)
        with_anchor_ids: True for motif_enrichment (needs all_anchor_ids for DE);
                         False for find_patterns / cell_ids
    """
    lines = [f"anchor_ct = {anchor_ct!r}"]
    if with_anchor_ids:
        lines.append("all_anchor_ids = np.where(np.array(sp.labels) == anchor_ct)[0]")
    return "\n".join(lines) + "\n"


def _render_missing_ctx(op_label: str, key_desc: str, prereq_hint: str) -> str:
    """Emit a `raise RuntimeError(...)` line for missing exec_ctx dependencies.

    Used at generation time when an operation's required ctx is absent — we
    emit a code line so the error surfaces at execution time with full context.
    """
    msg = f"{op_label} needs {key_desc} in exec_ctx — run {prereq_hint} first."
    return f"raise RuntimeError({msg!r})\n"


def _render_method_block(call_kind: str, method: str, params: dict, motif_cts) -> str:
    """Render the knn/dist method block for one of three call kinds.

    Args:
        call_kind: one of 'motif_enrichment', 'find_fp', 'get_anchor_motif_cell_ids'
        method:    'knn' or 'dist'
        params:    dict with (some of) 'k', 'max_dist', 'min_support'
        motif_cts: for motif_enrichment and cell_ids — a list of cell type names,
                   OR a raw Python expression string (e.g. 'significant_motifs.iloc[0]...'),
                   OR None. Ignored for find_fp.

    Returns the rendered code block (with trailing newline).
    """
    if method not in ('knn', 'dist'):
        raise ValueError(f"method must be 'knn' or 'dist', got {method!r}")
    if call_kind not in ('motif_enrichment', 'find_fp', 'get_anchor_motif_cell_ids'):
        raise ValueError(f"unknown call_kind: {call_kind!r}")

    # Render motif argument: list → rendered list literal; string → raw expression.
    if isinstance(motif_cts, str):
        motif_arg = motif_cts
    else:
        motif_arg = _render_motif_list(motif_cts)

    if call_kind == 'motif_enrichment':
        if method == 'knn':
            return (
                "neighborhood_method = 'knn'\n"
                f"neighborhood_k = {params['k']}\n"
                "motif_result = sp.motif_enrichment_knn(\n"
                f"    ct=anchor_ct, motifs={motif_arg},\n"
                f"    k={params['k']}, min_support={params['min_support']}, "
                f"max_dist={params['max_dist']},\n"
                "    return_cellID=True,\n"
                ")\n"
            )
        return (
            "neighborhood_method = 'dist'\n"
            f"neighborhood_max_dist = {params['max_dist']}\n"
            "motif_result = sp.motif_enrichment_dist(\n"
            f"    ct=anchor_ct, motifs={motif_arg},\n"
            f"    max_dist={params['max_dist']}, min_size=0, "
            f"min_support={params['min_support']},\n"
            "    return_cellID=True,\n"
            ")\n"
        )

    if call_kind == 'find_fp':
        if method == 'knn':
            return (
                "neighborhood_method = 'knn'\n"
                f"neighborhood_k = {params['k']}\n"
                "fp_result = sp.find_fp_knn(\n"
                f"    ct=anchor_ct, k={params['k']}, "
                f"min_support={params['min_support']}, max_dist={params['max_dist']})\n"
            )
        return (
            "neighborhood_method = 'dist'\n"
            f"neighborhood_max_dist = {params['max_dist']}\n"
            "fp_result = sp.find_fp_dist(\n"
            f"    ct=anchor_ct, max_dist={params['max_dist']}, min_size=0, "
            f"min_support={params['min_support']})\n"
        )

    # call_kind == 'get_anchor_motif_cell_ids'
    if method == 'knn':
        return (
            "ids = sp.get_anchor_motif_cell_ids(\n"
            f"    ct=anchor_ct, motif={motif_arg}, k={params['k']})\n"
            "sp.plot_all_center_motif(ct=anchor_ct, ids=ids)\n"
        )
    return (
        "ids = sp.get_anchor_motif_cell_ids(\n"
        f"    ct=anchor_ct, motif={motif_arg}, max_dist={params['max_dist']})\n"
        "sp.plot_all_center_motif(ct=anchor_ct, ids=ids)\n"
    )


def assemble_code(op_name: str, params: IntentParams, exec_ctx: dict = None) -> str:
    """Assemble executable code from structured templates + parsed parameters.

    Returns the assembled code string. No LLM call.
    """
    if exec_ctx is None:
        exec_ctx = {}
    # .get() because 'cell_ids' has no template entry; unknown ops fall through to the ValueError below.
    template = TEMPLATES.get(op_name)

    if op_name == 'motif_enrichment':
        code = _assemble_motif(template, params)
    elif op_name == 'de':
        code = _assemble_de(template, params)
    elif op_name == 'corr':
        code = _assemble_corr(template, params)
    elif op_name == 'find_patterns':
        code = _assemble_find_patterns(template, params)
    elif op_name == 'plot_fov':
        code = _assemble_plot_fov(template, params)
    elif op_name == 'plot_motif':
        code = _assemble_plot_motif(template, params, exec_ctx)
    elif op_name == 'cell_ids':
        code = _assemble_cell_ids(params, exec_ctx)
    elif op_name == 'de_custom':
        code = _assemble_de_custom(template, params, exec_ctx)
    elif op_name == 'sweep':
        code = _assemble_sweep(template, params)
    elif op_name == 'plot_gene_pair_spatial':
        code = _assemble_plot_gene_pair_spatial(template, params)
    elif op_name == 'find_patterns_grid':
        code = _assemble_find_patterns_grid(template, params)
    elif op_name == 'filter_result':
        code = _assemble_filter(params, exec_ctx)
    elif op_name == 'niche_freq':
        code = _assemble_niche_freq(params, exec_ctx)
    else:
        raise ValueError(f"No template for '{op_name}'")

    # Inject display hint so the display layer can read it from exec_ctx
    if params.display_n is not None:
        code = f"display_n = {params.display_n}\n" + code

    return code


def _assemble_motif(template: dict, params: IntentParams) -> str:
    anchor = params.anchor_cts[0] if params.anchor_cts else '_MISSING_'
    method = params.neighborhood_method or _DEFAULTS['method']
    method_params = {
        'k': params.neighborhood_k or _DEFAULTS['k'],
        'max_dist': params.neighborhood_max_dist or _DEFAULTS['max_dist'],
        'min_support': params.min_support if params.min_support is not None
                       else _DEFAULTS['min_support'],
    }
    parts = [
        _render_anchor_preamble(anchor, with_anchor_ids=True),
        _render_method_block('motif_enrichment', method, method_params, params.motif_cts),
        template['postprocess'],
    ]
    return "".join(parts)


def _motif_slice(params: IntentParams) -> str:
    """Build the pandas slice expression for iterating over significant motifs."""
    if params.motif_index is not None:
        i = params.motif_index
        return f"iloc[{i}:{i + 1}]"
    n = params.n_motifs or _DEFAULTS['n_motifs']
    return f"head({n})"


def _assemble_de(template: dict, params: IntentParams) -> str:
    return template['body'].format(
        motif_slice=_motif_slice(params),
        de_method=_DEFAULTS['de_method'],
    )


def _assemble_corr(template: dict, params: IntentParams) -> str:
    method = params.neighborhood_method or _DEFAULTS['method']
    if method == 'knn':
        neighborhood_arg = f"k=neighborhood_k"
    else:
        neighborhood_arg = f"max_dist=neighborhood_max_dist"
    return template['body'].format(motif_slice=_motif_slice(params),
                                   neighborhood_arg=neighborhood_arg)


def _assemble_find_patterns(template: dict, params: IntentParams) -> str:
    anchor = params.anchor_cts[0] if params.anchor_cts else '_MISSING_'
    method = params.neighborhood_method or _DEFAULTS['method']
    method_params = {
        'k': params.neighborhood_k or _DEFAULTS['k'],
        'max_dist': params.neighborhood_max_dist or _DEFAULTS['max_dist'],
        'min_support': params.min_support if params.min_support is not None
                       else _DEFAULTS['min_support'],
    }
    parts = [
        _render_anchor_preamble(anchor, with_anchor_ids=False),
        _render_method_block('find_fp', method, method_params, motif_cts=None),
        template['postprocess'],
    ]
    return "".join(parts)


def _assemble_plot_fov(template: dict, params: IntentParams) -> str:
    return template['body']


def _assemble_plot_motif(template: dict, params: IntentParams, exec_ctx: dict) -> str:
    max_dist = params.neighborhood_max_dist or _DEFAULTS['max_dist']
    if params.motif_cts and params.anchor_cts:
        motif_repr = "[" + ", ".join(f"'{ct}'" for ct in params.motif_cts) + "]"
        return template['specific'].format(anchor_ct=params.anchor_cts[0], motif=motif_repr, max_dist=max_dist)
    return template['from_enrichment'].format(motif_slice=_motif_slice(params), max_dist=max_dist)


def _assemble_cell_ids(params: IntentParams, exec_ctx: dict) -> str:
    anchor = params.anchor_cts[0] if params.anchor_cts else '_MISSING_'
    method = params.neighborhood_method or _DEFAULTS['method']
    method_params = {
        'k': params.neighborhood_k or _DEFAULTS['k'],
        'max_dist': params.neighborhood_max_dist or _DEFAULTS['max_dist'],
    }
    # Resolve motif: explicit list → list; else look up from significant_motifs by index.
    if params.motif_cts:
        motif_arg = params.motif_cts
    elif 'significant_motifs' in exec_ctx:
        idx = (params.n_motifs - 1) if params.n_motifs else 0
        motif_arg = f"sorted(list(significant_motifs.iloc[{idx}]['motifs']))"
    else:
        motif_arg = '_MISSING_MOTIF_'
    parts = [
        _render_anchor_preamble(anchor, with_anchor_ids=False),
        _render_method_block('get_anchor_motif_cell_ids', method, method_params, motif_arg),
    ]
    return "".join(parts)


def _assemble_de_custom(template: dict, params: IntentParams, exec_ctx: dict) -> str:
    de_method = _DEFAULTS['de_method']
    if 'ids' in exec_ctx:
        return template['from_ids'].format(de_method=de_method)
    ct1 = params.anchor_cts[0] if params.anchor_cts else '_MISSING_CT1_'
    ct2 = params.motif_cts[0] if params.motif_cts else '_MISSING_CT2_'
    return template['by_cell_type'].format(ct1=ct1, ct2=ct2, de_method=de_method)


def _assemble_sweep(template: dict, params: IntentParams) -> str:
    anchor = params.anchor_cts[0] if params.anchor_cts else '_MISSING_'
    sweep_param = params.sweep_param or _DEFAULTS['sweep_param']
    if params.sweep_values:
        sweep_values = [float(v) for v in params.sweep_values]
    else:
        sweep_values = _DEFAULTS[f'sweep_values_{sweep_param}']
    method = params.neighborhood_method or _DEFAULTS['method']
    motifs = _render_motif_list(params.motif_cts)

    preamble = template['preamble'].format(
        anchor_ct=anchor,
        sweep_param=sweep_param,
        sweep_values=sweep_values,
    )
    loop_key = 'loop_knn' if method == 'knn' else 'loop_dist'
    loop = template[loop_key].format(
        motifs=motifs,
        k_default=_DEFAULTS['k'],
        max_dist_default=_DEFAULTS['max_dist'],
        min_support_default=_DEFAULTS['min_support'],
    )
    return preamble + loop + template['postprocess']


def _assemble_find_patterns_grid(template: dict, params: IntentParams) -> str:
    max_dist = params.neighborhood_max_dist or _DEFAULTS['max_dist']
    min_support = (params.min_support
                   if params.min_support is not None
                   else _DEFAULTS['min_support'])
    return template['body'].format(max_dist=max_dist, min_support=min_support)


# filter_target → (exec_ctx key, flat-df local name, entity label, default sort col,
#                  default sort ascending, default threshold col)
_FLAT_FILTER_CONFIG = {
    'corr': {
        'ctx_key': 'corr_results',
        'local': '_all_corr',
        'entity': 'pairs',
        'topn_col': 'abs_combined_score',
        'topn_method': 'nlargest',
        'topn_label': 'gene pairs by abs_combined_score',
        'threshold_default_col': 'q_value_test1',
    },
    'de': {
        'ctx_key': 'de_results',
        'local': '_all_de',
        'entity': 'DE genes',
        'topn_col': 'adj-pval',
        'topn_method': 'nsmallest',
        'topn_label': 'DE genes by adj-pval',
        'threshold_default_col': 'adj-pval',
    },
}

_FILTER_CTX_KEY = {
    'motifs': 'significant_motifs',
    'corr': 'corr_results',
    'de': 'de_results',
}


def _assemble_filter(params: IntentParams, exec_ctx: dict) -> str:
    if not params.filter_target:
        raise ValueError("filter_result requires filter_target")
    if params.filter_target not in _FILTER_CTX_KEY:
        raise ValueError(f"unknown filter_target: {params.filter_target!r}")

    ctx_key = _FILTER_CTX_KEY[params.filter_target]
    kind = params.filter_kind or 'top_n'
    value = params.filter_value

    if ctx_key not in exec_ctx:
        return _render_missing_ctx(
            'filter_result',
            repr(ctx_key),
            f"the {params.filter_target} operation",
        )

    if params.filter_target == 'motifs':
        return _filter_motifs_body(kind, value)
    return _filter_flat_body(params.filter_target, kind, value, params.filter_sort_col)


def _filter_motifs_body(kind: str, value) -> str:
    if kind == 'contains':
        return (
            f"_target = {value!r}\n"
            "filtered_result = significant_motifs[\n"
            "    significant_motifs['motifs'].apply(lambda m: _target in set(m))\n"
            "]\n"
            "print(f'Filtered to {len(filtered_result)} motifs containing {_target!r}')\n"
            "print(filtered_result.head(10).to_string(index=False))\n"
        )
    if kind == 'top_n':
        n = int(value) if value is not None else 5
        return (
            f"filtered_result = significant_motifs.head({n})\n"
            f"print(f'Top {n} motifs:')\n"
            "print(filtered_result.to_string(index=False))\n"
        )
    raise ValueError(f"unsupported filter_kind for motifs: {kind}")


def _filter_flat_body(target: str, kind: str, value, sort_col) -> str:
    """Filter-body for dict-of-DataFrame targets (corr_results / de_results)."""
    cfg = _FLAT_FILTER_CONFIG[target]
    ctx_key = cfg['ctx_key']
    local = cfg['local']

    # Flatten dict-of-DataFrames into one df with an added 'motif' column.
    # .assign() returns a new frame without an explicit deep-copy step.
    header = (
        "_frames = [\n"
        f"    _df.assign(motif=_mk) for _mk, _df in {ctx_key}.items()\n"
        "]\n"
        f"{local} = pd.concat(_frames, ignore_index=True) if _frames else pd.DataFrame()\n"
    )

    if kind == 'contains':
        if target == 'corr':
            return header + (
                f"_gene = {value!r}\n"
                f"filtered_result = {local}[\n"
                f"    ({local}['gene_center'] == _gene) | ({local}['gene_motif'] == _gene)\n"
                "]\n"
                f"print(f'Filtered to {{len(filtered_result)}} pairs containing {{_gene!r}}')\n"
                "print(filtered_result.head(10).to_string(index=False))\n"
            )
        # de: substring match on gene column
        return header + (
            f"_pat = {value!r}\n"
            f"filtered_result = {local}[{local}['gene'].str.contains(_pat, na=False)]\n"
            f"print(f'Filtered to {{len(filtered_result)}} DE genes matching {{_pat!r}}')\n"
            "print(filtered_result.head(10).to_string(index=False))\n"
        )

    if kind == 'top_n':
        n = int(value) if value is not None else 5
        method = cfg['topn_method']  # 'nlargest' or 'nsmallest'
        col = cfg['topn_col']
        label = cfg['topn_label']
        return header + (
            f"filtered_result = {local}.{method}({n}, {col!r})\n"
            f"print(f'Top {n} {label}:')\n"
            "print(filtered_result.to_string(index=False))\n"
        )

    if kind == 'threshold':
        col = sort_col or cfg['threshold_default_col']
        entity = cfg['entity']
        return header + (
            f"filtered_result = {local}[{local}[{col!r}] < {float(value)}]\n"
            f"print(f'Filtered {{len(filtered_result)}} {entity} with {col} < {float(value)}')\n"
            "print(filtered_result.head(10).to_string(index=False))\n"
        )
    raise ValueError(f"unsupported filter_kind for {target}: {kind}")


def _assemble_niche_freq(params: IntentParams, exec_ctx: dict) -> str:
    if not params.anchor_cts:
        raise ValueError("niche_freq requires an anchor cell type")
    anchor = params.anchor_cts[0]
    max_dist = params.neighborhood_max_dist or _DEFAULTS['max_dist']

    if 'significant_motifs' in exec_ctx:
        pattern_var = 'significant_motifs'
    elif 'fp_result' in exec_ctx:
        pattern_var = 'fp_result'
    else:
        return _render_missing_ctx(
            'niche_freq',
            "'significant_motifs' or 'fp_result'",
            'motif_enrichment or find_patterns',
        )

    return (
        _render_anchor_preamble(anchor, with_anchor_ids=False)
        + "niche_freqs = retrieve_niche_pattern_freq(\n"
        f"    {pattern_var}, sp, anchor_ct, max_dist={max_dist},\n"
        ")\n"
        "plot_niche_pattern_freq(niche_freqs)\n"
    )


def _assemble_plot_gene_pair_spatial(template: dict, params: IntentParams) -> str:
    # Build the pair-selection code block based on precedence:
    # gene_pairs > pair_index > n_pairs > default top 3
    if params.gene_pairs:
        pairs_literal = "[" + ", ".join(
            f"({g1!r}, {g2!r})" for g1, g2 in params.gene_pairs
        ) + "]"
        pair_selection = (
            f"_pairs_to_plot = {pairs_literal}\n"
        )
    elif params.pair_index is not None:
        # 1-indexed → iloc[i-1:i]
        i = max(1, int(params.pair_index))
        pair_selection = (
            "_sig = corr_df[corr_df['if_significant']].sort_values("
            "'abs_combined_score', ascending=False)\n"
            f"_top = _sig.iloc[{i - 1}:{i}]\n"
            "_pairs_to_plot = list(zip(_top['gene_center'], _top['gene_motif']))\n"
        )
    elif params.n_pairs:
        n = int(params.n_pairs)
        pair_selection = (
            "_sig = corr_df[corr_df['if_significant']].sort_values("
            "'abs_combined_score', ascending=False)\n"
            f"_top = _sig.head({n})\n"
            "_pairs_to_plot = list(zip(_top['gene_center'], _top['gene_motif']))\n"
        )
    else:
        pair_selection = (
            "_sig = corr_df[corr_df['if_significant']].sort_values("
            "'abs_combined_score', ascending=False)\n"
            "_top = _sig.head(3)\n"
            "_pairs_to_plot = list(zip(_top['gene_center'], _top['gene_motif']))\n"
        )
    return template['body'].format(pair_selection=pair_selection)


# ---------------------------------------------------------------------------
# Code fix (shared across all tiers)
# ---------------------------------------------------------------------------

async def fix_code(original_code: str, error: str, data_summary: str, llm: Ollama) -> str:
    """Ask LLM to fix code that raised an error."""
    prompt = (
        "The following Python code raised an error. Fix it and return only the corrected code.\n\n"
        f"Dataset context:\n{data_summary}\n\n"
        f"Original code:\n{original_code}\n\n"
        f"Error:\n{error}\n\n"
        "Return only raw Python code, no markdown fences, no explanations."
    )
    response = await llm.acomplete(prompt)
    return strip_code_fences(response.text.strip())


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

def strip_code_fences(text: str) -> str:
    """Remove markdown code fences if the LLM wrapped its output in them."""
    lines = text.splitlines()
    if lines and lines[0].startswith("```"):
        lines = lines[1:]
    if lines and lines[-1].strip() == "```":
        lines = lines[:-1]
    return "\n".join(lines)
