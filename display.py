# display.py
"""Result display and diagnostics for each operation type."""
import re
from typing import Optional, List, Tuple

import chainlit as cl
import pandas as pd

from executor import ExecutionResult


# ---------------------------------------------------------------------------
# Per-operation display functions
# ---------------------------------------------------------------------------

async def display_motif(result: ExecutionResult, exec_ctx: dict):
    """Display motif enrichment results: selected columns top 5 + full heatmap."""
    parts = []
    elements = []

    sig = exec_ctx.get('significant_motifs')
    full = exec_ctx.get('motif_result')

    if isinstance(sig, pd.DataFrame) and len(sig) > 0:
        n_total = len(full) if isinstance(full, pd.DataFrame) else "?"
        n_sig = len(sig)
        display_cols = [c for c in ['motifs', 'n_center_motif', 'n_center', 'n_motif',
                                     'expectation', 'p-values', 'adj-pval', 'if_significant']
                        if c in sig.columns]
        preview = sig[display_cols].head(5).to_string(index=False)
        header = f"**Motif enrichment** ({n_total} motifs tested, {n_sig} significant)"
        suffix = f"\n_(showing top 5 of {n_sig})_" if n_sig > 5 else ""
        parts.append(f"{header}:\n```\n{preview}\n```{suffix}")

    for path in result.figure_paths:
        elements.append(cl.Image(path=path, name="figure", display="inline"))

    diag = _diag_motif(exec_ctx)
    if diag:
        parts.append(diag)

    content = "\n\n".join(parts) if parts else "Step completed (no text output)."
    await cl.Message(content=content, elements=elements).send()


async def display_de(result: ExecutionResult, exec_ctx: dict):
    """Display DE results per motif, split by direction."""
    parts = []
    elements = []

    display_n = exec_ctx.get('display_n', 5)
    de_results = exec_ctx.get('de_results')
    if isinstance(de_results, dict):
        for motif_key, df in de_results.items():
            parts.append(_summarize_de_df(motif_key, df, display_n))

    for path in result.figure_paths:
        elements.append(cl.Image(path=path, name="figure", display="inline"))

    diag = _diag_de(exec_ctx)
    if diag:
        parts.append(diag)

    content = "\n\n".join(parts) if parts else "Step completed (no text output)."
    await cl.Message(content=content, elements=elements).send()


async def display_corr(result: ExecutionResult, exec_ctx: dict):
    """Display co-variation results per motif, grouped by cell type."""
    parts = []
    elements = []

    corr_results = exec_ctx.get('corr_results')
    if isinstance(corr_results, dict):
        for motif_key, df in corr_results.items():
            parts.append(_summarize_corr_df(motif_key, df))

    for path in result.figure_paths:
        elements.append(cl.Image(path=path, name="figure", display="inline"))

    diag = _diag_corr(exec_ctx)
    if diag:
        parts.append(diag)

    content = "\n\n".join(parts) if parts else "Step completed (no text output)."
    await cl.Message(content=content, elements=elements).send()


async def display_generic(result: ExecutionResult, exec_ctx: dict):
    """Generic display: stdout + DataFrames + figures."""
    parts = []
    elements = []

    if result.stdout:
        parts.append(f"**Output:**\n```\n{result.stdout}\n```")

    for name, df in result.dataframes:
        parts.append(_summarize_df(name, df))

    for path in result.figure_paths:
        elements.append(cl.Image(path=path, name="figure", display="inline"))

    content = "\n\n".join(parts) if parts else "Step completed (no text output)."
    await cl.Message(content=content, elements=elements).send()


async def display_figures_only(result: ExecutionResult, exec_ctx: dict):
    """Display only figure images, no text."""
    elements = [cl.Image(path=p, name="figure", display="inline") for p in result.figure_paths]
    content = f"{len(elements)} figure(s) generated." if elements else "No figures generated."
    await cl.Message(content=content, elements=elements).send()


async def display_sweep(result: ExecutionResult, exec_ctx: dict):
    """Display sweep trajectory: summary table + trajectory figure."""
    parts = []
    elements = []

    summary = exec_ctx.get('sweep_summary')
    sweep_param = exec_ctx.get('sweep_param', 'parameter')
    if isinstance(summary, pd.DataFrame) and len(summary) > 0:
        parts.append(
            f"**Sweep of `{sweep_param}`** ({len(summary)} values tested):\n"
            f"```\n{summary.to_string(index=False)}\n```"
        )

    # Detail: also show significant motifs per sweep value
    sweep_results = exec_ctx.get('sweep_results')
    if isinstance(sweep_results, dict):
        detail_lines = []
        for sv, df in sweep_results.items():
            if isinstance(df, pd.DataFrame) and 'if_significant' in df.columns:
                sig = df[df['if_significant']]
                if len(sig) > 0:
                    top_motif = sorted(list(sig.iloc[0]['motifs'])) if 'motifs' in sig.columns else '?'
                    detail_lines.append(f"  {sweep_param}={sv}: {len(sig)} sig; top = {top_motif}")
        if detail_lines:
            parts.append("**Top motif per sweep value:**\n```\n" + "\n".join(detail_lines) + "\n```")

    for path in result.figure_paths:
        elements.append(cl.Image(path=path, name="figure", display="inline"))

    content = "\n\n".join(parts) if parts else "Sweep completed (no text output)."
    await cl.Message(content=content, elements=elements).send()


# Display function registry
DISPLAY_FNS = {
    'motif_enrichment': display_motif,
    'de': display_de,
    'corr': display_corr,
    'find_patterns': display_generic,
    'plot_fov': display_figures_only,
    'plot_motif': display_figures_only,
    'cell_ids': display_generic,
    'de_custom': display_de,
    'sweep': display_sweep,
}


async def display_result(op_name: str, result: ExecutionResult, exec_ctx: dict):
    """Dispatch to the correct display function for an operation."""
    fn = DISPLAY_FNS.get(op_name, display_generic)
    await fn(result, exec_ctx)


# ---------------------------------------------------------------------------
# Summarizers (internal)
# ---------------------------------------------------------------------------

def _summarize_df(name: str, df: pd.DataFrame) -> str:
    n_rows = len(df)
    summary_parts = [f"{n_rows} rows"]
    if 'if_significant' in df.columns:
        summary_parts.append(f"{int(df['if_significant'].sum())} significant")
    if 'adj-pval' in df.columns and n_rows > 0:
        summary_parts.append(f"min adj-pval={df['adj-pval'].min():.2e}")
    summary = ", ".join(summary_parts)
    preview = df.head(10).to_string(index=False)
    return f"**{name}** ({summary}):\n```\n{preview}\n```"


def _summarize_de_df(motif_key: str, df: pd.DataFrame, display_n: int = 5) -> str:
    n_tested = len(df)
    n_sig = int((df['adj-pval'] < 0.05).sum()) if 'adj-pval' in df.columns else 0
    header = f"**DE — {motif_key}** ({n_tested} genes tested, {n_sig} significant at adj-pval < 0.05)"

    if 'de_in' not in df.columns or 'adj-pval' not in df.columns:
        return header + "\n" + _summarize_df(f"DE — {motif_key}", df)

    cols = [c for c in ['gene', 'adj-pval', 'proportion_1', 'proportion_2', 'de_in']
            if c in df.columns]

    sig = df[df['adj-pval'] < 0.05].sort_values('adj-pval')
    up_g1 = sig[sig['de_in'] == 'group1'][cols].head(display_n)
    up_g2 = sig[sig['de_in'] == 'group2'][cols].head(display_n)

    parts = [header]
    if len(up_g1):
        parts.append("*Top up in motif+ cells (de\\_in = group1):*\n```\n"
                     + up_g1.to_string(index=False) + "\n```")
    else:
        parts.append("*No significant up-regulation in motif+ cells.*")

    if len(up_g2):
        parts.append("*Top up in motif− cells (de\\_in = group2):*\n```\n"
                     + up_g2.to_string(index=False) + "\n```")
    else:
        parts.append("*No significant up-regulation in motif− cells.*")

    return "\n\n".join(parts)


_CORR_TOP_PER_CT = 10


def _summarize_corr_df(motif_key: str, df: pd.DataFrame) -> str:
    n_tested = len(df)
    n_sig = int(df['if_significant'].sum()) if 'if_significant' in df.columns else 0
    header = f"**Covariation — {motif_key}** ({n_tested} gene pairs tested, {n_sig} significant)"

    score_col = 'abs_combined_score' if 'abs_combined_score' in df.columns else None
    sig_col = 'if_significant' if 'if_significant' in df.columns else None

    if score_col is None:
        return header + "\n" + _summarize_df(f"Covariation — {motif_key}", df)

    # Display columns (cell_type excluded — shown in group header instead)
    cols = [c for c in ['gene_center', 'gene_motif',
                         'delta_corr_test1', 'abs_combined_score', 'q_value_test1']
            if c in df.columns]

    has_ct = 'cell_type' in df.columns and df['cell_type'].nunique() > 1
    if has_ct:
        lines = [header]
        for ct, group in df.groupby('cell_type', sort=False):
            ct_sig = group[group[sig_col]] if sig_col else group
            n_ct_sig = len(ct_sig)
            if n_ct_sig > 0:
                top = ct_sig.nlargest(_CORR_TOP_PER_CT, score_col)
                suffix = f", showing top {_CORR_TOP_PER_CT}" if n_ct_sig > _CORR_TOP_PER_CT else ""
                lines.append(
                    f"\n**{ct}** ({n_ct_sig} significant{suffix}):\n"
                    f"```\n{top[cols].to_string(index=False)}\n```"
                )
            else:
                lines.append(f"\n**{ct}** — no significant gene pairs.")
        return "\n\n".join(lines)

    # Single cell type or no cell_type column
    if sig_col and n_sig > 0:
        display_df = df[df[sig_col]].nlargest(20, score_col)
        return (header + "\n\n*Top significant gene pairs (by abs\\_combined\\_score):*\n```\n"
                + display_df[cols].to_string(index=False) + "\n```")

    display_df = df.nlargest(5, score_col)
    return (header + "\n> No significant gene pairs found. Showing top 5 by score regardless.\n```\n"
            + display_df[cols].to_string(index=False) + "\n```")


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

def _diag_motif(exec_ctx: dict) -> Optional[str]:
    sig = exec_ctx.get('significant_motifs')
    full = exec_ctx.get('motif_result')
    if isinstance(sig, pd.DataFrame) and len(sig) == 0:
        parts = []
        # Show the full result so user can see p-values for non-significant motifs
        if isinstance(full, pd.DataFrame) and len(full) > 0:
            display_cols = [c for c in ['motifs', 'n_center_motif', 'n_center',
                                         'expectation', 'p-values', 'if_significant']
                            if c in full.columns]
            preview = full[display_cols].head(10).to_string(index=False)
            parts.append(f"> **No significant motifs found.** Full results:\n> ```\n{preview}\n> ```")
        parts.append(
            "> Suggestions:\n"
            "> - Lower `min_support` (e.g., 0.1 or 0.2)\n"
            "> - Increase `max_dist` (e.g., 30 or 50)\n"
            "> - Try a different anchor cell type or motif combination"
        )
        return "\n\n".join(parts)
    return None


def _diag_de(exec_ctx: dict) -> Optional[str]:
    de = exec_ctx.get('de_results')
    if isinstance(de, dict) and all(len(v) == 0 for v in de.values()):
        return "> **All DE results are empty.** Verify that center_ids and non_center_ids are non-empty."
    return None


def _diag_corr(exec_ctx: dict) -> Optional[str]:
    corr = exec_ctx.get('corr_results')
    if isinstance(corr, dict):
        all_empty = all(len(v) == 0 for v in corr.values())
        if all_empty:
            return "> **No gene co-variation results.** Check that enough cells are in the motif."
        none_sig = all(
            int(v['if_significant'].sum()) == 0 if 'if_significant' in v.columns else True
            for v in corr.values()
        )
        if none_sig:
            return "> **No significant gene pairs found.** All correlations below significance threshold."
    return None


# ---------------------------------------------------------------------------
# Error formatting
# ---------------------------------------------------------------------------

_ERROR_PATTERNS = [
    (r"KeyError: ['\"]center_id['\"]",
     "Motif result is missing cell ID columns (`center_id` / `neighbor_id`).",
     "Re-run **motif enrichment** — `return_cellID=True` is required for DE and co-variation."),
    (r"KeyError: ['\"]adj-pval['\"]",
     "Column name mismatch: expected `adj-pval` in the result DataFrame.",
     "Check that the installed SpatialQuery version uses `adj-pval` (hyphen, not underscore)."),
    (r"must contain at least 2",
     "Co-variation requires >= 2 non-center cell types in the motif.",
     "Choose a motif that includes at least 2 different neighbor cell types."),
    (r"(not found|KeyError:.*cell.?type)",
     "Cell type name not found in the dataset.",
     "Use only cell type names listed in the dataset context."),
    (r"NameError: name '(\w+)' is not defined",
     "A variable was used before it was assigned.",
     "Run the prerequisite analysis step first, then retry."),
]


def parse_error_message(error: str) -> str:
    """Map a traceback string to a user-friendly message with a fix hint."""
    for pattern, msg, fix in _ERROR_PATTERNS:
        if re.search(pattern, error, re.IGNORECASE):
            return f"> **Error:** {msg}\n> **Fix:** {fix}"
    last_line = next((l.strip() for l in reversed(error.splitlines()) if l.strip()), error)
    return f"> **Error:** {last_line}\n\n<details><summary>Full traceback</summary>\n\n```\n{error}\n```\n</details>"
