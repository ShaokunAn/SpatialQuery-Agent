# multi_anchor_helpers.py
"""Pure-logic helpers for multi-anchor operation dispatch.
Kept free of chainlit/ollama/LLM imports so tests can import cheaply."""
from typing import Optional

import pandas as pd


def snapshot_motif_summary(anchor_ct: str, exec_ctx: dict) -> dict:
    """Extract a compact per-anchor summary from exec_ctx after a motif run.

    Returns a dict with keys: anchor, n_total, n_significant, top_motif.
    Handles missing/empty exec_ctx gracefully.
    """
    motif_df = exec_ctx.get('motif_result')
    sig_df = exec_ctx.get('significant_motifs')
    summary = {
        'anchor': anchor_ct,
        'n_total': int(len(motif_df)) if isinstance(motif_df, pd.DataFrame) else 0,
        'n_significant': int(len(sig_df)) if isinstance(sig_df, pd.DataFrame) else 0,
        'top_motif': None,
    }
    if isinstance(sig_df, pd.DataFrame) and len(sig_df) > 0 and 'motifs' in sig_df.columns:
        summary['top_motif'] = sorted(list(sig_df.iloc[0]['motifs']))
    return summary


def format_summary_table(summaries: list) -> str:
    """Render a markdown table from a list of summary dicts."""
    if not summaries:
        return ""
    lines = [
        "| Anchor | Total | Significant | Top motif |",
        "| --- | --- | --- | --- |",
    ]
    for s in summaries:
        tm = ", ".join(s['top_motif']) if s['top_motif'] else "_none_"
        lines.append(f"| {s['anchor']} | {s['n_total']} | {s['n_significant']} | {tm} |")
    return "\n".join(lines)
