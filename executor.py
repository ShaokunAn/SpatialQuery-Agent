import io
import traceback
import tempfile
from dataclasses import dataclass, field
from typing import List, Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from SpatialQuery import (
    spatial_query, spatial_query_multi,
    retrieve_niche_pattern_freq, plot_niche_pattern_freq,
)

# Switch to non-interactive backend so figures can be saved without a display.
plt.switch_backend("agg")

# Names that are always in the initial exec_ctx and should never be overwritten.
_BUILTIN_NAMES = {
    'np', 'pd', 'plt', 'adata', 'sp', 'spatial_query', 'spatial_query_multi',
    'retrieve_niche_pattern_freq', 'plot_niche_pattern_freq',
}

# Types worth persisting across steps.
_PERSIST_TYPES = (pd.DataFrame, np.ndarray, list, dict, str, int, float, bool)


@dataclass
class ExecutionResult:
    stdout: str = ""
    dataframes: List[Tuple[str, pd.DataFrame]] = field(default_factory=list)
    figure_paths: List[str] = field(default_factory=list)
    error: Optional[str] = None


def build_exec_context(adata, summary: dict) -> dict:
    """Build the execution namespace for generated analysis code.

    Pre-initializes a spatial_query instance (sp) so generated code can
    call sp.method(...) without constructing the object itself.
    """
    spatial_key = summary.get("inferred_spatial_key") or "X_spatial"
    label_key = summary.get("inferred_label_key") or "cell_type"
    feature_name = summary.get("inferred_feature_name")

    sp = spatial_query(
        adata=adata,
        spatial_key=spatial_key,
        label_key=label_key,
        feature_name=feature_name,
    )

    return {
        "np": np,
        "pd": pd,
        "plt": plt,
        "adata": adata,
        "sp": sp,
        "spatial_query": spatial_query,
        "spatial_query_multi": spatial_query_multi,
        "retrieve_niche_pattern_freq": retrieve_niche_pattern_freq,
        "plot_niche_pattern_freq": plot_niche_pattern_freq,
    }


def execute_code(code: str, base_context: dict) -> ExecutionResult:
    """Execute a code string inside base_context and capture all outputs.

    After execution, all new user-defined variables of persitable types are
    written back into base_context so subsequent steps can access them regardless
    of what variable names the LLM chose to use.
    """
    result = ExecutionResult()
    stdout_buffer = io.StringIO()
    existing_figs = set(plt.get_fignums())
    original_ctx_keys = set(base_context.keys())
    namespace = dict(base_context)

    try:
        from contextlib import redirect_stdout
        with redirect_stdout(stdout_buffer):
            exec(compile(code, "<generated>", "exec"), namespace)  # noqa: S102

        result.stdout = stdout_buffer.getvalue()

        # Persist all new user-defined variables back to base_context.
        for name, val in namespace.items():
            if name.startswith("_") or name in _BUILTIN_NAMES:
                continue
            if name in original_ctx_keys:
                # Overwrite only if the value changed (e.g., updated list).
                if val is not base_context[name]:
                    base_context[name] = val
                continue
            if isinstance(val, _PERSIST_TYPES):
                base_context[name] = val

        # Collect new DataFrames for display.
        for name, val in namespace.items():
            if name.startswith("_") or name in original_ctx_keys:
                continue
            if isinstance(val, pd.DataFrame):
                result.dataframes.append((name, val))

        # Save any new matplotlib figures to temp PNG files.
        new_figs = set(plt.get_fignums()) - existing_figs
        for fig_num in sorted(new_figs):
            fig = plt.figure(fig_num)
            tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
            fig.savefig(tmp.name, bbox_inches="tight", dpi=150)
            tmp.close()
            result.figure_paths.append(tmp.name)
            plt.close(fig)

    except Exception:
        result.error = traceback.format_exc()

    return result


def clear_downstream(op_name: str, exec_ctx: dict):
    """Clear this operation's outputs and cascade to downstream operations."""
    from operations import OPERATIONS
    op = OPERATIONS.get(op_name)
    if not op:
        return
    for key in op.produces:
        exec_ctx.pop(key, None)
    # Cascade: find operations whose required_ctx includes any of op.produces
    for other_op in OPERATIONS.values():
        if other_op.name == op_name:
            continue
        if any(k in op.produces for k in other_op.required_ctx):
            clear_downstream(other_op.name, exec_ctx)
