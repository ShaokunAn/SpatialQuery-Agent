# router.py
"""Message router — fast-path routing for deterministic cases.

Non-deterministic routing (embedding shortlist → LLM selection) is handled by
intent.classify_and_parse() — see intent.py.
"""
import re
from dataclasses import dataclass, field
from typing import List, Optional

from operations import OPERATIONS, Operation
from session_state import SessionState

# ---------------------------------------------------------------------------
# Result type
# ---------------------------------------------------------------------------

@dataclass
class RouteResult:
    op: Optional[Operation]
    status: str  # 'matched' | 'menu' | 'needs_llm'
    suggestion_text: str = ""
    menu_ops: List[Operation] = field(default_factory=list)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CONTINUATION_PHRASES = {
    'yes', 'y', 'continue', 'proceed', 'run it', 'go ahead',
    'next', 'ok', 'sure', 'yep', 'yeah',
}

AUTOPILOT_TRIGGERS = {
    'full analysis', 'complete analysis', 'run all', 'end-to-end',
    'full pipeline', 'autopilot',
}

# Map natural-language phrases to operation names (for disambiguation)
_OP_NAME_PATTERNS = [
    (re.compile(r'\bmotif[_ ]?enrichment\b', re.I), 'motif_enrichment'),
    (re.compile(r'\bfind[_ ]?patterns?[_ ]?grid\b', re.I), 'find_patterns_grid'),
    (re.compile(r'\bfind[_ ]?patterns?\b', re.I), 'find_patterns'),
    (re.compile(r'\bplot[_ ]?fov\b', re.I), 'plot_fov'),
    (re.compile(r'\bplot[_ ]?motif\b', re.I), 'plot_motif'),
    (re.compile(r'\bgene[_ ]?pair[_ ]?spatial\b', re.I), 'plot_gene_pair_spatial'),
    (re.compile(r'\bniche[_ ]?freq\b', re.I), 'niche_freq'),
    (re.compile(r'\bsweep\b', re.I), 'sweep'),
    (re.compile(r'\bcell[_ ]?ids?\b', re.I), 'cell_ids'),
    (re.compile(r'\b(?:de|differential expression)\b', re.I), 'de'),
    (re.compile(r'\bcorr(?:elation)?\b', re.I), 'corr'),
]


def _explicit_op_name(msg: str) -> Optional[str]:
    """Return the operation name explicitly mentioned in *msg*, or None."""
    for pattern, op_name in _OP_NAME_PATTERNS:
        if pattern.search(msg):
            return op_name
    return None


# Follow-up detection patterns
_RERUN_PATTERNS = [
    re.compile(r'\b(re-?run|run again|try again|redo)\b', re.I),
    re.compile(r'(重新跑|再跑一次|重新分析)', re.I),
]
_PARAM_CHANGE_PATTERNS = [
    re.compile(r'\b(change|set|use|switch to|with|adjust)\b.*'
               r'\b(radius|dist|distance|knn|k\s*=|max_dist|support|method)\b', re.I),
    re.compile(r'(把.*改|换成|用.*方法)', re.I),
]
_PLOT_RESULT_PATTERNS = [
    re.compile(r'\b(plot|show|visualize|display|draw)\b.*'
               r'\b(it|this|result|the)\b', re.I),
    re.compile(r'(看一下|展示|画出来|可视化)', re.I),
]

# Ops excluded from menu
_EXCLUDED_OPS = {'question', 'pipeline'}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def route(message: str, state: SessionState) -> RouteResult:
    """Fast-path routing for deterministic cases only.

    Returns status='needs_llm' when no fast path matches — caller should
    use classify_and_parse() for LLM-based routing + intent extraction.

    Order: continuation → autopilot → specific follow-up (rerun/param/plot).
    """
    msg = message.strip()
    msg_lower = msg.lower()

    # Step 1: Continuation
    if msg_lower in CONTINUATION_PHRASES and state.suggested_ops:
        for op_name in state.suggested_ops:
            op = OPERATIONS.get(op_name)
            if op and _deps_satisfied(op, state.exec_ctx or {}):
                return RouteResult(op=op, status='matched')

    # Step 2: Autopilot
    if any(t in msg_lower for t in AUTOPILOT_TRIGGERS):
        return RouteResult(op=OPERATIONS['pipeline'], status='matched')

    # Step 3: Specific follow-up (rerun, param change, plot result)
    followup = _detect_followup_specific(msg, state)
    if followup:
        return RouteResult(op=followup, status='matched')

    # No fast path matched — needs LLM-based routing
    return RouteResult(op=None, status='needs_llm')


def build_menu(state: SessionState) -> List[Operation]:
    """Build the interactive operations menu, ordered by relevance."""
    ctx = state.exec_ctx or {}
    priority = []
    if state.last_op:
        last = OPERATIONS.get(state.last_op)
        if last:
            priority = list(last.next_suggestions)

    result = []
    added = set()

    for name in priority:
        op = OPERATIONS.get(name)
        if op and name not in _EXCLUDED_OPS:
            result.append(op)
            added.add(name)

    for name, op in OPERATIONS.items():
        if name in added or name in _EXCLUDED_OPS:
            continue
        if _deps_satisfied(op, ctx):
            result.append(op)
            added.add(name)

    for name, op in OPERATIONS.items():
        if name in added or name in _EXCLUDED_OPS:
            continue
        result.append(op)

    return result


def format_menu(menu_ops: List[Operation], exec_ctx: dict) -> str:
    """Format the menu for display to the user."""
    lines = ["**Available analyses:**\n"]
    for i, op in enumerate(menu_ops, 1):
        suffix = ""
        missing = [k for k in op.required_ctx if k not in (exec_ctx or {})]
        if missing:
            if 'significant_motifs' in missing or 'all_anchor_ids' in missing:
                suffix = " _(requires motif enrichment first)_"
            else:
                suffix = f" _(requires: {', '.join(missing)})_"
        lines.append(f"  {i}. **{op.name}** — {op.description}{suffix}")
    lines.append("\nType a number to select, or describe what you want to do.")
    return "\n".join(lines)


def check_deps(op: Operation, exec_ctx: dict) -> Optional[str]:
    """Return a user-facing message if op's dependencies are missing, else None."""
    missing = [k for k in op.required_ctx if k not in exec_ctx]
    if not missing:
        return None
    if 'significant_motifs' in missing or 'all_anchor_ids' in missing:
        return (
            "This analysis requires **motif enrichment** results first.\n"
            "Would you like to run motif enrichment now? Reply `yes` or specify parameters."
        )
    return f"Missing required data: {', '.join(missing)}. Please run the prerequisite step first."


def find_prereq(op: Operation) -> Optional[Operation]:
    """Return the operation that produces op's missing requirements, or None."""
    for req_key in op.required_ctx:
        for other in OPERATIONS.values():
            if req_key in other.produces:
                return other
    return None


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _deps_satisfied(op: Operation, exec_ctx: dict) -> bool:
    return all(k in exec_ctx for k in op.required_ctx)


def _detect_followup_specific(msg: str, state: SessionState) -> Optional[Operation]:
    """Detect specific follow-ups: rerun, param change, plot result.

    These have clear intent and should take priority over embedding classification.
    """
    if not state.last_op:
        return None
    last_op = OPERATIONS.get(state.last_op)
    if not last_op:
        return None

    # Rerun / parameter change → same operation, BUT only if the user
    # didn't explicitly name a different operation in the message.
    for p in _RERUN_PATTERNS + _PARAM_CHANGE_PATTERNS:
        if p.search(msg):
            named_op = _explicit_op_name(msg)
            if named_op and named_op != last_op.name:
                return OPERATIONS.get(named_op)
            return last_op

    # Plot result → find a plot op in next_suggestions
    for p in _PLOT_RESULT_PATTERNS:
        if p.search(msg):
            for name in last_op.next_suggestions:
                if name.startswith('plot_'):
                    return OPERATIONS.get(name)
            return OPERATIONS.get('plot_fov')

    return None


