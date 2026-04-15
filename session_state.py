# session_state.py
"""Unified session state — single context object shared across all pipeline stages."""
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class OpRecord:
    """Compact record of one completed operation."""
    op_name: str
    anchor_ct: Optional[str] = None
    method: Optional[str] = None
    max_dist: Optional[float] = None
    k: Optional[int] = None
    motifs: Optional[List[str]] = None


@dataclass
class SessionState:
    """All session context in one object. Passed to every pipeline stage."""

    # Data layer
    adata: Optional[Any] = None
    data_summary: Optional[str] = None
    exec_ctx: Optional[dict] = None
    cell_types: List[str] = field(default_factory=list)

    # Mode
    mode: str = "chatbot"

    # Conversation
    chat_history: List[dict] = field(default_factory=list)

    # Operation tracking
    op_history: List[OpRecord] = field(default_factory=list)
    suggested_ops: List[str] = field(default_factory=list)

    # Last operation parameter snapshot (auto-synced after each op)
    last_op: Optional[str] = None
    last_anchor_ct: Optional[str] = None
    last_method: Optional[str] = None
    last_max_dist: Optional[float] = None
    last_k: Optional[int] = None
    last_min_support: Optional[float] = None
    last_motifs: Optional[List[str]] = None

    @property
    def has_data(self) -> bool:
        return self.exec_ctx is not None and self.data_summary is not None


# Keys in exec_ctx that map to SessionState.last_* fields.
_SYNC_MAP = {
    'anchor_ct': 'last_anchor_ct',
    'neighborhood_method': 'last_method',
    'neighborhood_max_dist': 'last_max_dist',
    'neighborhood_k': 'last_k',
}


def sync_after_op(state: SessionState, op_name: str, intent) -> None:
    """Update state.last_* from exec_ctx + intent after a successful operation.

    Args:
        state: the session state to update
        op_name: name of the operation that just completed
        intent: IntentParams used for the operation
    """
    state.last_op = op_name

    # Sync from exec_ctx (variables the executed code stored)
    ctx = state.exec_ctx or {}
    for ctx_key, state_attr in _SYNC_MAP.items():
        if ctx_key in ctx:
            setattr(state, state_attr, ctx[ctx_key])

    # Sync from intent (covers params the template used but didn't store as variables)
    if intent.anchor_cts:
        state.last_anchor_ct = intent.anchor_cts[0]
    if intent.neighborhood_method:
        state.last_method = intent.neighborhood_method
    if intent.min_support is not None:
        state.last_min_support = intent.min_support
    if intent.motif_cts:
        state.last_motifs = intent.motif_cts

    # Append to history
    state.op_history.append(OpRecord(
        op_name=op_name,
        anchor_ct=state.last_anchor_ct,
        method=state.last_method,
        max_dist=state.last_max_dist,
        k=state.last_k,
        motifs=state.last_motifs,
    ))
