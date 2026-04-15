# intent.py
"""Combined operation classification + intent parsing via LLM."""
import difflib
import json
import re
from dataclasses import dataclass, field
from typing import Any, List, Optional, Tuple

from llama_index.llms.ollama import Ollama
from operations import OPERATIONS, Operation


@dataclass
class IntentParams:
    anchor_cts:            Optional[List[str]] = field(default=None)
    motif_cts:             Optional[List[str]] = field(default=None)
    neighborhood_method:   Optional[str]       = None
    neighborhood_k:        Optional[int]       = None
    neighborhood_max_dist: Optional[float]     = None
    min_support:           Optional[float]     = None
    n_motifs:              Optional[int]       = None
    motif_index:           Optional[int]       = None   # 0-indexed: "second motif" → 1

    # --- Phase B additions ---
    # Sweep (op: sweep)
    sweep_param:  Optional[str]         = None   # 'max_dist' | 'k' | 'min_support'
    sweep_values: Optional[List[float]] = None

    # plot_gene_pair_spatial
    gene_pairs: Optional[List[Tuple[str, str]]] = None
    pair_index: Optional[int]                   = None   # 1-indexed, "the 5th pair"
    n_pairs:    Optional[int]                   = None   # "top N pairs"

    # display hint (e.g. "show top 20 DE genes")
    display_n: Optional[int] = None

    # filter_result
    filter_target:   Optional[str] = None   # 'motifs' | 'corr' | 'de'
    filter_kind:     Optional[str] = None   # 'contains' | 'top_n' | 'threshold'
    filter_value:    Optional[Any] = None
    filter_sort_col: Optional[str] = None


# Patterns that indicate the user is referencing motifs from previous results,
# not specifying new cell types.
_MOTIF_REF_PATTERNS = [
    re.compile(r'\b(first|second|third|1st|2nd|3rd|\d+th)\s+(?:\w+\s+)?(motifs?|patterns?)', re.I),
    re.compile(r'\b(the|this|these|those)\s+(motifs?|patterns?|significant)', re.I),
    re.compile(r'\btop\s+\d+\s+(?:\w+\s+)?(motifs?|patterns?)', re.I),
    re.compile(r'\bfor\s+(the\s+)?(first|second|third|\d+)\s+(?:\w+\s+)?(motifs?|patterns?)', re.I),
    re.compile(r'(第[一二三四五\d]+个|前\d+个).*(motif|模式|pattern)', re.I),
]

_ORDINAL_MAP = {
    'first': 1, '1st': 1, 'second': 2, '2nd': 2,
    'third': 3, '3rd': 3, 'fourth': 4, '4th': 4, 'fifth': 5, '5th': 5,
}


def _has_ordinal_motif_ref(msg: str) -> bool:
    """Detect explicit ordinal or top-N motif references (e.g. 'first motif', 'top 3').

    These are unambiguous — the user is clearly referencing previous results
    by position. Used for post-processing when LLM didn't extract motif_cts.
    """
    # Only ordinal/top-N patterns — NOT vague "the motif" / "this motif"
    return bool(
        _MOTIF_REF_PATTERNS[0].search(msg)   # first/second/third motif
        or _MOTIF_REF_PATTERNS[2].search(msg)  # top N motifs
        or _MOTIF_REF_PATTERNS[3].search(msg)  # for the first/second motif
        or _MOTIF_REF_PATTERNS[4].search(msg)  # Chinese ordinals
    )


def _extract_motif_ref(msg: str):
    """Extract motif reference from message.

    Returns (number, is_top_n):
        is_top_n=True  → "top 3 motifs"  → n_motifs=3
        is_top_n=False → "second motif"  → motif_index=1 (0-indexed)
    """
    msg_lower = msg.lower()
    # "top N motifs/patterns" → process first N motifs
    m = re.search(r'\btop\s+(\d+)\s+(?:\w+\s+)?(?:motifs?|patterns?)', msg_lower)
    if m:
        return int(m.group(1)), True
    # Ordinals adjacent to motif/pattern (allow one optional word between)
    ordinal_alts = '|'.join(re.escape(w) for w in _ORDINAL_MAP)
    m = re.search(
        rf'\b({ordinal_alts})\s+(?:\w+\s+)?(?:motifs?|patterns?)', msg_lower,
    )
    if m:
        return _ORDINAL_MAP[m.group(1)], False
    # Numeric ordinals: "6th motif", "10th pattern"
    m = re.search(r'\b(\d+)(?:st|nd|rd|th)\s+(?:\w+\s+)?(?:motifs?|patterns?)', msg_lower)
    if m:
        return int(m.group(1)), False
    # Default: first motif
    return 1, False


def closest_cell_type(name: str, valid: List[str]) -> Optional[str]:
    """Return the closest valid cell type for *name*, or None if no good match."""
    if not valid:
        return None
    lower_map = {ct.lower(): ct for ct in valid}
    if name.lower() in lower_map:
        return lower_map[name.lower()]
    matches = difflib.get_close_matches(name.lower(), lower_map.keys(), n=1, cutoff=0.6)
    if matches:
        return lower_map[matches[0]]
    return None


async def parse_intent(user_message: str, state, llm: Ollama) -> IntentParams:
    """Extract structured analysis parameters from user_message.

    Uses state.op_history and state.last_* for context-aware extraction.
    Falls back to inheriting from state when LLM can't extract.
    """
    cell_types = state.cell_types or []
    ct_list_str = ", ".join(f'"{ct}"' for ct in cell_types[:30])

    # Build context from SessionState
    context_section = ""
    if state.op_history:
        history_lines = [
            f"  {h.op_name} | anchor={h.anchor_ct} | method={h.method}"
            + (f" | max_dist={h.max_dist}" if h.max_dist else "")
            + (f" | k={h.k}" if h.k else "")
            for h in state.op_history[-5:]
        ]
        context_section += (
            "Previous analysis steps (inherit parameters when user refers to earlier results):\n"
            + "\n".join(history_lines) + "\n\n"
        )

    # Add significant motifs context if available
    exec_ctx = state.exec_ctx or {}
    if 'significant_motifs' in exec_ctx:
        import pandas as pd
        sig = exec_ctx['significant_motifs']
        if isinstance(sig, pd.DataFrame) and len(sig) > 0:
            motif_lines = []
            for i, (_, row) in enumerate(sig.head(5).iterrows()):
                motif_lines.append(f"  Motif {i+1}: {sorted(list(row['motifs']))}")
            context_section += (
                "Previously found significant motifs (user may refer to these by position):\n"
                + "\n".join(motif_lines) + "\n\n"
            )

    ctx_params = []
    if state.last_anchor_ct:
        ctx_params.append(f"  anchor_ct = '{state.last_anchor_ct}'")
    if state.last_method:
        ctx_params.append(f"  neighborhood_method = '{state.last_method}'")
    if state.last_max_dist is not None:
        ctx_params.append(f"  neighborhood_max_dist = {state.last_max_dist}")
    if state.last_k is not None:
        ctx_params.append(f"  neighborhood_k = {state.last_k}")
    if ctx_params:
        context_section += (
            "Current session parameters (reuse unless user specifies differently):\n"
            + "\n".join(ctx_params) + "\n\n"
        )

    prompt = (
        "Extract analysis parameters from the user message below.\n"
        f"Valid cell types (use only these): [{ct_list_str}]\n\n"
        f"{context_section}"
        f'User message: "{user_message}"\n\n'
        "Return a JSON object. Include fields that are explicitly mentioned in the "
        "user message OR that can be inferred from the previous analysis context above. "
        "If the user refers to previous results (e.g. 'the pattern', 'test significance'), "
        "carry forward the anchor cell type and method from the context.\n\n"
        "IMPORTANT: If the user refers to motifs by position (e.g. 'first motif', "
        "'second motif', 'top 3 motifs', 'the motif'), set n_motifs and do NOT "
        "fill motif_cts. The motif will be resolved from previous results.\n"
        "Only set motif_cts when the user explicitly names specific cell types for the motif.\n\n"
        "Omit fields that the user did NOT mention and that cannot be inferred "
        "from context. If the user does not specify a value, omit the field — "
        "the system will use these defaults: neighborhood_method='dist', "
        "neighborhood_k=30, neighborhood_max_dist=20, min_support=0.5.\n"
        "{\n"
        '  "anchor_cts": ["<ct1>", "<ct2>"],\n'
        '  "motif_cts": ["<ct1>", "<ct2>"],\n'
        '  "neighborhood_method": "knn" or "dist" (default: "dist"),\n'
        '  "neighborhood_k": <integer> (default: 30),\n'
        '  "neighborhood_max_dist": <number> (default: 20),\n'
        '  "min_support": <number between 0 and 1> (default: 0.5),\n'
        '  "n_motifs": <integer, how many motifs to process — e.g. "first motif"=1, "top 3"=3>,\n'
        '  "display_n": <integer, how many results to show — e.g. "show top 20"=20>,\n'
        '  "sweep_param": "max_dist" or "k" or "min_support" (which parameter to sweep),\n'
        '  "sweep_values": [<number>, ...] (list of values to sweep over),\n'
        '  "filter_target": "motifs" or "corr" or "de" (what to filter),\n'
        '  "filter_kind": "contains" or "top_n" or "threshold",\n'
        '  "filter_value": <string or number depending on filter_kind>,\n'
        '  "n_pairs": <integer, how many gene pairs to plot>\n'
        "}\n"
        "Return only the JSON object, no explanation, no markdown."
    )
    try:
        response = await llm.acomplete(prompt)
        raw = response.text.strip()
        raw = re.sub(r"^```[a-z]*\n?", "", raw)
        raw = re.sub(r"\n?```$", "", raw)
        raw = re.sub(r"//[^\n]*", "", raw)
        data = json.loads(raw)
    except Exception:
        data = {}

    params = _parse_params_from_data(data, cell_types)

    # Post-processing: ordinal motif references ("first motif", "second motif", "top 3")
    # These are unambiguous — always override LLM, clear motif_cts.
    if _has_ordinal_motif_ref(user_message) and 'significant_motifs' in exec_ctx:
        params.motif_cts = None
        n, is_top_n = _extract_motif_ref(user_message)
        if is_top_n:
            params.n_motifs = n
            params.motif_index = None
        else:
            params.motif_index = n - 1  # convert to 0-indexed
            params.n_motifs = None

    # Safety net: if LLM dumped most cell types into motif_cts, clear it
    if (params.motif_cts and len(cell_types) > 5
            and len(params.motif_cts) > len(cell_types) // 2
            and 'significant_motifs' in exec_ctx):
        params.motif_cts = None

    # Fallback: inherit from state.last_* if LLM didn't extract
    if not params.anchor_cts and state.last_anchor_ct:
        ct = closest_cell_type(state.last_anchor_ct, cell_types)
        if ct:
            params.anchor_cts = [ct]
    if not params.neighborhood_method and state.last_method:
        params.neighborhood_method = state.last_method
    if params.neighborhood_max_dist is None and state.last_max_dist is not None:
        params.neighborhood_max_dist = state.last_max_dist
    if params.neighborhood_k is None and state.last_k is not None:
        params.neighborhood_k = state.last_k

    return params


async def classify_and_parse(
    user_message: str,
    state,
    candidates: List[Tuple[Operation, float]],
    llm: Ollama,
) -> Tuple[Optional[Operation], IntentParams]:
    """Combined operation selection + parameter extraction in one LLM call.

    Uses the embedding shortlist (candidates) to present the LLM with a small
    set of plausible operations. The LLM picks the best one AND extracts
    parameters simultaneously — no heuristics needed.

    Args:
        user_message: raw user message
        state: SessionState
        candidates: embedding shortlist from classifier.rank(), may be empty
        llm: the Ollama LLM instance

    Returns:
        (operation, intent_params) — operation may be None if LLM says 'question'
    """
    cell_types = state.cell_types or []
    ct_list_str = ", ".join(f'"{ct}"' for ct in cell_types[:30])
    exec_ctx = state.exec_ctx or {}

    # Build operation choices for the LLM
    if candidates:
        op_choices = []
        for op, score in candidates:
            deps_status = ""
            if op.required_ctx:
                missing = [k for k in op.required_ctx if k not in exec_ctx]
                if missing:
                    deps_status = f" [requires: {', '.join(missing)}]"
                else:
                    deps_status = " [dependencies satisfied]"
            op_choices.append(f"  - {op.name}: {op.description}{deps_status}")
        ops_section = "Candidate operations (ranked by relevance):\n" + "\n".join(op_choices)
    else:
        # No embedding matches — present all non-excluded ops
        op_choices = []
        for name, op in OPERATIONS.items():
            if name in ('question', 'pipeline'):
                continue
            op_choices.append(f"  - {op.name}: {op.description}")
        ops_section = "Available operations:\n" + "\n".join(op_choices)

    # Build context section (same as parse_intent)
    context_section = ""
    if state.op_history:
        history_lines = [
            f"  {h.op_name} | anchor={h.anchor_ct} | method={h.method}"
            + (f" | max_dist={h.max_dist}" if h.max_dist else "")
            + (f" | k={h.k}" if h.k else "")
            for h in state.op_history[-5:]
        ]
        context_section += (
            "Previous analysis steps:\n"
            + "\n".join(history_lines) + "\n\n"
        )

    if 'significant_motifs' in exec_ctx:
        import pandas as pd
        sig = exec_ctx['significant_motifs']
        if isinstance(sig, pd.DataFrame) and len(sig) > 0:
            motif_lines = []
            for i, (_, row) in enumerate(sig.head(5).iterrows()):
                motif_lines.append(f"  Motif {i+1}: {sorted(list(row['motifs']))}")
            context_section += (
                "Previously found significant motifs:\n"
                + "\n".join(motif_lines) + "\n\n"
            )

    if state.last_op:
        context_section += f"Last operation run: {state.last_op}\n"
    ctx_params = []
    if state.last_anchor_ct:
        ctx_params.append(f"  anchor_ct = '{state.last_anchor_ct}'")
    if state.last_method:
        ctx_params.append(f"  neighborhood_method = '{state.last_method}'")
    if state.last_max_dist is not None:
        ctx_params.append(f"  neighborhood_max_dist = {state.last_max_dist}")
    if state.last_k is not None:
        ctx_params.append(f"  neighborhood_k = {state.last_k}")
    if ctx_params:
        context_section += (
            "Current session parameters:\n"
            + "\n".join(ctx_params) + "\n\n"
        )

    prompt = (
        "You are a routing + parameter extraction system for SpatialQuery, "
        "a spatial transcriptomics analysis tool.\n\n"
        "Given the user's message, do TWO things:\n"
        "1. Pick the best operation from the candidates below\n"
        "2. Extract analysis parameters from the message\n\n"
        f"{ops_section}\n"
        "  - question: answer a general question about SpatialQuery (not an analysis)\n\n"
        f"Valid cell types: [{ct_list_str}]\n\n"
        f"{context_section}"
        f'User message: "{user_message}"\n\n'
        "IMPORTANT RULES:\n"
        "- 'motif_enrichment' = test whether a motif (cell type combination) is "
        "statistically significant around an anchor cell type. Use this when user "
        "says 'significant', 'enriched', 'explore', 'test', or specifies cell types "
        "to check as a motif. This is the FIRST analysis step.\n"
        "- 'de' = DE between motif-positive vs motif-negative ANCHOR cells "
        "(requires significant_motifs from prior motif enrichment). "
        "Pick 'de' whenever user says 'motif', 'first motif', 'significant motif', "
        "or refers to a previous motif enrichment result.\n"
        "- 'de_custom' = DE between two user-specified cell types "
        "(e.g. 'compare T cell vs B cell'). ONLY pick 'de_custom' when user "
        "explicitly names TWO cell types to compare and does NOT mention motifs.\n"
        "- 'corr' = gene co-variation / co-expression between anchor and "
        "neighboring motif cells (requires significant_motifs). Only pick corr when "
        "user explicitly asks for gene correlation, co-variation, or co-expression.\n"
        "- If user specifies motif cell types (e.g. 'motif=[A, B]') AND asks about "
        "significance/enrichment, pick 'motif_enrichment' and set motif_cts\n"
        "- If user says 'the motif' or 'first motif', they refer to previous "
        "motif enrichment results — set n_motifs (not motif_cts)\n"
        "- If an operation's dependencies are NOT satisfied, prefer an operation "
        "whose dependencies ARE satisfied, unless the user's intent is unambiguous\n"
        "- If the message is a question (what/how/why/explain), pick 'question'\n\n"
        "Return a JSON object with these fields:\n"
        "{\n"
        '  "operation": "<operation name from the list above>",\n'
        '  "anchor_cts": ["<ct1>"],\n'
        '  "motif_cts": ["<ct1>", "<ct2>"],\n'
        '  "neighborhood_method": "knn" or "dist" (default: "dist"),\n'
        '  "neighborhood_k": <integer> (default: 30),\n'
        '  "neighborhood_max_dist": <number> (default: 20),\n'
        '  "min_support": <number between 0 and 1> (default: 0.5),\n'
        '  "n_motifs": <integer>,\n'
        '  "display_n": <integer, how many results to show — e.g. "show top 20"=20>,\n'
        '  "sweep_param": "max_dist" or "k" or "min_support" (which parameter to sweep),\n'
        '  "sweep_values": [<number>, ...] (list of values to sweep over),\n'
        '  "filter_target": "motifs" or "corr" or "de" (what to filter),\n'
        '  "filter_kind": "contains" or "top_n" or "threshold",\n'
        '  "filter_value": <string or number depending on filter_kind>,\n'
        '  "n_pairs": <integer, how many gene pairs to plot>\n'
        "}\n"
        "Include only fields with values the user explicitly mentioned or that "
        "can be inferred from context. If the user does not specify a value, "
        "omit the field — the system will apply the defaults shown above. "
        "The 'operation' field is required.\n"
        "Return only the JSON object, no explanation, no markdown."
    )

    try:
        response = await llm.acomplete(prompt)
        raw = response.text.strip()
        raw = re.sub(r"^```[a-z]*\n?", "", raw)
        raw = re.sub(r"\n?```$", "", raw)
        raw = re.sub(r"//[^\n]*", "", raw)
        data = json.loads(raw)
    except Exception:
        data = {}

    # Extract operation
    op_name = data.get("operation", "")
    op = OPERATIONS.get(op_name)
    if not op:
        # Fallback: if candidates exist, use the top one
        if candidates:
            op = candidates[0][0]
        else:
            op = OPERATIONS.get("question")

    # Extract parameters (same logic as parse_intent)
    params = _parse_params_from_data(data, cell_types)

    # Post-processing: ordinal motif references ("first motif", "top 3")
    # Unambiguous — always override LLM, clear motif_cts.
    if _has_ordinal_motif_ref(user_message) and 'significant_motifs' in exec_ctx:
        params.motif_cts = None
        n, is_top_n = _extract_motif_ref(user_message)
        if is_top_n:
            params.n_motifs = n
            params.motif_index = None
        else:
            params.motif_index = n - 1
            params.n_motifs = None

    # Safety net: LLM dumped too many cell types
    if (params.motif_cts and len(cell_types) > 5
            and len(params.motif_cts) > len(cell_types) // 2
            and 'significant_motifs' in exec_ctx):
        params.motif_cts = None

    # Fallback: inherit from state
    if not params.anchor_cts and state.last_anchor_ct:
        ct = closest_cell_type(state.last_anchor_ct, cell_types)
        if ct:
            params.anchor_cts = [ct]
    if not params.neighborhood_method and state.last_method:
        params.neighborhood_method = state.last_method
    if params.neighborhood_max_dist is None and state.last_max_dist is not None:
        params.neighborhood_max_dist = state.last_max_dist
    if params.neighborhood_k is None and state.last_k is not None:
        params.neighborhood_k = state.last_k

    return op, params


def _parse_params_from_data(data: dict, cell_types: List[str]) -> IntentParams:
    """Parse IntentParams from a JSON dict (shared by parse_intent and classify_and_parse)."""
    params = IntentParams()

    raw_anchors = data.get("anchor_cts") or data.get("anchor_ct")
    if isinstance(raw_anchors, str):
        raw_anchors = [raw_anchors]
    if isinstance(raw_anchors, list):
        validated = [closest_cell_type(ct, cell_types) for ct in raw_anchors if isinstance(ct, str)]
        validated = [ct for ct in validated if ct is not None]
        if validated:
            params.anchor_cts = validated

    if "motif_cts" in data and isinstance(data["motif_cts"], list):
        validated = [closest_cell_type(ct, cell_types) for ct in data["motif_cts"] if isinstance(ct, str)]
        validated = [ct for ct in validated if ct is not None]
        if validated:
            params.motif_cts = validated

    if "neighborhood_method" in data and data["neighborhood_method"] in ("knn", "dist"):
        params.neighborhood_method = data["neighborhood_method"]
    if "neighborhood_k" in data:
        try:
            val = int(data["neighborhood_k"])
            if val >= 1:
                params.neighborhood_k = val
        except (TypeError, ValueError):
            pass
    if "neighborhood_max_dist" in data:
        try:
            val = float(data["neighborhood_max_dist"])
            if val > 0:
                params.neighborhood_max_dist = val
        except (TypeError, ValueError):
            pass
    if "min_support" in data:
        try:
            val = float(data["min_support"])
            params.min_support = max(0.0, min(val, 1.0))
        except (TypeError, ValueError):
            pass
    if "n_motifs" in data:
        try:
            val = int(data["n_motifs"])
            if val >= 1:
                params.n_motifs = val
        except (TypeError, ValueError):
            pass
    if "display_n" in data:
        try:
            val = int(data["display_n"])
            if val >= 1:
                params.display_n = val
        except (TypeError, ValueError):
            pass

    # --- Phase B fields ---
    _VALID_SWEEP_PARAMS = {'max_dist', 'k', 'min_support'}
    if data.get('sweep_param') in _VALID_SWEEP_PARAMS:
        params.sweep_param = data['sweep_param']
    if 'sweep_values' in data and isinstance(data['sweep_values'], list):
        try:
            vals = [float(v) for v in data['sweep_values']]
            # Drop non-positive values (k and max_dist must be > 0)
            vals = [v for v in vals if v > 0]
            if vals:
                params.sweep_values = vals
        except (TypeError, ValueError):
            pass

    if 'gene_pairs' in data and isinstance(data['gene_pairs'], list):
        pairs = []
        for pair in data['gene_pairs']:
            if isinstance(pair, (list, tuple)) and len(pair) == 2:
                pairs.append((str(pair[0]), str(pair[1])))
        if pairs:
            params.gene_pairs = pairs

    if 'pair_index' in data:
        try:
            val = int(data['pair_index'])
            if val >= 1:
                params.pair_index = val
        except (TypeError, ValueError):
            pass
    if 'n_pairs' in data:
        try:
            val = int(data['n_pairs'])
            if val >= 1:
                params.n_pairs = val
        except (TypeError, ValueError):
            pass

    _VALID_FILTER_TARGETS = {'motifs', 'corr', 'de'}
    _VALID_FILTER_KINDS = {'contains', 'top_n', 'threshold'}
    if data.get('filter_target') in _VALID_FILTER_TARGETS:
        params.filter_target = data['filter_target']
    if data.get('filter_kind') in _VALID_FILTER_KINDS:
        params.filter_kind = data['filter_kind']
    if 'filter_value' in data:
        params.filter_value = data['filter_value']
    if 'filter_sort_col' in data and isinstance(data['filter_sort_col'], str):
        params.filter_sort_col = data['filter_sort_col']

    return params


_CT_ARG_PATTERNS = [
    re.compile(r'\bct\s*=\s*["\']([^"\']+)["\']'),
    re.compile(r'\banchor_ct\s*=\s*["\']([^"\']+)["\']'),
    re.compile(r'\bmotif\s*=\s*\[([^\]]+)\]'),
    re.compile(r'\bmotif_cts\s*=\s*\[([^\]]+)\]'),
]


def validate_cell_types_in_code(code: str, valid_cell_types: List[str]) -> str:
    """Replace any invalid cell type string literals in *code* with closest valid match."""
    if not valid_cell_types:
        return code

    def fix_single(name: str) -> str:
        c = closest_cell_type(name, valid_cell_types)
        return c if c else name

    def fix_list_literal(match_str: str) -> str:
        def replacer(m):
            fixed = fix_single(m.group(1))
            return m.group(0).replace(m.group(1), fixed)
        return re.sub(r'["\']([^"\']+)["\']', replacer, match_str)

    result = code
    for pattern in _CT_ARG_PATTERNS[:2]:
        def scalar_replacer(m, _p=pattern):
            original = m.group(1)
            fixed = fix_single(original)
            if fixed != original:
                return m.group(0).replace(f"'{original}'", f"'{fixed}'").replace(
                    f'"{original}"', f'"{fixed}"'
                )
            return m.group(0)
        result = pattern.sub(scalar_replacer, result)

    for pattern in _CT_ARG_PATTERNS[2:]:
        def list_replacer(m):
            fixed_inner = fix_list_literal(m.group(1))
            if fixed_inner != m.group(1):
                return m.group(0).replace(m.group(1), fixed_inner)
            return m.group(0)
        result = pattern.sub(list_replacer, result)

    return result
