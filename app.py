# app.py
"""SpatialQuery Agent — Chainlit UI layer."""
import asyncio
import time
from typing import Optional

import chainlit as cl
import chromadb
import numpy as np
import pandas as pd
from llama_index.core import VectorStoreIndex, Settings
from llama_index.vector_stores.chroma import ChromaVectorStore
from llama_index.embeddings.huggingface import HuggingFaceEmbedding
from llama_index.llms.ollama import Ollama

from data_tools import (load_adata, inspect_adata, format_summary_for_llm,
                        format_params_for_confirmation, update_summary_params)
from session_state import SessionState, sync_after_op
from executor import build_exec_context, execute_code, ExecutionResult, clear_downstream
from intent import (IntentParams, parse_intent, classify_and_parse,
                     validate_cell_types_in_code, closest_cell_type)
from operations import OPERATIONS, Operation
from router import route, RouteResult, check_deps, find_prereq, build_menu, format_menu, CONTINUATION_PHRASES
from embed_classify import EmbeddingClassifier
from code_gen import assemble_code, fix_code
from display import display_result, parse_error_message
from multi_anchor_helpers import snapshot_motif_summary, format_summary_table

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
LLM_MODEL = "qwen2.5-coder:7b"
EMBED_MODEL = "BAAI/bge-base-en-v1.5"
CHROMA_PATH = "./chroma_db"
CHROMA_COLLECTION = "spatial_query"
LLM_TIMEOUT = 120.0
RAG_TOP_K = 3

# ---------------------------------------------------------------------------
# Model initialisation
# ---------------------------------------------------------------------------
print("Initializing models...")
embed_model = HuggingFaceEmbedding(model_name=EMBED_MODEL)
llm = Ollama(model=LLM_MODEL, request_timeout=LLM_TIMEOUT)
Settings.embed_model = embed_model
Settings.llm = llm

print("Building embedding classifier index...")
classifier = EmbeddingClassifier(embed_model)

HISTORY_WINDOW = 6


# ---------------------------------------------------------------------------
# Chat lifecycle
# ---------------------------------------------------------------------------

@cl.on_chat_start
async def start():
    msg = cl.Message(content="Starting SpatialQuery assistant, please wait...")
    await msg.send()

    db = chromadb.PersistentClient(path=CHROMA_PATH)
    chroma_collection = db.get_or_create_collection(CHROMA_COLLECTION)
    vector_store = ChromaVectorStore(chroma_collection=chroma_collection)
    index = VectorStoreIndex.from_vector_store(vector_store, embed_model=embed_model)
    query_engine = index.as_query_engine(similarity_top_k=RAG_TOP_K, streaming=True)

    state = SessionState()
    cl.user_session.set("query_engine", query_engine)
    cl.user_session.set("state", state)

    msg.content = "Hello! I am SpatialQuery Assistant. Connected to local knowledge base."
    await msg.update()

    res = await cl.AskUserMessage(
        content="Please enter the full path to your `.h5ad` file (or type `skip` for Q&A only).",
        timeout=120,
    ).send()

    path = res["output"].strip() if res else None

    if path and path.lower() != "skip":
        loading_msg = await cl.Message(content=f"Loading `{path}`...").send()
        try:
            adata = await asyncio.to_thread(load_adata, path)
            summary = inspect_adata(adata)
            loading_msg.content = (
                f"Data loaded: **{path.split('/')[-1]}** — "
                f"{summary['n_cells']} cells × {summary['n_genes']} genes"
            )
            await loading_msg.update()

            params_msg = format_params_for_confirmation(summary)
            confirm = await cl.AskUserMessage(
                content=(
                    f"{params_msg}\n\n"
                    "Type `ok` if correct. To fix, type only the parameters that need changing, e.g.:\n"
                    "`label_key=predicted_label feature_name=gene_name`"
                ),
                timeout=120,
            ).send()

            if confirm and confirm["output"].strip().lower() != "ok":
                summary = _parse_param_corrections(confirm["output"].strip(), summary, adata)

            summary_text = format_summary_for_llm(summary, filename=path.split("/")[-1])
            exec_ctx = await asyncio.to_thread(build_exec_context, adata, summary)
            exec_ctx['_cell_types'] = summary.get("cell_types", [])

            state.adata = adata
            state.data_summary = summary_text
            state.exec_ctx = exec_ctx
            state.cell_types = summary.get("cell_types", [])
            state.mode = "agent"

            await cl.Message(
                content=(
                    f"Data loaded.\n\n```\n{summary_text}\n```\n\n"
                    "**Agent mode.** Tell me what to analyse, or say "
                    "**\"run full analysis on [cell type]\"** for auto-pilot.\n\n"
                    "_Type `chat` to switch to chatbot mode, or `chat: <question>` to switch and ask in one step._"
                )
            ).send()

        except FileNotFoundError:
            await cl.Message(content=f"File not found: `{path}`").send()
        except Exception as e:
            await cl.Message(content=f"Failed to load file: {e}").send()
    else:
        await cl.Message(content="No file loaded. You can still ask general SpatialQuery questions.").send()


def _parse_param_corrections(user_input: str, summary: dict, adata) -> dict:
    spatial_key = summary["inferred_spatial_key"]
    label_key = summary["inferred_label_key"]
    feature_name = summary["inferred_feature_name"]
    for token in user_input.split():
        if "=" not in token:
            continue
        key, _, val = token.partition("=")
        key, val = key.strip(), val.strip()
        if key == "spatial_key":
            spatial_key = val
        elif key == "label_key":
            label_key = val
        elif key == "feature_name":
            feature_name = val
    return update_summary_params(summary, spatial_key, label_key, feature_name, adata)


# ---------------------------------------------------------------------------
# Message handler
# ---------------------------------------------------------------------------

@cl.on_message
async def main(message: cl.Message):
    state: SessionState = cl.user_session.get("state")
    query_engine = cl.user_session.get("query_engine")

    state.chat_history.append({"role": "user", "content": message.content})

    if not state.has_data:
        reply = await _answer_question(message.content, query_engine, state)
        _append_assistant(state.chat_history, reply)
        return

    # Mode prefix detection
    effective_message = await _parse_mode_prefix(message.content, state)
    if effective_message is None:  # bare "agent"/"chat" — already handled
        return

    if state.mode == "chatbot":
        reply = await _answer_question(effective_message, query_engine, state)
        _append_assistant(state.chat_history, reply)
        return

    # --- Agent mode ---
    route_result = route(effective_message, state)

    if route_result.status == 'needs_llm':
        candidates = classifier.rank(effective_message)

        if _is_clear_embedding_match(candidates):
            # Embedding is confident → use top match, extract params separately
            op = candidates[0][0]
            intent = await parse_intent(effective_message, state, llm)
        else:
            # Ambiguous or no candidates → LLM picks operation + params together
            op, intent = await classify_and_parse(
                effective_message, state, candidates, llm)

        if not op or op.name == 'question':
            reply = await _answer_question(effective_message, query_engine, state)
            _append_assistant(state.chat_history, reply)
            return

        if op.name == 'pipeline':
            await _run_pipeline(effective_message, query_engine, state)
            return

        # Dependency check
        dep_msg = check_deps(op, state.exec_ctx or {})
        if dep_msg:
            res = await cl.AskUserMessage(content=dep_msg, timeout=120).send()
            if res and res["output"].strip().lower() in CONTINUATION_PHRASES:
                prereq = find_prereq(op)
                if prereq:
                    intent = await _resolve_anchor_if_needed(prereq, intent, state)
                    if intent:
                        await _run_operation(prereq, intent, effective_message, state)
                        await _run_operation(op, intent, effective_message, state)
            return

        intent = await _resolve_anchor_if_needed(op, intent, state)
        if intent is None:
            return
        if intent.anchor_cts and len(intent.anchor_cts) > 1 and op.name in _NEEDS_ANCHOR:
            await _run_multi_anchor(op, intent, effective_message, state)
        else:
            await _run_operation(op, intent, effective_message, state)
        return

    # Fast-path matched
    op = route_result.op
    if not op:
        return

    if op.name == 'pipeline':
        await _run_pipeline(effective_message, query_engine, state)
        return

    # Parse parameters and run
    if effective_message.strip().lower() not in CONTINUATION_PHRASES:
        intent = await parse_intent(effective_message, state, llm)
    else:
        intent = IntentParams()

    # Dependency check
    dep_msg = check_deps(op, state.exec_ctx or {})
    if dep_msg:
        res = await cl.AskUserMessage(content=dep_msg, timeout=120).send()
        if res and res["output"].strip().lower() in CONTINUATION_PHRASES:
            prereq = find_prereq(op)
            if prereq:
                intent = await _resolve_anchor_if_needed(prereq, intent, state)
                if intent:
                    await _run_operation(prereq, intent, effective_message, state)
                    await _run_operation(op, intent, effective_message, state)
        return

    intent = await _resolve_anchor_if_needed(op, intent, state)
    if intent is None:
        return

    if intent.anchor_cts and len(intent.anchor_cts) > 1 and op.name in _NEEDS_ANCHOR:
        await _run_multi_anchor(op, intent, effective_message, state)
    else:
        await _run_operation(op, intent, effective_message, state)


def _append_assistant(chat_history: list, content: Optional[str]):
    if content:
        chat_history.append({"role": "assistant", "content": content})


# ---------------------------------------------------------------------------
# Mode prefix parsing
# ---------------------------------------------------------------------------

async def _parse_mode_prefix(content: str, state: SessionState) -> Optional[str]:
    """Parse mode prefix from message. Returns effective message, or None if bare mode switch."""
    msg_stripped = content.strip()
    msg_lower = msg_stripped.lower()

    for prefix in ('agent', 'chat'):
        if msg_lower == prefix:
            state.mode = "chatbot" if prefix == "chat" else "agent"
            await cl.Message(
                content=f"Switched to **{prefix} mode**. "
                        + ("Ask me anything." if prefix == "chat" else "Tell me what to analyse.")
            ).send()
            return None
        tag = prefix + ':'
        if msg_lower.startswith(tag):
            remainder = msg_stripped[len(tag):].strip()
            if remainder:
                state.mode = "chatbot" if prefix == "chat" else "agent"
                return remainder

    return content


def _parse_menu_selection(user_input: str, menu_ops: list) -> Optional[Operation]:
    text = user_input.strip()
    try:
        idx = int(text) - 1
        if 0 <= idx < len(menu_ops):
            return menu_ops[idx]
    except ValueError:
        pass
    return None


# ---------------------------------------------------------------------------
# Anchor resolution
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Embedding confidence check
# ---------------------------------------------------------------------------

_CLEAR_MATCH_THRESHOLD = 0.80   # top candidate must score at least this
_CLEAR_MATCH_GAP = 0.08         # gap to #2 must be at least this


def _is_clear_embedding_match(candidates) -> bool:
    """True when embedding has a single high-confidence match — no LLM needed."""
    if not candidates:
        return False
    top_score = candidates[0][1]
    if top_score < _CLEAR_MATCH_THRESHOLD:
        return False
    if len(candidates) >= 2 and (top_score - candidates[1][1]) < _CLEAR_MATCH_GAP:
        return False  # ambiguous — let LLM decide
    return True


_NEEDS_ANCHOR = {'motif_enrichment', 'find_patterns', 'cell_ids', 'de_custom'}


async def _resolve_anchor_if_needed(op, intent, state):
    if op.name in _NEEDS_ANCHOR and not intent.anchor_cts:
        return await _resolve_anchor(intent, state.cell_types)
    return intent


async def _resolve_anchor(intent: IntentParams, cell_types: list) -> Optional[IntentParams]:
    if intent.anchor_cts:
        return intent

    display_cts = cell_types[:30]
    ct_list_md = "\n".join(f"- {ct}" for ct in display_cts)
    suffix = f"\n_(showing first 30 of {len(cell_types)})_" if len(cell_types) > 30 else ""

    res = await cl.AskUserMessage(
        content=(
            "**Which cell type would you like to use as the anchor?**\n\n"
            f"Available cell types:\n{ct_list_md}{suffix}\n\n"
            "You can name multiple cell types separated by commas."
        ),
        timeout=180,
    ).send()

    if not res:
        await cl.Message(content="No anchor cell type provided. Analysis cancelled.").send()
        return None

    raw = res["output"].strip()
    candidates = [s.strip() for s in raw.replace(";", ",").split(",") if s.strip()]
    validated = [closest_cell_type(c, cell_types) for c in candidates]
    validated = [ct for ct in validated if ct is not None]

    if not validated:
        await cl.Message(content=f"Could not match '{raw}' to any known cell type.").send()
        return None

    intent.anchor_cts = validated
    return intent


# ---------------------------------------------------------------------------
# Core operation runner
# ---------------------------------------------------------------------------

async def _run_operation(op, intent, user_message, state):
    """Generate, execute, and display one operation."""
    params_line = _format_params_summary(intent)
    header = f"**{op.description}**"
    if params_line:
        header += f"\n{params_line}"
    await cl.Message(content=header).send()

    clear_downstream(op.name, state.exec_ctx)
    code = assemble_code(op.name, intent, state.exec_ctx)
    code = validate_cell_types_in_code(code, state.cell_types)
    await cl.Message(content=f"```python\n{code}\n```").send()

    # cl.Step shows a spinner at the bottom of the chat during execution
    async with cl.Step(name="Executing analysis", type="run") as step:
        t0 = time.time()
        result = await asyncio.to_thread(execute_code, code, state.exec_ctx)
        elapsed = time.time() - t0
        step.output = f"{'Completed' if not result.error else 'Failed'} in {elapsed:.1f}s"

    if result.error:
        await cl.Message(content=parse_error_message(result.error)).send()
        return

    if op.name == 'motif_enrichment':
        _fixup_motif_ctx(state.exec_ctx)

    sync_after_op(state, op.name, intent)
    await display_result(op.name, result, state.exec_ctx)
    await _suggest_next(op, state)


async def _run_multi_anchor(op, intent, user_message, state):
    """Run an operation for each anchor in intent.anchor_cts, then show a summary.

    Summary is only produced for motif_enrichment (the other _NEEDS_ANCHOR ops
    don't have the right shape for comparison).
    """
    anchors = intent.anchor_cts or []
    if len(anchors) <= 1:
        await _run_operation(op, intent, user_message, state)
        return

    await cl.Message(
        content=f"**Multi-anchor mode:** running {op.name} for {len(anchors)} anchors."
    ).send()

    summaries = []
    for anchor in anchors:
        per_intent = IntentParams(**{**intent.__dict__, 'anchor_cts': [anchor]})
        await _run_operation(op, per_intent, user_message, state)
        if op.name == 'motif_enrichment':
            summaries.append(snapshot_motif_summary(anchor, state.exec_ctx or {}))

    if summaries:
        table = format_summary_table(summaries)
        await cl.Message(
            content="**Multi-anchor motif summary:**\n\n" + table
        ).send()


def _format_params_summary(intent) -> str:
    """One-line summary of parsed parameters for display."""
    parts = []
    if intent.anchor_cts:
        parts.append(f"anchor: {', '.join(intent.anchor_cts)}")
    if intent.motif_cts:
        parts.append(f"motif: {', '.join(intent.motif_cts)}")
    if intent.motif_index is not None:
        parts.append(f"motif #{intent.motif_index + 1}")
    elif intent.n_motifs:
        parts.append(f"top {intent.n_motifs} motifs")
    if intent.neighborhood_method:
        method_label = "KNN" if intent.neighborhood_method == "knn" else "distance"
        parts.append(f"method: {method_label}")
    if intent.neighborhood_max_dist is not None:
        parts.append(f"max_dist: {intent.neighborhood_max_dist}")
    if intent.neighborhood_k is not None:
        parts.append(f"k: {intent.neighborhood_k}")
    if intent.min_support is not None:
        parts.append(f"min_support: {intent.min_support}")
    return " | ".join(parts) if parts else ""


def _fixup_motif_ctx(exec_ctx: dict):
    """Patch any variables the code may have omitted from the motif step."""
    if 'all_anchor_ids' not in exec_ctx and 'anchor_ct' in exec_ctx:
        exec_ctx['all_anchor_ids'] = np.where(
            np.array(exec_ctx['sp'].labels) == exec_ctx['anchor_ct']
        )[0]
    if 'neighborhood_method' not in exec_ctx:
        if 'neighborhood_k' in exec_ctx:
            exec_ctx['neighborhood_method'] = 'knn'
        elif 'neighborhood_max_dist' in exec_ctx:
            exec_ctx['neighborhood_method'] = 'dist'
    if 'significant_motifs' not in exec_ctx:
        _builtin = {'np', 'pd', 'plt', 'adata', 'sp', 'spatial_query', 'spatial_query_multi'}
        for k, v in exec_ctx.items():
            if k.startswith('_') or k in _builtin:
                continue
            if isinstance(v, pd.DataFrame) and 'if_significant' in v.columns:
                if len(v) > 0 and bool(v['if_significant'].all()):
                    exec_ctx['significant_motifs'] = v
                    break
        if 'significant_motifs' not in exec_ctx and 'motif_result' in exec_ctx:
            exec_ctx['significant_motifs'] = exec_ctx['motif_result'][
                exec_ctx['motif_result']['if_significant']
            ]


# ---------------------------------------------------------------------------
# Pipeline mode
# ---------------------------------------------------------------------------

_PIPELINE_STEPS = ['motif_enrichment', 'de', 'corr']
_PIPELINE_LABELS = {
    'motif_enrichment': 'Motif Enrichment',
    'de': 'Differential Expression',
    'corr': 'Co-variation',
}


async def _run_pipeline(user_message, query_engine, state):
    intent = await parse_intent(user_message, state, llm)
    intent = await _resolve_anchor_if_needed(OPERATIONS['motif_enrichment'], intent, state)
    if intent is None:
        return

    total = len(_PIPELINE_STEPS)
    pipeline_msg = cl.Message(
        content=f"**Auto-pilot mode** (0/{total}): starting pipeline..."
    )
    await pipeline_msg.send()

    for i, op_name in enumerate(_PIPELINE_STEPS, 1):
        op = OPERATIONS[op_name]
        pipeline_msg.content = (
            f"**Auto-pilot mode** ({i}/{total}): {_PIPELINE_LABELS[op_name]}"
        )
        await pipeline_msg.update()

        await _run_operation(op, intent, user_message, state)
        if op.produces and not any(k in state.exec_ctx for k in op.produces):
            pipeline_msg.content = (
                f"**Auto-pilot mode** — stopped at step {i}/{total}: "
                f"{_PIPELINE_LABELS[op_name]} produced no results."
            )
            await pipeline_msg.update()
            return

    pipeline_msg.content = f"**Auto-pilot mode** ({total}/{total}): complete."
    await pipeline_msg.update()
    await cl.Message(
        content="**Pipeline complete.** Describe your next analysis or ask a question."
    ).send()
    state.suggested_ops = []


# ---------------------------------------------------------------------------
# Suggestions
# ---------------------------------------------------------------------------

async def _suggest_next(op, state):
    available = []
    ctx = state.exec_ctx or {}
    for name in op.next_suggestions:
        next_op = OPERATIONS.get(name)
        if next_op and all(k in ctx for k in next_op.required_ctx):
            available.append(next_op)

    if not available:
        state.suggested_ops = []
        return

    lines = ["You can:"]
    for next_op in available:
        lines.append(f"- **{next_op.name}** — {next_op.description}")
    lines.append("- Adjust parameters and re-run")
    lines.append("- Or describe any other analysis")

    await cl.Message(content="\n".join(lines)).send()
    state.suggested_ops = [op.name for op in available]


# ---------------------------------------------------------------------------
# RAG Q&A
# ---------------------------------------------------------------------------

def _format_history(chat_history: list) -> str:
    recent = chat_history[-(HISTORY_WINDOW + 1):-1]
    if not recent:
        return ""
    lines = ["Recent conversation:"]
    for msg in recent:
        role = "User" if msg["role"] == "user" else "Assistant"
        lines.append(f"  {role}: {msg['content']}")
    return "\n".join(lines)


async def _answer_question(user_message, query_engine, state) -> str:
    data_summary = state.data_summary if state else None
    chat_history = state.chat_history if state else []
    data_context = f"\n\nLoaded dataset context:\n{data_summary}\n" if data_summary else ""
    history_block = _format_history(chat_history)
    history_section = f"\n{history_block}\n" if history_block else ""
    search_prompt = (
        f"You are an expert in SpatialQuery (a spatial transcriptomics library).{data_context}"
        f"{history_section}\n"
        f'The user asked: "{user_message}"\n\n'
        "Please rephrase this into a technical search query that likely matches API documentation.\n"
        "Search Query:"
    )
    parsed = await llm.acomplete(search_prompt)
    nodes = await query_engine.retriever.aretrieve(parsed.text)
    response = await query_engine.asynthesize(user_message, nodes)

    msg = cl.Message(content="")
    full_response = ""
    async for token in response.response_gen:
        await msg.stream_token(token)
        full_response += token
    await msg.send()
    return full_response
