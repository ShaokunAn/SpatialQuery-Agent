# embed_classify.py
"""Embedding-based operation classifier using cosine similarity."""
import numpy as np
from typing import Optional, Dict, List

from operations import OPERATIONS, Operation

# Example queries per operation — the embedding index.
EXAMPLE_QUERIES: Dict[str, List[str]] = {
    'motif_enrichment': [
        "explore the microenvironment of T cell",
        "explore micro environment of gut tube",
        "find enriched cell type neighborhoods",
        "run motif enrichment analysis",
        "what cell types surround the anchor cells",
        "analyze the spatial niche around macrophages",
        "test if the spatial patterns are significant",
        "find significant co-occurrence patterns",
        "neighborhood composition analysis",
        "which cell types are enriched around my anchor",
        "explore if this motif is significant surrounding the anchor",
        "is this cell type combination significantly enriched",
        "test significance of a specific motif around gut tube",
        "check if B cell and macrophage are enriched near T cell",
    ],
    'de': [
        "run differential expression analysis",
        "find differentially expressed genes",
        "perform DE analysis for the first motif",
        "what genes are upregulated in motif-positive cells",
        "differential analysis between motif positive and negative",
        "compare gene expression between groups with and without motif",
        "run DE between cells with and without the motif",
        "differential expression for the first significant motif",
        "DE on motif-positive vs motif-negative anchor cells",
        "find genes differentially expressed in motif neighborhoods",
    ],
    'corr': [
        "compute gene co-variation analysis",
        "gene-gene correlation between anchor and neighbor cells",
        "find co-expressed gene pairs",
        "analyze gene pair correlations in the niche",
        "run co-variation analysis for significant motifs",
        "perform covariation analysis between anchor and motif cells",
        "compute gene correlation for the first motif",
        "run covariation for this motif",
        "gene co-expression between center and neighbor cells",
    ],
    'find_patterns': [
        "find all frequent patterns without significance testing",
        "run find_fp without significance test",
        "discover co-occurrence patterns without enrichment",
        "find frequent cell type patterns only no significance",
        "list all patterns without statistical testing",
    ],
    'plot_fov': [
        "plot the tissue field of view",
        "show the tissue map with all cell types",
        "visualize the spatial layout of cells",
        "display the FOV",
        "draw the cell type map of the tissue",
    ],
    'plot_motif': [
        "plot the spatial distribution of the motif",
        "visualize motif location in tissue",
        "show where the motif appears spatially",
        "display the motif enrichment result on tissue",
        "plot the result of motif analysis",
    ],
    'cell_ids': [
        "get cell IDs for the anchor motif group",
        "extract cell indices for a specific motif",
        "retrieve anchor motif cell identifiers",
        "get the cell IDs without computing correlations",
    ],
    'de_custom': [
        "compare gene expression between T cell and B cell",
        "differential expression between two specific cell types",
        "custom DE analysis between user-defined groups",
        "run DE comparing two cell populations",
    ],
    'sweep': [
        "how does the niche around gut tube change with max_dist",
        "sweep max_dist from 5 to 20 and plot how many motifs remain significant",
        "sensitivity analysis of min_support for T cell neighborhood",
        "trajectory of significant motifs as k increases",
        "how sensitive is the result to max_dist choice",
        "scan parameter values and plot the sweep trajectory",
    ],
    'plot_gene_pair_spatial': [
        "plot the top 3 gene pairs spatially",
        "visualize the 5th significant gene pair",
        "show spatial distribution of HOXA1 and MEOX1",
        "plot gene pair spatial distribution from corr results",
        "display top N significant gene pairs on tissue map",
        "show where the top gene pairs are located",
    ],
    'find_patterns_grid': [
        "scan the whole tissue for frequent patterns",
        "find all co-occurrence patterns anywhere in the FOV",
        "discover patterns without specifying an anchor",
        "grid-based whole-tissue pattern scan",
        "find frequent cell type patterns across the entire tissue",
    ],
    'filter_result': [
        "filter the significant motifs to ones containing T cell",
        "from the motif enrichment result keep only motifs with B cell",
        "top 5 gene pairs by combined_score",
        "DE genes with adj-pval less than 0.01",
        "filter the results to only show pairs with HOXA1",
        "show me the top 10 DE genes",
    ],
    'niche_freq': [
        "show composition frequency around motif-positive anchors",
        "what's the cell type distribution in the motif neighborhoods",
        "niche composition breakdown",
        "retrieve niche pattern frequency",
        "plot niche composition heatmap",
        "cell type frequency around significant motifs",
    ],
}


class EmbeddingClassifier:
    """Classifies user messages to operations via embedding cosine similarity."""

    def __init__(self, embed_model):
        """Pre-compute embeddings for all example queries.

        Args:
            embed_model: Object with get_text_embedding(str) and optionally
                         get_text_embedding_batch(list) methods.
        """
        self._embed_model = embed_model
        self._op_names: List[str] = []
        self._op_embeddings: List[np.ndarray] = []
        self._build_index()

    def _build_index(self):
        for op_name, queries in EXAMPLE_QUERIES.items():
            if op_name not in OPERATIONS:
                continue
            try:
                embeddings = self._embed_model.get_text_embedding_batch(queries)
            except AttributeError:
                embeddings = [self._embed_model.get_text_embedding(q) for q in queries]
            emb_arr = np.array(embeddings, dtype=np.float32)
            # Pre-normalize so classify() needs only dot products, not per-call norms.
            norms = np.linalg.norm(emb_arr, axis=1, keepdims=True)
            norms = np.where(norms == 0, 1.0, norms)
            self._op_names.append(op_name)
            self._op_embeddings.append(emb_arr / norms)

    def classify(self, message: str, threshold: float = 0.7) -> Optional[Operation]:
        """Return the best-matching Operation if similarity >= threshold, else None."""
        results = self.rank(message, threshold=threshold)
        return results[0][0] if results else None

    def rank(self, message: str, threshold: float = 0.7) -> List[tuple]:
        """Return all operations above threshold as (Operation, score) pairs, sorted desc."""
        msg_emb = np.array(self._embed_model.get_text_embedding(message), dtype=np.float32)
        msg_norm = np.linalg.norm(msg_emb)
        if msg_norm == 0:
            return []
        msg_emb = msg_emb / msg_norm

        results = []
        for op_name, op_embs in zip(self._op_names, self._op_embeddings):
            max_sim = float(np.max(op_embs @ msg_emb))
            if max_sim >= threshold:
                op = OPERATIONS.get(op_name)
                if op:
                    results.append((op, max_sim))

        results.sort(key=lambda x: x[1], reverse=True)
        return results
