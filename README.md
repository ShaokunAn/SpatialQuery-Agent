# SpatialQuery-Agent

A natural-language agent for spatial transcriptomics analysis powered by [SpatialQuery](https://github.com/ShaokunAn/Spatial-Query). Translates conversational requests into executable SpatialQuery API calls through template-based code generation — no free-form LLM code generation.

## Architecture

```
User Message → Route → Classify → Parse → Generate → Execute → Display
```

- **Three-layer routing**: Regex fast-path → Embedding similarity → LLM arbitration
- **Template-based code generation**: All analysis code assembled from structured templates
- **Stateful sessions**: Execution context persists across operations with automatic dependency tracking

## Prerequisites

- Python 3.10+
- [Ollama](https://ollama.ai/) with `qwen2.5-coder:7b` model
- [SpatialQuery](https://github.com/ShaokunAn/Spatial-Query) package installed

## Installation

```bash
# 1. Clone this repository
git clone https://github.com/YOUR_USERNAME/SpatialQuery-Agent.git
cd SpatialQuery-Agent

# 2. Install dependencies
pip install -r requirements.txt

# 3. Install SpatialQuery (follow instructions at https://github.com/ShaokunAn/Spatial-Query)

# 4. Pull the LLM model
ollama pull qwen2.5-coder:7b

# 5. Build the knowledge base vector store
python build_db.py

# 6. Start the agent
chainlit run app.py
```

## Usage

Open the Chainlit UI in your browser (default: `http://localhost:8000`).

### Load data
```
/path/to/your/data.h5ad
```

### Run analyses via natural language
```
Find enriched neighborhoods around Endothelium cells
Run DE for the first significant motif, show top 20 genes
Sweep max_dist from 5 to 30
Plot the first motif
Run a full analysis pipeline for gut tube cells
```

### Mode switching
- `agent:` prefix — Analysis mode (default): executes SpatialQuery operations
- `chat:` prefix — Q&A mode: answers questions about the SpatialQuery API using RAG

## Supported Operations

| Category | Operations | Description |
|----------|-----------|-------------|
| Discovery | `motif_enrichment`, `find_patterns`, `find_patterns_grid` | Spatial co-occurrence patterns and significance testing |
| Downstream | `de`, `de_custom`, `corr` | Differential expression and gene co-variation |
| Visualization | `plot_fov`, `plot_motif`, `plot_gene_pair_spatial` | Spatial maps, motif distributions, gene pair heatmaps |
| Exploration | `sweep`, `filter_result`, `niche_freq`, `cell_ids` | Parameter sensitivity, result filtering, niche composition |
| Automation | `pipeline` | End-to-end: motif enrichment → DE → gene co-variation |

## Project Structure

```
SpatialQuery-Agent/
├── app.py                  # Chainlit UI and main orchestration
├── router.py               # Three-layer message routing
├── intent.py               # LLM-based intent classification and parameter extraction
├── embed_classify.py       # Embedding-based operation classifier
├── operations.py           # Operation registry and code templates
├── code_gen.py             # Template-based code assembly
├── executor.py             # Sandboxed code execution and state management
├── display.py              # Result rendering (tables, figures)
├── session_state.py        # Session state management
├── data_tools.py           # AnnData loading and inspection
├── multi_anchor_helpers.py # Multi-anchor analysis utilities
├── build_db.py             # ChromaDB vector store builder
├── knowledge_base/         # API documentation for RAG
│   ├── functions/          # Per-method documentation (34 files)
│   └── workflows/          # Workflow documentation (6 files)
├── chainlit.md             # Chainlit welcome page
├── requirements.txt        # Python dependencies
└── .chainlit/config.toml   # Chainlit configuration
```

## Citation

If you use this software in your research, please cite:

```
[Citation information to be added upon publication]
```
