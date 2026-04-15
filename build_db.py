import shutil
import os
from llama_index.core import VectorStoreIndex, SimpleDirectoryReader, StorageContext
from llama_index.vector_stores.chroma import ChromaVectorStore
from llama_index.embeddings.huggingface import HuggingFaceEmbedding
import chromadb

# 1. Set up the embedding model (using your specified BGE-base)
# The model will be downloaded automatically the first time it runs (approximately 500MB)
print("📥 Loading embedding model (BGE-base)...")
embed_model = HuggingFaceEmbedding(model_name="BAAI/bge-base-en-v1.5")

# 2. Read all files in the knowledge_base directory
print("📖 Reading documents...")
documents = SimpleDirectoryReader("./knowledge_base/functions", recursive=True).load_data()

# --- 新增调试代码 ---
print(f"📊 共加载了 {len(documents)} 个文档片段")
if len(documents) > 0:
    print(f"第一篇文档的内容预览: {documents[0].text[:100]}...")
else:
    print("❌ 警告：没有加载到任何文档！请检查路径是否正确！")
    exit() # 如果没读到文件，直接退出，不要建库
# ------------------

# 3. Initialize ChromaDB (local persistent storage)
# For demonstration purposes, delete the old database if it exists (do not do this in a production environment!)
if os.path.exists("./chroma_db"):
    shutil.rmtree("./chroma_db")

db = chromadb.PersistentClient(path="./chroma_db")
chroma_collection = db.get_or_create_collection("spatial_query")
vector_store = ChromaVectorStore(chroma_collection=chroma_collection)
storage_context = StorageContext.from_defaults(vector_store=vector_store)

# 4. Create the index (this is the most crucial step: Embedding -> Vector Store)
print("⚙️ Generating vectors and storing them in ChromaDB...")
index = VectorStoreIndex.from_documents(
    documents,
    storage_context=storage_context,
    embed_model=embed_model
)

print("✅ Database built successfully! Data is saved in the ./chroma_db folder.")