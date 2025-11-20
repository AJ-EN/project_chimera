import os
from dotenv import load_dotenv
from langchain_google_genai import GoogleGenerativeAIEmbeddings
from langchain_community.vectorstores import FAISS
from langchain_community.document_loaders import TextLoader, PyPDFLoader
from langchain_text_splitters import RecursiveCharacterTextSplitter

class Librarian:
    def __init__(self, data_dir="data"):
        load_dotenv()
        api_key = os.getenv("GOOGLE_API_KEY")
        if not api_key:
            raise ValueError("GOOGLE_API_KEY not found in environment variables")
            
        self.data_dir = data_dir
        self.embeddings = GoogleGenerativeAIEmbeddings(model="models/embedding-001", google_api_key=api_key)
        self.vector_store = None
        self.index_path = "faiss_index"
        
        # Initialize or load the index
        if os.path.exists(self.index_path):
            self.load_index()
        else:
            self.ingest_data()

    def ingest_data(self):
        """Loads data from the data directory and creates a FAISS index."""
        documents = []
        if not os.path.exists(self.data_dir):
            os.makedirs(self.data_dir)
            print(f"Created data directory: {self.data_dir}")
            return

        for filename in os.listdir(self.data_dir):
            file_path = os.path.join(self.data_dir, filename)
            if filename.endswith(".txt"):
                loader = TextLoader(file_path)
                documents.extend(loader.load())
            elif filename.endswith(".pdf"):
                loader = PyPDFLoader(file_path)
                documents.extend(loader.load())
        
        if not documents:
            print("No documents found to ingest.")
            return

        text_splitter = RecursiveCharacterTextSplitter(chunk_size=1000, chunk_overlap=200)
        texts = text_splitter.split_documents(documents)

        print(f"Ingesting {len(texts)} chunks...")
        self.vector_store = FAISS.from_documents(texts, self.embeddings)
        self.vector_store.save_local(self.index_path)
        print("Data ingested and index saved.")

    def load_index(self):
        """Loads the FAISS index from disk."""
        try:
            self.vector_store = FAISS.load_local(
                self.index_path, 
                self.embeddings,
                allow_dangerous_deserialization=True # Safe since we created it
            )
            print("FAISS index loaded.")
        except Exception as e:
            print(f"Error loading index: {e}. Re-ingesting...")
            self.ingest_data()

    def search(self, query):
        """Searches the vector store for relevant context."""
        if not self.vector_store:
            return "Librarian Error: No knowledge base available."
        
        docs = self.vector_store.similarity_search(query, k=2)
        if not docs:
            return "No relevant documents found."
            
        # Combine the content of the retrieved docs
        context = "\n\n".join([doc.page_content for doc in docs])
        return f"Librarian Findings:\n{context}"

# Test
if __name__ == "__main__":
    lib = Librarian()
    print(lib.search("E. Coli target"))
