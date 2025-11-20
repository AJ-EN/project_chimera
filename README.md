# 🧬 Project Chimera

Autonomous Drug Discovery Agent powered by Google Gemini, LangChain, and RDKit.

## 🚀 How to Run

1.  **Install Dependencies** (if you haven't already):
    ```bash
    pip install -r requirements.txt
    ```

2.  **Run the Application**:
    ```bash
    streamlit run main.py
    ```
    The app will open in your browser at `http://localhost:8501`.

## 🧪 Features

-   **Librarian Agent**: Searches scientific literature (PDFs in `data/`) for biological targets.
-   **Lab Rat Agent**: Simulates drug screening, calculates molecular properties (MW, LogP), and estimates binding affinity.
-   **Orchestrator**: The central AI brain that coordinates the agents.
-   **Visualization**: 2D molecule structures are automatically generated and displayed in the sidebar.

## 📂 Data Ingestion

To add new research papers:
1.  Place PDF files in the `data/` folder.
2.  Re-build the knowledge base:
    ```bash
    python agents/librarian.py
    ```

Find a cure for E. Coli.
