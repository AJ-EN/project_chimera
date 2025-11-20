# 🧬# 🧬 Project Chimera: The Autonomous Drug Discovery Agent

### 🏆 Submission for Kaggle Agents Intensive Capstone (Track: Agents for Good)

**Project Chimera** is a multi-agent system designed to democratize early-stage drug discovery. By combining Large Language Models (LLMs) with chemoinformatics tools (RDKit) and vector-based literature retrieval (RAG), Chimera acts as an "AI Laboratory Partner" that can research targets, screen molecules, and predict efficacy in seconds.

---

## 🌟 The Problem
Developing a new antibiotic takes **10-15 years** and costs **$1-2 billion**. A major bottleneck is the "Knowledge Gap": researchers must synthesize data from thousands of papers (unstructured data) and test millions of molecules (structured data).

## 💡 The Solution
Chimera bridges this gap using a **Multi-Agent Architecture**:
1. **The Librarian (RAG Agent):** Scans internal PDF repositories to identify biological targets (e.g., "Find the PDB ID for E. Coli Gyrase").
2. **The Lab Rat (Simulation Agent):** Uses Python libraries (`RDKit`) to validate chemical structures and simulate binding affinity.
3. **The Web Researcher:** Connects to Google Search for real-time clinical trial data.
4. **The Orchestrator:** A Gemini-powered brain that plans the workflow and synthesizes the final report.

---

## 🛠️ Tech Stack
* **Brain:** Google Gemini 2.5 Flash
* **Orchestration:** LangChain (Zero-Shot ReAct Agent)
* **Chemistry Engine:** RDKit (Open Source Cheminformatics)
* **Memory:** FAISS (Vector Database) & Google Generative Embeddings
* **Interface:** Streamlit

---

## ⚙️ Architecture
```
[User Query] --> [Orchestrator (Gemini 2.5)]
                        |
        -----------------------------------
        |               |                 |
   [Librarian]     [Lab Rat]      [Web Researcher]
   (Reads PDFs)    (Runs RDKit)   (Google Search)
```

---

## 🚀 How to Run

To add new research papers:
1.  Place PDF files in the `data/` folder.
2.  Re-build the knowledge base:
    ```bash
    python agents/librarian.py
    ```

Find a cure for E. Coli.
