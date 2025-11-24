import streamlit as st
from orchestrator import Orchestrator
from rdkit import Chem
from rdkit.Chem import Draw
import re


def render_molecule(smiles_string):
    """Draws the molecule to the sidebar or chat"""
    mol = Chem.MolFromSmiles(smiles_string)
    if mol:
        img = Draw.MolToImage(mol, size=(300, 300))
        return img
    return None


# Page Config (Make it look pro)
st.set_page_config(page_title="Project Chimera", page_icon="🧬", layout="wide")

# Custom CSS for that "Medical AI" look
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600&family=JetBrains+Mono:wght@400&display=swap');

    /* Global Reset & Typography */
    .stApp {
        font-family: 'Inter', sans-serif;
        /* Use Streamlit's native background variables to support both modes */
    }
    
    h1, h2, h3, h4, h5, h6 {
        font-weight: 600;
        letter-spacing: -0.02em;
        /* Color will adapt automatically */
    }
    
    /* Sidebar Styling */
    [data-testid="stSidebar"] {
        border-right: 1px solid var(--text-color-10); /* Subtle border */
    }
    
    /* Input Fields */
    .stTextInput > div > div > input {
        border-radius: 8px;
        padding: 12px;
        font-family: 'Inter', sans-serif;
        /* Let Streamlit handle colors, just add border/shape */
        border: 1px solid var(--text-color-20);
    }
    .stTextInput > div > div > input:focus {
        border-color: var(--primary-color);
        box-shadow: 0 0 0 3px rgba(47, 129, 247, 0.2);
    }
    
    /* Chat Interface */
    .stChatMessage {
        padding: 1.5rem;
        border-radius: 12px;
        margin-bottom: 1rem;
        border: 1px solid transparent;
    }
    
    [data-testid="stChatMessageUser"] {
        background-color: var(--secondary-background-color);
        border: 1px solid var(--text-color-10);
    }
    
    [data-testid="stChatMessageAvatarUser"] {
        background-color: #238636; /* Keep green for user */
    }
    
    [data-testid="stChatMessageAssistant"] {
        background-color: transparent;
        border: 1px solid var(--primary-color);
        background-color: rgba(47, 129, 247, 0.05); /* Subtle tint */
    }
    
    /* Custom Metric Cards */
    .metric-container {
        background-color: var(--secondary-background-color);
        border: 1px solid var(--text-color-10);
        border-radius: 8px;
        padding: 12px;
        margin-bottom: 8px;
    }
    .metric-label {
        font-size: 0.75rem;
        color: var(--text-color-60); /* 60% opacity text */
        text-transform: uppercase;
        letter-spacing: 0.05em;
    }
    .metric-value {
        font-family: 'JetBrains Mono', monospace;
        font-size: 1.1rem;
        color: var(--primary-color);
        font-weight: 500;
    }
    
    /* Header */
    .main-header {
        border-bottom: 1px solid var(--text-color-10);
        padding-bottom: 1rem;
        margin-bottom: 2rem;
    }
    .badge {
        background-color: var(--primary-color);
        color: white;
        padding: 4px 8px;
        border-radius: 4px;
        font-size: 0.75rem;
        font-weight: 600;
        margin-left: 10px;
        vertical-align: middle;
    }
    
    /* Helper classes for opacity */
    :root {
        --text-color-10: rgba(128, 128, 128, 0.1);
        --text-color-20: rgba(128, 128, 128, 0.2);
        --text-color-60: rgba(128, 128, 128, 0.6);
    }
</style>
""", unsafe_allow_html=True)

# Header Section
st.markdown("""
<div class="main-header">
    <h1>
        Project Chimera <span class="badge">ALPHA</span>
    </h1>
    <p style="color: var(--text-color-60); margin-top: -10px;">
        Autonomous Drug Discovery Protocol | Powered by Google Gemini
    </p>
</div>
""", unsafe_allow_html=True)

# Sidebar Header
with st.sidebar:
    st.image("https://cdn-icons-png.flaticon.com/512/2083/2083206.png", width=50) # Placeholder Icon
    st.markdown("### Molecule Analysis Unit")
    st.markdown("<div style='color: var(--text-color-60); font-size: 0.8rem;'>Real-time structure visualization and property prediction.</div>", unsafe_allow_html=True)
    st.markdown("---")

# Initialize the Brain
if "bot" not in st.session_state:
    st.session_state.bot = Orchestrator()

# Chat Interface
if "messages" not in st.session_state:
    st.session_state.messages = []

# Display Chat History
for message in st.session_state.messages:
    with st.chat_message(message["role"]):
        st.markdown(message["content"])

# User Input
if prompt := st.chat_input("Enter your research goal (e.g., 'Find inhibitors for Tuberculosis')..."):
    # Add user message to chat history
    st.session_state.messages.append({"role": "user", "content": prompt})
    with st.chat_message("user"):
        st.markdown(prompt)

    # Agent "Thinking" State
    with st.chat_message("assistant"):
        message_placeholder = st.empty()
        message_placeholder.markdown(
            """<span style='color: var(--primary-color); font-family: "JetBrains Mono", monospace;'>⚡ Dr. Chimera is analyzing protocol...</span>""", 
            unsafe_allow_html=True)

        try:
            # Run the Orchestrator
            response = st.session_state.bot.run(prompt)

            # Ensure response is a string (handle edge cases)
            if not isinstance(response, str):
                response = str(response)

            message_placeholder.markdown(response)

            # Save history
            st.session_state.messages.append(
                {"role": "assistant", "content": response})

            # --- VISUALIZATION MAGIC ---
            # Check for SMILES in the response (looking for code blocks or specific patterns)
            # Simple regex to find potential SMILES (this is a basic one, can be improved)
            # For now, we'll look for the specific format from our Lab Rat: "SMILES: `...`"
            smiles_match = re.search(r"SMILES: `([^`]+)`", response)

            if smiles_match:
                smiles_str = smiles_match.group(1)
                mol_img = render_molecule(smiles_str)

                if mol_img:
                    st.sidebar.image(
                        mol_img, caption=f"Molecule: {smiles_str}", use_container_width=True)

                    # Extract other metrics if available
                    mw_match = re.search(
                        r"Molecular Weight: ([\d.]+)", response)
                    ba_match = re.search(
                        r"Binding.*Score.*?: ([\d.]+)", response)

                    if mw_match:
                        st.sidebar.markdown(f"""
                        <div class="metric-container">
                            <div class="metric-label">Molecular Weight</div>
                            <div class="metric-value">{mw_match.group(1)} g/mol</div>
                        </div>
                        """, unsafe_allow_html=True)
                    if ba_match:
                        st.sidebar.markdown(f"""
                        <div class="metric-container">
                            <div class="metric-label">Binding Affinity</div>
                            <div class="metric-value">{ba_match.group(1)}</div>
                        </div>
                        """, unsafe_allow_html=True)

        except Exception as e:
            error_msg = f"System Failure: {str(e)}"
            message_placeholder.markdown(error_msg)
            st.session_state.messages.append(
                {"role": "assistant", "content": error_msg})
