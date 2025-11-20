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

# Custom CSS for that "Cyber-Lab" look
st.markdown("""
<style>
    .stApp { background-color: #0e1117; color: #c9d1d9; }
    .stTextInput input { color: #ffffff; background-color: #262730; }
</style>
""", unsafe_allow_html=True)

st.title("🧬 Project Chimera: Autonomous Drug Discovery")
st.markdown("### `Agents for Good` Track | Built with Google Gemini")

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
            "🧪 *Dr. Chimera is analyzing protocol...*")

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
                    st.sidebar.image(mol_img, caption=f"Molecule: {smiles_str}", use_container_width=True)
                    
                    # Extract other metrics if available
                    mw_match = re.search(r"Molecular Weight: ([\d.]+)", response)
                    ba_match = re.search(r"Binding.*Score.*?: ([\d.]+)", response)
                    
                    if mw_match:
                        st.sidebar.metric("Molecular Weight", f"{mw_match.group(1)} g/mol")
                    if ba_match:
                        st.sidebar.metric("Binding Affinity", ba_match.group(1))

        except Exception as e:
            error_msg = f"System Failure: {str(e)}"
            message_placeholder.markdown(error_msg)
            st.session_state.messages.append(
                {"role": "assistant", "content": error_msg})
