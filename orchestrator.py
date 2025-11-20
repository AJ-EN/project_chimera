import os
from dotenv import load_dotenv
from langchain_google_genai import ChatGoogleGenerativeAI
from langgraph.prebuilt import create_react_agent
from langchain_core.messages import SystemMessage
from agents.librarian import Librarian
from agents.lab_rat import LabRat
from agents.web_researcher import WebResearcher
from langchain_core.tools import tool

load_dotenv()

class Orchestrator:
    def __init__(self):
        self.llm = ChatGoogleGenerativeAI(
            model="gemini-2.5-flash",
            temperature=0.2,
            convert_system_message_to_human=True
        )
        
        self.librarian = Librarian()
        self.lab_rat = LabRat()
        self.web_researcher = WebResearcher()

        # Define tools using @tool decorator for LangGraph compatibility
        @tool
        def librarian_search(query: str) -> str:
            """Find biological targets (e.g., 'Find the target for E. Coli'). Returns text."""
            return self.librarian.search(query)
        
        @tool
        def lab_rat_simulation(molecule_name: str) -> str:
            """Screen drugs. Input: Molecule Name. Returns text report."""
            return self.lab_rat.run_simulation(molecule_name)
        
        @tool
        def web_search(query: str) -> str:
            """Search the web. Returns text."""
            return self.web_researcher.search(query)

        self.tools = [librarian_search, lab_rat_simulation, web_search]
        self.agent_graph = create_react_agent(self.llm, self.tools)

    def run(self, user_query):
        # --- THE FIX IS HERE: IMPROVED SYSTEM PROMPT ---
        system_instruction = """You are Dr. Chimera, an AI computational biologist.
        
        PROTOCOL:
        1. Receive the user's query.
        2. Select the correct tool (Librarian, Lab Rat, or Web Search).
        3. **CRITICAL:** When the tool returns a result, you MUST read it and then write a NATURAL LANGUAGE response.
        4. **NEVER** output the raw JSON or dictionary from the tool directly to the user.
        5. If the tool gives you JSON like {'text': '...'}, extract the 'text' part and summarize it.
        6. Always sound professional and scientific.
        
        Example of BAD Output: {'text': 'Cure found...'}
        Example of GOOD Output: "Based on my research, I have found a potential cure..."
        """
        
        inputs = {"messages": [("system", system_instruction), ("user", user_query)]}
        result = self.agent_graph.invoke(inputs)
        
        # --- DOUBLE SAFETY: CLEAN THE OUTPUT ---
        final_response = result["messages"][-1].content
        
        # Robust cleanup for Gemini's raw tool output format: "[{'type': 'text', ...}]"
        if isinstance(final_response, str):
            cleaned_response = final_response.strip()
            
            # Check if it looks like a list or dict
            if cleaned_response.startswith("[") or cleaned_response.startswith("{"):
                try:
                    import ast
                    parsed = ast.literal_eval(cleaned_response)
                    
                    # Case 1: List of dicts (e.g. [{'type': 'text', 'text': '...'}])
                    if isinstance(parsed, list) and len(parsed) > 0 and isinstance(parsed[0], dict):
                        # Look for 'text' key in the first element
                        if 'text' in parsed[0]:
                            return parsed[0]['text']
                            
                    # Case 2: Single dict (e.g. {'text': '...'})
                    elif isinstance(parsed, dict):
                        if 'text' in parsed:
                            return parsed['text']
                            
                except Exception:
                    # If parsing fails, return original text
                    pass
                 
        return final_response

if __name__ == "__main__":
    bot = Orchestrator()
    print(bot.run("Find a cure for E. Coli resistance"))