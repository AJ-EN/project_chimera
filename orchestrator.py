import os
from dotenv import load_dotenv
from langchain_google_genai import ChatGoogleGenerativeAI
from langgraph.prebuilt import create_react_agent # The Modern Way
from langchain_core.messages import SystemMessage

# Import your custom agents
from agents.librarian import Librarian
from agents.lab_rat import LabRat
from agents.web_researcher import WebResearcher
from langchain.tools import Tool

load_dotenv()

class Orchestrator:
    def __init__(self):
        # 1. Initialize the Brain (Gemini 2.5 Flash)
        self.llm = ChatGoogleGenerativeAI(
            model="gemini-2.5-flash",
            temperature=0.2,
            convert_system_message_to_human=True
        )
        
        # 2. Initialize your Specialist Agents
        self.librarian = Librarian()
        self.lab_rat = LabRat()
        self.web_researcher = WebResearcher()

        # 3. Define Tools (The "Hands" of the Agent)
        self.tools = [
            Tool(
                name="Librarian_Search",
                func=self.librarian.search,
                description="Use this to find biological targets (e.g., 'Find the target for E. Coli')."
            ),
            Tool(
                name="Lab_Rat_Simulation",
                func=self.lab_rat.run_simulation,
                description="Use this to screen drugs (e.g., 'Simulate Ciprofloxacin'). Input: Molecule Name."
            ),
            Tool(
                name="Web_Search",
                func=self.web_researcher.search,
                description="Use this for general questions or latest news (e.g., 'What is the flu?')."
            )
        ]

        # 4. Create the Graph (The "New" Way)
        # This automatically creates a 'Supervisor' workflow that the course teaches
        self.agent_graph = create_react_agent(self.llm, self.tools)

    def run(self, user_query):
        # 5. The System Prompt (The "Personality")
        system_instruction = """You are Dr. Chimera, a Nobel-prize winning computational biologist.
        Protocol:
        1. If the user asks for a cure/drug, FIRST use 'Librarian_Search' to find the biological target (PDB ID).
        2. THEN use 'Lab_Rat_Simulation' to screen candidate molecules against that target.
        3. ALWAYS explain the science behind your decision.
        4. If the user just says 'Hi' or asks general questions, use 'Web_Search'.
        """
        
        # 6. Run the Graph
        inputs = {"messages": [("system", system_instruction), ("user", user_query)]}
        
        # We stream the output to get the final response
        # In a real app, you might want to stream intermediate steps to the UI
        result = self.agent_graph.invoke(inputs)
        
        # Extract the final message content
        return result["messages"][-1].content

# Test block
if __name__ == "__main__":
    bot = Orchestrator()
    print(bot.run("Find a cure for E. Coli resistance"))
