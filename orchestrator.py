import os
from dotenv import load_dotenv
from langchain_google_genai import ChatGoogleGenerativeAI
from langchain.agents import initialize_agent, Tool, AgentType
from agents.librarian import Librarian
from agents.lab_rat import LabRat
from agents.web_researcher import WebResearcher

class Orchestrator:
    def __init__(self):
        self.llm = ChatGoogleGenerativeAI(
            model="gemini-2.5-flash",  # Use Flash for speed, Pro for reasoning
            temperature=0.2,
            convert_system_message_to_human=True
        )
        
        # Initialize Agents
        self.librarian = Librarian()
        self.lab_rat = LabRat()
        self.web_researcher = WebResearcher()

        # Define Tools
        self.tools = [
            Tool(
                name="Librarian_Search",
                func=self.librarian.search,
                description="Useful for finding specific biological targets and PDB IDs from the internal PDF library. Use this for specific drug discovery tasks."
            ),
            Tool(
                name="Lab_Rat_Simulation",
                func=self.lab_rat.run_simulation,
                description="Useful for screening drugs and calculating binding affinity scores. Input should be the molecule name."
            ),
            Tool(
                name="Web_Search",
                func=self.web_researcher.search,
                description="Useful for answering general questions, finding the latest research, or when the Librarian doesn't have the info. Use this for 'What is...', 'Latest research on...', or casual chat."
            )
        ]

        # Initialize the Agent
        self.agent = initialize_agent(
            self.tools,
            self.llm,
            agent=AgentType.ZERO_SHOT_REACT_DESCRIPTION,
            # This prints the "Thought Process" to the console (Crucial for debugging)
            verbose=True,
            handle_parsing_errors=True,
            max_iterations=5,
            early_stopping_method="generate"
        )


    def run(self, user_query):
        system_prompt = """
        You are Dr. Chimera, the Principal Investigator of an autonomous drug discovery lab.
        
        You have three modes of operation:
        1. **General Chat / Latest Research**: Use the 'Web_Search' tool.
           - Example: "Hi", "What can you do?", "Latest research on Alzheimer's".
           
        2. **Internal Knowledge Retrieval**: Use 'Librarian_Search'.
           - Use this when asked to find targets from your internal files.
           
        3. **Drug Discovery Protocol** (Strict):
           - If asked to "Find a cure" or "Screen drugs":
             1. Use 'Librarian_Search' to find a target (PDB ID).
             2. Use 'Lab_Rat_Simulation' to screen drugs.
             3. Synthesize results.
        
        CRITICAL INSTRUCTION: 
        Once you have a simulation result OR a web search answer, you MUST immediately provide a Final Answer. 
        Do not keep searching.
        
        Format:
        Final Answer: [Your recommendation here]
        """

        # We prepend the system prompt to the user query for the Zero Shot Agent
        full_prompt = f"{system_prompt}\n\nUser Query: {user_query}"
        return self.agent.run(full_prompt)


# Test block
if __name__ == "__main__":
    bot = Orchestrator()
    print(bot.run("Find a cure for E. Coli resistance"))
