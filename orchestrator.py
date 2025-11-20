import os
from dotenv import load_dotenv
from langchain_google_genai import ChatGoogleGenerativeAI
from langchain.agents import initialize_agent, AgentType
from langchain.tools import Tool

# Load environment variables
load_dotenv()


from agents.librarian import Librarian
from agents.lab_rat import LabRat

class Orchestrator:
    def __init__(self):
        self.llm = ChatGoogleGenerativeAI(
            model="gemini-2.5-flash",  # Use Flash for speed, Pro for reasoning
            temperature=0.2,
            convert_system_message_to_human=True
        )
        
        self.librarian = Librarian()
        self.lab_rat = LabRat()

        # We will add real tools later. For now, these are placeholders.
        self.tools = [
            Tool(
                name="Librarian_Search",
                func=self.librarian.search,
                description="Useful for finding biological targets and PDB IDs from scientific literature."
            ),
            Tool(
                name="Lab_Rat_Simulation",
                func=self.lab_rat.run_simulation,
                description="Useful for screening drugs and calculating binding affinity scores."
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
        Your goal is to find antibiotic candidates.
        
        Protocol:
        1. ALWAYS use the 'Librarian_Search' tool first to identify the biological target (PDB ID).
        2. Once you have a target, use the 'Lab_Rat_Simulation' tool to screen potential drugs.
        3. Synthesize the results into a final recommendation.
        
        CRITICAL INSTRUCTION: 
        Once you have a simulation result, you MUST immediately provide a Final Answer. 
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
