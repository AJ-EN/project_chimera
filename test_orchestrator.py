from orchestrator import Orchestrator
import sys

# Suppress warnings for cleaner output
import warnings
warnings.filterwarnings("ignore")

def test_orchestrator():
    print("Initializing Orchestrator...")
    try:
        bot = Orchestrator()
        print("Orchestrator initialized.")
        
        query = "Find a cure for E. Coli resistance"
        print(f"\nRunning query: '{query}'")
        
        # We'll run it and print the result
        # The verbose=True in orchestrator will show the thought process
        response = bot.run(query)
        
        print("\nFinal Response:")
        print(response)
        
    except Exception as e:
        print(f"\nError running Orchestrator: {e}")

if __name__ == "__main__":
    test_orchestrator()
