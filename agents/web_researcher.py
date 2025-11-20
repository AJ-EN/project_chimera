import os
from dotenv import load_dotenv
from google import genai
from google.genai import types

class WebResearcher:
    def __init__(self):
        load_dotenv()
        api_key = os.getenv("GOOGLE_API_KEY")
        if not api_key:
            raise ValueError("GOOGLE_API_KEY not found in environment variables")
        
        self.client = genai.Client(api_key=api_key)
        self.model = "gemini-2.5-flash"

    def search(self, query):
        """
        Performs a Google Search using Gemini Grounding to answer the query.
        """
        try:
            grounding_tool = types.Tool(
                google_search=types.GoogleSearch()
            )

            config = types.GenerateContentConfig(
                tools=[grounding_tool],
                system_instruction="You are a specialized medical research assistant. You ONLY answer questions related to medicine, biology, chemistry, drug discovery, and scientific research. If a user asks about unrelated topics (e.g., sports, politics, entertainment, general trivia), politely refuse and state that your purpose is strictly restricted to biomedical research."
            )

            response = self.client.models.generate_content(
                model=self.model,
                contents=query,
                config=config,
            )
            
            # Extract the text and any grounding metadata if needed
            # For now, we just return the text which should be grounded
            return response.text
            
        except Exception as e:
            return f"Web Search Error: {str(e)}"
