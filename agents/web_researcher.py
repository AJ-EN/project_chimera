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
            
            # Extract the text - handle different response formats
            result = response.text
            
            # parsing logic for when result is a string representation of a list/dict
            try:
                import ast
                # If it looks like a list or dict string, try to parse it
                if isinstance(result, str) and (result.strip().startswith('[') or result.strip().startswith('{')):
                    parsed = ast.literal_eval(result)
                    if isinstance(parsed, list) and len(parsed) > 0 and isinstance(parsed[0], dict):
                        return parsed[0].get('text', str(parsed))
                    elif isinstance(parsed, dict):
                        return parsed.get('text', str(parsed))
            except:
                # If parsing fails, just continue and return the original string
                pass

            # Sometimes the API returns a list of dicts (actual object), sometimes just text
            if isinstance(result, list):
                # Extract text from the first item if it's a list
                if result and isinstance(result[0], dict) and 'text' in result[0]:
                    return result[0]['text']
                return str(result)
            elif isinstance(result, dict):
                # Extract text if it's a dict
                return result.get('text', str(result))
            
            # Otherwise return as-is (should be a string)
            return result
            
        except Exception as e:
            return f"Web Search Error: {str(e)}"
