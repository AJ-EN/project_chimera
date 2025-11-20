import os
from dotenv import load_dotenv
from google import genai

load_dotenv()

# Check which key is available
api_key = os.getenv("GOOGLE_API_KEY") or os.getenv("GEMINI_API_KEY")

if not api_key:
    print("Error: No API key found in environment variables (GOOGLE_API_KEY or GEMINI_API_KEY).")
else:
    print(f"Found API Key: {api_key[:5]}...{api_key[-5:]}")
    
    try:
        # Initialize client
        client = genai.Client(api_key=api_key)
        
        print("Attempting to generate content with 'gemini-2.5-flash'...")
        response = client.models.generate_content(
            model="gemini-2.5-flash", 
            contents="Hello, are you working?"
        )
        print("\nSuccess! Response:")
        print(response.text)
    except Exception as e:
        print(f"\nError: {e}")
