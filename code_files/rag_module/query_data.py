# import argparse
# from langchain_community.vectorstores import Chroma
# from langchain_openai import OpenAIEmbeddings
# from langchain_openai import ChatOpenAI
# from langchain.prompts import ChatPromptTemplate
# import os #add import os

# CHROMA_PATH = "chroma"

# PROMPT_TEMPLATE = """
# Answer the question based only on the following context:

# {context}

# ---

# Answer the question based on the above context: {question}
# """

# def main():
#     print("--- STARTING QUERY_DATA.PY ---")
#     parser = argparse.ArgumentParser()
#     parser.add_argument("query_text", type=str, help="The query text.")
#     args = parser.parse_args()
#     query_text = args.query_text
#     print(f"Query Text: {query_text}")
#     print("Preparing the Chroma DB...")
#     api_key = os.environ.get("OPENAI_API_KEY") #get api key from environment.
#     embedding_function = OpenAIEmbeddings(openai_api_key=api_key) #pass api key.
#     db = Chroma(persist_directory=CHROMA_PATH, embedding_function=embedding_function)
#     print("Chroma DB prepared.")
#     print("Searching the Chroma DB...")
#     results = db.similarity_search_with_relevance_scores(query_text, k=3)
#     print(f"Search Results: {results}")
#     if len(results) == 0 or results[0][1] < 0.7:
#         print(f"Unable to find matching results.")
#         print("--- ENDING QUERY_DATA.PY ---")
#         return
#     context_text = "\n\n---\n\n".join([doc.page_content for doc, _score in results])
#     prompt_template = ChatPromptTemplate.from_template(PROMPT_TEMPLATE)
#     prompt = prompt_template.format(context=context_text, question=query_text)
#     print(f"Prompt: {prompt}")
#     print("Generating response from OpenAI...")
#     model = ChatOpenAI()
#     response_text = model.predict(prompt)
#     print("Response generated.")
#     sources = [doc.metadata.get("source", None) for doc, _score in results]
#     formatted_response = f"Response: {response_text}\nSources: {sources}"
#     print(formatted_response)
#     print("--- ENDING QUERY_DATA.PY ---")

# if __name__ == "__main__":
#     main()




###IMPORTANT: ABOVE IS THE ENTIRELY WORKING CODE WIHTOUT CSV. REVERT BACK TO THIS AT AN ERROR ^ 
import argparse
import os
from langchain_community.vectorstores import Chroma
from langchain_openai import OpenAIEmbeddings, ChatOpenAI
from langchain.prompts import ChatPromptTemplate

CHROMA_PATH = "chroma"

PROMPT_TEMPLATE = """
Answer the question based only on the following context and CSV data (if provided):

{context}

---

{csv_data_description}

---

Answer the question based on the above context and CSV data: {question}
"""

def process_csv_file(csv_file_path):
    print(f"Inside process_csv_file, csv_file_path: {csv_file_path}") #Debug print
    if not csv_file_path:
        print("csv_file_path is empty.") #debug print
        return ""

    try:
        df = pd.read_csv(csv_file_path)
        csv_description = f"CSV Data:\n{df.to_string()}"
        print(f"CSV processed successfully, description: {csv_description}") #debug print
        return csv_description
    except FileNotFoundError:
        error_message = f"Error: CSV file not found at {csv_file_path}"
        print(error_message) #debug print
        return error_message
    except Exception as e:
        error_message = f"Error processing CSV: {e}"
        print(error_message) #debug print
        return error_message

def main():
    print("--- STARTING QUERY_DATA.PY ---")
    parser = argparse.ArgumentParser()
    parser.add_argument("query_text", type=str, help="The query text.")
    parser.add_argument("--csv_file", type=str, help="Path to the CSV file (optional).")
    args = parser.parse_args()
    query_text = args.query_text
    csv_file_path = args.csv_file

    print(f"Query Text: {query_text}")
    print(f"CSV File: {csv_file_path}")

    csv_data_description = process_csv_file(csv_file_path)

    print(f"CSV Data Description: {csv_data_description}") #debug print

    print("Preparing the Chroma DB...")
    api_key = os.environ.get("OPENAI_API_KEY")
    embedding_function = OpenAIEmbeddings(openai_api_key=api_key)
    db = Chroma(persist_directory=CHROMA_PATH, embedding_function=embedding_function)
    print("Chroma DB prepared.")
    print("Searching the Chroma DB...")
    results = db.similarity_search_with_relevance_scores(query_text, k=3)
    print(f"Search Results: {results}")

    if len(results) == 0 or results[0][1] < 0.7:
        print("Unable to find matching results.")
        print("--- ENDING QUERY_DATA.PY ---")
        return

    context_text = "\n\n---\n\n".join([doc.page_content for doc, _score in results])
    prompt_template = ChatPromptTemplate.from_template(PROMPT_TEMPLATE)
    prompt = prompt_template.format(context=context_text, question=query_text, csv_data_description=csv_data_description)
    print(f"Prompt: {prompt}")
    print("Generating response from OpenAI...")
    model = ChatOpenAI()
    response_text = model.predict(prompt)
    print("Response generated.")
    sources = [doc.metadata.get("source", None) for doc, _score in results]
    formatted_response = f"Response: {response_text}\nSources: {sources}"
    print(formatted_response)
    print("--- ENDING QUERY_DATA.PY ---")

if __name__ == "__main__":
    main()