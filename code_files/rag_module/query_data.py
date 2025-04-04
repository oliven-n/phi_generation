import argparse
from langchain_community.vectorstores import Chroma
from langchain_openai import OpenAIEmbeddings
from langchain_openai import ChatOpenAI
from langchain.prompts import ChatPromptTemplate
import os #add import os

CHROMA_PATH = "chroma"

PROMPT_TEMPLATE = """
Answer the question based only on the following context:

{context}

---

Answer the question based on the above context: {question}
"""

def main():
    print("--- STARTING QUERY_DATA.PY ---")
    parser = argparse.ArgumentParser()
    parser.add_argument("query_text", type=str, help="The query text.")
    args = parser.parse_args()
    query_text = args.query_text
    print(f"Query Text: {query_text}")
    print("Preparing the Chroma DB...")
    api_key = os.environ.get("OPENAI_API_KEY") #get api key from environment.
    embedding_function = OpenAIEmbeddings(openai_api_key=api_key) #pass api key.
    db = Chroma(persist_directory=CHROMA_PATH, embedding_function=embedding_function)
    print("Chroma DB prepared.")
    print("Searching the Chroma DB...")
    results = db.similarity_search_with_relevance_scores(query_text, k=3)
    print(f"Search Results: {results}")
    if len(results) == 0 or results[0][1] < 0.7:
        print(f"Unable to find matching results.")
        print("--- ENDING QUERY_DATA.PY ---")
        return
    context_text = "\n\n---\n\n".join([doc.page_content for doc, _score in results])
    prompt_template = ChatPromptTemplate.from_template(PROMPT_TEMPLATE)
    prompt = prompt_template.format(context=context_text, question=query_text)
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