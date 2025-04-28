import os
from langchain_community.document_loaders import DirectoryLoader
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain_community.vectorstores import Chroma
from langchain_openai import OpenAIEmbeddings
import shutil
from dotenv import load_dotenv

load_dotenv()

DATA_PATH = r'C:\Users\noliv\phi_generation_repo_destination\phi_generation\code_files\rag_module\data'
CHROMA_PATH = "chroma"


def load_documents():
    loader = DirectoryLoader(DATA_PATH, glob="*.md")
    try:
        documents = loader.load()
        print(f"Loaded {len(documents)} documents.")
        return documents
    except Exception as e:
        print(f"Error loading documents: {e}")
        return []  # Return an empty list to prevent further errors


def split_text(documents):
    if not documents:  # Check if documents is empty
        return []  # Return an empty list if no documents

    text_splitter = RecursiveCharacterTextSplitter(
        chunk_size=2000,
        chunk_overlap=500,
        length_function=len,
        separators=["\n\n", "\n", " ", ""],
        is_separator_regex=False,
    )
    try:
        chunks = text_splitter.split_documents(documents)
        print(f"Split documents into {len(chunks)} chunks.")
        return chunks
    except Exception as e:
        print(f"Error splitting text: {e}")
        return []  # Return an empty list if there was an error


def save_to_chroma(chunks):
    if not chunks:  # Check if chunks is empty
        return

    api_key = os.environ.get("OPENAI_API_KEY")  # get api key from environment.
    embeddings = OpenAIEmbeddings(openai_api_key=api_key)
    try:
        if os.path.exists(CHROMA_PATH):
            shutil.rmtree(CHROMA_PATH)
        vector_store = Chroma.from_documents(
            chunks, embeddings, persist_directory=CHROMA_PATH
        )
        vector_store.persist()
        print(f"Saved {len(chunks)} chunks to {CHROMA_PATH}.")
    except Exception as e:
        print(f"Error saving to Chroma: {e}")


def generate_data_store():
    print("--- Generating Data Store ---")
    documents = load_documents()
    if documents:  # Only proceed if documents were loaded successfully
        chunks = split_text(documents)
        if chunks:  # Only proceed if chunks were split successfully
            save_to_chroma(chunks)
    print("--- Data Store Generation Complete ---")


def main():
    generate_data_store()


if __name__ == "__main__":
    main()


# import os
# from langchain_community.document_loaders import DirectoryLoader
# from langchain.text_splitter import RecursiveCharacterTextSplitter
# from langchain_community.vectorstores import Chroma
# from langchain_openai import OpenAIEmbeddings
# import openai
# from dotenv import load_dotenv
# import shutil

# load_dotenv()

# openai.api_key= ''
# DATA_PATH = ''
# CHROMA_PATH = "chroma"

# def load_documents():
#     loader = DirectoryLoader(DATA_PATH, glob="*.md")
#     documents = loader.load()
#     print(f"\n--- LOAD_DOCUMENTS ---")
#     print(f"Loaded {len(documents)} documents.")
#     for i, doc in enumerate(documents):
#         print(f"\n--- Document {i + 1} ---")
#         print(f"Source: {doc.metadata['source']}")
#         print(f"Content Length: {len(doc.page_content)}")
#         print(f"Content Type: {type(doc.page_content)}")
#         print(f"Content:\n{doc.page_content}\n--- END OF CONTENT ---\n")
#     return documents

# def split_text(documents):
#     text_splitter = RecursiveCharacterTextSplitter(
#         chunk_size=500,
#         chunk_overlap=50,
#         length_function=len,
#         separators=["\n\n", "\n", " ", ""],
#         is_separator_regex=False,
#     )
#     chunks = text_splitter.split_documents(documents)
#     print(f"\n--- SPLIT_TEXT ---")
#     print(f"Split {len(documents)} documents into {len(chunks)} chunks.")
#     for i, chunk in enumerate(chunks):
#         print(f"\n--- Chunk {i + 1} ---")
#         print(f"Chunk Length: {len(chunk.page_content)}")
#         print(f"Chunk Content:\n{chunk.page_content}\n--- END OF CHUNK ---\n")
#     if len(chunks) > 10:
#         document = chunks[10]
#         print(f"\n--- Chunk 11 Details ---")
#         print(f"Content: {document.page_content}")
#         print(f"Metadata: {document.metadata}")
#     else:
#         print("chunks list does not have 11 elements")
#     return chunks

# def save_to_chroma(chunks):
#     embeddings = OpenAIEmbeddings(openai_api_key=openai.api_key)
#     # if there's already stuff at the path chroma_path then delete it 
#     if os.path.exists(CHROMA_PATH):
#         shutil.rmtree(CHROMA_PATH) #delete the chroma db directory
    
#     vector_store = Chroma.from_documents(chunks, embeddings, persist_directory=CHROMA_PATH)
#     vector_store.persist()
#     print(f"\n--- SAVE_TO_CHROMA ---")
#     print(f"Saved {len(chunks)} chunks to {CHROMA_PATH}.")

# def generate_data_store():
#     print("\n--- GENERATE_DATA_STORE ---")
#     documents = load_documents()
#     chunks = split_text(documents)
#     save_to_chroma(chunks)

# def main():
#     print("\n--- MAIN ---")
#     generate_data_store()

# if __name__ == "__main__":
#     main()

# newer version than the one above. newest is top though, removing some print statements.
# import os
# from langchain_community.document_loaders import DirectoryLoader
# from langchain.text_splitter import RecursiveCharacterTextSplitter
# from langchain_community.vectorstores import Chroma
# from langchain_openai import OpenAIEmbeddings
# import shutil
# from dotenv import load_dotenv

# load_dotenv()


# DATA_PATH = r'C:\Users\noliv\phi-generation-repo-destination\phi-generation\code_files\rag_module\data'
# CHROMA_PATH = "chroma"

# def load_documents():
#     loader = DirectoryLoader(DATA_PATH, glob="*.md")
#     documents = loader.load()
#     print(f"\n--- LOAD_DOCUMENTS ---")
#     print(f"Loaded {len(documents)} documents.")
#     for i, doc in enumerate(documents):
#         print(f"\n--- Document {i + 1} ---")
#         print(f"Source: {doc.metadata['source']}")
#         print(f"Content Length: {len(doc.page_content)}")
#         print(f"Content Type: {type(doc.page_content)}")
#         print(f"Content:\n{doc.page_content}\n--- END OF CONTENT ---\n")
#     return documents

# def split_text(documents):
#     text_splitter = RecursiveCharacterTextSplitter(
#         chunk_size=500,
#         chunk_overlap=50,
#         length_function=len,
#         separators=["\n\n", "\n", " ", ""],
#         is_separator_regex=False,
#     )
#     chunks = text_splitter.split_documents(documents)
#     print(f"\n--- SPLIT_TEXT ---")
#     print(f"Split {len(documents)} documents into {len(chunks)} chunks.")
#     for i, chunk in enumerate(chunks):
#         print(f"\n--- Chunk {i + 1} ---")
#         print(f"Chunk Length: {len(chunk.page_content)}")
#         print(f"Chunk Content:\n{chunk.page_content}\n--- END OF CHUNK ---\n")
#     if len(chunks) > 10:
#         document = chunks[10]
#         print(f"\n--- Chunk 11 Details ---")
#         print(f"Content: {document.page_content}")
#         print(f"Metadata: {document.metadata}")
#     else:
#         print("chunks list does not have 11 elements")
#     return chunks

# def save_to_chroma(chunks):
#     api_key = os.environ.get("OPENAI_API_KEY") #get api key from environment.
#     embeddings = OpenAIEmbeddings(openai_api_key=api_key)
#     if os.path.exists(CHROMA_PATH):
#         shutil.rmtree(CHROMA_PATH)
#     vector_store = Chroma.from_documents(chunks, embeddings, persist_directory=CHROMA_PATH)
#     vector_store.persist()
#     print(f"\n--- SAVE_TO_CHROMA ---")
#     print(f"Saved {len(chunks)} chunks to {CHROMA_PATH}.")

# def generate_data_store():
#     print("\n--- GENERATE_DATA_STORE ---")
#     documents = load_documents()
#     chunks = split_text(documents)
#     save_to_chroma(chunks)

# def main():
#     print("\n--- MAIN ---")
#     generate_data_store()

# if __name__ == "__main__":
#     main()