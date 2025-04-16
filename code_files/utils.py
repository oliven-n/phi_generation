import pandas as pd
import os
import shutil

def csv_to_markdown_table(csv_file_path):
    """
    Converts a CSV file from a local file path into a markdown table.

    Args:
        csv_file_path: The local file path to the CSV file.

    Returns:
        A markdown table string, or an error message.
    """
    try:
        df = pd.read_csv(csv_file_path)
        markdown_table = df.to_markdown(index=False)
        return markdown_table
    except FileNotFoundError:
        return f"Error: CSV file not found at {csv_file_path}"
    except Exception as e:
        return f"Error converting CSV to markdown: {e}"

def var_markdown_to_data_dir(markdown_content, filename="data_table.md", destination_dir="code_files/data"):
    """
    Writes markdown content to a file in the specified directory, overwriting if the file exists.

    Args:
        markdown_content (str): The markdown content to write.
        filename (str): The filename to use for the markdown file.
        destination_dir (str, optional): The directory to write to.
                                        Defaults to "code_files/data".
                                        Can be a relative or absolute path.

    Returns:
        A message indicating success or an error message.
    """
    # Construct the absolute path to the destination directory
    abs_destination_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", destination_dir)

    if not os.path.exists(abs_destination_dir):
        try:
            os.makedirs(abs_destination_dir, exist_ok=True)
        except Exception as e:
            return f"Error: Could not create destination directory: {e}"

    file_path = os.path.join(abs_destination_dir, filename)

    try:
        with open(file_path, "w") as f:
            f.write(markdown_content)
        return f"Markdown written to {file_path}"
    except Exception as e:
        return f"Error writing markdown to file: {e}"

def copy_file_to_dir(local_file_path, destination_dir, enforce_md=False):
    """
    Copies a file from a local path to the specified destination directory.
    Overwrites the file in the destination directory if a file with the same name exists.
    Optionally enforces that the file being copied is a Markdown (.md) file.

    Args:
        local_file_path (str): The full path to the local file to copy.
        destination_dir (str): The full path to the destination directory.
        enforce_md (bool): If True, only copies files with the '.md' extension.

    Returns:
        A message indicating success or an error message.
    """
    if enforce_md and not local_file_path.lower().endswith(".md"):
        return "Error: Only files with the '.md' extension can be copied."

    if not os.path.exists(destination_dir):
        try:
            os.makedirs(destination_dir, exist_ok=True)
        except Exception as e:
            return f"Error: Could not create destination directory: {e}"

    try:
        file_name = os.path.basename(local_file_path)
        destination_path = os.path.join(destination_dir, file_name)
        shutil.copyfile(local_file_path, destination_path)
        return f"File '{file_name}' copied to '{destination_dir}'"
    except FileNotFoundError:
        return f"Error: File not found at {local_file_path}"
    except Exception as e:
        return f"An error occurred during file copy: {e}"

# Convenience function for copying to the rag_module data directory
def local_markdown_to_data_dir(local_file_path):
    """
    Copies a Markdown file from a local path to the rag_module's data directory.

    Args:
        local_file_path (str): Path to the Markdown file to copy.
    """
    rag_data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code_files", "rag_module", "data")
    result = copy_file_to_dir(local_file_path, rag_data_dir, enforce_md=True)
    print(result)

# Convenience function for adding any document content to the general data directory
def add_document_to_data_dir(document_content, filename):
    """
    Adds document content to a file in the 'data' directory. Overwrites if the file exists.

    Args:
        document_content (str): The content of the document to add.
        filename (str): The name of the file to create.

    Returns:
        A message indicating success or an error message.
    """
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
    file_path = os.path.join(data_dir, filename)
    try:
        with open(file_path, "w") as f:
            f.write(document_content)
        return f"Document written to {file_path}"
    except Exception as e:
        return f"Error writing document to file: {e}"