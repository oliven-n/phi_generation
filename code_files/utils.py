# utils.py
# Fix:  make sure you have pandas before running nat it's throwing erros
import pandas as pd
import os

# this function is for taking a local path csv and making it into a .md table for the database

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
    

# this function enables taking the resulting markdown table and throwing it into the data folder
def markdown_table_to_data_dir(markdown_table, filename="data_table.md"):
    """
    Inserts a markdown table into the 'data' directory.

    Args:
        markdown_table: The markdown table string.
        filename: The filename to use for the markdown file.

    Returns:
        A message indicating success or an error message.
    """
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
    if not os.path.exists(data_dir):
        try:
            os.makedirs(data_dir)
        except Exception as e:
            return f"Error: Could not create data directory: {e}"

    file_path = os.path.join(data_dir, filename)

    try:
        with open(file_path, "w") as f:
            f.write(markdown_table)
        return f"Markdown table written to {file_path}"
    except Exception as e:
        return f"Error writing markdown table to file: {e}"
    
    # code_files/utils.py

import shutil
import os

def copy_markdown_to_rag_data(local_file_path):
    """
    Copies a Markdown file from a local path to the rag_module's data directory.

    Args:
        local_file_path (str): Path to the Markdown file to copy.
    """

    rag_data_dir = "code_files/rag_module/data"  # Relative path to the data directory

    try:
        # 1. Ensure the destination directory exists
        os.makedirs(rag_data_dir, exist_ok=True)  # Create the directory if it doesn't exist

        # 2. Extract the filename
        file_name = os.path.basename(local_file_path)

        # 3. Construct the destination path
        destination_path = os.path.join(rag_data_dir, file_name)

        # 4. Copy the file
        shutil.copyfile(local_file_path, destination_path)

        print(f"Successfully copied '{file_name}' to '{rag_data_dir}'")

    except FileNotFoundError:
        print(f"Error: File not found at {local_file_path}")
    except Exception as e:
        print(f"An error occurred: {e}")