# phi-generation/tests/test_utils.py

import unittest
#import os
import pandas as pd

import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))  # Add parent directory to path
from code_files.utils import csv_to_markdown_table, markdown_table_to_data_dir
#from phi_generation.code_files.utils import csv_to_markdown_table, markdown_table_to_data_dir
class TestUtils(unittest.TestCase):

    def setUp(self):
        # Create a test CSV file
        self.csv_string = """Name,Age,City
Alice,30,New York
Bob,25,London
Charlie,35,Paris"""
        with open("test.csv", "w") as f:
            f.write(self.csv_string)

    def tearDown(self):
        # Clean up test files
        if os.path.exists("test.csv"):
            os.remove("test.csv")
        data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
        test_file = os.path.join(data_dir, "patient_data.md")
        if os.path.exists(test_file):
            os.remove(test_file)
        if os.path.exists(data_dir) and not os.listdir(data_dir): #remove data dir if empty.
            os.rmdir(data_dir)

    def test_csv_to_markdown_table_success(self):
        markdown_table = csv_to_markdown_table("test.csv")
        expected_output = """| Name    |   Age | City     |
|:--------|------:|:---------|
| Alice   |    30 | New York |
| Bob     |    25 | London   |
| Charlie |    35 | Paris    |"""
        self.assertEqual(markdown_table.strip(), expected_output.strip())

    def test_csv_to_markdown_table_file_not_found(self):
        error_message = csv_to_markdown_table("nonexistent.csv")
        self.assertEqual(error_message, "Error: CSV file not found at nonexistent.csv")

    def test_markdown_table_to_data_dir_success(self):
        markdown_table = csv_to_markdown_table("test.csv")
        result = markdown_table_to_data_dir(markdown_table, "patient_data.md")
        self.assertTrue("Markdown table written to" in result)

        data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
        test_file = os.path.join(data_dir, "patient_data.md")
        self.assertTrue(os.path.exists(test_file))
        with open(test_file, "r") as f:
            content = f.read().strip()
            expected_content = """| Name    |   Age | City     |
|:--------|------:|:---------|
| Alice   |    30 | New York |
| Bob     |    25 | London   |
| Charlie |    35 | Paris    |"""
            self.assertEqual(content, expected_content.strip())
    def test_markdown_table_to_data_dir_create_dir_failure(self):
        original_makedirs = os.makedirs
        def mock_makedirs(*args, **kwargs):
            raise OSError("Mock directory creation failure")
        os.makedirs = mock_makedirs
        markdown_table = csv_to_markdown_table("test.csv")
        result = markdown_table_to_data_dir(markdown_table, "patient_data.md")
        self.assertEqual(result, "Error: Could not create data directory: Mock directory creation failure")
        os.makedirs = original_makedirs #reset os.makedirs.