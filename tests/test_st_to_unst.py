# phi-generation/tests/test_st_to_unst.py

import unittest
import os
import sys
import pandas as pd
import numpy as np
import anndata as ad

# Ensure parent path is included
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# Correct import
from code_files.st_to_unst_module import st_to_unst as sun

class TestStToUnst(unittest.TestCase):

    def setUp(self):
        # Create a small test CSV file for patients
        self.test_csv_path = "test_patient_data.csv"
        self.csv_content = """Patient Name,Patient Age,Diabetes,Hypertension
Alice,30,1,0
Bob,45,0,1
Charlie,29,1,1"""
        with open(self.test_csv_path, "w") as f:
            f.write(self.csv_content)

    def tearDown(self):
        # Clean up test file
        if os.path.exists(self.test_csv_path):
            os.remove(self.test_csv_path)

    def test_obs_matches_csv(self):
        # Create AnnData
        adata = sun.process_csv_and_create_anndata(self.test_csv_path, patient_id_columns=["Patient Name"])

        # Load original CSV
        df_original = pd.read_csv(self.test_csv_path)

        # Check that .obs matches exactly
        pd.testing.assert_frame_equal(adata.obs.reset_index(drop=True), df_original.reset_index(drop=True))

    def test_patient_name_length_limit(self):
        # Create AnnData
        adata = sun.process_csv_and_create_anndata(self.test_csv_path, patient_id_columns=["Patient Name"])

        # Check that patient names (first column of adata.X) are less than 20 characters
        patient_names = adata.X[:, 0]
        for name in patient_names:
            self.assertTrue(len(str(name)) <= 20, f"Patient name too long: {name}")

    def test_unstructured_note_not_empty(self):
        # Ensure unstructured doctor notes are generated and not empty
        adata = sun.process_csv_and_create_anndata(self.test_csv_path, patient_id_columns=["Patient Name"])

        doctor_notes = adata.X[:, 1]
        for note in doctor_notes:
            self.assertTrue(isinstance(note, str))
            self.assertGreater(len(note.strip()), 0, "Generated doctor's note is empty!")

    def test_ann_data_var_names_correct(self):
        # Ensure var_names are set correctly
        adata = sun.process_csv_and_create_anndata(self.test_csv_path, patient_id_columns=["Patient Name"])

        self.assertEqual(list(adata.var_names), ["patient_name", "unstructured_note"])

    def test_uns_metadata_filename_set(self):
        # Check that the original CSV filename is stored in .uns
        adata = sun.process_csv_and_create_anndata(self.test_csv_path, patient_id_columns=["Patient Name"])
        
        self.assertIn('csv_filename', adata.uns)
        self.assertEqual(adata.uns['csv_filename'], os.path.basename(self.test_csv_path))
