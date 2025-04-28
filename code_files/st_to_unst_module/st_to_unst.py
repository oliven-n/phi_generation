
# gpt refactored to make it modular (loop over each patient and update accordingly) and then i debugged. Seems to work now!!

# code_files/st_to_unst_module/st_to_unst.py

import os
import numpy as np
import pandas as pd
import anndata as ad
from dotenv import load_dotenv
from langchain_openai import OpenAI  # UPDATED: new style import
from langchain.prompts import PromptTemplate

load_dotenv()

DEFAULT_PROMPT = """You are a doctor. Create an extensive, non-list public health report of the patient, using the csv data, in paragraph format. This should contain an introduction of 200 words, then additional 100 word paragraphs for each "relevant" illness and medication listed in the table, where "relevant" information is represented by a "1". Begin with listing the Patient's name, age, gender, and assigned number. The conditions and medications should not be separated from one another."""

def generate_doctors_note(patient_row, trait_columns, custom_prompt=None):
    """
    Generate a doctor's note using the patient information.
    """
    llm = OpenAI(temperature=0.7)

    patient_info = ", ".join([f"{col}: {patient_row[col]}" for col in patient_row.index if col not in trait_columns])
    relevant_traits = [col for col in trait_columns if str(patient_row[col]) == '1']

    if custom_prompt:
        prompt_template = PromptTemplate.from_template(custom_prompt)
    else:
        prompt_template = PromptTemplate.from_template(DEFAULT_PROMPT)

    prompt_text = prompt_template.format(patient_info=patient_info, relevant_traits=", ".join(relevant_traits))

    # UPDATED: use .invoke instead of .predict
    doctors_note = llm.invoke(prompt_text)
    return doctors_note

def create_initial_dataframe():
    """
    Creates an empty pandas DataFrame for storing patient name and doctor's note.
    """
    df = pd.DataFrame(columns=["patient_name", "unstructured_note"])
    return df

def update_dataframe(df, patient_name, doctors_note):
    """
    Appends a new patient entry into the DataFrame.
    """
    new_entry = {"patient_name": patient_name, "unstructured_note": doctors_note}
    df = pd.concat([df, pd.DataFrame([new_entry])], ignore_index=True)
    return df

def process_csv_and_create_anndata(csv_file_path, patient_id_columns=["Patient Name"], custom_prompt=None):
    """
    Processes the input CSV, generates doctor's notes, and builds AnnData.
    """
    df_csv = pd.read_csv(csv_file_path)

    if df_csv.empty:
        raise ValueError("Input CSV is empty.")

    trait_columns = [col for col in df_csv.columns if col not in patient_id_columns]

    # Step 1: Create empty working DataFrame
    working_df = create_initial_dataframe()

    # Step 2: Generate doctor's note per patient
    for idx, patient_row in df_csv.iterrows():
        patient_name = " ".join([str(patient_row[col]) for col in patient_id_columns])
        doctors_note = generate_doctors_note(patient_row, trait_columns, custom_prompt=custom_prompt)
        working_df = update_dataframe(working_df, patient_name, doctors_note)

    # Step 3: Create AnnData object properly
    adata = ad.AnnData(X=working_df.values)
    adata.obs = df_csv.copy()
    adata.var_names = working_df.columns.tolist()
    adata.uns['csv_filename'] = os.path.basename(csv_file_path)

    return adata
