# Super drastic changes impending

import argparse
import csv
import os
from langchain.llms import OpenAI
from langchain.prompts import PromptTemplate
import anndata as ad
import json

DEFAULT_PROMPT = """You are a doctor. Create an extensive, non-list public health report of the patient, using the csv data, in paragraph format. This should contain an introduction of 200 words, then additional 100 word paragraphs for each "relevant" illness and medication listed in the table, where "relevant" information is represented by a "1". Begin with listing the Patient's name, age, gender, and assigned number. The conditions and medications should not be separated from one another."""

def generate_patient_report(patient_data, trait_columns, custom_prompt=None):
    """
    Generates a doctor's note-style patient health record from structured patient data.

    Args:
        patient_data (dict): A dictionary representing a single patient's data.
                             Keys are column headers from the CSV, values are patient information.
        trait_columns (list): A list of column headers that represent traits (conditions, medications).
                              These columns are expected to contain binary values (1/0).
        custom_prompt (str, optional): A custom prompt for the LLM. Defaults to None.

    Returns:
        str: The generated patient health record.
    """
    try:
        patient_info = ""
        for key, value in patient_data.items():
            if key not in trait_columns:
                patient_info += f"{key}: {value}, "
        patient_info = patient_info.rstrip(', ') + "\n\n"

        relevant_traits = [trait for trait in trait_columns if patient_data.get(trait) == '1']

        if custom_prompt:
            prompt_template = PromptTemplate.from_template(custom_prompt)
        else:
            prompt_template = PromptTemplate.from_template(DEFAULT_PROMPT)

        llm = OpenAI(temperature=0.7)  # You can adjust temperature

        if relevant_traits:
            introduction_prompt = prompt_template.format(
                patient_info=patient_info,
                relevant_traits=", ".join(relevant_traits)
            )
            report = llm.predict(introduction_prompt)

            for trait in relevant_traits:
                specific_prompt_template = PromptTemplate.from_template(
                    f"Describe the condition/medication '{trait}' in a paragraph of approximately 100 words in the context of a patient health record."
                )
                specific_prompt = specific_prompt_template.format(trait=trait)
                report += "\n\n" + llm.predict(specific_prompt)
        else:
            basic_intro_template = PromptTemplate.from_template(
                "Generate a 200-word introduction for a patient health record with the following basic information: {patient_info}"
            )
            basic_intro_prompt = basic_intro_template.format(patient_info=patient_info)
            report = llm.predict(basic_intro_prompt)

        return report

    except Exception as e:
        return f"An error occurred: {e}"

def create_anndata_from_csv(csv_file, patient_id_columns=None):
    """
    Processes a CSV file, generates patient reports for each row, and creates an AnnData object.

    Args:
        csv_file (str): Path to the CSV file containing patient data.
        patient_id_columns (list, optional): List of column names to use as patient identifiers.
                                             Defaults to ['Name'].

    Returns:
        anndata.AnnData: An AnnData object containing patient reports and metadata.
    """
    try:
        with open(csv_file, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            patient_data = list(reader)

        if not patient_data:
            raise ValueError("No patient data found in the CSV file.")

        # Determine trait columns (columns with binary data)
        all_columns = patient_data[0].keys()
        if patient_id_columns is None:
            patient_id_columns = ['Name']
        trait_columns = [col for col in all_columns if col not in patient_id_columns]

        report_data = []
        patient_names = []
        for patient in patient_data:
            report = generate_patient_report(patient, trait_columns)
            report_data.append(report)
            patient_name = ", ".join([patient.get(col, 'N/A') for col in patient_id_columns])
            patient_names.append(patient_name)

        # Create AnnData object
        try:  # Nested try block for pandas and anndata imports
            import pandas as pd
            import anndata as ad
            adata = ad.AnnData(pd.DataFrame({'patient_name': patient_names, 'unstructured_phi': report_data}))
            adata.obs = pd.DataFrame(patient_data)
            adata.uns['csv_filename'] = os.path.basename(csv_file)  # Store CSV filename in metadata
            return adata
        except ImportError as e:
            raise ImportError(f"Error importing pandas or anndata: {e}")

    except FileNotFoundError:
        raise FileNotFoundError(f"CSV file not found at {csv_file}")
    except Exception as e:
        raise Exception(f"An error occurred during AnnData creation: {e}")
# def create_anndata_from_csv(csv_file, patient_id_columns=None):
#     """
#     Processes a CSV file, generates patient reports for each row, and creates an AnnData object.

#     Args:
#         csv_file (str): Path to the CSV file containing patient data.
#         patient_id_columns (list, optional): List of column names to use as patient identifiers.
#                                              These columns will be included in adata.X.
#                                              Defaults to ['Name'].

#     Returns:
#         anndata.AnnData: An AnnData object containing patient reports and metadata.
#     """
#     try:
#         with open(csv_file, 'r', newline='') as csvfile:
#             reader = csv.DictReader(csvfile)
#             patient_data = list(reader)

#         if not patient_data:
#             raise ValueError("No patient data found in the CSV file.")

#         # Determine trait columns (columns with binary data)
#         all_columns = patient_data[0].keys()
#         if patient_id_columns is None:
#             patient_id_columns = ['Name']
#         trait_columns = [col for col in all_columns if col not in patient_id_columns]

#         report_data = []
#         patient_names = []
#         for patient in patient_data:
#             report = generate_patient_report(patient, trait_columns)
#             report_data.append(report)
#             patient_name = ", ".join([patient.get(col, 'N/A') for col in patient_id_columns])
#             patient_names.append(patient_name)

#         # Create AnnData object
#         adata = ad.AnnData(pd.DataFrame({'patient_name': patient_names, 'unstructured_phi': report_data}))
#         adata.obs = pd.DataFrame(patient_data)
#         adata.uns['csv_filename'] = os.path.basename(csv_file)  # Store CSV filename in metadata

#         return adata

#     except FileNotFoundError:
#         raise FileNotFoundError(f"CSV file not found at {csv_file}")
#     except Exception as e:
#         raise Exception(f"An error occurred during AnnData creation: {e}")

# def create_anndata_from_csv(csv_file, patient_id_columns=None):
#     """
#     Processes a CSV file, generates patient reports, and creates an AnnData object.

#     Args:
#         csv_file (str): Path to the CSV file containing patient data.
#         patient_id_columns (list, optional): List of column names to use as patient identifiers.
#                                              These columns will be included in adata.X.
#                                              Defaults to ['Name'].

#     Returns:
#         anndata.AnnData: An AnnData object containing patient reports and metadata.
#     """
#     try:
#         with open(csv_file, 'r', newline='') as csvfile:
#             reader = csv.DictReader(csvfile)
#             patient_data = list(reader)

#         if not patient_data:
#             raise ValueError("No patient data found in the CSV file.")

#         # Determine trait columns (columns with binary data)
#         all_columns = patient_data[0].keys()
#         if patient_id_columns is None:
#             patient_id_columns = ['Name']
#         trait_columns = [col for col in all_columns if col not in patient_id_columns]

#         report_data = []
#         patient_names = []
#         for patient in patient_data:
#             report = generate_patient_report(patient, trait_columns)
#             report_data.append(report)
#             patient_name = ", ".join([patient.get(col, 'N/A') for col in patient_id_columns])
#             patient_names.append(patient_name)

#         # Create AnnData object
#         adata_X = pd.DataFrame({'patient_name': patient_names, 'unstructured_phi': report_data})
#         adata = ad.AnnData(adata_X)
#         adata.obs = pd.DataFrame(patient_data)
#         adata.uns['csv_filename'] = os.path.basename(csv_file)  # Store CSV filename in metadata

#         return adata

#     except FileNotFoundError:
#         raise FileNotFoundError(f"CSV file not found at {csv_file}")
#     except Exception as e:
#         raise Exception(f"An error occurred during AnnData creation: {e}")

def anndata_to_json(adata, filename="patient_reports.json"):
    """
    Converts the unstructured data from an AnnData object to a JSON file.

    Args:
        adata (anndata.AnnData): The AnnData object containing patient reports.
        filename (str, optional): The name of the JSON file to create.
                                 Defaults to "patient_reports.json".
    """
    try:
        json_data = {}
        for index, row in adata.X.iterrows():
            json_data[row['patient_name']] = row['unstructured_phi']

        with open(filename, 'w') as f:
            json.dump(json_data, f, indent=4)  # Use indent for readability

        if 'csv_filename' in adata.uns:
            json_metadata = {'csv_filename': adata.uns['csv_filename']}
            with open(filename.replace(".json", "_metadata.json"), 'w') as f:
              json.dump(json_metadata, f, indent=4)

        print(f"Patient reports written to {filename}")

    except Exception as e:
        print(f"An error occurred during JSON creation: {e}")

def main():
    parser = argparse.ArgumentParser(description="Generate patient health records and perform related operations.")
    parser.add_argument("csv_file", type=str, help="Path to the CSV file containing patient data.")
    parser.add_argument("--custom_prompt", type=str, help="Optional custom prompt for the LLM.", default=None)
    parser.add_argument("--output_json", action='store_true', help="Output reports to a JSON file.")
    args = parser.parse_args()

    try:
        adata = create_anndata_from_csv(args.csv_file)

        # Print a sample report
        print("\n--- Sample Patient Report ---")
        print(adata.X.iloc[0]['unstructured_phi'])

        if args.output_json:
            anndata_to_json(adata, filename="patient_reports.json")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()



# # code_files/st_to_unst_module/st_to_unst.py
# # This code file is VERY similar to the code for the rag_module, except is a vanilla llm instead of a RAG
# import argparse
# import csv
# from langchain.llms import OpenAI
# from langchain.prompts import PromptTemplate

# DEFAULT_PROMPT = """You are a doctor. Create an extensive, non-list public health report of the patient, using the csv data, in paragraph format. This should contain an introduction of 200 words, then additional 100 word paragraphs for each "relevant" illness and medication listed in the table, where "relevant" information is represented by a "1". Begin with listing the Patient's name, age, gender, and assigned number. The conditions and medications should not be separated from one another."""

# def generate_patient_report(csv_file, custom_prompt=None):
#     """
#     Generates a doctor's note-style patient health record from CSV data using an LLM.

#     Args:
#         csv_file (str): Path to the CSV file containing patient data.
#                         The first row should be headers (patient identifiers and traits).
#                         Subsequent rows represent patients with 1s (present) and 0s (absent) for traits.
#         custom_prompt (str, optional): A custom prompt for the LLM. Defaults to None.

#     Returns:
#         str: The generated patient health record.
#     """
#     try:
#         with open(csv_file, 'r', newline='') as csvfile:
#             reader = csv.DictReader(csvfile)
#             patient_data = list(reader)

#             if not patient_data:
#                 return "Error: No patient data found in the CSV file."

#             # Assuming only one patient per CSV for simplicity
#             patient = patient_data[0]

#             patient_info = f"Patient Name: {patient.get('Name', 'N/A')}, "
#             patient_info += f"Age: {patient.get('Age', 'N/A')}, "
#             patient_info += f"Gender: {patient.get('Gender', 'N/A')}, "
#             patient_info += f"Assigned Number: {patient.get('Number', 'N/A')}\n\n"

#             relevant_conditions_medications = []
#             for trait, value in patient.items():
#                 if value == '1' and trait not in ['Name', 'Age', 'Gender', 'Number']:
#                     relevant_conditions_medications.append(trait)

#             if custom_prompt:
#                 prompt_template = PromptTemplate.from_template(custom_prompt)
#             else:
#                 prompt_template = PromptTemplate.from_template(DEFAULT_PROMPT)

#             llm = OpenAI(temperature=0.7)  # You can adjust temperature and other LLM parameters

#             prompt_input = {
#                 "patient_info": patient_info,
#                 "relevant_traits": ", ".join(relevant_conditions_medications)
#             }

#             # Dynamically construct the prompt based on relevant traits
#             if relevant_conditions_medications:
#                 introduction_prompt = prompt_template.format(
#                     patient_info=prompt_input["patient_info"],
#                     relevant_traits=prompt_input["relevant_traits"]
#                 )
#                 report = llm.predict(introduction_prompt)

#                 for trait in relevant_conditions_medications:
#                     specific_prompt_template = PromptTemplate.from_template(f"Describe the condition/medication '{trait}' in a paragraph of approximately 100 words in the context of a patient health record.")
#                     specific_prompt = specific_prompt_template.format(trait=trait)
#                     report += "\n\n" + llm.predict(specific_prompt)
#             else:
#                 # If no relevant traits, just generate a basic introduction
#                 basic_intro_template = PromptTemplate.from_template("Generate a 200-word introduction for a patient health record with the following basic information: {patient_info}")
#                 basic_intro_prompt = basic_intro_template.format(patient_info=prompt_input["patient_info"])
#                 report = llm.predict(basic_intro_prompt)

#             return report

#     except FileNotFoundError:
#         return f"Error: CSV file not found at {csv_file}"
#     except Exception as e:
#         return f"An error occurred: {e}"

# def main():
#     parser = argparse.ArgumentParser(description="Generate a patient health record from a CSV file using an LLM.")
#     parser.add_argument("csv_file", type=str, help="Path to the CSV file containing patient data.")
#     parser.add_argument("--custom_prompt", type=str, help="Optional custom prompt for the LLM.", default=None)
#     args = parser.parse_args()

#     report = generate_patient_report(args.csv_file, args.custom_prompt)
#     print("\n--- Patient Health Record ---")
#     print(report)

# if __name__ == "__main__":
#     main()