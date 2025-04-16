# code_files/st_to_unst_module/st_to_unst.py

import argparse
import csv
from langchain.llms import OpenAI
from langchain.prompts import PromptTemplate

DEFAULT_PROMPT = """You are a doctor. Create an extensive, non-list public health report of the patient, using the csv data, in paragraph format. This should contain an introduction of 200 words, then additional 100 word paragraphs for each "relevant" illness and medication listed in the table, where "relevant" information is represented by a "1". Begin with listing the Patient's name, age, gender, and assigned number. The conditions and medications should not be separated from one another."""

def generate_patient_report(csv_file, custom_prompt=None):
    """
    Generates a doctor's note-style patient health record from CSV data using an LLM.

    Args:
        csv_file (str): Path to the CSV file containing patient data.
                        The first row should be headers (patient identifiers and traits).
                        Subsequent rows represent patients with 1s (present) and 0s (absent) for traits.
        custom_prompt (str, optional): A custom prompt for the LLM. Defaults to None.

    Returns:
        str: The generated patient health record.
    """
    try:
        with open(csv_file, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            patient_data = list(reader)

            if not patient_data:
                return "Error: No patient data found in the CSV file."

            # Assuming only one patient per CSV for simplicity
            patient = patient_data[0]

            patient_info = f"Patient Name: {patient.get('Name', 'N/A')}, "
            patient_info += f"Age: {patient.get('Age', 'N/A')}, "
            patient_info += f"Gender: {patient.get('Gender', 'N/A')}, "
            patient_info += f"Assigned Number: {patient.get('Number', 'N/A')}\n\n"

            relevant_conditions_medications = []
            for trait, value in patient.items():
                if value == '1' and trait not in ['Name', 'Age', 'Gender', 'Number']:
                    relevant_conditions_medications.append(trait)

            if custom_prompt:
                prompt_template = PromptTemplate.from_template(custom_prompt)
            else:
                prompt_template = PromptTemplate.from_template(DEFAULT_PROMPT)

            llm = OpenAI(temperature=0.7)  # You can adjust temperature and other LLM parameters

            prompt_input = {
                "patient_info": patient_info,
                "relevant_traits": ", ".join(relevant_conditions_medications)
            }

            # Dynamically construct the prompt based on relevant traits
            if relevant_conditions_medications:
                introduction_prompt = prompt_template.format(
                    patient_info=prompt_input["patient_info"],
                    relevant_traits=prompt_input["relevant_traits"]
                )
                report = llm.predict(introduction_prompt)

                for trait in relevant_conditions_medications:
                    specific_prompt_template = PromptTemplate.from_template(f"Describe the condition/medication '{trait}' in a paragraph of approximately 100 words in the context of a patient health record.")
                    specific_prompt = specific_prompt_template.format(trait=trait)
                    report += "\n\n" + llm.predict(specific_prompt)
            else:
                # If no relevant traits, just generate a basic introduction
                basic_intro_template = PromptTemplate.from_template("Generate a 200-word introduction for a patient health record with the following basic information: {patient_info}")
                basic_intro_prompt = basic_intro_template.format(patient_info=prompt_input["patient_info"])
                report = llm.predict(basic_intro_prompt)

            return report

    except FileNotFoundError:
        return f"Error: CSV file not found at {csv_file}"
    except Exception as e:
        return f"An error occurred: {e}"

def main():
    parser = argparse.ArgumentParser(description="Generate a patient health record from a CSV file using an LLM.")
    parser.add_argument("csv_file", type=str, help="Path to the CSV file containing patient data.")
    parser.add_argument("--custom_prompt", type=str, help="Optional custom prompt for the LLM.", default=None)
    args = parser.parse_args()

    report = generate_patient_report(args.csv_file, args.custom_prompt)
    print("\n--- Patient Health Record ---")
    print(report)

if __name__ == "__main__":
    main()