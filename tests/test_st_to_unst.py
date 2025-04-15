# tests/st_to_unst_tests.py

import unittest
from unittest.mock import patch, mock_open
from io import StringIO
import sys
import os
# it is better to import the py file and use the 
from code_files.st_to_unst_module import st_to_unst as sun

class TestSTToUNSTModule(unittest.TestCase):
    def test_generate_patient_report(self):
        # create a json with this information : ideas
        expected_type= str
        result = sun.generate_patient_report("dummy.csv")
        self.assertIsInstance(result,expected_type, "LLM output is not string.")

#     # Tests the generate_patient_report function when the patient data in the CSV
#     # contains relevant traits (value '1'). It mocks the file reading and the LLM's
#     # prediction to control the input and output for verification.
#     @patch('builtins.open', new_callable=mock_open, read_data="Name,Age,Gender,Number,Fever,Cough\nBob,42,Male,456,1,0\n")
#     @patch('langchain.llms.OpenAI.predict')
#     def test_generate_patient_report_with_relevant_traits(self, mock_predict, mock_file):
#         mock_predict.side_effect = [
#             "Patient Bob, age 42, male, assigned number 456 presents with initial observations warranting a comprehensive health assessment. Further investigation into the presence and duration of these indicators is crucial to formulating an accurate diagnosis and subsequent management plan. A thorough medical history, including any pre-existing conditions and relevant lifestyle factors, should be obtained to contextualize these findings. Additionally, a detailed physical examination, complemented by appropriate diagnostic testing, will be essential in elucidating the underlying etiology. Continuous monitoring of the patient's clinical course will guide the ongoing evaluation and refinement of the treatment strategy.",
#             "Fever, characterized by an elevation in body temperature beyond the normal range, can be indicative of various underlying conditions, most commonly infectious processes. It is essential to ascertain the onset, duration, pattern, and associated symptoms to narrow down the potential causes. Management strategies typically involve addressing the underlying etiology while providing symptomatic relief to ensure patient comfort and prevent complications.",
#         ]
#         report = sun.generate_patient_report("dummy.csv")
#         self.assertIn("Patient Bob, age 42, male, assigned number 456", report)
#         self.assertIn("Fever", report)
#         self.assertNotIn("Cough", report)
#         self.assertEqual(mock_predict.call_count, 2) # One for intro, one for fever

#     # Tests the generate_patient_report function when the patient data in the CSV
#     # does not contain any relevant traits (all trait values are '0'). It mocks
#     # the file reading and the LLM's prediction to ensure a basic introduction is generated.
#     @patch('builtins.open', new_callable=mock_open, read_data="Name,Age,Gender,Number,Headache\nCharlie,28,Other,789,0\n")
#     @patch('langchain.llms.OpenAI.predict')
#     def test_generate_patient_report_without_relevant_traits(self, mock_predict, mock_file):
#         mock_predict.side_effect = [
#             "Patient Charlie, age 28, other, assigned number 789 presents for a routine evaluation. At this time, there are no specific indicators or positive findings reported that necessitate immediate concern or further in-depth analysis. A general overview of the patient's well-being suggests a stable baseline, and no acute symptoms have been identified. Continued monitoring and adherence to recommended health maintenance guidelines will be advised to ensure the ongoing health and wellness of the individual. Should any new symptoms or concerns arise, prompt medical attention should be sought for appropriate assessment and management."
#         ]
#         report = sun.generate_patient_report("dummy.csv")
#         self.assertIn("Patient Charlie, age 28, other, assigned number 789", report)
#         self.assertNotIn("Headache", report)
#         self.assertEqual(mock_predict.call_count, 1) # Only for basic intro

#     # Tests the main function's behavior when run with the default prompt (no
#     # --custom_prompt argument provided). It mocks the command-line argument parsing,
#     # the generate_patient_report function, and the print function to verify the
#     # correct function call and output.
#     @patch('argparse.ArgumentParser.parse_args', return_value=unittest.mock.Mock(csv_file='test_data.csv', custom_prompt=None))
#     @patch('code_files.st_to_unst_module.st_to_unst.generate_patient_report')
#     @patch('builtins.print')
#     def test_main_with_default_prompt(self, mock_print, mock_generate, mock_args):
#         mock_generate.return_value = "Generated report with default prompt."
#         sun.main()
#         mock_generate.assert_called_once_with('test_data.csv', None)
#         mock_print.assert_called_with("\n--- Patient Health Record ---")
#         mock_print.assert_called_with("Generated report with default prompt.")

#     # Tests the main function's behavior when run with a custom prompt provided
#     # via the --custom_prompt argument. It mocks the command-line argument parsing,
#     # the generate_patient_report function, and the print function to verify the
#     # correct function call and output.
#     @patch('argparse.ArgumentParser.parse_args', return_value=unittest.mock.Mock(csv_file='custom.csv', custom_prompt='Custom prompt here'))
#     @patch('code_files.st_to_unst_module.st_to_unst.generate_patient_report')
#     @patch('builtins.print')
#     def test_main_with_custom_prompt(self, mock_print, mock_generate, mock_args):
#         mock_generate.return_value = "Generated report with custom prompt."
#         sun.main()
#         mock_generate.assert_called_once_with('custom.csv', 'Custom prompt here')
#         mock_print.assert_called_with("\n--- Patient Health Record ---")
#         mock_print.assert_called_with("Generated report with custom prompt.")

#     # Tests the scenario where the generate_patient_report function is called
#     # with a CSV file path that does not exist. It mocks the built-in open function
#     # to raise a FileNotFoundError and verifies that the correct error message is returned.
#     @patch('builtins.open', side_effect=FileNotFoundError)
#     def test_generate_patient_report_file_not_found(self, mock_file):
#         report = sun.generate_patient_report("nonexistent.csv")
#         self.assertEqual(report, "Error: CSV file not found at nonexistent.csv")

#     # Tests the robustness of the generate_patient_report function when the CSV
#     # data contains unexpected or invalid values (not '0' or '1') for the traits.
#     # It mocks the file reading and checks if the basic patient information is still
#     # present in the output. The LLM's handling of such input is not strictly tested here.
#     @patch('builtins.open', new_callable=mock_open, read_data="Name,Age,Gender,Number,Trait\nEve,60,Female,999,invalid\n")
#     def test_generate_patient_report_handles_invalid_data(self, mock_file):
#         report = sun.generate_patient_report("invalid_data.csv")
#         self.assertIn("Patient Eve, age 60, female, assigned number 999", report) # Basic info should still be present

if __name__ == '__main__':
    unittest.main()


