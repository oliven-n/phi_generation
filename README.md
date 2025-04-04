# phi-generation
This package (which is currently under development) aims to improve upon existing methods for generating simulated patient health records.

The package will have the following components, designed to be run separately or sequentially as a pipeline: 
1. Retrieval Augmented Generation Module 
> Qualitiatively improve style and structure of unstructured Patient Health Records (e.g., Doctor's notest) through use of Retrieval Augmented Generation (RAG) to supplement model records.
> *Note:* While quicker and with less data demands than fine-tuning a model for PHR production, it is possible that the.  For this reason and others, the error propagation module is included.
2. Structured to Unstructured Record Module
> Prompts LLM to infer a patient profile based on unstructured record.
3. Unstructured to Structured Record Module
> Converts written reports on patient to dataframes for diseases, symptoms, and treatments, with binary values reflecting if the patient presents with these.
4. Telephone Module: PHI Generation Error Propagation
> This module repeatedly converts between Unstructured and Structured PHI several times, to test the accuracy of the structured record when compared against a manually-verified ground truth. (Heh, get it, like the telephone game?)
> To test the accuracy of unstructured <-> structured record conversion and the accuracy of different prompting choices to the LLM, and to determine the future need for fine-tuning based on whether RAG inserts false positives based on documents.  Better design choices will result in records that are more robust against type conversions.
