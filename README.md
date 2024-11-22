# Osteoarthritis progression does not influence gut microbiota composition


We analyzed the gut microbiota of 1395 OA patients from 4 European cohorts, namely the Finrisk cohort, Estonia Biobank, Lifelines and UK twin study. We compared them with 1395 healthy control from a pool of total 17’641 participants. There were a total of 191 OA prevalent cases from a total of 8000 participants in Finrisk cohort, 590 out of 7019 in Lifelines, 557 out of 2509 in Estonia cohort and 57 out of 114 from the UK twin study. 
OA was defined with the help of ICD10 codes as M15-19, M16 was used for hip OA, M17 for knee OA, M18 for OA of the first carpometacarpal joint and M19 for other and unspecified OA. This was used for the Finrisk cohort, Estonia Biobank and the UK twin study, whereas for the Lifelines a self-assessment from a questionnaire was used. 

Participants with conditions known to affect the gut microbiome, namely inflammatory bowel syndrome, inflammatory bowel disease, coeliac disease, pregnancy and participants that were medicated for depression or diabetes were excluded. After exclusion 162 OA participants remained in Finrisk cohort, 482 in Estonia Cohort, 431 in Lifelines, and 57 in UK twin study. 

For the age test within the Finrisk cohort participants above the age of 55 were excluded for the “>55” group and participants below the age of 65 were excluded for the “<65” group. A total of 150 participants per group were selected to resemble the number of participants in the OA analysis.  

Next to the stringent exclusion criteria, we implemented a strict control matching for each cohort, to ensure that we limit the main confounding factors. Therein, for each case (OA positive), we selected a matching control (without replacement, so every patient can only be present once in the controls) following the criteria: (i) Gender: Exact match;  (ii) Age: Exact match, otherwise extend to +-1. (iii) BMI, Range of +- 1, otherwise extend to +-2 etc. (iv)  Activity Score. 
These variables were considered in a custom script to look for the best possible control. If a perfect match was found (same gender, same age, same BMI, same activity score), we chose this, otherwise we relaxed the requirements in age, BMI and activity score in a stepwise matter up to a maximum of +-5. In the end, we ensured that the parameters were not significantly different between control and cohort. These criteria are also written up in the "steps_to_select_controls.txt" file.


The quality of the reads was assessed using FastQC v0.11.9. 

Further steps (filtering and taxonomic profiling) were performed for each selected sample. An example worklfow is provided in "example_sample_pipeline.sh". 

Indicidual read count tables were in the end combined for both mapseq as well as mOTUs.

The resulting taxonomic profiles were combined with their metadata, and all statistical analyses (diversi9ties, differential abundance, machine learning..) were performed in R and are provided in "R_motus_final_OA_ctrl_final_clean.R".
