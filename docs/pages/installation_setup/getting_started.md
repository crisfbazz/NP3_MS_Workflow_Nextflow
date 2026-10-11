# Getting Started - Test Example

Let's use the LC-MS/MS raw samples from the L754_bacs dataset present in the 'test' folder to check the installation 
and run the NP³ MS workflow for the first time. 

A metadata table (Step 1) is already provided in the 'test' folder:
'test/L754_bacs/marine_bacteria_library_L754_metadata.csv'. 

Start by activating the NP³ MS workflow environment. Run the following command in the terminal: 

 
```  
mamba activate np3 
``` 

Now, we can run the NP³ MS workflow with the provided metadata table. 

The test dataset that we will use in this example gives good results with the pre-processing default parameters, 
so we will run the entire workflow with a single command. 

In the repository root folder run the following command in the terminal to execute all the workflow steps, using the 
provided metadata and data for the L754_bacs dataset: 


``` { .text .copy }
node np3_workflow.js run --output_name "L754_bacs_test" --output_path "test/L754_bacs" 
--metadata "test/L754_bacs/marine_bacteria_library_L754_metadata.csv" --raw_data_path "test/L754_bacs/mzxml" 
--tremolo_identification TRUE --verbose 1
``` 

If the workflow executed without any error message, the final results should be located in the folder 
'test/L754_bacs/L754_bacs_test'. The pre-processing result will be located at "test/L754_bacs/mzxml/processed_data". 

#### Results organization

The contents of the final result folder are explained below: 

~~~ 
L754_bacs_test 
│  
├── outs                                      <- the results from the clustering steps in separated folders, and inside them the results from the other workflow steps as described below 
│   │    
│   ├── L754_bacs_test                        <- the final clustering (Step 3) result folder - and final results from other steps!
│   │   │ 
|   |   ├── clust                             <- the folder with clusters membership files (which SCANS or msclusterID were joined in the final clustering step) (Step3) 
│   │   │ 
|   |   ├── count_tables                      <- the folder with the quantification tables from Step 4 named as "L754_bacs_test_(spectra|peak_area).csv" and "L754_bacs_test_peak_area_MS1.csv" 
│   │   │   |    
│   │   │   ├── clean                         <- the folder with the quantification tables from Steps 5, 7 and 9 
│   │   │   |  
│   │   │   └── merge                         <- the folder with the quantification tables from Steps 8 and 9 
│   │   │
|   |   ├── final_reports                     <- the folder with the final reports computed at the end of the processing based on the final clean counts and identifications
│   │   │   |    
│   │   │   ├── chemical_report               <- the folder with the chemical statistics and PCA plots (chemical_space_identifications subfolder)
│   │   │   |    
│   │   │   ├── molecular_networking_report   <- the folder with the molecular networks statistics
│   │   │   |  
│   │   │   └── quantification_report         <- the folder with the quantification statistics
│   │   │
|   |   ├── identifications                   <- the folder with the complete list of identifications from UNPD returned by tremolo and from GNPS2 Library Search              
│   │   │ 
|   |   ├── mgf                               <- the folder with the MGF files from the clustering Step 3 (named L754_bacs_test_all.mgf), containing the complete list of consensus spectra, and from the clean Step 5 (named L754_bacs_test_clean.mgf), containing the final list of clean consensus spectra. 
│   │   │    
|   |   └── molecular_networking              <- the folder with the spectral similarity molecular networks (SSMN) and the ionization variant annotation molecular network (IVAMN) in edge files (Steps 7 and 10). One table with the attributes of the IVAMN
│   │       |   
│   │       └── similarity_tables             <- the folder with the pairwise similarity tables (Step 5) 
│   │    
│   ├── B_X_Y                                 <- the clustering sub steps results folders, where X is the data collection batch number in the metadata table of each group of samples and Y is 0 if it is the result of a data clustering step or 1 if it is the result of a blank clustering step 
│   │ 
│   └── B_X                                   <- the data collection batch integration step results folders, where X is the data collection batch number 
│ 
├── np3_modifications.csv                     <- a copy of the ionization variant rules table used in Step 7
│ 
├── marine_bacteria_library_L754_metadata.csv <- a copy of the metadata table used in the job 
│ 
└── logRunParms                               <- the command line parameters values used in the run command, for reproducibility
~~~ 

For more details about what each file contains, its columns and values, and to better understand the results and how 
they are obtained see each step section for more details! 

## Recommendations to execute NP³ MS Workflow

The user can choose to run the entire workflow (Steps 2 to 10) with a single command 
**[run](../cmd_run/description.md)**, using the default parameters in 
the [pre-processing](../02_pre_process/description.md) Step 2 - as the example above -, or to run the pre-processing 
Step 2 separated, in order to optimize the LC-MS processing result, and then run the rest of the workflow with a single 
command *run* (Steps 3 to 10). 

We **highly recommend the users to run the pre-processing Step 2 separated** to optimize its most critical parameters: 

- The common MS1 peak width found in your data (peak_width); and 
- The expected deviation between the MS1 and the MS2 data (rt_tolerance).

The pre-processing results will directly impact the isomers definition throughout the rest of the 
workflow, a bad optimization of its critical parameters can, for example, lead to large MS1 peaks being split or 
adjacent peaks being joined, negatively impacting the final result. 

The users should **follow the [pre-processing 
optimization guide](../02_pre_process/optimization.md)** present in this command section when executing the workflow 
with a new dataset.
  


