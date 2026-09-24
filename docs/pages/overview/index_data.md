## Input
 
The main input for the NP³ MS Workflow (run command and major steps) are LC-MS/MS raw data files in mzXML or mzML formats (mzData format could also be used, but was not vastly tested). 
It fully supports the positive ([M+H]+) ion mode; and the negative ([M-H]-) ion mode may also be used, but the Steps 7 and 8 were not vastly tested for it and may lead to undesired results. 
The rest of the pipeline results can be used for the negative ion mode.

The NP³ pre-processed data and the clean data from the NP³ results are the inputs for the join_jobs command.

## Output 

The workflow output is organized in folders and depending on the command can contain: 
 
- The pre-processed spectra in MGF files, containing all detected and filtered MS2 ion enriched with MS1 information. (Step 2)
- The collection of consensus spectra from Step 3 and Step 5 in MGF files
- The consensus spectra quantification results in CSV files, referred to as count tables. 
The count tables contains a consensus spectra ID, named msclusterID (representing a m/z and a retention time interval), by row and its quantification in each sample and indicators by column. The following count tables are created:
    - Clustering count tables (Step 4)
    - Clean count tables, with the clustering count table redundancies aggregated - remove fragmented clusters (Step 5)
    - Merged count tables (Step 8)
- The count tables can also contain the following columns:
    - The libraries identification results (Step 6)
    - The ionization variants annotations grouped by type (Step 7)
    - The list of [M+H]+ representatives, suggesting the number of metabolites in the samples (Step 7)
    - The bioactivity correlation scores, that can be used to rank the consensus spectra. And the quantification grouping (Step 9)
- The molecular networks (MN) of ionization variants annotations (IVAMN) and of spectral similarity (SSMN) files from Steps 7 and 10:
    - IVAMN (Step 7)
    - SSMN (Step 10) 
    - SSMN filtered (Step 10)
    - IVAMN [M+H]+ (Step 10)
    - SSMN [M+H]+ filtered (Step 10)
- The final reports with quantification, chemical and network statistics and analysis of the final results
    - A chemical space of the identified spectra using PCA
 


## References
[^1]: C.F. Bazzano, et al. (2024). *NP³ MS Workflow*. Analytical Chemistry 96 DOI: 10.1021/acs.analchem.3c05829