## Command **gnps_result**
 
Join the GNPS library identification result from the Molecular Networking (download 'clustered spectra') or the 
Library Search (download 'all identifications') workflows to the count tables of the NP³ clustering or clean steps 
(Steps 3 or 5) from a run or a join_jobs result.

It computes the CDK top descriptors for the unique SMILES identified in GNPS and then make the curation of the 
identifications to score and filter the more reliable results. At the end, it also performs the final identification 
curation from UNPDxGNPS.

## Parameters

For the complete list of parameters run:
 
```{ .text .copy }
node np3_workflow.js gnps_result --help
```

Parameters for the command **gnps_result**:
 
- *\-i, \-\-cluster_info_path* <path\>       :  If joining the result of a Molecular Networking job, this should be the path to the file inside the folder named \'clusterinfo\' of the downloaded output from GNPS. Not used for results coming from the Library Search workflow. (default: "")
- *\-s, \-\-result_specnets_DB_path* <path\>       : If joining the result of a Molecular Networking job, this should be the path to the file inside the folder named \'result_specnets_DB\'; and if this is the result of a Library Search workflow, this should be the path to the file inside the downloaded folder. When using GNPS2, this should be the top 1 results.
- *\-c, \-\-count_file_path* <path\>       : Path to any of the count tables (peak_area or spectra) resulting from the NP³ MS workflow clustering or clean steps. If the peak_area is informed and the spectra table file exists in the same path (or the opposite), it will merge the GNPS results to both files
- *\-o, \-\-job_output_path* <path\>    :          path to the job final output data folder, inside the outs
                                        directory of the clustering result folder. It should contain
                                        the identifications folder, if not it will be created. The job
                                        name (output_name) will be extracted from here.
- *\-m, \-\-metadata* [file]   : path to the samples metadata table CSV file of the NP³ job. This is necessary to plot the distribution of the superclasses grouping by sample, it may be missing if this plot is not desired (leave as empty string).  (default: "")
- *\-h, \-\-help* : output usage information
 
## Results
 
The following columns with the GNPS results are added to the count tables: "gnps_SpectrumID", "gnps_Adduct", 
"gnps_IonMode", "gnps_Smiles", "gnps_InChIKey", "gnps_CAS_Number", "gnps_Compound_Name", "gnps_LibMZ", 
"gnps_MZErrorPPM", "gnps_MQScore", 
\newline "gnps_LibraryQualityString", "gnps_SharedPeaks", "gnps_Organism", "gnps_npclassifier_superclass", 
"gnps_npclassifier_class", "gnps_npclassifier_subclass" and "gnps_npclassifier_pathway".

See the [GNPS documentation](https://ccms-ucsd.github.io/GNPSDocumentation/spectrumcuration/#adding-single-spectra) 
for the description of these columns.
 
If there is more than one GNPS result for a single msclusterID the results are concatenated with a ';', except for 
the "gnps_Smiles" column which is concatenated with a ',' (ease visualization in cytoscape).

The Identification Curation results are concatenated to these created columns, detailed below.
Additional outputs are created in the job_output_path to the 'identifications' subfolder and to the final reports, 
also detailed below.

 
## Examples
 
Fake example to use the GNPS Identification Join command:
 
```{ .text .copy } 
node np3_workflow.js gnps_result --cluster_info_path 
"/path/to/the/output/dir/GNPS_result/clusterinfo/file" 
--result_specnets_DB_path 
"/path/to/the/output/dir/GNPS_result/result_specnets_DB/file.tsv" 
--ms_count_path "job_output_path/outs/output_name/count_files/clean/count_table.csv" 
--job_output_path "/path/to/the/output/NP3/output_name/"
```
 
## Details
 
The GNPS library identifications results should be obtained using the collection of consensus spectra present in the 
MGF files from the NP³ MS workflow output, from the clustering (Step 3) or the clean (Step 5) results. 
If joining the identifications to the clean table, the clean mgf must be used in the GNPS workflows. If joining the 
identifications to the clustering table, the all mgf must be used in the GNPS workflows. For GNPS2 workflows, the 
top 1 result is used for joining (not the top k) - only one identification (the best) is expected for each spectrum.
 
In the GNPS Molecular Networking (MN) workflow result, the 'clusterinfo' file contains the column 'SpecIdx' which is 
equal to the SCANS numbers present in the NP³ MS workflow MGF file used for the identification job. 
The SCANS numbers of the MGF files are equal to the msclursterIDs of the NP³ MS workflow count tables for the 
respective result (clean MGF for clean and all MGF for clustering). The column 'ClusterIdx' of the GNPS MN output 
'clusterinfo' file is equivalent to the column '#Scan' of the 'result_specnets_DB' file, and they are used to join the 
information of both files to obtain for each identification result the correct reference (the 'SpecIdx' column) to the 
msclusterIDs present in the NP³ MS workflow count tables.

In the GNPS or GNPS2 Library Search (LS) workflow result, the '#Scan' column present in the output file is equal to the 
SCANS numbers of the MGF files. These SCANS are equal to the msclusterIDs of the NP³ MS workflow count tables. 
In this case, the join is directly performed.

Finally, the joined GNPS results are also stored to a separated file called 'gnps_results_smiles.csv' inside the 
job_output_path/identifications folder. It will be used to compute the CDK descriptors for the unique SMILES from the 
GNPS result.

At the end, the GNPS identification curation is also performed to classify the retrieved library annotations 
(see respective [section](gnps_library_search.md/#gnps-identification-curation-and-unpdxgnps-final-curation-statistics) 
for more details). 
It is followed by the final identification curation (best from UNPD and GNPS) which is described in another 
[section](final_identification_curation.md).
