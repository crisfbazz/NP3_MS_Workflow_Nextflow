## Workflow Main Commands

The NP³ MS Workflow was developed in a Node.js CLI capable of automatically running the entire workflow, except for the first step. It can also execute each step of the workflow separately, plus the visualization commands and the command to include the GNPS library identification results in the results. 

The Steps 2 to 10 can be automatically executed with the NP³ MS Workflow command **run**. And different results from the NP³ MS Workflow may be joined using the command **join_jobs**, which automatically execute Steps 3 to 10 with adaptations for joining previous results in an incremental clustering approach.

**Final reports** are created at the end of the **run**, **join_jobs** and **clean** commands. These reports are separated in quantification, chemical and network statistics and analysis of the results. Detailed in the run command [Final Reports](../cmd_run/final_reports.md) section. The PCA plotting can be executed separated using the **pca_plot** command.
 
This workflow also contains two interactive commands for data visualization:
 
* **chr**: The first extracts chromatogram from raw MS1 data files and save the images to PNG files. Depending on the chosen parameters this can be a total ion chromatogram (TIC), a base peak chromatogram (BPC) or an extracted ion chromatogram (XIC), extracted from each sample/file with customized groups, and m/z and retention time windows. 
* **spectra_viewer**: The second visualizes and compares MS2 data from MGF files or peak lists. It can also save the images to PNG or SVG files. Currently, it is only supported for Unix OS.
 
Furthermore, the user may manually identify the resulting consensus spectra against the GNPS or GNPS2 [8] online spectral libraries and use another separated command of this workflow called **gnps_result** to join the GNPS identification results to the NP³ MS Workflow quantification tables. The GNPS2 Library Search Workflow may also be executed by this workflow offline using the command **gnps_library_search**, fully integrated with the NP³ results.

## Commands Overview

The NP³ MS Workflow executable is found in the root folder of the repository (named np3_workflow.js) and can be run with the command:
 
```{ .text .copy }
node np3_workflow.js [cmd] [options]
```
 
Where **cmd** are the available NP³ MS Workflow commands that the script can handle (the workflow steps) and **options** are the list of available options for each command.
 
The list of available commands is described below together with their respectively mandatory options (the full list of options is described in each command section). These lists can also be found by running the following command in the terminal:
 
```{ .text .copy }
node np3_workflow.js --help
```
 
Or simple:
 
```{ .text .copy } 
node np3_workflow -h
```
 
From this point on, when NP³ MS Workflow commands are described, the following conventions will be used: 

- Angled brackets (e.g., \<x\>) indicate a required input. 
- Square brackets (e.g. [y]) indicate an optional input. 
- The brackets should not be typed while running the command, they are only used to indicate the type of the option (see the examples of the commands).
 
Commands:
 
- **setup** : Check if the NP³ MS Workflow dependencies are installed, try to install missing R and python packages and compile the NP3_MSCluster algorithm. 
Also executes the setup of the libraries used for spectral searching (UNPD and GNPS2). 
 
- **run** [options] : Steps 2 to 10: Runs the entire NP³ MS Workflow.
    - List of mandatory options:
    - *\-n, \-\-output_name* \<name\>      : the job name. It will be used to name the output directory and the results of the final clustering integration step
    - *\-m, \-\-metadata* \<file\>          : path to the metadata table CSV file
    - *\-d, \-\-raw_data_path* \<path\>     : path to the folder containing the input LC-MS/MS raw spectra data files (mzXML format is recommended)
    - *\-o, \-\-output_path* \<path\>       : path to where the output directory will be created
 
- **pre_process** [options] : Step 2: This command runs the pre-process of the LC-MS/MS raw data. It extracts the list of MS1 peaks with their dimension information (minimum and maximum retention times, the peak area and ID) in each sample, matches the MS2 spectra retention time and precursor m/z against this list and assign to each MS2 spectra a MS1 peak that encompasses it. Additionally, a table with the MS1 peaks without any MS2 spectrum m/z and retention time match are stored in a count table of non-fragmented MS1 peaks.
    - List of mandatory options:
    -  *\-n, \-\-data_name* \<name\>  <x>      :   the data collection name for printting in the verbose messages
    -  *\-m, \-\-metadata* \<file\> <x>         :    path to the metadata table CSV file
    -  *\-d, \-\-raw_data_path* \<path\> <x>   :    path to the folder containing the input LC-MS/MS raw spectra data files
 
- **clustering** [options] : Steps 3 and 4: This command runs the NP3_MSCluster algorithm to perform the clustering of pre-processed MS/MS data into a collection of consensus spectra. Then, it runs the consensus spectra quantification to count the number of spectra and peak area by sample (clustering counts). If necessary, it runs Step 2. And it can also run the library spectra identifications (Step 6) for the collection of consensus spectra.
    - List of mandatory options:
    - *\-n, \-\-output_name* <name>   :       the job name. It will be used to name the output directory and the results of the final clustering integration step
    - *\-m, \-\-metadata*  <file>       :      path to the metadata table CSV file
    - *\-d, \-\-raw_data_path*  <path>    :     path to the folder containing the input LC-MS/MS raw spectra data files
    -  *\-o, \-\-output_path*  <path>     :     path to where the output directory will be created

- **clean** [options] : Step 5: This command runs the pairwise comparisons of the collection of consensus spectra (if not done yet) and then runs the cleaning of the clustering counts. It also runs Step 7 to annotate possible ion variants using the new clean count tables and to create the molecular network of annotations, and runs Step 10 to overwrite any old computation of the molecular network of similarities. It can also run the library spectra identifications (Step 6) for the collection of clean consensus spectra.
    - List of mandatory options:
    - *\-m, \-\-metadata* \<file\>          : path to the metadata table CSV file
    - *\-o, \-\-output_path* \<path\>       : path to the final output data folder, inside the 'outs' directory of the clustering result folder. It should contain the 'mgf' folder and the 'count_tables' folder with the peak area and spectra count tables in CSV files. The job name will be extracted from here
    - *\-y, \-\-processed_data_dir* \<path\>  : the path to the folder inside the raw data folder where the pre-processed data (MGFs) were stored.
 
- **tremolo** [options] : Step 6: (for Unix OS and positive ion mode only) This command runs the tremolo tool, a spectral identification tool used for spectra matching against In-Silico predicted MS/MS spectrum of Natural Products Database (ISDB) from the UNPD (Universal Natural Products Database). It also includes origin and class information of the compounds from NPClassifier, NPAtlas and ClassyFire.
    - List of mandatory options:
    - *\-o, \-\-output_path* <path>   :     path to where the spectral library search results will be stored
    - *\-g, \-\-mgf* <file>       :      path to the input MGF file with the MS/MS spectra data to be searched and identified
    
- **annotate_protonated** [options] : Step 7: (for positive ion mode only) This command runs the annotation of possible ionization variants in the clean count tables and creates the molecular network of annotations. It searches for adducts, neutral losses, multiple charges, dimers/trimers, isotopes and in-source fragmentation based on numerical equivalences and chemical rules. Finally, it runs a link analysis in the molecular network of annotations to assign some of the consensus spectra as putative [M+H]+ representatives.
    - List of mandatory options:
    - *\-m, \-\-metadata* \<file\>          : path to the metadata table CSV file
    - *\-o, \-\-output_path* \<path\>       : path to the output data folder, inside the 'outs' directory of the clustering result folder. It should contain the 'counts_table' folder and inside it the 'clean' subfolder with the clean count tables in CSV files. The name of the output folder will be used as the job name.
 
- **merge** [options] : Step 8: (for positive ion mode only) This command runs the merge of the clean count tables based on the annotated variants. It creates new symbolic spectra candidates representing the union of each spectra with its annotated variants. By default the merge is only performed for the consensus spectra assigned as a [M+H]+ representative ion, to better account for the quantifications of the true metabolites.
    - List of mandatory options:
    - *\-o, \-\-output_path* \<path\>       : path to the output data folder, inside the 'outs' directory of the clustering result folder. It should contain the 'counts_table' folder and inside it the 'clean' subfolder with the clean count tables in CSV files. The job name will be extracted from here
    - *\-y, \-\-processed_data_dir* \<path\>  : the path to the folder inside the raw data folder where the pre-processed data (MGFs) were stored.
    - *\-m, \-\-metadata* \<file\>          : path to the metadata table CSV file
 
- **corr** [options] : Step 9: This command runs the bioactivity correlation to rank the consensus spectra based on the scores computed for the selection of samples and bioactivity values present in the metadata table.
    - List of mandatory options:
    - *\-b, \-\-metadata* \<file\>\     : path to the metadata table CSV file. Used to retrieve the biocorrelation groups
    - *\-c, \-\-count_file_path* \<file\>\ \ \ \ \ \ \ \ \ \ : path to the count table CSV file
 
- **mn** [options] : Step 10: This command runs the creation of a molecular network of similarity based on the pairwise spectra similarity value above the given similarity cut-off. Then, a filter is applied on this network to limit the number of neighbors of each node (number of links) to the top K most similar ones and to limit the size of the components to a maximum number of nodes. The final filtered network contains components that represent the most analogous spectra, possible connecting spectra from similar chemical classes.
    - List of mandatory options:
    - *\-o, \-\-output_path* \<path\>       : path to the output data folder, inside the outs directory of the clustering result folder. It should contain the 'molecular_networking' folder and inside it the 'similarity_tables' folder. The job name will be extracted from here
 
- **gnps_library_search** [options] : Step 6.1: This command runs the GNPS2 Library Search workflow (offline) to identify the informed spectra against the ALL_GNPS_NO_PROPAGATED library (default for LC data). 
    - List of mandatory options:
    - *\-g, \-\-input_mgf_file* \<path\>      :  path to the input MGF file with the MS/MS spectra data to be searched and identified
    - *\-o, \-\-output_path* \<path\>         :  if the input is a NP3 result, the path to the final NP3 output data folder, inside the outs directory of the clustering result folder. It should contain the "identifications" folder, if not it will be created and the results will be stored in it. If the input is not a NP3 result, this should be a chosen result folder. The job name (output_name) will be extracted from here (basename).

- **gnps_result** [options] : This command join the GNPS library identification result from the Molecular Networking (download 'clustered spectra') or the Library Search (download 'All identifications') workflows to the count tables of the NP³ clustering or clean steps. 
    - List of mandatory options:
    - *\-i, \-\-cluster_info_path* \<path\>       :  If joining the result of a Molecular Networking job, this should be the path to the file inside the folder named \'clusterinfo\' of the downloaded output from GNPS. Not used for results coming from the Library Search workflow. (default: "")
    - *\-s, \-\-result_specnets_DB_path* \<path\>       : If joining the result of a Molecular Networking job, this should be the path to the file inside the folder named \'result_specnets_DB\'; and if this is the result of a Library Search workflow, this should be the path to the file inside the downloaded folder
    - *\-c, \-\-count_file_path* \<path\>       : Path to any of the count tables (peak_area or spectra) resulting from the NP³ MS Workflow clustering or clean steps. If the peak_area is informed and the spectra table file exists in the same path (or the opposite), it will merge the GNPS results to both files
    - *\-o, \-\-output_path* \<path\>    :          path to the final output data folder, inside the outs
                                        directory of the clustering result folder. It should contain
                                        the identifications folder, if not, it will be created. The job
                                        name (output_name) will be extracted from here.
   
- **pca_plot** [options] : This command creates a PCA plot of a new data in the NP³ reference chemical space, composed by UNPD+DrugBank+Allosteric. The procedure to create the PCA plots is to first compute the CDK descriptors of a provided list of SMILES and then to create the NP³ reference PCA and to transform these new data to the created chemical space.
    - List of mandatory options:
    - *\-t, \-\-table_path* \<path\>       :  The path to a table in CSV format containing SMILES string in a
                              column. Optionally, it may also contain the types/categories of
                              each SMILES entry in another column. It must be comma ","
                              separated.
    - *\-s, \-\-smiles_column* \<path\>       :  The name of the column in table_path containing the SMILES string.
                              
- **chr** [options] :        This command runs an interactive prompt to extract chromatogram(s) from raw MS1 data files (mzXML, mzData and mzML) and to save to PNG image files. Depending on the provided parameters this can be a total ion chromatogram (TIC - default), a base peak chromatogram (BPC) or an extracted ion chromatogram (XIC) extracted from each sample/file.
 
<!-- - **compare_spectra** [options] : An interactive prompt to compare two spectra from a MGF file, to plot them against it other and to save the image to a PNG file. -->
<!--     - List of mandatory options: -->
<!--     - *\-g, \-\-mgf* \<path\>        path to the input MGF file with the MS/MS spectra data to be compared -->
    
- **spectra_viewer** [options] : (for Unix OS only) This command runs an interactive Web App to visualize and compare MS2 spectra. It receives as input a MGF file or a peak list. It is also possible to manipulate, filter, calculate similarity of the spectra and save PNG or SVG plots. 
    
- **test** [options]      :       This command runs some use cases to test the NP³ MS Workflow consistency in all steps. This option is intended for debugging purposes, and is not a part of the analysis workflow.

## References
[^1]: C.F. Bazzano, et al. (2024). *NP³ MS Workflow*. Analytical Chemistry 96 DOI: 10.1021/acs.analchem.3c05829