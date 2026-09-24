## Steps 2 to 10: Command **run**
 
This command runs the entire NP³ MS Workflow pipeline - Steps 2 to 10.

## Parameters

For the complete list of parameters run:
 
```{ .text .copy }
node np3_workflow.js run --help
```
 
Options for the command *run*:
 
- *\-n, \-\-output_name* <name\>      : the job name. It will be used to name the output directory and the results of the final clustering integration step. It must have less than 80 characters.
- *\-m, \-\-metadata* <file\>          : path to the metadata table CSV file
- *\-d, \-\-raw_data_path* <path\>     : path to the folder containing the input LC-MS/MS raw spectra data files 
- *\-o, \-\-output_path* <path\>       : path of the output directory
- *\-f, \-\-fragment_tolerance* [x]  : the tolerance in Daltons for fragment peaks. Peaks in the original MS/MS spectra that are closer than this tolerance are merged in the clustering steps (Step 3). It is also used in the pre-processing (Step 2), in the spectra similarity comparisons and in the cleanning Step 5 (default: 0.05)
- *\-z, \-\-mz_tolerance* [x]              : this is the tolerance in Daltons for the m/z of the precursor that determines if two spectra will be compared and possibly joined. It is used in the clustering steps (Step 3), in the cleaning (Step 5), in the library identifications (Step 6) and in the annotation of ionization variants (Step 7) (default: 0.025)
- *\-p, \-\-ppm_tolerance* [x]        : the maximum tolerated m/z deviation in parts per million (ppm) to be used in the pre-processing (Step 2). Typically set to a generous multiple of the mass accuracy of the mass spectrometer. (default: 15)
- *\-a, \-\-ion_mode* [x]             :  the precursor ion mode. One of the following numeric values corresponding to an ion adduct type: '1' = Positive [M+H]+ or '2' = Negative [M-H]- (default: 1)
- *\-i, \-\-similarity_function* [x]  :    the similarity function to be used in the spectra comparison to create the pairwise similarity tables after clustering and clean steps. One of "np3_shifted_cosine" or "spec2vec". If "spec2vec" is selected, the model trained on UniqueInchikey subset (12,797 spectra) is used by spec2vec in the spectra comparison and the matchms library is used to compute the number of matched peaks between the compared spectra; otherwise, the NP³ shifted cosine function is used. (default: "np3_shifted_cosine")
- *\-s, \-\-similarity* [x]          : the minimum similarity to be considered in the hierarchical clustering of Step 3, starts in 0.70 and decrease to X in 15 rounds. Also used in the Step 5 quantification cleaning (default: 0.55)
- *\-g, \-\-similarity_blank* [x] :    the minimum similarity to be consider in the hierarchical clustering of the blank clustering steps, starts in 0.70 and decrease to X in 15 rounds. Only used in the clustering of blank samples (column SAMPLE_TYPE equals 'blank' in the metadata table) (default: 0.3)
- *\-w, \-\-similarity_mn* [x] :     the minimum similarity score that must occur between a pair of consensus spectra to connect them with a link in the molecular network of similarity. Lower values will increase the components sizes by inducing the connection of less related spectra; and higher values will limit the components sizes to the opposite (default: 0.6) 
- *\-t, \-\-rt_tolerance* [x,y]        : tolerances in seconds for the retention time width of the precursor that determines if two spectra will be compared and possibly joined. It is directly applied to the retention time minimum (subtracted) and maximum (added) of the spectra. It enlarges the peak boundaries to deal with misaligned samples or ionization variant spectra. The first tolerance [x] is used in the data, blank and batches integration steps from the clustering Step 3 and in Steps 2 (if no previous result is provided) and 7 (chemical annotations); and the tolerance [y] is used in the final integration step from the clustering Step 3 and in the clean Step 5 (default: 1,2)
- *\-\-min_matched_peaks* [x]           :  The minimum number of common peaks that two spectra must share to be connected by an edge in the filtered SSMN. Connections between spectra with less common peaks than this cutoff will be removed when filtering the SSMN. Except for when one of the spectra have a number of fragment peaks smaller than the given min_matched_peaks value, in this case the spectra must share at least 2 peaks. The fragment peaks count is performed after the spectra are normalized and cleaned. (default: 6)
- *\-k, \-\-net_top_k* [x]            :     the maximum number of connections for one single node in the molecular network of similarity. A link between two nodes is kept only if both nodes are within each other's X most similar nodes. Keeping this value low makes very large networks (many nodes) much easier to visualize (default: 15)
- *\-x, \-\-max_component_size* [x]   :     the maximum number of nodes that each component of the molecular network of similarity must have. The links of this network will be removed using an increasing cosine threshold until each component has at most X nodes. Keeping this value low makes very large networks (many nodes and links) much easier to visualize. (default: 200)
- *\-\-blank_expansion* [x]      :   the distance of neighborhood nodes from the blanks in IVAMN to be 
  					selected for removal in the final protonated networks. (0) to only remove blanks nodes,
  					(1) to remove nodes directly connected to a blank node, 
  					(2 or greater) to remove nodes in a distance equal to 2 or greater from a blank node, 
  					or (-1) to remove all possible neighbours and ancestors of a blank node (remove blank clusters) from IVAMN  (default: 0)
<!-- - \-m, \-\-model_dir [path]        : directory where model files are kept. If running MSCluster on Windows and not from the current directory you should specify the path to “Models Windows” (default: ./Models) -->
<!-- - \-r, \-\-num_rounds [n]          : determines how many rounds are used for the hierarchical clustering (default: 20) -->
<!-- - \-p, \-\-mixture_prob [x]        : the probability wrongfully adding a spectrum to a cluster (default: 0.4) -->
- *\-c, \-\-scale_factor* [x]          :  the scaling method to be used in the fragmented peak's intensities before any dot product comparison (Steps 3 and 6). Valid values are: 0 for the natural logarithm (ln) of the intensities; 1 for no scaling; and other values greater than zero for raising the fragment peaks intensities to the power of x (e.g. x = 0.5 is the square root scaling). [x] >= 0 (default: 0.5)
- *\-l, \-\-parallel_cores* [x]       :   the number of cores to be used for parallel processing in the spectra pairwise comparison (Step 5). x = 1 for disabling parallelization and x > 2 for enabling it. [x] >= 1 (default: 2)
- *\-y, \-\-processed_data_name* [x]  : the name of the folder inside the *raw_data_path* where the pre-processed data will be stored. If the given folder does not exist it will be created and the pre-process will be run using the default values of the missing options. Otherwise, it will depend on the *processed_data_overwrite* parameter value (default: "processed_data")
- *\-q, \-\-processed_data_overwrite* [x]  : A logical "TRUE" or "FALSE" indicating if the pre processed data present in the processed_data_name folder should be used (FALSE) and partially incremented if needed, in case it already exists, or overwritten and pre processed again (TRUE) with the default values of the missing Step 2 options (default: "FALSE")
- *\-\-bflag_cutoff* [x]            :      A positive numeric value to scale the interquartile range (IQR) of the blank spectra basePeakInt distribution and allow spectra with a basePeakInt value below this distribution median plus IQR\*bflag_cutoff to be joined with a blank spectrum without relying on the similarity value. Or FALSE to disable it. The IQR is the range between the 1st quartile (25th quantile) and the 3rd quartile (75th quantile). The spectra with a basePeakInt value <= median + IQR\*bflag_cutoff of the blank spectra basePeakInt distribution and BFLAG TRUE will be joined to a blank spectrum in the clean step. This cutoff will affect the spectra with BFLAG TRUE that would not get joined to a blank spectra when relying only on the similarity cutoff. This is a turn around to the fact that blank spectra have low quality spectra and thus can not fully rely on the similarity values. (default: 1.5)
- *\-\-noise_cutoff* [x]            :      A positive numeric value to scale the interquartile range (IQR) of the blank spectra basePeakInt distribution from the clustering Step 3 result and to remove the spectra with a basePeakInt value below this distribution median plus IQR*noise_cutoff after the clean Step 5. Or FALSE to disable it. The IQR is the range between the 1st quartile (25th quantile) and the 3rd quartile (75th quantile) of the distribution. When no blank sample is present in the metadata, the full distribution is used. This cutoff will affect the spectra with with a low basePeakInt value that probably are noise features.  If the clustering Step 3 results in more than 25000 spectra, the noise cutoff will be applied before the clean Step 5 to prevent a long processing time (default: "FALSE")
- *\-u, \-\-rules* [x]              :   path to the CSV file following the NP³ rules table format with the accepted ionization modification rules for detecting adducts, multiple charge and dimers/trimers variants, and their combination with neutral losses. To be used by the annotation algorithm (Step 7) (default: "rules/np3_modifications.csv")
- *\-r, \-\-trim_mz* [x]                : A logical "TRUE" or "FALSE" indicating if the spectra fragmented peaks around the precursor m/z +-20 Da should be deleted before the pairwise comparisons. If "TRUE" this removes the residual precursor ion, which is frequently observed in MS/MS spectra acquired on qTOFs. (default: "TRUE")
- *\-\-max_shift* [x]                   :  Maximum difference between precursor m/zs that will be used in the search of shifted m/z fragment ions in the NP³ shifted cosine function. Shifts greater than this value will be ignored and not used in the cosine computation. It can be useful to deal with local modifications of the same compound. (default: 200)
- *\-b, \-\-max_chunk_spectra* [x]      :   Maximum number of spectra (rows) to be loaded and processed in a chunk at the same time. In case of memory issues this value should be decreased. To be used in Steps 6, 7 and 10 (default: 3000)
- *\-e, \-\-method* [name]           : a character string indicating which correlation coefficient is to be computed. One of “pearson”, “kendall”, or “spearman” (Step 9) (default: spearman)
<!-- - *\-i, \-\-metfrag_identification* [x]  : a logical "TRUE" or "FALSE" indicating if the MetFrag tool should be used for spectra identification search against the PubChem database of the top correlated spectra. The Metfrag API is not very fast and could raise some errors, if this happens the user must stop the process (Step 5) (default: "FALSE") -->
- *\-j, \-\-tremolo_identification* [x]  : (not Windows OS's) A logical "TRUE" or "FALSE" indicating if the Tremolo tool should be used for the spectral matching against the ISDB from the UNPD (Step 5) (default: "TRUE")
- *--gnps_search_tool* [x]          :    the GNPS2 search tool to be used in the library searching against the ALL\_GNPS\_NO_PROPOGATED (Step 6.1). One of "gnps_indexed", "gnps", "gnps_new" or "" (disabled). The similarity function is hardcoded to be the cosine, the peak transformation function is the square root and top k equals 5. (default: "gnps_indexed")
- *--gnps_min_cosine* [x]         :      the similarity threshold for the GNPS2 library search that determines if two spectra are a match. The minimum cosine for a search match. Values greater or equal than 0.7 will retrieve more accurate results. (default: 0.7)
- *--gnps_min_matched_peaks* [x]     :   The minimum number of common peaks that the searched spectra must share with a library spectrum to be considered as a match. (default: 6)
- *--gnps_window_filter* [x]       :     If "TRUE", for each peak, it will check a window around that peak, if it is not one of the top peaks in terms of intensity in that window, it will be filtered out. This will speed up the library search and reduce the effect of noise. (default: "TRUE")
- *--gnps_analog_search* [x]      :      If "TRUE", also search for analog spectra in GNPS2 library search. This allows as a match the spectra with different precursor mass and similar fragmentation pattern.  (default: "FALSE")
- *--gnps_analog_max_shift* [x]    :     Maximum difference between precursor m/zs that will be used in the GNPS2 library search when analog_search is enabled. Only used when search_tool is "gnps_new". (default: 400)
- *--gnps_parallel_threads* [x]     :    the number of threads to be used for parallel processing in the search when search_tool is gnps_indexed. For parallelization set x >= 1 (default: 4)
- *\-v, \-\-verbose* [x]             : for values X\>0 show the scripts output information. For values greater or equal to 10 a consistency test of the results is also performed. (default: 0)
- *\-h, \-\-help*                    : output usage information
 
## Results
 
The **run** command result in a directory inside the `output_path` named with the `output_name` containing:
 
- A copy of the *metadata* and the *rules* files and the command line parameters values used in a file named 'logRunParms', for reproducibility
<!-- - a folder named 'spec_lists' containing the list of files used in each step run of the NP3_MSCluster algorithm;  -->
- A folder named 'outs' with the clustering steps results in separate folders containing: 
    - A subfolder named 'count_tables' with the Step 4 quantifications in CSV tables named as '<step_name\>_(spectra|peak_area).csv'.
    - Another subfolder named 'clust' with the clusters membership files (which SCANS or msclusterID were joined)
    - A third subfolder named 'mgf' with the resulting clusters consensus spectra in MGF files
    - A text file named 'logClusteringOutput' with the NP3_MSCluster log output. 
 
The clustering steps results folders are named as 'B\_<DATA\_COLLECTION_BATCH\>\_<X\>' where *DATA\_COLLECTION\_BATCH* 
is the data collection batch number in the metadata file of each group of samples and *X* is 0 if it is the result of a 
*data clustering step* or 1 if it is the result of a *blank clustering step*. 
The *data collection batch integration step* results are stored in folders named as 'B\_<DATA\_COLLECTION\_BATCH\>'.
 
The final integration step result is located inside the 'outs' directory in a folder named with the *output_name*. 
This folder contains the final quantification and is where the user will find the final results. 
It also contains the following data:
 
- Inside the 'count_tables' folder two subfolders named "clean" and "merge" containing CSV tables with the 
quantification and annotations from Steps 4 to 9, and the base peak intensity distribution plot of the clustering counts in a PNG image file;
- The 'identifications' folder with the tremolo and gnps library identification results;
- The 'molecular_networking' folder containing: 
    - One subfolder named "similarity_tables" with the pairwise similarity tables (Step 5);
    - Five molecular networks edge files (Steps 7 and 10): 
        - The ionization variant annoation molecular network named as: \newline
          '<*output_name*\>\_ivamn.selfloops'
        - The complete spectra similarity molecular network named as : \newline
          '<*output_name*\>\_ssmn_w\_<*similarity_mn*\>.selfloops', which contains all links with a similarity value above the cut-off
        - The filtered spectra similarity molecular network named as \newline
          '<*output_name*\>\_ssmn\_w\_<*similarity_mn*\>\_k\_<*net\_top\_k*\>\_x\_ \newline
           <*max_component_size*\>.selfloops'
        - The IVAMN [M+H]+ named as: \newline
          '<*output_name*\>\_ivamn\_protonated.selfloops'
        - The SSMN [M+H]+ filtered named as: \newline
          '<*output\_name*\>\_ssmn\_protonated\_w\_<*similarity_mn*\>\_k\_<*net_top_k*\>\_x\_ \newline 
          <*max\_component\_size*\>.selfloops'
    - One CSV table containing the molecular network of annotations attributes and the assigned [M+H]+ representatives, named as: \newline
      '<*output\_name*\>\_ivamn\_attributes.csv' (Step 7)
- The 'final_reports' folder containing statistics of the final result in separated subfolders for quantification, chemical and network analysis.

## Examples
 
Fake examples using the *run* command:

```{ .text .copy }
node np3_workflow.js run --output_name "test_np3" --output_path "/path/where/the/output/will/be/stored/" --metadata 
"/path/to/the/metadata/file/test_np3_metadata.csv" --raw_data_path "/path/to/the/raw/data/dir/" 
--fragment_tolerance 0.01
```
 
```{ .text .copy }
node np3_workflow run -n "test_np3_rt_tol" -o "/path/where/the/output/will/be/stored"  
-m "/path/to/the/metadata/file/test_np3_metadata.csv" -d "/path/to/the/raw/data/dir" 
-t 3.5,5
```