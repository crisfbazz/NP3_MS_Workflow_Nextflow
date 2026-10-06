# Step 5: Command **clean**
 
Compute the pairwise comparisons of the consensus spectra (if not done yet) and then clean the clustering counts. 

It also runs Step 7 to annotate possible ion variants using the new clean tables and to create the molecular network of 
annotations, and runs Step 10 to overwrite any old computation of the molecular network of similarities. 
It can also run the library spectra identifications (Step 6) for the new collection of clean consensus spectra. 
At the end, the final report is (re)computed.

## Parameters

For the complete list of parameters run:
 
```{ .text .copy }
node np3_workflow.js clean --help
```

Parameters for the command *clean*:
 
- *\-m, \-\-metadata* <file\>          : path to the metadata table CSV file
- *\-o, \-\-output_path* <path\>       : path to the final output data folder, inside the 'outs' directory of the clustering result folder. It should contain the 'mgf' folder and the 'count_tables' folder with the peak area and spectra count tables in CSV files. The job name will be extracted from here
- *\-y, \-\-processed_data_dir* <path\>  : the path to the folder inside the raw data folder where the pre-processed data (MGFs) were stored.
- *\-z, \-\-mz_tolerance* [x]              : the tolerance in Daltons that determines if two spectra will be compared and possibly joined. It is also used in Step 7 to detect possible ion variants (default: 0.025)
- *\-t, \-\-rt_tolerance* [x,y]        : tolerances in seconds for the retention time width of the precursor that determines if two spectra will be compared and possibly joined. It is directly applied to the retention time minimum (subtracted) and maximum (added) of the spectra. It enlarges the peak boundaries to deal with disaligned samples or ionization variant spectra. The first tolerance [x] is used in the annotation Step 7; and the tolerance [y] is used in the clean Step 5 (default: 1,2)
- *\-a, \-\-ion_mode* [x]             :  the precursor ion mode. One of the following numeric values corresponding to an ion adduct type: '1' = Positive [M+H]+ or '2' = Negative [M-H]- (default: 1)
- *\-i, \-\-similarity_function* [x]  :    the similarity function to be used in the spectra comparison to create the pairwise similarity tables after clustering and clean steps. One of "np3_shifted_cosine" or "spec2vec". If "spec2vec" is selected, the model trained on UniqueInchikey subset (12,797 spectra) is used by spec2vec in the spectra comparison and the matchms library is used to compute the number of matched peaks between the compared spectra; otherwise, the $NP^{3}$ shifted cosine function is used. (default: "np3_shifted_cosine")
- *\-s, \-\-similarity* [x]          : the similarity to be consider when joining clusters and merging their counts in not blank samples (default: 0.55)
- *\-g, \-\-similarity_blank* [x] :    the similarity to be consider when joining clusters and merging their counts in blank samples (column SAMPLE_TYPE equals 'blank' in the metadata table) (default: 0.3)
- *\-w, \-\-similarity_mn* [x] :     the minimum similarity score that must occur between a pair of consensus spectra to connect them with an edge in the molecular networking. Lower values will increase the component size of the clusters by inducing the connection of less related spectra; and higher values will limit the components sizes to the opposite (default: 0.6)
- *\-f, \-\-fragment_tolerance* [x]  : the tolerance in Daltons for fragment peaks. Peaks in a cluster spectrum that are closer than this are considered the same. (default: 0.05)
- *\-\-bflag_cutoff* [x]  :  A positive numeric value to scale the interquartile range (IQR) of the blank spectra basePeakInt distribution from the clustering result and to allow spectra with a basePeakInt value below this distribution median plus IQR*bflag_cutoff to be joined with a blank spectrum during the clean Step 5, without relying on the similarity value. Or FALSE to disable it. The IQR is the range between the 1st quartile (25th quantile) and the 3rd quartile (75th quantile) of the distribution. The spectra with a basePeakInt value <= median + IQR\*bflag_cutoff (from the blank spectra basePeakInt distribution) and BFLAG TRUE will be joined to a blank spectrum in the clean Step 5. This cutoff will affect the spectra with BFLAG TRUE that would not get joined to a blank spectra when relying only on the similarity cutoff. This is a turn around to the fact that blank spectra have low quality spectra and thus can not fully rely on the similarity values. (default: 1.5)
- *\-\-noise_cutoff* [x] :  A positive numeric value defining the minimum base peak intensity absolute value 
  					that a MS2 spectra must have to be kept after the clustering of Step 3. 
  					The MS2 spectra with a basePeakInt smaller than this value will be removed before clean Step 5.
  					 The default value is zero (disabled). Large values in this parameter may result in the loss of minority 
  					compounds together with noise spectra. (default: 0)
- *\-u, \-\-rules* [x]              :   path to the CSV file with the accepted ionization modification rules for detecting adducts, multiple charge and dimers/trimers variants, and their combination with neutral losses (Step 7). (default: "rules/np3_modifications.csv")
- *\-c, \-\-scale_factor* [x]          :  the scaling method to be used in the fragmented peak's intensities before any dot product comparison (Step 5). Valid values are: 0 for the natural logarithm (ln) of the intensities; 1 for no scaling; and other values greater than zero for raising the fragment peaks intensities to the power of x (e.g. x = 0.5 is the square root scaling). [x] >= 0 (default: 0.5)
- *\-\-min_matched_peaks* [x]           :  The minimum number of common peaks that two spectra must share to be connected by an edge in the filtered SSMN. Connections between spectra with less common peaks than this cutoff will be removed when filtering the SSMN. Except for when one of the spectra have a number of fragment peaks smaller than the given min_matched_peaks value, in this case the spectra must share at least 2 peaks. The fragment peaks count is performed after the spectra are normalized and cleaned. (default: 6)
- *\-k, \-\-net_top_k* [x]            :     the maximum number of connections for one single node in the similarity molecular networking. An edge between two nodes is kept only if both nodes are within each other's X most similar nodes. This restriction is applied to the spectra following the msclusterID's order, so smaller m/z's are limited first. Keeping this value low makes very large networks (many nodes) much easier to visualize (default: 10)
- *\-x, \-\-max_component_size* [x]    :    the maximum number of nodes that all component of 
                      the similarity molecular network must have. The edges of 
                      this network will be removed using an increasing cosine 
                      threshold until each network component has at most X nodes. 
                      Keeping this value low makes very large networks (many nodes 
                      and edges) much easier to visualize. (default: 200)
- *\-\-blank_expansion* [x]      :   the distance of neighborhood nodes from the blanks in IVAMN to be 
  					selected for removal in the final protonated networks. (0) to only remove blanks nodes,
  					(1) to remove nodes directly connected to a blank node, 
  					(2 or greater) to remove nodes in a distance equal to 2 or greater from a blank node, 
  					or (-1) to remove all possible neighbours and ancestors of a blank node (remove blank clusters) from IVAMN  (default: 0)
- *\-r, \-\-trim_mz* [x]                : A logical "TRUE" or "FALSE" indicating if the spectra fragmented peaks around the precursor m/z +-20 Da should be deleted before the pairwise comparisons. If "TRUE" this removes the residual precursor ion, which is frequently observed in MS/MS spectra acquired on qTOFs. (default: "TRUE")
- *\-\-max_shift* [x]                   :  Maximum difference between precursor m/zs that will be used in the search of shifted m/z fragment ions in the $NP^{3}$ shifted cosine function. Shifts greater than this value will be ignored and not used in the cosine computation. It can be useful to deal with local modifications of the same compound. (default: 200)
- *\-l, \-\-parallel_cores* [x]       :   the number of cores to be used for parallel processing (Step 5). x = 1 for disabling parallelization and x > 2 for enabling it. [x] >= 1 (default: 2)
- *\-e, \-\-method* [name]           : a character string indicating which correlation coefficient is to be computed. One of “pearson”, “kendall”, or “spearman” (default: spearman)
<!-- - *\-i, \-\-metfrag_identification* [x]  : a logical "TRUE" or "FALSE" indicating if the MetFrag tool should be used for spectra identification search against the PubChem database of the top correlated spectra. The Metfrag API is not very fast and could raise some errors, if this happens the user must stop the process. (default: "FALSE") -->
- *\-j, \-\-tremolo_identification* [x]  : (not Windows OS's) A logical "TRUE" or "FALSE" indicating if the Tremolo tool should be used for the spectral matching against the ISDB from the UNPD (default: "TRUE")
- *\-b, \-\-max_chunk_spectra* [x]      :   Maximum number of spectra (rows) to be loaded and processed in a chunk at the same time. In case of memory issues this value should be decreased (default: 3000)
- *\-v, \-\-verbose* [x]             : for values X\>0 show the scripts output information. (default: 0)
- *\-h, \-\-help*                    : output usage information
 
## Results
 
One subfolder inside the 'count\_tables' folder is created named 'clean' containing:
 
- A text file named 'analyseCountClusteringClean' with the clean count analyses;    
- Two CSV files with the clustering counts of spectra and peak area cleaned and annotated, named with the suffix '`output_name`\_clean_ann.csv';
- Copies of these CSV files with the correlation columns results are also created when there is a biocorrelation result, named with an additional suffix equals '\_corr_`method`.csv' and '\_corr_`method`\_bioAct.csv'.
 
The 'molecular_networking' folder is also created if not present yet, and inside it is created: 
 
- One subfolder named "similarity\_tables" with the clean version of the pairwise similarity table (n x n, where n is 
the final number of clean consensus spectra);
- Five molecular networks edge files (Steps 7 and 10): 
    - The ionization variant annotation molecular network named as: '`output_name`_ivamn.selfloops';
    - The complete spectra similarity molecular network named as : 
    '`output_name`\_ssmn_w\_`similarity_mn`.selfloops', which contains all links with a similarity value above the cut-off;
    - The filtered spectra similarity molecular network named as 
      '`output_name`\_ssmn_w_`similarity_mn`\_k_`net_top_k`\_x_ 
       `max_component_size`.selfloops';
    - The IVAMN [M+H]+ named as: 
      '`output_name`\_ivamn_protonated.selfloops';
    - The SSMN [M+H]+ filtered named as: 
      '`output_name`\_ssmn_protonated\_w\_`similarity_mn`\_k\_`net_top_k`\_x\_`max_component_size`.selfloops'.
- One CSV table containing the molecular network of annotations attributes and the assigned [M+H]+ representatives, named as 
    '`output_name`_ivamn_attributes.csv' (Step 7).
 
The 'final_reports' folder is also created if not present yet to store the recomputed statistics and analysis of the result.

Where the `output_name` is extracted from the `output_path`;
 
When running Tremolo the 'identification' folder is also created (if not present yet) with the identifications results 
inside it. These identifications are also added as new columns in the created clean count tables. 
 
## Examples
 
Fake example to execute the Cleaning command:

```{ .text .copy }  
node np3_workflow.js clean --metadata 
"/path/to/the/metadata/file/test_np3_metadata.csv" 
--output_path "/path/to/the/output/dir/test_np3/outs/test_np3" -b 1000
``` 
