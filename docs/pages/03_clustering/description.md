# Step 3 and 4: Command **clustering**
 
Runs the NP3_MSCluster algorithm to perform the **clustering** of pre-processed MS/MS data into a collection of consensus 
spectra. It relies on spectra similarity and chromatographic dimensions. 

And then, runs the consensus spectra **quantification**
(Step 4) to parse the clustering results and count the number of spectra and of peak area by sample. Finally, it computes 
appropriate indicators based on the sample's types. 

This command can also run the library spectra identifications (Step 6), and if 
necessary, it runs Step 2.

## Parameters

For the complete list of parameters run:
 
```{ .text .copy }
node np3_workflow.js clustering --help
```

Parameters for the command **clustering**:
 
- *\-n, \-\-output_name* \<name\>:       the job name. It will be used to name the output directory and the results from the final clustering integration step
- *\-m, \-\-metadata* \<file\>       :      path to the metadata table CSV file
- *\-d, \-\-raw_data_path* \<path\>   :     path to the folder containing the input LC-MS/MS raw spectra data files (mzXML format is recommended)
- *\-o, \-\-output_path* \<path\>     :     path to where the output directory will be created
- *\-f, \-\-fragment_tolerance* [x]  :    the tolerance in Daltons for fragment peaks. Peaks in the original spectra that are closer than this get merged by the NP3_MSCluster algorithm (default: 0.05)
- *\-z, \-\-mz_tolerance* [x]         :   this is the tolerance in Daltons for the m/z of the precursor that determines if two spectra will be compared and possibly joined. Used in the clustering job and in the library identifications (Step 6) (default: 0.025)
- *\-p, \-\-ppm_tolerance* [x]        :   the maximal tolerated m/z deviation in parts per million (ppm) to be used in the pre-processing step if needed (default: 5)
- *\-a, \-\-ion_mode* [x]            :    the precursor ion mode. One of the following numeric values corresponding to a ion adduct type: '1' = Positive [M+H]+ or '2' = Negative [M-H]- (default: 1)
- *\-s, \-\-similarity* [x]           :   the minimum similarity to be consider in the hierarchical clustering, starts in 0.70 and decrease to X in 15 rounds (default: 0.55)
- *\-g, \-\-similarity_blank* [x]      :  the minimum similarity to be consider in the hierarchical clustering of the blank clustering steps, starts in 0.70 and decrease to X in 15 rounds. Only used in the clustering of blank samples (column SAMPLE_TYPE equals 'blank' in the metadata table) (default: 0.3)
- *\-t, \-\-rt_tolerance* [x,y] :        tolerances in seconds for the retention time width of the precursor that determines if two spectra will be compared and possibly joined. The first tolerance [x] is used in the data, blank and batches integration clustering steps; and the second tolerance [y] is for the final integration step between all samples (default: 1,2)
- *\-y, \-\-processed_data_name* [x]  :   the name of the folder inside the raw_data_path where the pre-processed data was stored. If the given folder does not exist it will be created and the pre-process will be run using the default values of the missing options. Otherwise, it will depend on the processed_data_overwrite parameter value (default: "processed_data")
- *\-q, \-\-processed_data_overwrite* [x] :  A logical "TRUE" or "FALSE" indicating if the pre processed data present in the processed_data_name folder should be used (FALSE) and partially incremented if needed, in case it already exists, or overwritten and pre processed again (TRUE) with the default values of the missing Step 2 options (default: "FALSE")
- *\-c, \-\-scale_factor* [x]          :  the scaling method to be used in the fragmented peak's intensities before any dot product comparison (Step 3). Valid values are: 0 for the natural logarithm (ln) of the intensities; 1 for no scaling; and other values greater than zero for raising the fragment peaks intensities to the power of x (e.g., x = 0.5 is the square root scaling). [x] >= 0 (default: 0.5)
- *\-e, \-\-method* [name]             :  a character string indicating which correlation coefficient is to be computed. One of "pearson", "kendall", or "spearman" (default: "spearman")
- *\-x, \-\-min_peaks_output* [x]       : the minimum number of fragment peaks that a spectrum must have to be outputted after the final clustering step (Step 3). Spectra with less than X peaks will be discarded. x >= 1 (default: 5)
<!-- - *\-i, \-\-metfrag_identification* [x]  : a logical "TRUE" or "FALSE" indicating if MetFrag tool should be used for identification search against the PubChem database of the top correlated spectra. (default: "FALSE") -->
- *\-j, \-\-tremolo_identification* [x]  : (not Windows OS's) A logical "TRUE" or "FALSE" indicating if the Tremolo tool should be used for the spectral matching against the ISDB from the UNPD (default: "TRUE")
- *\-b, \-\-max_chunk_spectra* [x]      : Maximum number of spectra to be loaded and processed in a chunk at the same time. In case of memory issues this value should be decreased (default: 3000)
- *\-v, \-\-verbose* [x]              :   for values X>0 show the scripts output information (default: 0)
- *\-h, \-\-help*                     :   output usage information
 
 
## Results
 
A directory inside the `output_path` named with the `output_name` containing:
 
- A copy of the *metadata* file and the command line parameters values used in a file named 'logRunParms', for reproducibility
<!-- - a folder named 'spec_lists' containing the list of files used in each step run of the NP3_MSCluster algorithm;  -->
- A folder named 'outs' with the clustering steps results in separate folders containing: 
    - A subfolder named 'count\_tables' with the Step 4 quantification's in CSV tables named as 
  '<step_name>\_(spectra|peak\_area).csv'.
    - A subfolder named 'clust' with the cluster's membership files (which SCANS or msclusterID were joined)
    - A subfolder named 'mgf' with the resulting clusters consensus spectra in MGF files
    - A text file named 'logClusteringOutput' with the NP3_MSCluster log output. 
 
The clustering steps results folders are named as 'B\_\<DATA\_COLLECTION\_BATCH\>\_\<X\>' where *DATA\_COLLECTION\_BATCH\* 
is the data collection batch number in the metadata file of each group of samples and *X* is 0 if it is the result of a 
*data clustering step* or 1 if it is the result of a *blank clustering step*. 
The *data collection batch integration step* results are stored in folders named as 'B_\<DATA\_COLLECTION\_BATCH\>'.
 
The final integration step result is located inside the 'outs' directory in a folder named with the `output_name`. 
This folder contains the final counts and is where the user will find the final results. It also contains the tremolo 
identification results inside the 'identifications' folder.
 
## Examples
 
Fake examples to execute the Clustering command:
 
```{ .text .copy } 
node np3_workflow.js clustering --output_name "test_np3" --output_path "/path/where/the/output/will/be/stored" --metadata 
"/path/to/the/metadata/file/test_np3_metadata.csv" --raw_data_path 
"/path/to/the/raw/data/dir"
```
 
Same example with different retention time tolerances:

```{ .text .copy } 
node np3_workflow.js clustering -n "test_np3_rt_tol" -o "/path/where/the/output/will/be/stored" -m 
"/path/to/the/metadata/file/test_np3_metadata.csv" -d "/path/to/the/raw/data/dir" -t 3.5,5
```