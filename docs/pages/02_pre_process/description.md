## Step 2: Command **pre\_process**
 
Runs the pre-process of the LC-MS/MS raw data. 

It extracts the list of MS1 peaks with their dimension information (retention time minimum and maximum, the peak area 
and ID) in each sample, matches the MS2 spectra retention time and precursor m/z against this list and assign to each 
spectra a MS1 peak that encompasses it. Additionally, a table with the MS1 peaks without a MS2 spectrum m/z and retention 
time match are stored in a count table to account for the not fragmented MS1 peaks. Finally, a diagnostic is performed 
to evaluate the pre-process result followed by a suggestion for better *rt_tolerance* and *peak_width* parameters values, 
based on the list of fragmented MS1 peaks.
 
If *max_samples_batch_align* > 0, compute a suggestion for the retention time misalignment between samples of the same 
data collection batch and between all samples. 
This does not modify the samples retention time, it's only used to suggest a retention time tolerance value for the 
following commands (clustering, clean, annotate_protonated and merge).
 
It generates one MGF by sample containing all the detected MS2 spectra (high sensitivity) enriched with their 
respectively matched MS1 peak dimensions, e.g., retention time minimum, maximum, peak area and peak ID (given by the 
peak detection algorithm). When using the same `raw_data_path` to pre-process the LC-MS/MS raw data files using 
different metadata tables, a different *processed_data_name* must be used to avoid overwriting useful files and to 
store the created MGFs in different folders. If the new job intersects with the result present in the same 
`processed_data_name`, the user may choose to maintain this folder and the pre-processing will only process the missing 
samples (partial pre-processing) and reuse and increment the previous result. 

Running the pre-process step with a big data collection (number of samples above ~100) can be time consuming and 
impracticable to repeat this step more than one time to try to improve the parameters values based on the final 
diagnostic and suggestions. In this case, we recommend the user to build a metadata table with only a small selection 
of the samples, choosing the ones that could be more representative of the data collection and all blank samples (to 
correctly remove the blank m/zs). Then, iteratively pre-process these small selection of samples until better results 
are obtained based on the final diagnostic and suggestions for the parameters values. Finally, the user can run the 
pre-process for the complete data collection using the parameters that yielded the best result.
 
## Parameters

Parameters for the command *pre\_process*:
 
-  *\-n, \-\-data\_name* \<name\>  <x>      :   the data collection name for verbosity
-  *\-m, \-\-metadata* \<file\> <x>         :    path to the metadata table CSV file
-  *\-d, \-\-raw\_data\_path* \<path\> <x>   :    path to the folder containing the input LC-MS/MS raw spectra data files (mzXML format is recommended)
-  *\-y, \-\-processed\_data\_name* [x]       :       The name of the output directory that should be created inside the `raw_data_path` to store the processed data. When not using the default directory, this value should be informed in the following *run* or *clustering* commands. (default: processed_data)
-  *\-t, \-\-rt\_tolerance* [x]       :       the tolerance in seconds used to enlarge the MS1 peaks boundaries and accept as a match all MS2 ions that have a retention time value that is within a MS1 peak range. This value is applied to both sides of the MS1 peaks (RTmin - rt_tolerance and RTmax + rt_tolerance). Tries to overcome bad MS1 peak integrations. (default: 3)
-  *\-z, \-\-mz_tolerance* [x]       :      the tolerance in Daltons for matching a MS1 peak m/z with a MS2 spectrum precursor m/z. (default: 0.05)
-  *\-p, \-\-ppm_tolerance* [x]       :     the maximal tolerated m/z deviation in consecutive MS1 scans in parts per million (ppm) for the initial ROI definition of the R::xcms::centWave algorithm. Typically set to a generous multiple of the mass accuracy of the mass spectrometer. (default: 15)
-  *\-a, \-\-ion_mode* [x]              :    the precursor ion mode. One of the following numeric values corresponding to a ion adduct type: '1' = Positive [M+H]+ or '2' = Negative [M-H]- (default: 1)
-  *\-e, \-\-peak_width* [X,Y]      :       two numeric separated by comma ',' without spaces and using decimal point equals dot '.', containing the expected approximate peak width in chromatographic space. Given as a range (min, max) in seconds. The mean value will be used to simulate the width of the fake peaks (see Details section bellow) (default: 2,10)
-  *\-s, \-\-snthresh* [x]  :           a numeric defining the signal to noise ratio cutoff. (default: 0)
-  *\-f, \-\-pre_filter* [k,I]  :       two numeric separated by comma ',' without spaces and with decimal point equals dot '.' specifying the prefilter step for the first analysis step (ROI detection) of the R::xcms::centWave algorithm. Mass traces are only retained if they contain at least k peaks with intensity >= I. (default: 1,750)
-  *\-r, \-\-noise* [x]  :  a numeric allowing to set a minimum intensity required for centroids to be considered in the first analysis step (centroids with intensity < noise are omitted from ROI detection). (default: 500)
-  *\-c, \-\-mz_center_fun* [name]  :   Name of the function to calculate the m/z center of the chromatographic peak. Allowed are: "wMean": intensity weighted mean of the peak's m/z values, "mean": mean of the peak's m/z values, "apex": use the m/z value at the peak apex, "wMeanApex3": intensity weighted mean of the m/z value at the peak apex and the m/z values left and right of it and "meanApex3": mean of the m/z value of the peak apex and the m/z values left and right of it. (default: wMeanApex3)
-  *\-u, \-\-integrate_method* [x]  :   Integration method for the XCMS CentWave algorithm. For integrate = 1 peak limits are found through descent on the mexican hat filtered data, for integrate = 2 the descent is done on the real data. The latter method is more accurate but prone to noise, while the former is more robust, but less exact. (default: 2)
-  *\-g, \-\-fit_gauss* [x]  :          a logical "TRUE" or "FALSE" indicating whether or not a Gaussian should be fitted to each peak. Decreases performance. (default: FALSE)
-  *\-j, \-\-max_samples_batch_align* [x]  The maximum number of not blank samples to be selected in the metadata sequence for the batch alignment. Due to memory issues this value should not exceed 15 samples to avoid crashing the script. If x == 0 disable the alignment. (default: 5)
-  *\-i, \-\-min_fraction* [x]  :       a numeric defining the minimum fraction of samples in at least one sample group in which the peaks have to be present to be considered as a peak group (feature). Used in the alignment process to define hook peaks. (default: 0.3)
-  *\-w, \-\-bw* [x]  :     a numeric defining the bandwidth (standard deviation of the smoothing kernel) to be used. This option is passed to the density method used in the alignment process. (default: 2)
-  *\-b, \-\-bin_size* [x]  :           a numeric defining the size of the overlapping slices in mz dimension. This option is passed to the density method used in the alignment process. (default: 0.05)
-  *\-x, \-\-max_features* [x]  :       a numeric with the maximum number of peak groups to be identified in a single mz slice. This option is passed to the density method used in the alignment process. (default: 100)
-  *\-q, \-\-processed_data_overwrite* [x]  A logical "TRUE" or "FALSE" indicating if the pre-processed data present in the processed_data_name folder should be overwritten and pre-processed again in case it already exists. Otherwise a previous result may be used for a partial pre processing. (default: "FALSE")
 
-  *\-v, \-\-verbose* [x]  :            for values x>0 show the script output information (default: 0)
-  *\-h, \-\-help*                         output usage information
 
## Results
 
A folder named `processed_data_name` inside the `raw_data_path` containing:
 
  - One MGF per sample of the metadata table with the MS2 spectra enriched with their matched MS1 peak chromatographic dimensions. Each pre-processed file is named as "\<SAMPLE_CODE\>_peak_info.mgf", where \<SAMPLE_CODE\> is as defined in the metadata file located in the *metadata_path*.
  - A count file named "MS1_list_with_MS2.csv" with the quantification by sample of all the MS1 peaks that were assigned to a MS2 ion.
  - A count file named "MS1_list_no_MS2.csv" with the quantification by sample of all the MS1 peaks that were not assigned to a MS2 ion, e.g. probably not fragmented MS1 peaks, and thus these peaks will be missing in the final clustering counts. It will help to search for the ion's isotopic distributions.
  - A log file named "logPreProcessStatisticsWarning" or "logPreProcessStatistics" with the statistics of the rate of MS2 spectra without a MS1 peak correspondence and guidelines to improve this step result.
  - A table file named "log_MS2_no_MS1peak_match.csv" with the informations of the MS2 spectra that did not have a MS1 peak correspondence.
  - If `max_samples_batch_align` > 0, a CSV file named "samples_alignment.csv" is created with the maximum misalignment value in seconds for each data collection batch, as defined in the metadata table, and for all samples together. These values serve as a suggestion for the retention time tolerance to be used in the following steps of the workflow. This alignment process does not change the samples retention time. 
  - A file named "parameters.csv" with the parameter's values used, for reproducibility.
 
## Examples
 
Pre-processing with the default parameters values and compute the misalignment using at most 6 samples per data collection batch:
 
```{ .text .copy }
node np3_workflow pre_process --data_name "data_UHPLC_qTOF" --metadata "/path/to/the/metadata/file/test_np3_metadata.csv" 
--raw_data_path "/path/to/the/raw/data/dir" --max_samples_batch_align 6 --verbose 1
```
 
Pre-processing without computing the misalignment and with a different m/z tolerance value:
 
```{ .text .copy }
node np3_workflow pre_process --data_name "data_UHPLC_qTOF" -m 
"/path/to/the/metadata/file/test_np3_metadata.csv" -d "/path/to/the/raw/data/dir" -j 0 -z 0.01
```
 