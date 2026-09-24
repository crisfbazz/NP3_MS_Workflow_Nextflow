## *pre_process* Optimization Guide


The Pre-processing step performs the LC-MS data processing (feature finding) with XCMS. It extracts the list of MS¹ 
peaks for each sample, and then, directly impact in the isomers separation and resolution. The LC-MS processing is 
highly dependent on the parameters setup. The optimization of these parameters depends on the mass spectrometer equipment, 
the chromatographic method used and the data quality. With that pointed out, a diagnostic of Pre-process result and a 
suggestion of values for the **most critical parameters** were implemented to help the user in this optimization. 
The diagnostic and suggestion are still in a **Beta version**, and should be used as a guide and not as a rule.
 
In order to optimize the Pre-processing step, the user could follow two strategies: (1) run pre-processing with the 
default parameters and analyse the diagnostic result and suggestions; or (2) inspect the chromatograms of the dataset 
first (using a target m/z or randomly pick one) and run the pre-processing with a pre setup of parameters values that 
best represent what you observed in the real data. **The second option is more recommended**, because the suggestion of 
values for the most critical parameters will be biased by the parameters values used. The user can also choose to run 
the entire workflow with this first pre-processing result to examine the final clean quantification tables and decide if 
they are good enough for its purpose.
 
The critical parameters for the pre-processing step are: 

1. `ppm_tolerance` - the maximal tolerated m/z deviation in consecutive MS1 scans in parts per million (ppm). 
Typically set to a generous multiple of the mass accuracy of the mass spectrometer; 
2. `mz_tolerance` - the tolerance in Daltons for matching the m/z of a MS1 peak with the precursor m/z of a MS2 spectrum; 
3. `rt_tolerance` - the tolerance in seconds used to enlarge the MS1 peaks boundaries and accept as a match all MS2 
ions that have a retention time value within a listed MS1 peak range;
4. `peak_width` - the expected approximate peak width in chromatographic space (MS1 peaks). Given as a range (min, max) 
in seconds.
 
The first two critical parameters are more related to the mass spectrometrer equipment accuracy and the data resolution, 
 while the third and fourth critical parameters are also related to the chromatographic method and the data quality (well defined MS1 peaks). Then, to better define the last two critical parameters, the user must know your data and/or inspect its chromatograms to visualize the MS1 peaks' shape and definition. The MS1 peaks' shape and definition can vary between different datasets coming from the same equipment and can also vary between different m/zs from the same dataset. Because of that, the **last two parameters are the most critical ones**, with the `peak_width` being the most critical of all.

It is also important to notice that the noise removal parameters of Pre-processing, used to remove minority peaks, 
have very small default values close to the baseline of evaluated samples. For different samples they could also 
compromise the correspondence between MS¹ and MS². By default, Pre-processing set a minimum intensity equal to 500 to 
consider a MS¹ peak (parameter **noise**) and only consider MS1 peaks that have at least one ion with 750 of intensity 
(parameter **pre_filter**). If the minimum required MS¹ intensities are higher than the samples MS¹ baseline, real MS1 
peaks may be removed and the rate of non-correspondence between MS¹/MS² may increase (no matches in terms of m/z can be 
greatly affected here). The user must properly set the noise removal parameters values in Pre-processing to prevent this 
from happening (the minimum experimental intensity used to detect MS² ions is a good upper limit) and to better optimize 
the other critical parameters.

The user must take the following relations into account when choosing the critical parameters values: 
 
- ppm tolerance:
    - A small value can split MS1 peaks from the same ion into distinct peaks, with close but different m/zs
    - A big value can join MS1 peaks from different ions with different m/zs
- Retention time tolerance and m/z tolerance:  
    - A small value can prevent the match between a MS2 spectrum and its respective MS1 peak 
    - A big value can assign a MS2 spectrum to a wrong MS1 peak (this can hide a bad MS1 integration processing)
- Minimum peak width:  
    - A small value can split short MS1 peaks 
    - A big value can exclude very short peaks and/or join adjacent isomers 
- Maximum peak width:  
    - A small value can split large MS1 peaks 
    - A big value can join adjacent isomers 
  
As you can see by now, this optimization will be a balance between the amount of effort that the user is willing to 
give and how good is the MS1 peaks definition across the different m/zs of a dataset. 

At the end of the pre-processing step a diagnostic of its results is performed and stored in a log file that will be 
named 'logPreProcessStatisticsWarning' if a problem is detected in the diagnostic, or will be simple named 
'logPreProcessStatistics' if no problem is detected. This log is saved to the directory named with the 
'processed_data_name' inside the 'raw_data_path' folder, e.g. the pre-process output folder. 

The pre-process statistics log will start with the statistics of the percentage of MS2 spectra without a MS1 peak 
correspondence. This is presented using two variables, better described in the last section 4.3.3., which are: 

1. **rate\_MS2\_no\_MS1\_mz** - the percentage of MS2 spectra without a MS1 m/z range correspondence; and 
2. **rate\_MS2\_no\_MS1\_rt** - the percentage of MS2 spectra without a MS1 retention time range correspondence. 
These rates are computed for each sample and removing blank m/zs.
 
If any of the not blank samples have one of these rates above a cutoff of 5%, a report is printted informing the top 5 
problematic m/zs by sample. In these case the user should visualize the chromatogram of one or all (recommended) the 
printed m/zs to identify the problem(s). These could include bad defined peaks, for example jagged peaks or peaks with 
a long tail, or peaks with a width different from the values used in the `peak_width` parameter. The $NP^{3}$ MS 
workflow command *chr* can be used for this visualization or any other LC-MS/MS visualization tool. And then, the user 
have to decide to repeat the pre-process step with other parameters values or to accept this result. 
A descriptive statistics of these rates (mean, median and quartiles values) is also printted to help the user in this 
decision. When the job have a big number of samples it can be more reasonable to pay more attention to the mean/median 
values of these rates instead of expecting all samples to be below the cutoff (the best parameters for a given sample 
can give bad results in another sample). 

At the end of the diagnostic of MS1 and MS2 correspondence, a diagnostic to check if real peaks are being split is also 
performed. It will count the number of MS1 peaks (possible isomers) by m/z for each sample and aggregate this counts 
into a single distributions. Then, if more than 5% of the m/zs have more than 30 MS1 peaks a warnings will be printed 
indicating that real peaks could have being split into multiple smaller peaks. The user must evaluate if this is not 
expected for the dataset being used, and then choose to increase the maximum peak_width value to prevent this anomaly. 
The distribution of the number of m/zs by the number of MS1 peaks is stored in the file 'log_number_mzs_by_number_peaks.png' 
for further examination.

Following the diagnostic, a Beta suggestion for the `peak_width` and the `rt_tolerance` parameters values will be 
executed to help the user in this optimization. It starts by plotting the peak width distribution of the pre-process 
result, for all the MS1 peaks that had a MS2 correspondence and excluding blank m/zs. These distribution will show the 
common peak width returned by the LC-MS processing with the given parameters, then its important to notice that it will 
be biased by the parameters used, but it will also tend to go to the expected peak width for that dataset. For example, 
if the pre-process was executed with big values in the `peak_width` parameter, it will result in large MS1 peaks being 
detected even if that does not represent what is found in the dataset, because it will join adjacent peaks to fulfill 
the given `peak_width` values. And that's why we recommend the user to start with a pre setup of parameters that better 
represent your dataset.

Depending on the diagnostic result, the user should choose to follow or not the suggested values and rerun the pre-process 
if needed. The peak width distribution of the pre-process result is stored in the file 
'MS1_list_with_MS2_noBlank_peak_width_hist.png' to help the user in this decision.

The user can also look at the table with the list of detected MS1 peaks that had a MS2 correspondence to better 
understand the pre-process results. This table is named 'MS1_list_with_MS2.csv' and is located in the output folder of 
the pre-process step. When there is a target m/z, it is easier to look for it in this table to check if the results 
follow what is seeing in the chromatograms. Otherwise, the user could randomly choose a m/z to do it or use the 
problematic m/zs that appear in the diagnostic when a warning is emitted.
