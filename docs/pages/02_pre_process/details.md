## *pre_process* Details
 
The pre-process step use the R package XCMS CentWave algorithm [1] to process the LC-MS data and extract the samples' 
list of MS1 peaks, the XCMS::PeakGroup algorithm [2] to measure the disalignment between the samples, the MSnbase R 
package [3] to read the raw MS/MS data and a modified version of this package function to write to MGF files the 
pre-processed MS/MS data enriched with the MS1 peak dimensions information. Four fields are added to the MGF files 
header to store the MS1 information: RTMIN, RTMAX, PEAK_ID and PEAK_AREA.
 
The pre-process will enrich the MS2 spectra with a MS1 peak dimension. To do so, each detected MS2 ion's retention 
time and precursor m/z are matched against its sample list of MS1 peaks. The closest peak that encompasses the spectrum 
is assigned to it. This match reference is used by the write MGF function to write the peak dimensions information 
(retention time minimum and maximum, the peak area and the peak ID) in the header of each MS2 ion of the resulting 
pre-processed MGFs. The peak ID is unique for each MS1 peak, it is given by the CentWave algorithm and has the spectra's 
sample code of origin added as a suffix. The MS2 retention time is kept as the center of the peak.
 
Blank samples (in the metadata table SAMPLE_TYPE equals 'blank') receive a different treatment when enriching the spectra. 
To deal with the baseline problem, e.g., no defined MS1 peak for a given m/z, all MS2 ions of the blank samples are 
assigned with a baseline peak dimension which encompasses the entire chromatogram from retention time 0 to 1000000 seconds. 
These allows for more false positive spectra to be removed from the results.
 
In the matching peak step, depending on the parameter's values used, mostly the *peak_width*, some peaks can be incorrectly 
integrated by the CentWave algorithm generating false positive and false negative peaks. This can lead to no match between a 
MS2 ion and a MS1 peak. The reasons for not finding a peak match can be mostly because of three facts: (1) bad defined values 
for the *peak_width* parameter; (2) "bad" peaks that do not have a good gaussian shape (e.g., jagged peaks, very close 
intersecting isomers peaks) can make the CentWave algorithm not work properly and wrongly integrate the peaks (since it 
is an approximated solution) [4,5], or (3) a bad choice of the *mz_tolerance* or the *rt_tolerance* parameters values 
that could be too small to allow a match. The first and third cases can be modified as needed by the user for each job, 
and the second case depends on the data quality.
 
When no MS1 peak match is found for an MS2 spectra, a "fake" peak is created. The mean peak width, provided in the 
*peak_width* parameter, is used to simulate the fake peak width centered in the MS2 spectra retention time 
(RTmean +- mean(peak_width)/2) and the spectra intensity is used as the fake peak area (the clustering will simulate 
the peak area integration by summing the MS2 spectra intensities). When the peak match occurs because of the 
*rt_tolerance* value, e.g., the MS2 ion retention time is outside the integrated MS1 peak width but within a distance 
less than the *rt_tolerance* value, the peak width for that spectra is enlarged to encompass the spectra retention time 
correctly. When more than one match is found, the closest MS1 peak with respect to the peak retention time mean 
(center of the peak) to the spectra's retention time is chosen as the final match.
 
When no MS2 spectra receives an integrated MS1 peak, this peak dimensions and quantifications by sample are recorded in 
a file named "MS1_list_no_MS2.csv" inside the processed data folder. These records will generate a new count file with 
all MS1 peaks that were not fragmented and thus did not produce a MS2 spectra, and therefore, will be missing in the 
final clustering counts. This count is aggregated after a clustering job to produce the MS1 count file named with the 
suffix "_peak_area_MS1.csv", which will add up for the job completeness to try to guarantee that any ion that is 
detected in the data collection samples will be present in the final counts. These list of not fragmented MS1 peaks 
will also be used to search for isotopic patterns of multiple charged ions as part of the rules used in the annotation 
step (Step 7).
 
When *max_samples_batch_align* > 0, the XCMS::adjustRtime-peakGroups algorithm [2] is used to perform an alignment of 
the samples MS1 peaks and suggest a value to overcome the samples misalignment. The maximum deviation found between the 
sample's retention time is used to give a suggestion for the retention time tolerance value that should be used in the 
following workflow commands (in the clustering, clean and merge step, NOT in the pre-process step). It must be noted 
that the alignment does not change the samples retention time, it only gives a suggestion of the retention time 
misalignment between the samples. 
 
The samples alignment is a very computational consuming step and can possibly not work or crash in some computers with 
limited memory and/or processing capacity. To try to overcome this difficulty three programming decisions were made. 
First, the blank samples (metadata column SAMPLE_TYPE equals 'blank') are not used in the alignment, mainly because 
these samples usually have a lot of noise that could lead to a bad behavior of the algorithm. Second, the samples 
alignment of a same data collection batch can be made using only a given number of samples, defined by the parameter 
*max_samples_batch_align*, chosen following the not blank samples order in the provided metadata table (the first not 
blank sample of each data collection batch appearing in the metadata file is then always picked). And third, the 
alignment between all samples is always made using only one representative sample of each data collection batch, 
which is also hard coded to be the first not blank sample of each data collection batch appearing in the metadata file. 
These behaviors must be taken into account when setting the *max_samples_batch_align* value and when defining the samples 
order in the metadata file. The first sample of each data collection batch in the metadata file should be the more related 
samples between the data collections, e.g. the ones that have the same polarity or bigger chances of having more peaks in 
common. 

The parameters *min_fraction*, *bin_size*, *bw* and *max_features* should also be tuned to achieve a better alignment 
result. If the alignment persists to fail it is up to the user to guess a good value for the retention time tolerances 
of the following commands, and this value should represent the real misalignment between the samples. The alignment can 
fail if not enough common MS1 peaks are found between the samples, this is a limitation of the used algorithm.
 
If the alignment succeeds, then the user may evaluate to run the following $NP^{3}$ MS workflow commands with the 
suggested retention time tolerance values. The user can increase or decrease this tolerance, but we recommend to try to 
keep it in a range that will prevent joining close by peaks (isomers) or splitting a large peak depending on the common 
*peak_width* found in the raw data. It's up to the users to correctly define this tolerance values.

Finally, a diagnostic is performed to evaluate the pre-process result based on the rate of MS2 spectra without a MS1 peak 
correspondence (percentage of fake peaks). And then, a suggestion for better *peak_width* and *rt_tolerance* pre-process 
parameters values is performed to help the user improve this step result when needed.
The diagnostic computes the percentage of MS2 spectra without a MS1 m/z range correspondence (named rate_MS2_no_MS1_mz) 
and the percentage of MS2 spectra without a MS1 retention time range correspondence (named rate_MS2_no_MS1_rt) for each 
sample and removing blank m/z. If any of these rates is above a cutoff of 5% for a not blank sample, a report is 
printed informing the most problematic m/zs by sample to facilitate a manual evaluation of these cases by the user. 
The user should visualize the chromatogram of these m/zs to identify the problem(s), and then, decide to repeat the 
pre-process step with other parameters values or to accept the result. A descriptive statistics of these rates is also 
printed to help the user in this decision, when the job have a big number of samples it can be more reasonable to pay 
more attention to the mean/median values of these rates. Following this diagnostic, a suggestion for better *peak_width* 
and *rt_tolerance* parameters values is performed based on the peak width distribution of the MS1 list with a MS2 
correspondence (not fake peaks) and without blanks (table 'MS1_list_with_MS2.csv' in the pre-process result). 
The suggestion is computed as follow:

 - Retention time tolerance = median / 2
 - Peak width range (minimum and maximum) = Q25 - 1.5 * IQR, Q75 + 1.5 * IQR
   - The IQR (interquartile range) factors start at 1.5 and are decreased by 0.25 until the peak width suggestions are 
   - between the minimum and maximum peak width limit values.
   
Where the median, Q25, Q75 and IQR values corresponds, respectively, to the median, first quartile (percentile 25th), 
third quartile (percentile 75th) and interquartile range values of the peak width distribution of the not fake peaks 
without blanks. Depending on the shape of this distribution, the suggestions can be more or less accurate.

The user should use these suggestions as a guide for choosing better values for the parameters, depending on the 
diagnostic result and also based on your own analysis/knowlegde of the data (e.g. chromatograms inspection) and the 
peak width distribution of the real MS1 peaks.

At the end, a diagnostic of possible splitted peaks is also performed. It works by counting the number of MS1 peaks 
by m/z present in the MS2 data, without blanks. The number of MS1 peaks indicates the number of possible isomers by m/z. 
The distribution of the number of m/zs that had a given number of MS1 peaks is stored in the file 
'log_number_mzs_by_number_peaks.png' located in the pre-process result. This diagnostic will print a warning if there 
are more than 5% of the m/zs with more than 30 MS1 peaks (possible isomers), which can indicate that real peaks are 
being splitted in smaller peaks. The user should evaluate if this is really an issue for the dataset being used.
