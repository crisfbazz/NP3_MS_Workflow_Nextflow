## Step 4 Procedure Details
 
The spectra and peak area counts are performed by scanning and parsing the information present in the .clust files 
from the clustering output to compute the consensus spectra quantification by sample. It aggregates the information 
of the clusters members to obtain the total counts (number of spectra and peak area) by sample and the consensus peak 
profile information. A unique cluster ID is given to the consensus spectra, named msclusterIDs.
 
In the data and in the blank clustering steps the clusters members are the raw data SCANS, what gives for each cluster 
the number of spectra (SCANS's) and the peak profile information by sample. In the following integration steps the 
clusters members are the consensus spectra created in previous steps and the counts only merge and aggregate the values 
computed previously from the raw data.
 
The count of spectra quantification is obtained by summing for each cluster the number of MS2 SCANS detected in each 
sample, what gives the total number of times that a spectrum within a m/z tolerance was detected in a retention time 
range by sample. The count of peak area quantification is quite different, it is based on the peak's IDs of the cluster's 
members. The peak ID is used to detect the sample of origin and to parse the pre-processed data to extract the 
respective peak area, and then for each unique peak ID of the clusters members the count is obtained by summing the 
peak areas by sample.
 
The spectra and peak area count of all the blank samples (column SAMPLE_TYPE equals 'blank' in the metadata table) are 
summed for each m/z (row) and the result stored in a new column named *BLANKS\_TOTAL*. The same is performed for control 
and bed samples (SAMPLE_TYPE equals "control" and "bed", respectively), resulting in two new columns named 
*CONTROLS\_TOTAL* and *BEDS\_TOTAL*, respectively. These columns containing the total counts by sample type can help 
filter false positive candidates and are used to compute the flag indicators described below.
 
The counting script adds at most three new boolean (TRUE or FALSE) columns named *BFLAG*, *CFLAG* and *BEDFLAG*, if 
the respective sample type equals blank, control or bed is present in the metadata table. They are computed as follows: 

- if a msclusterID have the column *BLANKS\_TOTAL*, *CONTROLS\_TOTAL* or *BEDS\_TOTAL* greater than zero all spectra that 
have a m/z's within the given mass tolerance from it's m/z will have the column *BFLAG*, *CFLAG* or *BEDFLAG*, respectively, set to TRUE, otherwise it will be FALSE. 
- These column flags serve as a warning for false positive msclusterIDs, they signalize that there is a very close m/z that was detected in a blank, control or bed sample respectively, and is up to the user to check its correctness. 
 
The 'hit' sample type is used to dereplicate the m/z's that appeared in one of the samples that had a high 
bioactivity score or that has a m/z close enough to other m/z that appears in one of these samples. 
The counting scripts adds two new character columns named *DESREPLICATION* and *HFLAG*. 
The first stores all the hit samples codes in which each msclusterID have a count greater than zero. 
And the second column signalize all the hit samples codes that have a count greater than zero in any cluster ID with 
a m/z within the given mass tolerance of each msclusterID m/z (this column contains the first column results). 
The same way as the other flag columns, the *HFLAG* serve as a warning for false negative clusters, evidencing a 
spectrum that has a m/z very similar to a spectrum that appeared in a hit sample and is up to the user to check its 
correctness.
 
#### NP³ MS workflow Clustering Count Table Format

The resulting clustering count table has *m* rows and at most 20 + *n* columns, where *m* equals the number of 
consensus spectra in the end of the clustering job and *n* equals the number of processed samples present in the 
metadata file. A description of each column is presented below (as described above the weighted mean values are 
obtained using the cluster members MS2 summed intensities as weights).

|           Columns           |                                                                                                       Description                                                                                                        | Value Type |
|:---------------------------:|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:| :-------: |
|         msclusterID         |                                                                                         the consensus spectrum unique cluster ID                                                                                         | numeric |
|         numSpectra          |                                                                              total number of joined spectra (cluster members), cluster size                                                                              | numeric |
|         mzConsensus         |                                                                                         weighted mean of the cluster members m/z                                                                                         | numeric |
|           rtMean            |                                                                                weighted mean of the cluster members retention time center                                                                                | numeric |
|            rtMin            |                                                                               weighted mean of the cluster members retention time minimum                                                                                | numeric |
|            rtMax            |                                                                               weighted mean of the cluster members retention time maximum                                                                                | numeric |
|           peakIds           |                                                                           the peak IDs of the cluster members concatenated without duplicates                                                                            | string |
|            scans            |                                                               the MGF SCANS number with the sample of origin as suffix of the cluster members concatenated                                                               | string |
|         basePeakInt         |                                                                               the maximum base peak intensity value of the cluster members MS2                                                                              | numeric |
|           sumInts           |                                                                                        sum of the cluster members MS2 intensities                                                                                        | numeric |
| <X\_1\>\_<area\|spectra\> |            the count of area or spectra in the sample X_1, where 'spectra' is the number of spectra or 'area' is the unique sum of peak areas in the sample with SAMPLE_CODE equals X_1 in the metadata table            | numeric |
|             ...             |                                                                                                                                                                                                                          | |
|      <X\_*n*\>\_<area\|spectra\>      |                                                                                       the other samples counts for this mslusterID                                                                                       | numeric |
|        BLANKS_TOTAL         |                                                                           total number of spectra or peak area that appeared in blank samples                                                                            | numeric |
|       CONTROLS_TOTAL        |                                                                          total number of spectra or peak area that appeared in control samples                                                                           | numeric  |
|         BEDS_TOTAL          |                                                                            total number of spectra or peak area that appeared in bed samples                                                                             | numeric  |
|       DESREPLICATION        |                                                              concatenation of all the hit samples codes in which the spectra had a count greater than zero                                                               | string |
|            BFLAG            |                                                                 TRUE or FALSE indicating if there is a close m/z to the mzConsensus in any blank sample                                                                  | boolean |
|            CFLAG            |                                                                TRUE or FALSE indicating if there is a close m/z to the mzConsensus in any control sample                                                                 | boolean |
|           BEDFLAG           |                                                                  TRUE or FALSE indicating if there is a close m/z to the mzConsensus in any bed sample                                                                   | boolean |
|            HFLAG            | concatenation of all the hit samples codes that have a count greater than zero in any msclusterID with a mzConsensus within the given mass tolerance of the current cluster (contains the DESREPLICATION column results) | string |
|          peaksList          |                                                                the concatenated list of the fragmented peaks m/z's of the msclusterID consensus spectrum                                                                 | string |
|          peaksInt           |                                                             the concatenated list of the fragmented peaks intensities of the msclusterID consensus spectrum                                                              | string |

#### MS2 Base Peak Intensity Distribution 

At the end of the clustering step, the distribution of the MS2 base peak intensity values (column basePeakInt) is plotted 
and saved to a file named 'basePeakInt_distribution.png'. 

A descriptive analysis of the basePeakInt distribution is also computed and saved to a file named 
'basePeakInt_distribution_summary.txt'. 

If there is at least one blank sample in the metadata table, the base peak distribution plot will be colored to 
highlight the values from the consensus spectra that appear in blank samples and the descriptive analysis will be 
computed only for these consensus spectra from blank samples. 

The interquartile range (IQR) is also computed and printed at the end of the descriptive analysis file to help in 
the setup of the noise cutoff and of the bflag cutoff parameters (to be used in the clean Step 5). 
Also, vertical lines corresponding to different noise/bflag cutoff values are plotted in the base peak intensity 
distribution.