## Step 5 Details
 
The clustering (Step 3) identify multiple spectra of the same ion and replace them with a single representative 
consensus spectrum. The MS-Cluster algorithm is not an optimal solution, so even though it drastically reduces the 
data size keeping it's quality high, it is not able to remove all the repeated spectra of the same ion, leading to 
fragmented clusters. In the NP3_MSCluster algorithm, fragmented clusters are clusters with the same m/z within a m/z 
tolerance and that coeluted in a concurrent chromatographic peak (detected in the same retention time range) within a 
retention time tolerance.
 
The clean step is a second clustering strategy that was developed to reduce the remaining redundancies from the 
clustering result, and to overcome the fragmented clusters limitation. It is an optimal clustering with greedy 
heuristics that was implemented in R code with auxiliary compiled functions (dlls) in C++. 
 
The *clean* step starts by performing the pairwise comparisons of the consensus spectra. These comparisons result in 
two symmetric and quadratic matrix (n + 1)x(n + 1), where n is the number of consensus spectra from Step 3 and the 
remaining 1 row and column are the msclusterIDs (clusters IDs resulting from the clustering step). The upper triangular 
part of these matrix contains in one the similarity values and in the other the number of matched peaks between all 
consensus spectra. And their lower triangular part is empty to reduce repeating the same computations. These similarity 
values are computed using a shifted version of the cosine function [8,9], which we call NP³ shifted cosine function, 
or spec2vec [24], using the model trained on UniqueInchikey subset (12,797 spectra). If spec2vec is used, the 
matchms [25] library is used to compute the number of matched peaks between the compared spectra.
 
The NP³ shifted cosine algorithm first matches the same m/z's from fragmented peaks list of the spectra being 
compared within a mass tolerance; then, for the remaining not matched fragmented peaks of the spectra it adds a shift 
equals to the difference of the spectra precursor m/z in all fragmented peaks' m/z of the spectrum with the smaller 
precursor m/z. And finally, it tries to match again the remaining peaks of the spectra within the mass tolerance. 
The shift is used to overcome a small motif difference that could be present in the ions structure, it will produce 
a shift in some fragmented peaks m/z's equal to the motif mass. Only motifs greater than the m/z tolerance and smaller 
than the max_shift parameter value are considered, a precursor m/z difference outside this range will be set to zero, 
and then, the shift won't be applied. Using the precursor m/z difference in the cosine computation is a heuristic that 
can increase the similarity score of analogue structures and possible ionization variants (e.g., different adducts, 
neutral losses, etc). This implementation was inspired by the modified cosine implemented in the GNPS workflow [8].
 
Next, the *clean* step use the computed similarity table of the resulting consensus spectra to reduce the clustering 
count tables. It works in the following steps: 
 
1. Data cleaning step: the fragmented clusters are joined if their consensus spectra have a similarity value above the 
cutoff, or if they share at least one peakId or if they pass the bflag cutoff criteria. The similarity table is reduced 
based on the performed joins always keeping the maximum similarity value between the aggregated rows and columns (greedy 
heuristic). This procedure is repeated at most 10 times to take into account the results of the performed aggregations 
of previous steps, it stops earlier if no join is performed. At the first time of this step, only blank spectra are 
considered for the joining and teh bflag cutoff is applied. Finally, a noise cutoff is applied to remove spectra with 
a low base peak intensity value;
2. Indicators step: the samples types indicators are computed for each row of the just obtained clean counts. 
At the end, all the indicators of Step 4 are computed plus another indicator showing the distance of each spectra to a 
blank spectrum (when blank samples are present). Additionally, a fragmented cluster indicator is computed to check the 
remaining not joined clean consensus spectra. The joined clusters have their fragmented peak list and intensities 
aggregated in a consensus that is stored in the clean count table, and saved to the clean MGF file. 
The IDs of the joined clusters are stored in the clean count tables as well, in a column named 'joinedIds'; 
 
In the Data cleaning step, when searching for fragmented clusters an additional condition was added in the criteria to 
join two spectra. It checks if the peak center deviation between the spectra is less than 4 times the retention time 
tolerance or if their peak boundaries deviation is less than 2 times the retention time tolerance. 
If both of these conditions are not satisfied, the spectra are kept separated. By doing this, we try to prevent joining 
adjacent isomers that have a high similarity value. 

Also in the Data cleaning step, a bflag cutoff was implemented to allow joining spectra from blank samples that do not 
have a similarity value above the cutoff. This helps to reduce the fragmented clusters with bflag TRUE and a low base 
peak intensity, that could not fully rely in the similarity values. The bflag cutoff is computed as the median value of 
the base peak intensity distribution of blank spectra plus the factor informed by the user times the interquartile 
range (IQR) of this distribution. The IQR is the range between the 1st quartile (25th quantile) and the 3rd quartile 
(75th quantile) of a distribution. The consensus spectra with a basePeakInt value <= median + IQR*bflag_factor 
(from the blank spectra basePeakInt distribution) and BFLAG TRUE will be joined to a blank spectrum independent of 
their similarity values.

At the end of the clean step a noise cutoff is applied to remove spectra with a low base peak intensity value. 
The noise cutoff is computed as the the median value of the base peak intensity distribution of blank spectra plus 
the factor informed by the user times the interquartile range (IQR) of this distribution. When no blank sample is 
present in the metadata, the full distribution of the base peak intensity from the clustering counts is used. 
This cutoff will affect the consensus spectra with a low basePeakInt value that probably are noise features. 
If the clustering Step 3 resulted in more than 15000 consensus spectra, the noise cutoff will be applied before the 
clean step to prevent a long processing time.

The base peak intensity distribution plotted at the end of the clustering step helps the users to better define these 
cutoffs' factors and to better understand its effect on their dataset.
 
The count tables, in terms of the number of spectra and of peak area, of the joined clusters are aggregated following 
the same rules applied in Step 4 and the joined clusters peak information (m/z, retention time mean, minimum and maximum) 
are computed as the average of the joining clusters values weighted by their intensities (sumInts).
 
The following columns are added to the count tables in the *clean* step:
 
**Table: Clean Count Tables New Columns from the Data Cleaning Step**

| Columns | Description | Value Type |
| :--------------- | --------------------------------------- | :-------: |
| joinedIDs | the msclusterIDs of the joined spectra in the cleaning step, separated by a ';' | character with concatenated numbers |
| numJoins | the number of spectra that were joined in the cleaning step | numeric |
| BLANK_DIST | the distance in the molecular network of similarity from the current spectra to a spectrum of a blank sample. This distance ranges from 0 (if the current spectra appear in a blank sample) to 3 (if there are at least 3 links between the current spectra and a spectrum of a blank sample in the molecular network of similarity). If this value is NA it means that this distance is at least greater than 3, and thus was not computed. | numeric |
| fragmented_cluster | the number of clean consensus spectra identified as a fragmented cluster of the current spectrum or -1 to indicate that the current spectrum is a fragmented cluster of another more intense spectrum. The clean consensus spectrum that receives a positive value greater or equal than one in this column is the candidate with the highest basePeakInt among the other spectra signalized as fragmented clusters within the same m/z and retention time interval, which receives a -1 value. The spectrum without fragmented clusters receives a value equals to 0. | numeric |

