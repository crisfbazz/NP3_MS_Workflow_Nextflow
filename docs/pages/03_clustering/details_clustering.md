## Step 3 Procedure Details

Relying on the principle of similarity between experimental spectra[^1][^2] the samples clustering is performed in 
steps, based on the defined data collection batches and the samples types, in the metadata table. Thus, the more 
similar spectra, guaranteed of being detected in the same conditions, are forced to be joined previous and more 
information be aggregated to them before clustering all the samples spectra together.
 
The NP³ MS workflow spectra clustering is performed in four steps:
 
1. The **data clustering step** - all the not blank samples (column SAMPLE_TYPE equals 'sample', 'bed', 'control' or 
'hit' in the metadata table) of each data collection batch are clustered together; 
2. The **blank clustering step** - all the blank samples (column SAMPLE_TYPE equals 'blank') of each data collection 
batch are clustered together, to try to avoid the noise impact of these samples in the remaining clustering steps; 
3. The **data collection batch integration step** - the results of the data and the blank clustering steps of each data 
collection batch are clustered together; 
4. And the **final integration step** - the results of all data collection batch integration steps are clustered together 
to integrate all results.
 
The NP3_MSCluster results are the MGF with the collection of consensus spectra from the outputted clusters and .clust 
files (in CSV format) with the cluster's information describing which scans were joined to make up the consensus spectra. 
In the data and in the blank clustering steps the resulting .clust files contains the spectra detected in the raw files 
(real SCANS), with the matched peak information, and in the subsequent clustering steps the SCANS are the clusters IDs 
(named msclusterID) of the consensus spectra created in the previous steps by the workflow Step 4. 
The info present in the .clust files may be used for debugging.

## Step 3 Methodology Details
 
Tandem mass spectrometry (MS/MS) experiments often generate redundant data sets containing multiple spectra of the same 
molecule. Clustering of MS/MS spectra takes advantage of this redundancy by identifying multiple spectra of the same ion 
and replacing them with a single representative spectrum[^3].
 
The NP3_MSCluster is a modified version of the MS\-Cluster algorithm[^1]. The MS\-Cluster is a simple and effective 
algorithm designed to rapidly process large MS/MS data sets (even in the excess of 10 million spectra), while ensuring 
the high quality of the resulting clusters. Its optimizations were made for proteomics focusing on processing peptides 
data. This open-source software was developed in C++ and is available for download at http://peptide.ucsd.edu. 
It is part of the GNPS workflow[^4].
 
To use the MS\-Cluster in the NP³ MS workflow the bias toward proteomics had to be overcome with modifications toward 
metabolites processing. The major limitation of its use is that stereo isomers are not resolved as separate clusters, 
because the original algorithm did not use the retention time for joining the spectra. Another limitation is that this 
algorithm is not an optimal solution, it approximates a hierarchical clustering leading to fragmented clusters, that is, 
several distinct clusters containing spectra of a same ion. These fragmented clusters can lead to a wrong count of the 
detected ions from the same metabolite and the solution for this problem is presented in the section of Step 5 clean.
 
An ion in a LC\-MS/MS experiment is detected in a chromatographic peak (MS1 peak), which is characterized by a m/z range, 
a retention time range and the total intensity (peak area). The MS2 data are the fragmented ions (spectra) of the MS1 peaks. 
The *pre_process* step is capable of matching the peak dimension information of the MS1 peaks with all detected MS2 spectra, 
enriching the spectra with chromatographic information and making possible to use this information in the clustering job.
 
The NP3_MSCluster algorithm was developed to support the clustering of metabolites and to separate stereo isomers. 
It parses the MS1 peak dimensions information from the NP³ MS workflow pre-processed MS/MS data and incorporates it in 
the cluster's attributes, by taking the average of its members peak dimensions. With the use of this additional 
information, it prevents the join of spectra that have the same m/z but different retention times (descendant from 
different MS1 peaks, not concurrent), what characterizes stereo isomers. 
 
A tolerance in the spectra retention time range is used to deal with the calibration problem when multiple MS runs are 
being clustered and the runs were not carefully aligned or can't fully rely on the alignment. Due to the retention time 
tolerance very close by peaks can still be joined and the user must take this into account when setting this tolerance 
value. The default parameters values were refined for a UHPLC-MS/MS-qTOF equipment, the UHPLC Acquity HClass Waters and 
the spectrometer ESI\-QqTOF Impact II Bruker.
 
The NP3_MSCluster algorithm only considerate two clusters within a given mass tolerance as joining candidates if the 
retention time mean (center of the peak) of one cluster is contained in the retention time range of the other cluster, 
within a retention time tolerance applied in the peak boundaries. Then the clustering happens as in the original 
MS\-Cluster algorithm. If the cosine similarity of the joining candidates is above the threshold, they will be joined 
in a consensus spectrum, called the cluster representative. The consensus spectrum peak dimensions attributes are 
computed by summing the total MS2 intensities of the cluster members and using it to calculate the weighted average of 
the members peak dimensions values, with the precursor intensity of each cluster member as weight. 
At the end of the clustering process each consensus cluster will represent the spectrum of one MS1 peak, except for 
the fragmented clusters limitation and the possible join of very close by peaks (very close by isomers).
 
The MS\-Cluster configurations and some hard coded parameters had to be modified to eliminate its bias toward proteomics. 
The filtering (*sqs*) and the annotation parameters (*assign\-charges* and *correct\-pm*) were disabled from the default 
settings following the authors recommendation[^5]. To reduce the loss of information from the detected MS2 spectra, 
the hardcoded filter of spectra with too few numbers of fragmented peaks was modified in the entire code to only filter 
spectra with no fragmented peaks. The original algorithm filtered spectra with less than 7 peaks in the input and with 
less than 15 peaks in the output. A parameter was added to allow the user to control this filtering in the final 
clustering step, by choosing the minimum number of fragmented peaks that a spectrum must have to be outputted 
(default to 5). Some other hard coded heuristic parameters were relaxed to allow more comparisons, aiming to reduce 
the fragmented clusters occurrence with exchange in performance (see file NP3_MSCluster/src/MsCluster/MsClusterIncludes.h 
for the list of changed parameters).
 
Scaling peak intensities has been shown to improve the quality of the similarity computations[^2]. 
The scaling method adopted by the MS\-Cluster algorithm is the natural logarithm of the peaks' intensities, which they 
found to be the most suitable for their data[^1]. But for our metabolomic data we found that the square root scaling 
method was the most suitable one and it is the recommended one for dot-product algorithms (cosine)[^2][^6]. 
To overcome this limitation a new parameter was added to allow the user to choose which scaling method is to be used 
in the entire workflow, and the available options are: no scaling; natural logarithm scaling; and scaling to the power 
of x, where x can assume any value greater than zero (by default set to square root scaling x = 0.5). 

## References
[^1]: Frank AM, Bandeira N, Shen Z, Tanner S, Briggs SP, Smith RD, Pevzner PA. Clustering millions of tandem mass spectra. 
J Proteome Res. 2008 Jan;7(1):113-22. doi: 10.1021/pr070361e. Epub 2007 Dec 8. PMID: 18067247; PMCID: PMC2533155. (MS-Clsuter)
[^2]: Stein SE, Scott DR. Optimization and testing of mass spectral library search algorithms for compound 
identification. J Am Soc Mass Spectrom. 1994 Sep;5(9):859-66. doi: 10.1016/1044-0305(94)87009-8. PMID: 24222034.
[^3]: Ramsay, D.. “Applications of Clustering Algorithms in the Analysis of Mass Spectrometry Data.” (2017).
[^4]: Wang, M., Carver, J., Phelan, V. et al. Sharing and community curation of mass spectrometry data with Global 
Natural Products Social Molecular Networking. Nat Biotechnol 34, 828–837 (2016). https://doi.org/10.1038/nbt.3597 (GNPS).
[^5]: Watrous, J., Roach, P., Alexandrov, T., Heath, B., Yang, J., Kersten, R., Voort, M., Pogliano, K., Gross, H., 
Raaijmakers, J., Moore, B., Laskin, J., Bandeira, N., & Dorrestein, P. (2012). Mass spectral molecular networking of 
living microbial colonies. Proceedings of the National Academy of Sciences, 109(26), E1743–E1752.
[^6]: van den Berg, R.A., Hoefsloot, H.C., Westerhuis, J.A. et al. Centering, scaling, and transformations: 
improving the biological information content of metabolomics data. BMC Genomics 7, 142 (2006). 
https://doi.org/10.1186/1471-2164-7-142