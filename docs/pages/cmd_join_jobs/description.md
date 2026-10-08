# Command *join_jobs* 

Command to join different NP³ jobs into a single united job using an incremental clustering approach. 
It concatenates different results from the **run** or the **join_jobs** commands without the need of running them all 
together again in a single **run** result, and always uses the first appearing job in the `metadata_join` as the 
reference for keeping the same msclusterIDs. 

This command uses a different metadata, called `metadata_join`, defining the jobs to be joined and their reference 
codes. It uses the clean results from the provided NP³ jobs and execute the main pipeline from Steps 3 to 10 with some 
modifications and adaptations in an incremental clustering manner, except for Step 8 which is skipped 
(it may be executed a posteriori if needed). At the end, the final reports are created for the joined results!

The next [subsection](metadata_join.md) define how to create the `metadata_join` table, its expected format, and 
the data organization to execute the **join_jobs** command. More details on its [methodology](details.md) in the 
following subsection.

The **join_jobs** command may be useful for processing growing libraries, which will have new datasets being included 
from time to time and need to keep integrating new results without loosing previous reference IDs from the main result; 
or for processing very large jobs, which may be divided into smaller jobs and then joined by chunks with a smaller 
memory footprint (divide and conquer strategy). 

## Parameters

For the complete list of parameters run:
 
```{ .text .copy }
node np3_workflow.js join_jobs --help
```

Parameters for the command **join_jobs**:
 
- *-n, --output_name* <name\> : the job name. It will be used to name the output directory and 
  					the results from joining the jobs. It must have less than 80 characters.
  
- *-m, --metadata_join* <file\> : path to the metadata table CSV file defining the jobs to be joined. Different format, see manual.
- *-d, --jobs_data_path* <path\> : path to the folder containing the input jobs result to be joined, 
  this should contain their previous NP³ result, named accordingly to what is specified in the metadata_join. 
  Their clean mgf and quantification tables will be used.
- *-y, --pre_processed_dir_path* <path\> : path to the folder containing the input jobs pre processing result, 
  this should contain all the original jobs previous NP³ pre processing result in separated folders named 
  accordingly to what is specified in the metadata_join.
- *-o, --output_path* <path\> : path to where the output directory will be created
- *-f, --fragment_tolerance* [x] : the tolerance in Daltons for fragment peaks. Peaks in the
  					original spectra that are closer than this get merged by
  					the NP3_MSCluster algorithm. Also used in the pre process (Step 2)
   (default: 0.05)
- *-z, --mz_tolerance* [x] : this is the tolerance in Daltons for the m/z of the
  					precursor that determines if two spectra will be compared
  					and possibly joined. Used in the clustering job and
  					in the library identifications (Step 6) (default: 0.025)
- *-p, --ppm_tolerance* [x] : the maximal tolerated m/z deviation in parts per million (ppm)
  					to be used in the pre-processing step if ran
  					 (default: 5)
- *-t, --rt_tolerance* [x] : tolerance in seconds for the retention time width of the
  					precursor that determines if two spectra will be compared
  					and possibly joined. It is directly applied to the retention
  					time minimum (subtracted) and maximum (added) of the spectra.
  					 (default: 2)
- *-a, --ion_mode* [x] : the precursor ion mode. One of the following numeric values
  					corresponding to a ion adduct type: '1' = [M+H]+ or 
  					'2' = [M-H]- (default: 1)
- *-i, --similarity_function* [x] : the similarity function to be used in the spectra comparison 
  					to create the pairwise similarity table after clustering and clean steps. 
  					If "spec2vec" is selected, the model trained on UniqueInchikey subset (12,797 spectra) 
  					is used by spec2vec in the spectra comparison and the matchms library is used to compute 
  					the number of matched peaks between the compared spectra; otherwise, the NP³ shifted cosine 
  					function is used. One of "np3_shifted_cosine" or "spec2vec". (default: "np3_shifted_cosine")
- *-s, --similarity* [x] : the minimum similarity to be consider in the hierarchical
  					clustering, starts in 0.70 and decrease to X in 15 rounds
  					 (default: 0.6)
- *-g, --similarity_blank* [x] : the minimum similarity to be consider in the hierarchical
  					clustering of the blank clustering steps, starts in 0.70
  					and decrease to X in 15 rounds. Only used in the 
  					clustering of blank samples (column SAMPLE\_TYPE equals
  					'blank' in the metadata table) (default: 0.3)
- *--bflag_cutoff* [x] : A positive numeric value to scale the interquartile range (IQR)
  					of the blank spectra basePeakInt distribution from the clustering result and to allow spectra with a basePeakInt
  					value below this distribution median plus IQR*bflag_cutoff to be
  					joined with a blank spectrum during the clean Step 5, without relying on the similarity value.
  					Or FALSE to disable it.
  					The IQR is the range between the 1st quartile (25th quantile) and the
  					3rd quartile (75th quantile) of the distribution. The spectra with a
  					basePeakInt value <= median + IQR*bflag_cutoff (from the
  					blank spectra basePeakInt distribution) and BFLAG TRUE will be joined to a blank spectrum
  					in the clean Step 5. This cutoff will affect the spectra with BFLAG TRUE
  					that would not get joined to a blank spectra when relying only on the
  					similarity cutoff. This is a turn around to the fact that blank spectra
  					have low quality spectra and thus can not fully rely on the similarity values. (default: 1.5)
- *--noise_cutoff* [x] : A positive numeric value defining the minimum base peak intensity absolute value 
  					that a MS2 spectra must have to be kept after the clustering of Step 3. 
  					The MS2 spectra with a basePeakInt smaller than this value will be removed before clean Step 5.
  					 The default value is zero (disabled). Large values in this parameter may result in the loss of minority 
  					compounds together with noise spectra. (default: 0)
- *-c, --scale_factor* [x] : the scaling method to be used in the fragmented peak's
  					intensities before any dot product comparison (Step 3).
  					Valid values are: 0 for the natural logarithm (ln) of the
  					intensities; 1 for no scaling; and other values greater
  					than zero for raising the fragment peaks intensities to
  					the power of x (e.g. x = 0.5 is the square root scaling).
  					[x] >= 0 (default: 0.5)
- *-e, --method* [name] : a character string indicating which correlation coefficient is
  					to be computed. One of "pearson", "kendall", or "spearman"
  					 (default: "spearman")
- *-x, --min_peaks_output* [x] : the minimum number of fragment peaks that a spectrum must have
  					to be outputted after the final clustering step. Spectra
  					with less than x fragmented peaks will be discarded. x >= 1
  					 (default: 5)
- *-j, --tremolo_identification* [x] : (not Windows OS's) A logical "TRUE" or "FALSE" indicating if
  					the Tremolo tool should be used for the spectral matching
  					against the ISDB from the UNPD (default: "TRUE")
- *-r, --trim_mz* [x] : A logical "TRUE" or "FALSE" indicating if the spectra fragmented 
  					peaks around the precursor m/z +-20 Da should be deleted 
  					before the pairwise comparisons. If "TRUE" this removes the 
  					residual precursor ion, which is frequently observed in MS/MS 
  					spectra acquired on qTOFs. (default: "TRUE")
- *--max_shift* [x] : Maximum difference between precursor m/zs that will be used in the search of shifted m/z fragment ions in the NP³ shifted cosine function. Shifts
                                       greater than this value will be ignored and not used in the cosine computation. It can be useful to deal with local modifications of the same
                                       compound. (default: 200)
- *-l, --parallel_cores* [x] : the number of cores to be used for parallel processing
  					in Step 5 spectra comparison. x = 1 for disabling parallelization and x > 2
  					for enabling it. x >= 1 (default: 2)
- *-w, --similarity_mn* [x] : the minimum similarity score that must occur between a pair of consensus MS/MS spectra in order to create an edge in the molecular networking. Lower
                                       values will increase the component size of the clusters by inducing the connection of less related MS/MS spectra; and higher values will  limit the
                                       components sizes to the opposite (default: 0.6)
- *-k, --net_top_k* [x] : the maximum number of connection for one single node in the
  					similarity molecular networking. An edge between two nodes
  					is kept only if both nodes are within each other's [x]
  					most similar nodes. Keeping this value low makes 
  					very large networks (many nodes) much easier to visualize (default: 15)
- *-x, --max_component_size* [x] : the maximum number of nodes that all component of 
  					the similarity molecular network must have. The edges of 
  					this network will be removed using an increasing cosine 
  					threshold until each network component has at most X nodes. 
  					Keeping this value low makes very large networks (many nodes 
  					and edges) much easier to visualize. (default: 200)
- *--min_matched_peaks* [x] : The minimum number of common peaks that two spectra must share to be connected by an edge in the filtered SSMN. Connections between spectra with less
                                       common peaks than this cutoff will be removed when filtering the SSMN. Except for when one of the spectra have a number of fragment peaks smaller than
                                       the given min_matched_peaks value, in this case the spectra must share at least 2 peaks. The fragment peaks count is performed after the spectra are
                                       normalized and cleaned. (default: 6)
- *--blank_expansion* [x] : the distance of neighborhood nodes from the blanks in IVAMN to be 
  					selected for removal in the final protonated networks. 
  					(0) to only remove blanks nodes,  
  					(1) to remove nodes directly connected to a blank node, 
  					(2 or greater) to remove nodes in a distance equal to 2 or greater from a blank node, 
  					or (-1) to remove all possible neighbors and ancestors of a blank node (remove blank clusters) from IVAMN  (default: 0)
- *-b, --max_chunk_spectra* [x] : Maximum number of spectra to be loaded and processed in a
  					chunk at the same time. In case of memory issues this
  					value should be decreased (default: 3000)
- *-v, --verbose* [x] : for values X>0 show the scripts output information
  					 (default: 0)
- *-h, --help* : display help for command


## Results

A directory inside the `output_path` named with the `output_name` containing:

- A copy of the `metadata_join` file and the command line parameters values used in a file named 'logRunParms', for 
reproducibility.
- Two automatically created metadata tables containing: 
    - One the original samples concatenated in a single file named 'original_samples_METADATA.csv'; 
    - And another with the original NP³ jobs that were joined in this process and in any provided joined job in a 
  file named 'original_jobs_METADATA_JOIN.csv'. For reproducibility and future joins with this result.
- A folder named 'outs' with the clustering result of the joined jobs in a single sub folder named with the 
`output_name` containing:
  - A sub folder named 'count_tables' with the Step 4 quantification in CSV tables named as 
  '`output_name_<spectra|peak_area>.csv'. And inside it the clean tables in a folder named 'clean'.
  - Another sub folder named 'clust' with the clusters membership files (which SCANS or msclusterID were joined).
  - A third sub folder named 'mgf' with the resulting clean consensus spectra in MGF files.
  - A fourth sub folder named 'identifications' with the tremolo identification results in a csv table.
  - A fifth sub folder named 'molecular_networking' with the molecular networking of this joined job, both SSMN and IVAMN.
  - A text file named 'logClusteringOutput' with the NP3_MSCluster log output.

## Examples

Join original jobs A and B (results from the run command).

```{ .text .copy }
node np3_workflow.js join_jobs --output_name "test_join_a_b" --output_path 
"/path/where/the/output/will/be/stored" --metadata_join 
"/path/to/the/metadata_join/file/test_np3_join_a_b_metadata.csv" 
--pre_processed_dir_path "/path/to/the/dir/with/joining/jobs/pre/process/results" 
--jobs_data_path "/path/where/the/joining/jobs/output/is/stored"
```

Join the joined job from A and B (result from the join_jobs command) with original job C (result from the run command).
```{ .text .copy }
node np3_workflow.js join_jobs -n "test_join_ab_c" -o 
"/path/where/the/output/will/be/stored" -m 
"/path/to/the/metadata_join/file/test_np3_join_ab_c_metadata.csv" -y 
"/path/to/the/dir/with/joining/jobs/pre/process/results" -d 
"/path/where/the/joining/jobs/output/is/stored" -t 5 -v 10
```

## *join_jobs* Workflow

This command performs:

1. The setup of the `metadata_join` information of the defined NP³ jobs; 
    - Concatenates the samples metadata of all original jobs (from the original **run** results); 
2. Fix of the original samples SAMPLE_CODE (resolve duplicates);
3. Executes the Clustering Step 3 and 4 of all the clean data from the provided jobs together; 
    - Keep as reference the msclusterID from the first appearing job in the provided `metadata_join`, additional m/z will 
   have an incremental ID starting from the biggest reference ID
4. Executes the Clean step 5 and compares the final joined clean consensus spectra pairwise;
5. Executes the Spectral Identification Step 6 (UNPD) and Step 6.1 (GNPS2) to retrieve library annotations;
6. Executes Step 7 adapted for joining jobs: updates the original jobs ionization variants annotations with the final joined 
   clean consensus spectra and merge the updated IVAMNS, recomputes the [M+H]+ and keeps the previous protonated 
   representatives that remain valid;
7. Executes Step 9 Biocorrelation and groupings using the original samples bioactivity and grouping, if present;
8. Executes Step 10: creates the spectra similarity molecular networking (SSMN) with the final joined result and the 
corresponding protonated networks.

The **join_jobs** can be used to join the results from multiple original jobs and 
also from previous joined jobs with a new original or joined job. There is no limit for
joining NP³ results, coming from the **run** or the **join_job** commands, the user may continuously join results from 
different jobs in a final incremental result without loosing reference from a previous main result. 

If there is any duplicated sample code among different original jobs being joined (SAMPLE_CODE column from the original 
samples' metadata), these codes are automatically updated with a numerical and incremental suffix. The original sample 
codes are maintained in a separated column (see details in next [subsections](details.md)).