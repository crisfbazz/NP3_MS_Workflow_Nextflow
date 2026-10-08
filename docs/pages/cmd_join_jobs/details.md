# Command *join_jobs* Details

The join_jobs command uses the main NP³ MS Workflow steps (from command **run**) to integrate different jobs by uniting 
their NP³ results. Some steps needed to be adapted (Steps 4 and 5) or replaced (step 7) to allow joining original NP³ 
results using an incremental approach.

The inputs to the join_jobs command are the pre processed data and the final results from other NP³ executions (the clean 
data is used), coming from a **run** execution (called original job) or from a **join_jobs** execution (called joined job). 

The clean data used are the clean MGF file, the clean count tables and the annotations of 
ionization variants. The ionization variants are not recomputed here because they are directly related to the original 
samples which were collected together, and instead, the original annotations are maintained from the original jobs, 
joined and updated here. 

## **join_jobs** Setup

The **join_jobs** command starts with a setup process to retrieve the metadata information of the NP³ jobs being joined 
and to automatically create two default metadata:

1. The first is the **original samples metadata**, named "original_samples_METADATA.csv", which concatenates the original 
samples of all original jobs being joined by concatenating all original metadata. The format of the samples metadata 
follow the format of the original metadata presented in section 4.1 which is Step 1 of the NP³ MS workflow. 
It uses the column METADATA_NAME from the `metadata_join` to retrieve all the samples metadata from the original jobs.
If there is a joined job in the current join_jobs execution (JOINED_JOB == 1 in the provided `metadata_join`), 
its original samples metadata is retrieved and concatenated with the list of samples from the new original jobs or from 
other joined job present in the provided `metadata_join`.
    - Two additional columns are created here to allow duplicated SAMPLE_CODE among different jobs: 
   (1) the "SAMPLE_CODE_ORIGINAL": which stores the original sample codes from the original jobs, as it is found in 
   the *run* results of the original jobs; and (2) the "SAMPLE_CODE_JOINED_JOBS": which stores the SAMPLE_CODE used in 
   the last result of the provided jobs (useful when joining jobs recursively/continuously).
    - The SAMPLE_CODE column is thus modified to automatically remove duplicates. To do this, the *join_jobs* setup adds 
   a numerical and incremental suffix to all duplicated columns equals to a '_' followed by a number starting in 1. 
   This number is incremented until all duplicates are removed.
2. And the other is the **original joined jobs metadata**, named "original_jobs_METADATA_JOIN.csv", which concatenates the 
list of all original NP³ jobs being joined. Its format follows the `metadata_join` format. If there is a joined job in 
the current join_jobs execution (JOINED_JOB == 1 in the provided `metadata_join`), its original joined jobs metadata is 
retrieved and concatenated with the list of new original jobs or other joined jobs metadata present in the provided 
`metadata_join`. 

With this two metadata, the reference to the original jobs is automatically kept throughout different executions of 
the *join_jobs* commands, even when the executions recursively call the join of previous joined jobs. 
With the original samples metadata created, a check is performed to guarantee that all jobs' code are unique among the 
jobs being joined, otherwise the join will fail and the user must correct this to proceed.

## **join_jobs** Workflow

The **join_jobs** processing workflow starts with the clustering (Step 3) of the clean MGF of all jobs being joined. The 
clustering step reduces all the clean MGF of the provided jobs to a single clustered MGF and the command proceed to the 
quantification by original samples (Step 4) of the clustered consensus spectra that were joined. 

The quantification step was adapted here to use the clean tables from the jobs being joined and to merge their counts 
for the new clustered spectra. It starts by updating the SAMPLE_CODE from the quantification columns using the last 
used value present in SAMPLE_CODE_JOINED_JOBS column and replacing it by the list with no duplicates present in the 
SAMPLE_CODE column (if there is any). It also updates the SAMPLE_CODE in the scans and peakIds columns, 
guaranteeing that the correct reference to the original samples is not lost due to any resolved duplicates. It creates 
two columns: 

- One called "joinedOriginJobsID" to keep the reference to the original msclusterIDs of the original jobs 
that may be clustered here, it set this column equals to the original msclusterID concatenated with the original 
JOB_CODE for all clean consensus spectra from a original job and then concatenate these values for all clustered spectra; 
- And another called "joinedJobsID" to keep the reference to the last msclusterIDs from the joining jobs that were 
clustered, similarly this column is set equals to the msclusterID of the joining results concatenated with the provided 
and respective JOB_CODE. For an already joined job being joined again, the column "joinedOriginJobsID" is kept for the 
reference. 

In the resulting quantification table, a new column named "msclusterID_integrative" is created to store the
msclusterID of the reference job, which is the first appearing job in the metadata_join. This column will keep the 
reference to a previous result, and will allow keeping the same msclusterID for the same spectra in the joined jobs result. 
Here the "msclusterID_integrative" will concatenate all msclusterID coming from the reference job in each joined 
clustered spectra and will give a new id to the new spectra that did not appear in the reference job equal to the 
maximum ID present in the reference job plus one and incrementally until all joined clustered spectra are named.

The quantification step performs the count by number of spectra as a simple sum of the clean consensus spectra that 
were joined in the clustering step and performs the count by peak area using the pre-processed data of all original 
samples being joined. It uses the "scans" column to retrieve the reference to the original peaks information and the 
original samples metadata to retrieve all the original SAMPLES_CODEs present in the pre-processed result, and then, 
correctly compute the peak area of the joined clustered consensus spectra. The peak area computation was adapted to 
allow retrieving the peaks information from multiple pre-processed data, the PRE_PROCESSED_DATA_NAME from the 
metadata_join is used here for reference.

Next, the join_jobs workflow proceeds to the cleaning (Step 5), which is executed similar to the main workflow. The 
main difference here is that at the end of the process the final msclusterID are saved in a new column called 
"cleanClustID" and the msclusterIDs are set to be equal to the minimum ID present in the mscluster_integrative column - 
this way the IDs from the reference job are kept in the final joined clean result. Another difference is that the peak 
area computation was adapted for multiple pre-processed data, but works similarly. Then, the clean joined data are 
compared pairwise and the joined clean MGF is identified against UNPD using tremolo (Step 6) and against GNPS2 (Step 6.1).

Following is the annotate_protonated (Step 7) which works in an incremental approach instead of recomputing the ionization 
variants. The annotate_protonated step was replaced here by a procedure to merge the ionization annotations of the jobs 
being joined and then to recompute the [M+H]+ while keeping some of the original protonated representatives. 

The *join_jobs* flow of Step 7 is as follows: 
1. first it retrieves the IVAMNs of all the joining jobs and map their 
msclusterIDs to the new joined IDs (uses column joinedJobsID for reference), selfloops are removed and duplicated 
edges are merged. 
   - The edges attributes are recomputed with the new similarity values, the mzError receives the mean 
value of the merged edges and the rtError is recomputed again for the joined msclusterIDs. 
   - The annotations rules are applied again and removes any invalid annotation (similarity cutoff), 
   any resulting empty edge is removed. 
   - This result in the final joined IVAMN, selfloops are added again at the end for the missing nodes. 
2. Then, the find protonated script is executed for the final joined IVAMN, resulting in a new list of [M+H]+ ions 
(column protonated_representatives = 1). 
   - This list is further merged with the previous protonated_representatives that have an in-degree > 0 in the final joined 
   IVAMN (keep [M+H]+ with valid annotations), this prevents loosing relevant information from the previous jobs. 
3. Finally, a new script parses the information present in the final joined IVAMN and writes them as annotation columns 
to the joined clean tables, with a similar format to the main flow. 
   - The only missing column here is the 'analogs' one, which is not 
      recomputed in the new procedure. The rest of the ionization annotation columns and the protonated_representative are 
      written to the clean tables and to the joined IVAMN attribute table.

Next, the _join_jobs_ workflow skips Step 8, the merging may result in too much additional data and may be executed 
separated by the user at the end of this procedure if needed. To execute the merge with the _join_jobs_ result the user 
must inform the created original samples metadata as the metadata parameter. 

Then, the _join_jobs_ command proceeds to 
the biocorrelation (Step 9), using the automatically concatenated original samples metadata to perform the 
correlations and any grouping present in the original metadata. 
At the end of this command, the user may add new correlations, 
bioactivities or quantification groupings to the original samples metadata and execute this step again (using command 
**corr**).

The *join_jobs* workflow ends with the molecular networking (Step 10) similar to the main flow. The data organization of 
the output is similar to the main flow and the additional steps may be executed the same way.