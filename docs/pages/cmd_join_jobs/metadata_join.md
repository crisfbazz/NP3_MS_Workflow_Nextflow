# **join_jobs** Setup

## Metadata_join Format

The `metadata_join` table defines the list of NP³ jobs to be united, it defines the joining jobs codes, original 
samples' metadata name and original pre processing result. It also defines if the joining jobs came from a **run** 
(called here original job) or a **join_jobs** (called here joined job) command result. The order of the jobs defined in 
the metadata_join matters, the *first* defined job is used as the *reference* job, from which the msclusterIDs' will be 
maintained in the final clean result (incremental clustering).

This table must be a CSV file with columns separated by comma 
',' and with the dot '.' as decimal point character. It can be created and edited in EXCEL, LibreOffice CALC or any 
other spreadsheet program, but must be saved in CSV (UTF-8) format. If any text name in the metadata_join file contains 
special characters, the file must be saved with the parameter to quote text cells set to TRUE. After the *metadata_join* 
file was finished it is recommended to open it with a simple text viewer, e.g., notepad, to make sure that the column 
separator character and the decimal point character were properly set. 
 
A *metadata_join* template file is available in the NP³ MS workflow repository, named 
'[METADATA_JOIN_TEMPLATE.csv](https://github.com/danielatrivella/NP3_MS_Workflow/blob/master/METADATA_JOIN_template.csv)'. 
It uses the Bra346 and the MA9 datasets as examples, provided in the NP³ paper.
 
The `metadata_join` table must contain the following mandatory columns (in uppercase):
 
- JOBNAME \- the name of the job to be joined. This name must be equal to the *output_name* used to obtain this job 
result, it must match the job final results naming. This name will be used to retrieve the job data to perform the joining.
- JOB_CODE \- unique and syntactically valid (see below) name for each job. The use of the same JOB_CODE in more than 
one entry will generate an error. A syntactically valid name must start with a letter and consist of letters, numbers, 
and underscore characters. R reserved words (https://cran.r-project.org/doc/manuals/r-release/R-lang.html#Reserved-words) 
are not syntactically valid names. The JOB_CODE name will be used to reference the original NP³ jobs and to keep the 
joined msclusterIDs track.
- METADATA_NAME \- the name of the metadata table used to execute the job, if this job is a result of the *run* command 
(JOINED_JOB column equals 0) this is the original samples metadata name, or if this job is a result of the **join_jobs** 
command (JOINED_JOB column equals 1), this is the used metadata_join name. This name must be the complete file name with 
extension, and it must point to the metadata used to obtain the provided results, automatically stored in the root 
folder of the jobs results. It will be used to automatically retrieve the data and information from the jobs results.  
- PRE_PROCESSED_DATA_NAME \- the name of the folder with the pre processing result of the original job. It will be used 
to retrieve the pre processed data of the job. If this is a joined job (JOINED_JOB column equals 1), this column is skipped.
- JOINED_JOB \- A 0 or 1 column value to define if this job is the result of a original job coming from the *run* 
command execution (0) or if it is the result of a joined job coming from another **join_jobs** command execution (1).

The user is free to add any additional column to the metadata_join file, for example to add relevant descriptions for 
each job. These additional columns will be ignored by the workflow commands, but can be useful to add information 
related to the jobs being joined. It is important that additional column names are not equal or a prefix of any of the 
mandatory or optional columns of the metadata_join file.

The samples' quantification groupings, correlations and bioactivities will be retrieved from the samples' metadata of 
the original jobs. These information are not informed again, the user may modify the samples metadata of the jobs being 
joined to automatically use it in the **join_jobs**. Anyway, in the **join_jobs** setup, the original samples metadata is 
automatically created (more details in next [subsection](details.md)) and may be used a posteriori to add new biocorrelations or quantification groupings to the 
joining results, and the command **corr** may be used to recompute them.

## Data Organization

The data organization to execute the **join_jobs** command must contain three things (defined by mandatory parameters):

1. A metadata_join table defining the jobs to be united and following the format and mandatory columns defined above
(parameter `metadata_join`).
2. A folder containing the results of the NP³ jobs to be united, these may be the results from a **run** or a previous 
**join_jobs** execution, in all cases only the results of the jobs being joined are needed. All jobs defined in the 
metadata_join must be placed inside this folder. The data to be united will be retrieved from here (parameter
`jobs_data_path`).
3. And a folder containing the results of the pre-processing of the original NP³ jobs to be united or united by 
any included joined job. These data will be used to correctly compute the peak area of the final joined clean consensus 
spectra (paramater `pre_processed_dir_path`).

In a project intended to continually grow with new jobs being collected and processed, this structure must be kept to 
continuous join a new processed dataset.
