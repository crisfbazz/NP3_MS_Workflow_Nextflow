## Step 1: Creation of a Metadata Table
 
The metadata table must be a **CSV** file with columns separated by comma ',' and with the dot '.' as decimal point character. 
It can be created and edited in EXCEL, LibreOffice CALC or any other spreadsheet program, but must be saved in CSV (UTF-8) format. 
If any text name in the metadata file contains special characters, the file must be saved with the parameter to quote text cells set to TRUE. 

After the metadata file was finished it is recommended to open it with a simple text viewer, e.g., notepad, to make sure that the column separator character and the decimal point character were properly set. 
 
A metadata template file is available in the NP³ MS workflow repository: '[METADATA_TEMPLATE.csv](https://github.com/danielatrivella/NP3_MS_Workflow/blob/master/METADATA_template.csv)'.
 
## Format 

The metadata table must contain the following **mandatory columns** (in uppercase):
 
- _FILENAME_ \- the raw data file names (including extensions but not including complete paths). 
        - The accepted file format is mzXML or mzML (mzData format could also be used, but was not vastly tested).
- SAMPLE_CODE \- unique and syntactically valid (see below) names for each file. 
The use of the same SAMPLE_CODE in more than one entry will generate an error.
    - A syntactically valid name must start with a letter and consist of letters, numbers, and underscore characters. 
  [R reserved words](https://cran.r-project.org/doc/manuals/r-release/R-lang.html#Reserved-words) are not syntactically valid names.
    - The SAMPLE_CODE names will be used to name and reference the pre-processed files in Step 2 and to name the output 
count table's columns of Steps 3 to 9.
- _DATA_COLLECTION_BATCH_ \- a number starting with 1 and progressing in ascending order (e.g., 1,2,3,...) indicating 
to which data collection batch group each sample file belongs. 
    - Samples from the same data collection batch should be chemically and physically related 
(e.g. collected in the same LC-MS/MS run, from the same screening, from the same origin or species), 
ensuring closely related experimental conditions. 
    - See below about the batches impact and the order of the samples in each batch.
- _SAMPLE_TYPE_ \- one among "sample", "blank", "bed", "control" or "hit". 
    - If the sample file is a chromatographic blank this column must be set to "blank" (e.g. DMSO, solvents, etc); 
    - If the sample is the culture media (containing molecules that should be identified for later exclusion or separation) set it to "bed"; 
    - If the sample is a control (e.g. the culture media + the non-elicited microorganism) set it to "control"; 
    - If the sample is bioactive and is desired to dereplicate its m/z's, set it to "hit"; 
    - And if the file is not a blank sample nor any of the above set it to "sample" - ordinary sample. 
    - Spectra from blank, control, bed and hit samples will be used to signal possible false positive spectra and to 
dereplicate m/z's that appear in hit samples (see Step 4 details).
 
The metadata table may contain the following **optional columns**:
 
- *BIOACTIVITY_<NAME\>* - zero or more columns starting with the prefix "BIOACTIVITY_" followed by a unique name <NAME\> defining different bioactivities values for each sample. 
    - Set it to a positive number (greater or equal than zero) indicating the sample bioactivity score, it could be an inhibition or an activation of a target. 
    - These values should be normalized, e.g., from 0 to 100, to be used in the correlation computation.
- *COR_<NAME\>* - zero or more columns starting with the prefix "COR_" followed by a unique name <NAME\> defining the correlation group.
Each one of these columns must have values 0 or 1. 
    - If the sample spectra or peak area count is to be used in the COR_<NAME\> correlation with each defined bioactivity scores (BIOACTIVITY_<NAME\> columns), set it to 1. Otherwise, set it to 0. 
    - In Step 9 each one of these columns will result in a new column for each available bioactivity score in the output count tables, 
containing the bioactivity correlation scores computed with the respectively selected samples. 
    - For better correlation results the selected groups of samples should belong to the same bioactivity peak, e.g., a selection that contains very active, middle active and not active samples of a same bioactivity peak (a 'hit'). 
    - We strongly recommend to use the blank, bed and control samples in this selection to eliminate false positives. 
    - Bioactivity data and LC-MS/MS data must be collected from chemical samples at the same concentration. 
  This will guarantee real correlations, especially if the samples and bioactivity data come from a screening of natural product samples.
- *GR_<NAME\>* - zero or more columns starting with the prefix "GR_" followed by a unique name <NAME\> defining a set of grouping for the spectra quantification. 
    - Each of those columns may have one or more tags to label a group of samples. 
    - A column "GR_<NAME\>" with groups "<g1\>", "<g2\>" and "<g3\>" will bind to the final count tables the columns "<NAME_g1\>", 
  "<NAME_g2\>" and "<NAME_g3\>" containing the sum of peak areas or number of spectra of the samples that belongs to each group for each m/z row. 
    - For example, to create a set called *colors* with the groups *red*, *green* and *blue*, add to the metadata the 
  column "GR_colors", tag the samples that belongs to the groups *red*, *green* and *blue* by filling the rows 
  (you may have empty observations if it is necessary) and it will bind the columns *colors_red*, *colors_green* and *colors_blue* to the final count tables.
 
The user is free to add any additional column to the metadata file, for example to add relevant descriptions for each file. 
These additional columns will be ignored by the workflow commands, but can be useful to add information related to the job samples. 

It's important that additional column names are not equal or a prefix of any of the mandatory or optional columns of the metadata file.
 
## Data Collection Batch Impacts

#### Misalignment suggestion

The first files (rows) of each data collection batch number of a metadata table, excluding blanks, 
are used at the end of **Step 2** - pre_process to perform the suggestion of misalignment between the samples 
(this suggestion may help setting the parameter rt_tolerance value in Steps 3 to 7). 

This information must be taken into account when defining the files order in the metadata table. 
It is expected that the first file of each data collection batch number is the more related sample between the 
different data collection batches, e.g., 
those that have the same polarity or other characteristic that provides them bigger chances of having more MS1 peaks in common. 
Files from extracts are a good choice to go first. This will only influence in the proposed suggestion for 
the retention time tolerance value, it will **not modify the data** (see Step 2 pre_process Details).
 
#### Clustering 

The **Step 3** clustering will be performed in batches using the DATA_COLLECTION_BATCH number and the 
SAMPLE_TYPE information to group the samples in the following criteria:

1. First, all the samples, excluding blanks, (column "SAMPLE_TYPE" equal to "sample", "hit", "bed" or "control") 
from the same batch (which have the same DATA_COLLECTION_BATCH number) are clustered in the *data clustering step*. 
2. Next, all the blank samples (column "SAMPLE_TYPE" equal to "blank") from the same batch are clustered in the 
*blank clustering step*, in which the retention time is not used to better deal with baseline blanks. 
3. Subsequently, all sub-batches (data and blank steps results) from the same data collection batch are clustered 
together in the *data collection batch integration step*.
4. Finally, all batches are clustered together in the *final integration step*. 

This way, spectra that were detected in the same conditions, and thus tend to be more similar, 
are enriched first before being clustered with the spectra from a different batch (less related samples). 
And at the end, all spectra are clustered together.

#### Ionization Variant Annotation and [M+H]+ choice

During **Step 7** Annotation [M+H]+ only the consensus spectra *m/z* that appear in at least one data collection batch 
in common may be annotated as a pair of ionization variants. 

This means that the ionization variants are only searched within each data collection batch, 
and again the data collection batch should group the related samples (e.g. same extract, same origin, same species) 
from which it's expected and possible (physically or chemically) to find ionization variants.

Thus, the DATA_COLLECTION_BATCH grouping will directly impact the number of possible ionization variants annotations and 
this will indirectly impact the choice of the putative [M+H]+, which are selected based on the set of ionization variants annotations.
