# Step 8: Command *merge*
 
Runs the merge of the clean count tables based on the annotated ionization variants and putative [M+H]+ 
(for positive ion mode only). 

It creates new symbolic spectra candidates representing the union of each [M+H]+ representative 
consensus spectra with its annotated ionization variants. 
This union is performed for each type of annotation (adducts + neutral losses, multiple charges, dimers/trimers, 
isotopes and in-source fragments) and by combining all of them together, what can lead to at most 31 new symbolic 
spectra by consensus spectra. 
 
By default, the merge is only performed for the consensus spectra assigned as a [M+H]+ representative, 
to better account for the quantifications of the putative true metabolites. 
To merge using all consensus spectra set the parameter `merge_protonated` to FALSE.
 
## Parameters

For the complete list of parameters run:
 
```{ .text .copy }
node np3_workflow.js merge --help
```

Parameters for the command **merge**:

- *\-o, \-\-output_path* \<path\>       : path to the output data folder, inside the 'outs' directory of the clustering result folder. It should contain the 'counts_table' folder and inside it the 'clean' subfolder with the clean count tables in CSV files. The job name will be extracted from here
- *\-y, \-\-processed_data_dir* \<path\>  : the path to the folder inside the raw data folder where the pre-processed data (MGFs) were stored.
- *\-m, \-\-metadata* \<file\>          : path to the metadata table CSV file
- *\-p \-\-merge_protonated* [x]          : A boolean TRUE or FALSE indicating if only the [M+H]+ representative consensus spectra should be merged. If FALSE merge all msclusterID's (default: "TRUE")
- *\-e, \-\-method* [name]           : a character string indicating which correlation coefficient is to be computed. One of “pearson”, “kendall”, or “spearman” (default: spearman)
- *\-v, \-\-verbose* [x]             : for values X\>0 show the scripts output information. (default: 0)
- *\-h, \-\-help*                    : output usage information
 
## Results
 
One subfolder inside the 'count_tables' folder is created named 'merge' containing:
 
- Two CSV files with the cleaned and annotated counts of spectra and peak area merged and new symbolic clusters added 
as new rows, named with the suffix '_merged_ann.csv';
- CSV files with the correlation columns added are also included when there is a biocorrelation result.
 
## Examples
 
Fake example to execute the Merging command:

```{ .text .copy }
node np3_workflow.js merge --output_path "/path/to/the/output/dir/test_np3/outs/test_np3"
```
