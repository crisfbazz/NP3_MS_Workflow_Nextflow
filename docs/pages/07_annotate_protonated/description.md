# Step 7: Command *annotate_protonated*
 
Annotate possible ionization variants in the clean count tables and creates the **Ionization Variant Annotation Molecular 
Network (IVAMN)** - for positive ion mode only. It searches for adducts, neutral losses, multiple charge, dimers/trimers, 
isotopes and in-source fragments based on numerical equivalences and chemical rules. Finally, it runs a link analysis 
in the IVAMN to assign some consensus spectra as putative [M+H]+ representatives. 
 
The number of putative [M+H]+ representatives offer a suggestion to the number of real metabolites present in the samples.
 
The nodes of the IVAMN are the set of clean consensus spectra and the links of this network connects two consensus 
spectra that have a chemical annotation. The links have a direction, pointing from the consensus spectra considered as 
an ion variant to the consensus spectra considered as a putative [M+H]+ candidate in the respective annotation (e.g., 
[M+Na]+ -> [M+H]+).
 
## Parameters

For the complete list of parameters run:
 
```{ .text .copy }
node np3_workflow.js annotate_protonated --help
```

Parameters for the command *annotate_protonated*:
 
- *\-m, \-\-metadata* <file\>          : path to the metadata table CSV file
- *\-o, \-\-output_path* <path\>       : path to the output data folder, inside the 'outs' directory of the clustering result folder. It should contain the 'counts_table' folder and inside it the 'clean' subfolder with the clean count tables in CSV files. The job name will be extracted from here
- *\-z, \-\-mz_tolerance* [x]              : the tolerance in Daltons for matching the numerical rules of detecting adducts, neutral losses, multiple charge, dimers/trimers and isotopes variants (default: 0.025)
- *\-f, \-\-fragment_tolerance* [x]  : the tolerance in Daltons for matching the numerical rules of in-source fragments and multiple charge isotopic patterns. (default: 0.025)
- *\-t, \-\-rt_tolerance* [x]        : tolerance in seconds to enlarge the peak boundaries for detecting concurrent variant spectra. (default: "2")
- *\-i, \-\-absolute_ms2_int_cutoff* [x] :    The absolute intensity cutoff for fragmented MS2 peaks in the interval of 0 to 1000 (default: 15)
- *\-a, \-\-ion_mode* [x]             :  the precursor ion mode. One of the following numeric values corresponding to an ion adduct type: '1' = Positive [M+H]+ or '2' = Negative [M-H]- (default: 1)
- *\-u, \-\-rules* [x]              :   path to the CSV file with the accepted ionization modification rules for detecting adducts, multiple charge and dimers/trimers variants, and their combination with neutral losses. (default: "rules/np3_modifications.csv")
- *\-c, \-\-scale_factor* [x]          :  the scaling method to be used in the fragmented peak's intensities before any dot product comparison (Step 5). Valid values are: 0 for the natural logarithm (ln) of the intensities; 1 for no scaling; and other values greater than zero for raising the fragment peaks intensities to the power of x (e.g. x = 0.5 is the square root scaling). [x] >= 0 (default: 0.5)
- *\-b, \-\-max_chunk_spectra* [x]      :   Maximum number of spectra (rows) to be loaded and processed in a chunk at the same time. In case of memory issues this value should be decreased (default: 3000)
- *\-v, \-\-verbose* [x]             : for values X\>0 show the scripts output information. (default: 0)
- *\-h, \-\-help*                    : output usage information
 
## Results
 
It creates inside the 'count_tables/clean' folder:
 
- Two CSV files with the clean counts of spectra and peak area concatenated with the annotations result as new columns, 
named with the suffix '_ann.csv'
 
The 'molecular_networking' folder is created if not present yet, and inside it is created: 
 
- The IVAMN edge file named as '`output_name`_ivamn.selfloops';
- One CSV table containing the IVAMN attributes and the assigned [M+H]+ representatives, named as 
    '`output_name`_ivamn_attributes.csv'
 
Where the `output_name` is extracted from the `output_path`;
 
## Examples
 
Fake example to execute the Annotation [M+H]+ command:

```{ .text .copy }
node np3_workflow.js annotate_protonated --metadata 
"/path/to/the/metadata/file/test_np3_metadata.csv" 
--output_path "/path/to/the/output/dir/test_np3/outs/test_np3" -i 5
```
