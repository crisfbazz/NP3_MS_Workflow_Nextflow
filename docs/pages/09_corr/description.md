# Step 9: Command *corr*
 
The biocorrelation command computes for each consensus spectra present in the provided count table the 
**biocorrelation score** between its count of spectra or of peak area by sample (quantification) and the samples 
bioactivity score. This command also computes the **grouping** of 
the consensus spectra quantification by sample if defined in the `metadata` table (columns named with prefix "GR_"). 

The biocorrelation is performed for each correlation group pairwise combined with each bioactivity score 
defined in the `metadata` table (columns named with the prefix "BIOACTIVITY_" and "COR_"). 
The correlation may be performed using different `method`s: 'pearson', 'kendall' or 'spearman' (default).
 
The quantification grouping may be used to aggregate any characteristics of the samples. And the
biocorrelation score can be used to rank the consensus spectra and to select the top possible candidates 
responsible for the observed hits in bioactivity experiments. 
These information may also be visualized in the molecular networks and applied for the discovery of active compounds, 
as done in[^1][^2].
 
## Parameters

For the complete list of parameters run:
 
```{ .text .copy }
node np3_workflow.js corr --help
```

Parameters for the command **corr**:
 
- *\-b, \-\-metadata* <file\>     : path to the metadata table CSV file. Used to retrieve the biocorrelation groups and quantification grouping. For a joined job, this must be the original samples metadata table.
- *\-c, \-\-count_file_path* <file\>  : path to the count table CSV file
- *\-e, \-\-method* [name]  : a character string indicating which correlation coefficient is to be computed. One of “pearson”, “kendall”, or “spearman” (default: spearman)
- *\-b, \-\-bio_cutoff* [name]  : a bioactivity cutoff value. Bioactivities in the metadata table that are less than the bio_cutoff value will be set to zero. (default: 0)
- *\-v, \-\-verbose* [x] : for values X>0 show the scripts output information (default: 0)
- *\-h, \-\-help* : output usage information
 
## Results
 
A CSV table in the count_file directory named as:

- '<count\_file\_name\>\_corr\_<method\>.csv' 

Where the *count\_file\_name* is the name of the provided `count_file_path` 
  table. The new tables contain the original count tables data plus new columns with the correlation scores of each 
  consensus spectrum for each correlation group and bioactivity score, as defined in the `metadata` table. 
  It may also contain new columns for the quantification grouping, if defined in the `metadata`.

Another table is also created as a copy of the first table, named with the suffix 'bioAct' as:

- '<count\_file\_name\>\_corr\_`method`\_bioAct.csv'

Containing new rows at the beginning of the table with the bioactivities scores of each sample placed above the 
  original counts table header (first rows) - as provided in the `metadata` table. 
This may ease the visualization of the quantification alongside with the bioactivity scores.
 
## Examples
 
Fake example to execute the Biocorrelation command:

```{ .text .copy } 
node np3_workflow.js corr --metadata "/path/to/the/metadata/file/test_np3_metadata.csv" 
--count_file "/path/to/the/metadata/file/np3_job_spectra.csv"
```

## References
[^1]: C.F. Bazzano, et al. (2024). *NP³ MS Workflow*. Analytical Chemistry 96 DOI: 10.1021/acs.analchem.3c05829
[^2]: Louis-Félix Nothias, Mélissa Nothias-Esposito, Ricardo da Silva, Mingxun Wang, Ivan Protsyuk, Zheng Zhang, Abi Sarvepalli, Pieter Leyssen, David Touboul, Jean Costa, Julien Paolini, Theodore Alexandrov, Marc Litaudon, and Pieter C. Dorrestein. Bioactivity-Based Molecular Networking for the Discovery of Drug Leads in Natural Product Bioassay-Guided Fractionation Journal of Natural Products 2018 81 (4), 758-767 DOI: 10.1021/acs.jnatprod.7b00737