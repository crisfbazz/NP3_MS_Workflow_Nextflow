## Step 6.1: Command **gnps_library_search**

This command runs the GNPS2 Library Search[^1] workflow (offline) to identify the informed spectra against the 
ALL\_GNPS\_NO\_PROPAGATED library (default for LC data). Library annotations data was enriched in April 2026.

The GNPS identification curation is executed at the end to score and classify the retrieved library annotations 
(more details in next subsections below).

## Parameters

Parameters for the command *gnps_library_search*:

- *\-g, \-\-input\_mgf\_file* <path\> : path to the input MGF file with the MS/MS spectra data
                                     to be searched and identified
- *\-o, \-\-output\_path* <path> : if the input is a NP3 result, the path to the final NP3
                                     output data folder, inside the outs directory of the
                                     clustering result folder. It should contain the
                                     identifications folder, if not it will be created and
                                     the results will be stored in it. If the input is not a
                                     NP3 result, this should be a chosen result folder. The
                                     job name (output_name) will be extracted from here
                                     (basename).
- *\-z, \-\-mz\_tolerance* [x]  :       the tolerance for parent mass search in Daltons,
                                     depending on the instrument and data accuracy. (default:
                                     0.025)
- *\-f, \-\-fragment\_tolerance* [x]  :     the tolerance in Daltons for fragment peaks. Used for
                                     comparing the mass of the peaks during the search.
                                     (default: 0.05)
- *\-a, \-\-ion_mode* [x]             :  the precursor ion mode. One of the following numeric values corresponding to an ion adduct type: '1' = Positive [M+H]+ or '2' = Negative [M-H]-.  This will be used to select the library MGF file. (default: 1)
- *\-i, \-\-search\_tool* [x]  :        the GNPS2 search tool to be used in the library
                                     searching against the ALL_GNPS_NO_PROPAGATED. One of
                                     "gnps_indexed", "gnps", "gnps_new" or "" (disabled). The
                                     similarity function is hardcoded to be the cosine and
                                     the peak transformation function is the square root.
                                     (default: "gnps_indexed")
- *\-s, \-\-search\_min\_cosine* [x]  :  the similarity threshold that determines if two spectra
                                     are a match. The minimum cosine for a search match.
                                     Values greater or equal than 0.7 will retrieve more
                                     accurate results. (default: 0.7)
- *-p \-\-search\_min\_matched\_peaks* [x] : The minimum number of common peaks that the searched
                                     spectra must share with a library spectrum to be
                                     considered as a match. (default: 6)
- *\-k, \-\-top\_k* [x]  :  defines the maximal number of results returned by the
                                     GNPS Search Tool for each input spectrum. Additionally,
                                     the top 1 result is always computed at the end of the
                                     search. (default: 5)
- *\-r, \-\-trim\_mz* [x]  :      Filter precursor peaks and peaks around precursor m/z.A logical "TRUE" or "FALSE" indicating if the spectra fragmented 
  					peaks around the precursor m/z +-20 Da should be deleted 
  					before the search comparisons. If "TRUE" this removes the 
  					residual precursor ion, which is frequently observed in MS/MS 
  					spectra acquired on qTOFs. (default: "TRUE")
- *\-w, \-\-window\_filter* [x]  :      If "TRUE", for each peak, it will check a window around
                                     that peak, if it is not one of the top peaks in terms of
                                     intensity in that window, it will be filtered out. This
                                     will speed up the search and reduce the effect of noise.
                                     (default: "TRUE")
- *\-\-analog\_search* [x]  :    If "TRUE", also search for analog spectra. This allows
                                     as a match the spectra with different precursor mass and
                                     similar fragmentation pattern.  (default: "FALSE")
- *\-\-analog\_max\_shift* [x]  :       Maximum difference between precursor m/zs that will be
                                     used in the search when analog_search is enabled. Only
                                     used when search_tool is "gnps_new". (default: 400)
- *\-l, \-\-parallel\_threads* [x]  :   the number of threads to be used for parallel processing
                                     in the search when search_tool is gnps_indexed. For
                                     parallelization set x >= 1 (default: 8)
- *\-c, \-\-count\_file\_path* [path]   :    Path to any of the count tables (peak_area or spectra)
                                     resulting from the NP3 clustering or clean steps (which
                                     matches with the input_mgf_file used). If the peak_area
                                     is informed and the spectra table file exists in the
                                     same path (or the opposite), it will merge the GNPS
                                     results to both files. Ignore if this is not a NP3
                                     result (leave as empty string). (default: "")
- *\-m, \-\-metadata* [file]  :  path to the metadata table CSV file of the NP3 job. This is necessary to plot the distribution of the superclasses grouping by sample, it may be missing if this plot is not desired (leave as empty string).
   (default: "")
   
## Results

Creates a directory named "identifications", inside the provided `output_path`, to store the results. 

Three tables will be stored in this directory: 

1. One containing all the `top_k` library matches with no annotations named with the `input_mgf_file` basename 
concatenated with the library file name (ALL\_GNPS\_NO\_PROPAGATED\_IonMode\_<ion\_mode\>.mgf) and the search tool used;
2. Another containing table 1 enriched with GNPS2 annotations named with the `output_path` basename plus the tag 
"library_search" and the search tool used;
3. And the final result table with the top 1 result named with table 2 name plus a suffix equals "top1".

A log file with the searching outputs is also created, named "logGNPS2LibrarySearch".

## Examples

Example identifying the result of one of the test cases using the gnps_new search tool.

```{ .text .copy }
node np3_workflow.js gnps_library_search -g 
test/L754_bacs/L754_bacs_blanks_one_sample/outs/L754_bacs_blanks_one_sample/mgf/L754_bacs_blanks_one_sample_clean.mgf 
-o test/L754_bacs/L754_bacs_blanks_one_sample/outs/L754_bacs_blanks_one_sample/ -i gnps_new
```

Now, example of identifying a MGF data in the negative ion mode:

```{r, eval=F} 
node np3_workflow.js gnps_library_search -a 2 -g data_in_negative_ion_mode.mgf 
-o output_path/data/library_search_result/ -i gnps_new
```

## Details

The GNPS2 Library Search workflow was adapted and integrated to the NP³ MS Workflow. All GNPS2 LC libraries, 
aggregated in the [ALL_GNPS_NO_PROPOGATED](https://external.gnps2.org/gnpslibrary) library, are used in this routine. 
It's MGF is retrieved locally by the **setup** command and is separated in two MGF files, one containing data in positive 
ion mode and another in negative ion mode only. Then, the setup creates a table with this MGF library ion headers 
information and enriches this table with GNPS2 online annotations previous organized from April 2026 (all offline). 
This enrichment is followed by a validation routine of the retrieved SMILES to maintain only the valid ones that can 
be processed by Cytoscape using ChemViz2 and to remove double quotes from the compound names. 
If a SMILES present in the MGF is validated it will be used in the entries where the online SMILES was invalid. 

The GNPS2 Library Search annotations will be updated with the GNPS2 online data from time to time by the admins of 
this repository and at that moment the *setup* command will need to be executed again to retrieve the new updates. 

Three search tools from [GNPS2](https://gnps2.org/homepage) are used here for the local identifications: gnps_indexed (current default in GNPS2),
gnps and gnps_new. The following filters that are applied in GNPS2 Library Search workflow were also adapted here to 
use the offline annotations enrichment. The top1 selection was adapted to resolve ties in the MQScore by using the 
number of matched peaks and the library quality string of the retrieved matches (different from GNPS2 online which only 
uses the MQScore).

The ALL_GNPS_NO_PROPAGATED library contains 956358 experimental spectra MS/MS with distinct SpectrumID's, from which 
862848 (90%) spectra have a valid SMILES from 90253 unique molecules, which counts the unique SMILES present in the 
annotations. The data in positive ion mode (ALL_GNPS_NO_PROPAGATED_IonMode_Positive) contains 713858 (75%) spectra and 
the data in negative ion mode (ALL_GNPS_NO_PROPAGATED_IonMode_Negative) contains 242500 (25%) spectra. 
Additional disk space is necessary to store this library offline, around 4Gbs are used.

The GNPS2 data and the codes adapted and modified from GNPS2 are located in: '[src/GNPS_LibrarySearch](https://github.com/danielatrivella/NP3_MS_Workflow/tree/master/src/GNPS_LibrarySearch)'; 

The enriched annotations are located in:
'[GNPS_LibrarySearch/data/library_summary_ALL_GNPS_NO_PROPAGATED_annotations.tsv.tar.gz](https://github.com/danielatrivella/NP3_MS_Workflow/blob/master/src/GNPS_LibrarySearch/data/library_summary_ALL_GNPS_NO_PROPOGATED_annotations.tsv.tar.gz)'. 

The retrieved libraries MGF files will be stored in:
'[GNPS_LibrarySearch/data/](https://github.com/danielatrivella/NP3_MS_Workflow/tree/master/src/GNPS_LibrarySearch/data)libraries'.

## GNPS Identification Curation and UNPDxGNPS Final Curation Statistics

A curation of the GNPS identification results is also performed for the gnps_library_search and gnps_result in order 
to score and rank the best library annotations from the GNPS libraries. The curation is executed after the gnps result 
is merged to the NP³ count tables.

The curation of the GNPS identification result start by comparing the SMILES from the tremolo_SMILES_best and the 
gnps_Smiles using Tanimoto, and stores the similarity scores in the 'tanimoto_unpd_gnps' column.

Next, it creates the column 'gnps_GoldCategory' and set it as 'GOLD' for the spectra that have gnps_MZErrorPPM <= 20 
and gnps_MQScore >= 0.9, otherwise set it as 'out'. The spectra set as 'GOLD' contains a very reliable identification.

Then, after the gold results were selected, the curation proceeds to a categorization and scoring of the GNPS result, 
similarly to the UNPD Identification Curation. Ten categories and scores were defined to group the best results 
according to the reliability of their matches or to put them out of the analysis, using the criteria presented below:

| GNPS Category | Score | Criteria | 
| :--------------- | :---: | :------------------------------------: | 
| Top1\|mqs>=0.9 | 100 | mzE <= 20 and nsP >= 6 |
| Top2\|0.8<=mqs<0.9 | 80 | mzE <= 20 and nsP >= 6 |
| Top3\|sp<6 mqs>=0.9 | 60 | mzE <= 20 and nsP < 6 |
| Top4\|sp<6 0.8<=mqs<0.9 | 50| mzE <= 20 and nsP < 6 |
| Top5\|mqs>=0.7 | 45 | mzE <= 20 and nsP >= 6 and mqs<0.8 |
| Analog1\|mqs>=0.9 | 40 | mzE > 20  and nsP >= 6 |
| Analog2\|0.8<=mqs<0.9 | 30 | mzE > 20 and nsP >= 6 |
| Analog3\|sp<6 mqs>=0.9 |  20 | mzE > 20 and nsP < 6 |
| Analog4\|sp<6 0.8<=mqs<0.9 | 10 | mzE > 20 and nsP < 6 |
| Analog5\|mqs>=0.7 | 5 | mzE > 20 and nsP >= 6 and mqs<0.8 |
| out | 0 | did not passed in any category - unreliable identification |
  
* Where mqs is the gnps_MQScore, nsP is the gnps_SharedPeaks and mzE is the gnps_MZErrorPPM values. 

The defined category for each GNPS identification is stored in the column 'gnps_category' and its corresponding score 
is stored in the column 'gnps_score'.

The superclass of the GNPS result is cleaned to only retrieve the NPClassifier or first result present (before the 
pipe '|') and its stored in the column 'gnps_npclassifier_superclass_clean'. This clean result is used to create the 
'gnps_curated_superclass' column, which receives the gnps_npclassifier_superclass_clean where gnps_score > 0 (or 
gnps_category != 'out'). Following, it is created the 'gnps_curated_superclass_grouping' column with the grouping of the 
gnps_curated_superclass and the 'gnps_curated_superclass_GR_<superclass group\>' columns checking the occurrence of the 
spectra in each superclass group. The superclass grouping follow the table defined in the [Identification Curation of 
UNPD-Tremolo result of Step 6](tremolo.md/#superclass-grouping).

Finally, the final identification curation is performed following [next section](final_identification_curation.md) 
to select the best library annotation from UNPDxGNPS, if the UNPD identification result is present, 
otherwise only the curated library annotation from GNPS is 
used. The columns starting with the prefix 'curated_lib_annotation_*' will store the final curated library annotation, 
their origin (UNPD or GNPS), score, quality, SMILES and superclass. 

#### UNPDxGNPS Chemical Statistics

At the end of the GNPS curation step, the final NP³ result will have the best curated identification from UNPD and GNPS. 
If the `job_output_path` is informed, the commands further creates two chemical reports with the GNPS and UNPDxGNPS 
identification statistics and stores them to the `job_output_path`/final_reports/chemical_reports folder. 
In the chemical report, this command also plots the NP³ chemical space with the curated GNPS identification and the 
origin of the final curated identification (UNPDxGNPS). 

For the chemical space plotting using PCA, the CDK descriptors of the GNPS result are previous computed and stored in 
the `job_output_path`/identifications folder in a CSV file named 'gnps_results_smiles_descriptorsCDK.csv'. The file 
'gnps_results_smiles.csv' is created from the GNPS joined identifications and is used for the CDK descriptors computation. 
The top 24 descriptors are computed here for each unique SMILES present in the gnps_Smiles column. 
Only the identifications from not blank and not bed m/zs (BLANKS_TOTAL == 0 and BEDS_TOTAL == 0) are used in the PCA 
plot, and a separated plot is created filtering only the protonated m/zs.

If the metadata path is informed, the plot with the distribution of the superclass grouping are also created for GNPS 
and the final curated identification (UNPDxGNPS) using columns 'gnps_curated_superclass_grouping' and 
'curated_lib_annotation_superclass_grouping', respectively, and are stored in the chemical report folder.

## References

[^1]: Wang, M., Carver, J., Phelan, V. et al. Sharing and community curation of mass spectrometry data with Global Natural Products Social Molecular Networking. Nat Biotechnol 34, 828–837 (2016). https://doi.org/10.1038/nbt.3597 (GNPS)