## Details

The final identification curation is performed at the end of the UNPD curation and at the end of the GNPS curation. 
It selects the origin (UNDP x GNPS) of the best library annotation present in the current NP³ result, apply a 
filter based on their curated score by origin and also set a quality class based on this score. 

In the final reports, the final curated library annotations are used for the final statistics and compared with the 
best identifications from UNPD and GNPS not curated (complete identification result). The final curated library 
annotations are stored in different columns with the prefix "curated_lib_annotation_*". 

The final identification curation starts by creating the column "curated_lib_annotation_origin" with the library origin 
(UNPD or GNPS) of the best library annotation with a score greater than 0. 
This column is set to 'GNPS' if the gnps_score > 0 and gnps_score >= tremolo_UNPD_score_best. If the tremolo result 
is not present, uses only gnps_score > 0 to set the best origin as 'GNPS', the rest is left empty. 
Otherwise, if the GNPS result is not present, set it to 'UNPD' where tremolo_UNPD_score_best > 0. 
The preference here goes to GNPS, which is an experimental spectra library.

Then, using the selected origin of each final curated identification, the procedure retrieves its respective 
score (from 'gnps_score' for 'GNPS' and from 'tremolo_UNPD_score_best' for 'UNPD') and stores it in the column 
"curated_lib_annotation_score". When no origin was selected, the score is set to 0. 
These scores are defined in Section 4.6.4 for UNPD and Section 4.12.4 for GNPS.

Using the score of the final curated identification, a **quality grouping** is defined and assigned to column 
"curated_lib_annotation_quality" using three quality groups: 

| Quality Group | Criteria | Description |
| :--: | :---------: | :---------: |
| 1 | curated_lib_annotation_score >= 50 | representing good/top quality library annotations |
| 2 | curated_lib_annotation_score > 0 and curated_lib_annotation_score < 50 | less reliable library annotations, mostly analogs or low similarity score |
| 0 | curated_lib_annotation_score == 0 | category 'out' - not annotated |

Next, the SMILES of the final curated identification is set to the column "curated_lib_annotation_SMILES", 
which receives the value from the column "gnps_Smiles" where origin equals to "GNPS" or the value from the 
column "tremolo_SMILES_best" where origin equals to "UNPD", otherwise leave it empty. Similarly, 
the superclass of the final curated identification is set to the column "curated_lib_annotation_superclass", 
which receives the value from the columns "gnps_curated_superclass" or "tremolo_curated_superclass"; 
the final curated ID is set to the column "curated_lib_annotation_ID", which receives the value from the 
columns "gnps_SpectrumID" or "tremolo_UNPD_IDs_best"; and the final curated compound name is set to the 
column "curated_lib_annotation_compoundName", which receives values from the columns "gnps_Compound_Name" or 
"tremolo_chemicalNames_best".

Finally, using the values present in the superclass of the final curated library annotation, the columns with 
the superclass grouping are created accordingly named "curated_lib_annotation_superclass_grouping" and 
"curated_lib_annotation_superclass_GR_<superclass_group\>", where superclass_group are the 10 superclass 
groups defined in the previous section plus the not annotated group.

The final curated identifications present in the columns "curated_lib_annotation_*" are defined as the more 
reliable library annotations of each consensus spectra from a NP³ result.
