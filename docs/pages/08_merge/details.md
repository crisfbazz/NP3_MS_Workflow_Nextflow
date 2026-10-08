# Step 8: *merge* Details
 
The merge step extends the clean count tables with new rows containing symbolic clusters resulting from the aggregation 
by row of each spectra with its annotated variants. It contains all the clean results plus new symbolic spectra created 
from the merged annotations, with a quantification that could be more representative of a true metabolite (merged count 
tables). By default, the merge is only performed for the consensus spectra assigned as a [M+H]+ representative.
 
## Methodology

The merge algorithm is performed by consensus spectra, following these steps: first all annotated isotopes variants are
merged and results in new symbolic spectra, with their counts aggregated as in Step 4 (summing the number of spectra and 
summing the unique peak areas); then all annotated adducts and neutral losses are merged, the new symbolic spectra from 
the previous merge (isotopes) are also considered; next the same is done for dimers/trimers, multiple charge and finally 
for fragments. The merge can result in too many new symbolic clusters, for a single spectrum the maximum is 31 new 
symbolic clusters. We recommend to use the correlation score to filter the relevant new symbolic clusters and to set the 
parameter merge_protonated to TRUE to only apply the merge algorithm to a subset of the consensus spectra. 
The merge is not applied to spectra that appears in blank samples.
 
It's up to the user to identify relevant symbolic clusters. Since the annotations are based on evidences of the 
considered chemical rules, some symbolic clusters may include the aggregation of wrong annotations and must be 
evaluated carefully.
 
The 'peakLists' and 'peakInts' columns are also recomputed aggregating the merged spectra fragmented peaks lists 
similar to how it is done in the clean Step 5. The new symbolic spectra are not exported to a MGF file because their 
lists of fragmented peaks could be very noisy and were not vastly tested. In the next releases of the workflow, we plan 
to develop a more sophisticated heuristic to perform the merge and to provide a more relevant and easier to use result. 
Still, the new symbolic spectra can improve the biocorrelation score of metabolites that were detected as different 
types of ionization variants.
 
#### Merge Count Tables New Columns

The following columns are added to the count tables in the *merge* step:
 
| Columns | Description | Value Type |
| :--------------------- | --------------------------------------- | :-------: |
| mergedIDs_all | the msclusterIDs of all merged spectra concatenated by a ';'. The first msclusterID of this list is from the spectrum that received the merge | character |
| mergedIDs_adducts | the msclusterIDs of all adducts and neutral losses variants merged concatenated by a ';' | character |
| precursorMz_adducts | the mzConsensus of all adducts and neutral losses variants merged concatenated by a ';' | character |
| mergedIDs_dimers | the msclusterIDs of all dimers/trimers variants merged concatenated by a ';' | character |
| precursorMz_dimers | the mzConsensus of all dimers/trimers variants merged concatenated by a ';' | character |
| mergedIDs_multiCharges | the msclusterIDs of all multiple charge and multiple charge isotopic variants merged concatenated by a ';' | character |
| precursorMz_multiCharges | the mzConsensus of all multiple charge and multiple charge isotopic variants merged concatenated by a ';' | character |
| mergedIDs_fragments | the msclusterIDs of all fragments variants merged concatenated by a ';' | character |
| precursorMz_fragments | the mzConsensus of all fragments variants merged concatenated by a ';' | character |