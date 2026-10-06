# Step 7: *annotate_protonated* Details

## Ionization Annotation Algorithm
 
During the detection of a given metabolite in a MS experiment some ionization variants can occur, what adds more 
redundancy and ambiguity to the list of consensus spectra. To try to overcome this limitation of the technique a 
chemical annotation based on chemical rules and numerical equivalences is performed to detected some of the variations 
in the clean set of consensus spectra and creates the IVAMN with the results.
 
The ionization variants considered are adducts, neutral losses, multiple charge, multiple charge isotopic patterns, 
dimers/trimers, isotopes and in-source fragmentation. Adducts and dimers/trimers combined with neutral losses are also 
considered. Multiple charge isotopic pattern is considered looking at the MS2 and the MS1 list that do not have a 
consensus spectrum assigned to it, in the last case the detection is only signalize with a flag added to clean count tables.
 
The annotation algorithm considers every consensus spectrum as a putative [M+H]+ candidate and parses the spectra that 
are in concurrent chromatographic peaks (same retention time interval with a tolerance applied in the peak boundaries) 
in search of m/z shifts corresponding to known annotations (numerical rules with a tolerance), and with the use of 
additional chemical rules in some cases (see table below for more details). Only the consensus spectra that appears in 
at least one data collection batch in common can be annotated. This algorithm also computes the m/z error of the 
numerical rules (within the tolerance) and the retention time error (between the peak boundaries – retention time minimum 
and maximum) associated to each annotation. And two flags named 'multicharge_ion' and 'isotope_ion' are created. 
The 'multicharge_ion' flag signalize potential multiple charge ions without the need of the mono charge ion detection 
(which can be absent), it indicates if the ion has a multiple charge isotopic pattern (a m/z shift of 0.5 for double 
charge and of 0.33 for triple charge variants) detected in the MS1 or MS2 lists. The 'multicharge_ion' flag receives a 
value of 2 if the isotopic pattern and the mono charge variant were detected in the MS1 or MS2 list, 1 if the isotopic 
pattern was detected but the mono charge is absent, and 0 otherwise. The 'isotope_ion' flag signalize potential isotope 
ions, it receives a value of 1 if the respective spectrum was annotated as a isotope ion, and 0 otherwise.
 
The annotation algorithm guides its behavior using the provided rules table, where the user specifies in the correct 
format (described below) all ionization variants that should be detected. Except for in-source fragmentation and multiple 
charge isotopic patterns that are always considered. It's important to notice that the annotations are based on chemical 
evidence and numerical equivalence of the used rules, and thus they can be wrong in some cases or some true annotations
that are not covered could be missing. A fine tune of the tolerance parameters (for precursor m/z, retention time deviation 
and similarity cutoff) and of the defined annotations in the rules table should be manually performed by the user to reduce 
false positive and false negative rates, according to the experimental conditions.

#### Ionization Variants Rules Table Format

The NP³ MS workflow default rules table of accepted ionization modifications can be found inside the 'rules' folder 
in the '[np3_modifications.csv](https://github.com/danielatrivella/NP3_MS_Workflow/blob/master/rules/np3_modifications.csv)' CSV file. It must contain the adducts, neutral losses, multiple charge, dimers/trimers 
and carbon isotopes variants to be considered and which adducts or dimers/trimers should also be considered in 
combination with a neutral loss of water or ammonia. A copy of this file should be used in order to select a subset of 
the default rules and/or to add new ones. These new created file must follow the rules format described in the two 
tables below and be passed to the workflow commands to be used in place of the default rules table.
 
The NP³ MS workflow rules table must be defined following the format of the table below with the same column 
names, where *X*, *Y* and *Z* are values provided by the user:
 

|ion|mzdiff| charge | neutral_loss_h2o | neutral_loss_nh3 | neutral_loss |                                                                                                                                           sim_cutoff                                                                                                                                           |
|:--------------------:|::|:-----:| :--------: | :--------: | :--------: |:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
| [**Z**M+**X**]**Y**+ |   The m/z difference to the monoisotopic m/z ([M]+) after taking *Y* and *Z* into account   |+Y| 0 or 1, if 1 and the ion is an adduct or dimer/trimer variant also considerates it combined with a neutral loss of water (i.e., *X*-*H2O*), else nothing happen | 0 or 1, if 1 and the ion is an adduct or dimer/trimer variant also considerates it combined with a neutral loss of ammonia (i.e., *X*-*NH3*), else nothing happen | 0 or 1, if the ion modification is a neutral loss set it as 1, otherwise as 0 (if it is an adduct or dimer/trimer combined with a neutral loss set as 0) | A numeric ranging from 0 to 1. The similarity cut-off value to accept the annotation, for adducts, neutral losses and isotopes. Multiple charge variants do not rely on the similarity value. The spectra must have a similarity value greater or equal to this value to be annotated \newline |
| [**Z**M-**X**]**Y**- | The m/z difference to the the monoisotopic m/z ([M]-) after taking *Y* and *Z* into account |-Y| 0 or 1 | 0 or 1 | 0 or 1 |                                                                                                                                 A numeric ranging from 0 to 1                                                                                                                                  |

 
The 'ion' column will be used to label the detected annotations, spaces are removed and ignored. *X* must be a string 
with the ionization modification to be considered (e.g. for adducts Na, NH4, 2K-H or for neutral losses H-H2O, H-NH3, 
H-H2O-NH3). *Y* is the ion charge, it is used to define double or triple charge variants and when *Y=1* it must be 
omitted in the 'ion' column. *Z* is the number of repetitions of a monomer, used to define dimers or trimers variants, 
otherwise it must be omitted in the 'ion' column. The mzdiff column is defined in the same way as the 'Mass' column 
present in the Fiehn Lab Mass Spectrometry Adduct Calculator [12] and the numerical rules follow their standard. 
The difference in the NP³ MS workflow annotations is that the calculations are done assuming that the current spectra 
being annotated is a [M+H]+ or a [M-H]- ion, and so the standard ionization mode is hard coded to be subtracted from 
the spectra m/z, but the rules definition in the table are the same.
 
The following table describes the chemical and numerical rules used by the annotation algorithm to detect each type of 
possible ionization variant and the correct format of their annotation label (to be placed in the 'ion' column).
 
Where:
 
- M : is the m/z of the current spectrum being annotated (considered as a [M+H]+ or a [M-H]- ion)
- m : is the m/z of any other spectra that appears in a concurrent retention time interval and in at least one data
collection batch in common with the current spectrum M 
- ion_mode : is the ionization mode informed by the user (+1 for [M+H]+ and -1 for [M-H]-) which is converted to the 
hydrogen mass equals 1.00783 or -1.00783 Da
- mz_tol : is the precursor m/z tolerance in Daltons informed by the user (default 0.025)
- fragment_tol : is the MS2 fragmented peaks m/z tolerance in Daltons informed by the user (default 0.025)
- m_H2O : is the mass of the water molecule equals 18.01056 Da
- m_NH3 : is the mass of the ammonia molecule equals 17.03052 Da
- m_iso : is the mass difference between the carbon isotopes 13C - 12C equals 1.0033 Da
- abs(**expression**) : is the function that returns the absolute value of the expression inside the parentheses, e.g., 
the positive value
 
#### Chemical and Numerical Rules for Ionization Variants Annotation

| Variant Type | Annotation Format | Rules | Description |
| :------: | :-------------: | :--------------------------: | :---------------: |
| Adduct | [M+X]+ or [M+X-H2O]+ or [M+X-NH3]+ | abs(M-ion_mode-m-m_X) <= mz_tol OR abs(M-ion_mode-m-m_X-m_H2O) <= mz_tol OR abs(M-ion_mode-m-m_X-m_NH3) <= mz_tol AND the spectra have a similarity value greater or equal than the 'sim_cutoff' value of the respective annotation \newline | m_X is the mass of the adduct X |
| Multiple Charge Isotopic Pattern | [M+X]Y+ isotopic | abs((M-ion_mode) - m) - 1/Y <= fragment_tol AND the peak area of the greater m/z is less than 2/3 of the other spectra peak area | Y is the charge of the multiple charge variant, only Y = 2 or 3 are considered. Multiple charge variants do not rely on the similarity value; This annotation is always considered \newline |
| Multiple Charge | [M+X]Y+ | ((M-ion_mode)/Y + m_X) - m <= mz_tol AND the presence of the isotopic pattern where an m/z equals m+1/Y must be present in the count tables (the not fragmented MS1 peaks count table is also checked) | m_X is the mass of the adduct X and Y is the charge of the multiple charge variant, only Y = 2 or 3 are considered. Multiple charge variants do not rely on the similarity value. \newline|
| Dimer or Trimer | [Z\*M+X]+ or [Z\*M+X-H2O]+ or [Z\*M+X-NH3]+| (Z\*(M-ion_mode) + m_X) - m <= mz_tol AND the monomer must not have fragmented peaks greater than the precursor m/z plus a gap for neutral losses AND one of the following is true (1) the dimer or trimer do not have fragmented peaks between the m/z of the monomer and its precursor m/z plus a gap for neutral losses and it has no fragmented peaks greater than its precursor m/z plus a gap for neutral losses OR (2) the candidates have a spectra similarity greater or equal than the 'sim_cutoff' value of the respective annotation \newline | m_X is the mass of the adduct X (less m_h2o or m_nh3 when X-H2O or X-NH3) and Z is 2 if it is a dimer and 3 if it is a trimer; this annotation is only considered if no multiple charge variant was detected between spectra M and m |
| In-source Fragmentation | fragment | m appears in the list of fragmented peaks of M, within a tolerance equals 'fragment_tol' AND m < M - 2\*iso_mass AND with an intensity greater than a MS2 baseline cut-off (15 is used by default) AND the spectra have a trimmed similarity value greater or equal than 0.2 | the fragmented peaks intensity are normalized from 0 to 1000 and then scaled. The scale is also applied to the baseline cut-off before using it; A trimmed similarity is computed by first trimming the fragmented peaks of both spectra using the smaller precursor m/z less 2*m_iso and then computing their cosine similarity, the maximum similarity value between the shifted cosine and the trimmed cosine is used; This annotation is always considered \newline |
| Neutral Loss | [M+H-X]+ | abs(M-m-m_X) <= mz_tol AND the spectra have a similarity value greater or equal than the 'sim_cutoff' value of the respective annotation \newline | m_X is the mass of the neutral loss X, e.g. if X=H2O then m_X = m_h2o. The 'neutral_loss' column must be set to 1 |
| Isotope | [M+X]+ | (M-m) - m_iso * X <= mz_tol AND the spectra have a similarity value greater or equal than the 'sim_cutoff' value of the respective annotation | X is an integer representing the carbon isotope deviation to be considered in the precursor m/z, e.g. X=1 for C13-C12 (m_iso difference) and X=2 for C14-C12 (2\*m_iso difference). We recommend setting X <= 2 and a high 'sim_cutoff' |

 
Only the rules for the positive ion mode are present in the above table, but the equivalent rules are also applied when the negative ion mode is selected but they were not vastly tested for this ion mode and thus could lead to undesired results.
 
The annotations of ionization variants detected for each msclusterID (row) being considered as a [M+H]+ ion is stored in the count tables in 6 new columns, one for each type of annotation, except for neutral losses that are grouped with adducts, named as: 'adducts', 'isotopes', 'dimers', 'multiCharges', 'fragments' and 'analogs.' Where 'analogs' annotate the concurrent spectra that have a similarity value greater or equal than 0.7 and did not receive any annotation from the spectra of the respective row.  
 
The annotations are stored in those columns in the following format:
 
>  *<ann\>* (sim *<sim_value\>* - mzE *<mz_error\>* - rtE *<rt_error\>*)[*<variantID\>*]\{*<#samples\>*\}
 
Where,
 
> **ann** is the detected annotation, following the format present in Table 4
 
> **sim_value** is the spectra similarity value. If the annotation is 'fragment' the similarity value is the 
trimmed spectra similarity.
 
> **mz_error** is the m/z error related with the numerical rule of the respective annotation
 
> **rt_error** is the retention time error between the spectra peak boundaries
 
> **variantID** is the variant spectra msclusterID in the clean count table
 
> **#samples** is the number of common samples between the spectra, the number of samples that both spectra appear
 
If more than one annotation of a same type is present for a given spectra, they are concatenated using a ';' character.
 
## IVAMN Creation and [M+H]+ Assignment Algorithms 
 
After the variants annotation algorithm the ionization variant annotation molecular network (IVAMN) is created based on the detected chemical annotations. The nodes of this network are the clean consensus spectra and the links of this network connects two consensus spectra that have a chemical annotation. The links have a direction, pointing from the consensus spectra considered as an ion variant to the consensus spectra considered as a putative [M+H]+ candidate in the respective annotation (e.g., [M+Na]+ -> [M+H]+). 
 
The IVAMN .selfloop file have the following columns, comma separated:
 
- "msclusterID_source" : the number of the msclusterID of the source node (ionization variant)
- "msclusterID_target" : the number of the msclusterID of the target node (considered [M+H]+ ion)
- "cosine" : the similarity value between the source and the target node
- "annotation" : the annotation detected between the source and the target node
- "mzError" : the m/z numerical rule error associated with the detected annotation
- "rtError" : the retention time error between the peak boundaries of the source and target nodes
- "numCommonSamples" : the number of samples that both the source and the target node appear
- "componentIndex" : the component index of the source and target nodes
 
Finally, an algorithm to assign some of the clean consensus spectra as the most likely [M+H]+ (protonated) candidates is applied in the IVAMN. This algorithm performs a link analysis in IVAMN to assign some nodes as representative [M+H]+ in each component (set of interconnected nodes) of this network. It works as follows: 
 
1. First it runs the PageRank [13], a link analysis algorithm; 
2. Then, for each component of IVAMN it selects the node with the highest score in the PageRank algorithm as a putative [M+H]+ representant. 
    - The nodes signalized as possible multiple charge ions are excluded from this selection, as well as ions annotated as neutral losses of a multiple charge node. 
    - In the case of a tie in the PageRank score between nodes that are in a same cycle (e.g., that have links pointing to each other), the node pointed by the link of the annotation with the lowest m/z error is used to select the putative [M+H]+ representant of the current iteration. If the tie is between nodes that are not in a same cycle, the node with the lowest ID is selected first; 
3. Next it removes all the nodes ancestors to the node selected as representant (e.g., all nodes that have a path to this node) from the current component 
4. Repeats steps 2. and 3. until all nodes of the current component are removed (have a putative [M+H]+ representant that covers it). 
 
A validation metric is computed for each node assigned as a representative [M+H]+ candidate, equals the sum of the numerical rules' errors of the annotations of all its ancestors. The result of this algorithm is a list of consensus spectra that represent the most likely [M+H]+ candidates of the collection, and thus, can be used as a suggestion of the number of real metabolites present in the samples. Some ambiguities may occur in the assignment of the putative [M+H]+ representants in the components of the molecular network of annotations, due mostly to possible missing/mistaken annotations and the PageRank solution, and a more robust solution to this limitation will be addressed in the next version of this workflow. 
 
#### New Columns from the Protonated Assignment

The following columns are added to the count tables in the *annotate_protonated* step:
 
| Columns | Description | Value Type |
| :--------------------- | --------------------------------------- | :-------: |
| protonated_representative | 1 if the respective spectrum was assigned as a [M+H]+ representative, 0 otherwise | numeric |
| protonated_mzError_sum | if protonated_representative is 1, equals the sum of the numerical rules' errors of the annotations of all its ancestors in the molecular networking of annotations, NA otherwise | numeric |
| protonated_rtError_sum | if protonated_representative is 1, equals the sum of the retention times errors of the annotations of all its ancestors in the molecular networking of annotations, NA otherwise | numeric |

 
An attributes table is also created for the molecular network of annotations containing the nodes information 
(m/z and retention times), the 'multicharge_ion' flag, the above columns and the following columns:
 
- "in_degree" : the number of incoming connections of the node
- "total_degree" : the number of incoming and outgoing connections of the node 
<!-- - "hub" : the HITS algorithm [14] score, it estimates the node value based on outgoing links -->
<!-- - "authority" : the HITS algorithm [14] score, it estimates the node value based on the incoming links -->
- "pagerank" : the PageRank algorithm scores. It computes a ranking of the nodes in the network based on the structure of the incoming links
- "number_ancestors" : the number of ancestor nodes of the current node
- "protonated_num_ancestors_edges" : the number of connections between all the ancestor's nodes of the current node