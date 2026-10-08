# Step 10: *mn* Details
 
The molecular networks are a great way to visualize the clustering results. 
The nodes of the SSMN are the clean consensus spectra and the edges/links of this network connect two consensus spectra 
that have a similarity score greater or equal to the provided cosine similarity cut-off (default to 0.6 is used). 

The resulting .selfloops edges files can be 
opened in graph visualization tools, such as Cytoscape[^1] and Gephi[^2]. 
Then the count tables can also be loaded to the networks to add the other nodes information present in these tables, 
by matching their msclusterIDs.
 
## Methodology

The creation of the spectra similarity molecular network (SSMN) is based on the pairwise similarity comparison of the 
clean consensus spectra using the NP³ shifted cosine similarity score or the spec2vec scores (see Step 5 details), 
ranging from 0 (totally dissimilar) to 1 (completely identical). 
 
The SSMN is further filtered by removing the links between spectra with a number of common peaks that do not respect the 
minimum number of matched peaks cut-off (default to at least 6 matched peaks), 
limiting the number of edges/links of each node to the top K most similar neighbors (default K equals 15 was used) 
and pruning the remaining links with an increasing cosine similarity cut-off to limit the components sizes to a 
provided maximum number of nodes (default of 200 nodes is used). 
These filtering was based in the GNPS molecular networking analysis code[^3] with some adaptations.
 
If the resulting SSMN is too dense (many nodes and links forming a giant single component), the user can try to 
restrict the filters with a bigger similarity and minimum number of common peaks cut-offs 
(`similarity_mn` and `min_matched_peaks` parameters) and a smaller top K neighbors and component size cut-offs 
(`net_top_k` and `max_component_size` parameters).
 
#### SSMN columns

The SSMN .selfloop files have the following columns, comma separated:
 
- "msclusterID_source" : the number of the msclusterID of the source node spectrum
- "msclusterID_target" : the number of the msclusterID of the target node spectrum
- "cosine" : the similarity value between the source and the target node
- "num_matched_peaks": the number of common peaks between the source and target node spectra
- "annotation" : the existing annotation between the source and the target node (this can be missing for some node pairs)
- "num_peaks_source": the number of fragment peaks in the source node spectrum after normalization and cleaning
- "num_peaks_target": the number of fragment peaks in the target node spectrum after normalization and cleaning
- "componentIndex" : the component index of the source and target nodes
 
The SSMN is an undirected graph. The source and the target nodes could be shifted, their order do not matter but the 
source node has always a smaller msclusterID than the target node due to the implemented algorithm.

#### [M+H]+ Networks

At the end, if the IVAMN is present, a [M+H]+ analysis is executed and results in the protonated networks. 
It creates de SSMN [M+H]+ filtered and the IVAMN [M+H]+ networks without blanks and without culture media m/z - 
remove nodes that appear in a *blank* or *bed* sample.

The [M+H]+ analysis uses as input the complete SSMN, the IVAMN, the IVAMN attribute table (protonated information) and 
the clean table. 
It works as follows: 

1. In the IVAMN, the blank nodes (BLANKS_TOTAL > 0) are removed; if argument blank_expansion > 0, the blank neighbors 
and ancestors (or successors, ignore the links direction) are also removed together; 
otherwise only remove the blank nodes (default behavior);
1. Next, the culture media nodes (BEDS_TOTAL > 0) are removed from the IVAMN;
2. Then, the [M+H]+ (protonated_representative == 1) are selected in the remaining IVAMN, 
resulting in the **IVAMN [M+H]+ network**;
3. And the nodes from IVAMN [M+H]+ are used to select the final nodes from the complete SSMN, 
resulting in the SSMN [M+H]+ network;
4. Finally, the SSMN [M+H]+ is filtered using the `min_matched_peaks`, `net_top_k` and `max_component_size` parameters, 
resulting in the **SSMN [M+H]+ filtered**.

## References
[^1]: Lopes CT, Franz M, Kazi F, Donaldson SL, Morris Q, Bader GD. Cytoscape Web: an interactive web-based network browser. Bioinformatics. 2010 Sep 15;26(18):2347-8. doi: 10.1093/bioinformatics/btq430. Epub 2010 Jul 23. PMID: 20656902; PMCID: PMC2935447
[^2]: Bastian M., Heymann S., Jacomy M. (2009). Gephi: an open source software for exploring and manipulating networks. International AAAI Conference on Weblogs and Social Media
[^3]: Wang, M., Carver, J., Phelan, V. et al. Sharing and community curation of mass spectrometry data with Global Natural Products Social Molecular Networking. Nat Biotechnol 34, 828–837 (2016). https://doi.org/10.1038/nbt.3597 (GNPS)