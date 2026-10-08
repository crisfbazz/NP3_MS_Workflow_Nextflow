# Step 10: Command *mn* (Molecular Networking)
 
The mn command creates a spectra similarity molecular network (SSMN) which connects clean consensus spectra (the nodes) 
based on the pairwise spectra similarity value above a given similarity cut-off (edges/links defined by `similarity_mn`). 
And applies a filter in this network to keep only the most 
strongly connected spectra (higher similarity values) and this creates the SSMN filtered. 
The filtered SSMN contains components that represent the most analogous spectra, possible connecting spectra 
from similar chemical classes.

If the IVAMN is present, the [M+H]+ analysis is executed and results in the protonated IVAMN and protonated 
SSMN filtered. The ionization variants annotations from Step 7 are also included in the SSMN as labels when present.

The filters applied in the SSMN are used to remove edges/links between spectra that 
have less peaks in common than the minimum number of matched peaks (`min_matched_peaks`), 
to limit the number of neighbors of each node (number of edges/links) to the top K most similar ones (`net_top_k`) 
and to limit the size of the components to a maximum number of nodes (`max_component_size`).
 
## Parameters

For the complete list of parameters run:
 
```{ .text .copy }
node np3_workflow.js mn --help
```

Parameters for the command **mn**:
 
- *\-o, \-\-output_path* <path\>       : path to the output data folder, inside the outs directory of the clustering result folder. It should contain the 'molecular_networking' folder and inside it the 'similarity_tables' folder. The job name will be extracted from here
- *\-w, \-\-similarity_mn* [x]  : the minimum similarity score that must occur between a pair of consensus MS/MS spectra in order to create an edge in the molecular networking. Lower values will increase the component size of the clusters by inducing the connection of less related MS/MS spectra; and higher values will limit the components sizes to the opposite (default: 0.6)
- *\-\-min_matched_peaks* [x]           :  The minimum number of common peaks that two spectra must share to be connected by an edge in the filtered SSMN. Connections between spectra with less common peaks than this cutoff will be removed when filtering the SSMN. Except for when one of the spectra have a number of fragment peaks smaller than the given min_matched_peaks value, in this case the spectra must share at least 2 peaks. The fragment peaks count is performed after the spectra are normalized and cleaned. (default: 6)
- *\-k, \-\-net_top_k* [x]          : the maximum number of connections for one single node in the similarity molecular networking. An edge between two nodes is kept only if both nodes are within each other's [x] most similar nodes. This restriction is applied to the spectra following the msclusterID's order, so smaller m/z's are limited first. Keeping this value low makes very large networks (many nodes) much easier to visualize (default: 10)
- *\-x, \-\-max_component_size* [x]     :   the maximum number of nodes that all component of 
                      the similarity molecular network must have. The edges of 
                      this network will be removed using an increasing cosine 
                      threshold until each network component has at most X nodes. 
                      Keeping this value low makes very large networks (many nodes 
                      and edges) much easier to visualize. (default: 200)
- *\-\-blank_expansion* [x]      :   the distance of neighborhood nodes from the blanks in IVAMN to be 
  					selected for removal in the final protonated networks. (0) to only remove blanks nodes,
  					(1) to remove nodes directly connected to a blank node, 
  					(2 or greater) to remove nodes in a distance equal to 2 or greater from a blank node, 
  					or (-1) to remove all possible neighbours and ancestors of a blank node (remove blank clusters) from IVAMN  (default: 0)
- *\-b, \-\-max_chunk_spectra* [name]           : Maximum number of spectra (rows) to be loaded and processed in a chunk at the same time. In case of memory issues this value should be decreased (default: 3000)
- *\-v, \-\-verbose* [x]             : for values X\>0 show the scripts output information. (default: 0)
- *\-h, \-\-help*                    : output usage information
 
## Results
 
Four files are created inside the 'molecular_networking' folder with the SSMN, the filtered SSMN, the protonated 
IVAMN and the protonated SSMN filtered:
 
- The complete spectra similarity molecular network named as : \newline
'\<*output_name*\>\_ssmn_w\_\<*similarity_mn*\>.selfloops', which contains all links with a similarity value above the cut-off
- The filtered spectra similarity molecular network named as \newline
  '\<*output_name*\>\_ssmn_w_\<*similarity_mn*\>\_k_\<*net_top_k*\>\_x_ \newline
   \<*max_component_size*\>.selfloops'
- The IVAMN [M+H]+ named as: \newline
  '\<*output_name*\>\_ivamn_protonated.selfloops'
- The SSMN [M+H]+ filtered named as: \newline
  '\<*output_name*\>\_ssmn_protonated\_w\_\<*similarity_mn*\>\_k\_\<*net_top_k*\>\_x\_ \newline 
  \<*max_component_size*\>.selfloops'
    
Where the 'output_name' is extracted from the 'output_path';
 
## Examples
 
Fake example to execute the Molecular Network of spectra similarity command:

```{ .text .copy }
node np3_workflow.js mn --output_path "/path/to/the/output/dir/test_np3/outs/test_np3" 
--similarity_mn 0.7 --verbose 1
```
 
