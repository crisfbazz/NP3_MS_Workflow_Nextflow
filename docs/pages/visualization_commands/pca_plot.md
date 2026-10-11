# Command *pca_plot*
 
This command creates a PCA plot of a new data in the NP³ reference chemical space, composed by 
UNPD+DrugBank+Allosteric datasets (detailed in the [Final Reports](../cmd_run/final_reports.md#np3-chemical-space)). 

The procedure to create the PCA plots is to first compute the CDK descriptors of a 
provided list of SMILES, then to create the NP³ reference PCA and to transform these new data to the created 
chemical space. The resulting PCA is exported to a PNG image.

The provided new data table may contain a column defining the types/categories of each SMILES entry to be used to 
label the points in the PCA plot.

For processing very big SMILES, the user may manually increase the underline java maximum memory limit size to 4Gb or 
more (default to 2Gb). In the case of an error "Java Exception <no description because toString() failed>", before 
executing this command the user may run the following code in the command line: `export _JAVA_OPTIONS="-Xmx4g"`.
  
## Parameters

For the complete list of parameters run:
 
```{ .text .copy }
node np3_workflow.js pca_plot --help
```

Parameters for the command **pca_plot**:
 
- *\-t, \-\-table_path* <path\>       :  The path to a table in CSV format containing SMILES string in a column. Optionally, it may also contain the types/categories of each SMILES entry in another column. It must be comma "," separated.
- *\-s, \-\-smiles_column* <name\>       :  The name of the column in table_path containing the SMILES string.
- *\-d, \-\-data_type* [name]       : A name defining the type of the data to be used to label the
                              points or the name of a column from table_path containing the
                              types/categories names of each SMILES entry to be used to label
                              the points in the final PCA plot. This data_type must be
                              different from "UNPD" or "GNPS" and, in the case of a column
                              name, it must contain at most 8 different classes, otherwise
                              only the provided data_type name will be used to label the
                              points. These column values will go to the legend of the PCA plot.
                              (default: "new_data")
- *\-o, \-\-output_path* [path]    :          The path to the output directory where the final PCA plots will be stored, together with the descriptors table. If empty string
                              "" (default), the output_path is set to the table_path base
                              directory. (default: "")
- *\-n, \-\-output_name* [file]   : The name to be used to name the final plots, used as a prefix in the PCA plots naming. (default: "new_data")
- *\-h, \-\-help* : output usage information
 
## Results
 
In the output_path directory four files are created: 

1. One table with the computed CDK descriptors for each valid SMILES present in the provided table_path, named with the 
table_path name plus the suffix "_descriptorsCDK.csv";
2. Two PCA plots containing the reference NP³ dataset used to create the reference chemical space and the new 
transformed data. These plots are named as: "`output_name`_chemical_space_NP3_`data_type`_PCA_scores.png" and 
"`output_name`_chemical_space_NP3_`data_type`_PCA_biplot_components.png", this second plot is similar to the first plus 
the principal components;
3. The PCA quality of representation circle plot with the reference components named as 
"pca_quality_representation_cos2_NP3_reference.png".

The third plot result is always the same and depends only on the NP³ reference dataset, used to create the PCA 
chemical space.
 
## Examples
 
Example calling the creation of a PCA plot for an external data:
 
```{ .text .copy }
node np3_workflow.js pca_plot --table_path "/path/to/the/table/test_pca.csv" 
--smiles_column SMILES --data_type Libraries 
--output_path "/path/to/the/output/directory/" --output_name test_libraries
```
