# Command *chr*
 
Runs the R interactive script to extract chromatogram(s) from raw MS1 data files (mzXML, mzData and mzML) and saves to 
PNG image(s) file(s). 

Depending on the provided parameters this can be a total ion chromatogram (TIC), a base peak chromatogram (BPC) or an 
extracted ion chromatogram (XIC) extracted from each sample/file. 

In the interactive prompt is possible to select a mz and/or retention time window to better visualize interesting parts 
of the sample's chromatograms. The plots can be grouped by data collection batch, by a selection of samples, or all the 
samples can be plotted together.
 
## Parameters

For the complete list of parameters run:
 
```{ .text .copy }
node np3_workflow.js chr --help
```

Parameters for the command *chr*:
 
- *\-n, \-\-data_name* [name] : the data collection name. Used for verbosity (default: -)
- *\-b, \-\-metadata* <file\>  : path to the metadata table CSV file
- *\-d, \-\-raw_data_path* <path\>   : path to the folder containing the input LC-MS/MS raw spectra data in mzXML, mzData and mzML format
- *\-h, \-\-help*                   output usage information
 
## Results
 
PNG image(s) file(s) with the extracted chromatogram(s) of the selected sample(s) and grouped as specified in the options.
 
## Examples
 
Fake example to execute the MS¹ Viewer command:

```{ .text .copy }
node np3_workflow.js chr --metadata "/path/to/the/metadata/file/test_np3_metadata.csv" 
--data_name "test_np3" --raw_data_path "/path/to/the/raw/data/directory"
```

```{ .text .copy } 
node np3_workflow.js chr -b "/path/to/the/metadata/file/test_np3_metadata.csv" 
-d "/path/to/the/raw/data/directory"
```