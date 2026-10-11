# Command *spectra_viewer*

This command runs the Streamlit interactive Web App to visualize and compare MS2 spectra. 

It receives as input a MGF file or a peak list. Currently, it is only supported for Unix OS.

It is also possible to manipulate, filter, calculate similarity of the spectra and save to PNG or SVG plots.

## Parameters

For the complete list of parameters run:
 
```{ .text .copy }
node np3_workflow.js spectra_viewer --help
```

Parameters for the command *spectra_viewer*:

- *\-p, \-\-port* [port_number]  localhost server port number (default: "8501")
- *\-h, \-\-help*              output usage information

## Examples:

Example to execute the MS² spectra viewer application:

```{ .text .copy }
node np3_workflow.js spectra_viewer
```

Example to execute it in a different port:

```{ .text .copy }
node np3_workflow.js spectra_viewer --port 8080
```
