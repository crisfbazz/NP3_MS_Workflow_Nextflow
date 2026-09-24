## Setup Command

The NP³ MS Workflow **setup** command is used to **automatically**:

- Check packages installation
- Retrieve remaining dependencies
- Compile the NP3_MSCluster algorithm
- Unzip and configure UNPD in-silico library
- Retrieve and configure GNPS2_ALL LC libraries
- Retrieve spec2vec models

To execute it go to the repository folder, open a terminal window and run the following command:
 
```{ .text .copy } 
node np3_workflow.js setup
```
 
  + If it runs without any ERROR message you are good to go! Otherwise, look for dependencies problems and retry. If it persists to fail follow the **Manual installation** or contact the dev team.

## Manual Compilation

To manually **compile** the NP3_MSCluster algorithm, inside the NP³ MS workflow repository go to the *NP3_MSCluster* folder and in the terminal run:
 
```{ .text .copy }  
make clean
make
```
 
If it runs without any ERROR message you are good to go! Otherwise, look for dependencies problems and retry with superuser privileges.
 
--------------------------------------------------------------------------------
 
## Updates
 
After any update in the NP³ MS workflow, it will be signalized in the repository if the **setup** command needs to be executed again.
This may be necessary to recompile the NP3_MSCluster algorithm or configure the libraries or models to allow the new updates to work.