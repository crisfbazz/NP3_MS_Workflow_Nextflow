`Currently the *windows installation is not working* in all PCs, we recommend users to use Linux instead.`

## Installation with Conda or Mamba

 `We recommend using Mamba because it's considerably faster, but if you want to use Conda, just replace commands starting with 'mamba' by 'conda'.`

`System requirements: Packages tzdata, zip, wget and which must be installed.`

NP³ MS Workflow repository includes a conda/mamba environment file for Unix and Windows to help the users install package dependencies.  



First, download or clone the workflow repository. 

**Download** the repository:

- Click on the green button named '[Code](https://github.com/danielatrivella/NP3_MS_Workflow/archive/refs/heads/master.zip)' and then click on the 'Download ZIP' option.
    - The repository contains the entire Universal Natural Products Database (UNPD) and the GNPS_ALL LC libraries, so this download can take a while to finish (around 2 Gb). 
- Extract the zip file 

Or **clone** the repository:

- Make sure you have Git installed (instruction in https://github.com/git-guides/install-git)
- Open the command line in the folder you want to place this repository and then type:

```{ .text .copy }  
git clone https://github.com/danielatrivella/NP3_MS_Workflow.git
``` 

Then download and install **miniforge** to install the conda environment with the mamba package included from the following link:

- [Miniforge for Linux amd64](https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh)
- [Miniforge for Windows](https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe)
- [More OS and architectures](https://github.com/conda-forge/miniforge/releases)
    - During installation on Windows OS check the option: 'Add to my PATH environment variable'  
        - If you do not check this option you have to manually add the conda anda mamba executables to your PATH environment variables 
    - To install on Unix OS run the following terminal command in the folder where the Anaconda installer was downloaded 
        ```{ .text .copy }    
        sh Miniforge3-Linux-x86_64.sh
        ```      
        - During installation when asked 'Do you wish the installer to initialize Anaconda3 by running conda init?' answer 'yes' 
        - If you answer 'no', then you need to manually type the command `conda init` before using Mamba.
    - If you already have a conda environment installed on your device, you can install the mamba package exclusively in the base environment with the command below, but it is recommended that you do a fresh install removing the old conda environment.
        ```{ .text .copy } 
        conda install conda-forge::mamba -y
        ```

For more detailed instructions how to setup conda, visit the [conda-forge/miniforge's github](https://github.com/conda-forge/miniforge)
                   
To verify if the Conda and Mamba installation was successful, open a **new** terminal window and type:   

```{ .text .copy }   
mamba --version
``` 

If it outputs the current version of `mamba` or `conda` you are good to go.    

Now, go the NP³ MS Workflow repository folder that you have just downloaded and extracted, 
and open a terminal window there.  

Then, create the NP³ MS Workflow conda environment to automatically install almost all the workflow programs and 
required packages using the following terminal command:

Unix OS: 

```{ .text .copy }    
mamba env create -f environment_np3_unix.yml 
``` 

Windows OS: 

```{ .text .copy }   
mamba env create -f environment_np3_win.yml 
``` 

The NP³ MS Workflow conda environment must be created and activated before running the NP³ MS Workflow commands. If the NP³ MS Workflow environment was created successfully, activate it with the following command: 

```{ .text .copy }  
mamba activate np3 
```   
 
You should see the '**(np3)**' tag as the first thing in your terminal line. 
Every time you open a new terminal you must execute this command once again in order to activate the NP³ MS Workflow environment 
before executing the NP³ MS Workflow commands.
 
Now install the libraries required by the node.js program:
 
In the terminal execute the following command:
```{ .text .copy }   
npm install shelljs@0.8.4 commander@5.1.0
```
 
Make sure all programs and OS libraries were properly installed and continue to the remaining installation setup. 
If you could not install the environment or failed in any of the above steps go to the **manual** installation below.

-------------------------------------------------------

## Manual installation
For **manually** installing the NP³ MS Workflow dependencies, use the following links to download the required programs and install them. Then, go to the repository folder and run the following commands in the terminal to install the programs' required packages:
 
- *node.js* (LTS version) - https://nodejs.org/en/download/
    + Also need to install the *npm* package manager (automatically installed with node.js \- need to check the option in the installation setup for Windows OS)
    + Terminal command to install the *node.js* packages with *npm*:
    ```{ .text .copy } 
    npm install shelljs@0.8.4 commander@5.1.0
    ```
- *make* - https://www.gnu.org/software/make/
- Compilers:
    + *gcc* - https://gcc.gnu.org/
    + *g++* - https://gcc.gnu.org/projects/cxx-status.html
        + These compilers are usually installed for most Unix distribution.
- *R >= 3.6.3* - https://www.r-project.org/
    + Terminal command to install the required *R* packages:
    ```{ .text .copy }
    sudo Rscript src/R_requiriments.R
    ```
   <!-- $ R CMD javareconf -->
- *Python 3.7* - https://www.python.org/download/releases/3.7/
    + Also install *pip* - https://pypi.org/project/pip/
    + Terminal command to install the required *python 3* packages with *pip*:
    ```{ .text .copy } 
    pip install -r src/python_requirements
    ```
       
If any installation fails, look for dependencies problems and retry.

---------------------------------------------------

## Setup Overview 
 
The workflow command **setup** (more instructions in next subsection) checks the programs and packages installation and 
tries to automatically install the missing ones (except for the *node.js* packages). Then it configure the used libraries and models, ang compile the NP3_MSCluster algorithm.

If the workflow command **setup** fails the programs and packages must be manually installed, 
and then this command must be used to configure the libraries and to automatically compile the NP3_MSCluster algorithm. 
This compilation may also be executed manually (more instructions in the next subsection).
