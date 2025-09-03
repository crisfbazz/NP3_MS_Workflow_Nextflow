#  NP³ MS Workflow Nextflow 

To run the workflow to test simply do

```
make run
```

The NP3 testing dataset is used by default.

To learn more about the NP3 pipeline, please checkout the [repository](https://github.com/danielatrivella/NP3_MS_Workflow).

To update the NP³ code with the latest commit and setup do

```
make update_np3
```

To learn NextFlow checkout this documentation:

https://www.nextflow.io/docs/latest/index.html


## Installation NP³ MS Workflow with Nextflow

You will need to have conda and mamba installed to run things locally. 

Create the np3_nextflow conda environment stored in: bin/environment_np3_nextflow_unix.yml

```
mamba env create -f bin/environment_np3_nextflow_unix.yml 
```

Then, install further nodejs dependencies and execute the NP³ setup automatically by running:

```
sh bin/install_setup.sh
```

You are good to go!

Or, the last installation step can be executed manually with the following steps:

Activate the np3_nextflow environment:

```
mamba activate np3_nextflow 
```

Then, install some additional nodejs packages and audit the installation: 

```
npm install shelljs@0.8.4 commander@5.1.0
npm audit fix
```

Finally, go to the 'bin/NP3_MS_Workflow' directory and run the NP³ MS Workflow setup (check installed packages, compile some code, aggregate some data from UNPD-ISDB and download large data from spec2vec):

```
cd bin/NP3_MS_Workflow
node np3_workflow.js setup 
```

You are good to go!

## Deployment to GNPS2

Check [Nexftlow template instructions from Mingxun Wang](https://github.com/Wang-Bioinformatics-Lab/Nextflow_Workflow_Template).

