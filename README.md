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

Create the np3_nextflow environment stored in: bin/environment_np3_nextflow_unix.yml

```
mamba env create -f bin/environment_np3_nextflow_unix.yml 
```

And activate the np3_nextflow environment:

```
mamba activate np3_nextflow 
```

Then, install some additional nodejs packages: 

```
npm install shelljs@0.8.4 commander@5.1.0 
```

Finally, go to the 'bin/NP3_MS_Workflow' directory and run the NP³ MS Workflow setup (check installed packages, compile some code, aggregate some data from UNPD and download large data from spec2vec):

```
node np3_workflow.js setup 
```

You are good to go!


## Deployment to GNPS2

Check [Nexftlow template instructions from Mingxun Wang](https://github.com/Wang-Bioinformatics-Lab/Nextflow_Workflow_Template).

