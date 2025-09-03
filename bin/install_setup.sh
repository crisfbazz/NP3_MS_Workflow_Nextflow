conda activate np3_nextflow
# install nodejs packages not present in conda and audit the installation
npm install shelljs@0.8.4 commander@5.1.0 
npm audit fix
# execute the np3 update and setup run
cd bin/NP3_MS_Workflow && git pull origin master &&	node np3_workflow.js setup && cd -
