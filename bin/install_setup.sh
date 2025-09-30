# install nodejs packages not present in conda and audit the installation
npm install shelljs@0.8.4 commander@5.1.0 
# install tzdata - internal R dependency
apt-get install -y tzdata 
# execute the np3 update and setup run
cd NP3_MS_Workflow && node np3_workflow.js setup && cd -
