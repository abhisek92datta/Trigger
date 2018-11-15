# Trigger Efficiency Analysis for 2018 Data

# Installation :

source /cvmfs/cms.cern.ch/cmsset_default.sh

export SCRAM_ARCH="slc6_amd64_gcc630"

export CMSSW_VERSION="CMSSW_10_2_1"

cmsrel $CMSSW_VERSION

cd ${CMSSW_VERSION}/src

voms-proxy-init -voms cms

cmsenv

git clone https://gitlab.cern.ch/abdatta/Trigger.git

cd Trigger/

git checkout Data_2018

cd ..

scram b

cd Trigger/TriggerAnalyzer/


# Step 1 : Running analyzer on MiniAOD files (Data or MC) to produce ntuples

# Setting the Global Tag

For Prompt-Reco 2018 Data : in test/trigger_data_PromptReco.py

process.GlobalTag.globaltag = '101X_dataRun2_Prompt_v9'

For Re-Reco 2018 Data : 

process.GlobalTag.globaltag =

For 2018 MC : 

process.GlobalTag.globaltag = 

# To run locally :

Put desired filename in python config file , eg. test/trigger_data_PromptReco.py

cmsRun test/trigger_data_PromptReco.py > output_log.txt

# To run on the grid using CRAB

Set the request name, name of destination folder and the desired Dataset and JSON file

1. config.General.requestName = '< request name >'

2. config.Data.outputDatasetTag = '< output dataset tag >'

3. config.Data.inputDataset = '< Dataset Name >'

4. config.Data.lumiMask = 'data/< JSON filename >'  (only for DATA)

Modify these lines according to your storage area in the crabConfig_*.py files

1. config.Data.outLFNDirBase = '< output location path >'

2. config.Site.storageSite = 'T3_US_FNALLPC'

crab submit -c < crab config filename >  (for eg. < crab config filename > : crabConfig_Data_jetht_A_v1_promptreco.py)


# Step 2 : Running macros to obtain Trigger Efficiencies using the ntuples to produce root files with plots

Repeat the Installation steps (in FNAL account if ntuples stored in FNAL storage area)

cd macros/Trigger_Efficiency_Codes/

(These examples are only for 2017C_promptreco, you can do similarly for other Data or MC)

# For L1

cd L1/

cp ../../../interface/TriggerStudyEventVars.h .

First, put the correct path and number of paths for the relevant ntuples in input_L1_Global_SE_2017C_promptreco.txt and input_L1_Global_SM_2017C_promptreco.txt. Then,

cp input_L1_Global_SE_2017C_promptreco.txt input.txt

root -l -q Trigger_Efficiency_calc_L1_Global_SE.C > out_l1_global_se_2017C_promptreco.txt

cp input_L1_Global_SM_2017C_promptreco.txt input.txt

root -l -q Trigger_Efficiency_calc_L1_Global_SM.C > out_l1_global_sm_2017C_promptreco.txt

rm -rf input.txt

Repeat this for all the different Datasets 

cd ..

# For HLT

cd HLT/

# For Global Efficiency

cd Global/

cp ../../../../interface/TriggerStudyEventVars.h .

First, put the correct path and number of paths for the relevant ntuples in input_HLT_Global_MET_2017C_promptreco.txt and input_HLT_Global_SM_2017C_promptreco.txt. Then,

cp input_HLT_Global_MET_2017C_promptreco.txt input.txt

root -l -q Trigger_Efficiency_calc_HLT_Global_MET.C > out_hlt_global_met_2017C_promptreco.txt

cp input_HLT_Global_SM_2017C_promptreco.txt input.txt

root -l -q Trigger_Efficiency_calc_HLT_Global_SM.C > out_hlt_global_sm_2017C_promptreco.txt

rm -rf input.txt

Repeat this for all the different Datasets 

cd ..

# For Jet/HT Leg Efficiency

cd Jet_HT_Leg/

cp ../../../../interface/TriggerStudyEventVars.h .

First, put the correct path and number of paths for the relevant ntuples in input_HLT_JetHTLeg_SE_2017C_promptreco.txt. Then,

cp input_HLT_JetHTLeg_SE_2017C_promptreco.txt input.txt

root -l -q Trigger_Efficiency_calc_HLT_JetHTLeg_SE.C > out_hlt_jethtleg_se_2017C_promptreco.txt

rm -rf input.txt

Repeat this for all the different Datasets 

cd ..

# For Electron Leg Efficiency

cd Ele_leg/

# For Jet Path

cd Jet_Path/

cp ../../../../../interface/TriggerStudyEventVars.h .

First, put the correct path and number of paths for the relevant ntuples in input_HLT_EleLeg_JET_2017C_promptreco.txt. Then,

cp input_HLT_EleLeg_JET_2017C_promptreco.txt input.txt

root -l -q Trigger_Efficiency_calc_HLT_EleLeg_JET.C> out_hlt_eleleg_jet_2017C_promptreco.txt

rm -rf input.txt

Repeat this for all the different Datasets 

cd ..

# For HT Path

cd HT_Path/

cp ../../../../../interface/TriggerStudyEventVars.h .

copy TriggerStudyEventVars.h from interface/ to this location

First, put the correct path and number of paths for the relevant ntuples in input_HLT_EleLeg_HT_2017C_promptreco.txt. Then,

cp input_HLT_EleLeg_HT_2017C_promptreco.txt input.txt

root -l -q Trigger_Efficiency_calc_HLT_EleLeg_HT.C> out_hlt_eleleg_ht_2017C_promptreco.txt

rm -rf input.txt

Repeat this for all the different Datasets 

cd ..

cd ..

Shortcut : You can use run_macros_promptreco.sh, run_macros_rereco.sh and run_macro_mc.sh (in macros/Trigger_Efficiency_Codes/) to run all the above Trigger Efficiency macros together instead one at a time

# Step 3 : Running macros on the root plots to produce pdf/png plots

cd ../../Trigger_Plot_Codes/

Collect all the root files produced by the macros to this location and then the following macros :

# For L1 plots

root -l -q Compare_Efficiencies_L1.C

# For HLT plots

root -l -q Compare_Efficiencies_HLT.C

# For combining relevant plots in one root file

root -l -q Combine_plots_root.C



# Additional CRAB info

# to check status of CRAB jobs
crab status -d < crab_folder >

# to kill CRAB jobs
crab kill -d < crab_folder >



# Additonal GIT commands

# to clone repo and swtich to desired branch
git clone < GITLAB repo link >

git checkout < branch name >

# to create a new branch
git checkout -b < branch name >

# to check current branch
git branch

# to pull changes from GITLAB repo to local 
git pull

# to push changes to GITLAB repo (always pull before pushing)
git add .

git commit -m “Comment”

git push origin < branch name >


# Adding together histograms/Tefficiencies in root files
hadd Combined_Root_file.root Individual_Root_files_*.root


