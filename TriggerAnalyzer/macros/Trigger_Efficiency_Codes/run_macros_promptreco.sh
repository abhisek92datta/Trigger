cd L1/
cp ../../../interface/TriggerStudyEventVars.h .

cp input_L1_Global_SM_2018A_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_L1_Global_SM.C > out_l1_global_sm_2018A_promptreco.txt

cp input_L1_Global_SM_2018B_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_L1_Global_SM.C > out_l1_global_sm_2018B_promptreco.txt

cp input_L1_Global_SM_2018C_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_L1_Global_SM.C > out_l1_global_sm_2018C_promptreco.txt

cp input_L1_Global_SM_2018D_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_L1_Global_SM.C > out_l1_global_sm_2018D_promptreco.txt

hadd Efficiency_Trigger_histos_L1_Global_SM_2018ABCD_Promptreco.root Efficiency_Trigger_histos_L1_Global_SM_2018*_Promptreco.root

echo "L1 SM Done"

cp input_L1_Global_MET_2018A_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_L1_Global_MET.C > out_l1_global_met_2018A_promptreco.txt

cp input_L1_Global_MET_2018B_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_L1_Global_MET.C > out_l1_global_met_2018B_promptreco.txt

cp input_L1_Global_MET_2018C_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_L1_Global_MET.C > out_l1_global_met_2018C_promptreco.txt

cp input_L1_Global_MET_2018D_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_L1_Global_MET.C > out_l1_global_met_2018D_promptreco.txt

hadd Efficiency_Trigger_histos_L1_Global_MET_2018ABCD_Promptreco.root Efficiency_Trigger_histos_L1_Global_MET_2018*_Promptreco.root

echo "L1 MET Done"

rm -rf input.txt
cd ../HLT/Global/
cp ../../../../interface/TriggerStudyEventVars.h .

cp input_HLT_Global_SM_2018A_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_HLT_Global_SM.C > out_hlt_global_sm_2018A_promptreco.txt

cp input_HLT_Global_SM_2018B_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_HLT_Global_SM.C > out_hlt_global_sm_2018B_promptreco.txt

cp input_HLT_Global_SM_2018C_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_HLT_Global_SM.C > out_hlt_global_sm_2018C_promptreco.txt

cp input_HLT_Global_SM_2018D_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_HLT_Global_SM.C > out_hlt_global_sm_2018D_promptreco.txt

hadd Efficiency_Trigger_histos_HLT_Global_SM_2018ABCD_Promptreco.root Efficiency_Trigger_histos_HLT_Global_SM_2018*_Promptreco.root

echo "HLT Global SM Done"

cp input_HLT_Global_MET_2018A_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_HLT_Global_MET.C > out_hlt_global_met_2018A_promptreco.txt

cp input_HLT_Global_MET_2018B_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_HLT_Global_MET.C > out_hlt_global_met_2018B_promptreco.txt

cp input_HLT_Global_MET_2018C_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_HLT_Global_MET.C > out_hlt_global_met_2018C_promptreco.txt

cp input_HLT_Global_MET_2018D_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_HLT_Global_MET.C > out_hlt_global_met_2018D_promptreco.txt

hadd Efficiency_Trigger_histos_HLT_Global_MET_2018ABCD_Promptreco.root Efficiency_Trigger_histos_HLT_Global_MET_2018*_Promptreco.root

echo "HLT Global MET Done"

rm -rf input.txt
cd ../Jet_HT_Leg/
cp ../../../../interface/TriggerStudyEventVars.h .

cp input_HLT_JetHTLeg_EG_2018A_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_HLT_JetHTLeg_EG.C > out_hlt_jethtleg_eg_2018A_promptreco.txt

cp input_HLT_JetHTLeg_EG_2018B_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_HLT_JetHTLeg_EG.C > out_hlt_jethtleg_eg_2018B_promptreco.txt

cp input_HLT_JetHTLeg_EG_2018C_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_HLT_JetHTLeg_EG.C > out_hlt_jethtleg_eg_2018C_promptreco.txt

cp input_HLT_JetHTLeg_EG_2018D_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_HLT_JetHTLeg_EG.C > out_hlt_jethtleg_eg_2018D_promptreco.txt

hadd Efficiency_Trigger_histos_HLT_JetHTLeg_EG_2018ABCD_Promptreco.root Efficiency_Trigger_histos_HLT_JetHTLeg_EG_2018*_Promptreco.root

echo "HLT JetHT Leg Done"

rm -rf input.txt
cd ../Ele_Leg/HT_Path/
cp ../../../../../interface/TriggerStudyEventVars.h .

cp input_HLT_EleLeg_HT_2018A_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_HLT_EleLeg_HT.C > out_hlt_eleleg_ht_2018A_promptreco.txt

cp input_HLT_EleLeg_HT_2018B_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_HLT_EleLeg_HT.C > out_hlt_eleleg_ht_2018B_promptreco.txt

cp input_HLT_EleLeg_HT_2018C_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_HLT_EleLeg_HT.C > out_hlt_eleleg_ht_2018C_promptreco.txt

cp input_HLT_EleLeg_HT_2018D_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_HLT_EleLeg_HT.C > out_hlt_eleleg_ht_2018D_promptreco.txt

hadd Efficiency_Trigger_histos_HLT_EleLeg_HT_2018ABCD_Promptreco.root Efficiency_Trigger_histos_HLT_EleLeg_HT_2018*_Promptreco.root

echo "HLT Ele Leg HT path Done"

rm -rf input.txt
cd ../Jet_Path/
cp ../../../../../interface/TriggerStudyEventVars.h .

cp input_HLT_EleLeg_JET_2018A_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_HLT_EleLeg_JET.C > out_hlt_eleleg_jet_2018A_promptreco.txt

cp input_HLT_EleLeg_JET_2018B_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_HLT_EleLeg_JET.C > out_hlt_eleleg_jet_2018B_promptreco.txt

cp input_HLT_EleLeg_JET_2018C_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_HLT_EleLeg_JET.C > out_hlt_eleleg_jet_2018C_promptreco.txt

cp input_HLT_EleLeg_JET_2018D_promptreco.txt input.txt
root -l -q Trigger_Efficiency_calc_HLT_EleLeg_JET.C > out_hlt_eleleg_jet_2018D_promptreco.txt

hadd Efficiency_Trigger_histos_HLT_EleLeg_JET_2018ABCD_Promptreco.root Efficiency_Trigger_histos_HLT_EleLeg_JET_2018*_Promptreco.root

echo "HLT Ele Leg JET path Done"

rm -rf input.txt
cd ../../../











