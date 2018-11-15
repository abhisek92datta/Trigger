#ifndef Trigger_TriggerAnalyzer_TriggerStudyEventVars_h
#define Trigger_TriggerAnalyzer_TriggerStudyEventVars_h

//
// Dependencies (#includes)
//
#include <iostream>
#include <vector>
#include "TLorentzVector.h"

//#ifdef __MAKECINT__
//#pragma link C++ class std::vector< TLorentzVector >+;
//#endif

using namespace std;



typedef std::vector<std::vector<double> > vvdouble;
typedef std::vector<std::vector<std::string> > vvstring;
typedef std::vector<double> vdouble;
typedef std::vector<float> vfloat;
typedef std::vector<string> vstring;
typedef std::vector<bool> vbool;
typedef std::vector<int> vint;
typedef std::vector< TLorentzVector > vecTLorentzVector;

//
// Utility Class for Handling Event Variables
//

struct triggerStudyEventVars{


  //////////////////////////////////////////////////////////////////////////
  ///  Tree branches/leaves
  //////////////////////////////////////////////////////////////////////////

  explicit triggerStudyEventVars() { }

  /////

  int event_nr_;
  int run_nr_;
  int lumi_nr_;

  double rho_;

  // L1
  int pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_;
  int pass_L1_LooseIsoEG28er2p1_HTT100er_;
  int pass_L1_Mu_all_;
  int pass_L1_EG_all_;
  int pass_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3_;
  int pass_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3_;
  int pass_L1_LooseIsoEG24er2p1_HTT100er_;
  int pass_L1_LooseIsoEG26er2p1_HTT100er_;
  int pass_L1_LooseIsoEG30er2p1_HTT100er_;

  // HLT

  int pass_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_;
  int pass_HLT_Ele28_eta2p1_WPTight_Gsf_HT150_;
  int pass_HLT_Ele27_WPTight_Gsf_;
  int pass_HLT_Ele32_WPTight_Gsf_;
  int pass_HLT_Ele32_WPTight_Gsf_L1DoubleEG_;
  int pass_HLT_Ele35_WPTight_Gsf_;
  int pass_HLT_Ele35_WPTight_Gsf_L1EGMT_;
  int pass_HLT_Ele38_WPTight_Gsf_;
  int pass_HLT_Ele40_WPTight_Gsf_;
  int pass_HLT_IsoMu30_;
  int pass_HLT_IsoMu27_;
  int pass_HLT_IsoMu24_eta2p1_;
  int pass_HLT_IsoMu24_;
  int pass_HLT_PFMET110_PFMHT110_IDTight_;
  int pass_HLT_PFMET120_PFMHT120_IDTight_;
  int pass_HLT_PFMET130_PFMHT130_IDTight_;
  int pass_HLT_PFMET140_PFMHT140_IDTight_;
  int pass_HLT_PFHT180_;
  int pass_HLT_PFHT250_;
  int pass_HLT_PFHT370_;
  int pass_HLT_PFHT430_;
  int pass_HLT_PFHT510_;
  int pass_HLT_PFHT590_;
  int pass_HLT_PFHT680_;
  int pass_HLT_PFHT780_;
  int pass_HLT_PFHT890_;
  int pass_HLT_PFHT1050_;
  int pass_HLT_PFJet40_;
  int pass_HLT_PFJet60_;
  int pass_HLT_PFJet80_;
  int pass_HLT_PFJet140_;
  int pass_HLT_PFJet200_;
  int pass_HLT_PFJet260_;
  int pass_HLT_PFJet320_;
  int pass_HLT_PFJet400_;
  int pass_HLT_PFJet450_;
  int pass_HLT_PFJet500_;
  int pass_HLT_PFJet550_;

  int prescale_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_;
  int prescale_HLT_Ele28_eta2p1_WPTight_Gsf_HT150_;
  int prescale_HLT_Ele27_WPTight_Gsf_;
  int prescale_HLT_Ele32_WPTight_Gsf_;
  int prescale_HLT_Ele32_WPTight_Gsf_L1DoubleEG_;
  int prescale_HLT_Ele35_WPTight_Gsf_;
  int prescale_HLT_Ele35_WPTight_Gsf_L1EGMT_;
  int prescale_HLT_Ele38_WPTight_Gsf_;
  int prescale_HLT_Ele40_WPTight_Gsf_;
  int prescale_HLT_IsoMu30_;
  int prescale_HLT_IsoMu27_;
  int prescale_HLT_IsoMu24_eta2p1_;
  int prescale_HLT_IsoMu24_;
  int prescale_HLT_PFMET110_PFMHT110_IDTight_;
  int prescale_HLT_PFMET120_PFMHT120_IDTight_;
  int prescale_HLT_PFMET130_PFMHT130_IDTight_;
  int prescale_HLT_PFMET140_PFMHT140_IDTight_;
  int prescale_HLT_PFHT180_;
  int prescale_HLT_PFHT250_;
  int prescale_HLT_PFHT370_;
  int prescale_HLT_PFHT430_;
  int prescale_HLT_PFHT510_;
  int prescale_HLT_PFHT590_;
  int prescale_HLT_PFHT680_;
  int prescale_HLT_PFHT780_;
  int prescale_HLT_PFHT890_;
  int prescale_HLT_PFHT1050_;
  int prescale_HLT_PFJet40_;
  int prescale_HLT_PFJet60_;
  int prescale_HLT_PFJet80_;
  int prescale_HLT_PFJet140_;
  int prescale_HLT_PFJet200_;
  int prescale_HLT_PFJet260_;
  int prescale_HLT_PFJet320_;
  int prescale_HLT_PFJet400_;
  int prescale_HLT_PFJet450_;
  int prescale_HLT_PFJet500_;
  int prescale_HLT_PFJet550_;

  // HLT Objects

  std::vector<double> pt_trigger_object_;
  std::vector<double> eta_trigger_object_;
  std::vector<double> phi_trigger_object_;

  std::vector<std::vector<string>> filter_trigger_object_;

  // Kinematics

  vdouble jet_pt_;
  vdouble jet_eta_;
  vdouble jet_phi_;
  vdouble jet_energy_;
  vdouble jet_csv_;
  vint jet_is_ID_Tight_;

  vdouble ele_pt_;
  vdouble ele_eta_;
  vdouble ele_sceta_;
  vdouble ele_phi_;
  vdouble ele_energy_;
  vdouble ele_iso_;
  vint ele_is_ID_Tight_;
  vint ele_is_ID_Loose_;

  vdouble mu_pt_;
  vdouble mu_eta_;
  vdouble mu_phi_;
  vdouble mu_energy_;
  vdouble mu_iso_;
  vint mu_is_ID_Tight_;

  double pfMET_pt_;
  double pfMET_phi_;

};


typedef std::vector<triggerStudyEventVars> vtriggerStudyEventVars;

#endif
