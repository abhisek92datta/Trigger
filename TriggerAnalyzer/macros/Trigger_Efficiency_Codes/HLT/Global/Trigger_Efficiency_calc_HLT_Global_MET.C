#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH3.h"
#include "TH2F.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TAxis.h"
#include "TMath.h"
#include "TRandom3.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <exception>
#include <cmath> 
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TGraphAsymmErrors.h"
#include "TVector.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "Math/Interpolator.h"
#include "TEfficiency.h"
#include <time.h>

#ifdef __MAKECINT__
#pragma link C++ class std::vector< TLorentzVector >+; 
#endif

#if !defined(__CINT__) && !defined(__MAKECINT__)

#include "TriggerStudyEventVars.h"

#endif

typedef std::vector< TLorentzVector >          vecTLorentzVector;
typedef std::vector<int>                       vint;
typedef std::vector<double>                    vdouble;
typedef std::vector<std::vector<double> >      vvdouble;
typedef std::vector<std::vector<int> >         vvint;

double dR_calc( double eta1, double phi1, double eta2, double phi2){

	double deta = eta1 - eta2;
	double dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
	double dR = sqrt(deta*deta + dphi*dphi);
	return dR;
}


int reco_hlt_dR_match( double reco_eta, double reco_phi, std::vector<double> hlt_eta, std::vector<double> hlt_phi, double dR_cut ){

	double dR_min = 9999;
	int match_index = -1;

	for(unsigned int i=0; i<hlt_eta.size(); i++){

		double deta = reco_eta - hlt_eta[i];
		double dphi = TVector2::Phi_mpi_pi(reco_phi - hlt_phi[i]);
		double dR = sqrt(deta*deta + dphi*dphi);

		if(dR < dR_cut){
			if(dR < dR_min){
				match_index = i;
				dR_min = dR;
			}
		}

	}
	std::cout<<dR_min<<"  "<<match_index<<"\n";
	return match_index;
}


int check_filter_label( std::vector<string> filter_label_list , string filter ){

	int filter_match = 0;
	for(unsigned int i=0; i<filter_label_list.size(); i++){
		if( filter_label_list[i].compare(filter) == 0 ){
			filter_match = 1;
			break;
		}
	}
	return filter_match;
}

void Trigger_Efficiency_calc_HLT_Global_MET( int maxNentries=-1, int Njobs=1, int jobN=1 ) {

  clock_t t;
  t = clock();

  ifstream fin;
  int n_input_files;
  string file;
  std::vector<std::string> treefilenames;
  std::string output_filename;

  fin.open("input.txt");
  fin>>n_input_files;
  for(int i=0; i<n_input_files; i++){
	fin>>file;
	treefilenames.push_back(file);
  }
  fin>>output_filename;
  fin.close();

  int N_all_total = 0;
  int N_total = 0;
  int N_e = 0;
  int N_e_den1 = 0;
  int N_e_den2 = 0;
  int N_e_den11 = 0;
  int N_e_den22 = 0;
  int N_e_den33 = 0;
  int N_e_den44 = 0;
  int N_e_ele27 = 0;
  int N_e_ele32 = 0;
  int N_e_ele32_l1eg = 0;
  int N_e_ele35 = 0;
  int N_e_ele30_jet35 = 0;
  int N_e_ele28_ht150 = 0;
  double effi_e_ele27 = 0;
  double effi_e_ele35 = 0;
  double effi_e_ele32 = 0;
  double effi_e_ele32_l1eg = 0;
  double effi_e_ele30_jet35 = 0;
  double effi_e_ele28_ht150 = 0;
  
  std::string str_jobN;
  std::stringstream stream;
  stream << jobN;
  str_jobN = stream.str();
  
  std::string histofilename = output_filename;

  TChain *chain = new TChain("trigger/triggerTree");
 
  for(int i=0; i<n_input_files; i++)
	chain->Add(treefilenames[i].c_str());

  triggerStudyEventVars *eve=0;
  chain->SetBranchAddress("eve", &eve );

  TFile histofile(histofilename.c_str(),"recreate");
  histofile.cd();

  //Histograms
	 
  TH1::SetDefaultSumw2();

	int n_ptBins;
	double x_ptBins[200];
	int s=0;
	int n_bin=0;
	while(s<100){
		x_ptBins[n_bin++]=s;
		s=s+10;
	}
	while(s<200){
		x_ptBins[n_bin++]=s;
		s=s+25;
	}
	while(s<=400){
		x_ptBins[n_bin++]=s;
		s=s+50;
	}
	n_ptBins = n_bin-1;

	s=n_bin=0;
	int n_htBins;
	double x_htBins[200];
	while(s<300){
		x_htBins[n_bin++]=s;
		s=s+25;
	}
	while(s<600){
		x_htBins[n_bin++]=s;
		s=s+50;
	}
	while(s<=1000){
		x_htBins[n_bin++]=s;
		s=s+100;
	}
	n_htBins = n_bin-1;

	TEfficiency* Eff_ele_pt_ele27 = new TEfficiency("Eff_ele_pt_ele27","Efficiency vs Electron pT for Ele27;Electron pT (GeV);Efficiency",n_ptBins,x_ptBins);
	TEfficiency* Eff_ele_eta_ele27 = new TEfficiency("Eff_ele_eta_ele27","Efficiency vs Electron eta for Ele27;Electron eta;Efficiency",10,-3,3);
	TEfficiency* Eff_njets_ele27 = new TEfficiency("Eff_njets_ele27","Efficiency vs Nr. of Jets for Ele27;Nr. of Jets;Efficiency",15,0,15);
    TEfficiency* Eff_ele_pt_ele32 = new TEfficiency("Eff_ele_pt_ele32","Efficiency vs Electron pT for Ele32;Electron pT (GeV);Efficiency",n_ptBins,x_ptBins);
    TEfficiency* Eff_ele_eta_ele32 = new TEfficiency("Eff_ele_eta_ele32","Efficiency vs Electron eta for Ele32;Electron eta;Efficiency",10,-3,3);
    TEfficiency* Eff_njets_ele32 = new TEfficiency("Eff_njets_ele32","Efficiency vs Nr. of Jets for Ele32;Nr. of Jets;Efficiency",15,0,15);
    TEfficiency* Eff_ele_pt_ele32_l1eg = new TEfficiency("Eff_ele_pt_ele32_l1eg","Efficiency vs Electron pT for Ele32_l1eg;Electron pT (GeV);Efficiency",n_ptBins,x_ptBins);
    TEfficiency* Eff_ele_eta_ele32_l1eg = new TEfficiency("Eff_ele_eta_ele32_l1eg","Efficiency vs Electron eta for Ele32_l1eg;Electron eta;Efficiency",10,-3,3);
    TEfficiency* Eff_njets_ele32_l1eg = new TEfficiency("Eff_njets_ele32_l1eg","Efficiency vs Nr. of Jets for Ele32_l1eg;Nr. of Jets;Efficiency",15,0,15);
	TEfficiency* Eff_ele_pt_ele35 = new TEfficiency("Eff_ele_pt_ele35","Efficiency vs Electron pT for Ele35;Electron pT (GeV);Efficiency",n_ptBins,x_ptBins);
	TEfficiency* Eff_ele_eta_ele35 = new TEfficiency("Eff_ele_eta_ele35","Efficiency vs Electron eta for Ele35;Electron eta;Efficiency",10,-3,3);
	TEfficiency* Eff_njets_ele35 = new TEfficiency("Eff_njets_ele35","Efficiency vs Nr. of Jets for Ele35;Nr. of Jets;Efficiency",15,0,15);
	TEfficiency* Eff_ele_pt_ele30_jet35 = new TEfficiency("Eff_ele_pt_ele30_jet35","Efficiency vs Electron pT for Ele30_Jet35;Electron pT (GeV);Efficiency",n_ptBins,x_ptBins);
	TEfficiency* Eff_ele_eta_ele30_jet35 = new TEfficiency("Eff_ele_eta_ele30_jet35","Efficiency vs Electron eta for Ele30_Jet35;Electron eta;Efficiency",10,-3,3);
	TEfficiency* Eff_jet_pt_ele30_jet35 = new TEfficiency("Eff_jet_pt_ele30_jet35","Efficiency vs Leading Jet pT for Ele30_Jet35;Leading Jet pT (GeV);Efficiency",n_ptBins,x_ptBins);
	TEfficiency* Eff_jet_eta_ele30_jet35 = new TEfficiency("Eff_jet_eta_ele30_jet35","Efficiency vs Leading Jet eta for Ele30_Jet35;Leading Jet eta;Efficiency",10,-3,3);
	TEfficiency* Eff_ht_ele30_jet35 = new TEfficiency("Eff_ht_ele30_jet35","Efficiency vs HT for Ele30_Jet35;HT (GeV);Efficiency",n_htBins,x_htBins);
	TEfficiency* Eff_njets_ele30_jet35 = new TEfficiency("Eff_njets_ele30_jet35","Efficiency vs Nr. of Jets for Ele30_Jet35;Nr. of Jets;Efficiency",15,0,15);
	TEfficiency* Eff_ele_pt_ele28_ht150 = new TEfficiency("Eff_ele_pt_ele28_ht150","Efficiency vs Electron pT for Ele28_HT150;Electron pT (GeV);Efficiency",n_ptBins,x_ptBins);
	TEfficiency* Eff_ele_eta_ele28_ht150 = new TEfficiency("Eff_ele_eta_ele28_ht150","Efficiency vs Electron eta for Ele28_HT150;Electron eta;Efficiency",10,-3,3);
	TEfficiency* Eff_jet_pt_ele28_ht150 = new TEfficiency("Eff_jet_pt_ele28_ht150","Efficiency vs Leading Jet pT for Ele28_HT150;Leading Jet pT (GeV);Efficiency",n_ptBins,x_ptBins);
	TEfficiency* Eff_jet_eta_ele28_ht150 = new TEfficiency("Eff_jet_eta_ele28_ht150","Efficiency vs Leading Jet eta for Ele28_HT150;Leading Jet eta;Efficiency",10,-3,3);
	TEfficiency* Eff_ht_ele28_ht150 = new TEfficiency("Eff_ht_ele28_ht150","Efficiency vs HT for Ele28_HT150;HT (GeV);Efficiency",n_htBins,x_htBins);
	TEfficiency* Eff_njets_ele28_ht150 = new TEfficiency("Eff_njets_ele28_ht150","Efficiency vs Nr. of Jets for Ele28_HT150;Nr. of Jets;Efficiency",15,0,15);
	
  //Event Loop

  int nentries = chain->GetEntries();
  std::cout << "\n\t Number of entries = " << nentries << std::endl;
  std::cout << "\t Max number of entries = " << maxNentries << std::endl;
  std::cout << "\n" << std::endl;

  int use_nentries = std::max( maxNentries, nentries);

  int NeventsPerJob = int( double(use_nentries)/double(Njobs) + 0.000001 ) + 1;

  int firstEvent = (jobN-1)*NeventsPerJob + 1;
  int lastEvent  = firstEvent + NeventsPerJob;
  if( jobN==Njobs ) lastEvent = -1;
  if( jobN==1 ) firstEvent = 0;


  std::cout << "========  Starting Event Loop  ========" << std::endl;
  for (Long64_t ievt=0; ievt<chain->GetEntries();ievt++) {    //Long64_t

        N_all_total++;

	    if( ievt<firstEvent ) continue;
	    if( ievt==lastEvent ) break;

      	if( ievt==0 )        std::cout << "     Event " << ievt+1 << std::endl;
      	if( ievt%10000==0 && ievt!=1) std::cout << "           " << ievt << "\t" 
  					     << int(double(ievt-1)/double(nentries)*100) << "% done" << std::endl;
	  
	 	if( ievt==(maxNentries+1) && ievt!=0 ) break;

	  
	  	chain->GetEntry(ievt);

        int run_nr = eve->run_nr_;
        int event_nr = eve->event_nr_;
        int lumi_nr = eve->lumi_nr_;

        N_total++;

        double rho = eve->rho_;
		
		//Grab Specific Lepton information from trees

		vdouble ele_pt = eve->ele_pt_;
        vdouble ele_eta = eve->ele_eta_;
        vdouble ele_sceta = eve->ele_sceta_;
        vdouble ele_phi = eve->ele_phi_;
        vdouble ele_energy = eve->ele_energy_;
        vdouble ele_iso = eve->ele_iso_;
        vint ele_is_ID_Tight = eve->ele_is_ID_Tight_;
        vint ele_is_ID_Loose = eve->ele_is_ID_Loose_;

        vdouble mu_pt = eve->mu_pt_;
        vdouble mu_eta = eve->mu_eta_;
        vdouble mu_phi = eve->mu_phi_;
        vdouble mu_energy = eve->mu_energy_;
        vdouble mu_iso = eve->mu_iso_;
        vint mu_is_ID_Tight = eve->mu_is_ID_Tight_;

        vdouble jet_pt = eve->jet_pt_;
        vdouble jet_eta = eve->jet_eta_;
        vdouble jet_phi = eve->jet_phi_;
        vdouble jet_energy = eve->jet_energy_;
		vint jet_is_ID_Tight = eve->jet_is_ID_Tight_;

        double met_pt = eve->pfMET_pt_;
        double met_phi = eve->pfMET_phi_;

	    int pass_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned = eve->pass_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_;
	    int pass_HLT_Ele28_eta2p1_WPTight_Gsf_HT150 = eve->pass_HLT_Ele28_eta2p1_WPTight_Gsf_HT150_;
	    int pass_HLT_Ele27_WPTight_Gsf = eve->pass_HLT_Ele27_WPTight_Gsf_;
		int pass_HLT_Ele32_WPTight_Gsf = eve->pass_HLT_Ele32_WPTight_Gsf_;
        int pass_HLT_Ele32_WPTight_Gsf_L1DoubleEG = eve->pass_HLT_Ele32_WPTight_Gsf_L1DoubleEG_;
		int pass_HLT_Ele35_WPTight_Gsf = eve->pass_HLT_Ele35_WPTight_Gsf_;
		int pass_HLT_Ele35_WPTight_Gsf_L1EGMT = eve->pass_HLT_Ele35_WPTight_Gsf_L1EGMT_;
	    int pass_HLT_Ele38_WPTight_Gsf = eve->pass_HLT_Ele38_WPTight_Gsf_;
	    int pass_HLT_Ele40_WPTight_Gsf = eve->pass_HLT_Ele40_WPTight_Gsf_;
	    int pass_HLT_IsoMu30 = eve->pass_HLT_IsoMu30_;
	    int pass_HLT_IsoMu27 = eve->pass_HLT_IsoMu27_;
	    int pass_HLT_IsoMu24_eta2p1 = eve->pass_HLT_IsoMu24_eta2p1_;
	    int pass_HLT_IsoMu24 = eve->pass_HLT_IsoMu24_;
	    int pass_HLT_PFMET110_PFMHT110_IDTight = eve->pass_HLT_PFMET110_PFMHT110_IDTight_;
	    int pass_HLT_PFMET120_PFMHT120_IDTight = eve->pass_HLT_PFMET120_PFMHT120_IDTight_;
	    int pass_HLT_PFMET130_PFMHT130_IDTight = eve->pass_HLT_PFMET130_PFMHT130_IDTight_;
	    int pass_HLT_PFMET140_PFMHT140_IDTight = eve->pass_HLT_PFMET140_PFMHT140_IDTight_;
		int pass_HLT_PFHT180 = eve->pass_HLT_PFHT180_;
	    int pass_HLT_PFHT250 = eve->pass_HLT_PFHT250_;
		int pass_HLT_PFHT370 = eve->pass_HLT_PFHT370_;
	    int pass_HLT_PFHT430 = eve->pass_HLT_PFHT430_;
	    int pass_HLT_PFHT510 = eve->pass_HLT_PFHT510_;
	    int pass_HLT_PFHT590 = eve->pass_HLT_PFHT590_;
	    int pass_HLT_PFHT680 = eve->pass_HLT_PFHT680_;
	    int pass_HLT_PFHT780 = eve->pass_HLT_PFHT780_;
	    int pass_HLT_PFHT890 = eve->pass_HLT_PFHT890_;
	    int pass_HLT_PFHT1050 = eve->pass_HLT_PFHT1050_;
		int pass_HLT_PFJet40 = eve->pass_HLT_PFJet40_;
		int pass_HLT_PFJet60 = eve->pass_HLT_PFJet60_;
		int pass_HLT_PFJet80 = eve->pass_HLT_PFJet80_;
	    int pass_HLT_PFJet140 = eve->pass_HLT_PFJet140_;
	    int pass_HLT_PFJet200 = eve->pass_HLT_PFJet200_;
	    int pass_HLT_PFJet260 = eve->pass_HLT_PFJet260_;
	    int pass_HLT_PFJet320 = eve->pass_HLT_PFJet320_;
	    int pass_HLT_PFJet400 = eve->pass_HLT_PFJet400_;
	    int pass_HLT_PFJet450 = eve->pass_HLT_PFJet450_;
	    int pass_HLT_PFJet500 = eve->pass_HLT_PFJet500_;
		int pass_HLT_PFJet550 = eve->pass_HLT_PFJet550_;

        int prescale_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned = eve->prescale_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_;
        int prescale_HLT_Ele28_eta2p1_WPTight_Gsf_HT150 = eve->prescale_HLT_Ele28_eta2p1_WPTight_Gsf_HT150_;
        int prescale_HLT_Ele27_WPTight_Gsf = eve->prescale_HLT_Ele27_WPTight_Gsf_;
        int prescale_HLT_Ele32_WPTight_Gsf = eve->prescale_HLT_Ele32_WPTight_Gsf_;
        int prescale_HLT_Ele32_WPTight_Gsf_L1DoubleEG = eve->prescale_HLT_Ele32_WPTight_Gsf_L1DoubleEG_;
        int prescale_HLT_Ele35_WPTight_Gsf = eve->prescale_HLT_Ele35_WPTight_Gsf_;
        int prescale_HLT_Ele35_WPTight_Gsf_L1EGMT = eve->prescale_HLT_Ele35_WPTight_Gsf_L1EGMT_;
        int prescale_HLT_Ele38_WPTight_Gsf = eve->prescale_HLT_Ele38_WPTight_Gsf_;
        int prescale_HLT_Ele40_WPTight_Gsf = eve->prescale_HLT_Ele40_WPTight_Gsf_;
        int prescale_HLT_IsoMu30 = eve->prescale_HLT_IsoMu30_;
        int prescale_HLT_IsoMu27 = eve->prescale_HLT_IsoMu27_;
        int prescale_HLT_IsoMu24_eta2p1 = eve->prescale_HLT_IsoMu24_eta2p1_;
        int prescale_HLT_IsoMu24 = eve->prescale_HLT_IsoMu24_;
        int prescale_HLT_PFMET110_PFMHT110_IDTight = eve->prescale_HLT_PFMET110_PFMHT110_IDTight_;
        int prescale_HLT_PFMET120_PFMHT120_IDTight = eve->prescale_HLT_PFMET120_PFMHT120_IDTight_;
        int prescale_HLT_PFMET130_PFMHT130_IDTight = eve->prescale_HLT_PFMET130_PFMHT130_IDTight_;
        int prescale_HLT_PFMET140_PFMHT140_IDTight = eve->prescale_HLT_PFMET140_PFMHT140_IDTight_;
        int prescale_HLT_PFHT180 = eve->prescale_HLT_PFHT180_;
        int prescale_HLT_PFHT250 = eve->prescale_HLT_PFHT250_;
        int prescale_HLT_PFHT370 = eve->prescale_HLT_PFHT370_;
        int prescale_HLT_PFHT430 = eve->prescale_HLT_PFHT430_;
        int prescale_HLT_PFHT510 = eve->prescale_HLT_PFHT510_;
        int prescale_HLT_PFHT590 = eve->prescale_HLT_PFHT590_;
        int prescale_HLT_PFHT680 = eve->prescale_HLT_PFHT680_;
        int prescale_HLT_PFHT780 = eve->prescale_HLT_PFHT780_;
        int prescale_HLT_PFHT890 = eve->prescale_HLT_PFHT890_;
        int prescale_HLT_PFHT1050 = eve->prescale_HLT_PFHT1050_;
        int prescale_HLT_PFJet40 = eve->prescale_HLT_PFJet40_;
        int prescale_HLT_PFJet60 = eve->prescale_HLT_PFJet60_;
        int prescale_HLT_PFJet80 = eve->prescale_HLT_PFJet80_;
        int prescale_HLT_PFJet140 = eve->prescale_HLT_PFJet140_;
        int prescale_HLT_PFJet200 = eve->prescale_HLT_PFJet200_;
        int prescale_HLT_PFJet260 = eve->prescale_HLT_PFJet260_;
        int prescale_HLT_PFJet320 = eve->prescale_HLT_PFJet320_;
        int prescale_HLT_PFJet400 = eve->prescale_HLT_PFJet400_;
        int prescale_HLT_PFJet450 = eve->prescale_HLT_PFJet450_;
        int prescale_HLT_PFJet500 = eve->prescale_HLT_PFJet500_;
        int prescale_HLT_PFJet550 = eve->prescale_HLT_PFJet550_;

		std::vector<double> pt_trigger_object = eve->pt_trigger_object_;
		std::vector<double> eta_trigger_object = eve->eta_trigger_object_;
		std::vector<double> phi_trigger_object = eve->phi_trigger_object_;
		std::vector<std::vector<string>> filter_trigger_object = eve->filter_trigger_object_;

        int nele_tight = 0;
        int nele_loose = 0;
        int nmu_tight = 0;
        int nmu_loose = 0;
        double electron_pt = -9999;
        double electron_eta = -9999;
        double electron_phi = -9999;
        double muon_pt = -9999;

		double HT = 0;
        int njets = 0;
        double leading_jet_pt = -9999;
		double leading_jet_eta = -9999;
		double leading_jet_phi = -9999;

		// Reference Trigger pass selection
	    //if( pass_HLT_PFMET110_PFMHT110_IDTight != 1 && pass_HLT_PFMET120_PFMHT120_IDTight != 1 && pass_HLT_PFMET130_PFMHT130_IDTight != 1 && pass_HLT_PFMET140_PFMHT140_IDTight != 1 )
		if( pass_HLT_PFMET110_PFMHT110_IDTight != 1 && pass_HLT_PFMET120_PFMHT120_IDTight != 1 )
		  continue;

		// Selecting Electrons
        for(unsigned int i=0; i<ele_pt.size(); i++) {
            double absEta = fabs(ele_sceta[i]);
            bool inCrack = absEta >= 1.4442 && absEta <= 1.5660;
            if(ele_pt[i] > 15 && fabs(ele_eta[i]) < 2.4 && !inCrack && ele_is_ID_Tight[i] == 1) {
                nele_loose++;
                if(ele_pt[i] > 30 && fabs(ele_eta[i]) < 2.1) {
                    nele_tight++;
                    if(ele_pt[i] > electron_pt){
                        electron_pt = ele_pt[i];
                        electron_eta = ele_eta[i];
                        electron_phi = ele_phi[i];
                    }
                }
            }
        }

	  	// Selecting Muons
        for(unsigned int i=0; i<mu_pt.size(); i++) {
          if(mu_pt[i] > 15 && fabs(mu_eta[i]) < 2.4 && mu_is_ID_Tight[i] == 1 && mu_iso[i] < 0.25) {
              nmu_loose++;
              if(mu_pt[i] > 26 && fabs(mu_eta[i]) < 2.1 && mu_iso[i] < 0.15) {
                  nmu_tight++;
                  if(mu_pt[i] > muon_pt)
                      muon_pt = mu_pt[i];
              }
          }
        }

	    // Selecting Jets
    	for(unsigned int i=0; i<jet_pt.size(); i++) {
           if(jet_is_ID_Tight[i]!=1)
				continue;
       	   if(jet_pt[i] > 30) {
         	     if(fabs(jet_eta[i]) < 2.4){
                    double dR = dR_calc(jet_eta[i], jet_phi[i], electron_eta, electron_phi);
          	        if (dR > 0.4){
           	           njets++;
           	           HT = HT + jet_pt[i];
						  if(jet_pt[i] > leading_jet_pt){
            	              leading_jet_pt = jet_pt[i];
							  leading_jet_eta = jet_eta[i];
							  leading_jet_phi = jet_phi[i];
						  }
				    }
              	}
		   }
        }

        // Single Electron event selection (requires MET > 100 GeV)
        if(nele_tight == 1){
            if( nele_loose == 1 && met_pt > 100){
                if(njets >= 1){
                    N_e++;

					// Matching of RECO objects with HLT Objects for signal trigger and reference triggers

					double dR_min1;
					double dR_min2;
					double dR_min11;
					double dR_min22;
                    double dR_min33;
                    double dR_min44;
					double dR_min3;
					double dR_min4;
					double dR;

					dR_min1 = 9999;
					dR_min2 = 9999;
					dR_min11 = 9999;
					dR_min22 = 9999;
                    dR_min33 = 9999;
                    dR_min44 = 9999;
					dR_min3 = 9999;
					dR_min4 = 9999;

					int ele_label1_match = 0;
					int ele_label2_match = 0;
					int ele_label11_match = 0;
					int ele_label22_match = 0;
                    int ele_label33_match = 0;
                    int ele_label44_match = 0;
					int ele_dR1_match_index = 0;
					int ele_dR2_match_index = 0;
					int ele_dR11_match_index = 0;
					int ele_dR22_match_index = 0;
                    int ele_dR33_match_index = 0;
                    int ele_dR44_match_index = 0;
					int ele_sig1_hlt_match = 0;
					int ele_sig2_hlt_match = 0;
					int ele_sig11_hlt_match = 0;
					int ele_sig22_hlt_match = 0;
                    int ele_sig33_hlt_match = 0;
                    int ele_sig44_hlt_match = 0;

					int jet_label1_match = 0;
					int jet_label2_match = 0;
					int jet_dR1_match_index = 0;
					int jet_dR2_match_index = 0;
					int jet_sig1_hlt_match = 0;
					int jet_sig2_hlt_match = 0;

					for(unsigned int l=0; l<eta_trigger_object.size(); l++){

						// Matching Electron

						ele_label1_match = check_filter_label( filter_trigger_object[l] , "hltEle30erJetC34WPTightGsfTrackIsoFilter" );
						if (ele_label1_match == 1){
							dR = dR_calc(electron_eta, electron_phi, eta_trigger_object[l], phi_trigger_object[l]);
							if( dR < 0.3 ){
								if( dR < dR_min1){
									ele_dR1_match_index = l;
									ele_sig1_hlt_match = 1;
									dR_min1 = dR;
								}
							}
						}

						ele_label2_match = check_filter_label( filter_trigger_object[l] , "hltEle28erHTT100WPTightGsfTrackIsoFilter" );
						if (ele_label2_match == 1){
							dR = dR_calc(electron_eta, electron_phi, eta_trigger_object[l], phi_trigger_object[l]);
							if( dR < 0.3 ){
								if( dR < dR_min2){
									ele_dR2_match_index = l;
									ele_sig2_hlt_match = 1;
									dR_min2 = dR;
								}
							}
						}

						ele_label11_match = check_filter_label( filter_trigger_object[l] , "hltEle27WPTightGsfTrackIsoFilter" );
						if (ele_label11_match == 1){
						 	dR = dR_calc(electron_eta, electron_phi, eta_trigger_object[l], phi_trigger_object[l]);
							if( dR < 0.3 ){
						 		if( dR < dR_min11){
						 			ele_dR11_match_index = l;
						 			ele_sig11_hlt_match = 1;
						 			dR_min11 = dR;
						 		}
						 	}
						}

                        ele_label33_match = check_filter_label( filter_trigger_object[l] , "hltEle32WPTightGsfTrackIsoFilter" );
                        if (ele_label33_match == 1){
                            dR = dR_calc(electron_eta, electron_phi, eta_trigger_object[l], phi_trigger_object[l]);
                            if( dR < 0.3 ){
                                if( dR < dR_min33){
                                    ele_dR33_match_index = l;
                                    ele_sig33_hlt_match = 1;
                                    dR_min33 = dR;
                                }
                            }
                        }

                        ele_label44_match = check_filter_label( filter_trigger_object[l] , "hltEle32L1DoubleEGWPTightGsfTrackIsoFilter" );
                        if (ele_label44_match == 1){
                            dR = dR_calc(electron_eta, electron_phi, eta_trigger_object[l], phi_trigger_object[l]);
                            if( dR < 0.3 ){
                                if( dR < dR_min44){
                                    ele_dR44_match_index = l;
                                    ele_sig44_hlt_match = 1;
                                    dR_min44 = dR;
                                }
                            }
                        }

						ele_label22_match = check_filter_label( filter_trigger_object[l] , "hltEle35noerWPTightGsfTrackIsoFilter" );
						if (ele_label22_match == 1){
							dR = dR_calc(electron_eta, electron_phi, eta_trigger_object[l], phi_trigger_object[l]);
							if( dR < 0.3 ){
								if( dR < dR_min22){
									ele_dR22_match_index = l;
									ele_sig22_hlt_match = 1;
									dR_min22 = dR;
								}
							}
						}

						// Matching Leading Jet

						jet_label1_match = check_filter_label( filter_trigger_object[l] , "hltEle30PFJet35EleCleaned" );
						if (jet_label1_match == 1){
							dR = dR_calc(leading_jet_eta, leading_jet_phi, eta_trigger_object[l], phi_trigger_object[l]);
							if( dR < 0.3 ){
								if( dR < dR_min3){
									jet_dR1_match_index = l;
									jet_sig1_hlt_match = 1;
									dR_min3 = dR;
								}
							}
						}
					}

					jet_label2_match = 1;
					jet_dR2_match_index = -1;
					jet_sig2_hlt_match = 1;

					int no_overlap = 0;

					if(pass_HLT_Ele27_WPTight_Gsf != -1 && prescale_HLT_Ele27_WPTight_Gsf != 0){
						N_e_den11++;
						no_overlap = 1;

						int pass_ele27 = pass_HLT_Ele27_WPTight_Gsf * ele_sig11_hlt_match * no_overlap;

						if(pass_ele27 == 1)
							N_e_ele27++;
						Eff_ele_pt_ele27->Fill(pass_ele27,electron_pt);
						Eff_ele_eta_ele27->Fill(pass_ele27,electron_eta);
						Eff_njets_ele27->Fill(pass_ele27,njets);
					}

                    if(pass_HLT_Ele32_WPTight_Gsf != -1 && prescale_HLT_Ele32_WPTight_Gsf != 0){
                        N_e_den33++;
                        no_overlap = 1;

                        int pass_ele32 = pass_HLT_Ele32_WPTight_Gsf * ele_sig33_hlt_match * no_overlap;

                        if(pass_ele32 == 1)
                            N_e_ele32++;
                        Eff_ele_pt_ele32->Fill(pass_ele32,electron_pt);
                        Eff_ele_eta_ele32->Fill(pass_ele32,electron_eta);
                        Eff_njets_ele32->Fill(pass_ele32,njets);
                    }

                    if(pass_HLT_Ele32_WPTight_Gsf_L1DoubleEG != -1 && prescale_HLT_Ele32_WPTight_Gsf_L1DoubleEG != 0){
                        N_e_den44++;
                        no_overlap = 1;

                        int pass_ele32_l1eg = pass_HLT_Ele32_WPTight_Gsf_L1DoubleEG * ele_sig44_hlt_match * no_overlap;

                        if(pass_ele32_l1eg == 1)
                            N_e_ele32_l1eg++;
                        Eff_ele_pt_ele32_l1eg->Fill(pass_ele32_l1eg,electron_pt);
                        Eff_ele_eta_ele32_l1eg->Fill(pass_ele32_l1eg,electron_eta);
                        Eff_njets_ele32_l1eg->Fill(pass_ele32_l1eg,njets);
                    }

					if(pass_HLT_Ele35_WPTight_Gsf != -1 && prescale_HLT_Ele35_WPTight_Gsf != 0){
						N_e_den22++;
						no_overlap = 1;

						int pass_ele35 = pass_HLT_Ele35_WPTight_Gsf * ele_sig22_hlt_match * no_overlap;

						if(pass_ele35 == 1)
							N_e_ele35++;
						Eff_ele_pt_ele35->Fill(pass_ele35,electron_pt);
						Eff_ele_eta_ele35->Fill(pass_ele35,electron_eta);
						Eff_njets_ele35->Fill(pass_ele35,njets);
					}

					if(pass_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned != -1 && prescale_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned != 0){

						N_e_den1++;
						/*
						if ( dR_calc( leading_jet_eta, leading_jet_phi, eta_trigger_object[ele_dR1_match_index], phi_trigger_object[ele_dR1_match_index]  ) > 0.3 ) {
							if ( dR_calc( electron_eta, electron_phi, eta_trigger_object[jet_dR1_match_index], phi_trigger_object[jet_dR1_match_index]  ) > 0.3 )
								no_overlap = 1;
						}
						*/
						no_overlap = 1;

						int pass_ele30_jet35 = pass_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned * ele_sig1_hlt_match * jet_sig1_hlt_match * no_overlap;

						if(pass_ele30_jet35 == 1)
							N_e_ele30_jet35++;
						Eff_ele_pt_ele30_jet35->Fill(pass_ele30_jet35,electron_pt);
						Eff_ele_eta_ele30_jet35->Fill(pass_ele30_jet35,electron_eta);
						Eff_jet_pt_ele30_jet35->Fill(pass_ele30_jet35,leading_jet_pt);
						Eff_jet_eta_ele30_jet35->Fill(pass_ele30_jet35,leading_jet_eta);
						Eff_ht_ele30_jet35->Fill(pass_ele30_jet35,HT);
						Eff_njets_ele30_jet35->Fill(pass_ele30_jet35,njets);
					}

					if(pass_HLT_Ele28_eta2p1_WPTight_Gsf_HT150 != -1 && prescale_HLT_Ele28_eta2p1_WPTight_Gsf_HT150 != 0){
						N_e_den2++;

						int pass_ele28_ht150 = pass_HLT_Ele28_eta2p1_WPTight_Gsf_HT150 * ele_sig2_hlt_match * jet_sig2_hlt_match;
						if(pass_ele28_ht150 == 1)
							N_e_ele28_ht150++;
						Eff_ele_pt_ele28_ht150->Fill(pass_ele28_ht150,electron_pt);
						Eff_ele_eta_ele28_ht150->Fill(pass_ele28_ht150,electron_eta);
						Eff_jet_pt_ele28_ht150->Fill(pass_ele28_ht150,leading_jet_pt);
						Eff_jet_eta_ele28_ht150->Fill(pass_ele28_ht150,leading_jet_eta);
						Eff_ht_ele28_ht150->Fill(pass_ele28_ht150,HT);
						Eff_njets_ele28_ht150->Fill(pass_ele28_ht150,njets);
					}
                }
            }
        }
   }

	effi_e_ele27 = double(N_e_ele27) / N_e_den11;
    effi_e_ele32 = double(N_e_ele32) / N_e_den33;
    effi_e_ele32_l1eg = double(N_e_ele32_l1eg) / N_e_den44;
    effi_e_ele35 = double(N_e_ele35) / N_e_den22;
	effi_e_ele30_jet35 = double(N_e_ele30_jet35) / N_e_den1;
	effi_e_ele28_ht150 = double(N_e_ele28_ht150) / N_e_den2;

	t=clock()-t;

	std::cout << " Done! " <<((float)t)/CLOCKS_PER_SEC<< std::endl;
	std::cout<<"**********************************************************************************************\n";
	std::cout<<"Total No. of events in all runs : "<<N_all_total<<"\n";
	std::cout<<"Total No. of events in selected runs : "<<N_total<<"\n";
	std::cout<<"No. of events passing Single Electron event selection + reference triggers only : "<<N_e<<"\n";
	std::cout<<"Denominator for Ele27 HLT Path : "<<N_e_den11<<"\n";
	std::cout<<"No. of events and efficiency passing Single Electron event selection + reference triggers + Ele27 HLT Path : "<<N_e_ele27<<" ,  efficiency = "<<effi_e_ele27<<"\n";
    std::cout<<"Denominator for Ele32 HLT Path : "<<N_e_den33<<"\n";
    std::cout<<"No. of events and efficiency passing Single Electron event selection + reference triggers + Ele32 HLT Path : "<<N_e_ele32<<" ,  efficiency = "<<effi_e_ele32<<"\n";
    std::cout<<"Denominator for Ele32_L1DoubleEG HLT Path : "<<N_e_den44<<"\n";
    std::cout<<"No. of events and efficiency passing Single Electron event selection + reference triggers + Ele32_L1DoubleEG HLT Path : "<<N_e_ele32_l1eg<<" ,  efficiency = "<<effi_e_ele32_l1eg<<"\n";
	std::cout<<"Denominator for Ele35 HLT Path : "<<N_e_den22<<"\n";
	std::cout<<"No. of events and efficiency passing Single Electron event selection + reference triggers + Ele35 HLT Path : "<<N_e_ele35<<" ,  efficiency = "<<effi_e_ele35<<"\n";
	std::cout<<"Denominator for Ele30_Jet35 HLT Path : "<<N_e_den1<<"\n";
	std::cout<<"No. of events and efficiency passing Single Electron event selection + reference triggers + Ele30_Jet35 HLT Path : "<<N_e_ele30_jet35<<" ,  efficiency = "<<effi_e_ele30_jet35<<"\n";
	std::cout<<"Denominator for Ele28_HTT150 HLT Path : "<<N_e_den2<<"\n";
	std::cout<<"No. of events and efficiency passing Single Electron event selection + reference triggers + Ele28_HTT150 HLT Path : "<<N_e_ele28_ht150<<" , efficiency = "<<effi_e_ele28_ht150<<"\n";
	std::cout<<"**********************************************************************************************\n";
	histofile.Write();
	histofile.Close();

}


