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

void Trigger_Efficiency_calc_L1_Global_MET( int maxNentries=-1, int Njobs=1, int jobN=1 ) {

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
  int N_e_EG_all = 0;
  int N_e_ele30_jet34 = 0;
  int N_e_ele28_ht100 = 0;
  int N_e_ele30_jet34_OR = 0;
  int N_e_ele28_ht100_OR = 0;
  int N_e_ele30_jet34_OR_prescaled = 0;
  int N_e_ele28_ht100_OR_prescaled = 0;
  double effi_e_EG_all = 0;
  double effi_e_ele30_jet34 = 0;
  double effi_e_ele28_ht100 = 0;
  double effi_e_ele30_jet34_OR = 0;
  double effi_e_ele28_ht100_OR = 0;
  double effi_e_ele30_jet34_OR_prescaled = 0;
  double effi_e_ele28_ht100_OR_prescaled = 0;
  
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

  TEfficiency* Eff_ele_pt_EG_all = new TEfficiency("Eff_ele_pt_EG_all","Efficiency vs Electron pT for OR of all single EG;Electron pT (GeV);Efficiency", n_ptBins,x_ptBins);
  TEfficiency* Eff_ele_eta_EG_all = new TEfficiency("Eff_ele_eta_EG_all","Efficiency vs Electron Eta for OR of all single EG;Electron Eta;Efficiency", 10,-3,3);
  TEfficiency* Eff_njets_EG_all = new TEfficiency("Eff_njets_EG_all","Efficiency vs Nr. of Jets for OR of all single EG;Nr. of Jets;Efficiency",15,0,15);
  TEfficiency* Eff_ele_pt_ele30_jet34 = new TEfficiency("Eff_ele_pt_ele30_jet34","Efficiency vs Electron pT for Ele30_Jet34;Electron pT (GeV);Efficiency", n_ptBins,x_ptBins);
  TEfficiency* Eff_ele_eta_ele30_jet34 = new TEfficiency("Eff_ele_eta_ele30_jet34","Efficiency vs Electron Eta for Ele30_Jet34;Electron Eta;Efficiency", 10,-3,3);
  TEfficiency* Eff_jet_pt_ele30_jet34 = new TEfficiency("Eff_jet_pt_ele30_jet34","Efficiency vs Leading Jet pT for Ele30_Jet34;Leading Jet pT (GeV);Efficiency", n_ptBins,x_ptBins);
  TEfficiency* Eff_jet_eta_ele30_jet34 = new TEfficiency("Eff_jet_eta_ele30_jet34","Efficiency vs Leading Jet Eta for Ele30_Jet34;Leading Jet Eta;Efficiency", 10,-3,3);
  TEfficiency* Eff_njets_ele30_jet34 = new TEfficiency("Eff_njets_ele30_jet34","Efficiency vs Nr. of Jets for Ele30_Jet34;Nr. of Jets;Efficiency",15,0,15);
  TEfficiency* Eff_ele_pt_ele28_ht100 = new TEfficiency("Eff_ele_pt_ele28_ht100","Efficiency vs Electron pT for Ele28_HT100;Electron pT (GeV);Efficiency", n_ptBins,x_ptBins);
  TEfficiency* Eff_ele_eta_ele28_ht100 = new TEfficiency("Eff_ele_eta_ele28_ht100","Efficiency vs Electron Eta for Ele28_HT100;Electron Eta;Efficiency", 10,-3,3);
  TEfficiency* Eff_ht_ele28_ht100 = new TEfficiency("Eff_ht_ele28_ht100","Efficiency vs HT for Ele28_HT100;HT (GeV);Efficiency", n_htBins,x_htBins);
  TEfficiency* Eff_njets_ele28_ht100 = new TEfficiency("Eff_njets_ele28_ht100","Efficiency vs Nr. of Jets for Ele28_HT100;Nr. of Jets;Efficiency",15,0,15);
  TEfficiency* Eff_ele_pt_ele30_jet34_OR = new TEfficiency("Eff_ele_pt_ele30_jet34_OR","Efficiency vs Electron pT for Ele30_Jet34 OR Ele40 OR IsoEle38 OR IsoEle36er;Electron pT (GeV);Efficiency", n_ptBins,x_ptBins);
  TEfficiency* Eff_ele_eta_ele30_jet34_OR = new TEfficiency("Eff_ele_eta_ele30_jet34_OR","Efficiency vs Electron Eta for Ele30_Jet34 OR Ele40 OR IsoEle38 OR IsoEle36er;Electron Eta;Efficiency", 10,-3,3);
  TEfficiency* Eff_jet_pt_ele30_jet34_OR = new TEfficiency("Eff_jet_pt_ele30_jet34_OR","Efficiency vs Leading Jet pT for Ele30_Jet34 OR Ele40 OR IsoEle38 OR IsoEle36er;Leading Jet pT (GeV);Efficiency", n_ptBins,x_ptBins);
  TEfficiency* Eff_jet_eta_ele30_jet34_OR = new TEfficiency("Eff_jet_eta_ele30_jet34_OR","Efficiency vs Leading Jet Eta for Ele30_Jet34 OR Ele40 OR IsoEle38 OR IsoEle36er;Leading Jet Eta;Efficiency", 10,-3,3);
  TEfficiency* Eff_njets_ele30_jet34_OR = new TEfficiency("Eff_njets_ele30_jet34_OR","Efficiency vs Nr. of Jets for Ele30_Jet34 OR Ele40 OR IsoEle38 OR IsoEle36er;Nr. of Jets;Efficiency",15,0,15);
  TEfficiency* Eff_ele_pt_ele28_ht100_OR = new TEfficiency("Eff_ele_pt_ele28_ht100_OR","Efficiency vs Electron pT for Ele28_HTT100 OR Ele40 OR IsoEle38 OR IsoEle36er;Electron pT (GeV);Efficiency", n_ptBins,x_ptBins);
  TEfficiency* Eff_ele_eta_ele28_ht100_OR = new TEfficiency("Eff_ele_eta_ele28_ht100_OR","Efficiency vs Electron Eta for Ele28_HTT100 OR Ele40 OR IsoEle38 OR IsoEle36er;Electron Eta;Efficiency", 10,-3,3);
  TEfficiency* Eff_ht_ele28_ht100_OR = new TEfficiency("Eff_ht_ele28_ht100_OR","Efficiency vs HT for Ele28_HTT100 OR Ele40 OR IsoEle38 OR IsoEle36er;HT (GeV);Efficiency", n_htBins,x_htBins);
  TEfficiency* Eff_njets_ele28_ht100_OR = new TEfficiency("Eff_njets_ele28_ht100_OR","Efficiency vs Nr. of Jets for Ele28_HTT100 OR Ele40 OR IsoEle38 OR IsoEle36er;Nr. of Jets;Efficiency",15,0,15);
  TEfficiency* Eff_ele_pt_ele30_jet34_OR_prescaled = new TEfficiency("Eff_ele_pt_ele30_jet34_OR_prescaled","Efficiency vs Electron pT for Ele30_Jet34 OR Ele26_Jet34 OR Ele28_Jet34;Electron pT (GeV);Efficiency", n_ptBins,x_ptBins);
  TEfficiency* Eff_ele_eta_ele30_jet34_OR_prescaled = new TEfficiency("Eff_ele_eta_ele30_jet34_OR_prescaled","Efficiency vs Electron Eta for Ele30_Jet34 OR Ele26_Jet34 OR Ele28_Jet34;Electron Eta;Efficiency", 10,-3,3);
  TEfficiency* Eff_jet_pt_ele30_jet34_OR_prescaled = new TEfficiency("Eff_jet_pt_ele30_jet34_OR_prescaled","Efficiency vs Leading Jet pT for Ele30_Jet34 OR Ele26_Jet34 OR Ele28_Jet34;Leading Jet pT (GeV);Efficiency", n_ptBins,x_ptBins);
  TEfficiency* Eff_jet_eta_ele30_jet34_OR_prescaled = new TEfficiency("Eff_jet_eta_ele30_jet34_OR_prescaled","Efficiency vs Leading Jet Eta for Ele30_Jet34 OR Ele26_Jet34 OR Ele28_Jet34;Leading Jet Eta;Efficiency", 10,-3,3);
  TEfficiency* Eff_njets_ele30_jet34_OR_prescaled = new TEfficiency("Eff_njets_ele30_jet34_OR_prescaled","Efficiency vs Nr. of Jets for Ele30_Jet34 OR Ele26_Jet34 OR Ele28_Jet34;Nr. of Jets;Efficiency",15,0,15);
  TEfficiency* Eff_ele_pt_ele28_ht100_OR_prescaled = new TEfficiency("Eff_ele_pt_ele28_ht100_OR_prescaled","Efficiency vs Electron pT for Ele28_HTT100 OR Ele24_HTT100 OR Ele26_HTT100;Electron pT (GeV);Efficiency", n_ptBins,x_ptBins);
  TEfficiency* Eff_ele_eta_ele28_ht100_OR_prescaled = new TEfficiency("Eff_ele_eta_ele28_ht100_OR_prescaled","Efficiency vs Electron Eta for Ele28_HTT100 OR Ele24_HTT100 OR Ele26_HTT100;Electron Eta;Efficiency", 10,-3,3);
  TEfficiency* Eff_ht_ele28_ht100_OR_prescaled = new TEfficiency("Eff_ht_ele28_ht100_OR_prescaled","Efficiency vs HT for Ele28_HTT100 OR Ele24_HTT100 OR Ele26_HTT100;HT (GeV);Efficiency", n_htBins,x_htBins);
  TEfficiency* Eff_njets_ele28_ht100_OR_prescaled = new TEfficiency("Eff_njets_ele28_ht100_OR_prescaled","Efficiency vs Nr. of Jets for Ele28_HTT100 OR Ele24_HTT100 OR Ele26_HTT100;Nr. of Jets;Efficiency",15,0,15);

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

        int pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3 = eve->pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_;
        int pass_L1_LooseIsoEG28er2p1_HTT100er = eve->pass_L1_LooseIsoEG28er2p1_HTT100er_;
		int pass_L1_EG_all = eve->pass_L1_EG_all_;
		int pass_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3 = eve->pass_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3_;
        int pass_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3 = eve->pass_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3_;
        int pass_L1_LooseIsoEG24er2p1_HTT100er = eve->pass_L1_LooseIsoEG24er2p1_HTT100er_;
        int pass_L1_LooseIsoEG26er2p1_HTT100er = eve->pass_L1_LooseIsoEG26er2p1_HTT100er_;
        int pass_L1_LooseIsoEG30er2p1_HTT100er = eve->pass_L1_LooseIsoEG30er2p1_HTT100er_;

        int pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_OR = 0;
        int pass_L1_LooseIsoEG28er2p1_HTT100er_OR = 0;
        int pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_OR_prescaled = 0;
        int pass_L1_LooseIsoEG28er2p1_HTT100er_OR_prescaled = 0;

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

        // Checking whether reference triggers are fired
        //int pass_L1_Mu_all = eve->pass_L1_Mu_all_;
        //if( pass_L1_Mu_all == 0 )
        //    continue;

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
        	          double deta = jet_eta[i] - electron_eta;
            	      double dphi = TVector2::Phi_mpi_pi(jet_phi[i] - electron_phi);
                	  double dR = sqrt( deta*deta + dphi*dphi );
                	  if (dR > 0.4){
                    	  njets++;
                      	HT = HT + jet_pt[i];
                        if(jet_pt[i] > leading_jet_pt){
						  leading_jet_pt = jet_pt[i];
                          leading_jet_eta = jet_eta[i];
                        }
					  }
              	}
         	 }
		}

        // Single Electron Event Selection (require MET > 100 GeV)
        if(nele_tight == 1){
            if( nele_loose == 1 && met_pt > 100 ){
                if(njets >= 1){
                    N_e++;
					if(pass_L1_EG_all == 1)
						N_e_EG_all++;
                    if(pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3 == 1)
                        N_e_ele30_jet34++;
                    if((pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3 == 1) || (pass_L1_EG_all == 1) || (pass_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3 == 1) || (pass_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3 == 1)) {
                        N_e_ele30_jet34_OR++;
                        pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_OR = 1;
                    }
                    if((pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3 == 1) || (pass_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3 == 1) || (pass_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3 == 1)) {
                        N_e_ele30_jet34_OR_prescaled++;
                        pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_OR_prescaled = 1;
                    }
                    if(pass_L1_LooseIsoEG28er2p1_HTT100er == 1)
                        N_e_ele28_ht100++;
                    if((pass_L1_LooseIsoEG28er2p1_HTT100er == 1) || (pass_L1_EG_all == 1) || (pass_L1_LooseIsoEG24er2p1_HTT100er == 1) || (pass_L1_LooseIsoEG26er2p1_HTT100er == 1) ) {
                        N_e_ele28_ht100_OR++;
                        pass_L1_LooseIsoEG28er2p1_HTT100er_OR = 1;
                    }
                    if((pass_L1_LooseIsoEG28er2p1_HTT100er == 1) || (pass_L1_LooseIsoEG24er2p1_HTT100er == 1) || (pass_L1_LooseIsoEG26er2p1_HTT100er == 1) ) {
                        N_e_ele28_ht100_OR_prescaled++;
                        pass_L1_LooseIsoEG28er2p1_HTT100er_OR_prescaled = 1;
                    }

					Eff_ele_pt_EG_all->Fill(pass_L1_EG_all,electron_pt);
                    Eff_ele_eta_EG_all->Fill(pass_L1_EG_all,electron_eta);
                    Eff_njets_EG_all->Fill(pass_L1_EG_all,njets);
                    Eff_ele_pt_ele30_jet34->Fill(pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3,electron_pt);
                    Eff_ele_eta_ele30_jet34->Fill(pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3,electron_eta);
                    Eff_jet_pt_ele30_jet34->Fill(pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3,leading_jet_pt);
                    Eff_jet_eta_ele30_jet34->Fill(pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3,leading_jet_eta);
                    Eff_njets_ele30_jet34->Fill(pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3,njets);
                    Eff_ele_pt_ele28_ht100->Fill(pass_L1_LooseIsoEG28er2p1_HTT100er,electron_pt);
                    Eff_ele_eta_ele28_ht100->Fill(pass_L1_LooseIsoEG28er2p1_HTT100er,electron_eta);
                    Eff_ht_ele28_ht100->Fill(pass_L1_LooseIsoEG28er2p1_HTT100er,HT);
                    Eff_njets_ele28_ht100->Fill(pass_L1_LooseIsoEG28er2p1_HTT100er,njets);
                    Eff_ele_pt_ele30_jet34_OR->Fill(pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_OR,electron_pt);
                    Eff_ele_eta_ele30_jet34_OR->Fill(pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_OR,electron_eta);
                    Eff_jet_pt_ele30_jet34_OR->Fill(pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_OR,leading_jet_pt);
                    Eff_jet_eta_ele30_jet34_OR->Fill(pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_OR,leading_jet_eta);
                    Eff_njets_ele30_jet34_OR->Fill(pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_OR,njets);
                    Eff_ele_pt_ele28_ht100_OR->Fill(pass_L1_LooseIsoEG28er2p1_HTT100er_OR,electron_pt);
                    Eff_ele_eta_ele28_ht100_OR->Fill(pass_L1_LooseIsoEG28er2p1_HTT100er_OR,electron_eta);
                    Eff_ht_ele28_ht100_OR->Fill(pass_L1_LooseIsoEG28er2p1_HTT100er_OR,HT);
                    Eff_njets_ele28_ht100_OR->Fill(pass_L1_LooseIsoEG28er2p1_HTT100er_OR,njets);
                    Eff_ele_pt_ele30_jet34_OR_prescaled->Fill(pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_OR_prescaled,electron_pt);
                    Eff_ele_eta_ele30_jet34_OR_prescaled->Fill(pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_OR_prescaled,electron_eta);
                    Eff_jet_pt_ele30_jet34_OR_prescaled->Fill(pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_OR_prescaled,leading_jet_pt);
                    Eff_jet_eta_ele30_jet34_OR_prescaled->Fill(pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_OR_prescaled,leading_jet_eta);
                    Eff_njets_ele30_jet34_OR_prescaled->Fill(pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_OR_prescaled,njets);
                    Eff_ele_pt_ele28_ht100_OR_prescaled->Fill(pass_L1_LooseIsoEG28er2p1_HTT100er_OR_prescaled,electron_pt);
                    Eff_ele_eta_ele28_ht100_OR_prescaled->Fill(pass_L1_LooseIsoEG28er2p1_HTT100er_OR_prescaled,electron_eta);
                    Eff_ht_ele28_ht100_OR_prescaled->Fill(pass_L1_LooseIsoEG28er2p1_HTT100er_OR_prescaled,HT);
                    Eff_njets_ele28_ht100_OR_prescaled->Fill(pass_L1_LooseIsoEG28er2p1_HTT100er_OR_prescaled,njets);
                }
            }
        }

   }

   effi_e_EG_all = double(N_e_EG_all) / N_e;
   effi_e_ele30_jet34 = double(N_e_ele30_jet34) / N_e;
   effi_e_ele28_ht100 = double(N_e_ele28_ht100) / N_e;
   effi_e_ele30_jet34_OR = double(N_e_ele30_jet34_OR) / N_e;
   effi_e_ele28_ht100_OR = double(N_e_ele28_ht100_OR) / N_e;
   effi_e_ele30_jet34_OR_prescaled = double(N_e_ele30_jet34_OR_prescaled) / N_e;
   effi_e_ele28_ht100_OR_prescaled = double(N_e_ele28_ht100_OR_prescaled) / N_e;

   t=clock()-t;
	
   std::cout << " Done! " <<((float)t)/CLOCKS_PER_SEC<< std::endl;
   std::cout<<"**********************************************************************************************\n";
   std::cout<<"Total No. of events in all runs : "<<N_all_total<<"\n";
   std::cout<<"Total No. of events in selected runs : "<<N_total<<"\n";
   std::cout<<"No. of events and efficiency passing Single Electron event selection only : "<<N_e<<"\n";
   std::cout<<"No. of events passing Single Electron event selection + OR of all singleEG L1 seeds : "<<N_e_EG_all<<",  efficiency = "<<effi_e_EG_all<<"\n";
   std::cout<<"No. of events passing Single Electron event selection + Ele30_Jet34 L1 seed : "<<N_e_ele30_jet34<<" ,  efficiency = "<<effi_e_ele30_jet34<<"\n";
   std::cout<<"No. of events passing Single Electron event selection + Ele28_HTT100 L1 seed : "<<N_e_ele28_ht100<<" , efficiency = "<<effi_e_ele28_ht100<<"\n";
   std::cout<<"No. of events passing Single Electron event selection + Ele30_Jet34 L1 seed OR Single EG seeds OR Ele26_Jet34 OR Ele28_Jet34 : efficiency = "<<N_e_ele30_jet34_OR<<",  "<<effi_e_ele30_jet34_OR<<"\n";
   std::cout<<"No. of events passing Single Electron event selection + Ele28_HTT100 L1 seed OR Single EG seeds OR Ele24_HTT100 OR Ele26_HTT100 : efficiency = "<<N_e_ele28_ht100_OR<<",  "<<effi_e_ele28_ht100_OR<<"\n";
   std::cout<<"No. of events passing Single Electron event selection + Ele30_Jet34 L1 seed OR Ele26_Jet34 OR Ele28_Jet34  : efficiency = "<<N_e_ele30_jet34_OR_prescaled<<",  "<<effi_e_ele30_jet34_OR_prescaled<<"\n";
   std::cout<<"No. of events passing Single Electron event selection + Ele28_HTT100 L1 seed OR Ele24_HTT100 OR Ele26_HTT100 : efficiency = "<<N_e_ele28_ht100_OR_prescaled<<",  "<<effi_e_ele28_ht100_OR_prescaled<<"\n";
   std::cout<<"**********************************************************************************************\n";
   histofile.Write();
   histofile.Close();
  
}
