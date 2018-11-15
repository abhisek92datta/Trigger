#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH1D.h"
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
#include <cstdio>
#include <sstream>
#include "TGraphAsymmErrors.h"
#include "TVector.h"
#include "TLorentzVector.h"
#include "TEfficiency.h"
#include "THStack.h"


void Compare_Efficiencies_HLT_all_runs() {

    TH1::SetDefaultSumw2();
    TFile *fileA;
    TFile *fileB;
	TFile *fileC;
	TFile *fileD;

    // List of Files

    ifstream fin;
    fin.open("files_HLT_all_runs.txt");
    char file_runs[200][200];
    int nfiles = 0;

    while(!fin.eof()){
        fin>>file_runs[nfiles];
        nfiles++;
    }
    nfiles--;
    fin.close();

    // List of Histos
    fin.open("histos_HLT.txt");
    char histonames[200][200];
    int nhistos = 0;

    while(!fin.eof()){
        fin>>histonames[nhistos];
        nhistos++;
    }
    fin.close();

    //Looping over files
    for(int i=0; i<nfiles; i++ ){

        char fnameA[200];
		char fnameB[200];
		char fnameC[200];
		char fnameD[200];
        std::sprintf(fnameA, "%sA_Promptreco.root", file_runs[i]);
		std::sprintf(fnameB, "%sB_Promptreco.root", file_runs[i]);
		std::sprintf(fnameC, "%sC_Promptreco.root", file_runs[i]);
		std::sprintf(fnameD, "%sD_Promptreco.root", file_runs[i]);
        fileA = new TFile(fnameA);
		fileB = new TFile(fnameB);
		fileC = new TFile(fnameC);
		fileD = new TFile(fnameD);
		std::string type1, type2, type3;
		std::string fname = fnameA;

        if(fname.find("Global_MET") != std::string::npos)
            type1 = "Global_MET";
		else if(fname.find("Global_SM") != std::string::npos)
			type1 = "Global_SM";
		else if(fname.find("Global_MC") != std::string::npos)
			type1 = "Global_MC";
		else if(fname.find("JetHTLeg_EG") != std::string::npos)
			type1 = "JetHTLeg_EG";
		else if(fname.find("JetHTLeg_MC") != std::string::npos)
			type1 = "JetHTLeg_MC";
		else if(fname.find("EleLeg_JET") != std::string::npos)
			type1 = "EleLeg_JET";
		else if(fname.find("EleLeg_MC_JET") != std::string::npos)
			type1 = "EleLeg_MC_JET";
		else if(fname.find("EleLeg_HT") != std::string::npos)
			type1 = "EleLeg_HT";
		else if(fname.find("EleLeg_MC_HT") != std::string::npos)
			type1 = "EleLeg_MC_HT";
		else
			type1 = "";

		if(fname.find("Promptreco") != std::string::npos)
			type2 = "Promptreco";
		else if(fname.find("Rereco") != std::string::npos)
			type2 = "Rereco";
		else if(fname.find("AODSIM") != std::string::npos)
			type2 = "AODSIM";
		else
			type2 = "";

		if(fname.find("ttjets") != std::string::npos)
			type3 = "ttjets";
		else
			type3 = "2018_AtoD";

		std::vector<TEfficiency*> effi_ele_pt_ele27;
		std::vector<TEfficiency*> effi_ele_eta_ele27;
		std::vector<TEfficiency*> effi_njets_ele27;
		std::vector<std::string> name_ele_pt_ele27;
		std::vector<std::string> name_ele_eta_ele27;
		std::vector<std::string> name_njets_ele27;
        std::vector<TEfficiency*> effi_ele_pt_ele32;
        std::vector<TEfficiency*> effi_ele_eta_ele32;
        std::vector<TEfficiency*> effi_njets_ele32;
        std::vector<std::string> name_ele_pt_ele32;
        std::vector<std::string> name_ele_eta_ele32;
        std::vector<std::string> name_njets_ele32;
		std::vector<TEfficiency*> effi_ele_pt_ele35;
		std::vector<TEfficiency*> effi_ele_eta_ele35;
		std::vector<TEfficiency*> effi_njets_ele35;
		std::vector<std::string> name_ele_pt_ele35;
		std::vector<std::string> name_ele_eta_ele35;
		std::vector<std::string> name_njets_ele35;
        std::vector<TEfficiency*> effi_ele_pt_elejet;
		std::vector<TEfficiency*> effi_ele_eta_elejet;
        std::vector<TEfficiency*> effi_jet_pt_elejet;
		std::vector<TEfficiency*> effi_jet_eta_elejet;
        std::vector<TEfficiency*> effi_ht_elejet;
        std::vector<TEfficiency*> effi_njets_elejet;
        std::vector<std::string> name_ele_pt_elejet;
		std::vector<std::string> name_ele_eta_elejet;
        std::vector<std::string> name_jet_pt_elejet;
		std::vector<std::string> name_jet_eta_elejet;
        std::vector<std::string> name_ht_elejet;
        std::vector<std::string> name_njets_elejet;
		std::vector<TEfficiency*> effi_ele_pt_eleht;
		std::vector<TEfficiency*> effi_ele_eta_eleht;
		std::vector<TEfficiency*> effi_jet_pt_eleht;
		std::vector<TEfficiency*> effi_jet_eta_eleht;
		std::vector<TEfficiency*> effi_ht_eleht;
		std::vector<TEfficiency*> effi_njets_eleht;
		std::vector<std::string> name_ele_pt_eleht;
		std::vector<std::string> name_ele_eta_eleht;
		std::vector<std::string> name_jet_pt_eleht;
		std::vector<std::string> name_jet_eta_eleht;
		std::vector<std::string> name_ht_eleht;
		std::vector<std::string> name_njets_eleht;

        std::string name = "";

        for(int j=0; j<nhistos; j++){

            std::string hname = histonames[j];

            TEfficiency *hA = (TEfficiency*)fileA->Get((hname).c_str());
            TEfficiency *hB = (TEfficiency*)fileB->Get((hname).c_str());
			TEfficiency *hC = (TEfficiency*)fileC->Get((hname).c_str());
			TEfficiency *hD = (TEfficiency*)fileD->Get((hname).c_str());
            hA->SetTitle((hname).c_str());
            hA->SetName((hname).c_str());
            hB->SetTitle((hname).c_str());
            hB->SetName((hname).c_str());
			hC->SetTitle((hname).c_str());
			hC->SetName((hname).c_str());
			hD->SetTitle((hname).c_str());
			hD->SetName((hname).c_str());

            if(hname.find("_ele_pt_") != std::string::npos){
                if( (hname.find("jet35") != std::string::npos) && (type1 != "EleLeg_HT") && (type1 != "EleLeg_MC_HT") ){
					//name = "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned";
                    effi_ele_pt_elejet.push_back(hA);
                    effi_ele_pt_elejet.push_back(hB);
					effi_ele_pt_elejet.push_back(hC);
					effi_ele_pt_elejet.push_back(hD);
                    name_ele_pt_elejet.push_back("Run 2018A");
                    name_ele_pt_elejet.push_back("Run 2018B");
					name_ele_pt_elejet.push_back("Run 2018C");
					name_ele_pt_elejet.push_back("Run 2018D");
                }
                else if( (hname.find("ht150") != std::string::npos) && (type1 != "EleLeg_JET") && (type1 != "EleLeg_MC_JET") ){
					//name = "HLT_Ele28_eta2p1_WPTight_Gsf_HT150";
                    effi_ele_pt_eleht.push_back(hA);
					effi_ele_pt_eleht.push_back(hB);
					effi_ele_pt_eleht.push_back(hC);
					effi_ele_pt_eleht.push_back(hD);
                    name_ele_pt_eleht.push_back("Run 2018A");
                    name_ele_pt_eleht.push_back("Run 2018B");
					name_ele_pt_eleht.push_back("Run 2018C");
					name_ele_pt_eleht.push_back("Run 2018D");
                }
				else if( (hname.find("ele32") != std::string::npos) ){
					//name = "Ele32_WPTight_Gsf";
                    effi_ele_pt_ele32.push_back(hA);
					effi_ele_pt_ele32.push_back(hB);
					effi_ele_pt_ele32.push_back(hC);
					effi_ele_pt_ele32.push_back(hD);
                    name_ele_pt_ele32.push_back("Run 2018A");
					name_ele_pt_ele32.push_back("Run 2018B");
					name_ele_pt_ele32.push_back("Run 2018C");
					name_ele_pt_ele32.push_back("Run 2018D");
				}
				else if( (hname.find("ele35") != std::string::npos) ){
					//name = "Ele35_WPTight_Gsf";
                    effi_ele_pt_ele35.push_back(hA);
					effi_ele_pt_ele35.push_back(hB);
					effi_ele_pt_ele35.push_back(hC);
					effi_ele_pt_ele35.push_back(hD);
                    name_ele_pt_ele35.push_back("Run 2018A");
					name_ele_pt_ele35.push_back("Run 2018B");
					name_ele_pt_ele35.push_back("Run 2018C");
					name_ele_pt_ele35.push_back("Run 2018D");
				}
            }

			if(hname.find("_ele_eta_") != std::string::npos){
				if( (hname.find("jet35") != std::string::npos) && (type1 != "EleLeg_HT") && (type1 != "EleLeg_MC_HT") ){
					//name = "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned";
                    effi_ele_eta_elejet.push_back(hA);
					effi_ele_eta_elejet.push_back(hB);
					effi_ele_eta_elejet.push_back(hC);
					effi_ele_eta_elejet.push_back(hD);
                    name_ele_eta_elejet.push_back("Run 2018A");
					name_ele_eta_elejet.push_back("Run 2018B");
					name_ele_eta_elejet.push_back("Run 2018C");
					name_ele_eta_elejet.push_back("Run 2018D");
				}
				else if( (hname.find("ht150") != std::string::npos) && (type1 != "EleLeg_JET") && (type1 != "EleLeg_MC_JET") ){
					//name = "HLT_Ele28_eta2p1_WPTight_Gsf_HT150";
                    effi_ele_eta_eleht.push_back(hA);
					effi_ele_eta_eleht.push_back(hB);
					effi_ele_eta_eleht.push_back(hC);
					effi_ele_eta_eleht.push_back(hD);
                    name_ele_eta_eleht.push_back("Run 2018A");
					name_ele_eta_eleht.push_back("Run 2018B");
					name_ele_eta_eleht.push_back("Run 2018C");
					name_ele_eta_eleht.push_back("Run 2018D");
				}
				else if( (hname.find("ele32") != std::string::npos) ){
					//name = "Ele32_WPTight_Gsf";
                    effi_ele_eta_ele32.push_back(hA);
					effi_ele_eta_ele32.push_back(hB);
					effi_ele_eta_ele32.push_back(hC);
					effi_ele_eta_ele32.push_back(hD);
                    name_ele_eta_ele32.push_back("Run 2018A");
					name_ele_eta_ele32.push_back("Run 2018B");
					name_ele_eta_ele32.push_back("Run 2018C");
					name_ele_eta_ele32.push_back("Run 2018D");
				}
				else if( (hname.find("ele35") != std::string::npos) ){
					//name = "Ele35_WPTight_Gsf";
                    effi_ele_eta_ele35.push_back(hA);
					effi_ele_eta_ele35.push_back(hB);
					effi_ele_eta_ele35.push_back(hC);
					effi_ele_eta_ele35.push_back(hD);
                    name_ele_eta_ele35.push_back("Run 2018A");
					name_ele_eta_ele35.push_back("Run 2018B");
					name_ele_eta_ele35.push_back("Run 2018C");
					name_ele_eta_ele35.push_back("Run 2018D");
				}
			}

			if(hname.find("_jet_pt_") != std::string::npos){
				if( (hname.find("jet35") != std::string::npos) && (type1 != "EleLeg_HT") && (type1 != "EleLeg_MC_HT") ){
					//name = "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned";
                    effi_jet_pt_elejet.push_back(hA);
					effi_jet_pt_elejet.push_back(hB);
					effi_jet_pt_elejet.push_back(hC);
					effi_jet_pt_elejet.push_back(hD);
                    name_jet_pt_elejet.push_back("Run 2018A");
					name_jet_pt_elejet.push_back("Run 2018B");
					name_jet_pt_elejet.push_back("Run 2018C");
					name_jet_pt_elejet.push_back("Run 2018D");
				}
				else if( (hname.find("ht150") != std::string::npos) && (type1 != "EleLeg_JET") && (type1 != "EleLeg_MC_JET") ){
					//name = "HLT_Ele28_eta2p1_WPTight_Gsf_HT150";
                    effi_jet_pt_eleht.push_back(hA);
					effi_jet_pt_eleht.push_back(hB);
					effi_jet_pt_eleht.push_back(hC);
					effi_jet_pt_eleht.push_back(hD);
                    name_jet_pt_eleht.push_back("Run 2018A");
					name_jet_pt_eleht.push_back("Run 2018B");
					name_jet_pt_eleht.push_back("Run 2018C");
					name_jet_pt_eleht.push_back("Run 2018D");
				}
			}

			if(hname.find("_jet_eta_") != std::string::npos){
				if( (hname.find("jet35") != std::string::npos) && (type1 != "EleLeg_HT") && (type1 != "EleLeg_MC_HT") ){
					//name = "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned";
                    effi_jet_eta_elejet.push_back(hA);
					effi_jet_eta_elejet.push_back(hB);
					effi_jet_eta_elejet.push_back(hC);
					effi_jet_eta_elejet.push_back(hD);
                    name_jet_eta_elejet.push_back("Run 2018A");
                    name_jet_eta_elejet.push_back("Run 2018B");
					name_jet_eta_elejet.push_back("Run 2018C");
					name_jet_eta_elejet.push_back("Run 2018D");
				}
				else if( (hname.find("ht150") != std::string::npos) && (type1 != "EleLeg_JET") && (type1 != "EleLeg_MC_JET") ){
					//name = "HLT_Ele28_eta2p1_WPTight_Gsf_HT150";
                    effi_jet_eta_eleht.push_back(hA);
					effi_jet_eta_eleht.push_back(hB);
					effi_jet_eta_eleht.push_back(hC);
					effi_jet_eta_eleht.push_back(hD);
                    name_jet_eta_eleht.push_back("Run 2018A");
					name_jet_eta_eleht.push_back("Run 2018B");
					name_jet_eta_eleht.push_back("Run 2018C");
					name_jet_eta_eleht.push_back("Run 2018D");
				}
			}

			if(hname.find("_ht_") != std::string::npos){
				if( (hname.find("jet35") != std::string::npos) && (type1 != "EleLeg_HT") && (type1 != "EleLeg_MC_HT") ){
					//name = "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned";
					effi_ht_elejet.push_back(hA);
                    effi_ht_elejet.push_back(hB);
					effi_ht_elejet.push_back(hC);
					effi_ht_elejet.push_back(hD);
                    name_ht_elejet.push_back("Run 2018A");
					name_ht_elejet.push_back("Run 2018B");
					name_ht_elejet.push_back("Run 2018C");
					name_ht_elejet.push_back("Run 2018D");
				}
				else if( (hname.find("ht150") != std::string::npos) && (type1 != "EleLeg_JET") && (type1 != "EleLeg_MC_JET") ){
					//name = "HLT_Ele28_eta2p1_WPTight_Gsf_HT150";
                    effi_ht_eleht.push_back(hA);
					effi_ht_eleht.push_back(hB);
					effi_ht_eleht.push_back(hC);
					effi_ht_eleht.push_back(hD);
                    name_ht_eleht.push_back("Run 2018A");
					name_ht_eleht.push_back("Run 2018B");
					name_ht_eleht.push_back("Run 2018C");
					name_ht_eleht.push_back("Run 2018D");
				}
			}

            if(hname.find("_njets_") != std::string::npos){
				if( (hname.find("jet35") != std::string::npos) && (type1 != "EleLeg_HT") && (type1 != "EleLeg_MC_HT") ){
					//name = "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned";
                    effi_njets_elejet.push_back(hA);
					effi_njets_elejet.push_back(hB);
					effi_njets_elejet.push_back(hC);
					effi_njets_elejet.push_back(hD);
                    name_njets_elejet.push_back("Run 2018A");
					name_njets_elejet.push_back("Run 2018B");
					name_njets_elejet.push_back("Run 2018C");
					name_njets_elejet.push_back("Run 2018D");
				}
				else if( (hname.find("ht150") != std::string::npos) && (type1 != "EleLeg_JET") && (type1 != "EleLeg_MC_JET") ){
					//name = "HLT_Ele28_eta2p1_WPTight_Gsf_HT150";
                    effi_njets_eleht.push_back(hA);
					effi_njets_eleht.push_back(hB);
					effi_njets_eleht.push_back(hC);
					effi_njets_eleht.push_back(hD);
                    name_njets_eleht.push_back("Run 2018A");
					name_njets_eleht.push_back("Run 2018B");
					name_njets_eleht.push_back("Run 2018C");
					name_njets_eleht.push_back("Run 2018D");
				}
				else if( (hname.find("ele32") != std::string::npos) ){
					//name = "Ele32_WPTight_Gsf";
                    effi_njets_ele32.push_back(hA);
					effi_njets_ele32.push_back(hB);
					effi_njets_ele32.push_back(hC);
					effi_njets_ele32.push_back(hD);
                    name_njets_ele32.push_back("Run 2018A");
					name_njets_ele32.push_back("Run 2018B");
					name_njets_ele32.push_back("Run 2018C");
					name_njets_ele32.push_back("Run 2018D");
				}
				else if( (hname.find("ele35") != std::string::npos) ){
					//name = "Ele35_WPTight_Gsf";
                    effi_njets_ele35.push_back(hA);
					effi_njets_ele35.push_back(hB);
					effi_njets_ele35.push_back(hC);
					effi_njets_ele35.push_back(hD);
                    name_njets_ele35.push_back("Run 2018A");
					name_njets_ele35.push_back("Run 2018B");
					name_njets_ele35.push_back("Run 2018C");
					name_njets_ele35.push_back("Run 2018D");
				}
            }

        }

        gStyle->SetLegendBorderSize(1);
        gStyle->SetLegendTextSize(0.027);

		/*
		// Efficiency vs Electron Pt for Ele27 Trigger
		if(effi_ele_pt_ele27.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(20,0,300,1.03,("DATA " + type1 + " : HLT_Ele27_WPTight_Gsf Efficiency vs Electron pT; pT (Electron) [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.60,0.15,0.80,0.35);
			TLine *line1 = new TLine(20,1,300,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_ele_pt_ele27.size(); u++){
				leg1->AddEntry(effi_ele_pt_ele27[u],(name_ele_pt_ele27[u]).c_str(),"L");
				effi_ele_pt_ele27[u]->SetLineWidth(1);
				effi_ele_pt_ele27[u]->SetLineColor(u+4);
				effi_ele_pt_ele27[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_Ele27_" + type1 + "_" + type2 + "_" + type3 + "_" + "ele_pt_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

		// Efficiency vs Electron Eta for Ele27 Trigger
		if(effi_ele_eta_ele27.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(-3,0,3,1.03,("DATA " + type1 + " : HLT_Ele27_WPTight_Gsf Efficiency vs Electron eta; eta (Electron) ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.60,0.15,0.80,0.35);
			TLine *line1 = new TLine(-3,1,3,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_ele_eta_ele27.size(); u++){
				leg1->AddEntry(effi_ele_eta_ele27[u],(name_ele_eta_ele27[u]).c_str(),"L");
				effi_ele_eta_ele27[u]->SetLineWidth(1);
				effi_ele_eta_ele27[u]->SetLineColor(u+4);
				effi_ele_eta_ele27[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_Ele27_" + type1 + "_" + type2 + "_" + type3 + "_" + "ele_eta_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

		 // Efficiency vs Nr. of jets for Ele27 Trigger
		 if(effi_njets_ele27.size() != 0){
		 TCanvas *c1 = new TCanvas("c1","test",1100,650);
		 c1->DrawFrame(0,0,10,1.03,("DATA " + type1 + " : HLT_Ele27_WPTight_Gsf Efficiency vs Jet Multiplicity; Nr. of Jets ; Efficiency").c_str());
		 c1->SetGrid();
		 TLegend* leg1;
		 leg1 = new TLegend(0.60,0.15,0.80,0.35);
		 TLine *line1 = new TLine(0,1,10,1);
		 line1->SetLineStyle(kDashed);
		 line1->SetLineWidth(1);
		 leg1->SetFillColor(kWhite);
		 leg1->SetFillStyle(1001);
		 for(unsigned int u=0; u<effi_njets_ele27.size(); u++){
		 leg1->AddEntry(effi_njets_ele27[u],(name_njets_ele27[u]).c_str(),"L");
		 effi_njets_ele27[u]->SetLineWidth(1);
		 effi_njets_ele27[u]->SetLineColor(u+4);
		 effi_njets_ele27[u]->Draw("L same");
		 }
		 leg1->Draw("same");
		 line1->Draw("same");
		 c1->Print(("HLT_Ele27_" + type1 + "_" + type2 + "_" + type3 + "_" + "njets_turn_on.pdf").c_str());
		 delete c1;
		 delete leg1;
		 delete line1;
		 }
		*/

        // Efficiency vs Electron Pt for Ele32 Trigger
        if(effi_ele_pt_ele32.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(20,0,300,1.03,("DATA " + type1 + " : HLT_Ele32_WPTight_Gsf Efficiency vs Electron pT; pT (Electron) [GeV] ; Efficiency").c_str());
            c1->SetGrid();
            TLegend* leg1;
            leg1 = new TLegend(0.60,0.15,0.80,0.35);
            TLine *line1 = new TLine(20,1,300,1);
            line1->SetLineStyle(kDashed);
            line1->SetLineWidth(1);
            leg1->SetFillColor(kWhite);
            leg1->SetFillStyle(1001);
            for(unsigned int u=0; u<effi_ele_pt_ele32.size(); u++){
                leg1->AddEntry(effi_ele_pt_ele32[u],(name_ele_pt_ele32[u]).c_str(),"L");
                effi_ele_pt_ele32[u]->SetLineWidth(1);
                effi_ele_pt_ele32[u]->SetLineColor(u+4);
                effi_ele_pt_ele32[u]->Draw("L same");
            }
            leg1->Draw("same");
            line1->Draw("same");
            c1->Print(("HLT_Ele32_" + type1 + "_" + type2 + "_" + type3 + "_" + "ele_pt_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

        // Efficiency vs Electron Eta for Ele32 Trigger
        if(effi_ele_eta_ele32.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(-3,0,3,1.03,("DATA " + type1 + " : HLT_Ele32_WPTight_Gsf Efficiency vs Electron eta; eta (Electron) ; Efficiency").c_str());
            c1->SetGrid();
            TLegend* leg1;
            leg1 = new TLegend(0.60,0.15,0.80,0.35);
            TLine *line1 = new TLine(-3,1,3,1);
            line1->SetLineStyle(kDashed);
            line1->SetLineWidth(1);
            leg1->SetFillColor(kWhite);
            leg1->SetFillStyle(1001);
            for(unsigned int u=0; u<effi_ele_eta_ele32.size(); u++){
                leg1->AddEntry(effi_ele_eta_ele32[u],(name_ele_eta_ele32[u]).c_str(),"L");
                effi_ele_eta_ele32[u]->SetLineWidth(1);
                effi_ele_eta_ele32[u]->SetLineColor(u+4);
                effi_ele_eta_ele32[u]->Draw("L same");
            }
            leg1->Draw("same");
            line1->Draw("same");
            c1->Print(("HLT_Ele32_" + type1 + "_" + type2 + "_" + type3 + "_" + "ele_eta_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

        // Efficiency vs Nr. of jets for Ele32 Trigger
        if(effi_njets_ele32.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(0,0,10,1.03,("DATA " + type1 + " : HLT_Ele32_WPTight_Gsf Efficiency vs Jet Multiplicity; Nr. of Jets ; Efficiency").c_str());
            c1->SetGrid();
            TLegend* leg1;
            leg1 = new TLegend(0.60,0.15,0.80,0.35);
            TLine *line1 = new TLine(0,1,10,1);
            line1->SetLineStyle(kDashed);
            line1->SetLineWidth(1);
            leg1->SetFillColor(kWhite);
            leg1->SetFillStyle(1001);
            for(unsigned int u=0; u<effi_njets_ele32.size(); u++){
                leg1->AddEntry(effi_njets_ele32[u],(name_njets_ele32[u]).c_str(),"L");
                effi_njets_ele32[u]->SetLineWidth(1);
                effi_njets_ele32[u]->SetLineColor(u+4);
                effi_njets_ele32[u]->Draw("L same");
            }
            leg1->Draw("same");
            line1->Draw("same");
            c1->Print(("HLT_Ele32_" + type1 + "_" + type2 + "_" + type3 + "_" + "njets_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

		// Efficiency vs Electron Pt for Ele35 Trigger
		if(effi_ele_pt_ele35.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(20,0,300,1.03,("DATA " + type1 + " : HLT_Ele35_WPTight_Gsf Efficiency vs Electron pT; pT (Electron) [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.60,0.15,0.80,0.35);
			TLine *line1 = new TLine(20,1,300,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_ele_pt_ele35.size(); u++){
				leg1->AddEntry(effi_ele_pt_ele35[u],(name_ele_pt_ele35[u]).c_str(),"L");
				effi_ele_pt_ele35[u]->SetLineWidth(1);
				effi_ele_pt_ele35[u]->SetLineColor(u+4);
				effi_ele_pt_ele35[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_Ele35_" + type1 + "_" + type2 + "_" + type3 + "_" + "ele_pt_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

		// Efficiency vs Electron Eta for Ele35 Trigger
		if(effi_ele_eta_ele35.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(-3,0,3,1.03,("DATA " + type1 + " : HLT_Ele35_WPTight_Gsf Efficiency vs Electron eta; eta (Electron) ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.60,0.15,0.80,0.35);
			TLine *line1 = new TLine(-3,1,3,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_ele_eta_ele35.size(); u++){
				leg1->AddEntry(effi_ele_eta_ele35[u],(name_ele_eta_ele35[u]).c_str(),"L");
				effi_ele_eta_ele35[u]->SetLineWidth(1);
				effi_ele_eta_ele35[u]->SetLineColor(u+4);
				effi_ele_eta_ele35[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_Ele35_" + type1 + "_" + type2 + "_" + type3 + "_" + "ele_eta_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

		// Efficiency vs Nr. of jets for Ele35 Trigger
		if(effi_njets_ele35.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(0,0,10,1.03,("DATA " + type1 + " : HLT_Ele35_WPTight_Gsf Efficiency vs Jet Multiplicity; Nr. of Jets ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.60,0.15,0.80,0.35);
			TLine *line1 = new TLine(0,1,10,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_njets_ele35.size(); u++){
				leg1->AddEntry(effi_njets_ele35[u],(name_njets_ele35[u]).c_str(),"L");
				effi_njets_ele35[u]->SetLineWidth(1);
				effi_njets_ele35[u]->SetLineColor(u+4);
				effi_njets_ele35[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_Ele35_" + type1 + "_" + type2 + "_" + type3 + "_" + "njets_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

        // Efficiency vs Electron Pt for Ele+Jet Trigger
        if(effi_ele_pt_elejet.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(20,0,300,1.03,("DATA " + type1 + " : HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned Efficiency vs Electron pT; pT (Electron) [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.60,0.15,0.80,0.35);
            TLine *line1 = new TLine(20,1,300,1);
            line1->SetLineStyle(kDashed);
            line1->SetLineWidth(1);
            leg1->SetFillColor(kWhite);
            leg1->SetFillStyle(1001);
            for(unsigned int u=0; u<effi_ele_pt_elejet.size(); u++){
                leg1->AddEntry(effi_ele_pt_elejet[u],(name_ele_pt_elejet[u]).c_str(),"L");
                effi_ele_pt_elejet[u]->SetLineWidth(1);
                effi_ele_pt_elejet[u]->SetLineColor(u+4);
                effi_ele_pt_elejet[u]->Draw("L same");
            }
            leg1->Draw("same");
            line1->Draw("same");
            c1->Print(("HLT_EleJet_" + type1 + "_" + type2 + "_" + type3 + "_" + "ele_pt_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

		// Efficiency vs Electron Eta for Ele+Jet Trigger
		if(effi_ele_eta_elejet.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(-3,0,3,1.03,("DATA " + type1 + " : HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned Efficiency vs Electron eta; eta (Electron) ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.60,0.15,0.80,0.35);
			TLine *line1 = new TLine(-3,1,3,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_ele_eta_elejet.size(); u++){
				leg1->AddEntry(effi_ele_eta_elejet[u],(name_ele_eta_elejet[u]).c_str(),"L");
				effi_ele_eta_elejet[u]->SetLineWidth(1);
				effi_ele_eta_elejet[u]->SetLineColor(u+4);
				effi_ele_eta_elejet[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_EleJet_" + type1 + "_" + type2 + "_" + type3 + "_" + "ele_eta_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

        // Efficiency vs Jet Pt for Ele+Jet Trigger
        if(effi_jet_pt_elejet.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(20,0,300,1.03,("DATA " + type1 + " : HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned Efficiency vs Leading Jet pT; pT (Leading Jet) [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.60,0.15,0.80,0.35);
            TLine *line1 = new TLine(20,1,300,1);
            line1->SetLineStyle(kDashed);
            line1->SetLineWidth(1);
            leg1->SetFillColor(kWhite);
            leg1->SetFillStyle(1001);
            for(unsigned int u=0; u<effi_jet_pt_elejet.size(); u++){
                leg1->AddEntry(effi_jet_pt_elejet[u],(name_jet_pt_elejet[u]).c_str(),"L");
                effi_jet_pt_elejet[u]->SetLineWidth(1);
                effi_jet_pt_elejet[u]->SetLineColor(u+2);
                effi_jet_pt_elejet[u]->Draw("L same");
            }
            leg1->Draw("same");
            line1->Draw("same");
            c1->Print(("HLT_EleJet_" + type1 + "_" + type2 + "_" + type3 + "_" + "jet_pt_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

		// Efficiency vs Jet Eta for Ele+Jet Trigger
		if(effi_jet_eta_elejet.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(-3,0,3,1.03,("DATA " + type1 + " : HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned Efficiency vs Leading Jet eta; eta (Leading Jet) ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.60,0.15,0.80,0.35);
			TLine *line1 = new TLine(-3,1,3,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_jet_eta_elejet.size(); u++){
				leg1->AddEntry(effi_jet_eta_elejet[u],(name_jet_eta_elejet[u]).c_str(),"L");
				effi_jet_eta_elejet[u]->SetLineWidth(1);
				effi_jet_eta_elejet[u]->SetLineColor(u+4);
				effi_jet_eta_elejet[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_EleJet_" + type1 + "_" + type2 + "_" + type3 + "_" + "jet_eta_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

        // Efficiency vs HT for Ele+Jet Trigger
        if(effi_ht_elejet.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(0,0,1000,1.03,("DATA " + type1 + " : HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned Efficiency vs HT; HT [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.60,0.15,0.80,0.35);
            TLine *line1 = new TLine(0,1,1000,1);
            line1->SetLineStyle(kDashed);
            line1->SetLineWidth(1);
            leg1->SetFillColor(kWhite);
            leg1->SetFillStyle(1001);
            for(unsigned int u=0; u<effi_ht_elejet.size(); u++){
                leg1->AddEntry(effi_ht_elejet[u],(name_ht_elejet[u]).c_str(),"L");
                effi_ht_elejet[u]->SetLineWidth(1);
                effi_ht_elejet[u]->SetLineColor(u+2);
                effi_ht_elejet[u]->Draw("L same");
            }
            leg1->Draw("same");
            line1->Draw("same");
            c1->Print(("HLT_EleJet_" + type1 + "_" + type2 + "_" + type3 + "_" + "ht_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

        // Efficiency vs Nr. of jets for Ele+Jet Trigger
        if(effi_njets_elejet.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(0,0,10,1.03,("DATA " + type1 + " : HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned Efficiency vs Jet Multiplicity; Nr. of Jets ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.60,0.15,0.80,0.35);
            TLine *line1 = new TLine(0,1,10,1);
            line1->SetLineStyle(kDashed);
            line1->SetLineWidth(1);
            leg1->SetFillColor(kWhite);
            leg1->SetFillStyle(1001);
            for(unsigned int u=0; u<effi_njets_elejet.size(); u++){
                leg1->AddEntry(effi_njets_elejet[u],(name_njets_elejet[u]).c_str(),"L");
                effi_njets_elejet[u]->SetLineWidth(1);
                effi_njets_elejet[u]->SetLineColor(u+4);
                effi_njets_elejet[u]->Draw("L same");
            }
            leg1->Draw("same");
            line1->Draw("same");
            c1->Print(("HLT_EleJet_" + type1 + "_" + type2 + "_" + type3 + "_" + "njets_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

		// Efficiency vs Electron Pt for Ele+HT Trigger
		if(effi_ele_pt_eleht.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(20,0,300,1.03,("DATA " + type1 + " : HLT_Ele28_eta2p1_WPTight_Gsf_HT150 Efficiency vs Electron pT; pT (Electron) [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.60,0.15,0.80,0.35);
			TLine *line1 = new TLine(20,1,300,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_ele_pt_eleht.size(); u++){
				leg1->AddEntry(effi_ele_pt_eleht[u],(name_ele_pt_eleht[u]).c_str(),"L");
				effi_ele_pt_eleht[u]->SetLineWidth(1);
				effi_ele_pt_eleht[u]->SetLineColor(u+4);
				effi_ele_pt_eleht[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_EleHT_" + type1 + "_" + type2 + "_" + type3 + "_" + "ele_pt_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

		// Efficiency vs Electron Eta for Ele+HT Trigger
		if(effi_ele_eta_eleht.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(-3,0,3,1.03,("DATA " + type1 + " : HLT_Ele28_eta2p1_WPTight_Gsf_HT150 Efficiency vs Electron eta; eta (Electron) ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.60,0.15,0.80,0.35);
			TLine *line1 = new TLine(-3,1,3,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_ele_eta_eleht.size(); u++){
				leg1->AddEntry(effi_ele_eta_eleht[u],(name_ele_eta_eleht[u]).c_str(),"L");
				effi_ele_eta_eleht[u]->SetLineWidth(1);
				effi_ele_eta_eleht[u]->SetLineColor(u+4);
				effi_ele_eta_eleht[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_EleHT_" + type1 + "_" + type2 + "_" + type3 + "_" + "ele_eta_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

		// Efficiency vs Jet Pt for Ele+HT Trigger
		if(effi_jet_pt_eleht.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(20,0,300,1.03,("DATA " + type1 + " : HLT_Ele28_eta2p1_WPTight_Gsf_HT150 Efficiency vs Leading Jet pT; pT (Leading Jet) [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.60,0.15,0.80,0.35);
			TLine *line1 = new TLine(20,1,300,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_jet_pt_eleht.size(); u++){
				leg1->AddEntry(effi_jet_pt_eleht[u],(name_jet_pt_eleht[u]).c_str(),"L");
				effi_jet_pt_eleht[u]->SetLineWidth(1);
				effi_jet_pt_eleht[u]->SetLineColor(u+2);
				effi_jet_pt_eleht[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_EleHT_" + type1 + "_" + type2 + "_" + type3 + "_" + "jet_pt_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

		// Efficiency vs Jet Eta for Ele+HT Trigger
		if(effi_jet_eta_eleht.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(-3,0,3,1.03,("DATA " + type1 + " : HLT_Ele28_eta2p1_WPTight_Gsf_HT150 Efficiency vs Leading Jet eta; eta (Leading Jet) ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.60,0.15,0.80,0.35);
			TLine *line1 = new TLine(-3,1,3,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_jet_eta_eleht.size(); u++){
				leg1->AddEntry(effi_jet_eta_eleht[u],(name_jet_eta_eleht[u]).c_str(),"L");
				effi_jet_eta_eleht[u]->SetLineWidth(1);
				effi_jet_eta_eleht[u]->SetLineColor(u+4);
				effi_jet_eta_eleht[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_EleHT_" + type1 + "_" + type2 + "_" + type3 + "_" + "jet_eta_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

		// Efficiency vs HT for Ele+HT Trigger
		if(effi_ht_eleht.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(0,0,1000,1.03,("DATA " + type1 + " : HLT_Ele28_eta2p1_WPTight_Gsf_HT150 Efficiency vs HT; HT [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.60,0.15,0.80,0.35);
			TLine *line1 = new TLine(0,1,1000,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_ht_eleht.size(); u++){
				leg1->AddEntry(effi_ht_eleht[u],(name_ht_eleht[u]).c_str(),"L");
				effi_ht_eleht[u]->SetLineWidth(1);
				effi_ht_eleht[u]->SetLineColor(u+2);
				effi_ht_eleht[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_EleHT_" + type1 + "_" + type2 + "_" + type3 + "_" + "ht_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

		// Efficiency vs Nr. of jets for Ele+HT Trigger
		if(effi_njets_eleht.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(0,0,10,1.03,("DATA " + type1 + " : HLT_Ele28_eta2p1_WPTight_Gsf_HT150 Efficiency vs Jet Multiplicity; Nr. of Jets ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.60,0.15,0.80,0.35);
			TLine *line1 = new TLine(0,1,10,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_njets_eleht.size(); u++){
				leg1->AddEntry(effi_njets_eleht[u],(name_njets_eleht[u]).c_str(),"L");
				effi_njets_eleht[u]->SetLineWidth(1);
				effi_njets_eleht[u]->SetLineColor(u+4);
				effi_njets_eleht[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_EleHT_" + type1 + "_" + type2 + "_" + type3 + "_" + "njets_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

    }

    delete fileA;
    delete fileB;
	delete fileC;
	delete fileD;

	return;
}


