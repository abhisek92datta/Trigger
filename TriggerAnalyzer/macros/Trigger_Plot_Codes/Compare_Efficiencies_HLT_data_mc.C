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


void Compare_Efficiencies_HLT_data_mc() {

    TH1::SetDefaultSumw2();
    TFile *file_data;
	TFile *file_mc;

    // List of Files

    ifstream fin;
    fin.open("files_HLT_data_mc.txt");
    char filenames_data[200][200];
	char filenames_mc[200][200];
    int nfiles = 0;

    while(!fin.eof()){
        fin>>filenames_data[nfiles];
		fin>>filenames_mc[nfiles];
        nfiles++;
    }
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

		std::string fname_data = filenames_data[i];
		std::string fname_mc = filenames_mc[i];
        file_data = new TFile(filenames_data[i]);
		file_mc = new TFile(filenames_mc[i]);
        std::string type1, type2, leg_name;

        if(fname_data.find("Global_MET") != std::string::npos)
            type1 = "Global_MET";
		if(fname_data.find("Global_SM") != std::string::npos)
			type1 = "Global_SM";
		else if(fname_data.find("JetHTLeg_EG") != std::string::npos)
			type1 = "JetHTLeg_EG";
		else if(fname_data.find("EleLeg_JET") != std::string::npos)
			type1 = "EleLeg_JET";
		else if(fname_data.find("EleLeg_HT") != std::string::npos)
			type1 = "EleLeg_HT";
		else
			type1 = "";

		if(fname_data.find("2018ABCD_Promptreco") != std::string::npos){
			type2 = "2018_DataABCD_Promptreco_MC";
			leg_name = "Data Runs 2018A to D";
		}
		else{
			type2 = "";
			leg_name = "";
		}

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

            TEfficiency *h_data = (TEfficiency*)file_data->Get((hname).c_str());
			TEfficiency *h_mc = (TEfficiency*)file_mc->Get((hname).c_str());
            h_data->SetTitle((hname).c_str());
            h_data->SetName((hname).c_str());
			h_mc->SetTitle((hname).c_str());
			h_mc->SetName((hname).c_str());

            if(hname.find("_ele_pt_") != std::string::npos){
                if( (hname.find("jet35") != std::string::npos) && (type1 != "EleLeg_HT") && (type1 != "EleLeg_MC_HT") ){
					//name = "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned";
                    effi_ele_pt_elejet.push_back(h_data);
					effi_ele_pt_elejet.push_back(h_mc);
                    name_ele_pt_elejet.push_back(leg_name);
					name_ele_pt_elejet.push_back("MC ttbar 2018");
                }
                else if( (hname.find("ht150") != std::string::npos) && (type1 != "EleLeg_JET") && (type1 != "EleLeg_MC_JET") ){
					//name = "HLT_Ele28_eta2p1_WPTight_Gsf_HT150";
					effi_ele_pt_eleht.push_back(h_data);
					effi_ele_pt_eleht.push_back(h_mc);
					name_ele_pt_eleht.push_back(leg_name);
					name_ele_pt_eleht.push_back("MC ttbar 2018");
                }
				else if( (hname.find("ele32") != std::string::npos) ){
					effi_ele_pt_ele32.push_back(h_data);
					effi_ele_pt_ele32.push_back(h_mc);
					name_ele_pt_ele32.push_back(leg_name);
					name_ele_pt_ele32.push_back("MC ttbar 2018");
				}
				else if( (hname.find("ele35") != std::string::npos) ){
					effi_ele_pt_ele35.push_back(h_data);
					effi_ele_pt_ele35.push_back(h_mc);
					name_ele_pt_ele35.push_back(leg_name);
					name_ele_pt_ele35.push_back("MC ttbar 2018");
				}
            }

			if(hname.find("_ele_eta_") != std::string::npos){
				if( (hname.find("jet35") != std::string::npos) && (type1 != "EleLeg_HT") && (type1 != "EleLeg_MC_HT") ){
					//name = "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned";
					effi_ele_eta_elejet.push_back(h_data);
					effi_ele_eta_elejet.push_back(h_mc);
					name_ele_eta_elejet.push_back(leg_name);
					name_ele_eta_elejet.push_back("MC ttbar 2018");
				}
				else if( (hname.find("ht150") != std::string::npos) && (type1 != "EleLeg_JET") && (type1 != "EleLeg_MC_JET") ){
					//name = "HLT_Ele28_eta2p1_WPTight_Gsf_HT150";
					effi_ele_eta_eleht.push_back(h_data);
					effi_ele_eta_eleht.push_back(h_mc);
					name_ele_eta_eleht.push_back(leg_name);
					name_ele_eta_eleht.push_back("MC ttbar 2018");
				}
				else if( (hname.find("ele32") != std::string::npos) ){
					effi_ele_eta_ele32.push_back(h_data);
					effi_ele_eta_ele32.push_back(h_mc);
					name_ele_eta_ele32.push_back(leg_name);
					name_ele_eta_ele32.push_back("MC ttbar 2018");
				}
				else if( (hname.find("ele35") != std::string::npos) ){
					effi_ele_eta_ele35.push_back(h_data);
					effi_ele_eta_ele35.push_back(h_mc);
					name_ele_eta_ele35.push_back(leg_name);
					name_ele_eta_ele35.push_back("MC ttbar 2018");
				}
			}

			if(hname.find("_jet_pt_") != std::string::npos){
				if( (hname.find("jet35") != std::string::npos) && (type1 != "EleLeg_HT") && (type1 != "EleLeg_MC_HT") ){
					//name = "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned";
					effi_jet_pt_elejet.push_back(h_data);
					effi_jet_pt_elejet.push_back(h_mc);
					name_jet_pt_elejet.push_back(leg_name);
					name_jet_pt_elejet.push_back("MC ttbar 2018");
				}
				else if( (hname.find("ht150") != std::string::npos) && (type1 != "EleLeg_JET") && (type1 != "EleLeg_MC_JET") ){
					//name = "HLT_Ele28_eta2p1_WPTight_Gsf_HT150";
					effi_jet_pt_eleht.push_back(h_data);
					effi_jet_pt_eleht.push_back(h_mc);
					name_jet_pt_eleht.push_back(leg_name);
					name_jet_pt_eleht.push_back("MC ttbar 2018");
				}
			}

			if(hname.find("_jet_eta_") != std::string::npos){
				if( (hname.find("jet35") != std::string::npos) && (type1 != "EleLeg_HT") && (type1 != "EleLeg_MC_HT") ){
					//name = "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned";
					effi_jet_eta_elejet.push_back(h_data);
					effi_jet_eta_elejet.push_back(h_mc);
					name_jet_eta_elejet.push_back(leg_name);
					name_jet_eta_elejet.push_back("MC ttbar 2018");
				}
				else if( (hname.find("ht150") != std::string::npos) && (type1 != "EleLeg_JET") && (type1 != "EleLeg_MC_JET") ){
					//name = "HLT_Ele28_eta2p1_WPTight_Gsf_HT150";
					effi_jet_eta_eleht.push_back(h_data);
					effi_jet_eta_eleht.push_back(h_mc);
					name_jet_eta_eleht.push_back(leg_name);
					name_jet_eta_eleht.push_back("MC ttbar 2018");
				}
			}

			if(hname.find("_ht_") != std::string::npos){
				if( (hname.find("jet35") != std::string::npos) && (type1 != "EleLeg_HT") && (type1 != "EleLeg_MC_HT") ){
					//name = "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned";
					effi_ht_elejet.push_back(h_data);
					effi_ht_elejet.push_back(h_mc);
					name_ht_elejet.push_back(leg_name);
					name_ht_elejet.push_back("MC ttbar 2018");
				}
				else if( (hname.find("ht150") != std::string::npos) && (type1 != "EleLeg_JET") && (type1 != "EleLeg_MC_JET") ){
					//name = "HLT_Ele28_eta2p1_WPTight_Gsf_HT150";
					effi_ht_eleht.push_back(h_data);
					effi_ht_eleht.push_back(h_mc);
					name_ht_eleht.push_back(leg_name);
					name_ht_eleht.push_back("MC ttbar 2018");
				}
			}

            if(hname.find("_njets_") != std::string::npos){
				if( (hname.find("jet35") != std::string::npos) && (type1 != "EleLeg_HT") && (type1 != "EleLeg_MC_HT") ){
					//name = "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned";
					effi_njets_elejet.push_back(h_data);
					effi_njets_elejet.push_back(h_mc);
					name_njets_elejet.push_back(leg_name);
					name_njets_elejet.push_back("MC ttbar 2018");
				}
				else if( (hname.find("ht150") != std::string::npos) && (type1 != "EleLeg_JET") && (type1 != "EleLeg_MC_JET") ){
					//name = "HLT_Ele28_eta2p1_WPTight_Gsf_HT150";
					effi_njets_eleht.push_back(h_data);
					effi_njets_eleht.push_back(h_mc);
					name_njets_eleht.push_back(leg_name);
					name_njets_eleht.push_back("MC ttbar 2018");
				}
				else if( (hname.find("ele32") != std::string::npos) ){
					effi_njets_ele32.push_back(h_data);
					effi_njets_ele32.push_back(h_mc);
					name_njets_ele32.push_back(leg_name);
					name_njets_ele32.push_back("MC ttbar 2018");
				}
				else if( (hname.find("ele35") != std::string::npos) ){
					effi_njets_ele35.push_back(h_data);
					effi_njets_ele35.push_back(h_mc);
					name_njets_ele35.push_back(leg_name);
					name_njets_ele35.push_back("MC ttbar 2018");
				}
            }

        }

        gStyle->SetLegendBorderSize(1);
        gStyle->SetLegendTextSize(0.027);

		/*
		// Efficiency vs Electron Pt for Ele27 Trigger
		if(effi_ele_pt_ele27.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(20,0,300,1.03,("DATA_MC " + type1 + " : Ele27_WPTight_Gsf Efficiency vs Electron pT; pT (Electron) [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.55,0.15,0.80,0.30);
			TLine *line1 = new TLine(20,1,300,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_ele_pt_ele27.size(); u++){
				leg1->AddEntry(effi_ele_pt_ele27[u],(name_ele_pt_ele27[u]).c_str(),"L");
				effi_ele_pt_ele27[u]->SetLineWidth(1);
				effi_ele_pt_ele27[u]->SetLineColor(u+6);
				effi_ele_pt_ele27[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_Ele27_" + type1 + "_" + type2 + "_ele_pt_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

		// Efficiency vs Electron Eta for Ele27 Trigger
		if(effi_ele_eta_ele27.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(-3,0,3,1.03,("DATA_MC " + type1 + " : Ele27_WPTight_Gsf Efficiency vs Electron eta; eta (Electron) ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.55,0.15,0.80,0.30);
			TLine *line1 = new TLine(-3,1,3,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_ele_eta_ele27.size(); u++){
				leg1->AddEntry(effi_ele_eta_ele27[u],(name_ele_eta_ele27[u]).c_str(),"L");
				effi_ele_eta_ele27[u]->SetLineWidth(1);
				effi_ele_eta_ele27[u]->SetLineColor(u+6);
				effi_ele_eta_ele27[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_Ele27_" + type1 + "_" + type2 + "_ele_eta_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

		 // Efficiency vs Nr. of jets for Ele27 Trigger
		 if(effi_njets_ele27.size() != 0){
		 TCanvas *c1 = new TCanvas("c1","test",1100,650);
		 c1->DrawFrame(0,0,10,1.03,("DATA_MC " + type1 + " : Ele27_WPTight_Gsf Efficiency vs Jet Multiplicity; Nr. of Jets ; Efficiency").c_str());
		 c1->SetGrid();
		 TLegend* leg1;
		 leg1 = new TLegend(0.55,0.15,0.80,0.30);
		 TLine *line1 = new TLine(0,1,10,1);
		 line1->SetLineStyle(kDashed);
		 line1->SetLineWidth(1);
		 leg1->SetFillColor(kWhite);
		 leg1->SetFillStyle(1001);
		 for(unsigned int u=0; u<effi_njets_ele27.size(); u++){
		 leg1->AddEntry(effi_njets_ele27[u],(name_njets_ele27[u]).c_str(),"L");
		 effi_njets_ele27[u]->SetLineWidth(1);
		 effi_njets_ele27[u]->SetLineColor(u+6);
		 effi_njets_ele27[u]->Draw("L same");
		 }
		 leg1->Draw("same");
		 line1->Draw("same");
		 c1->Print(("HLT_Ele27_" + type1 + "_" + type2 + "_njets_turn_on.pdf").c_str());
		 delete c1;
		 delete leg1;
		 delete line1;
		 }
		*/

        // Efficiency vs Electron pT for Ele32 Trigger

        if(effi_ele_pt_ele32.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(20,0,300,1.03,("DATA_MC " + type1 + " : Ele32_WPTight_Gsf Efficiency vs Electron pT; pT (Electron) [GeV] ; Efficiency").c_str());
            c1->SetGrid();
            TLegend* leg1;
            leg1 = new TLegend(0.55,0.15,0.80,0.30);
            TLine *line1 = new TLine(20,1,300,1);
            line1->SetLineStyle(kDashed);
            line1->SetLineWidth(1);
            leg1->SetFillColor(kWhite);
            leg1->SetFillStyle(1001);
            for(unsigned int u=0; u<effi_ele_pt_ele32.size(); u++){
                leg1->AddEntry(effi_ele_pt_ele32[u],(name_ele_pt_ele32[u]).c_str(),"L");
                effi_ele_pt_ele32[u]->SetLineWidth(1);
                effi_ele_pt_ele32[u]->SetLineColor(u+6);
                effi_ele_pt_ele32[u]->Draw("L same");
            }
            leg1->Draw("same");
            line1->Draw("same");
            c1->Print(("HLT_Ele32_" + type1 + "_" + type2 + "_ele_pt_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

        // Efficiency vs Electron Eta for Ele32 Trigger
        if(effi_ele_eta_ele32.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(-3,0,3,1.03,("DATA_MC " + type1 + " : Ele32_WPTight_Gsf Efficiency vs Electron eta; eta (Electron) ; Efficiency").c_str());
            c1->SetGrid();
            TLegend* leg1;
            leg1 = new TLegend(0.55,0.15,0.80,0.30);
            TLine *line1 = new TLine(-3,1,3,1);
            line1->SetLineStyle(kDashed);
            line1->SetLineWidth(1);
            leg1->SetFillColor(kWhite);
            leg1->SetFillStyle(1001);
            for(unsigned int u=0; u<effi_ele_eta_ele32.size(); u++){
                leg1->AddEntry(effi_ele_eta_ele32[u],(name_ele_eta_ele32[u]).c_str(),"L");
                effi_ele_eta_ele32[u]->SetLineWidth(1);
                effi_ele_eta_ele32[u]->SetLineColor(u+6);
                effi_ele_eta_ele32[u]->Draw("L same");
            }
            leg1->Draw("same");
            line1->Draw("same");
            c1->Print(("HLT_Ele32_" + type1 + "_" + type2 + "_ele_eta_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

        // Efficiency vs Nr. of jets for Ele32 Trigger
        if(effi_njets_ele32.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(0,0,10,1.03,("DATA_MC " + type1 + " : Ele32_WPTight_Gsf Efficiency vs Jet Multiplicity; Nr. of Jets ; Efficiency").c_str());
            c1->SetGrid();
            TLegend* leg1;
            leg1 = new TLegend(0.55,0.15,0.80,0.30);
            TLine *line1 = new TLine(0,1,10,1);
            line1->SetLineStyle(kDashed);
            line1->SetLineWidth(1);
            leg1->SetFillColor(kWhite);
            leg1->SetFillStyle(1001);
            for(unsigned int u=0; u<effi_njets_ele32.size(); u++){
                leg1->AddEntry(effi_njets_ele32[u],(name_njets_ele32[u]).c_str(),"L");
                effi_njets_ele32[u]->SetLineWidth(1);
                effi_njets_ele32[u]->SetLineColor(u+6);
                effi_njets_ele32[u]->Draw("L same");
            }
            leg1->Draw("same");
            line1->Draw("same");
            c1->Print(("HLT_Ele32_" + type1 + "_" + type2 + "_njets_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

		 // Efficiency vs Electron Pt for Ele35 Trigger
		 if(effi_ele_pt_ele35.size() != 0){
		 TCanvas *c1 = new TCanvas("c1","test",1100,650);
		 c1->DrawFrame(20,0,300,1.03,("DATA_MC " + type1 + " : Ele35_WPTight_Gsf Efficiency vs Electron pT; pT (Electron) [GeV] ; Efficiency").c_str());
		 c1->SetGrid();
		 TLegend* leg1;
		 leg1 = new TLegend(0.55,0.15,0.80,0.30);
		 TLine *line1 = new TLine(20,1,300,1);
		 line1->SetLineStyle(kDashed);
		 line1->SetLineWidth(1);
		 leg1->SetFillColor(kWhite);
		 leg1->SetFillStyle(1001);
		 for(unsigned int u=0; u<effi_ele_pt_ele35.size(); u++){
		 leg1->AddEntry(effi_ele_pt_ele35[u],(name_ele_pt_ele35[u]).c_str(),"L");
		 effi_ele_pt_ele35[u]->SetLineWidth(1);
		 effi_ele_pt_ele35[u]->SetLineColor(u+6);
		 effi_ele_pt_ele35[u]->Draw("L same");
		 }
		 leg1->Draw("same");
		 line1->Draw("same");
		 c1->Print(("HLT_Ele35_" + type1 + "_" + type2 + "_ele_pt_turn_on.pdf").c_str());
		 delete c1;
		 delete leg1;
		 delete line1;
		 }

		 // Efficiency vs Electron Eta for Ele35 Trigger
		 if(effi_ele_eta_ele35.size() != 0){
		 TCanvas *c1 = new TCanvas("c1","test",1100,650);
		 c1->DrawFrame(-3,0,3,1.03,("DATA_MC " + type1 + " : Ele35_WPTight_Gsf Efficiency vs Electron eta; eta (Electron) ; Efficiency").c_str());
		 c1->SetGrid();
		 TLegend* leg1;
		 leg1 = new TLegend(0.55,0.15,0.80,0.30);
		 TLine *line1 = new TLine(-3,1,3,1);
		 line1->SetLineStyle(kDashed);
		 line1->SetLineWidth(1);
		 leg1->SetFillColor(kWhite);
		 leg1->SetFillStyle(1001);
		 for(unsigned int u=0; u<effi_ele_eta_ele35.size(); u++){
		 leg1->AddEntry(effi_ele_eta_ele35[u],(name_ele_eta_ele35[u]).c_str(),"L");
		 effi_ele_eta_ele35[u]->SetLineWidth(1);
		 effi_ele_eta_ele35[u]->SetLineColor(u+6);
		 effi_ele_eta_ele35[u]->Draw("L same");
		 }
		 leg1->Draw("same");
		 line1->Draw("same");
		 c1->Print(("HLT_Ele35_" + type1 + "_" + type2 + "_ele_eta_turn_on.pdf").c_str());
		 delete c1;
		 delete leg1;
		 delete line1;
		 }

		// Efficiency vs Nr. of jets for Ele35 Trigger
		if(effi_njets_ele35.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(0,0,10,1.03,("DATA_MC " + type1 + " : Ele35_WPTight_Gsf Efficiency vs Jet Multiplicity; Nr. of Jets ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.55,0.15,0.80,0.30);
			TLine *line1 = new TLine(0,1,10,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_njets_ele35.size(); u++){
				leg1->AddEntry(effi_njets_ele35[u],(name_njets_ele35[u]).c_str(),"L");
				effi_njets_ele35[u]->SetLineWidth(1);
				effi_njets_ele35[u]->SetLineColor(u+6);
				effi_njets_ele35[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_Ele35_" + type1 + "_" + type2 + "_njets_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

        // Efficiency vs Electron Pt for Ele+Jet Trigger
        if(effi_ele_pt_elejet.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(20,0,300,1.03,("DATA_MC " + type1 + " : HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned Efficiency vs Electron pT; pT (Electron) [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.55,0.15,0.80,0.30);
            TLine *line1 = new TLine(20,1,300,1);
            line1->SetLineStyle(kDashed);
            line1->SetLineWidth(1);
            leg1->SetFillColor(kWhite);
            leg1->SetFillStyle(1001);
            for(unsigned int u=0; u<effi_ele_pt_elejet.size(); u++){
                leg1->AddEntry(effi_ele_pt_elejet[u],(name_ele_pt_elejet[u]).c_str(),"L");
                effi_ele_pt_elejet[u]->SetLineWidth(1);
                effi_ele_pt_elejet[u]->SetLineColor(u+6);
                effi_ele_pt_elejet[u]->Draw("L same");
            }
            leg1->Draw("same");
            line1->Draw("same");
            c1->Print(("HLT_EleJet_" + type1 + "_" + type2 + "_ele_pt_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

		// Efficiency vs Electron Eta for Ele+Jet Trigger
		if(effi_ele_eta_elejet.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(-3,0,3,1.03,("DATA_MC " + type1 + " : HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned Efficiency vs Electron eta; eta (Electron) ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.55,0.15,0.80,0.30);
			TLine *line1 = new TLine(-3,1,3,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_ele_eta_elejet.size(); u++){
				leg1->AddEntry(effi_ele_eta_elejet[u],(name_ele_eta_elejet[u]).c_str(),"L");
				effi_ele_eta_elejet[u]->SetLineWidth(1);
				effi_ele_eta_elejet[u]->SetLineColor(u+6);
				effi_ele_eta_elejet[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_EleJet_" + type1 + "_" + type2 + "_ele_eta_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

        // Efficiency vs Jet Pt for Ele+Jet Trigger
        if(effi_jet_pt_elejet.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(20,0,300,1.03,("DATA_MC " + type1 + " : HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned Efficiency vs Leading Jet pT; pT (Leading Jet) [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.55,0.15,0.80,0.30);
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
            c1->Print(("HLT_EleJet_" + type1 + "_" + type2 + "_jet_pt_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

		// Efficiency vs Jet Eta for Ele+Jet Trigger
		if(effi_jet_eta_elejet.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(-3,0,3,1.03,("DATA_MC " + type1 + " : HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned Efficiency vs Leading Jet eta; eta (Leading Jet) ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.55,0.15,0.80,0.30);
			TLine *line1 = new TLine(-3,1,3,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_jet_eta_elejet.size(); u++){
				leg1->AddEntry(effi_jet_eta_elejet[u],(name_jet_eta_elejet[u]).c_str(),"L");
				effi_jet_eta_elejet[u]->SetLineWidth(1);
				effi_jet_eta_elejet[u]->SetLineColor(u+6);
				effi_jet_eta_elejet[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_EleJet_" + type1 + "_" + type2 + "_jet_eta_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

        // Efficiency vs HT for Ele+Jet Trigger
        if(effi_ht_elejet.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(0,0,1000,1.03,("DATA_MC " + type1 + " : HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned Efficiency vs HT; HT [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.55,0.15,0.80,0.30);
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
            c1->Print(("HLT_EleJet_" + type1 + "_" + type2 + "_ht_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

        // Efficiency vs Nr. of jets for Ele+Jet Trigger
        if(effi_njets_elejet.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(0,0,10,1.03,("DATA_MC " + type1 + " : HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned Efficiency vs Jet Multiplicity; Nr. of Jets ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.55,0.15,0.80,0.30);
            TLine *line1 = new TLine(0,1,10,1);
            line1->SetLineStyle(kDashed);
            line1->SetLineWidth(1);
            leg1->SetFillColor(kWhite);
            leg1->SetFillStyle(1001);
            for(unsigned int u=0; u<effi_njets_elejet.size(); u++){
                leg1->AddEntry(effi_njets_elejet[u],(name_njets_elejet[u]).c_str(),"L");
                effi_njets_elejet[u]->SetLineWidth(1);
                effi_njets_elejet[u]->SetLineColor(u+6);
                effi_njets_elejet[u]->Draw("L same");
            }
            leg1->Draw("same");
            line1->Draw("same");
            c1->Print(("HLT_EleJet_" + type1 + "_" + type2 + "_njets_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

		// Efficiency vs Electron Pt for Ele+HT Trigger
		if(effi_ele_pt_eleht.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(20,0,300,1.03,("DATA_MC " + type1 + " : HLT_Ele28_eta2p1_WPTight_Gsf_HT150 Efficiency vs Electron pT; pT (Electron) [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.55,0.15,0.80,0.30);
			TLine *line1 = new TLine(20,1,300,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_ele_pt_eleht.size(); u++){
				leg1->AddEntry(effi_ele_pt_eleht[u],(name_ele_pt_eleht[u]).c_str(),"L");
				effi_ele_pt_eleht[u]->SetLineWidth(1);
				effi_ele_pt_eleht[u]->SetLineColor(u+6);
				effi_ele_pt_eleht[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_EleHT_" + type1 + "_" + type2 + "_ele_pt_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

		// Efficiency vs Electron Eta for Ele+HT Trigger
		if(effi_ele_eta_eleht.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(-3,0,3,1.03,("DATA_MC " + type1 + " : HLT_Ele28_eta2p1_WPTight_Gsf_HT150 Efficiency vs Electron eta; eta (Electron) ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.55,0.15,0.80,0.30);
			TLine *line1 = new TLine(-3,1,3,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_ele_eta_eleht.size(); u++){
				leg1->AddEntry(effi_ele_eta_eleht[u],(name_ele_eta_eleht[u]).c_str(),"L");
				effi_ele_eta_eleht[u]->SetLineWidth(1);
				effi_ele_eta_eleht[u]->SetLineColor(u+6);
				effi_ele_eta_eleht[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_EleHT_" + type1 + "_" + type2 + "_ele_eta_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

		// Efficiency vs Jet Pt for Ele+HT Trigger
		if(effi_jet_pt_eleht.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(20,0,300,1.03,("DATA_MC " + type1 + " : HLT_Ele28_eta2p1_WPTight_Gsf_HT150 Efficiency vs Leading Jet pT; pT (Leading Jet) [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.55,0.15,0.80,0.30);
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
			c1->Print(("HLT_EleHT_" + type1 + "_" + type2 + "_jet_pt_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

		// Efficiency vs Jet Eta for Ele+HT Trigger
		if(effi_jet_eta_eleht.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(-3,0,3,1.03,("DATA_MC " + type1 + " : HLT_Ele28_eta2p1_WPTight_Gsf_HT150 Efficiency vs Leading Jet eta; eta (Leading Jet) ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.55,0.15,0.80,0.30);
			TLine *line1 = new TLine(-3,1,3,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_jet_eta_eleht.size(); u++){
				leg1->AddEntry(effi_jet_eta_eleht[u],(name_jet_eta_eleht[u]).c_str(),"L");
				effi_jet_eta_eleht[u]->SetLineWidth(1);
				effi_jet_eta_eleht[u]->SetLineColor(u+6);
				effi_jet_eta_eleht[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_EleHT_" + type1 + "_" + type2 + "_jet_eta_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

		// Efficiency vs HT for Ele+HT Trigger
		if(effi_ht_eleht.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(0,0,1000,1.03,("DATA_MC " + type1 + " : HLT_Ele28_eta2p1_WPTight_Gsf_HT150 Efficiency vs HT; HT [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.55,0.15,0.80,0.30);
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
			c1->Print(("HLT_EleHT_" + type1 + "_" + type2 + "_ht_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

		// Efficiency vs Nr. of jets for Ele+HT Trigger
		if(effi_njets_eleht.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(0,0,10,1.03,("DATA_MC " + type1 + " : HLT_Ele28_eta2p1_WPTight_Gsf_HT150 Efficiency vs Jet Multiplicity; Nr. of Jets ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.55,0.15,0.80,0.30);
			TLine *line1 = new TLine(0,1,10,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_njets_eleht.size(); u++){
				leg1->AddEntry(effi_njets_eleht[u],(name_njets_eleht[u]).c_str(),"L");
				effi_njets_eleht[u]->SetLineWidth(1);
				effi_njets_eleht[u]->SetLineColor(u+6);
				effi_njets_eleht[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_EleHT_" + type1 + "_" + type2 + "_njets_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

    }

    delete file_data;
	delete file_mc;

	return;
}


