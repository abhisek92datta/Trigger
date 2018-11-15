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


void Compare_Efficiencies_L1_data_mc() {

    TH1::SetDefaultSumw2();
	TFile *file_data;
	TFile *file_mc;

    // List of Files

    ifstream fin;
    fin.open("files_L1_data_mc.txt");
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

    fin.open("histos_L1.txt");
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

        if(fname_data.find("Global_SM") != std::string::npos)
            type1 = "Global_SM";
        else if(fname_data.find("Global_MET") != std::string::npos)
            type1 = "Global_MET";
		else
			type1 = "";

		if(fname_data.find("2018ABCD_Promptreco") != std::string::npos){
			type2 = "2018ABCD_Promptreco";
			leg_name = "Data Runs 2018A to D : ";
		}
        else if(fname_data.find("2018ABCD_Rereco") != std::string::npos){
            type2 = "2018ABCD_Rereco";
            leg_name = "Data Runs 2018A to D : ";
        }
		else{
			type2 = "";
			leg_name = "";
		}

		std::vector<TEfficiency*> effi_ele_pt_EGall;
        std::vector<TEfficiency*> effi_ele_eta_EGall;
		std::vector<std::string> name_ele_pt_EGall;
        std::vector<std::string> name_ele_eta_EGall;
        std::vector<TEfficiency*> effi_ele_pt_elejet;
        std::vector<TEfficiency*> effi_ele_eta_elejet;
        std::vector<TEfficiency*> effi_jet_pt_elejet;
        std::vector<TEfficiency*> effi_jet_eta_elejet;
        std::vector<TEfficiency*> effi_njets_elejet;
        std::vector<std::string> name_ele_pt_elejet;
        std::vector<std::string> name_ele_eta_elejet;
        std::vector<std::string> name_jet_pt_elejet;
        std::vector<std::string> name_jet_eta_elejet;
        std::vector<std::string> name_njets_elejet;
		std::vector<TEfficiency*> effi_ele_pt_eleht;
        std::vector<TEfficiency*> effi_ele_eta_eleht;
		std::vector<TEfficiency*> effi_ht_eleht;
		std::vector<TEfficiency*> effi_njets_eleht;
		std::vector<std::string> name_ele_pt_eleht;
        std::vector<std::string> name_ele_eta_eleht;
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
                if(hname.find("jet34") != std::string::npos){
                    if(hname.find("OR") != std::string::npos)
                        name = "IsoEG30er_Jet34 OR all SingleEG";
                    else
                        name = "IsoEG30er_Jet34";
					effi_ele_pt_elejet.push_back(h_data);
					effi_ele_pt_elejet.push_back(h_mc);
					name_ele_pt_elejet.push_back((leg_name + name).c_str());
					name_ele_pt_elejet.push_back(("MC ttbar 2018 : " + name).c_str());
                }
                if(hname.find("ht100") != std::string::npos){
                    if(hname.find("OR") != std::string::npos)
                        name = "IsoEG28er_HTT100 OR all SingleEG";
                    else
                        name = "IsoEG28er_HTT100";
					effi_ele_pt_eleht.push_back(h_data);
					effi_ele_pt_eleht.push_back(h_mc);
					name_ele_pt_eleht.push_back((leg_name + name).c_str());
					name_ele_pt_eleht.push_back(("MC ttbar 2018 : " + name).c_str());
                }
				if(hname.find("EG_all") != std::string::npos){
					name = "OR all SingleEG";
					effi_ele_pt_EGall.push_back(h_data);
					effi_ele_pt_EGall.push_back(h_mc);
					name_ele_pt_EGall.push_back((leg_name + name).c_str());
					name_ele_pt_EGall.push_back(("MC ttbar 2018 : " + name).c_str());
				}
            }

            if(hname.find("_ele_eta_") != std::string::npos){
                if(hname.find("jet34") != std::string::npos){
                    if(hname.find("OR") != std::string::npos)
                        name = "IsoEG30er_Jet34 OR all SingleEG";
                    else
                        name = "IsoEG30er_Jet34";
                    effi_ele_eta_elejet.push_back(h_data);
                    effi_ele_eta_elejet.push_back(h_mc);
                    name_ele_eta_elejet.push_back((leg_name + name).c_str());
                    name_ele_eta_elejet.push_back(("MC ttbar 2018 : " + name).c_str());
                }
                if(hname.find("ht100") != std::string::npos){
                    if(hname.find("OR") != std::string::npos)
                        name = "IsoEG28er_HTT100 OR all SingleEG";
                    else
                        name = "IsoEG28er_HTT100";
                    effi_ele_eta_eleht.push_back(h_data);
                    effi_ele_eta_eleht.push_back(h_mc);
                    name_ele_eta_eleht.push_back((leg_name + name).c_str());
                    name_ele_eta_eleht.push_back(("MC ttbar 2018 : " + name).c_str());
                }
                if(hname.find("EG_all") != std::string::npos){
                    name = "OR all SingleEG";
                    effi_ele_eta_EGall.push_back(h_data);
                    effi_ele_eta_EGall.push_back(h_mc);
                    name_ele_eta_EGall.push_back((leg_name + name).c_str());
                    name_ele_eta_EGall.push_back(("MC ttbar 2018 : " + name).c_str());
                }
            }

            if(hname.find("_jet_pt_") != std::string::npos){
                if(hname.find("OR") != std::string::npos)
                    name = "IsoEG30er_Jet34 OR all SingleEG";
                else
                    name = "IsoEG30er_Jet34";
				effi_jet_pt_elejet.push_back(h_data);
				effi_jet_pt_elejet.push_back(h_mc);
				name_jet_pt_elejet.push_back((leg_name + name).c_str());
				name_jet_pt_elejet.push_back(("MC ttbar 2018 : " + name).c_str());
            }

            if(hname.find("_jet_eta_") != std::string::npos){
                if(hname.find("OR") != std::string::npos)
                    name = "IsoEG30er_Jet34 OR all SingleEG";
                else
                    name = "IsoEG30er_Jet34";
                effi_jet_eta_elejet.push_back(h_data);
                effi_jet_eta_elejet.push_back(h_mc);
                name_jet_eta_elejet.push_back((leg_name + name).c_str());
                name_jet_eta_elejet.push_back(("MC ttbar 2018 : " + name).c_str());
            }

            if(hname.find("_ht_") != std::string::npos){
                if(hname.find("OR") != std::string::npos)
                    name = "IsoEG28er_HTT100 OR all SingleEG";
                else
                    name = "IsoEG28er_HTT100";
				effi_ht_eleht.push_back(h_data);
				effi_ht_eleht.push_back(h_mc);
				name_ht_eleht.push_back((leg_name + name).c_str());
				name_ht_eleht.push_back(("MC ttbar 2018 : " + name).c_str());
            }

            if(hname.find("_njets_") != std::string::npos){
                if(hname.find("jet34") != std::string::npos){
                    if(hname.find("OR") != std::string::npos)
                        name = "IsoEG30er_Jet34 OR all SingleEG";
                    else
                        name = "IsoEG30er_Jet34";
					effi_njets_elejet.push_back(h_data);
					effi_njets_elejet.push_back(h_mc);
					name_njets_elejet.push_back((leg_name + name).c_str());
					name_njets_elejet.push_back(("MC ttbar 2018 : " + name).c_str());
                }
                if(hname.find("ht100") != std::string::npos){
                    if(hname.find("OR") != std::string::npos)
                        name = "IsoEG28er_HTT100 OR all SingleEG";
                    else
                        name = "IsoEG28er_HTT100";
					effi_njets_eleht.push_back(h_data);
					effi_njets_eleht.push_back(h_mc);
					name_njets_eleht.push_back((leg_name + name).c_str());
					name_njets_eleht.push_back(("MC ttbar 2018 : " + name).c_str());
                }
            }

        }

        gStyle->SetLegendBorderSize(1);
        gStyle->SetLegendTextSize(0.027);

		// Efficiency vs Electron Pt for OR of all single EG Triggers
		if(effi_ele_pt_EGall.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(20,0.30,150,1.03,("DATA_MC " + type1 + " : OR of all SingleEG Efficiency vs Electron pT; pT (Electron) [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1 = new TLegend(0.30,0.15,0.85,0.35);
			TLine *line1 = new TLine(20,1,150,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_ele_pt_EGall.size(); u++){
				leg1->AddEntry(effi_ele_pt_EGall[u],(name_ele_pt_EGall[u]).c_str(),"L");
				effi_ele_pt_EGall[u]->SetLineWidth(1);
				effi_ele_pt_EGall[u]->SetLineColor(u+6);
				effi_ele_pt_EGall[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("L1_EGall_" + type1 + "_" + type2 + "_ele_pt_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

        // Efficiency vs Electron Eta for OR of all single EG Triggers
        if(effi_ele_eta_EGall.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(-3,0.30,3,1.03,("DATA_MC " + type1 + " : OR of all SingleEG Efficiency vs Electron Eta; Eta (Electron) ; Efficiency").c_str());
            c1->SetGrid();
            TLegend* leg1 = new TLegend(0.30,0.15,0.85,0.35);
            TLine *line1 = new TLine(20,1,150,1);
            line1->SetLineStyle(kDashed);
            line1->SetLineWidth(1);
            leg1->SetFillColor(kWhite);
            leg1->SetFillStyle(1001);
            for(unsigned int u=0; u<effi_ele_eta_EGall.size(); u++){
                leg1->AddEntry(effi_ele_eta_EGall[u],(name_ele_eta_EGall[u]).c_str(),"L");
                effi_ele_eta_EGall[u]->SetLineWidth(1);
                effi_ele_eta_EGall[u]->SetLineColor(u+6);
                effi_ele_eta_EGall[u]->Draw("L same");
            }
            leg1->Draw("same");
            line1->Draw("same");
            c1->Print(("L1_EGall_" + type1 + "_" + type2 + "_ele_eta_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

        // Efficiency vs Electron Pt for Ele+Jet Trigger
        if(effi_ele_pt_elejet.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(20,0.30,150,1.03,("DATA_MC " + type1 + " : L1_IsoEG30er_Jet34 Efficiency vs Electron pT; pT (Electron) [GeV] ; Efficiency").c_str());
			c1->SetGrid();
            TLegend* leg1 = new TLegend(0.30,0.15,0.85,0.35);
            TLine *line1 = new TLine(20,1,150,1);
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
            c1->Print(("L1_EleJet_" + type1 + "_" + type2 + "_ele_pt_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

        // Efficiency vs Electron Eta for Ele+Jet Trigger
        if(effi_ele_pt_elejet.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(-3,0.30,3,1.03,("DATA_MC " + type1 + " : L1_IsoEG30er_Jet34 Efficiency vs Electron Eta; Eta (Electron) ; Efficiency").c_str());
            c1->SetGrid();
            TLegend* leg1 = new TLegend(0.30,0.15,0.85,0.35);
            TLine *line1 = new TLine(20,1,150,1);
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
            c1->Print(("L1_EleJet_" + type1 + "_" + type2 + "_ele_eta_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

		// Efficiency vs Electron Pt for Ele+HT Trigger
		if(effi_ele_pt_eleht.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(20,0.30,150,1.03,("DATA_MC " + type1 + " : L1_IsoEG28er_HTT100 vs Electron pT; pT (Electron) [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1 = new TLegend(0.30,0.15,0.85,0.35);
			TLine *line1 = new TLine(20,1,150,1);
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
			c1->Print(("L1_EleHT_" + type1 + "_" + type2 + "_ele_pt_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

        // Efficiency vs Electron Eta for Ele+HT Trigger
        if(effi_ele_eta_eleht.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(-3,0.30,3,1.03,("DATA_MC " + type1 + " : L1_IsoEG28er_HTT100 vs Electron Eta; Eta (Electron) ; Efficiency").c_str());
            c1->SetGrid();
            TLegend* leg1 = new TLegend(0.30,0.15,0.85,0.35);
            TLine *line1 = new TLine(20,1,150,1);
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
            c1->Print(("L1_EleHT_" + type1 + "_" + type2 + "_ele_eta_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

        // Efficiency vs Jet Pt for Ele+Jet Trigger
        if(effi_jet_pt_elejet.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(20,0.10,150,1.03,("DATA_MC " + type1 + " : L1_IsoEG30er_Jet34 Efficiency vs Leading Jet pT; pT (Leading Jet) [GeV] ; Efficiency").c_str());
			c1->SetGrid();
            TLegend* leg1 = new TLegend(0.30,0.15,0.85,0.35);
            TLine *line1 = new TLine(20,1,150,1);
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
            c1->Print(("L1_EleJet_" + type1 + "_" + type2 + "_jet_pt_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

        // Efficiency vs Jet Eta for Ele+Jet Trigger
        if(effi_jet_eta_elejet.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(-3,0.30,3,1.03,("DATA_MC " + type1 + " : L1_IsoEG30er_Jet34 Efficiency vs Leading Jet Eta; Eta (Leading Jet) ; Efficiency").c_str());
            c1->SetGrid();
            TLegend* leg1 = new TLegend(0.30,0.15,0.85,0.35);
            TLine *line1 = new TLine(20,1,150,1);
            line1->SetLineStyle(kDashed);
            line1->SetLineWidth(1);
            leg1->SetFillColor(kWhite);
            leg1->SetFillStyle(1001);
            for(unsigned int u=0; u<effi_jet_eta_elejet.size(); u++){
                leg1->AddEntry(effi_jet_eta_elejet[u],(name_jet_eta_elejet[u]).c_str(),"L");
                effi_jet_eta_elejet[u]->SetLineWidth(1);
                effi_jet_eta_elejet[u]->SetLineColor(u+2);
                effi_jet_eta_elejet[u]->Draw("L same");
            }
            leg1->Draw("same");
            line1->Draw("same");
            c1->Print(("L1_EleJet_" + type1 + "_" + type2 + "_jet_eta_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

        // Efficiency vs HT for Ele+HT Trigger
        if(effi_ht_eleht.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(0,0.20,500,1.03,("DATA_MC " + type1 + " : L1_IsoEG28er_HTT100 Efficiency vs HT; HT [GeV] ; Efficiency").c_str());
			c1->SetGrid();
            TLegend* leg1 = new TLegend(0.30,0.15,0.85,0.35);
            TLine *line1 = new TLine(0,1,500,1);
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
            c1->Print(("L1_EleHT_" + type1 + "_" + type2 + "_ht_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

        // Efficiency vs Nr. of jets for Ele+Jet Trigger
        if(effi_njets_elejet.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(0,0.40,10,1.03,("DATA_MC " + type1 + " : L1_IsoEG30er_Jet34 Efficiency vs Jet Multiplicity; Nr. of Jets ; Efficiency").c_str());
  			c1->SetGrid();
			TLegend* leg1 = new TLegend(0.30,0.15,0.85,0.35);
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
            c1->Print(("L1_EleJet_" + type1 + "_" + type2 + "_njets_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

		// Efficiency vs Nr. of jets for Ele+HT Trigger
		if(effi_njets_eleht.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(0,0.40,10,1.03,("DATA_MC " + type1 + " : L1_IsoEG28er_HTT100 Efficiency vs Jet Multiplicity; Nr. of Jets ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1 = new TLegend(0.30,0.15,0.85,0.35);
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
			c1->Print(("L1_EleHT_" + type1 + "_" + type2 + "_njets_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

    }

	delete file_data;
	delete file_mc;
	
	return;
}


