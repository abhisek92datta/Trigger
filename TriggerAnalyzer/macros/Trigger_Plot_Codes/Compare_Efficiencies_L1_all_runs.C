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


void Compare_Efficiencies_L1_all_runs() {

    TH1::SetDefaultSumw2();
    TFile *fileA;
	TFile *fileB;
	TFile *fileC;
	TFile *fileD;

    // List of Files

    ifstream fin;
    fin.open("files_L1_all_runs.txt");
    char file_runs[200][200];
    int nfiles = 0;

    while(!fin.eof()){
        fin>>file_runs[nfiles];
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
		std::string fname = fnameB;

        if(fname.find("Global_SM") != std::string::npos)
            type1 = "Global_SM";
        else if(fname.find("Global_MET") != std::string::npos)
            type1 = "Global_MET";
		else if(fname.find("Global_MC") != std::string::npos)
			type1 = "Global_MC";
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
                if(hname.find("jet34") != std::string::npos){
					if(hname.find("OR") != std::string::npos){
                        name = "IsoEG30er_Jet34 OR all SingleEG";
                        effi_ele_pt_elejet.push_back(hA);
                        effi_ele_pt_elejet.push_back(hB);
						effi_ele_pt_elejet.push_back(hC);
						effi_ele_pt_elejet.push_back(hD);
						name_ele_pt_elejet.push_back("Run 2018A");
						name_ele_pt_elejet.push_back("Run 2018B");
						name_ele_pt_elejet.push_back("Run 2018C");
						name_ele_pt_elejet.push_back("Run 2018D");

					}
                }
                if(hname.find("ht100") != std::string::npos){
					if(hname.find("OR") != std::string::npos){
                        name = "IsoEG28er_HTT100 OR all SingleEG";
                        effi_ele_pt_eleht.push_back(hA);
                        effi_ele_pt_eleht.push_back(hB);
						effi_ele_pt_eleht.push_back(hC);
						effi_ele_pt_eleht.push_back(hD);
                        name_ele_pt_eleht.push_back("Run 2018A");
						name_ele_pt_eleht.push_back("Run 2018B");
						name_ele_pt_eleht.push_back("Run 2018C");
						name_ele_pt_eleht.push_back("Run 2018D");
					}
                }
				if(hname.find("EG_all") != std::string::npos){
						name = "OR all SingleEG";
                        effi_ele_pt_EGall.push_back(hA);
                        effi_ele_pt_EGall.push_back(hB);
						effi_ele_pt_EGall.push_back(hC);
						effi_ele_pt_EGall.push_back(hD);
                        name_ele_pt_EGall.push_back("Run 2018A");
						name_ele_pt_EGall.push_back("Run 2018B");
						name_ele_pt_EGall.push_back("Run 2018C");
						name_ele_pt_EGall.push_back("Run 2018D");
				}
            }

            if(hname.find("_ele_eta_") != std::string::npos){
                if(hname.find("jet34") != std::string::npos){
                    if(hname.find("OR") != std::string::npos){
                        name = "IsoEG30er_Jet34 OR all SingleEG";
                        effi_ele_eta_elejet.push_back(hA);
                        effi_ele_eta_elejet.push_back(hB);
                        effi_ele_eta_elejet.push_back(hC);
                        effi_ele_eta_elejet.push_back(hD);
                        name_ele_eta_elejet.push_back("Run 2018A");
                        name_ele_eta_elejet.push_back("Run 2018B");
                        name_ele_eta_elejet.push_back("Run 2018C");
                        name_ele_eta_elejet.push_back("Run 2018D");

                    }
                }
                if(hname.find("ht100") != std::string::npos){
                    if(hname.find("OR") != std::string::npos){
                        name = "IsoEG28er_HTT100 OR all SingleEG";
                        effi_ele_eta_eleht.push_back(hA);
                        effi_ele_eta_eleht.push_back(hB);
                        effi_ele_eta_eleht.push_back(hC);
                        effi_ele_eta_eleht.push_back(hD);
                        name_ele_eta_eleht.push_back("Run 2018A");
                        name_ele_eta_eleht.push_back("Run 2018B");
                        name_ele_eta_eleht.push_back("Run 2018C");
                        name_ele_eta_eleht.push_back("Run 2018D");
                    }
                }
                if(hname.find("EG_all") != std::string::npos){
                    name = "OR all SingleEG";
                    effi_ele_eta_EGall.push_back(hA);
                    effi_ele_eta_EGall.push_back(hB);
                    effi_ele_eta_EGall.push_back(hC);
                    effi_ele_eta_EGall.push_back(hD);
                    name_ele_eta_EGall.push_back("Run 2018A");
                    name_ele_eta_EGall.push_back("Run 2018B");
                    name_ele_eta_EGall.push_back("Run 2018C");
                    name_ele_eta_EGall.push_back("Run 2018D");
                }
            }

            if(hname.find("_jet_pt_") != std::string::npos){
				if(hname.find("OR") != std::string::npos){
                    name = "IsoEG30er_Jet34 OR all SingleEG";
                    effi_jet_pt_elejet.push_back(hA);
                    effi_jet_pt_elejet.push_back(hB);
					effi_jet_pt_elejet.push_back(hC);
					effi_jet_pt_elejet.push_back(hD);
                    name_jet_pt_elejet.push_back("Run 2018A");
					name_jet_pt_elejet.push_back("Run 2018B");
					name_jet_pt_elejet.push_back("Run 2018C");
					name_jet_pt_elejet.push_back("Run 2018D");
				}
            }

            if(hname.find("_jet_eta_") != std::string::npos){
                if(hname.find("OR") != std::string::npos){
                    name = "IsoEG30er_Jet34 OR all SingleEG";
                    effi_jet_eta_elejet.push_back(hA);
                    effi_jet_eta_elejet.push_back(hB);
                    effi_jet_eta_elejet.push_back(hC);
                    effi_jet_eta_elejet.push_back(hD);
                    name_jet_eta_elejet.push_back("Run 2018A");
                    name_jet_eta_elejet.push_back("Run 2018B");
                    name_jet_eta_elejet.push_back("Run 2018C");
                    name_jet_eta_elejet.push_back("Run 2018D");
                }
            }

            if(hname.find("_ht_") != std::string::npos){
				if(hname.find("OR") != std::string::npos){
                    name = "IsoEG28er_HTT100 OR all SingleEG";
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
                if(hname.find("jet34") != std::string::npos){
					if(hname.find("OR") != std::string::npos){
                        name = "IsoEG30er_Jet34 OR all SingleEG";
                        effi_njets_elejet.push_back(hA);
                        effi_njets_elejet.push_back(hB);
						effi_njets_elejet.push_back(hC);
						effi_njets_elejet.push_back(hD);
                        name_njets_elejet.push_back("Run 2018A");
                        name_njets_elejet.push_back("Run 2018B");
						name_njets_elejet.push_back("Run 2018C");
						name_njets_elejet.push_back("Run 2018D");
					}
                }
                if(hname.find("ht100") != std::string::npos){
					if(hname.find("OR") != std::string::npos){
                        name = "IsoEG28er_HTT100 OR all SingleEG";
                        effi_njets_eleht.push_back(hA);
						effi_njets_eleht.push_back(hB);
						effi_njets_eleht.push_back(hC);
						effi_njets_eleht.push_back(hD);
                        name_njets_eleht.push_back("Run 2018A");
						name_njets_eleht.push_back("Run 2018B");
						name_njets_eleht.push_back("Run 2018C");
						name_njets_eleht.push_back("Run 2018D");
					}
                }
            }

        }

        gStyle->SetLegendBorderSize(1);
        gStyle->SetLegendTextSize(0.027);

		// Efficiency vs Electron Pt for OR of all SingleEG Triggers
		if(effi_ele_pt_EGall.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(20,0.30,150,1.03,("DATA " + type1 + " : OR all SingleEG Efficiency vs Electron pT; pT (Electron) [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1 = new TLegend(0.60,0.15,0.80,0.35);
			TLine *line1 = new TLine(20,1,150,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_ele_pt_EGall.size(); u++){
				leg1->AddEntry(effi_ele_pt_EGall[u],(name_ele_pt_EGall[u]).c_str(),"L");
				effi_ele_pt_EGall[u]->SetLineWidth(1);
				effi_ele_pt_EGall[u]->SetLineColor(u+4);
				effi_ele_pt_EGall[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("L1_EGall_" + type1 + "_" + type2 + "_" + type3 + "_ele_pt_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

        // Efficiency vs Electron Eta for OR of all SingleEG Triggers
        if(effi_ele_eta_EGall.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(-3,0.30,3,1.03,("DATA " + type1 + " : OR all SingleEG Efficiency vs Electron Eta; Eta (Electron) ; Efficiency").c_str());
            c1->SetGrid();
            TLegend* leg1 = new TLegend(0.60,0.15,0.80,0.35);
            TLine *line1 = new TLine(20,1,150,1);
            line1->SetLineStyle(kDashed);
            line1->SetLineWidth(1);
            leg1->SetFillColor(kWhite);
            leg1->SetFillStyle(1001);
            for(unsigned int u=0; u<effi_ele_eta_EGall.size(); u++){
                leg1->AddEntry(effi_ele_eta_EGall[u],(name_ele_eta_EGall[u]).c_str(),"L");
                effi_ele_eta_EGall[u]->SetLineWidth(1);
                effi_ele_eta_EGall[u]->SetLineColor(u+4);
                effi_ele_eta_EGall[u]->Draw("L same");
            }
            leg1->Draw("same");
            line1->Draw("same");
            c1->Print(("L1_EGall_" + type1 + "_" + type2 + "_" + type3 + "_ele_eta_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

        // Efficiency vs Electron Pt for Ele+Jet Trigger
        if(effi_ele_pt_elejet.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(20,0.30,150,1.03,("DATA " + type1 + " : L1_IsoEG30er_Jet34 OR all SingleEG Efficiency vs Electron pT; pT (Electron) [GeV] ; Efficiency").c_str());
			c1->SetGrid();
            TLegend* leg1 = new TLegend(0.60,0.15,0.80,0.35);
            TLine *line1 = new TLine(20,1,150,1);
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
            c1->Print(("L1_EleJet_" + type1 + "_" + type2 + "_" + type3 + "_ele_pt_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

        // Efficiency vs Electron Eta for Ele+Jet Trigger
        if(effi_ele_eta_elejet.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(-3,0.30,3,1.03,("DATA " + type1 + " : L1_IsoEG30er_Jet34 OR all SingleEG Efficiency vs Electron Eta; Eta (Electron) ; Efficiency").c_str());
            c1->SetGrid();
            TLegend* leg1 = new TLegend(0.60,0.15,0.80,0.35);
            TLine *line1 = new TLine(20,1,150,1);
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
            c1->Print(("L1_EleJet_" + type1 + "_" + type2 + "_" + type3 + "_ele_eta_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

		// Efficiency vs Electron Pt for Ele+HT Trigger
		if(effi_ele_pt_eleht.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(20,0.30,150,1.03,("DATA " + type1 + " : L1_IsoEG28er_HTT100 OR all SingleEG vs Electron pT; pT (Electron) [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1 = new TLegend(0.60,0.15,0.80,0.35);
			TLine *line1 = new TLine(20,1,150,1);
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
			c1->Print(("L1_EleHT_" + type1 + "_" + type2 + "_" + type3 + "_ele_pt_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

        // Efficiency vs Electron Eta for Ele+HT Trigger
        if(effi_ele_eta_eleht.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(-3,0.30,3,1.03,("DATA " + type1 + " : L1_IsoEG28er_HTT100 OR all SingleEG vs Electron Eta; Eta (Electron) ; Efficiency").c_str());
            c1->SetGrid();
            TLegend* leg1 = new TLegend(0.60,0.15,0.80,0.35);
            TLine *line1 = new TLine(20,1,150,1);
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
            c1->Print(("L1_EleHT_" + type1 + "_" + type2 + "_" + type3 + "_ele_eta_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

        // Efficiency vs Jet Pt for Ele+Jet Trigger
        if(effi_jet_pt_elejet.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(20,0.10,150,1.03,("DATA " + type1 + " : L1_IsoEG30er_Jet34 OR all SingleEG Efficiency vs Leading Jet pT; pT (Leading Jet) [GeV] ; Efficiency").c_str());
			c1->SetGrid();
            TLegend* leg1 = new TLegend(0.60,0.15,0.80,0.35);
            TLine *line1 = new TLine(20,1,150,1);
            line1->SetLineStyle(kDashed);
            line1->SetLineWidth(1);
            leg1->SetFillColor(kWhite);
            leg1->SetFillStyle(1001);
            for(unsigned int u=0; u<effi_jet_pt_elejet.size(); u++){
                leg1->AddEntry(effi_jet_pt_elejet[u],(name_jet_pt_elejet[u]).c_str(),"L");
                effi_jet_pt_elejet[u]->SetLineWidth(1);
                effi_jet_pt_elejet[u]->SetLineColor(u+4);
                effi_jet_pt_elejet[u]->Draw("L same");
            }
            leg1->Draw("same");
            line1->Draw("same");
            c1->Print(("L1_EleJet_" + type1 + "_" + type2 + "_" + type3 + "_jet_pt_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

        // Efficiency vs Jet Eta for Ele+Jet Trigger
        if(effi_jet_eta_elejet.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(-3,0.30,3,1.03,("DATA " + type1 + " : L1_IsoEG30er_Jet34 OR all SingleEG Efficiency vs Leading Jet Eta; Eta (Leading Jet) ; Efficiency").c_str());
            c1->SetGrid();
            TLegend* leg1 = new TLegend(0.60,0.15,0.80,0.35);
            TLine *line1 = new TLine(20,1,150,1);
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
            c1->Print(("L1_EleJet_" + type1 + "_" + type2 + "_" + type3 + "_jet_eta_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

        // Efficiency vs HT for Ele+HT Trigger
        if(effi_ht_eleht.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(0,0.20,500,1.03,("DATA " + type1 + " : L1_IsoEG28er_HTT100 OR all SingleEG Efficiency vs HT; HT [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1 = new TLegend(0.60,0.15,0.80,0.35);
            TLine *line1 = new TLine(0,1,500,1);
            line1->SetLineStyle(kDashed);
            line1->SetLineWidth(1);
            leg1->SetFillColor(kWhite);
            leg1->SetFillStyle(1001);
            for(unsigned int u=0; u<effi_ht_eleht.size(); u++){
                leg1->AddEntry(effi_ht_eleht[u],(name_ht_eleht[u]).c_str(),"L");
                effi_ht_eleht[u]->SetLineWidth(1);
                effi_ht_eleht[u]->SetLineColor(u+4);
                effi_ht_eleht[u]->Draw("L same");
            }
            leg1->Draw("same");
            line1->Draw("same");
            c1->Print(("L1_EleHT_" + type1 + "_" + type2 + "_" + type3 + "_ht_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

        // Efficiency vs Nr. of jets for Ele+Jet Trigger
        if(effi_njets_elejet.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(0,0.40,10,1.03,("DATA " + type1 + " : L1_IsoEG30er_Jet34 OR all SingleEG Efficiency vs Jet Multiplicity; Nr. of Jets ; Efficiency").c_str());
			c1->SetGrid();
            TLegend* leg1 = new TLegend(0.60,0.15,0.80,0.35);
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
            c1->Print(("L1_EleJet_" + type1 + "_" + type2 + "_" + type3 + "_njets_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

		// Efficiency vs Nr. of jets for Ele+HT Trigger
		if(effi_njets_eleht.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(0,0.40,10,1.03,("DATA " + type1 + " : L1_IsoEG28er_HTT100 OR all SingleEG Efficiency vs Jet Multiplicity; Nr. of Jets ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1 = new TLegend(0.60,0.15,0.80,0.35);
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
			c1->Print(("L1_EleHT_" + type1 + "_" + type2 + "_" + type3 + "_njets_turn_on.pdf").c_str());
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


