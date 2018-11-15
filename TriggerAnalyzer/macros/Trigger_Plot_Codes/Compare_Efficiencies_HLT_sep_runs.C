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
#include <sstream>
#include "TGraphAsymmErrors.h"
#include "TVector.h"
#include "TLorentzVector.h"
#include "TEfficiency.h"
#include "THStack.h"


void Compare_Efficiencies_HLT_sep_runs() {

    TH1::SetDefaultSumw2();
    TFile *file;

    // List of Files

    ifstream fin;
    fin.open("files_HLT_sep_runs.txt");
    char filenames[200][200];
    int nfiles = 0;

    while(!fin.eof()){
        fin>>filenames[nfiles];
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

        file = new TFile(filenames[i]);
        std::string fname = filenames[i];
        std::string type0, type1, type2, type3;

		if(fname.find("Global_MET") != std::string::npos){
			type0 = "DATA";
            type1 = "Global_MET";
		}
		else if(fname.find("Global_SM") != std::string::npos){
			type0 = "DATA";
			type1 = "Global_SM";
		}
		else if(fname.find("Global_MC_MET") != std::string::npos){
			type0 = "MC";
			type1 = "Global_MC_MET";
		}
		else if(fname.find("Global_MC_SM") != std::string::npos){
			type0 = "MC";
			type1 = "Global_MC_SM";
		}
		else if(fname.find("JetHTLeg_EG") != std::string::npos){
			type0 = "DATA";
			type1 = "JetHTLeg_EG";
		}
		else if(fname.find("JetHTLeg_MC_EG") != std::string::npos){
			type0 = "MC";
			type1 = "JetHTLeg_MC_EG";
		}
		else if(fname.find("EleLeg_JET") != std::string::npos){
			type0 = "DATA";
			type1 = "EleLeg_JET";
		}
		else if(fname.find("EleLeg_MC_JET") != std::string::npos){
			type0 = "MC";
			type1 = "EleLeg_MC_JET";
		}
		else if(fname.find("EleLeg_HT") != std::string::npos){
			type0 = "DATA";
			type1 = "EleLeg_HT";
		}
		else if(fname.find("EleLeg_MC_HT") != std::string::npos){
			type0 = "MC";
			type1 = "EleLeg_MC_HT";
		}
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

        if(fname.find("2018ABCD") != std::string::npos)
            type3 = "2018ABCD";
		else if(fname.find("2018A") != std::string::npos)
			type3 = "2018A";
		else if(fname.find("2018B") != std::string::npos)
			type3 = "2018B";
		else if(fname.find("2018C") != std::string::npos)
			type3 = "2018C";
		else if(fname.find("2018D") != std::string::npos)
			type3 = "2018D";
		else if(fname.find("ttjets") != std::string::npos)
			type3 = "ttjets";
		else
			type3 = "";

        std::vector<TEfficiency*> effi_ele_pt;
		std::vector<TEfficiency*> effi_ele_eta;
        std::vector<TEfficiency*> effi_jet_pt;
		std::vector<TEfficiency*> effi_jet_eta;
        std::vector<TEfficiency*> effi_ht;
        std::vector<TEfficiency*> effi_njets;
        std::vector<std::string> name_ele_pt;
		std::vector<std::string> name_ele_eta;
        std::vector<std::string> name_jet_pt;
		std::vector<std::string> name_jet_eta;
        std::vector<std::string> name_ht;
        std::vector<std::string> name_njets;

        std::string name = "";

        for(int j=0; j<nhistos; j++){

            std::string hname = histonames[j];

            TEfficiency *h = (TEfficiency*)file->Get((hname).c_str());
            h->SetTitle((hname).c_str());
            h->SetName((hname).c_str());

            if(hname.find("_ele_pt_") != std::string::npos){
                if( (hname.find("jet35") != std::string::npos) && (type1 != "EleLeg_HT") && (type1 != "EleLeg_MC_HT") ){
					name = "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned";
                    effi_ele_pt.push_back(h);
                    name_ele_pt.push_back((name).c_str());
                }
                else if( (hname.find("ht150") != std::string::npos) && (type1 != "EleLeg_JET") && (type1 != "EleLeg_MC_JET") ){
					name = "HLT_Ele28_eta2p1_WPTight_Gsf_HT150";
					effi_ele_pt.push_back(h);
                    name_ele_pt.push_back((name).c_str());
                }
				else if( (hname.find("ele32") != std::string::npos)){
					name = "HLT_Ele32_WPTight_Gsf";
					effi_ele_pt.push_back(h);
					name_ele_pt.push_back((name).c_str());
				}
				else if( (hname.find("ele35") != std::string::npos)){
					name = "HLT_Ele35_WPTight_Gsf";
					effi_ele_pt.push_back(h);
					name_ele_pt.push_back((name).c_str());
				}
            }

			if(hname.find("_ele_eta_") != std::string::npos){
				if( (hname.find("jet35") != std::string::npos) && (type1 != "EleLeg_HT") && (type1 != "EleLeg_MC_HT") ){
					name = "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned";
					effi_ele_eta.push_back(h);
					name_ele_eta.push_back((name).c_str());
				}
				else if( (hname.find("ht150") != std::string::npos) && (type1 != "EleLeg_JET") && (type1 != "EleLeg_MC_JET") ){
					name = "HLT_Ele28_eta2p1_WPTight_Gsf_HT150";
					effi_ele_eta.push_back(h);
					name_ele_eta.push_back((name).c_str());
				}
				else if( (hname.find("ele32") != std::string::npos)){
				 	name = "HLT_Ele32_WPTight_Gsf";
				 	effi_ele_eta.push_back(h);
				 	name_ele_eta.push_back((name).c_str());
                }
                else if( (hname.find("ele35") != std::string::npos)){
					name = "HLT_Ele35_WPTight_Gsf";
					effi_ele_eta.push_back(h);
					name_ele_eta.push_back((name).c_str());
				}
			}

			if(hname.find("_jet_pt_") != std::string::npos){
				if( (hname.find("jet35") != std::string::npos) && (type1 != "EleLeg_HT") && (type1 != "EleLeg_MC_HT") ){
					name = "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned";
					effi_jet_pt.push_back(h);
					name_jet_pt.push_back((name).c_str());
				}
			}

			if(hname.find("_jet_eta_") != std::string::npos){
				if( (hname.find("jet35") != std::string::npos) && (type1 != "EleLeg_HT") && (type1 != "EleLeg_MC_HT") ){
					name = "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned";
					effi_jet_eta.push_back(h);
					name_jet_eta.push_back((name).c_str());
				}
				else if( (hname.find("ht150") != std::string::npos) && (type1 != "EleLeg_JET") && (type1 != "EleLeg_MC_JET") ){
					name = "HLT_Ele28_eta2p1_WPTight_Gsf_HT150";
					effi_jet_eta.push_back(h);
					name_jet_eta.push_back((name).c_str());
				}
			}

			if(hname.find("_ht_") != std::string::npos){
				if( (hname.find("ht150") != std::string::npos) && (type1 != "EleLeg_JET") && (type1 != "EleLeg_MC_JET") ){
					name = "HLT_Ele28_eta2p1_WPTight_Gsf_HT150";
					effi_ht.push_back(h);
					name_ht.push_back((name).c_str());
				}
			}

            if(hname.find("_njets_") != std::string::npos){
				if( (hname.find("jet35") != std::string::npos) && (type1 != "EleLeg_HT") && (type1 != "EleLeg_MC_HT") ){
					name = "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned";
					effi_njets.push_back(h);
					name_njets.push_back((name).c_str());
				}
				else if( (hname.find("ht150") != std::string::npos) && (type1 != "EleLeg_JET") && (type1 != "EleLeg_MC_JET") ){
					name = "HLT_Ele28_eta2p1_WPTight_Gsf_HT150";
					effi_njets.push_back(h);
					name_njets.push_back((name).c_str());
				}
				else if( (hname.find("ele32") != std::string::npos)){
				 	name = "HLT_Ele32_WPTight_Gsf";
				 	effi_njets.push_back(h);
				 	name_njets.push_back((name).c_str());
                }
                else if( (hname.find("ele35") != std::string::npos)){
					name = "HLT_Ele35_WPTight_Gsf";
					effi_njets.push_back(h);
					name_njets.push_back((name).c_str());
				}
            }

        }

        gStyle->SetLegendBorderSize(1);
        gStyle->SetLegendTextSize(0.025);

        // Efficiency vs Electron Pt
        if(effi_ele_pt.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(20,0,400,1.03,(type0 + " " + type1 + "_" + type2 + "_" + type3 + " : HLT Efficiency vs Electron pT; pT (Electron) [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.3,0.15,0.85,0.35);
            TLine *line1 = new TLine(20,1,400,1);
            line1->SetLineStyle(kDashed);
            line1->SetLineWidth(1);
            leg1->SetFillColor(kWhite);
            leg1->SetFillStyle(1001);
            for(unsigned int u=0; u<effi_ele_pt.size(); u++){
                leg1->AddEntry(effi_ele_pt[u],(name_ele_pt[u]).c_str(),"L");
                effi_ele_pt[u]->SetLineWidth(1);
                effi_ele_pt[u]->SetLineColor(u+6);
                effi_ele_pt[u]->Draw("L same");
            }
            leg1->Draw("same");
            line1->Draw("same");
            c1->Print(("HLT_" + type1 + "_" + type2 + "_" + type3 + "_" + "ele_pt_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

		// Efficiency vs Electron Eta
		if(effi_ele_eta.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(-3,0,3,1.03,(type0 + " " + type1 + "_" + type2 + "_" + type3 + " : HLT Efficiency vs Electron eta; eta (Electron) ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.3,0.15,0.85,0.35);
			TLine *line1 = new TLine(-3,1,3,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_ele_eta.size(); u++){
				leg1->AddEntry(effi_ele_eta[u],(name_ele_eta[u]).c_str(),"L");
				effi_ele_eta[u]->SetLineWidth(1);
				effi_ele_eta[u]->SetLineColor(u+6);
				effi_ele_eta[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_" + type1 + "_" + type2 + "_" + type3 + "_" + "ele_eta_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

        // Efficiency vs Jet Pt
        if(effi_jet_pt.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(20,0,400,1.03,(type0 + " " + type1 + "_" + type2 + "_" + type3 + " : HLT Efficiency vs Leading Jet pT; pT (Leading Jet) [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.3,0.15,0.85,0.35);
            TLine *line1 = new TLine(20,1,400,1);
            line1->SetLineStyle(kDashed);
            line1->SetLineWidth(1);
            leg1->SetFillColor(kWhite);
            leg1->SetFillStyle(1001);
            for(unsigned int u=0; u<effi_jet_pt.size(); u++){
                leg1->AddEntry(effi_jet_pt[u],(name_jet_pt[u]).c_str(),"L");
                effi_jet_pt[u]->SetLineWidth(1);
                effi_jet_pt[u]->SetLineColor(u+2);
                effi_jet_pt[u]->Draw("L same");
            }
            leg1->Draw("same");
            line1->Draw("same");
            c1->Print(("HLT_" + type1 + "_" + type2 + "_" + type3 + "_" + "jet_pt_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

		// Efficiency vs Jet Eta
		if(effi_jet_eta.size() != 0){
			TCanvas *c1 = new TCanvas("c1","test",1100,650);
			c1->DrawFrame(-3,0,3,1.03,(type0 + " " + type1 + "_" + type2 + "_" + type3 + " : HLT Efficiency vs Leading Jet eta; eta (Leading Jet) ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.3,0.15,0.85,0.35);
			TLine *line1 = new TLine(-3,1,3,1);
			line1->SetLineStyle(kDashed);
			line1->SetLineWidth(1);
			leg1->SetFillColor(kWhite);
			leg1->SetFillStyle(1001);
			for(unsigned int u=0; u<effi_jet_eta.size(); u++){
				leg1->AddEntry(effi_jet_eta[u],(name_jet_eta[u]).c_str(),"L");
				effi_jet_eta[u]->SetLineWidth(1);
				effi_jet_eta[u]->SetLineColor(u+6);
				effi_jet_eta[u]->Draw("L same");
			}
			leg1->Draw("same");
			line1->Draw("same");
			c1->Print(("HLT_" + type1 + "_" + type2 + "_" + type3 + "_" + "jet_eta_turn_on.pdf").c_str());
			delete c1;
			delete leg1;
			delete line1;
		}

        // Efficiency vs HT
        if(effi_ht.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(0,0,1000,1.03,(type0 + " " + type1 + "_" + type2 + "_" + type3 + " : HLT Efficiency vs HT; HT [GeV] ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.3,0.15,0.85,0.35);
            TLine *line1 = new TLine(0,1,1000,1);
            line1->SetLineStyle(kDashed);
            line1->SetLineWidth(1);
            leg1->SetFillColor(kWhite);
            leg1->SetFillStyle(1001);
            for(unsigned int u=0; u<effi_ht.size(); u++){
                leg1->AddEntry(effi_ht[u],(name_ht[u]).c_str(),"L");
                effi_ht[u]->SetLineWidth(1);
                effi_ht[u]->SetLineColor(u+2);
                effi_ht[u]->Draw("L same");
            }
            leg1->Draw("same");
            line1->Draw("same");
            c1->Print(("HLT_" + type1 + "_" + type2 + "_" + type3 + "_" + "ht_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

        // Efficiency vs Nr. of jets
        if(effi_njets.size() != 0){
            TCanvas *c1 = new TCanvas("c1","test",1100,650);
            c1->DrawFrame(0,0,10,1.03,(type0 + " " + type1 + "_" + type2 + "_" + type3 + " : HLT Efficiency vs Jet Multiplicity; Nr. of Jets ; Efficiency").c_str());
			c1->SetGrid();
			TLegend* leg1;
			leg1 = new TLegend(0.3,0.15,0.85,0.35);
            TLine *line1 = new TLine(0,1,10,1);
            line1->SetLineStyle(kDashed);
            line1->SetLineWidth(1);
            leg1->SetFillColor(kWhite);
            leg1->SetFillStyle(1001);
            for(unsigned int u=0; u<effi_njets.size(); u++){
                leg1->AddEntry(effi_njets[u],(name_njets[u]).c_str(),"L");
                effi_njets[u]->SetLineWidth(1);
                effi_njets[u]->SetLineColor(u+6);
                effi_njets[u]->Draw("L same");
            }
            leg1->Draw("same");
            line1->Draw("same");
            c1->Print(("HLT_" + type1 + "_" + type2 + "_" + type3 + "_" + "njets_turn_on.pdf").c_str());
            delete c1;
            delete leg1;
            delete line1;
        }

    }

    delete file;
	return;
}


