#ifndef Trigger_data_analyzer_cc
#define Trigger_data_analyzer_cc


/// Includes
#include "Trigger_data_analyzer.h"

/// Constructor
Trigger_data_analyzer::Trigger_data_analyzer(const edm::ParameterSet &iConfig)
{

    hltTag = iConfig.getUntrackedParameter<string>("HLTsource");
    filterTag = iConfig.getUntrackedParameter<string>("PATsource");

    triggerResultsToken = consumes <edm::TriggerResults> (edm::InputTag(std::string("TriggerResults"), std::string(""), hltTag));
    triggerPrescalesToken = consumes <pat::PackedTriggerPrescales> (edm::InputTag(std::string("patTrigger")));
    filterResultsToken  = consumes <edm::TriggerResults> (edm::InputTag(std::string("TriggerResults"), std::string(""), filterTag));
    gtStage2DigisToken  = consumes <BXVector<GlobalAlgBlk>> (edm::InputTag(std::string("gtStage2Digis"), std::string(""), filterTag));
	triggerEventToken = consumes <trigger::TriggerEvent> (edm::InputTag(std::string("hltTriggerSummaryAOD"), std::string(""), hltTag));
	//triggerObjectsToken = consumes <pat::TriggerObjectStandAloneCollection> (edm::InputTag(std::string("selectedPatTrigger")));
	triggerObjectsToken = consumes <pat::TriggerObjectStandAloneCollection> (edm::InputTag(std::string("slimmedPatTrigger")));

    RhoToken = consumes<double>(edm::InputTag(std::string("fixedGridRhoFastjetAll")));
    vertexToken = consumes <vector<reco::Vertex>> (edm::InputTag(std::string("offlineSlimmedPrimaryVertices")));
    electronToken = consumes <vector<pat::Electron>> (edm::InputTag(std::string("slimmedElectrons")));
    muonToken = consumes <vector<pat::Muon>> (edm::InputTag(std::string("slimmedMuons")));
    jetToken = consumes <vector<pat::Jet>> (edm::InputTag(std::string("slimmedJets")));
    pfMetToken = consumes <vector<pat::MET>> (edm::InputTag(std::string("slimmedMETs")));

    m_ttree = fs_->make<TTree>("triggerTree", "triggerTree");
	eve=0;
	m_ttree->Branch("eve", "triggerStudyEventVars", &eve, 8000, 1);

}

/// Destructor
Trigger_data_analyzer::~Trigger_data_analyzer()
{

}

double Trigger_data_analyzer::get_ele_iso(const pat::Electron &iElectron, double & rho) {

    double pfIsoCharged = iElectron.pfIsolationVariables().sumChargedHadronPt;
    double pfIsoNeutral = iElectron.pfIsolationVariables().sumNeutralHadronEt + iElectron.pfIsolationVariables().sumPhotonEt;
    double Eta = fabs(iElectron.superCluster()->position().eta());
    double EffArea = 0;

    if (Eta >= 0. && Eta < 1.0) EffArea = 0.1566;
    else if (Eta >= 1.0 && Eta < 1.479) EffArea = 0.1626;
    else if (Eta >= 1.479 && Eta < 2.0) EffArea = 0.1073;
    else if (Eta >= 2.0 && Eta < 2.2) EffArea = 0.0854;
    else if (Eta >= 2.2 && Eta < 2.3) EffArea = 0.1051;
    else if (Eta >= 2.3 && Eta < 2.4) EffArea = 0.1204;
    else if (Eta >= 2.4 && Eta < 5.0) EffArea = 0.1524;


    double correction = rho*EffArea;
    double pfIsoPUSubtracted = std::max( 0.0, pfIsoNeutral - correction );
    return ((pfIsoCharged + pfIsoPUSubtracted)/iElectron.pt());

}


double Trigger_data_analyzer::get_mu_iso(const pat::Muon &iMuon) {

    double pfIsoCharged = iMuon.pfIsolationR04().sumChargedHadronPt;
    double pfIsoNeutral = iMuon.pfIsolationR04().sumNeutralHadronEt + iMuon.pfIsolationR04().sumPhotonEt;
    double correction =  0.5*iMuon.pfIsolationR04().sumPUPt;
    double pfIsoPUSubtracted = std::max( 0.0, pfIsoNeutral - correction );
    return ((pfIsoCharged + pfIsoPUSubtracted)/iMuon.pt());

}

int Trigger_data_analyzer::is_ele_tightid(const pat::Electron& iElectron, double & rho) {

    double SCeta = (iElectron.superCluster().isAvailable()) ? iElectron.superCluster()->position().eta() : -99;
    double absSCeta = fabs(SCeta);
    double sc_energy = iElectron.superCluster()->energy();
    double relIso = get_ele_iso(iElectron, rho);

    bool isEB = ( absSCeta <= 1.479 );

    double full5x5_sigmaIetaIeta = iElectron.full5x5_sigmaIetaIeta();
    double dEtaInSeed = iElectron.superCluster().isNonnull() && iElectron.superCluster()->seed().isNonnull() ? iElectron.deltaEtaSuperClusterTrackAtVtx() - iElectron.superCluster()->eta() + iElectron.superCluster()->seed()->eta() : std::numeric_limits<float>::max();
    double fabsdEtaInSeed=fabs(dEtaInSeed);
    double dPhiIn = fabs( iElectron.deltaPhiSuperClusterTrackAtVtx() );
    double hOverE = iElectron.hcalOverEcal();

    double ooEmooP = -999;
    if( iElectron.ecalEnergy() == 0 ) ooEmooP = 1e30;
    else if( !std::isfinite(iElectron.ecalEnergy()) ) ooEmooP = 1e30;
    else ooEmooP = fabs(1.0/iElectron.ecalEnergy() - iElectron.eSuperClusterOverP()/iElectron.ecalEnergy() );

    double d0 = -999;
    double dZ = -999;
    double expectedMissingInnerHits = 999;
    if( iElectron.gsfTrack().isAvailable() ){
        d0 = fabs(iElectron.gsfTrack()->dxy(vertex.position()));
        dZ = fabs(iElectron.gsfTrack()->dz(vertex.position()));
        expectedMissingInnerHits = iElectron.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
    }

    bool passConversionVeto = ( iElectron.passConversionVeto() );

    bool pass = false;
    if( isEB ){
                pass = ( full5x5_sigmaIetaIeta < 0.0104 &&
                        fabsdEtaInSeed < 0.00353 &&
                        dPhiIn < 0.0499 &&
                        hOverE < ( 0.026 + (0.5/sc_energy)+ (0.201*rho/sc_energy) ) &&
                        ooEmooP < 0.0278 &&
                        d0 < 0.05 &&
                        dZ < 0.1 &&
                        expectedMissingInnerHits <= 1 &&
                        passConversionVeto &&
                        relIso < 0.0361
                        );
    }
    else{
                pass = ( full5x5_sigmaIetaIeta < 0.0305 &&
                        fabsdEtaInSeed < 0.00567 &&
                        dPhiIn < 0.0165 &&
                        hOverE < ( 0.026 + (0.5/sc_energy)+ (0.201*rho/sc_energy) ) &&
                        ooEmooP < 0.0158 &&
                        d0 < 0.1 &&
                        dZ < 0.2 &&
                        expectedMissingInnerHits <= 1 &&
                        passConversionVeto &&
                        relIso < 0.094
                        );
    }

	if (pass == true)
		return 1;
	else
		return 0;
}

int Trigger_data_analyzer::is_ele_looseid(const pat::Electron& iElectron, double & rho) {

    double SCeta = (iElectron.superCluster().isAvailable()) ? iElectron.superCluster()->position().eta() : -99;
    double absSCeta = fabs(SCeta);
    double sc_energy = iElectron.superCluster()->energy();
    double relIso = get_ele_iso(iElectron, rho);

    bool isEB = ( absSCeta <= 1.479 );

    double full5x5_sigmaIetaIeta = iElectron.full5x5_sigmaIetaIeta();
    double dEtaInSeed = iElectron.superCluster().isNonnull() && iElectron.superCluster()->seed().isNonnull() ? iElectron.deltaEtaSuperClusterTrackAtVtx() - iElectron.superCluster()->eta() + iElectron.superCluster()->seed()->eta() : std::numeric_limits<float>::max();
    double fabsdEtaInSeed=fabs(dEtaInSeed);
    double dPhiIn = fabs( iElectron.deltaPhiSuperClusterTrackAtVtx() );
    double hOverE = iElectron.hcalOverEcal();

    double ooEmooP = -999;
    if( iElectron.ecalEnergy() == 0 ) ooEmooP = 1e30;
    else if( !std::isfinite(iElectron.ecalEnergy()) ) ooEmooP = 1e30;
    else ooEmooP = fabs(1.0/iElectron.ecalEnergy() - iElectron.eSuperClusterOverP()/iElectron.ecalEnergy() );

    double d0 = -999;
    double dZ = -999;
    double expectedMissingInnerHits = 999;
    if( iElectron.gsfTrack().isAvailable() ){
        d0 = fabs(iElectron.gsfTrack()->dxy(vertex.position()));
        dZ = fabs(iElectron.gsfTrack()->dz(vertex.position()));
        expectedMissingInnerHits = iElectron.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
    }

    bool passConversionVeto = ( iElectron.passConversionVeto() );

    bool pass = false;
    if( isEB ){
        pass = ( full5x5_sigmaIetaIeta < 0.0105 &&
                fabsdEtaInSeed < 0.00387 &&
                dPhiIn < 0.0716 &&
                hOverE < ( 0.05 + (0.5/sc_energy)+ (0.201*rho/sc_energy) ) &&
                ooEmooP < 0.129 &&
                d0 < 0.05 &&
                dZ < 0.1 &&
                expectedMissingInnerHits <= 1 &&
                passConversionVeto &&
                relIso < 0.133
                );
    }
    else{
        pass = ( full5x5_sigmaIetaIeta < 0.0356 &&
                fabsdEtaInSeed < 0.0072 &&
                dPhiIn < 0.147 &&
                hOverE < ( 0.0414 + (0.5/sc_energy)+ (0.201*rho/sc_energy) ) &&
                ooEmooP < 0.0875 &&
                d0 < 0.1 &&
                dZ < 0.2 &&
                expectedMissingInnerHits <= 1 &&
                passConversionVeto &&
                relIso < 0.146
                );
    }

    if (pass == true)
        return 1;
    else
        return 0;
}


int Trigger_data_analyzer::is_mu_tightid(const pat::Muon& iMuon){

    if( !iMuon.globalTrack().isAvailable() ) return 0;

    bool passesGlobalTrackID = ( (iMuon.globalTrack()->normalizedChi2() < 10.)
                                && (iMuon.globalTrack()->hitPattern().numberOfValidMuonHits() > 0)
                                );
    if(!passesGlobalTrackID) return 0;


    if(!iMuon.muonBestTrack().isAvailable() )return 0;
    bool passesMuonBestTrackID = ( (fabs(iMuon.muonBestTrack()->dxy(vertex.position())) < 0.2)
                                  && (fabs(iMuon.muonBestTrack()->dz(vertex.position())) < 0.5)
                                  );
    if(!passesMuonBestTrackID) return 0;

    if(!iMuon.innerTrack().isAvailable() ) return 0;
    bool passesInnerTrackID = (iMuon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0);
    if(!passesInnerTrackID) return 0;

    if(!iMuon.track().isAvailable() ) return 0;
    bool passesTrackID = (iMuon.track()->hitPattern().trackerLayersWithMeasurement() > 5);
    if(!passesTrackID) return 0;
    
    if(iMuon.numberOfMatchedStations() <= 1) return 0;
    
    if(!iMuon.isPFMuon()) return 0;
    
    return 1;
    
}

int Trigger_data_analyzer::is_jet_tightid(const pat::Jet &iJet) {

	bool tight = false;

    if( fabs(iJet.eta())<=2.7 ){
        tight = (
                 iJet.neutralHadronEnergyFraction() < 0.90 &&
                 //iJet.chargedEmEnergyFraction() < 0.90 &&
                 iJet.neutralEmEnergyFraction() < 0.90 &&
                 iJet.numberOfDaughters() > 1
                );

        if( fabs(iJet.eta())<=2.4 ){
            tight = ( tight &&
                     iJet.chargedHadronEnergyFraction() > 0.0 &&
                     iJet.chargedMultiplicity() > 0
                     );
        }
    }
    else if( fabs(iJet.eta())>2.7 && fabs(iJet.eta())<=3.0 ){
        tight = (
                 iJet.neutralEmEnergyFraction() < 0.99 &&
                 iJet.neutralEmEnergyFraction() > 0.02 &&
                 iJet.neutralMultiplicity() > 2
                 );

    }
    else if( fabs(iJet.eta())>3.0 ){
        tight = (
                 iJet.neutralHadronEnergyFraction() > 0.02 &&
                 iJet.neutralEmEnergyFraction() < 0.90 &&
                 iJet.neutralMultiplicity() > 10
                 );
    }

	if(tight == true)
		return 1;
	else
		return 0;
}

double Trigger_data_analyzer::get_jet_csv(const pat::Jet& jet, const std::string taggername ){

	double defaultFailure = -0.1;

	double bTagVal = jet.bDiscriminator(taggername);

	if(std::isnan(bTagVal))
		return defaultFailure;

	return bTagVal;
}

int Trigger_data_analyzer::Check_PV(reco::Vertex & vtx) {

    if (vtx.isFake() || vtx.ndof() < 4.0 || abs(vtx.z()) > 24.0 ||
        abs(vtx.position().Rho()) > 2.0)
        return 0;
    else
        return 1;
}

void Trigger_data_analyzer::ntuple_initialize(){

	eve->event_nr_ = -99;
	eve->run_nr_ = -99;
	eve->lumi_nr_ = -99;

	eve->rho_ = -99;

	eve->pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_ = -99;
	eve->pass_L1_LooseIsoEG28er2p1_HTT100er_ = -99;
	eve->pass_L1_Mu_all_ = -99;
	eve->pass_L1_EG_all_ = -99;
    eve->pass_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3_ = -99;
    eve->pass_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3_ = -99;
    eve->pass_L1_LooseIsoEG24er2p1_HTT100er_ = -99;
    eve->pass_L1_LooseIsoEG26er2p1_HTT100er_ = -99;
    eve->pass_L1_LooseIsoEG30er2p1_HTT100er_ = -99;
	eve->pass_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_ = -99;
	eve->pass_HLT_Ele28_eta2p1_WPTight_Gsf_HT150_ = -99;
	eve->pass_HLT_Ele27_WPTight_Gsf_ = -99;
    eve->pass_HLT_Ele32_WPTight_Gsf_ = -99;
    eve->pass_HLT_Ele32_WPTight_Gsf_L1DoubleEG_ = -99;
	eve->pass_HLT_Ele35_WPTight_Gsf_ = -99;
    eve->pass_HLT_Ele35_WPTight_Gsf_L1EGMT_ = -99;
	eve->pass_HLT_Ele38_WPTight_Gsf_ = -99;
	eve->pass_HLT_Ele40_WPTight_Gsf_ = -99;
    eve->pass_HLT_IsoMu30_ = -99;
	eve->pass_HLT_IsoMu27_ = -99;
	eve->pass_HLT_IsoMu24_eta2p1_ = -99;
	eve->pass_HLT_IsoMu24_ = -99;
	eve->pass_HLT_PFMET110_PFMHT110_IDTight_ = -99;
	eve->pass_HLT_PFMET120_PFMHT120_IDTight_ = -99;
	eve->pass_HLT_PFMET130_PFMHT130_IDTight_ = -99;
	eve->pass_HLT_PFMET140_PFMHT140_IDTight_ = -99;
	eve->pass_HLT_PFHT180_ = -99;
	eve->pass_HLT_PFHT250_ = -99;
	eve->pass_HLT_PFHT370_ = -99;
	eve->pass_HLT_PFHT430_ = -99;
	eve->pass_HLT_PFHT510_ = -99;
	eve->pass_HLT_PFHT590_ = -99;
	eve->pass_HLT_PFHT680_ = -99;
	eve->pass_HLT_PFHT780_ = -99;
	eve->pass_HLT_PFHT890_ = -99;
	eve->pass_HLT_PFHT1050_ = -99;
	eve->pass_HLT_PFJet40_ = -99;
	eve->pass_HLT_PFJet60_ = -99;
	eve->pass_HLT_PFJet80_ = -99;
	eve->pass_HLT_PFJet140_ = -99;
	eve->pass_HLT_PFJet200_ = -99;
	eve->pass_HLT_PFJet260_ = -99;
	eve->pass_HLT_PFJet320_ = -99;
	eve->pass_HLT_PFJet400_ = -99;
	eve->pass_HLT_PFJet450_ = -99;
	eve->pass_HLT_PFJet500_ = -99;
	eve->pass_HLT_PFJet550_ = -99;

    eve->prescale_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_ = -99;
    eve->prescale_HLT_Ele28_eta2p1_WPTight_Gsf_HT150_ = -99;
    eve->prescale_HLT_Ele27_WPTight_Gsf_ = -99;
    eve->prescale_HLT_Ele32_WPTight_Gsf_ = -99;
    eve->prescale_HLT_Ele32_WPTight_Gsf_L1DoubleEG_ = -99;
    eve->prescale_HLT_Ele35_WPTight_Gsf_ = -99;
    eve->prescale_HLT_Ele35_WPTight_Gsf_L1EGMT_ = -99;
    eve->prescale_HLT_Ele38_WPTight_Gsf_ = -99;
    eve->prescale_HLT_Ele40_WPTight_Gsf_ = -99;
    eve->prescale_HLT_IsoMu30_ = -99;
    eve->prescale_HLT_IsoMu27_ = -99;
    eve->prescale_HLT_IsoMu24_eta2p1_ = -99;
    eve->prescale_HLT_IsoMu24_ = -99;
    eve->prescale_HLT_PFMET110_PFMHT110_IDTight_ = -99;
    eve->prescale_HLT_PFMET120_PFMHT120_IDTight_ = -99;
    eve->prescale_HLT_PFMET130_PFMHT130_IDTight_ = -99;
    eve->prescale_HLT_PFMET140_PFMHT140_IDTight_ = -99;
    eve->prescale_HLT_PFHT180_ = -99;
    eve->prescale_HLT_PFHT250_ = -99;
    eve->prescale_HLT_PFHT370_ = -99;
    eve->prescale_HLT_PFHT430_ = -99;
    eve->prescale_HLT_PFHT510_ = -99;
    eve->prescale_HLT_PFHT590_ = -99;
    eve->prescale_HLT_PFHT680_ = -99;
    eve->prescale_HLT_PFHT780_ = -99;
    eve->prescale_HLT_PFHT890_ = -99;
    eve->prescale_HLT_PFHT1050_ = -99;
    eve->prescale_HLT_PFJet40_ = -99;
    eve->prescale_HLT_PFJet60_ = -99;
    eve->prescale_HLT_PFJet80_ = -99;
    eve->prescale_HLT_PFJet140_ = -99;
    eve->prescale_HLT_PFJet200_ = -99;
    eve->prescale_HLT_PFJet260_ = -99;
    eve->prescale_HLT_PFJet320_ = -99;
    eve->prescale_HLT_PFJet400_ = -99;
    eve->prescale_HLT_PFJet450_ = -99;
    eve->prescale_HLT_PFJet500_ = -99;
    eve->prescale_HLT_PFJet550_ = -99;


	eve->pt_trigger_object_.clear();
	eve->eta_trigger_object_.clear();
	eve->phi_trigger_object_.clear();

	eve->filter_trigger_object_.clear();

	eve->pfMET_pt_ = -99;
	eve->pfMET_phi_ = -99;

	eve->jet_pt_.clear();
	eve->jet_eta_.clear();
	eve->jet_phi_.clear();
	eve->jet_energy_.clear();
	eve->jet_csv_.clear();
	eve->jet_is_ID_Tight_.clear();

	eve->ele_pt_.clear();
	eve->ele_eta_.clear();
	eve->ele_sceta_.clear();
	eve->ele_phi_.clear();
	eve->ele_energy_.clear();
	eve->ele_iso_.clear();
	eve->ele_is_ID_Tight_.clear();
    eve->ele_is_ID_Loose_.clear();

	eve->mu_pt_.clear();
	eve->mu_eta_.clear();
	eve->mu_phi_.clear();
	eve->mu_energy_.clear();
	eve->mu_iso_.clear();
	eve->mu_is_ID_Tight_.clear();

	return;
}


// ------------ method called for each event  ------------
void Trigger_data_analyzer::analyze(const edm::Event &iEvent,
                         const edm::EventSetup &iSetup)
{

    //std::cout<<"Step 1 \n";
    ntuple_initialize();

    int run_nr = iEvent.id().run();
    int event_nr = iEvent.id().event();
    int lumi_nr = iEvent.id().luminosityBlock();

    eve->run_nr_ = run_nr;
    eve->event_nr_ = event_nr;
    eve->lumi_nr_ = lumi_nr;

    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByToken(triggerResultsToken, triggerResults);

    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    iEvent.getByToken(triggerPrescalesToken, triggerPrescales);

    edm::Handle<BXVector<GlobalAlgBlk>> gtStage2Digis_;
    iEvent.getByToken(gtStage2DigisToken, gtStage2Digis_);
    BXVector<GlobalAlgBlk> gtStage2Digis = *gtStage2Digis_;

	edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
	iEvent.getByToken(triggerObjectsToken, triggerObjects);

    //std::cout<<"Step 2 \n";

	// L1 pass info

    int pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3 = -1;
    int pass_L1_LooseIsoEG28er2p1_HTT100er = -1;
    int pass_L1_Mu_all = -1;
	int pass_L1_EG_all = -1;
    int pass_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3 = -1;
    int pass_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3 = -1;
    int pass_L1_LooseIsoEG24er2p1_HTT100er = -1;
    int pass_L1_LooseIsoEG26er2p1_HTT100er = -1;
    int pass_L1_LooseIsoEG30er2p1_HTT100er = -1;

    if (run_nr < 316058){

		pass_L1_Mu_all = 0;
		for(int i=0; i<=23; i++){
            if ((i<=11) || (i>=13 && i<=15) || (i>=17)){
                if(gtStage2Digis.at(0,0).getAlgoDecisionFinal(i)==1){
                    pass_L1_Mu_all = 1;
                    break;
                }
            }
		}

        pass_L1_EG_all = 0;
        for(int i=55; i<=90; i++){
            if ((i==55) || (i>=59 && i<=68) || (i>=73 && i<=77) || (i>=86)){
                if(gtStage2Digis.at(0,0).getAlgoDecisionFinal(i)==1){
                    pass_L1_EG_all = 1;
                    break;
                }
            }
        }

        pass_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3 = gtStage2Digis.at(0,0).getAlgoDecisionFinal(226);
        pass_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3 = gtStage2Digis.at(0,0).getAlgoDecisionFinal(227);
        pass_L1_LooseIsoEG24er2p1_HTT100er = gtStage2Digis.at(0,0).getAlgoDecisionFinal(229);
        pass_L1_LooseIsoEG26er2p1_HTT100er = gtStage2Digis.at(0,0).getAlgoDecisionFinal(230);
        pass_L1_LooseIsoEG30er2p1_HTT100er = 0;

		pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3 = gtStage2Digis.at(0,0).getAlgoDecisionFinal(228);
		pass_L1_LooseIsoEG28er2p1_HTT100er = gtStage2Digis.at(0,0).getAlgoDecisionFinal(231);
	}
	else if (run_nr >= 316058 && run_nr < 319913 && run_nr != 319506 && run_nr != 319523 && run_nr != 319527 && run_nr != 319528 && run_nr != 319606 && run_nr != 319623 && run_nr != 319625 && run_nr != 319800 && run_nr != 319856){

        pass_L1_Mu_all = 0;
        for(int i=0; i<=33; i++){
            if ((i<=23) || (i>=25)){
                if(gtStage2Digis.at(0,0).getAlgoDecisionFinal(i)==1){
                    pass_L1_Mu_all = 1;
                    break;
                }
            }
        }

        pass_L1_EG_all = 0;
        for(int i=166; i<=195; i++){
            if ((i<=173) || (i==183) || (i==185) || (i==186) || (i==188) || (i==189) || (i>=191)){
                if(gtStage2Digis.at(0,0).getAlgoDecisionFinal(i)==1){
                    pass_L1_EG_all = 1;
                    break;
                }
            }
        }

        pass_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3 = gtStage2Digis.at(0,0).getAlgoDecisionFinal(234);
        pass_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3 = gtStage2Digis.at(0,0).getAlgoDecisionFinal(235);
        pass_L1_LooseIsoEG24er2p1_HTT100er = gtStage2Digis.at(0,0).getAlgoDecisionFinal(238);
        pass_L1_LooseIsoEG26er2p1_HTT100er = gtStage2Digis.at(0,0).getAlgoDecisionFinal(239);
        pass_L1_LooseIsoEG30er2p1_HTT100er = gtStage2Digis.at(0,0).getAlgoDecisionFinal(241);

        pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3 = gtStage2Digis.at(0,0).getAlgoDecisionFinal(236);
        pass_L1_LooseIsoEG28er2p1_HTT100er = gtStage2Digis.at(0,0).getAlgoDecisionFinal(240);
	}
    else {

        pass_L1_Mu_all = 0;
        for(int i=0; i<=33; i++){
            if ((i<=23) || (i>=25)){
                if(gtStage2Digis.at(0,0).getAlgoDecisionFinal(i)==1){
                    pass_L1_Mu_all = 1;
                    break;
                }
            }
        }

        pass_L1_EG_all = 0;
        for(int i=162; i<=196; i++){
            if ((i==162) || (i==167) || (i==168) || (i==169) || (i==170) || (i==171) || (i==172) || (i==173) || (i==183) || (i==185) || (i==186) || (i==189) || (i==190) || (i>=192)){
                if(gtStage2Digis.at(0,0).getAlgoDecisionFinal(i)==1){
                    pass_L1_EG_all = 1;
                    break;
                }
            }
        }

        pass_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3 = gtStage2Digis.at(0,0).getAlgoDecisionFinal(234);
        pass_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3 = gtStage2Digis.at(0,0).getAlgoDecisionFinal(235);
        pass_L1_LooseIsoEG24er2p1_HTT100er = gtStage2Digis.at(0,0).getAlgoDecisionFinal(238);
        pass_L1_LooseIsoEG26er2p1_HTT100er = gtStage2Digis.at(0,0).getAlgoDecisionFinal(239);
        pass_L1_LooseIsoEG30er2p1_HTT100er = gtStage2Digis.at(0,0).getAlgoDecisionFinal(241);

        pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3 = gtStage2Digis.at(0,0).getAlgoDecisionFinal(236);
        pass_L1_LooseIsoEG28er2p1_HTT100er = gtStage2Digis.at(0,0).getAlgoDecisionFinal(240);
    }

	//std::cout<<"Step 3 \n";

	eve->pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_ = pass_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3;
	eve->pass_L1_LooseIsoEG28er2p1_HTT100er_ = pass_L1_LooseIsoEG28er2p1_HTT100er;
    eve->pass_L1_Mu_all_ = pass_L1_Mu_all;
    eve->pass_L1_EG_all_ = pass_L1_EG_all;
    eve->pass_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3_ = pass_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3;
    eve->pass_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3_ = pass_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3;
    eve->pass_L1_LooseIsoEG24er2p1_HTT100er_ = pass_L1_LooseIsoEG24er2p1_HTT100er;
    eve->pass_L1_LooseIsoEG26er2p1_HTT100er_ = pass_L1_LooseIsoEG26er2p1_HTT100er;
    eve->pass_L1_LooseIsoEG30er2p1_HTT100er_ = pass_L1_LooseIsoEG30er2p1_HTT100er;
    
	// HLT pass info

	int pass_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned = -1;
	int pass_HLT_Ele28_eta2p1_WPTight_Gsf_HT150 = -1;
	int pass_HLT_Ele27_WPTight_Gsf = -1;
    int pass_HLT_Ele32_WPTight_Gsf = -1;
    int pass_HLT_Ele32_WPTight_Gsf_L1DoubleEG = -1;
	int pass_HLT_Ele35_WPTight_Gsf = -1;
    int pass_HLT_Ele35_WPTight_Gsf_L1EGMT = -1;
	int pass_HLT_Ele38_WPTight_Gsf = -1;
	int pass_HLT_Ele40_WPTight_Gsf = -1;
    int pass_HLT_IsoMu30 = -1;
	int pass_HLT_IsoMu27 = -1;
	int pass_HLT_IsoMu24_eta2p1 = -1;
	int pass_HLT_IsoMu24 = -1;
	int pass_HLT_PFMET110_PFMHT110_IDTight = -1;
	int pass_HLT_PFMET120_PFMHT120_IDTight = -1;
	int pass_HLT_PFMET130_PFMHT130_IDTight = -1;
	int pass_HLT_PFMET140_PFMHT140_IDTight = -1;
	int pass_HLT_PFHT180 = -1;
	int pass_HLT_PFHT250 = -1;
	int pass_HLT_PFHT370 = -1;
	int pass_HLT_PFHT430 = -1;
	int pass_HLT_PFHT510 = -1;
	int pass_HLT_PFHT590 = -1;
	int pass_HLT_PFHT680 = -1;
	int pass_HLT_PFHT780 = -1;
	int pass_HLT_PFHT890 = -1;
	int pass_HLT_PFHT1050 = -1;
	int pass_HLT_PFJet40 = -1;
	int pass_HLT_PFJet60 = -1;
	int pass_HLT_PFJet80 = -1;
	int pass_HLT_PFJet140 = -1;
	int pass_HLT_PFJet200 = -1;
	int pass_HLT_PFJet260 = -1;
	int pass_HLT_PFJet320 = -1;
	int pass_HLT_PFJet400 = -1;
	int pass_HLT_PFJet450 = -1;
	int pass_HLT_PFJet500 = -1;
	int pass_HLT_PFJet550 = -1;

    int prescale_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned = -1;
    int prescale_HLT_Ele28_eta2p1_WPTight_Gsf_HT150 = -1;
    int prescale_HLT_Ele27_WPTight_Gsf = -1;
    int prescale_HLT_Ele32_WPTight_Gsf = -1;
    int prescale_HLT_Ele32_WPTight_Gsf_L1DoubleEG = -1;
    int prescale_HLT_Ele35_WPTight_Gsf = -1;
    int prescale_HLT_Ele35_WPTight_Gsf_L1EGMT = -1;
    int prescale_HLT_Ele38_WPTight_Gsf = -1;
    int prescale_HLT_Ele40_WPTight_Gsf = -1;
    int prescale_HLT_IsoMu30 = -1;
    int prescale_HLT_IsoMu27 = -1;
    int prescale_HLT_IsoMu24_eta2p1 = -1;
    int prescale_HLT_IsoMu24 = -1;
    int prescale_HLT_PFMET110_PFMHT110_IDTight = -1;
    int prescale_HLT_PFMET120_PFMHT120_IDTight = -1;
    int prescale_HLT_PFMET130_PFMHT130_IDTight = -1;
    int prescale_HLT_PFMET140_PFMHT140_IDTight = -1;
    int prescale_HLT_PFHT180 = -1;
    int prescale_HLT_PFHT250 = -1;
    int prescale_HLT_PFHT370 = -1;
    int prescale_HLT_PFHT430 = -1;
    int prescale_HLT_PFHT510 = -1;
    int prescale_HLT_PFHT590 = -1;
    int prescale_HLT_PFHT680 = -1;
    int prescale_HLT_PFHT780 = -1;
    int prescale_HLT_PFHT890 = -1;
    int prescale_HLT_PFHT1050 = -1;
    int prescale_HLT_PFJet40 = -1;
    int prescale_HLT_PFJet60 = -1;
    int prescale_HLT_PFJet80 = -1;
    int prescale_HLT_PFJet140 = -1;
    int prescale_HLT_PFJet200 = -1;
    int prescale_HLT_PFJet260 = -1;
    int prescale_HLT_PFJet320 = -1;
    int prescale_HLT_PFJet400 = -1;
    int prescale_HLT_PFJet450 = -1;
    int prescale_HLT_PFJet500 = -1;
    int prescale_HLT_PFJet550 = -1;

	if( triggerResults.isValid() ){

		std::vector<std::string> triggerNames = hlt_config_.triggerNames();

		for( unsigned int iPath=0; iPath<triggerNames.size(); iPath++ ){

			std::string pathName = triggerNames[iPath];
			unsigned int hltIndex = hlt_config_.triggerIndex(pathName);
            //int prescale = hlt_config_.prescaleValue(iEvent, iSetup, pathName);
            int prescale = triggerPrescales->getPrescaleForIndex(hltIndex);

			if( hltIndex >= triggerResults->size() ) continue;

			int accept = triggerResults->accept(hltIndex);

			std::string pathNameNoVer = hlt_config_.removeVersion(pathName);

			//if( accept ) hlt_cppath_[pathNameNoVer]+=1;

            if( pathName.find("HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v")!=std::string::npos ){
				pass_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned = (accept) ? 1 : 0;
                prescale_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned = prescale;
            }

            if( pathName.find("HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v")!=std::string::npos ){
				pass_HLT_Ele28_eta2p1_WPTight_Gsf_HT150 = (accept) ? 1 : 0;
                prescale_HLT_Ele28_eta2p1_WPTight_Gsf_HT150 = prescale;
            }

			if( pathName.find("HLT_Ele27_WPTight_Gsf_v")!=std::string::npos ){
				pass_HLT_Ele27_WPTight_Gsf = (accept) ? 1 : 0;
                prescale_HLT_Ele27_WPTight_Gsf = prescale;
            }

            if( pathName.find("HLT_Ele32_WPTight_Gsf_v")!=std::string::npos ){
                pass_HLT_Ele32_WPTight_Gsf = (accept) ? 1 : 0;
                prescale_HLT_Ele32_WPTight_Gsf = prescale;
            }

            if( pathName.find("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v")!=std::string::npos ){
                pass_HLT_Ele32_WPTight_Gsf_L1DoubleEG = (accept) ? 1 : 0;
                prescale_HLT_Ele32_WPTight_Gsf_L1DoubleEG = prescale;
            }

			if( pathName.find("HLT_Ele35_WPTight_Gsf_v")!=std::string::npos ){
				pass_HLT_Ele35_WPTight_Gsf = (accept) ? 1 : 0;
                prescale_HLT_Ele35_WPTight_Gsf = prescale;
            }

            if( pathName.find("HLT_Ele35_WPTight_Gsf_L1EGMT_v")!=std::string::npos ){
                pass_HLT_Ele35_WPTight_Gsf_L1EGMT = (accept) ? 1 : 0;
                prescale_HLT_Ele35_WPTight_Gsf_L1EGMT = prescale;
            }

			if( pathName.find("HLT_Ele38_WPTight_Gsf_v")!=std::string::npos ){
				pass_HLT_Ele38_WPTight_Gsf = (accept) ? 1 : 0;
                prescale_HLT_Ele38_WPTight_Gsf = prescale;
            }

			if( pathName.find("HLT_Ele40_WPTight_Gsf_v")!=std::string::npos ){
				pass_HLT_Ele40_WPTight_Gsf = (accept) ? 1 : 0;
                prescale_HLT_Ele40_WPTight_Gsf = prescale;
            }

			if( pathName.find("HLT_IsoMu30_v")!=std::string::npos ){
				pass_HLT_IsoMu30 = (accept) ? 1 : 0;
                prescale_HLT_IsoMu30 = prescale;
            }

            if( pathName.find("HLT_IsoMu27_v")!=std::string::npos ){
                pass_HLT_IsoMu27 = (accept) ? 1 : 0;
                prescale_HLT_IsoMu27 = prescale;
            }

			if( pathName.find("HLT_IsoMu24_eta2p1_v")!=std::string::npos ){
				pass_HLT_IsoMu24_eta2p1 = (accept) ? 1 : 0;
                prescale_HLT_IsoMu24_eta2p1 = prescale;
            }

			if( pathName.find("HLT_IsoMu24_v")!=std::string::npos ){
				pass_HLT_IsoMu24 = (accept) ? 1 : 0;
                prescale_HLT_IsoMu24 = prescale;
            }

			if( pathName.find("HLT_PFMET110_PFMHT110_IDTight_v")!=std::string::npos ){
				pass_HLT_PFMET110_PFMHT110_IDTight = (accept) ? 1 : 0;
                prescale_HLT_PFMET110_PFMHT110_IDTight = prescale;
            }

			if( pathName.find("HLT_PFMET120_PFMHT120_IDTight_v")!=std::string::npos ){
				pass_HLT_PFMET120_PFMHT120_IDTight = (accept) ? 1 : 0;
                prescale_HLT_PFMET120_PFMHT120_IDTight = prescale;
            }

			if( pathName.find("HLT_PFMET130_PFMHT130_IDTight_v")!=std::string::npos ){
				pass_HLT_PFMET130_PFMHT130_IDTight = (accept) ? 1 : 0;
                prescale_HLT_PFMET130_PFMHT130_IDTight = prescale;
            }

			if( pathName.find("HLT_PFMET140_PFMHT140_IDTight_v")!=std::string::npos ){
				pass_HLT_PFMET140_PFMHT140_IDTight = (accept) ? 1 : 0;
                prescale_HLT_PFMET140_PFMHT140_IDTight = prescale;
            }

			if( pathName.find("HLT_PFHT180_v")!=std::string::npos ){
				pass_HLT_PFHT180 = (accept) ? 1 : 0;
                prescale_HLT_PFHT180 = prescale;
            }

			if( pathName.find("HLT_PFHT250_v")!=std::string::npos ){
				pass_HLT_PFHT250 = (accept) ? 1 : 0;
                prescale_HLT_PFHT250 = prescale;
            }

			if( pathName.find("HLT_PFHT370_v")!=std::string::npos ){
				pass_HLT_PFHT370 = (accept) ? 1 : 0;
                prescale_HLT_PFHT370 = prescale;
            }

			if( pathName.find("HLT_PFHT430_v")!=std::string::npos ){
				pass_HLT_PFHT430 = (accept) ? 1 : 0;
                prescale_HLT_PFHT430 = prescale;
            }

			if( pathName.find("HLT_PFHT510_v")!=std::string::npos ){
				pass_HLT_PFHT510 = (accept) ? 1 : 0;
                prescale_HLT_PFHT510 = prescale;
            }

			if( pathName.find("HLT_PFHT590_v")!=std::string::npos ){
				pass_HLT_PFHT590 = (accept) ? 1 : 0;
                prescale_HLT_PFHT590 = prescale;
            }

			if( pathName.find("HLT_PFHT680_v")!=std::string::npos ){
				pass_HLT_PFHT680 = (accept) ? 1 : 0;
                prescale_HLT_PFHT680 = prescale;
            }

			if( pathName.find("HLT_PFHT780_v")!=std::string::npos ){
				pass_HLT_PFHT780 = (accept) ? 1 : 0;
                prescale_HLT_PFHT780 = prescale;
            }

			if( pathName.find("HLT_PFHT890_v")!=std::string::npos ){
				pass_HLT_PFHT890 = (accept) ? 1 : 0;
                prescale_HLT_PFHT890 = prescale;
            }

			if( pathName.find("HLT_PFHT1050_v")!=std::string::npos ){
				pass_HLT_PFHT1050 = (accept) ? 1 : 0;
                prescale_HLT_PFHT1050 = prescale;
            }

			if( pathName.find("HLT_PFJet40_v")!=std::string::npos ){
				pass_HLT_PFJet40 = (accept) ? 1 : 0;
                prescale_HLT_PFJet40 = prescale;
            }

			if( pathName.find("HLT_PFJet60_v")!=std::string::npos ){
				pass_HLT_PFJet60 = (accept) ? 1 : 0;
                prescale_HLT_PFJet60 = prescale;
            }

			if( pathName.find("HLT_PFJet80_v")!=std::string::npos ){
				pass_HLT_PFJet80 = (accept) ? 1 : 0;
                prescale_HLT_PFJet80 = prescale;
            }

			if( pathName.find("HLT_PFJet140_v")!=std::string::npos ){
				pass_HLT_PFJet140 = (accept) ? 1 : 0;
                prescale_HLT_PFJet140 = prescale;
            }

			if( pathName.find("HLT_PFJet200_v")!=std::string::npos ){
				pass_HLT_PFJet200 = (accept) ? 1 : 0;
                prescale_HLT_PFJet200 = prescale;
            }

			if( pathName.find("HLT_PFJet260_v")!=std::string::npos ){
				pass_HLT_PFJet260 = (accept) ? 1 : 0;
                prescale_HLT_PFJet260 = prescale;
            }

			if( pathName.find("HLT_PFJet320_v")!=std::string::npos ){
				pass_HLT_PFJet320 = (accept) ? 1 : 0;
                prescale_HLT_PFJet320 = prescale;
            }

			if( pathName.find("HLT_PFJet400_v")!=std::string::npos ){
				pass_HLT_PFJet400 = (accept) ? 1 : 0;
                prescale_HLT_PFJet400 = prescale;
            }

			if( pathName.find("HLT_PFJet450_v")!=std::string::npos ){
				pass_HLT_PFJet450 = (accept) ? 1 : 0;
                prescale_HLT_PFJet450 = prescale;
            }

            if( pathName.find("HLT_PFJet500_v")!=std::string::npos ){
				pass_HLT_PFJet500 = (accept) ? 1 : 0;
                prescale_HLT_PFJet500 = prescale;
            }

            if( pathName.find("HLT_PFJet550_v")!=std::string::npos ){
				pass_HLT_PFJet550 = (accept) ? 1 : 0;
                prescale_HLT_PFJet550 = prescale;
            }
		}
	}

	eve->pass_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_ = pass_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned;
	eve->pass_HLT_Ele28_eta2p1_WPTight_Gsf_HT150_ = pass_HLT_Ele28_eta2p1_WPTight_Gsf_HT150;
	eve->pass_HLT_Ele27_WPTight_Gsf_ = pass_HLT_Ele27_WPTight_Gsf;
    eve->pass_HLT_Ele32_WPTight_Gsf_ = pass_HLT_Ele32_WPTight_Gsf;
    eve->pass_HLT_Ele32_WPTight_Gsf_L1DoubleEG_ = pass_HLT_Ele32_WPTight_Gsf_L1DoubleEG;
	eve->pass_HLT_Ele35_WPTight_Gsf_ = pass_HLT_Ele35_WPTight_Gsf;
    eve->pass_HLT_Ele35_WPTight_Gsf_L1EGMT_ = pass_HLT_Ele35_WPTight_Gsf_L1EGMT;
	eve->pass_HLT_Ele38_WPTight_Gsf_ = pass_HLT_Ele38_WPTight_Gsf;
	eve->pass_HLT_Ele40_WPTight_Gsf_ = pass_HLT_Ele40_WPTight_Gsf;
	eve->pass_HLT_IsoMu30_ = pass_HLT_IsoMu30;
    eve->pass_HLT_IsoMu27_ = pass_HLT_IsoMu27;
	eve->pass_HLT_IsoMu24_eta2p1_ = pass_HLT_IsoMu24_eta2p1;
	eve->pass_HLT_IsoMu24_ = pass_HLT_IsoMu24;
	eve->pass_HLT_PFMET110_PFMHT110_IDTight_ = pass_HLT_PFMET110_PFMHT110_IDTight;
	eve->pass_HLT_PFMET120_PFMHT120_IDTight_ = pass_HLT_PFMET120_PFMHT120_IDTight;
	eve->pass_HLT_PFMET130_PFMHT130_IDTight_ = pass_HLT_PFMET130_PFMHT130_IDTight;
	eve->pass_HLT_PFMET140_PFMHT140_IDTight_ = pass_HLT_PFMET140_PFMHT140_IDTight;
	eve->pass_HLT_PFHT180_ = pass_HLT_PFHT180;
	eve->pass_HLT_PFHT250_ = pass_HLT_PFHT250;
	eve->pass_HLT_PFHT370_ = pass_HLT_PFHT370;
	eve->pass_HLT_PFHT430_ = pass_HLT_PFHT430;
	eve->pass_HLT_PFHT510_ = pass_HLT_PFHT510;
	eve->pass_HLT_PFHT590_ = pass_HLT_PFHT590;
	eve->pass_HLT_PFHT680_ = pass_HLT_PFHT680;
	eve->pass_HLT_PFHT780_ = pass_HLT_PFHT780;
	eve->pass_HLT_PFHT890_ = pass_HLT_PFHT890;
	eve->pass_HLT_PFHT1050_ = pass_HLT_PFHT1050;
	eve->pass_HLT_PFJet40_ = pass_HLT_PFJet40;
	eve->pass_HLT_PFJet60_ = pass_HLT_PFJet60;
	eve->pass_HLT_PFJet80_ = pass_HLT_PFJet80;
	eve->pass_HLT_PFJet140_ = pass_HLT_PFJet140;
	eve->pass_HLT_PFJet200_ = pass_HLT_PFJet200;
	eve->pass_HLT_PFJet260_ = pass_HLT_PFJet260;
	eve->pass_HLT_PFJet320_ = pass_HLT_PFJet320;
	eve->pass_HLT_PFJet400_ = pass_HLT_PFJet400;
	eve->pass_HLT_PFJet450_ = pass_HLT_PFJet450;
	eve->pass_HLT_PFJet500_ = pass_HLT_PFJet500;
	eve->pass_HLT_PFJet550_ = pass_HLT_PFJet550;

    eve->prescale_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_ = prescale_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned;
    eve->prescale_HLT_Ele28_eta2p1_WPTight_Gsf_HT150_ = prescale_HLT_Ele28_eta2p1_WPTight_Gsf_HT150;
    eve->prescale_HLT_Ele27_WPTight_Gsf_ = prescale_HLT_Ele27_WPTight_Gsf;
    eve->prescale_HLT_Ele32_WPTight_Gsf_ = prescale_HLT_Ele32_WPTight_Gsf;
    eve->prescale_HLT_Ele32_WPTight_Gsf_L1DoubleEG_ = prescale_HLT_Ele32_WPTight_Gsf_L1DoubleEG;
    eve->prescale_HLT_Ele35_WPTight_Gsf_ = prescale_HLT_Ele35_WPTight_Gsf;
    eve->prescale_HLT_Ele35_WPTight_Gsf_L1EGMT_ = prescale_HLT_Ele35_WPTight_Gsf_L1EGMT;
    eve->prescale_HLT_Ele38_WPTight_Gsf_ = prescale_HLT_Ele38_WPTight_Gsf;
    eve->prescale_HLT_Ele40_WPTight_Gsf_ = prescale_HLT_Ele40_WPTight_Gsf;
    eve->prescale_HLT_IsoMu30_ = prescale_HLT_IsoMu30;
    eve->prescale_HLT_IsoMu27_ = prescale_HLT_IsoMu27;
    eve->prescale_HLT_IsoMu24_eta2p1_ = prescale_HLT_IsoMu24_eta2p1;
    eve->prescale_HLT_IsoMu24_ = prescale_HLT_IsoMu24;
    eve->prescale_HLT_PFMET110_PFMHT110_IDTight_ = prescale_HLT_PFMET110_PFMHT110_IDTight;
    eve->prescale_HLT_PFMET120_PFMHT120_IDTight_ = prescale_HLT_PFMET120_PFMHT120_IDTight;
    eve->prescale_HLT_PFMET130_PFMHT130_IDTight_ = prescale_HLT_PFMET130_PFMHT130_IDTight;
    eve->prescale_HLT_PFMET140_PFMHT140_IDTight_ = prescale_HLT_PFMET140_PFMHT140_IDTight;
    eve->prescale_HLT_PFHT180_ = prescale_HLT_PFHT180;
    eve->prescale_HLT_PFHT250_ = prescale_HLT_PFHT250;
    eve->prescale_HLT_PFHT370_ = prescale_HLT_PFHT370;
    eve->prescale_HLT_PFHT430_ = prescale_HLT_PFHT430;
    eve->prescale_HLT_PFHT510_ = prescale_HLT_PFHT510;
    eve->prescale_HLT_PFHT590_ = prescale_HLT_PFHT590;
    eve->prescale_HLT_PFHT680_ = prescale_HLT_PFHT680;
    eve->prescale_HLT_PFHT780_ = prescale_HLT_PFHT780;
    eve->prescale_HLT_PFHT890_ = prescale_HLT_PFHT890;
    eve->prescale_HLT_PFHT1050_ = prescale_HLT_PFHT1050;
    eve->prescale_HLT_PFJet40_ = prescale_HLT_PFJet40;
    eve->prescale_HLT_PFJet60_ = prescale_HLT_PFJet60;
    eve->prescale_HLT_PFJet80_ = prescale_HLT_PFJet80;
    eve->prescale_HLT_PFJet140_ = prescale_HLT_PFJet140;
    eve->prescale_HLT_PFJet200_ = prescale_HLT_PFJet200;
    eve->prescale_HLT_PFJet260_ = prescale_HLT_PFJet260;
    eve->prescale_HLT_PFJet320_ = prescale_HLT_PFJet320;
    eve->prescale_HLT_PFJet400_ = prescale_HLT_PFJet400;
    eve->prescale_HLT_PFJet450_ = prescale_HLT_PFJet450;
    eve->prescale_HLT_PFJet500_ = prescale_HLT_PFJet500;
    eve->prescale_HLT_PFJet550_ = prescale_HLT_PFJet550;


	// HLT Object Info

	std::vector<double> pt_trigger_object;
	std::vector<double> eta_trigger_object;
	std::vector<double> phi_trigger_object;

	std::vector<std::vector<string>> filter_trigger_object;

	pt_trigger_object.clear();
	eta_trigger_object.clear();
	phi_trigger_object.clear();

	filter_trigger_object.clear();

	std::vector<string> filter_list;
	filter_list.clear();
	filter_list.push_back("hltEle30erJetC34WPTightGsfTrackIsoFilter");
	filter_list.push_back("hltEle30PFJet35EleCleaned");
	filter_list.push_back("hltEle28erHTT100WPTightGsfTrackIsoFilter");
	filter_list.push_back("hltPFHTJet30");
	filter_list.push_back("hltPFHT150Jet30");
	filter_list.push_back("hltEle27WPTightGsfTrackIsoFilter");
    filter_list.push_back("hltEle32WPTightGsfTrackIsoFilter");
    filter_list.push_back("hltEle32L1DoubleEGWPTightGsfTrackIsoFilter");
	filter_list.push_back("hltEle35noerWPTightGsfTrackIsoFilter");
    filter_list.push_back("hltSingleEle35WPTightGsfL1EGMTTrackIsoFilter");
	filter_list.push_back("hltEle38noerWPTightGsfTrackIsoFilter");
	filter_list.push_back("hltEle40noerWPTightGsfTrackIsoFilter");
	filter_list.push_back("hltSinglePFJet40");
	filter_list.push_back("hltSinglePFJet60");
	filter_list.push_back("hltSinglePFJet80");
	filter_list.push_back("hltSinglePFJet140");
	filter_list.push_back("hltSinglePFJet200");
	filter_list.push_back("hltSinglePFJet260");
	filter_list.push_back("hltSinglePFJet320");
	filter_list.push_back("hltSinglePFJet400");
	filter_list.push_back("hltSinglePFJet450");
	filter_list.push_back("hltSinglePFJet500");
	filter_list.push_back("hltSinglePFJet550");

	if( triggerObjects.isValid() && triggerResults.isValid() ){

		const edm::TriggerNames &names = iEvent.triggerNames(*triggerResults);

		for (pat::TriggerObjectStandAlone obj : *triggerObjects) {

			obj.unpackPathNames(names);
			obj.unpackFilterLabels(iEvent, *triggerResults);

			std::vector<string> labels;
			labels.clear();

			for (unsigned h = 0; h < obj.filterLabels().size(); h++) {
				string filter_label = obj.filterLabels()[h];

				for(unsigned int u=0 ; u < filter_list.size(); u++){
					if( filter_label.compare(filter_list[u]) == 0 ){
						labels.push_back(filter_label);
						break;
					}
				}
			}

			if(labels.size() > 0){
				pt_trigger_object.push_back(obj.pt());
				eta_trigger_object.push_back(obj.eta());
				phi_trigger_object.push_back(obj.phi());
				filter_trigger_object.push_back(labels);
			}
		}
	}

    eve->pt_trigger_object_ = pt_trigger_object;
	eve->eta_trigger_object_ = eta_trigger_object;
	eve->phi_trigger_object_ = phi_trigger_object;

	eve->filter_trigger_object_ = filter_trigger_object;


	// Kinematic info

    int check_pv;

    std::vector<double> ele_pt;
    std::vector<double> ele_eta;
    std::vector<double> ele_sceta;
    std::vector<double> ele_phi;
    std::vector<double> ele_energy;
    std::vector<double> ele_iso;
    std::vector<int> ele_is_ID_Tight;
    std::vector<int> ele_is_ID_Loose;

    std::vector<double> mu_pt;
    std::vector<double> mu_eta;
    std::vector<double> mu_phi;
    std::vector<double> mu_energy;
    std::vector<double> mu_iso;
    std::vector<int> mu_is_ID_Tight;

    std::vector<double> jet_pt;
    std::vector<double> jet_eta;
    std::vector<double> jet_phi;
    std::vector<double> jet_energy;
	std::vector<double> jet_csv;
	std::vector<int> jet_is_ID_Tight;

    double met_pt;
    double met_phi;

    edm::Handle<double> Rho;
    iEvent.getByToken(RhoToken, Rho);
    double rho = *Rho;

    eve->rho_ = rho;

    edm::Handle<vector<reco::Vertex>> vtxHandle;
    iEvent.getByToken(vertexToken,vtxHandle);

    vertex = *(vtxHandle->begin());

    edm::Handle<vector<pat::Electron>> helectrons;
    iEvent.getByToken(electronToken,helectrons);
    vector<pat::Electron> electrons = *helectrons;

    edm::Handle<vector<pat::Muon>> hmuons;
    iEvent.getByToken(muonToken,hmuons);
    vector<pat::Muon> muons = *hmuons;

    edm::Handle<vector<pat::Jet>> hjets;
    iEvent.getByToken(jetToken,hjets);
    vector<pat::Jet> jets = *hjets;

    edm::Handle<vector<pat::MET>> hmet;
    iEvent.getByToken(pfMetToken,hmet);
    vector<pat::MET> met = *hmet;

    //std::cout<<"Step 4 \n";

    check_pv = Check_PV(vertex);

    int n_ele_loose = 0;
    int n_mu_loose = 0;
    int n_jet = 0;

    for( std::vector<pat::Electron>::const_iterator iEle = electrons.begin(); iEle != electrons.end(); iEle++ ){
        pat::Electron ele = *iEle;
        double absEta = fabs(iEle->superCluster()->eta());
        bool inCrack = absEta >= 1.4442 && absEta <= 1.5660;
		if( (iEle->pt() > 10) && (fabs(iEle->eta()) < 2.4) && (is_ele_looseid(ele, rho) == 1) && !inCrack ){
			ele_pt.push_back(iEle->pt());
			ele_eta.push_back(iEle->eta());
			ele_sceta.push_back(iEle->superCluster()->eta());
			ele_phi.push_back(iEle->phi());
			ele_energy.push_back(iEle->energy());
			ele_iso.push_back(get_ele_iso(ele,rho));
			ele_is_ID_Tight.push_back(is_ele_tightid(ele, rho));
            ele_is_ID_Loose.push_back(is_ele_looseid(ele, rho));
            n_ele_loose++;
		}
    }

    for( std::vector<pat::Muon>::const_iterator iMu = muons.begin(); iMu != muons.end(); iMu++ ){
        pat::Muon mu = *iMu;
		if( (iMu->pt() > 10) && (fabs(iMu->eta()) < 2.4) && (is_mu_tightid(mu) == 1) && (get_mu_iso(mu) < 0.25) ){
			mu_pt.push_back(iMu->pt());
			mu_eta.push_back(iMu->eta());
			mu_phi.push_back(iMu->phi());
			mu_energy.push_back(iMu->energy());
			mu_iso.push_back(get_mu_iso(mu));
			mu_is_ID_Tight.push_back(is_mu_tightid(mu));
            n_mu_loose++;
		}
    }

    for( std::vector<pat::Jet>::const_iterator iJet = jets.begin(); iJet != jets.end(); iJet++ ){
		pat::Jet jet = *iJet;
		if( (iJet->pt() > 10) && (fabs(iJet->eta()) < 2.4) && (is_jet_tightid(jet) == 1) ){
			jet_pt.push_back(iJet->pt());
			jet_eta.push_back(iJet->eta());
			jet_phi.push_back(iJet->phi());
			jet_energy.push_back(iJet->energy());
			jet_csv.push_back(get_jet_csv(jet,"pfDeepCSVJetTags:probb") + get_jet_csv(jet,"pfDeepCSVJetTags:probbb"));
			jet_is_ID_Tight.push_back(is_jet_tightid(jet));
            n_jet++;
		}
    }

    met_pt = met.front().pt();
    met_phi = met.front().phi();

    eve->ele_pt_ = ele_pt;
    eve->ele_eta_ = ele_eta;
    eve->ele_sceta_ = ele_sceta;
    eve->ele_phi_ = ele_phi;
    eve->ele_energy_ = ele_energy;
    eve->ele_iso_ = ele_iso;
    eve->ele_is_ID_Tight_ = ele_is_ID_Tight;
    eve->ele_is_ID_Loose_ = ele_is_ID_Loose;

    eve->mu_pt_ = mu_pt;
    eve->mu_eta_ = mu_eta;
    eve->mu_phi_ = mu_phi;
    eve->mu_energy_ = mu_energy;
    eve->mu_iso_ = mu_iso;
    eve->mu_is_ID_Tight_ = mu_is_ID_Tight;

    eve->jet_pt_ = jet_pt;
    eve->jet_eta_ = jet_eta;
    eve->jet_phi_ = jet_phi;
    eve->jet_energy_ = jet_energy;
	eve->jet_csv_ = jet_csv;
	eve->jet_is_ID_Tight_ = jet_is_ID_Tight;

    eve->pfMET_pt_ = met_pt;
    eve->pfMET_phi_ = met_phi;

    //std::cout<<"Step 5 \n";

    // Fill tree for selected events only
    bool event_selection = false;

    if(check_pv==1){
        if(n_ele_loose >= 1){
            if(n_jet >= 1)
                event_selection = true;
        }
    }

    if (event_selection)
        m_ttree->Fill();

}

// ------------ method called once each job just before starting event loop
// ------------
void Trigger_data_analyzer::beginJob()
{

}

// ------------ method called once each job just after ending the event loop
// ------------
void Trigger_data_analyzer::endJob() { return; }

// ------------ method called when starting to processes a run  ------------
void Trigger_data_analyzer::beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup)
{

	bool hltchanged = true;
	if (!hlt_config_.init(iRun, iSetup, hltTag, hltchanged)) {
		std::cout << "Warning, didn't find trigger process HLT with input tag:" << hltTag << std::endl;
		return;
	}

    return;
}

// ------------ method called when ending the processing of a run  ------------
void Trigger_data_analyzer::endRun(const edm::Run &, const edm::EventSetup &)
{
    return;
}

// ------------ method fills 'descriptions' with the allowed parameters for the
// module  ------------
void Trigger_data_analyzer::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
    return;
}

// define this as a CMSSW plugin
DEFINE_FWK_MODULE(Trigger_data_analyzer);

#endif
