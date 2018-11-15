#ifndef Trigger_data_analyzer_h
#define Trigger_data_analyzer_h

/// Core libraries
#include <cmath>  // arctan
#include <cstdio> // printf, fprintf
#include <fstream>
#include <memory>
#include <stdexcept> // standard exceptions
#include <vector>
#include <map>

/// CMSSW user include files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidateFwd.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1HFRings.h"
#include "DataFormats/L1Trigger/interface/L1HFRingsFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerEvmReadoutRecord.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h"
#include "DataFormats/L1Trigger/interface/L1HFRings.h"
#include "DataFormats/L1Trigger/interface/L1HFRingsFwd.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/strbitset.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoBTag/BTagTools/interface/SignedTransverseImpactParameter.h"

#include "Trigger/TriggerAnalyzer/interface/TriggerStudyEventVars.h"

#include "LHAPDF/LHAPDF.h"

/// ROOT includes
#include "TH1.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TMath.h"


class Trigger_data_analyzer : public edm::EDAnalyzer
{
  public:
    explicit Trigger_data_analyzer(const edm::ParameterSet &);
    ~Trigger_data_analyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions &);

  private:
    /*
     * Function/method section
    */

    /// Standard EDAnalyzer functions
    void analyze(const edm::Event &, const edm::EventSetup &) override;

    void beginJob() override;
    void endJob() override;

    void beginRun(const edm::Run &, const edm::EventSetup &) override;
    void endRun(const edm::Run &, const edm::EventSetup &) override;

	void ntuple_initialize();
    double get_ele_iso(const pat::Electron &, double &);
    int is_ele_tightid(const pat::Electron &, double &);
    int is_ele_looseid(const pat::Electron &, double &);
    double get_mu_iso(const pat::Muon &);
    int is_mu_tightid(const pat::Muon &);
	int is_jet_tightid(const pat::Jet &);
	double get_jet_csv(const pat::Jet& , const std::string );
    int Check_PV(reco::Vertex &);

    edm::EDGetTokenT <edm::TriggerResults> triggerResultsToken;
    edm::EDGetTokenT <pat::PackedTriggerPrescales> triggerPrescalesToken;
    edm::EDGetTokenT <edm::TriggerResults> filterResultsToken;
    edm::EDGetTokenT <BXVector<GlobalAlgBlk>> gtStage2DigisToken;
	edm::EDGetTokenT <trigger::TriggerEvent> triggerEventToken;
	edm::EDGetTokenT <pat::TriggerObjectStandAloneCollection> triggerObjectsToken;

    edm::EDGetTokenT<double> RhoToken;
    edm::EDGetTokenT <vector<reco::Vertex>> vertexToken;
    edm::EDGetTokenT <vector<pat::Electron>> electronToken;
    edm::EDGetTokenT <vector<pat::Muon>> muonToken;
    edm::EDGetTokenT <vector<pat::Jet>> jetToken;
    edm::EDGetTokenT <vector<pat::MET>> pfMetToken;

    reco::Vertex vertex;

	HLTConfigProvider hlt_config_;
	HLTConfigProvider filter_config_;

    std::string hltTag;
    std::string filterTag;

    std::map<std::string, int> l1talgo_cppath_;

    std::map<std::string, int> hlt_cppath_;
    std::map<std::string, int> flt_cppath_;

    std::vector<std::string> hlt_triggerNames_;
    std::vector<std::string> flt_filterNames_;

    edm::Service<TFileService> fs_;

	triggerStudyEventVars *eve;
    TTree *m_ttree;

};

#endif
