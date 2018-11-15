import os
import FWCore.ParameterSet.Config as cms

from RecoJets.Configuration.RecoJets_cff import *
from RecoJets.Configuration.RecoPFJets_cff import *
from JetMETCorrections.Configuration.JetCorrectionProducersAllAlgos_cff import *
from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

process = cms.Process("MAOD")

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())

#process.source = cms.Source("PoolSource",
#     fileNames = cms.untracked.vstring(
#     '/store/data/Run2018A/SingleMuon/MINIAOD/PromptReco-v1/000/315/257/00000/1C3A2CCE-604B-E811-99DA-FA163E1FBD5A.root',
#	  '/store/data/Run2018A/EGamma/MINIAOD/PromptReco-v1/000/315/255/00000/2C278415-5B4B-E811-9447-FA163E62E0F4.root'
#     )
#)

#import FWCore.PythonUtilities.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename = 'data/Cert_314472-316271_13TeV_PromptReco_Collisions18_JSON.txt').getVLuminosityBlockRange()

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
#process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")

process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_cff" )
#For 2018 Prompt-Reco
process.GlobalTag.globaltag = '101X_dataRun2_Prompt_v11'
#For 2018 Re-Reco
#process.GlobalTag.globaltag = ''

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(True)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))


process.trigger = cms.EDAnalyzer('Trigger_data_analyzer',
                                 HLTsource = cms.untracked.string("HLT"),
                                 PATsource = cms.untracked.string("RECO")
                                 )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("trigger_ntuple.root")
                                   )

process.p = cms.Path(process.trigger)




