from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'trigger_data_2018_SingleMuon_2018A_v2_promptreco'

config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'test/trigger_data_PromptReco.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = ['data']
config.JobType.maxJobRuntimeMin = 3600
config.JobType.maxMemoryMB = 2500

config.section_("Data")
config.Data.inputDataset = '/SingleMuon/Run2018A-PromptReco-v2/MINIAOD'
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.allowNonValidInputDataset = True
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = 'data/Cert_314472-323523_13TeV_PromptReco_Collisions18_JSON.txt'
config.Data.unitsPerJob = 20
config.Data.publication = False
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.outputDatasetTag = 'trigger_data_2018_SingleMuon_2018A_v2_promptreco'
config.Data.outLFNDirBase = '/store/user/abdatta/Trigger_2018/Data_Prompt_Reco_2018/'
config.Data.ignoreLocality = False

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
