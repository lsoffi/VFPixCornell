from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'GGH_HToZZTo4L_FPix_1212'

config.section_("JobType")
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ConfFile_cfg.py'

config.section_("Data")
#config.Data.inputDataset = '/VBF_HToZZTo4L_M-125_14TeV-powheg-pythia6/Summer12-PU50bx25_POSTLS161_V12-v1/GEN-SIM-RAW-RECO'
#config.Data.inputDataset = '/VBF_HToMuMu_M-125_14TeV-powheg-pythia6/TP2023SHCALDR-SHCALMar26_PU140BX25_PH2_1K_FB_V6-v1/GEN-SIM-RECO'
#config.Data.inputDataset = '/GluGluToHToZZTo4L_M-125_14TeV-powheg-pythia6/Summer12-PU50bx25_POSTLS161_V12-v1/GEN-SIM-RAW-RECO'
config.Data.inputDataset = '/GluGluToHToZZTo4L_M-125_14TeV-powheg-pythia6/TP2023SHCALDR-SHCALJan23_PU140BX25_PH2_1K_FB_V6-v1/GEN-SIM-RECO'
#'/GluGluToHToZZTo4L_M-125_14TeV-powheg-pythia6/TP2023SHCALDR-SHCALMar26_PU140BX25_PH2_1K_FB_V6-v1/GEN-SIM-RECO'
#config.Data.inputDataset ='/VBF_HToZZTo4L_M-125_14TeV-powheg-pythia6/TP2023SHCALDR-SHCALMar26_PU140BX25_PH2_1K_FB_V6-v3/GEN-SIM-RECO'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

#config.Data.lumiMask = 'http://soffi.web.cern.ch/soffi/jsonGood_0716.txt'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/DCSOnly/json_DCSONLY_Run2015B.txt'
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/soffi/'
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'GGH_HToZZTo4L_FPix_1212'
config.section_('User')
config.section_("Site")
config.Site.storageSite = 'T2_IT_Rome'

