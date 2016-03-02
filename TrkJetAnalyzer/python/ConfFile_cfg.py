import FWCore.ParameterSet.Config as cms
from RecoJets.JetProducers.TrackJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *

process = cms.Process("TrkJetPhase2")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(150))

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#VBF invisible
#'/store/mc/TP2023SHCALDR/VBF_HToInv_M-125_14TeV-powheg-pythia6/GEN-SIM-RECO/SHCALMar26_PU140BX25_PH2_1K_FB_V6-v1/50000/FA1FEEE5-AEDB-E411-81E1-002618FDA237.root',
#'/store/mc/TP2023SHCALDR/VBF_HToInv_M-125_14TeV-powheg-pythia6/GEN-SIM-RECO/SHCALMar26_PU140BX25_PH2_1K_FB_V6-v1/50000/FAF1FE37-7ADB-E411-B660-00259059642A.root',
#'/store/mc/TP2023SHCALDR/VBF_HToInv_M-125_14TeV-powheg-pythia6/GEN-SIM-RECO/SHCALMar26_PU140BX25_PH2_1K_FB_V6-v1/50000/FC525A03-55DB-E411-8398-0025905B8596.root',
#'/store/mc/TP2023SHCALDR/VBF_HToInv_M-125_14TeV-powheg-pythia6/GEN-SIM-RECO/SHCALMar26_PU140BX25_PH2_1K_FB_V6-v1/50000/FC86F1EB-28DB-E411-85E1-00261894396D.root',
'/store/mc/TP2023SHCALDR/VBF_HToInv_M-125_14TeV-powheg-pythia6/GEN-SIM-RECO/SHCALMar26_PU140BX25_PH2_1K_FB_V6-v1/50000/FE5A38D0-81DB-E411-8B8B-0025905A608A.root')
#VBF 4leptons
#'/store/mc/TP2023SHCALDR/VBF_HToZZTo4L_M-125_14TeV-powheg-pythia6/GEN-SIM-RECO/SHCALJan23_PU140BX25_PH2_1K_FB_V6-v2/10000/30E51CAA-C6A8-E411-97ED-C4346BC7EDD8.root')
#root://xrootd.unl.edu//store/mc/TP2023SHCALDR/VBF_HToZZTo4L_M-125_14TeV-powheg-pythia6/GEN-SIM-RECO/SHCALMar26_PU140BX25_PH2_1K_FB_V6-v3/00000/02449722-9DE5-E411-838A-E0CB4E1A1147.root')
#'/store/mc/TP2023SHCALDR/VBF_HToZZTo4L_M-125_14TeV-powheg-pythia6/GEN-SIM-RECO/SHCALMar26_PU140BX25_PH2_1K_FB_V6-v3/00000/02449722-9DE5-E411-838A-E0CB4E1A1147.root',
#'/store/mc/TP2023SHCALDR/VBF_HToZZTo4L_M-125_14TeV-powheg-pythia6/GEN-SIM-RECO/SHCALMar26_PU140BX25_PH2_1K_FB_V6-v3/00000/040E818A-ACE5-E411-BBA0-20CF30561711.root',
#'/store/mc/TP2023SHCALDR/VBF_HToZZTo4L_M-125_14TeV-powheg-pythia6/GEN-SIM-RECO/SHCALMar26_PU140BX25_PH2_1K_FB_V6-v3/00000/0449CCF7-AEE5-E411-9510-00259073E500.root',
#'/store/mc/TP2023SHCALDR/VBF_HToZZTo4L_M-125_14TeV-powheg-pythia6/GEN-SIM-RECO/SHCALMar26_PU140BX25_PH2_1K_FB_V6-v3/00000/062285A5-95E5-E411-BAA6-E0CB4E29C4E6.root')
#'/store/mc/TP2023SHCALDR/VBF_HToZZTo4L_M-125_14TeV-powheg-pythia6/GEN-SIM-RECO/SHCALMar26_PU140BX25_PH2_1K_FB_V6-v3/00000/0A2A9592-9BE5-E411-B238-90E6BAE8CC18.root' )
#GGH 4leptons
#'/store/mc/TP2023SHCALDR/GluGluToHToZZTo4L_M-125_14TeV-powheg-pythia6/GEN-SIM-RECO/SHCALMar26_PU140BX25_PH2_1K_FB_V6-v1/10000/009B3F0C-EAD9-E411-B790-E0CB4E19F967.root',
#'/store/mc/TP2023SHCALDR/GluGluToHToZZTo4L_M-125_14TeV-powheg-pythia6/GEN-SIM-RECO/SHCALMar26_PU140BX25_PH2_1K_FB_V6-v1/10000/00A16110-F5D9-E411-9858-00259074AE54.root',
#'/store/mc/TP2023SHCALDR/GluGluToHToZZTo4L_M-125_14TeV-powheg-pythia6/GEN-SIM-RECO/SHCALMar26_PU140BX25_PH2_1K_FB_V6-v1/10000/0450F75E-EAD9-E411-B350-0025907B4EF2.root',
#'/store/mc/TP2023SHCALDR/GluGluToHToZZTo4L_M-125_14TeV-powheg-pythia6/GEN-SIM-RECO/SHCALMar26_PU140BX25_PH2_1K_FB_V6-v1/10000/0826A6E8-A2D7-E411-B32F-E0CB4E29C4C4.root',
#'/store/mc/TP2023SHCALDR/GluGluToHToZZTo4L_M-125_14TeV-powheg-pythia6/GEN-SIM-RECO/SHCALMar26_PU140BX25_PH2_1K_FB_V6-v1/10000/0AE4A0E7-A2D7-E411-8CA9-E0CB4EA0A904.root' )
#'/store/mc/RunIISpring15MiniAODv2/GluGluHToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/00289B81-C56D-E511-8C21-0025905C2CBC.root')
)

#produce trk jet
TrackJetParameters = cms.PSet(
    src            = cms.InputTag('trackRefsForJets'),
    srcPVs         = cms.InputTag('offlinePrimaryVertices'),
    jetType        = cms.string('TrackJet'),
    doOutputJets   = cms.bool(True),
    jetPtMin       = cms.double(3.0),
    inputEMin      = cms.double(0.0),
    inputEtMin     = cms.double(0.0),
    doPVCorrection = cms.bool(False),
    # pileup with offset correction
    doPUOffsetCorr = cms.bool(False),
    # if pileup is false, these are not read:
    nSigmaPU = cms.double(1.0),
    radiusPU = cms.double(0.5),  
    # fastjet-style pileup     
    doAreaFastjet   = cms.bool(False),
    doRhoFastjet    = cms.bool(False),
    doAreaDiskApprox= cms.bool( False),
    voronoiRfact    = cms.double(-0.9),
    # if doPU is false, these are not read:
    Active_Area_Repeats = cms.int32(1),
    GhostArea = cms.double(0.01),
    Ghost_EtaMax = cms.double(5.0),
    # only use the tracks that were used to fit the vertex
    UseOnlyVertexTracks = cms.bool(False),
    # only consider the highest-sum-pT PV for clustering
    UseOnlyOnePV        = cms.bool(False),
    # maximum z-distance between track and vertex for association (in cm)
    DzTrVtxMax          = cms.double(10),
    # maximum xy-distance between track and vertex for association (in cm)
    DxyTrVtxMax         = cms.double(0.2),
    # minimum number of degrees of freedom to call a PV a good vertex
    MinVtxNdof          = cms.int32(5),
    # maximum z distance to origin to c  pvTrueX_ = genParticles->at (2).vx ();
   
    MaxVtxZ             = cms.double(15.),
    useDeterministicSeed= cms.bool( True ),
    minSeed             = cms.uint32( 14327 )
)


process.ak5TrackJets = cms.EDProducer(
    "FastjetJetProducer",
    TrackJetParameters,
   # DzTrVtxMax          = cms.double(1),  
    AnomalousCellParameters,
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.5)
    )




process.TrkJetAnalyzer = cms.EDAnalyzer('TrkJetAnalyzer',
                                      GenJetCollection       = cms.InputTag("ak5GenJets"),
                                      GenParticlesCollection = cms.InputTag("genParticles"),
                                      TrackJetCollection    = cms.InputTag("ak5TrackJets"),
                                      TrackCollection    = cms.InputTag("generalTracks"),
                                    #  SimTracksCollection = cms.InputTag ("g4SimHits", ""),
                                      VertexCollection = cms.InputTag("offlinePrimaryVertices"),
                                      PUSummaryInfoCollection = cms.InputTag("addPileupInfo"),
                                    #  PfCandidatesCollection = cms.InputTag ("particleFlow"),
                                      minTracks = cms.untracked.uint32(1000)
         )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('histoOUTPUT_VBF.root')
                                   )
#dump event content
process.dump=cms.EDAnalyzer('EventContentAnalyzer')


process.p = cms.Path(process.ak5TrackJets*process.TrkJetAnalyzer)
