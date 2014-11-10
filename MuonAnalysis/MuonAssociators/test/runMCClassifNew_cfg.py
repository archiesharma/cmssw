import FWCore.ParameterSet.Config as cms
process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.Geometry.GeometryExtended2023Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.GlobalTag.globaltag = 'DES23_62_V1::All'
process.maxEvents = cms.untracked.PSet(
     input = cms.untracked.int32(10000) 
)

#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
#process.Tracer = cms.Service("Tracer")

readFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource",fileNames = readFiles)
readFiles.extend( [


      '/store/user/calabria/MinBias_TuneZ2star_14TeV-pythia6/calabria_MinBias3_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023ScenarioNoUpg_3Step/7e82c0353afc1625567b3844e4c947b8/step3_100_1_q5W.root'
#'file:/tmp/archie/step3_1_1_Kc7.root'
# '/store/user/archie/MinBias_TuneZ2star_14TeV-pythia6/CMSSW_620_Patch1_SLHC12_punchthrough_MinBiasNoBkg_Upg_Idcuts_alleta_v1/e2ee4bb6f39f0809b677b3ae987fefab/myTestMuons_New_MC_9_1_9jD.root'
#'/store/user/archie/calabria_MinBias_GEN-SIM_CMSSW_6_2_0_SLHC12_2017Scenario_1Step/CMSSW_620_Patch1_SLHC12_2017Scenario_3Step_MinBias_GEN-SIM-DIGI-RECO_new/18a63ba141c3c40ec5c6c092f4e65730/step3_1_1_Kc7.root'
##with background:
#'file:/afs/cern.ch/work/m/mileva/punchThrough/CMSSW_6_2_0_SLHC12/src/files_bkg/step3_1000_1_m8t.root',
#'file:/afs/cern.ch/work/m/mileva/punchThrough/CMSSW_6_2_0_SLHC12/src/files_bkg/step3_1001_1_KZe.root',
#'file:/afs/cern.ch/work/m/mileva/punchThrough/CMSSW_6_2_0_SLHC12/src/files_bkg/step3_100_1_5UV.root',
#'file:/afs/cern.ch/work/m/mileva/punchThrough/CMSSW_6_2_0_SLHC12/src/files_bkg/step3_101_1_E8S.root',
#'file:/afs/cern.ch/work/m/mileva/punchThrough/CMSSW_6_2_0_SLHC12/src/files_bkg/step3_102_1_jUo.root',
#'file:/afs/cern.ch/work/m/mileva/punchThrough/CMSSW_6_2_0_SLHC12/src/files_bkg/step3_103_1_MNx.root'


##read files from eos:
#'root://eoscms//eos/cms/store/user/mileva/gemPunchThrough/step3_1000_1_m8t.root',
#'root://eoscms//eos/cms/store/user/mileva/gemPunchThrough/step3_1001_1_KZe.root',
#'root://eoscms//eos/cms/store/user/mileva/gemPunchThrough/step3_100_1_5UV.root'


##without background:
#'file:/afs/cern.ch/work/m/mileva/punchThrough/CMSSW_6_2_0_SLHC12/src/files_Nobkg/step3_1000_1_o1X.root',
#'file:/afs/cern.ch/work/m/mileva/punchThrough/CMSSW_6_2_0_SLHC12/src/files_Nobkg/step3_1001_1_At2.root',
#'file:/afs/cern.ch/work/m/mileva/punchThrough/CMSSW_6_2_0_SLHC12/src/files_Nobkg/step3_100_1_j9X.root',
#'file:/afs/cern.ch/work/m/mileva/punchThrough/CMSSW_6_2_0_SLHC12/src/files_Nobkg/step3_101_1_ENu.root',
#'file:/afs/cern.ch/work/m/mileva/punchThrough/CMSSW_6_2_0_SLHC12/src/files_Nobkg/step3_102_1_UIO.root'

]);

process.checkMuSize = cms.EDFilter("MuonsizeFilter"
                                            )
## require scraping filter
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                     applyfilter = cms.untracked.bool(True),
                                     debugOn = cms.untracked.bool(False),
                                     numtrack = cms.untracked.uint32(10),
                                     thresh = cms.untracked.double(0.2)
                                     )

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                            vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                            minimumNDOF = cms.uint32(4) ,
                                            maxAbsZ = cms.double(24),
                                            maxd0 = cms.double(2)
                                            )


#from SimGeneral.MixingModule.mixNoPU_cfi                          import *
#from SimGeneral.TrackingAnalysis.trackingParticlesNoSimHits_cfi   import * 
#from SimMuon.MCTruth.MuonAssociatorByHitsESProducer_NoSimHits_cfi import *

process.load('MuonAnalysis.MuonAssociators.muonClassificationByHitsNew_cfi')

##triggers
#from HLTrigger.HLTfilters.hltHighLevel_cfi import *
#process.triggerSelection = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::HLT", HLTPaths = ["HLT_Mu40_eta2p1_v*"])


process.p = cms.Path(
                     process.checkMuSize*
#                     process.scrapingVeto*
        	     process.primaryVertexFilter*
                     (
#                     process.mix +
#                     process.trackingParticlesNoSimHits +
	             process.classByHits  
#                     process.classByHitsE +
#                     process.classByHitsO +
#                     process.classByHitsU 
                     )

		 ) 

#process.MessageLogger.categories += [ 'MuonMCClassifier' ]
#process.MessageLogger.cerr.MuonMCClassifier = cms.untracked.PSet(
#    optionalPSet = cms.untracked.bool(False),
#    limit = cms.untracked.int32(10000000)
#)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("myTestMuons_New_MC.root"),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep recoMuons_*_*_*",
        "keep *_*_*_Demo",
        "keep *_genParticles_*_*",
        "keep recoTrackExtras_standAloneMuons_*_*",          ## track states at the muon system, used both by patMuons and standAloneMuons
        "keep recoTracks_standAloneMuons_*_*",                ## bare standalone muon tracks, using standalone muon momentum (without BS constraint)
        "keep recoTracks_refittedStandAloneMuons_*_*",                ## bare standalone muon tracks, using standalone muon momentum (without BS constraint)
        "keep recoTrackExtras_refittedStandAloneMuons_*_*",
        "keep *_generalTracks_*_*",
#        "keep edmTriggerResults_*_*_HLT",                    
##        "keep l1extraL1MuonParticles_l1extraParticles_*_*",  
        "keep *_offlinePrimaryVertices_*_*",                 
        "keep *_offlineBeamSpot_*_*"                     
    ),
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring("p")),
)
process.end = cms.EndPath(process.out)


