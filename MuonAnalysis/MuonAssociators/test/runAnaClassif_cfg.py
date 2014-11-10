import FWCore.ParameterSet.Config as cms
process = cms.Process("anaClas")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.Geometry.GeometryExtended2023Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.GlobalTag.globaltag = 'DES23_62_V1::All'
process.maxEvents = cms.untracked.PSet(
     input = cms.untracked.int32(-1) 
)


readFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource",fileNames = readFiles)
readFiles.extend( [
#'file:myTestMuons_New_MC.root'


       '/store/user/archie/Muminus_Pt100-gun/CMSSW_620_Patch1_SLHC12_punchthrough_SingleMu2023_NoUpg_alleta_v1/e2ee4bb6f39f0809b677b3ae987fefab/myTestMuons_New_MC_10_2_eUC.root',
       '/store/user/archie/Muminus_Pt100-gun/CMSSW_620_Patch1_SLHC12_punchthrough_SingleMu2023_NoUpg_alleta_v1/e2ee4bb6f39f0809b677b3ae987fefab/myTestMuons_New_MC_11_2_PZE.root',
       '/store/user/archie/Muminus_Pt100-gun/CMSSW_620_Patch1_SLHC12_punchthrough_SingleMu2023_NoUpg_alleta_v1/e2ee4bb6f39f0809b677b3ae987fefab/myTestMuons_New_MC_12_2_dry.root',
       '/store/user/archie/Muminus_Pt100-gun/CMSSW_620_Patch1_SLHC12_punchthrough_SingleMu2023_NoUpg_alleta_v1/e2ee4bb6f39f0809b677b3ae987fefab/myTestMuons_New_MC_13_2_NHG.root',
       '/store/user/archie/Muminus_Pt100-gun/CMSSW_620_Patch1_SLHC12_punchthrough_SingleMu2023_NoUpg_alleta_v1/e2ee4bb6f39f0809b677b3ae987fefab/myTestMuons_New_MC_14_2_lSU.root',
       '/store/user/archie/Muminus_Pt100-gun/CMSSW_620_Patch1_SLHC12_punchthrough_SingleMu2023_NoUpg_alleta_v1/e2ee4bb6f39f0809b677b3ae987fefab/myTestMuons_New_MC_15_2_FEF.root',
       '/store/user/archie/Muminus_Pt100-gun/CMSSW_620_Patch1_SLHC12_punchthrough_SingleMu2023_NoUpg_alleta_v1/e2ee4bb6f39f0809b677b3ae987fefab/myTestMuons_New_MC_16_2_CRm.root',
       '/store/user/archie/Muminus_Pt100-gun/CMSSW_620_Patch1_SLHC12_punchthrough_SingleMu2023_NoUpg_alleta_v1/e2ee4bb6f39f0809b677b3ae987fefab/myTestMuons_New_MC_17_2_Z4N.root',
       '/store/user/archie/Muminus_Pt100-gun/CMSSW_620_Patch1_SLHC12_punchthrough_SingleMu2023_NoUpg_alleta_v1/e2ee4bb6f39f0809b677b3ae987fefab/myTestMuons_New_MC_18_2_pZw.root',
       '/store/user/archie/Muminus_Pt100-gun/CMSSW_620_Patch1_SLHC12_punchthrough_SingleMu2023_NoUpg_alleta_v1/e2ee4bb6f39f0809b677b3ae987fefab/myTestMuons_New_MC_19_2_p0J.root',
       '/store/user/archie/Muminus_Pt100-gun/CMSSW_620_Patch1_SLHC12_punchthrough_SingleMu2023_NoUpg_alleta_v1/e2ee4bb6f39f0809b677b3ae987fefab/myTestMuons_New_MC_1_2_rqt.root',
       '/store/user/archie/Muminus_Pt100-gun/CMSSW_620_Patch1_SLHC12_punchthrough_SingleMu2023_NoUpg_alleta_v1/e2ee4bb6f39f0809b677b3ae987fefab/myTestMuons_New_MC_20_2_VVV.root',
       '/store/user/archie/Muminus_Pt100-gun/CMSSW_620_Patch1_SLHC12_punchthrough_SingleMu2023_NoUpg_alleta_v1/e2ee4bb6f39f0809b677b3ae987fefab/myTestMuons_New_MC_2_2_Us3.root',
       '/store/user/archie/Muminus_Pt100-gun/CMSSW_620_Patch1_SLHC12_punchthrough_SingleMu2023_NoUpg_alleta_v1/e2ee4bb6f39f0809b677b3ae987fefab/myTestMuons_New_MC_3_2_ao1.root',
       '/store/user/archie/Muminus_Pt100-gun/CMSSW_620_Patch1_SLHC12_punchthrough_SingleMu2023_NoUpg_alleta_v1/e2ee4bb6f39f0809b677b3ae987fefab/myTestMuons_New_MC_4_2_6xI.root',
       '/store/user/archie/Muminus_Pt100-gun/CMSSW_620_Patch1_SLHC12_punchthrough_SingleMu2023_NoUpg_alleta_v1/e2ee4bb6f39f0809b677b3ae987fefab/myTestMuons_New_MC_5_2_Okk.root',
       '/store/user/archie/Muminus_Pt100-gun/CMSSW_620_Patch1_SLHC12_punchthrough_SingleMu2023_NoUpg_alleta_v1/e2ee4bb6f39f0809b677b3ae987fefab/myTestMuons_New_MC_6_2_OHG.root',
       '/store/user/archie/Muminus_Pt100-gun/CMSSW_620_Patch1_SLHC12_punchthrough_SingleMu2023_NoUpg_alleta_v1/e2ee4bb6f39f0809b677b3ae987fefab/myTestMuons_New_MC_7_2_ZCT.root',
       '/store/user/archie/Muminus_Pt100-gun/CMSSW_620_Patch1_SLHC12_punchthrough_SingleMu2023_NoUpg_alleta_v1/e2ee4bb6f39f0809b677b3ae987fefab/myTestMuons_New_MC_8_2_y6k.root',
       '/store/user/archie/Muminus_Pt100-gun/CMSSW_620_Patch1_SLHC12_punchthrough_SingleMu2023_NoUpg_alleta_v1/e2ee4bb6f39f0809b677b3ae987fefab/myTestMuons_New_MC_9_2_ct2.root'

]);

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

process.anaclas = cms.EDAnalyzer('AnaClassif',
                               muonPtCut = cms.double( 5 ),	## put this 45 for tight selection
                               isoCut = cms.double( 10 )	## put this 0.05 (loose) or 0.10 (tight) isolation in cone 0.3
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('test_SingleMu23NoUpg_allEta.root')
)

##triggers
#from HLTrigger.HLTfilters.hltHighLevel_cfi import *
#process.triggerSelection = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::HLT", HLTPaths = ["HLT_Mu40_eta2p1_v*"])


process.p = cms.Path(
                     #process.scrapingVeto*
        	     #process.primaryVertexFilter*
		     process.anaclas
		 ) 

