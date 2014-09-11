import FWCore.ParameterSet.Config as cms

process = cms.Process("EDMtoMEConvert")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("DQMServices.Components.EDMtoMEConverter_cff")
process.load("DQMServices.Components.DQMEnvironment_cfi")
#process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

process.load("Validation.Configuration.postValidation_cff")
process.load("Validation.RecoMuon.PostProcessorHLT_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'file:validationEDM.root'

'file:/tmp/archie/validationEDM_10_1_L6s.root',
'file:/tmp/archie/validationEDM_11_1_ozA.root',
'file:/tmp/archie/validationEDM_12_1_ubU.root',
'file:/tmp/archie/validationEDM_13_1_zQG.root',
'file:/tmp/archie/validationEDM_14_1_JDl.root',
'file:/tmp/archie/validationEDM_15_1_oBY.root',
'file:/tmp/archie/validationEDM_16_1_Wf2.root',
'file:/tmp/archie/validationEDM_17_1_1E3.root',
'file:/tmp/archie/validationEDM_18_1_wYH.root',
'file:/tmp/archie/validationEDM_19_1_w3i.root',
'file:/tmp/archie/validationEDM_1_1_qrI.root',
'file:/tmp/archie/validationEDM_20_1_U6t.root',
'file:/tmp/archie/validationEDM_21_1_1iV.root',
'file:/tmp/archie/validationEDM_2_1_0cj.root',
'file:/tmp/archie/validationEDM_4_1_ows.root',
'file:/tmp/archie/validationEDM_5_1_bec.root',
'file:/tmp/archie/validationEDM_7_1_t6d.root',
'file:/tmp/archie/validationEDM_9_1_2bG.root',
'file:/tmp/archie/validationEDM_8_1_tCV.root'




)
)

process.DQMStore.referenceFileName = ""
process.DQMStore.collateHistograms = False

process.dqmSaver.convention = "Offline"
#Settings equivalent to 'RelVal' convention:
process.dqmSaver.saveByRun = cms.untracked.int32(-1)
process.dqmSaver.saveAtJobEnd = cms.untracked.bool(True)
process.dqmSaver.forceRunNumber = cms.untracked.int32(1)
#End of 'RelVal convention settings
process.dqmSaver.workflow = "/GlobalValidation/Test/RECO"


from Validation.RecoMuon.PostProcessor_cff import *
process.p1 = cms.Path(process.EDMtoMEConverter*
#                      process.postValidation*
                       process.recoMuonPostProcessors* 
#                      process.recoMuonPostProcessorsHLT*
                      process.dqmSaver)
