import FWCore.ParameterSet.Config as cms

process = cms.Process("Dump")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D12Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D12_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '90X_upgrade2023_realistic_v1', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(	
    'file:/tmp/archie/step2.root'
#   '/store/user/rgupta/Signal/CMSSW_9_1_X_SingleNuE10_Digi_8Apr/170408_135442/0000/step2_SingleNuE10_1.root' 
    )
)

process.dumper = cms.EDAnalyzer("MuonNewGEMDigis",
    simInputLabel = cms.InputTag("g4SimHits","MuonME0Hits"), 
    simTrackLabel = cms.InputTag("g4SimHits"),
    digiInputLabel = cms.InputTag("simMuonGEMDigis"), 
    
    debug = cms.bool(False)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('testGEMDigiReader.root')
)

process.p = cms.Path(process.dumper)
