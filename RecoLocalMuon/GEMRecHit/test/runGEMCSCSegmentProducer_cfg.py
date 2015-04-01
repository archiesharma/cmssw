import FWCore.ParameterSet.Config as cms

process = cms.Process("RERECO")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2019Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2019_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2019', '')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.contentAna = cms.EDAnalyzer("EventContentAnalyzer")


process.load('RecoLocalMuon.GEMRecHit.gemcscSegments_cfi')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
                                      # 'file:out_reco10_paranormal.root',
                                      'file:out_reco10_overlap.root'
   

    )
)

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(
        'file:out_rereco_overlap.root'
    ),
    outputCommands = cms.untracked.vstring(
        'keep  *_gemcscSegments_*_*',
    )
)


process.contentAna = cms.EDAnalyzer("EventContentAnalyzer")
process.reco_step    = cms.Path(process.gemcscSegments)
process.endjob_step  = cms.Path(process.endOfProcess)
process.out_step     = cms.EndPath(process.output)

process.schedule = cms.Schedule(
    process.reco_step,
    process.endjob_step,
    process.out_step
)
