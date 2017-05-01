import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.Geometry.GeometryExtended2023D4_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = '90X_upgrade2023_realistic_v4'
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:0E45AC9C-FB06-E711-ACE4-0CC47A4D762A.root'
    )
)

process.demo = cms.EDAnalyzer('GEMSegmentStudies'
)


process.p = cms.Path(process.demo)
