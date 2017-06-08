import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
process = cms.Process("Demo",eras.Phase2C2_timing)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.Geometry.GeometryExtended2023D4_cff')
process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

#process.GlobalTag.globaltag = '90X_upgrade2023_realistic_v4'
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:0E45AC9C-FB06-E711-ACE4-0CC47A4D762A.root'
         '/store/relval/CMSSW_9_0_0_pre5/RelValZMM_14/GEN-SIM-RECO/PU25ns_90X_upgrade2023_realistic_v4_D4TPU200-v1/00000/0033F30E-7D05-E711-AAD6-0CC47A4C8F18.root'
    )
)

process.demo = cms.EDAnalyzer('GEMSegmentStudies'
)


process.p = cms.Path(process.demo)

                           
