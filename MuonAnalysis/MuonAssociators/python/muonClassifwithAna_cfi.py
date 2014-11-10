import FWCore.ParameterSet.Config as cms
#from RecoMuon.TrackingTools.MuonServiceProxy_cff import MuonServiceProxy
#from SimGeneral.MixingModule.mixNoPU_cfi                          import *
#from SimGeneral.TrackingAnalysis.trackingParticlesNoSimHits_cfi   import * 
from Validation.RecoMuon.selectors_cff import muonTPSet
from SimMuon.MCTruth.MuonAssociatorByHitsESProducer_cfi import *

classByHits = cms.EDAnalyzer("MuonMCClassifAndAna",
    simLabel = cms.InputTag("mix","MergedTrackTruth"),
    muonLabel = cms.InputTag("muons"),
    linkToGenParticles = cms.bool(False),
    genParticles = cms.InputTag("genParticles"),
    muAssocLabel = cms.string("muonAssociatorByHits"),
    decayRho  = cms.double(200),
    decayAbsZ = cms.double(400),
    muonPtCut = cms.double(0.),
    muonMinEtaCut = cms.double(1.6),
    muonMaxEtaCut = cms.double(2.1),
    isoCut = cms.double(100)             #tight 0.05, loose 0.1

)


    
