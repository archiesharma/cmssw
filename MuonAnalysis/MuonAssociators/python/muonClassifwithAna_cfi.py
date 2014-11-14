import FWCore.ParameterSet.Config as cms
#from RecoMuon.TrackingTools.MuonServiceProxy_cff import MuonServiceProxy
#from SimGeneral.MixingModule.mixNoPU_cfi                          import *
#from SimGeneral.TrackingAnalysis.trackingParticlesNoSimHits_cfi   import * 


from Validation.RecoMuon.selectors_cff import muonTPSet
from SimMuon.MCTruth.MuonAssociatorByHitsESProducer_cfi import *

selectedVertices = cms.EDFilter("VertexSelector",
   src = cms.InputTag('offlinePrimaryVertices'),
   cut = cms.string("isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0 & abs(z) < 24 & abs(position.Rho) < 2."),
   filter = cms.bool(False)
)


classByHits = cms.EDAnalyzer("MuonMCClassifAndAna",
    simLabel = cms.InputTag("mix","MergedTrackTruth"),
    muonLabel = cms.InputTag("muons"),
    linkToGenParticles = cms.bool(False),
    genParticles = cms.InputTag("genParticles"),
    muAssocLabel = cms.string("muonAssociatorByHits"),
    decayRho  = cms.double(200),
    decayAbsZ = cms.double(400),
    muonPtCut = cms.double(0.),
    muonMinEtaCut = cms.double(1.5),
    muonMaxEtaCut = cms.double(2.1),
    isoCut = cms.double(100),             #tight 0.05, loose 0.1
    vxtTag = cms.InputTag("selectedVertices")
)

classByHitsPtCut15 = classByHits.clone(
 muonPtCut = cms.double(15.),
 muonMinEtaCut = cms.double(1.5),
 muonMaxEtaCut = cms.double(2.1),
)

classByHitsPtCut3 = classByHits.clone(
 muonPtCut = cms.double(3.),
 muonMinEtaCut = cms.double(1.5),
 muonMaxEtaCut = cms.double(2.1),
)

classByHitsPtCut5 = classByHits.clone(
 muonPtCut = cms.double(5.),
 muonMinEtaCut = cms.double(1.5),
 muonMaxEtaCut = cms.double(2.1),
) 

muonClassificationByHits = cms.Sequence(
 classByHits +
 classByHitsPtCut15 +
 classByHitsPtCut3 +
 classByHitsPtCut5
)   
