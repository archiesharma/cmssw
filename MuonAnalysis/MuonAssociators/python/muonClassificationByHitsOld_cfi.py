import FWCore.ParameterSet.Config as cms
#from RecoMuon.TrackingTools.MuonServiceProxy_cff import MuonServiceProxy
#from SimGeneral.MixingModule.mixNoPU_cfi                          import *
#from SimGeneral.TrackingAnalysis.trackingParticlesNoSimHits_cfi   import * 
from Validation.RecoMuon.selectors_cff import muonTPSet
from SimMuon.MCTruth.MuonAssociatorByHitsESProducer_cfi import *		#muonAssociatorByHits

classByHitsB = cms.EDProducer("MuonMCClassifier",
    simLabel = cms.InputTag("mix","MergedTrackTruth"),		##ok
    muonLabel = cms.InputTag("muons"),				##ok
    linkToGenParticles = cms.bool(False),          ##ok  produce also a collection of GenParticles for secondary muons
    genParticles = cms.InputTag("genParticles"),  ##ok  and associations to primary and secondaries
    muAssocLabel = cms.string("muonAssociatorByHits"),	##ok
    decayRho  = cms.double(200), # to classifiy differently decay muons included in ppMuX
    decayAbsZ = cms.double(400), # and decay muons that could not be in ppMuX
    muonPtCut = cms.double(1.),
    muonMinEtaCut = cms.double(0.),
    muonMaxEtaCut = cms.double(1.2)
)

classByHitsE = classByHitsB.clone(
    muonMinEtaCut = cms.double(1.2),
    muonMaxEtaCut = cms.double(3.)
)

classByHitsO = classByHitsB.clone(
    muonMinEtaCut = cms.double(0.9),
    muonMaxEtaCut = cms.double(1.2)
)

classByHitsU = classByHitsB.clone(
    muonMinEtaCut = cms.double(1.6),
    muonMaxEtaCut = cms.double(3.)
)

muonClassificationByHits = cms.Sequence(
      classByHitsB +
      classByHitsE +
      classByHitsO +
      classByHitsU 

)



    
