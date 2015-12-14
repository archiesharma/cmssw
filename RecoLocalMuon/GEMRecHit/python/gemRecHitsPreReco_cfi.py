import FWCore.ParameterSet.Config as cms

gemRecHitsPreReco = cms.EDProducer("GEMRecHitProducerPreReco",
    recAlgoConfig = cms.PSet(
                    recOnlyMuons = cms.bool(True),
    ),
    prerecAlgo = cms.string('GEMRecHitStandardAlgoPreReco'),
    gemPreDigiLabel = cms.InputTag("simMuonGEMDigiswithTiming"),
)


