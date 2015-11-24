import FWCore.ParameterSet.Config as cms

# Module to create simulated GEM Pre Reco digis.
simMuonGEMDigiswithTiming = cms.EDProducer("GEMDigiPreRecoProducer",
    inputCollection = cms.string('g4SimHitsMuonGEMHits'),
    digiPreRecoModelString = cms.string('PreRecoGaussian'),
    timeResolution = cms.double(0.001), # in ns
    phiResolution = cms.double(0.05),   # in cm average resolution along local x in case of no correlation
    etaResolution = cms.double(1.),     # in cm average resolution along local y in case of no correlation
    useCorrelation = cms.bool(False),
    useEtaProjectiveGEO = cms.bool(False),
    averageEfficiency = cms.double(0.98),
    digitizeOnlyMuons = cms.bool(True),
    # simulateIntrinsicNoise = cms.bool(False), # intrinsic noise --> not implemented
    # averageNoiseRate = cms.double(0.001),     # intrinsic noise --> not implemented
    simulateElectronBkg = cms.bool(False),       # True - will simulate electron background
    simulateNeutralBkg  = cms.bool(False),       # True - will simulate neutral (n+g)  background
    minBunch = cms.int32(-5),                   # [x 25 ns], forms the readout window together with maxBunch,
    maxBunch = cms.int32(3),                    # we should think of shrinking this window ...
)
