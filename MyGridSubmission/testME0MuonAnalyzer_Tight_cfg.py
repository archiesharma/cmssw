import FWCore.ParameterSet.Config as cms

process = cms.Process("Test")


process.load('Configuration.Geometry.GeometryExtended2023MuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023Muon_cff')


process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')

process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")

process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")

process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")

process.load( "DQMServices/Core/DQMStore_cfg" )

process.load('Validation.RecoMuon.associators_cff')
process.load('Validation.RecoMuon.selectors_cff')
process.load('Validation.RecoMuon.MuonTrackValidator_cfi')
process.load('Validation.RecoMuon.RecoMuonValidator_cfi')



process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2019', '')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50)
)



process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # 'file:/store/group/upgrade/muon/ME0GlobalReco/ME0MuonReRun_DY_SLHC23patch1/M-20_TuneZ2star_14TeV_6_2_0_SLHC23patch1_2023/output_9_1_hVY.root'
        'file:/cmshome/piet/SLC6/ME0_Studies/SLHC26_InTimeOutOfTimePU/CMSSW_6_2_0_SLHC26_patch3/src/MyAnalyzers/MyME0InTimePUAnalyzer/DYToMuMu_M-20_HGCALGS_PU140_ME0_RECO_100ps_amandeep_100_v2.root'
        )


)



from Validation.RecoMuon.selectors_cff import *
from Validation.RecoMuon.associators_cff import *
# Configurations for MuonTrackValidators
import Validation.RecoMuon.MuonTrackValidator_cfi


# Configurations for RecoMuonValidators
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from Validation.RecoMuon.RecoMuonValidator_cfi import *

#import SimGeneral.MixingModule.mixNoPU_cfi
from SimMuon.MCTruth.MuonAssociatorByHitsESProducer_NoSimHits_cfi import *
from SimMuon.MCTruth.MuonAssociatorByHits_cfi import muonAssociatorByHitsCommonParameters

process.TrackAssociatorByChi2ESProducer = Validation.RecoMuon.associators_cff.TrackAssociatorByChi2ESProducer.clone(chi2cut = 6.0,ComponentName = 'TrackAssociatorByChi2')

process.recoMuonValidation = cms.Sequence(#probeTracks_seq*
    #(selectedVertices * selectedFirstPrimaryVertex) * 
    #bestMuonTuneP_seq*
    #muonColl_seq*trackColl_seq*extractedMuonTracks_seq*bestMuon_seq*trackerMuon_seq*
    me0muonColl_seq
    #((process.muonValidation_seq))
    )

process.Test = cms.EDAnalyzer("ME0MuonAnalyzer",
                              HistoFile = cms.string('ME0MuonAnalyzerOutput.root'),
                              HistoFolder = cms.string('ME0MuonAnalyzerOutput'),
                              FakeRatePtCut = cms.double(5.0),
                              MatchingWindowDelR = cms.double (.15),
                              RejectEndcapMuons = cms.bool(False),
                              UseAssociators = cms.bool(True),
                              associators = cms.vstring('TrackAssociatorByChi2'),
                              #label = cms.VInputTag('bestMuon'),
                              label = cms.VInputTag('me0muon'),
                              #associatormap = cms.InputTag("tpToMuonTrackAssociation"),

                              # selection of GP for evaluation of efficiency
                              ptMinGP = cms.double(0.9),
                              minRapidityGP = cms.double(-2.5),
                              maxRapidityGP = cms.double(2.5),
                              tipGP = cms.double(3.5),
                              lipGP = cms.double(30.0),
                              chargedOnlyGP = cms.bool(True),
                              statusGP = cms.int32(1),
                              pdgIdGP = cms.vint32(13, -13),
                              #parametersDefiner = cms.string('LhcParametersDefinerForTP'),
                              ME0MuonSelectionType = cms.string('Tight'),
)

process.p = cms.Path(process.recoMuonValidation*process.Test)



# process.PoolSource.fileNames = [
# 'file:out_me0_test.root'
# ]
