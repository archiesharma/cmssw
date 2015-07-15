import FWCore.ParameterSet.Config as cms

#Track selector
from Validation.RecoMuon.selectors_cff import *

#TrackAssociation
from SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi import *
import SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi
import SimTracker.TrackAssociation.TrackAssociatorByPosition_cfi

TrackAssociatorByHits = SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi.quickTrackAssociatorByHits.clone( ComponentName = 'TrackAssociatorByHits' )

OnlineTrackAssociatorByHits = SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi.quickTrackAssociatorByHits.clone()
OnlineTrackAssociatorByHits.ComponentName = 'OnlineTrackAssociatorByHits'
OnlineTrackAssociatorByHits.UseGrouped = cms.bool(False)
OnlineTrackAssociatorByHits.UseSplitting = cms.bool(False)
OnlineTrackAssociatorByHits.ThreeHitTracksAreSpecial = False

TrackAssociatorByPosDeltaR = SimTracker.TrackAssociation.TrackAssociatorByPosition_cfi.TrackAssociatorByPosition.clone()
TrackAssociatorByPosDeltaR.ComponentName = 'TrackAssociatorByDeltaR'
TrackAssociatorByPosDeltaR.method = cms.string('momdr')
TrackAssociatorByPosDeltaR.QCut = cms.double(0.5)
TrackAssociatorByPosDeltaR.ConsiderAllSimHits = cms.bool(True)

#
# Configuration for Muon track extractor
#

selectedVertices = cms.EDFilter("VertexSelector",
    src = cms.InputTag('offlinePrimaryVertices'),
    cut = cms.string("isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0 & abs(z) < 24 & abs(position.Rho) < 2."),
    filter = cms.bool(False)
)

selectedFirstPrimaryVertex = cms.EDFilter("PATSingleVertexSelector",
    mode = cms.string('firstVertex'),
    vertices = cms.InputTag('selectedVertices'),
    filter = cms.bool(False)
)

trackWithVertexSelector = cms.EDProducer("TrackWithVertexSelector",
    # -- module configuration --
    src = cms.InputTag('generalTracks'),
    quality = cms.string("highPurity"),
    useVtx = cms.bool(False),
    vertexTag = cms.InputTag('selectedFirstPrimaryVertex'),
    nVertices = cms.uint32(1),
    vtxFallback = cms.bool(True),
    copyExtras = cms.untracked.bool(False),
    copyTrajectories = cms.untracked.bool(False),
    # --------------------------
    # -- these are the vertex compatibility cuts --
    zetaVtx = cms.double(0.2),
    rhoVtx = cms.double(0.1),
    # ---------------------------------------------
    # -- dummy selection on tracks --
    etaMin = cms.double(0.0),
    etaMax = cms.double(5.0),
    ptMin = cms.double(0.00001),
    ptMax = cms.double(999999.),
    d0Max = cms.double(999999.),
    dzMax = cms.double(999999.),
    normalizedChi2 = cms.double(999999.),
    numberOfValidHits = cms.uint32(0),
    numberOfLostHits = cms.uint32(999),
    numberOfValidPixelHits = cms.uint32(0),
    ptErrorCut = cms.double(999999.)
    # ------------------------------
)

muonTRKSTA = cms.EDFilter("MuonSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("isTrackerMuon && isStandAloneMuon"),
    filter = cms.bool(False)
)

muonTRKSTAGLB = cms.EDFilter("MuonSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("isTrackerMuon && isStandAloneMuon && isGlobalMuon"),
    filter = cms.bool(False)
)

muonTRKSTA20 = cms.EDFilter("MuonSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("pt > 20 && isTrackerMuon && isStandAloneMuon"),
    filter = cms.bool(False)
)

muonPt3 = cms.EDFilter("MuonSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("pt > 3."),
    filter = cms.bool(False)
)

muonPt5 = cms.EDFilter("MuonSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("pt > 5."),
    filter = cms.bool(False)
)

muonPt10 = cms.EDFilter("MuonSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("pt > 10."),
    filter = cms.bool(False)
)

muonPt15 = cms.EDFilter("MuonSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("pt > 15."),
    filter = cms.bool(False)
)

muonPt20 = cms.EDFilter("MuonSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("pt > 20."),
    filter = cms.bool(False)
)

#-----------------------------------------------------------------------------------------------------------------------

import PhysicsTools.RecoAlgos.recoTrackSelector_cfi
staMuonsPt3 = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
staMuonsPt3.ptMin = cms.double(3.0)
staMuonsPt3.quality = cms.vstring('')
staMuonsPt3.minHit = cms.int32(0)
staMuonsPt3.src = cms.InputTag("standAloneMuons:UpdatedAtVtx")

staMuonsPt5 = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
staMuonsPt5.ptMin = cms.double(5.0)
staMuonsPt5.quality = cms.vstring('')
staMuonsPt5.minHit = cms.int32(0)
staMuonsPt5.src = cms.InputTag("standAloneMuons:UpdatedAtVtx")

staMuonsPt10 = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
staMuonsPt10.ptMin = cms.double(10.0)
staMuonsPt10.quality = cms.vstring('')
staMuonsPt10.minHit = cms.int32(0)
staMuonsPt10.src = cms.InputTag("standAloneMuons:UpdatedAtVtx")

staMuonsPt15 = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
staMuonsPt15.ptMin = cms.double(15.0)
staMuonsPt15.quality = cms.vstring('')
staMuonsPt15.minHit = cms.int32(0)
staMuonsPt15.src = cms.InputTag("standAloneMuons:UpdatedAtVtx")

staMuonsPt20 = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
staMuonsPt20.ptMin = cms.double(20.0)
staMuonsPt20.quality = cms.vstring('')
staMuonsPt20.minHit = cms.int32(0)
staMuonsPt20.src = cms.InputTag("standAloneMuons:UpdatedAtVtx")

#-----------------------------------------------------------------------------------------------------------------------

tevMuonsFirstHit20 = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
tevMuonsFirstHit20.ptMin = cms.double(20.0)
tevMuonsFirstHit20.quality = cms.vstring('')
tevMuonsFirstHit20.minHit = cms.int32(0)
tevMuonsFirstHit20.src = cms.InputTag("tevMuons:firstHit")

tevMuonsPicky20 = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
tevMuonsPicky20.ptMin = cms.double(20.0)
tevMuonsPicky20.quality = cms.vstring('')
tevMuonsPicky20.minHit = cms.int32(0)
tevMuonsPicky20.src = cms.InputTag("tevMuons:picky")

#-----------------------------------------------------------------------------------------------------------------------

#import SimMuon.MCTruth.MuonTrackProducer_cfi
#extractedGlobalMuons = SimMuon.MCTruth.MuonTrackProducer_cfi.muonTrackProducer.clone()
#extractedGlobalMuons.selectionTags = ('AllGlobalMuons',)
#extractedGlobalMuons.trackType = "globalTrack"

extractedTRKSTAMuons = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonTRKSTA"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuon')
)

extractedTRKSTAGLBMuons = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonTRKSTAGLB"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuon')
)

extractedTRKSTAMuons20 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonTRKSTA20"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuon')
)

extractedGlobalMuons = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muons"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('AllGlobalMuons'),
   trackType = cms.string('globalTrack')
)

extractedGlobalMuons3 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt3"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('AllGlobalMuons'),
   trackType = cms.string('globalTrack')
)

extractedGlobalMuons5 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt5"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('AllGlobalMuons'),
   trackType = cms.string('globalTrack')
)

extractedGlobalMuons10 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt10"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('AllGlobalMuons'),
   trackType = cms.string('globalTrack')
)

extractedGlobalMuons15 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt15"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('AllGlobalMuons'),
   trackType = cms.string('globalTrack')
)

extractedGlobalMuons20 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt20"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('AllGlobalMuons'),
   trackType = cms.string('globalTrack')
)

#-----------------------------------------------------------------------------------------------------------------------

bestMuon = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muons"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuon')
)

bestMuon3 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt3"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuon')
)

bestMuon5 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt5"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuon')
)

bestMuon10 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt10"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuon')
)

bestMuon15 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt15"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuon')
)

bestMuon20 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt20"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuon')
)

#-----------------------------------------------------------------------------------------------------------------------

bestMuonLoose = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muons"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonLoose')
)

bestMuonLoose3 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt3"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonLoose')
)

bestMuonLoose5 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt5"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonLoose')
)

bestMuonLoose10 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt10"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonLoose')
)

bestMuonLoose15 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt15"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonLoose')
)

bestMuonLoose20 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt20"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonLoose')
)

#-----------------------------------------------------------------------------------------------------------------------

bestMuonLoose2 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muons"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonLoose2')
)

bestMuonLoose32 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt3"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonLoose2')
)

bestMuonLoose52 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt5"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonLoose2')
)

bestMuonLoose102 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt10"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonLoose2')
)

bestMuonLoose152 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt15"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonLoose2')
)

bestMuonLoose202 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt20"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonLoose2')
)

#-----------------------------------------------------------------------------------------------------------------------

bestMuonTight = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muons"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonTight')
)

bestMuonTight3 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt3"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonTight')
)

bestMuonTight5 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt5"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonTight')
)

bestMuonTight10 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt10"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonTight')
)

bestMuonTight15 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt15"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonTight')
)

bestMuonTight20 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt20"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonTight')
)

#-----------------------------------------------------------------------------------------------------------------------

bestMuonTightNoIPz = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muons"),
   vxtTag = cms.InputTag("selectedVertices"),
   useIPxy = cms.untracked.bool(True),
   useIPz = cms.untracked.bool(False),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonTight')
)

bestMuonTightNoIPz3 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt3"),
   vxtTag = cms.InputTag("selectedVertices"),
   useIPxy = cms.untracked.bool(True),
   useIPz = cms.untracked.bool(False),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonTight')
)

bestMuonTightNoIPz5 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt5"),
   vxtTag = cms.InputTag("selectedVertices"),
   useIPxy = cms.untracked.bool(True),
   useIPz = cms.untracked.bool(False),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonTight')
)

bestMuonTightNoIPz10 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt10"),
   vxtTag = cms.InputTag("selectedVertices"),
   useIPxy = cms.untracked.bool(True),
   useIPz = cms.untracked.bool(False),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonTight')
)

bestMuonTightNoIPz15 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt15"),
   vxtTag = cms.InputTag("selectedVertices"),
   useIPxy = cms.untracked.bool(True),
   useIPz = cms.untracked.bool(False),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonTight')
)

bestMuonTightNoIPz20 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt20"),
   vxtTag = cms.InputTag("selectedVertices"),
   useIPxy = cms.untracked.bool(True),
   useIPz = cms.untracked.bool(False),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonTight')
)

#-----------------------------------------------------------------------------------------------------------------------

bestMuonTuneP = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muons"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('bestMuonTuneP')
)

#-----------------------------------------------------------------------------------------------------------------------

trackerMuons = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muons"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('AllTrackerMuons'),
   trackType = cms.string('bestMuon')
)

TMOneStationTight = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muons"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('TMOneStationTight'),
   trackType = cms.string('bestMuon')
)

TMLastStationAngTight = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muons"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('TMLastStationAngTight'),
   trackType = cms.string('bestMuon')
)

TrackerMuonArbitrated = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muons"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('TrackerMuonArbitrated'),
   trackType = cms.string('bestMuon')
)

trackerMuons5 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt5"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('AllTrackerMuons'),
   trackType = cms.string('bestMuon')
)

TMOneStationTight5 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt5"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('TMOneStationTight'),
   trackType = cms.string('bestMuon')
)

TMLastStationAngTight5 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt5"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('TMLastStationAngTight'),
   trackType = cms.string('bestMuon')
)

TrackerMuonArbitrated5 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt5"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('TrackerMuonArbitrated'),
   trackType = cms.string('bestMuon')
)

trackerMuons10 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt10"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('AllTrackerMuons'),
   trackType = cms.string('bestMuon')
)

TMOneStationTight10 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt10"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('TMOneStationTight'),
   trackType = cms.string('bestMuon')
)

TMLastStationAngTight10 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt10"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('TMLastStationAngTight'),
   trackType = cms.string('bestMuon')
)

TrackerMuonArbitrated10 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt10"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('TrackerMuonArbitrated'),
   trackType = cms.string('bestMuon')
)

trackerMuons20 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt20"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('AllTrackerMuons'),
   trackType = cms.string('bestMuon')
)

TMOneStationTight20 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt20"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('TMOneStationTight'),
   trackType = cms.string('bestMuon')
)

TMLastStationAngTight20 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt20"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('TMLastStationAngTight'),
   trackType = cms.string('bestMuon')
)

TrackerMuonArbitrated20 = cms.EDProducer("MuonTrackCollProducer",
   muonsTag = cms.InputTag("muonPt20"),
   vxtTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('TrackerMuonArbitrated'),
   trackType = cms.string('bestMuon')
)

#-----------------------------------------------------------------------------------------------------------------------

muonColl_seq = cms.Sequence(
			     #muonTRKSTA * 
			     #muonTRKSTAGLB * 
			     #muonTRKSTA20 *
			     muonPt5 * 
			     muonPt10 *
			     muonPt20)
trackColl_seq = cms.Sequence(staMuonsPt20)
extractedMuonTracks_seq = cms.Sequence( #extractedTRKSTAMuons * 
				  	#extractedTRKSTAGLBMuons *
					#extractedTRKSTAMuons20 *
				        extractedGlobalMuons * 
					extractedGlobalMuons5 * 
					extractedGlobalMuons10 * 
					extractedGlobalMuons20)
bestMuon_seq = cms.Sequence( 
			     #bestMuon * bestMuonLoose * bestMuonTight * bestMuonTightNoIPz * 
			     #bestMuon5 * bestMuonLoose5 * bestMuonTight5 * bestMuonTightNoIPz5 *
			     #bestMuon10 * bestMuonLoose10 * bestMuonTight10 * bestMuonTightNoIPz10 *
			     #bestMuon20 * bestMuonLoose20 * bestMuonTight20 * bestMuonTightNoIPz20
			     bestMuonLoose * bestMuonLoose5 *  bestMuonLoose10 * bestMuonLoose20 *
			     bestMuonTight * bestMuonTight5 * bestMuonTight10 * bestMuonTight20 *
			     bestMuonLoose2 * bestMuonLoose52  * bestMuonLoose102  * bestMuonLoose202 
			   )
bestMuonTuneP_seq = cms.Sequence( bestMuonTuneP )
trackerMuon_seq = cms.Sequence(
			       #trackerMuons * TrackerMuonArbitrated * TMOneStationTight * TMLastStationAngTight *
			       #trackerMuons5 * TrackerMuonArbitrated5 * TMOneStationTight5 * TMLastStationAngTight5 *
			       #trackerMuons10 * TrackerMuonArbitrated10 * TMOneStationTight10 * TMLastStationAngTight10 *
			       #trackerMuons20 * TrackerMuonArbitrated20 * TMOneStationTight20 * TMLastStationAngTight20
			       TrackerMuonArbitrated * TrackerMuonArbitrated5 * TrackerMuonArbitrated10 * TrackerMuonArbitrated20
		   	      )
#
# Configuration for Seed track extractor
#
#
#import SimMuon.MCTruth.SeedToTrackProducer_cfi
#seedsOfSTAmuons = SimMuon.MCTruth.SeedToTrackProducer_cfi.SeedToTrackProducer.clone()
#seedsOfSTAmuons.L2seedsCollection = cms.InputTag("ancientMuonSeed")
#seedsOfSTAmuons_seq = cms.Sequence( seedsOfSTAmuons )

# select probe tracks
#import PhysicsTools.RecoAlgos.recoTrackSelector_cfi
probeTracks = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
probeTracks.quality = cms.vstring('highPurity')
probeTracks.tip = cms.double(3.5)
probeTracks.lip = cms.double(30.)
probeTracks.ptMin = cms.double(4.0)
probeTracks.minRapidity = cms.double(-2.5)
probeTracks.maxRapidity = cms.double(2.5)
probeTracks_seq = cms.Sequence( probeTracks )

#
# Associators for Full Sim + Reco:
#

tpToTkmuTrackAssociation = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByHits'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
#    label_tr = cms.InputTag('generalTracks')
    label_tr = cms.InputTag('probeTracks')
)

tpToStaTrackAssociation = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('standAloneMuons','')
)

tpToStaUpdTrackAssociation = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('standAloneMuons','UpdatedAtVtx')
)

tpToGlbTrackAssociation = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('extractedGlobalMuons')
)

tpToStaSETTrackAssociation = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('standAloneSETMuons','')
)

tpToStaSETUpdTrackAssociation = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('standAloneSETMuons','UpdatedAtVtx')
)

tpToGlbSETTrackAssociation = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('globalSETMuons')
)

tpToTevFirstTrackAssociation = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('tevMuons','firstHit')
)

tpToTevPickyTrackAssociation = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('tevMuons','picky')
)
tpToTevDytTrackAssociation = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('tevMuons','dyt')
)

tpToL2TrackAssociation = cms.EDProducer('TrackAssociatorEDProducer',
    ignoremissingtrackcollection=cms.untracked.bool(True),
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('hltL2Muons','')
)

tpToL2UpdTrackAssociation = cms.EDProducer('TrackAssociatorEDProducer',
    ignoremissingtrackcollection=cms.untracked.bool(True),
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('hltL2Muons','UpdatedAtVtx')
)

tpToL3TrackAssociation = cms.EDProducer("TrackAssociatorEDProducer",
    ignoremissingtrackcollection=cms.untracked.bool(True),
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('hltL3Muons')
)

tpToL3TkTrackTrackAssociation = cms.EDProducer("TrackAssociatorEDProducer",
    ignoremissingtrackcollection=cms.untracked.bool(True),
    associator = cms.string('OnlineTrackAssociatorByHits'),
    label_tp = cms.InputTag('mix','MergedTrackTruth'),
    label_tr = cms.InputTag('hltL3TkTracksFromL2','')
)

tpToL3L2TrackTrackAssociation = cms.EDProducer("TrackAssociatorEDProducer",
    ignoremissingtrackcollection=cms.untracked.bool(True),
    associator = cms.string('OnlineTrackAssociatorByHits'),
    label_tp = cms.InputTag('mix','MergedTrackTruth'),
    label_tr = cms.InputTag('hltL3Muons:L2Seeded')
)



#MuonAssociation
import SimMuon.MCTruth.MuonAssociatorByHits_cfi

tpToTkMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
#tpToStaSeedAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
tpToStaMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
tpToStaUpdMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
tpToGlbMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
tpToStaRefitMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
tpToStaRefitUpdMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
tpToStaSETMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
tpToStaSETUpdMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
tpToGlbSETMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
tpToTevFirstMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
tpToTevPickyMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
tpToTevDytMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
tpToL3TkMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
tpToL2MuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
tpToL2UpdMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
tpToL3MuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()

# few more association modules usable for the Upgrade TP studies 
tpToTkSelMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
tpToStaUpdSelMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
tpToGlbSelMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()

tpToTkMuonAssociation.tpTag = 'mix:MergedTrackTruth'
#tpToTkMuonAssociation.tracksTag = 'generalTracks'
tpToTkMuonAssociation.tracksTag = 'probeTracks'
tpToTkMuonAssociation.UseTracker = True
tpToTkMuonAssociation.UseMuon = False

#tpToStaSeedAssociation.tpTag = 'mix:MergedTrackTruth'
#tpToStaSeedAssociation.tracksTag = 'seedsOfSTAmuons'
#tpToStaSeedAssociation.UseTracker = False
#tpToStaSeedAssociation.UseMuon = True
#

tpToStaMuonAssociation.tpTag = 'mix:MergedTrackTruth'
tpToStaMuonAssociation.tracksTag = 'standAloneMuons'
tpToStaMuonAssociation.UseTracker = False
tpToStaMuonAssociation.UseMuon = True

tpToStaUpdMuonAssociation.tpTag = 'mix:MergedTrackTruth'
tpToStaUpdMuonAssociation.tracksTag = 'standAloneMuons:UpdatedAtVtx'
tpToStaUpdMuonAssociation.UseTracker = False
tpToStaUpdMuonAssociation.UseMuon = True

tpToGlbMuonAssociation.tpTag = 'mix:MergedTrackTruth'
tpToGlbMuonAssociation.tracksTag = 'extractedGlobalMuons'
tpToGlbMuonAssociation.UseTracker = True
tpToGlbMuonAssociation.UseMuon = True

tpToStaRefitMuonAssociation.tpTag = 'mix:MergedTrackTruth'
tpToStaRefitMuonAssociation.tracksTag = 'refittedStandAloneMuons'
tpToStaRefitMuonAssociation.UseTracker = False
tpToStaRefitMuonAssociation.UseMuon = True

tpToStaRefitUpdMuonAssociation.tpTag = 'mix:MergedTrackTruth'
tpToStaRefitUpdMuonAssociation.tracksTag = 'refittedStandAloneMuons:UpdatedAtVtx'
tpToStaRefitUpdMuonAssociation.UseTracker = False
tpToStaRefitUpdMuonAssociation.UseMuon = True

tpToStaSETMuonAssociation.tpTag = 'mix:MergedTrackTruth'
tpToStaSETMuonAssociation.tracksTag = 'standAloneSETMuons'
tpToStaSETMuonAssociation.UseTracker = False
tpToStaSETMuonAssociation.UseMuon = True

tpToStaSETUpdMuonAssociation.tpTag = 'mix:MergedTrackTruth'
tpToStaSETUpdMuonAssociation.tracksTag = 'standAloneSETMuons:UpdatedAtVtx'
tpToStaSETUpdMuonAssociation.UseTracker = False
tpToStaSETUpdMuonAssociation.UseMuon = True

tpToGlbSETMuonAssociation.tpTag = 'mix:MergedTrackTruth'
tpToGlbSETMuonAssociation.tracksTag = 'globalSETMuons'
tpToGlbSETMuonAssociation.UseTracker = True
tpToGlbSETMuonAssociation.UseMuon = True

tpToTevFirstMuonAssociation.tpTag = 'mix:MergedTrackTruth'
tpToTevFirstMuonAssociation.tracksTag = 'tevMuons:firstHit'
tpToTevFirstMuonAssociation.UseTracker = True
tpToTevFirstMuonAssociation.UseMuon = True

tpToTevPickyMuonAssociation.tpTag = 'mix:MergedTrackTruth'
tpToTevPickyMuonAssociation.tracksTag = 'tevMuons:picky'
tpToTevPickyMuonAssociation.UseTracker = True
tpToTevPickyMuonAssociation.UseMuon = True

tpToTevDytMuonAssociation.tpTag = 'mix:MergedTrackTruth'
tpToTevDytMuonAssociation.tracksTag = 'tevMuons:dyt'
tpToTevDytMuonAssociation.UseTracker = True
tpToTevDytMuonAssociation.UseMuon = True

tpToL3TkMuonAssociation.tpTag = 'mix:MergedTrackTruth'
tpToL3TkMuonAssociation.tracksTag = 'hltL3TkTracksFromL2'
tpToL3TkMuonAssociation.DTrechitTag = 'hltDt1DRecHits'
tpToL3TkMuonAssociation.UseTracker = True
tpToL3TkMuonAssociation.UseMuon = False
tpToL3TkMuonAssociation.ignoreMissingTrackCollection = True
tpToL3TkMuonAssociation.UseSplitting = False
tpToL3TkMuonAssociation.UseGrouped = False

tpToL2MuonAssociation.tpTag = 'mix:MergedTrackTruth'
tpToL2MuonAssociation.tracksTag = 'hltL2Muons'
tpToL2MuonAssociation.DTrechitTag = 'hltDt1DRecHits'
tpToL2MuonAssociation.UseTracker = False
tpToL2MuonAssociation.UseMuon = True
tpToL2MuonAssociation.ignoreMissingTrackCollection = True

tpToL2UpdMuonAssociation.tpTag = 'mix:MergedTrackTruth'
tpToL2UpdMuonAssociation.tracksTag = 'hltL2Muons:UpdatedAtVtx'
tpToL2UpdMuonAssociation.DTrechitTag = 'hltDt1DRecHits'
tpToL2UpdMuonAssociation.UseTracker = False
tpToL2UpdMuonAssociation.UseMuon = True
tpToL2UpdMuonAssociation.ignoreMissingTrackCollection = True

tpToL3MuonAssociation.tpTag = 'mix:MergedTrackTruth'
tpToL3MuonAssociation.tracksTag = 'hltL3Muons'
tpToL3MuonAssociation.DTrechitTag = 'hltDt1DRecHits'
tpToL3MuonAssociation.UseTracker = True
tpToL3MuonAssociation.UseMuon = True
tpToL3MuonAssociation.ignoreMissingTrackCollection = True
tpToL3MuonAssociation.UseSplitting = False
tpToL3MuonAssociation.UseGrouped = False

# few more association modules usable for the Upgrade TP studies 
#tpToTkSelMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
#tpToStaUpdSelMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
#tpToGlbSelMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()

tpToTkSelMuonAssociation.tpTag = 'mix:MergedTrackTruth'
tpToTkSelMuonAssociation.tracksTag = 'probeTracks'
tpToTkSelMuonAssociation.UseTracker = True
tpToTkSelMuonAssociation.UseMuon = False
tpToTkSelMuonAssociation.EfficiencyCut_track = cms.double(0.5)
tpToTkSelMuonAssociation.PurityCut_track = cms.double(0.75)

tpToStaUpdSelMuonAssociation.tpTag = 'mix:MergedTrackTruth'
tpToStaUpdSelMuonAssociation.tracksTag = 'standAloneMuons:UpdatedAtVtx' 
tpToStaUpdSelMuonAssociation.UseTracker = False
tpToStaUpdSelMuonAssociation.UseMuon = True
tpToStaUpdSelMuonAssociation.includeZeroHitMuons = False

tpToGlbSelMuonAssociation.tpTag = 'mix:MergedTrackTruth'
tpToGlbSelMuonAssociation.tracksTag = 'extractedGlobalMuons' 
tpToGlbSelMuonAssociation.UseTracker = True
tpToGlbSelMuonAssociation.UseMuon = True
tpToGlbSelMuonAssociation.EfficiencyCut_track = cms.double(0.5)
tpToGlbSelMuonAssociation.PurityCut_track = cms.double(0.75)
tpToGlbSelMuonAssociation.acceptOneStubMatchings = False
tpToGlbSelMuonAssociation.includeZeroHitMuons = False
 
#
# Associators for cosmics:
#

tpToTkCosmicTrackAssociation = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByHits'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('ctfWithMaterialTracksP5LHCNavigation')
)

tpToStaCosmicTrackAssociation = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('cosmicMuons')
)

tpToGlbCosmicTrackAssociation = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('globalCosmicMuons')
)

tpToTkCosmicMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
tpToStaCosmicMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
tpToGlbCosmicMuonAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()

tpToTkCosmicMuonAssociation.tpTag = 'mix:MergedTrackTruth'
tpToTkCosmicMuonAssociation.tracksTag = 'ctfWithMaterialTracksP5LHCNavigation'
tpToTkCosmicMuonAssociation.UseTracker = True
tpToTkCosmicMuonAssociation.UseMuon = False

tpToStaCosmicMuonAssociation.tpTag = 'mix:MergedTrackTruth'
tpToStaCosmicMuonAssociation.tracksTag = 'cosmicMuons'
tpToStaCosmicMuonAssociation.UseTracker = False
tpToStaCosmicMuonAssociation.UseMuon = True

tpToGlbCosmicMuonAssociation.tpTag = 'mix:MergedTrackTruth'
tpToGlbCosmicMuonAssociation.tracksTag = 'globalCosmicMuons'
tpToGlbCosmicMuonAssociation.UseTracker = True
tpToGlbCosmicMuonAssociation.UseMuon = True


#
# The full-sim association sequences
#

muonAssociation_seq = cms.Sequence(
    extractedMuonTracks_seq
#    +seedsOfSTAmuons_seq
    +probeTracks_seq+(tpToTkMuonAssociation+tpToTkmuTrackAssociation)
#    +(tpToStaSeedAssociation+tpToStaMuonAssociation+tpToStaUpdMuonAssociation+tpToGlbMuonAssociation)
    +(tpToStaMuonAssociation+tpToStaUpdMuonAssociation+tpToGlbMuonAssociation)
#   +(tpToStaTrackAssociation+tpToStaUpdTrackAssociation+tpToGlbTrackAssociation)
#    
# few more association modules usable for the Upgrade TP studies 
    +(tpToTkSelMuonAssociation+tpToStaUpdSelMuonAssociation+tpToGlbSelMuonAssociation) 
)
muonAssociationTEV_seq = cms.Sequence(
    (tpToTevFirstMuonAssociation+tpToTevPickyMuonAssociation+tpToTevDytMuonAssociation)
#    +(tpToTevFirstTrackAssociation+tpToTevPickyTrackAssociation)
)
muonAssociationRefit_seq = cms.Sequence(
    (tpToStaRefitMuonAssociation+tpToStaRefitUpdMuonAssociation)
)
muonAssociationSET_seq = cms.Sequence(
    (tpToStaSETMuonAssociation+tpToStaSETUpdMuonAssociation+tpToGlbSETMuonAssociation)
#    +(tpToStaSETTrackAssociation+tpToStaSETUpdTrackAssociation+tpToGlbSETTrackAssociation)
)
muonAssociationCosmic_seq = cms.Sequence(
    (tpToTkCosmicMuonAssociation+tpToStaCosmicMuonAssociation+tpToGlbCosmicMuonAssociation)
#    +(tpToTkCosmicTrackAssociation+tpToStaCosmicTrackAssociation+tpToGlbCosmicTrackAssociation)
)
muonAssociationHLT_seq = cms.Sequence(
    (tpToL2MuonAssociation+tpToL2UpdMuonAssociation+tpToL3MuonAssociation+tpToL3TkMuonAssociation)
#    +(tpToL2TrackAssociation+tpToL2UpdTrackAssociation+tpToL3TrackAssociation+tpToL3TkTrackTrackAssociation)
)


#
# Associators for Fast Sim
#

tpToTkmuTrackAssociationFS = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByHits'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
#    label_tr = cms.InputTag('generalTracks')
    label_tr = cms.InputTag('probeTracks')
)

tpToStaTrackAssociationFS = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('standAloneMuons','')
)

tpToStaUpdTrackAssociationFS = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('standAloneMuons','UpdatedAtVtx')
)

tpToGlbTrackAssociationFS = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('extractedGlobalMuons')
)

tpToTevFirstTrackAssociationFS = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('tevMuons','firstHit')
)

tpToTevPickyTrackAssociationFS = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('tevMuons','picky')
)

tpToTevDytTrackAssociationFS = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('tevMuons','dyt')
)

tpToL2TrackAssociationFS = cms.EDProducer('TrackAssociatorEDProducer',
    ignoremissingtrackcollection=cms.untracked.bool(True),
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('hltL2Muons','')
)

tpToL2UpdTrackAssociationFS = cms.EDProducer('TrackAssociatorEDProducer',
    ignoremissingtrackcollection=cms.untracked.bool(True),
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('hltL2Muons','UpdatedAtVtx')
)

tpToL3TrackAssociationFS = cms.EDProducer("TrackAssociatorEDProducer",
    ignoremissingtrackcollection=cms.untracked.bool(True),
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mix', 'MergedTrackTruth'),
    label_tr = cms.InputTag('hltL3Muons')
)

tpToL3TkTrackTrackAssociationFS = cms.EDProducer("TrackAssociatorEDProducer",
    ignoremissingtrackcollection=cms.untracked.bool(True),
    associator = cms.string('OnlineTrackAssociatorByHits'),
    label_tp = cms.InputTag('mix','MergedTrackTruth'),
    label_tr = cms.InputTag('hltL3TkTracksFromL2','')
)

tpToL3L2TrackTrackAssociationFS = cms.EDProducer("TrackAssociatorEDProducer",
    ignoremissingtrackcollection=cms.untracked.bool(True),
    associator = cms.string('OnlineTrackAssociatorByHits'),
    label_tp = cms.InputTag('mix','MergedTrackTruth'),
    label_tr = cms.InputTag('hltL3Muons:L2Seeded')
)


#MuonAssociation
import SimMuon.MCTruth.MuonAssociatorByHits_cfi

baseMuonAssociatorFS = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
baseMuonAssociatorFS.tpTag = 'mix:MergedTrackTruth'
baseMuonAssociatorFS.UseTracker = True
baseMuonAssociatorFS.UseMuon = True
baseMuonAssociatorFS.simtracksTag = "famosSimHits"
baseMuonAssociatorFS.DTsimhitsTag  = "MuonSimHits:MuonDTHits"
baseMuonAssociatorFS.CSCsimHitsTag = "MuonSimHits:MuonCSCHits"
baseMuonAssociatorFS.RPCsimhitsTag = "MuonSimHits:MuonRPCHits"
baseMuonAssociatorFS.GEMsimhitsTag = "MuonSimHits:MuonGEMHits"
baseMuonAssociatorFS.simtracksXFTag = "mix:famosSimHits"
baseMuonAssociatorFS.DTsimhitsXFTag  = "mix:MuonSimHitsMuonDTHits"
baseMuonAssociatorFS.CSCsimHitsXFTag = "mix:MuonSimHitsMuonCSCHits"
baseMuonAssociatorFS.RPCsimhitsXFTag = "mix:MuonSimHitsMuonRPCHits"
baseMuonAssociatorFS.GEMsimhitsXFTag = "mix:MuonSimHitsMuonGEMHits"
baseMuonAssociatorFS.ROUList = ['famosSimHitsTrackerHits']


tpToTkMuonAssociationFS  = baseMuonAssociatorFS.clone()
#tpToStaSeedAssociationFS = baseMuonAssociatorFS.clone()
tpToStaMuonAssociationFS  = baseMuonAssociatorFS.clone()
tpToStaUpdMuonAssociationFS  = baseMuonAssociatorFS.clone()
tpToStaRefitMuonAssociationFS  = baseMuonAssociatorFS.clone()
tpToStaRefitUpdMuonAssociationFS  = baseMuonAssociatorFS.clone()
tpToGlbMuonAssociationFS  = baseMuonAssociatorFS.clone()
tpToTevFirstMuonAssociationFS = baseMuonAssociatorFS.clone()
tpToTevPickyMuonAssociationFS = baseMuonAssociatorFS.clone()
tpToTevDytMuonAssociationFS = baseMuonAssociatorFS.clone()
tpToL3TkMuonAssociationFS = baseMuonAssociatorFS.clone()
tpToL2MuonAssociationFS   = baseMuonAssociatorFS.clone()
tpToL2UpdMuonAssociationFS   = baseMuonAssociatorFS.clone()
tpToL3MuonAssociationFS   = baseMuonAssociatorFS.clone()

#tpToTkMuonAssociationFS.tracksTag = 'generalTracks'
tpToTkMuonAssociationFS.tracksTag = 'probeTracks'
tpToTkMuonAssociationFS.tpTag = 'mix:MergedTrackTruth'
tpToTkMuonAssociationFS.UseTracker = True
tpToTkMuonAssociationFS.UseMuon = False

#tpToStaSeedAssociationFS.tpTag = 'mix:MergedTrackTruth'
#tpToStaSeedAssociationFS.tracksTag = 'seedsOfSTAmuons'
#tpToStaSeedAssociationFS.UseTracker = False
#tpToStaSeedAssociationFS.UseMuon = True

tpToStaMuonAssociationFS.tpTag = 'mix:MergedTrackTruth'
tpToStaMuonAssociationFS.tracksTag = 'standAloneMuons'
tpToStaMuonAssociationFS.UseTracker = False
tpToStaMuonAssociationFS.UseMuon = True

tpToStaUpdMuonAssociationFS.tracksTag = 'standAloneMuons:UpdatedAtVtx'
tpToStaUpdMuonAssociationFS.UseTracker = False
tpToStaUpdMuonAssociationFS.UseMuon = True

tpToStaRefitMuonAssociationFS.tracksTag = 'refittedStandAloneMuons'
tpToStaRefitMuonAssociationFS.UseTracker = False
tpToStaRefitMuonAssociationFS.UseMuon = True

tpToStaRefitUpdMuonAssociationFS.tracksTag = 'refittedStandAloneMuons:UpdatedAtVtx'
tpToStaRefitUpdMuonAssociationFS.UseTracker = False
tpToStaRefitUpdMuonAssociationFS.UseMuon = True

tpToGlbMuonAssociationFS.tracksTag = 'extractedGlobalMuons'
tpToGlbMuonAssociationFS.tpTag = 'mix:MergedTrackTruth'
tpToGlbMuonAssociationFS.UseTracker = True
tpToGlbMuonAssociationFS.UseMuon = True

tpToTevFirstMuonAssociationFS.tracksTag = 'tevMuons:firstHit'
tpToTevFirstMuonAssociationFS.tpTag = 'mix:MergedTrackTruth'
tpToTevFirstMuonAssociationFS.UseTracker = True
tpToTevFirstMuonAssociationFS.UseMuon = True

tpToTevPickyMuonAssociationFS.tracksTag = 'tevMuons:picky'
tpToTevPickyMuonAssociationFS.tpTag = 'mix:MergedTrackTruth'
tpToTevPickyMuonAssociationFS.UseTracker = True
tpToTevPickyMuonAssociationFS.UseMuon = True

tpToTevDytMuonAssociationFS.tracksTag = 'tevMuons:dyt'
tpToTevDytMuonAssociationFS.tpTag = 'mix:MergedTrackTruth'
tpToTevDytMuonAssociationFS.UseTracker = True
tpToTevDytMuonAssociationFS.UseMuon = True

tpToL3TkMuonAssociationFS.tracksTag = 'hltL3TkTracksFromL2'
tpToL3TkMuonAssociationFS.tpTag = 'mix:MergedTrackTruth'
tpToL3TkMuonAssociationFS.UseTracker = True
tpToL3TkMuonAssociationFS.UseMuon = False
tpToL3TkMuonAssociationFS.ignoreMissingTrackCollection = True
tpToL3TkMuonAssociationFS.UseSplitting = False
tpToL3TkMuonAssociationFS.UseGrouped = False

tpToL2MuonAssociationFS.tracksTag = 'hltL2Muons'
tpToL2MuonAssociationFS.UseTracker = False
tpToL2MuonAssociationFS.UseMuon = True
tpToL2MuonAssociationFS.ignoreMissingTrackCollection = True

tpToL2UpdMuonAssociationFS.tracksTag = 'hltL2Muons:UpdatedAtVtx'
tpToL2UpdMuonAssociationFS.UseTracker = False
tpToL2UpdMuonAssociationFS.UseMuon = True
tpToL2UpdMuonAssociationFS.ignoreMissingTrackCollection = True

tpToL3MuonAssociationFS.tracksTag = 'hltL3Muons'
tpToL3MuonAssociationFS.tpTag = 'mix:MergedTrackTruth'
tpToL3MuonAssociationFS.UseTracker = True
tpToL3MuonAssociationFS.UseMuon = True
tpToL3MuonAssociationFS.ignoreMissingTrackCollection = True
tpToL3MuonAssociationFS.UseSplitting = False
tpToL3MuonAssociationFS.UseGrouped = False



muonAssociationFastSim_seq = cms.Sequence(
        extractedMuonTracks_seq
#        +seedsOfSTAmuons_seq
        +probeTracks+(tpToTkMuonAssociationFS+tpToTkmuTrackAssociationFS) 
#        +(tpToStaSeedAssociationFS+tpToStaMuonAssociationFS+tpToStaUpdMuonAssociationFS+tpToGlbMuonAssociationFS)
        +(tpToStaMuonAssociationFS+tpToStaUpdMuonAssociationFS+tpToGlbMuonAssociationFS)
        +(tpToStaRefitMuonAssociationFS+tpToStaRefitUpdMuonAssociationFS)
        +(tpToTevFirstMuonAssociationFS+tpToTevPickyMuonAssociationFS+tpToTevDytMuonAssociationFS)
#        +tpToStaTrackAssociationFS+tpToStaUpdTrackAssociationFS+tpToGlbTrackAssociationFS
#        +tpToTevFirstTrackAssociationFS+tpToTevPickyTrackAssociationFS
        )
muonAssociationHLTFastSim_seq = cms.Sequence(
    tpToL2MuonAssociationFS+tpToL2UpdMuonAssociationFS+tpToL3MuonAssociationFS+tpToL3TkMuonAssociationFS
#    +tpToL2TrackAssociationFS+tpToL2UpdTrackAssociationFS+tpToL3TrackAssociationFS+tpToL3TkTrackTrackAssociationFS
    )

