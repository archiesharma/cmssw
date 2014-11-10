// -*- C++ -*-
//
// Package:    MuonMCClassifier
// Class:      MuonMCClassifier
//
/**\class MuonMCClassifier MuonMCClassifier.cc MuonAnalysis/MuonAssociators/src/MuonMCClassifier.cc


 CLASSIFICATION: For each RECO Muon, match to SIM particle, and then:
 - If the SIM is not a Muon, label as Punchthrough (1)
 - If the SIM is a Muon, then look at it's provenance.
 A) the SIM muon is also a GEN muon, whose parent is NOT A HADRON AND NOT A TAU
 -> classify as "primary" (4).
 B) the SIM muon is also a GEN muon, whose parent is HEAVY FLAVOURED HADRON OR A TAU
 -> classify as "heavy flavour" (3)
 C) classify as "light flavour/decay" (2)

 In any case, if the TP is not preferentially matched back to the same RECO muon,
 label as Ghost (flip the classification)


 FLAVOUR:
 - for non-muons: 0
 - for primary muons: 13
 - for non primary muons: flavour of the mother: abs(pdgId) of heaviest quark, or 15 for tau

 */
//
// Original Author:  Nov 16 16:12 (lxplus231.cern.ch)
//         Created:  Sun Nov 16 16:14:09 CET 2008
// $Id: MuonMCClassifier.cc,v 1.9 2013/06/24 12:53:19 speer Exp $
//
//
// system include files
#include <memory>
#include <set>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/Association.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimMuon/MCTruth/interface/MuonAssociatorByHits.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include <SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h>

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

//rumi
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

//
// class decleration
class MuonMCClassifier: public edm::EDProducer
{
public:
  explicit MuonMCClassifier(const edm::ParameterSet&);
//  std::vector<double> findSimVtx(edm::Event& iEvent);
//  bool isSoft(edm::Event& iEvent, reco::Muon myMuon, bool useIPxy, bool useIPz);
//  bool isTight(edm::Event& iEvent, reco::Muon myMuon, bool useIPxy, bool useIPz);
  ~MuonMCClassifier();
  void myClassification(size_t myNmu, edm::Handle<edm::View<reco::Muon> > myMuons_handle,
      edm::RefToBaseVector<reco::Muon> selMuons, MuonAssociatorByHits::MuonToSimCollection myRecSimColl,
      MuonAssociatorByHits::SimToMuonCollection mySimRecColl, std::vector<int> &myClassif, std::vector<int> &myExt,
      std::vector<int> &myHitsPdgId, std::vector<int> &myMomPdgId, std::vector<int> &myGmomPdgId,
      std::vector<int> &myMomStatus, std::vector<int> &myFlav, std::vector<int> &myMomFlav,
      std::vector<int> &myGmomFlav, std::vector<int> &myHmomFlav, std::vector<int> &myTpId,
      std::vector<float> &myProdRho, std::vector<float> &myProdZ, std::vector<float> &myMomRho,
      std::vector<float> &myMomZ, std::vector<float> &myTpAssoQuality,
      std::auto_ptr<reco::GenParticleCollection> &secondaries, std::vector<int> &muToPrimary,
      std::vector<int> &muToSecondary, edm::Handle < reco::GenParticleCollection > &genParticles);

private:
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  bool debug;

  /// The RECO objects
  edm::InputTag trackingParticles_;	 /// The TrackingParticle objects
  edm::InputTag muons_;
  bool linkToGenParticles_;		/// Create a link to the generator level particles
  edm::InputTag genParticles_;
  std::string associatorLabel_;		 /// The Associations
//not used  bool hasMuonCut_;			/// A preselection cut for the muon
//not used  StringCutObjectSelector<reco::Muon> muonCut_;
  double decayRho_, decayAbsZ_;		/// Cylinder to use to decide if a decay is early or late

  /// Track to use
//  MuonAssociatorByHits::MuonTrackType trackType_;
//  bool useIPxy, useIPz;
  /// Returns the flavour given a pdg id code
  int flavour(int pdgId) const;
  double muonPtCut_;
  double muonMinEtaCut_;
  double muonMaxEtaCut_;

  /// Write a ValueMap<int> in the event
  template<typename T>
  void writeValueMap(edm::Event &iEvent, const edm::Handle<edm::View<reco::Muon> > & handle,
      const std::vector<T> & values, const std::string & label) const;

//Get Mother of the tracking particle
  TrackingParticleRef getTpMother(TrackingParticleRef tp)
  {
    if (tp.isNonnull() && tp->parentVertex().isNonnull() && !tp->parentVertex()->sourceTracks().empty())
    {
      return tp->parentVertex()->sourceTracks()[0];
    }
    else
    {
      return TrackingParticleRef();
    }
  }

  /// Convert TrackingParticle into GenParticle, save into output collection,
  /// if mother is primary set reference to it,
  /// return index in output collection
  int convertAndPush(const TrackingParticle &tp, reco::GenParticleCollection &out, const TrackingParticleRef &momRef,
      const edm::Handle<reco::GenParticleCollection> & genParticles) const;

};

/*
std::vector<double> MuonMCClassifier::findSimVtx(edm::Event& iEvent)
{
  edm::Handle < reco::GenParticleCollection > genParticles;
  iEvent.getByLabel("genParticles", genParticles);
  std::vector<double> vtxCoord;
  vtxCoord.push_back(0);
  vtxCoord.push_back(0);
  vtxCoord.push_back(0);
  vtxCoord.push_back(0);
  vtxCoord.push_back(0);
  vtxCoord.push_back(0);
  vtxCoord.push_back(0);

  if (genParticles.isValid())
  {
    for (reco::GenParticleCollection::const_iterator itg = genParticles->begin(); itg != genParticles->end(); ++itg)
    {
      int id = itg->pdgId();
      int status = itg->status();

      if (abs(id) == 23 && status == 3)
      {
        vtxCoord[0] = 1;
        vtxCoord[4] = (double) (itg->vx());
        vtxCoord[5] = (double) (itg->vy());
        vtxCoord[6] = (double) (itg->vz());
      }

      if (abs(id) == 13 && status == 1)
      {
        vtxCoord[1] = (double) (itg->vx());
        vtxCoord[2] = (double) (itg->vy());
        vtxCoord[3] = (double) (itg->vz());
      }
    }
  }

  return vtxCoord;
}

//for standalone mouns - it is enough to ask whether the muon is standalone, no special function is needed

bool MuonMCClassifier::isSoft(edm::Event& iEvent, reco::Muon myMuon, bool useIPxy, bool useIPz)
{
  bool result = false;
  if (myMuon.muonBestTrack().isNonnull() && myMuon.innerTrack().isNonnull())
  {

    edm::Handle < reco::VertexCollection > vertexHandle;
    if(vertexHandle.isValid())
    {
    iEvent.getByLabel("selectedVertices", vertexHandle);
    const reco::VertexCollection* vertices = vertexHandle.product();

    bool isGood = muon::isGoodMuon(myMuon, muon::TMOneStationTight);
    bool trkLayMeas = myMuon.muonBestTrack()->hitPattern().trackerLayersWithMeasurement() > 5;
    bool pxlLayMeas = myMuon.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0;
    bool quality = myMuon.innerTrack()->quality(reco::Track::highPurity);
    bool ipxy = false;
    bool ipz = false;
    if (vertices->size() != 0 && useIPxy)
      ipxy = fabs(myMuon.muonBestTrack()->dxy((*vertices)[0].position())) < 0.2;
    else
      ipxy = true;
    if (vertices->size() != 0 && useIPz)
      ipz = fabs(myMuon.muonBestTrack()->dz((*vertices)[0].position())) < 0.5;
    else
      ipz = true;
//    if (isGood && trkLayMeas && pxlLayMeas && quality && ipxy && ipz)
    if (isGood && trkLayMeas && pxlLayMeas && quality)
      result = true;
    if(isGood && trkLayMeas && pxlLayMeas && quality && ipxy && ipz) std::cout << "old is good" << std::endl;
    else std::cout << "soft is " << result << std::endl;
  }
  }
  return result;
}

bool MuonMCClassifier::isTight(edm::Event& iEvent, reco::Muon myMuon, bool useIPxy, bool useIPz)
{
  bool result = false;

  if (myMuon.muonBestTrack().isNonnull() && myMuon.innerTrack().isNonnull() && myMuon.globalTrack().isNonnull())
  {
    std::vector<double> vtxCoord = findSimVtx(iEvent);
    GlobalPoint point(vtxCoord[1], vtxCoord[2], vtxCoord[3]);
    GlobalPoint pointDY(vtxCoord[4], vtxCoord[5], vtxCoord[6]);

    double muonZ = pointDY.z();
    edm::Handle < reco::VertexCollection > vertexHandle;
    if(vertexHandle.isValid())
    {
    iEvent.getByLabel("selectedVertices", vertexHandle);
    const reco::VertexCollection* vertices = vertexHandle.product();

    double distInit = 24;
    int indexFinal = 0;
    for (int i = 0; i < (int) vertices->size(); i++)
    {
      double vtxZ = (*vertices)[i].z();
      double dist = fabs(muonZ - vtxZ);
      if (dist < distInit)
      {
        distInit = dist;
        indexFinal = i;
      }
    }

    double ipxySim = 999;
    double ipzSim = 999;

    if (vtxCoord[0] < 0.5)
    {
      ipxySim = fabs(myMuon.muonBestTrack()->dxy(math::XYZPoint(point.x(), point.y(), point.z())));
      ipzSim = fabs(myMuon.muonBestTrack()->dz(math::XYZPoint(point.x(), point.y(), point.z())));
    }
    else if (vtxCoord[0] > 0.5)
    {
      ipxySim = fabs(myMuon.muonBestTrack()->dxy(math::XYZPoint(pointDY.x(), pointDY.y(), pointDY.z())));
      ipzSim = fabs(myMuon.muonBestTrack()->dz(math::XYZPoint(pointDY.x(), pointDY.y(), pointDY.z())));
    }

    bool ipxySimBool = ipxySim < 0.2;
    bool ipzSimBool = ipzSim < 0.5;

    bool trkLayMeas = myMuon.muonBestTrack()->hitPattern().trackerLayersWithMeasurement() > 5;
    bool isGlb = myMuon.isGlobalMuon();
    bool isPF = myMuon.isPFMuon();
    bool chi2 = myMuon.globalTrack()->normalizedChi2() < 10.;
    bool validHits = myMuon.globalTrack()->hitPattern().numberOfValidMuonHits() > 0;
    bool matchedSt = myMuon.numberOfMatchedStations() > 1;
    bool ipxy = false;
    bool ipz = false;

    if (vertices->size() != 0 && useIPxy == true)
    {

      if (vtxCoord[0] > 0.5)
        ipxy = fabs(myMuon.muonBestTrack()->dxy((*vertices)[indexFinal].position())) < 0.2;
      else
        ipxy = ipxySimBool;
    }

    else if (vertices->size() == 0 && useIPxy == true)
      ipxy = false;
    else if (useIPxy == false)
      ipxy = true;

    if (vertices->size() != 0 && useIPz == true)
    {

      if (vtxCoord[0] > 0.5)
        ipz = fabs(myMuon.muonBestTrack()->dz((*vertices)[indexFinal].position())) < 0.5;
      else
        ipz = ipzSimBool;

    }
    else if (vertices->size() == 0 && useIPz == true)
      ipz = false;
    else if (useIPz == false)
      ipz = true;

    bool validPxlHit = myMuon.innerTrack()->hitPattern().pixelLayersWithMeasurement(3, 2) > 0;

//    if (trkLayMeas && isGlb && isPF && chi2 && validHits && matchedSt && ipxy && ipz && validPxlHit)
     if (trkLayMeas && isGlb && isPF && chi2 && validHits && matchedSt && validPxlHit)
      result = true;
     if(trkLayMeas && isGlb && isPF && chi2 && validHits && matchedSt && ipxy && ipz && validPxlHit) std::cout << "old is good" << std::endl;
     else std::cout << "soft is " << result << std::endl;
  }
  }
  return result;
}
*/

MuonMCClassifier::MuonMCClassifier(const edm::ParameterSet &iConfig) :
    trackingParticles_(iConfig.getParameter < edm::InputTag > ("simLabel")) //
        , muons_(iConfig.getParameter < edm::InputTag > ("muonLabel"))  //
        , linkToGenParticles_(iConfig.getParameter<bool>("linkToGenParticles"))  //
        , genParticles_(linkToGenParticles_ ? iConfig.getParameter < edm::InputTag > ("genParticles") : edm::InputTag("NONE"))  //
        , associatorLabel_(iConfig.getParameter < std::string > ("muAssocLabel")) //
//        , hasMuonCut_(iConfig.existsAs < std::string > ("muonPreselection"))  //
//        , muonCut_(hasMuonCut_ ? iConfig.getParameter < std::string > ("muonPreselection") : "")  //
        , decayRho_(iConfig.getParameter<double>("decayRho")) //
        , decayAbsZ_(iConfig.getParameter<double>("decayAbsZ"))   //
//        , useIPxy(iConfig.getUntrackedParameter<bool>("useIPxy", true)) //
//        , useIPz(iConfig.getUntrackedParameter<bool>("useIPz", true)) //
        , muonPtCut_(iConfig.getParameter<double>("muonPtCut"))  //
        , muonMinEtaCut_(iConfig.getParameter<double>("muonMinEtaCut")) //
        , muonMaxEtaCut_(iConfig.getParameter<double>("muonMaxEtaCut"))
{
/*
  std::string trackType = iConfig.getParameter < std::string > ("trackType");
  if (trackType == "innerTrack")
    trackType_ = MuonAssociatorByHits::InnerTk;
  else if (trackType == "outerTrack")
    trackType_ = MuonAssociatorByHits::OuterTk;
  else if (trackType == "globalTrack")
    trackType_ = MuonAssociatorByHits::GlobalTk;
  else if (trackType == "innerTrackPlusSegments")
    trackType_ = MuonAssociatorByHits::Segments;
  else
    throw cms::Exception("Configuration") << "Track type '" << trackType << "' not supported.\n";
*/

  debug = false;
  //glob tight
  produces < edm::ValueMap<int> > ("glbT");
  produces < edm::ValueMap<int> > ("extglbT");
  produces < edm::ValueMap<int> > ("flavglbT");
  produces < edm::ValueMap<int> > ("hitsPdgIdglbT");
  produces < edm::ValueMap<int> > ("momPdgIdglbT");
  produces < edm::ValueMap<int> > ("momFlavglbT");
  produces < edm::ValueMap<int> > ("momStatusglbT");
  produces < edm::ValueMap<int> > ("gmomPdgIdglbT");
  produces < edm::ValueMap<int> > ("gmomFlavglbT");
  produces < edm::ValueMap<int> > ("hmomFlavglbT"); // heaviest mother flavour
  produces < edm::ValueMap<int> > ("tpIdglbT");
  produces < edm::ValueMap<float> > ("prodRhoglbT");
  produces < edm::ValueMap<float> > ("prodZglbT");
  produces < edm::ValueMap<float> > ("momRhoglbT");
  produces < edm::ValueMap<float> > ("momZglbT");
  produces < edm::ValueMap<float> > ("tpAssoQualityglbT");

  produces < edm::ValueMap<int> > ("glbL");
  produces < edm::ValueMap<int> > ("extglbL");
  produces < edm::ValueMap<int> > ("flavglbL");
  produces < edm::ValueMap<int> > ("hitsPdgIdglbL");
  produces < edm::ValueMap<int> > ("momPdgIdglbL");
  produces < edm::ValueMap<int> > ("momFlavglbL");
  produces < edm::ValueMap<int> > ("momStatusglbL");
  produces < edm::ValueMap<int> > ("gmomPdgIdglbL");
  produces < edm::ValueMap<int> > ("gmomFlavglbL");
  produces < edm::ValueMap<int> > ("hmomFlavglbL"); // heaviest mother flavour
  produces < edm::ValueMap<int> > ("tpIdglbL");
  produces < edm::ValueMap<float> > ("prodRhoglbL");
  produces < edm::ValueMap<float> > ("prodZglbL");
  produces < edm::ValueMap<float> > ("momRhoglbL");
  produces < edm::ValueMap<float> > ("momZglbL");
  produces < edm::ValueMap<float> > ("tpAssoQualityglbL");

  produces < edm::ValueMap<int> > ("sta");
  produces < edm::ValueMap<int> > ("extsta");
  produces < edm::ValueMap<int> > ("flavsta");
  produces < edm::ValueMap<int> > ("hitsPdgIdsta");
  produces < edm::ValueMap<int> > ("momPdgIdsta");
  produces < edm::ValueMap<int> > ("momFlavsta");
  produces < edm::ValueMap<int> > ("momStatussta");
  produces < edm::ValueMap<int> > ("gmomPdgIdsta");
  produces < edm::ValueMap<int> > ("gmomFlavsta");
  produces < edm::ValueMap<int> > ("hmomFlavsta"); // heaviest mother flavour
  produces < edm::ValueMap<int> > ("tpIdsta");
  produces < edm::ValueMap<float> > ("prodRhosta");
  produces < edm::ValueMap<float> > ("prodZsta");
  produces < edm::ValueMap<float> > ("momRhosta");
  produces < edm::ValueMap<float> > ("momZsta");
  produces < edm::ValueMap<float> > ("tpAssoQualitysta");

  if (linkToGenParticles_)
  {
//    produces < reco::GenParticleCollection > ("secondaries");
    produces < edm::Association<reco::GenParticleCollection> > ("toPrimaries");
    produces < edm::Association<reco::GenParticleCollection> > ("toSecondaries");

/*
//for glb tight
    produces < reco::GenParticleCollection > ("secondariesforGlbT");
    produces < edm::Association<reco::GenParticleCollection> > ("toPrimariesforGlbT");
    produces < edm::Association<reco::GenParticleCollection> > ("toSecondariesforGlbT");
//for glb soft
    produces < reco::GenParticleCollection > ("secondariesforGlbL");
    produces < edm::Association<reco::GenParticleCollection> > ("toPrimariesforGlbL");
    produces < edm::Association<reco::GenParticleCollection> > ("toSecondariesforGlbL");
//for standalone
    produces < reco::GenParticleCollection > ("secondariesforSta");
    produces < edm::Association<reco::GenParticleCollection> > ("toPrimariesforSta");
    produces < edm::Association<reco::GenParticleCollection> > ("toSecondariesforSta");
*/
  }
}

MuonMCClassifier::~MuonMCClassifier()
{
}

void MuonMCClassifier::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace reco;

std::cout << "file without ID cuts!" << std::endl;

  if (debug)
    edm::LogVerbatim("MuonMCClassifier") << "\n sono in MuonMCClassifier !";

  edm::Handle < edm::View<reco::Muon> > myMuons;
  iEvent.getByLabel(muons_, myMuons);

  edm::Handle < TrackingParticleCollection > trackingParticles;
  iEvent.getByLabel(trackingParticles_, trackingParticles);

  edm::Handle < reco::GenParticleCollection > genParticlesforGlbT;
  edm::Handle < reco::GenParticleCollection > genParticlesforGlbL;
  edm::Handle < reco::GenParticleCollection > genParticlesforSta;
//  edm::Handle < std::vector<int> > genBarcodes;
  if (linkToGenParticles_)
  {
    iEvent.getByLabel(genParticles_, genParticlesforGlbT);
    iEvent.getByLabel(genParticles_, genParticlesforGlbL);
    iEvent.getByLabel(genParticles_, genParticlesforSta);
    //   iEvent.getByLabel(genParticles_, genBarcodes);
  }

  edm::ESHandle < TrackAssociatorBase > associatorBase;
  iSetup.get<TrackAssociatorRecord>().get(associatorLabel_, associatorBase);
  const MuonAssociatorByHits * assoByHits = dynamic_cast<const MuonAssociatorByHits *>(associatorBase.product());
  if (assoByHits == 0)
    throw cms::Exception("Configuration") << "The Track Associator with label '" << associatorLabel_
        << "' is not a MuonAssociatorByHits.\n";

  MuonAssociatorByHits::MuonToSimCollection recSimColl;
  MuonAssociatorByHits::SimToMuonCollection simRecColl;
/*
  if (debug)
    edm::LogVerbatim("MuonMCClassifier") << "\n ***************************************************************** ";
  if (debug)
    edm::LogVerbatim("MuonMCClassifier") << " RECO MUON association, type:  " << trackType_;
  if (debug)
    edm::LogVerbatim("MuonMCClassifier") << " ***************************************************************** \n";
*/
  edm::RefToBaseVector < reco::Muon > selTightMuons;
  edm::RefToBaseVector < reco::Muon > selSoftMuons;
  edm::RefToBaseVector < reco::Muon > selStaMuons;

    for (size_t i = 0, n = myMuons->size(); i < n; ++i)
    {

     if ((*myMuons)[i].pt() > muonPtCut_ && fabs((*myMuons)[i].eta()) < muonMaxEtaCut_)

      //if ((*myMuons)[i].pt() > muonPtCut_ && fabs((*myMuons)[i].eta()) >= muonMinEtaCut_
        //  && fabs((*myMuons)[i].eta()) <= muonMaxEtaCut_)
      {
        edm::RefToBase < reco::Muon > rmu = myMuons->refAt(i);
        //if (muonCut_((*myMuons)[i]))
       // selMuons.push_back(rmu);
        if (debug)
          std::cout << "event: run\t" << (iEvent.id()).run() << "\tevent #\t" << (iEvent.id()).event() << std::endl; //
        if (debug)
          std::cout << "muon[" << i << "] with " << "muon pt = " << (*myMuons)[i].pt() << "\t and eta = "
              << (*myMuons)[i].eta() << std::endl; //
        if (debug)
          std::cout << "muon[" << i << "] isTracker = " << (*myMuons)[i].isTrackerMuon() << "\t isGlobal = "
              << (*myMuons)[i].isGlobalMuon() << "\t isStandalone = " << (*myMuons)[i].isStandAloneMuon() << std::endl;

//        bool soft = isSoft(iEvent, (*myMuons)[i], useIPxy, useIPz);
//        if (soft)
        if(muon::isGoodMuon((*myMuons)[i], muon::TMOneStationTight))
        {
          selSoftMuons.push_back(rmu);
          std::cout << "I am Good muon!" << std::endl;
        }

//        bool tight = isTight(iEvent, (*myMuons)[i], useIPxy, useIPz);
//        if (tight)
        if((*myMuons)[i].isGlobalMuon())
        {
          selTightMuons.push_back(rmu);
          std::cout << "I am Global muon!" << std::endl;
        }
        if ((*myMuons)[i].isStandAloneMuon())
          selStaMuons.push_back(rmu);
      }
    }

//create tracking particle collection for all the muons
  edm::RefVector < TrackingParticleCollection > allTPs;
  for (size_t i = 0, n = trackingParticles->size(); i < n; ++i)
  {
    allTPs.push_back(TrackingParticleRef(trackingParticles, i));
    if (debug)
      std::cout << "In trackingParticles" << std::endl;
  }

//create reco to sim and sim to reco associations
//global tight muons assotiated by global track
  MuonAssociatorByHits::MuonToSimCollection recSimColl_glbTight;
  MuonAssociatorByHits::SimToMuonCollection simRecColl_glbTight;
  edm::LogVerbatim("MuonMCClassifier") << "\n Global Tigh Muon association by global track";  //
  assoByHits->associateMuons(recSimColl_glbTight, simRecColl_glbTight, selTightMuons, MuonAssociatorByHits::GlobalTk, allTPs, &iEvent, &iSetup);  //

  //global loose muon assosiated by TMOneStationTight - require one well matched segment
  MuonAssociatorByHits::MuonToSimCollection recSimColl_glbloose;
  MuonAssociatorByHits::SimToMuonCollection simRecColl_glbloose;
  edm::LogVerbatim("MuonMCClassifier") << "\n Global Loose Muon association by inner track ";  //
  assoByHits->associateMuons(recSimColl_glbloose, simRecColl_glbloose, selSoftMuons, MuonAssociatorByHits::Segments, allTPs, &iEvent, &iSetup);  //

  //standalone loose muon assosiated by outher track
  MuonAssociatorByHits::MuonToSimCollection recSimColl_sta;
  MuonAssociatorByHits::SimToMuonCollection simRecColl_sta;
  edm::LogVerbatim("MuonMCClassifier") << "\n Global Loose Muon association by inner track ";  //
  assoByHits->associateMuons(recSimColl_sta, simRecColl_sta, selStaMuons, MuonAssociatorByHits::OuterTk, allTPs, &iEvent, &iSetup);  //

//  typedef MuonAssociatorByHits::MuonToSimCollection::const_iterator r2s_it;
//  typedef MuonAssociatorByHits::SimToMuonCollection::const_iterator s2r_it;
  size_t nmu = myMuons->size();
  if(debug) std::cout << "there are " << nmu << "  reco::Muons" << std::endl;
  if(debug) std::cout << "=================================" << std::endl;
  edm::LogVerbatim("MuonMCClassifier") << "\n There are " << nmu << " reco::Muons.";
/*
  //clear and initialise the std::vectors in myClassification()
  std::vector<int> classif, ext;
  std::vector<int> hitsPdgId, momPdgId, gmomPdgId, momStatus;
  std::vector<int> flav, momFlav, gmomFlav, hmomFlav;
  std::vector<int> tpId;
  std::vector<float> prodRho, prodZ, momRho, momZ;
  std::vector<float> tpAssoQuality;
*/
  std::vector<int> classift(nmu, 0), extt(nmu, 0);
  std::vector<int> hitsPdgIdt(nmu, 0), momPdgIdt(nmu, 0), gmomPdgIdt(nmu, 0), momStatust(nmu, 0);
  std::vector<int> flavt(nmu, 0), momFlavt(nmu, 0), gmomFlavt(nmu, 0), hmomFlavt(nmu, 0);
  std::vector<int> tpIdt(nmu, -1);
  std::vector<float> prodRhot(nmu, 0.0), prodZt(nmu, 0.0), momRhot(nmu, 0.0), momZt(nmu, 0.0);
  std::vector<float> tpAssoQualityt(nmu, -1);


  std::auto_ptr<reco::GenParticleCollection> secondariesforGlbT;     // output collection of secondary muons
  std::auto_ptr<reco::GenParticleCollection> secondariesforGlbL;     // output collection of secondary muons
  std::auto_ptr<reco::GenParticleCollection> secondariesforSta;     // output collection of secondary muons

  std::vector<int> muToPrimaryforGlbT(nmu, -1), muToSecondaryforGlbT(nmu, -1); // map from input into (index) in output, -1 for null
  std::vector<int> muToPrimaryforGlbL(nmu, -1), muToSecondaryforGlbL(nmu, -1); // map from input into (index) in output, -1 for null
  std::vector<int> muToPrimaryforSta(nmu, -1), muToSecondaryforSta(nmu, -1); // map from input into (index) in output, -1 for null

  if (linkToGenParticles_)
  {
    secondariesforGlbT.reset(new reco::GenParticleCollection());
    secondariesforGlbL.reset(new reco::GenParticleCollection());
    secondariesforSta.reset(new reco::GenParticleCollection());
//    edm::RefProd < reco::GenParticleCollection > priRP(edm::Handle < reco::GenParticleCollection > genParticles);
//    edm::RefProd < reco::GenParticleCollection > secRP(secHandle);
//    std::auto_ptr<edm::Association<reco::GenParticleCollection> > outPri(new edm::Association<reco::GenParticleCollection>(priRP)); //
//    std::auto_ptr<edm::Association<reco::GenParticleCollection> > outSec(new edm::Association<reco::GenParticleCollection>(secRP)); //
//    edm::Association<reco::GenParticleCollection>::Filler fillPri(*outPri), fillSec(*outSec);

  }

  //for GLOBAL TIGHT MUONS
  myClassification(nmu, myMuons, selTightMuons, recSimColl_glbTight, simRecColl_glbTight, classift, extt, hitsPdgIdt, momPdgIdt,
      gmomPdgIdt, momStatust, flavt, momFlavt, gmomFlavt, hmomFlavt, tpIdt, prodRhot, prodZt, momRhot, momZt, tpAssoQualityt,
      secondariesforGlbT, muToPrimaryforGlbT, muToSecondaryforGlbT, genParticlesforGlbT);

  writeValueMap(iEvent, myMuons, classift, "glbT");
  writeValueMap(iEvent, myMuons, extt, "extglbT");
  writeValueMap(iEvent, myMuons, flavt, "flavglbT");
  writeValueMap(iEvent, myMuons, tpIdt, "tpIdglbT");
  writeValueMap(iEvent, myMuons, hitsPdgIdt, "hitsPdgIdglbT");
  writeValueMap(iEvent, myMuons, momPdgIdt, "momPdgIdglbT");
  writeValueMap(iEvent, myMuons, momStatust, "momStatusglbT");
  writeValueMap(iEvent, myMuons, momFlavt, "momFlavglbT");
  writeValueMap(iEvent, myMuons, gmomPdgIdt, "gmomPdgIdglbT");
  writeValueMap(iEvent, myMuons, gmomFlavt, "gmomFlavglbT");
  writeValueMap(iEvent, myMuons, hmomFlavt, "hmomFlavglbT");
  writeValueMap(iEvent, myMuons, prodRhot, "prodRhoglbT");
  writeValueMap(iEvent, myMuons, prodZt, "prodZglbT");
  writeValueMap(iEvent, myMuons, momRhot, "momRhoglbT");
  writeValueMap(iEvent, myMuons, momZt, "momZglbT");
  writeValueMap(iEvent, myMuons, tpAssoQualityt, "tpAssoQualityglbT");

  if (linkToGenParticles_)
  {
//   edm::OrphanHandle < reco::GenParticleCollection > secHandle = iEvent.put(secondariesforGlbT, "secondariesforGlbT");
   edm::OrphanHandle < reco::GenParticleCollection > secHandleforGlbT = iEvent.put(secondariesforGlbT, "secondariesforGlbT");

//    edm::RefProd < reco::GenParticleCollection > priRP(genParticles);
//    edm::RefProd < reco::GenParticleCollection > secRP(secHandle);
//    std::auto_ptr<edm::Association<reco::GenParticleCollection> > outPri(new edm::Association<reco::GenParticleCollection>(priRP)); //
//    std::auto_ptr<edm::Association<reco::GenParticleCollection> > outSec(new edm::Association<reco::GenParticleCollection>(secRP)); //
//    edm::Association<reco::GenParticleCollection>::Filler fillPri(*outPri), fillSec(*outSec);

    edm::RefProd < reco::GenParticleCollection > priRPforGlbT(genParticlesforGlbT);
    edm::RefProd < reco::GenParticleCollection > secRPforGlbT(secHandleforGlbT);
    std::auto_ptr<edm::Association<reco::GenParticleCollection> > outPriforGlbT(new edm::Association<reco::GenParticleCollection>(priRPforGlbT)); //
    std::auto_ptr<edm::Association<reco::GenParticleCollection> > outSecforGlbT(new edm::Association<reco::GenParticleCollection>(secRPforGlbT)); //
    edm::Association<reco::GenParticleCollection>::Filler fillPriforGlbT(*outPriforGlbT), fillSecforGlbT(*outSecforGlbT);

    fillPriforGlbT.insert(myMuons, muToPrimaryforGlbT.begin(), muToPrimaryforGlbT.end());
    fillSecforGlbT.insert(myMuons, muToSecondaryforGlbT.begin(), muToSecondaryforGlbT.end());
    fillPriforGlbT.fill();
    fillSecforGlbT.fill();
//    iEvent.put(outPriforGlbT, "toPrimariesGlbT");
//    iEvent.put(outSecforGlbT, "toSecondariesGlbT");
    iEvent.put(outPriforGlbT, "toPrimaries");
    iEvent.put(outSecforGlbT, "toSecondaries");
  }


  //for GLOBAL SOFT MUONS
  std::vector<int> classifl(nmu, 0), extl(nmu, 0);
  std::vector<int> hitsPdgIdl(nmu, 0), momPdgIdl(nmu, 0), gmomPdgIdl(nmu, 0), momStatusl(nmu, 0);
  std::vector<int> flavl(nmu, 0), momFlavl(nmu, 0), gmomFlavl(nmu, 0), hmomFlavl(nmu, 0);
  std::vector<int> tpIdl(nmu, -1);
  std::vector<float> prodRhol(nmu, 0.0), prodZl(nmu, 0.0), momRhol(nmu, 0.0), momZl(nmu, 0.0);
  std::vector<float> tpAssoQualityl(nmu, -1);

  myClassification(nmu, myMuons, selSoftMuons, recSimColl_glbloose, simRecColl_glbloose, classifl, extl, hitsPdgIdl, momPdgIdl,
      gmomPdgIdl, momStatusl, flavl, momFlavl, gmomFlavl, hmomFlavl, tpIdl, prodRhol, prodZl, momRhol, momZl, tpAssoQualityl,
      secondariesforGlbL, muToPrimaryforGlbL, muToSecondaryforGlbL, genParticlesforGlbL);

  writeValueMap(iEvent, myMuons, classifl, "glbL");
  writeValueMap(iEvent, myMuons, extl, "extglbL");
  writeValueMap(iEvent, myMuons, flavl, "flavglbL");
  writeValueMap(iEvent, myMuons, tpIdl, "tpIdglbL");
  writeValueMap(iEvent, myMuons, hitsPdgIdl, "hitsPdgIdglbL");
  writeValueMap(iEvent, myMuons, momPdgIdl, "momPdgIdglbL");
  writeValueMap(iEvent, myMuons, momStatusl, "momStatusglbL");
  writeValueMap(iEvent, myMuons, momFlavl, "momFlavglbL");
  writeValueMap(iEvent, myMuons, gmomPdgIdl, "gmomPdgIdglbL");
  writeValueMap(iEvent, myMuons, gmomFlavl, "gmomFlavglbL");
  writeValueMap(iEvent, myMuons, hmomFlavl, "hmomFlavglbL");
  writeValueMap(iEvent, myMuons, prodRhol, "prodRhoglbL");
  writeValueMap(iEvent, myMuons, prodZl, "prodZglbL");
  writeValueMap(iEvent, myMuons, momRhol, "momRhoglbL");
  writeValueMap(iEvent, myMuons, momZl, "momZglbL");
  writeValueMap(iEvent, myMuons, tpAssoQualityl, "tpAssoQualityglbL");

  if (linkToGenParticles_)
  {
    edm::OrphanHandle < reco::GenParticleCollection > secHandleforGlbL = iEvent.put(secondariesforGlbL, "secondariesforGlbL");

    edm::RefProd < reco::GenParticleCollection > priRPforGlbL(genParticlesforGlbL);
    edm::RefProd < reco::GenParticleCollection > secRPforGlbL(secHandleforGlbL);
    std::auto_ptr<edm::Association<reco::GenParticleCollection> > outPriforGlbL(new edm::Association<reco::GenParticleCollection>(priRPforGlbL)); //
    std::auto_ptr<edm::Association<reco::GenParticleCollection> > outSecforGlbL(new edm::Association<reco::GenParticleCollection>(secRPforGlbL)); //
    edm::Association<reco::GenParticleCollection>::Filler fillPriforGlbL(*outPriforGlbL), fillSecforGlbL(*outSecforGlbL);

    fillPriforGlbL.insert(myMuons, muToPrimaryforGlbL.begin(), muToPrimaryforGlbL.end());
    fillSecforGlbL.insert(myMuons, muToSecondaryforGlbL.begin(), muToSecondaryforGlbL.end());
    fillPriforGlbL.fill();
    fillSecforGlbL.fill();
//    iEvent.put(outPriforGlbL, "toPrimariesGlbL");
//    iEvent.put(outSecforGlbL, "toSecondariesGlbL");
    iEvent.put(outPriforGlbL, "toPrimaries");
    iEvent.put(outSecforGlbL, "toSecondaries");

  }


  //for STANDALONE MUONS
  std::vector<int> classifs(nmu, 0), exts(nmu, 0);
  std::vector<int> hitsPdgIds(nmu, 0), momPdgIds(nmu, 0), gmomPdgIds(nmu, 0), momStatuss(nmu, 0);
  std::vector<int> flavs(nmu, 0), momFlavs(nmu, 0), gmomFlavs(nmu, 0), hmomFlavs(nmu, 0);
  std::vector<int> tpIds(nmu, -1);
  std::vector<float> prodRhos(nmu, 0.0), prodZs(nmu, 0.0), momRhos(nmu, 0.0), momZs(nmu, 0.0);
  std::vector<float> tpAssoQualitys(nmu, -1);

  myClassification(nmu, myMuons, selStaMuons, recSimColl_sta, simRecColl_sta, classifs, exts, hitsPdgIds, momPdgIds,
      gmomPdgIds, momStatuss, flavs, momFlavs, gmomFlavs, hmomFlavs, tpIds, prodRhos, prodZs, momRhos, momZs, tpAssoQualitys,
      secondariesforSta, muToPrimaryforSta, muToSecondaryforSta, genParticlesforSta);

  writeValueMap(iEvent, myMuons, classifs, "sta");
  writeValueMap(iEvent, myMuons, exts, "extsta");
  writeValueMap(iEvent, myMuons, flavs, "flavsta");
  writeValueMap(iEvent, myMuons, tpIds, "tpIdsta");
  writeValueMap(iEvent, myMuons, hitsPdgIds, "hitsPdgIdsta");
  writeValueMap(iEvent, myMuons, momPdgIds, "momPdgIdsta");
  writeValueMap(iEvent, myMuons, momStatuss, "momStatussta");
  writeValueMap(iEvent, myMuons, momFlavs, "momFlavsta");
  writeValueMap(iEvent, myMuons, gmomPdgIds, "gmomPdgIdsta");
  writeValueMap(iEvent, myMuons, gmomFlavs, "gmomFlavsta");
  writeValueMap(iEvent, myMuons, hmomFlavs, "hmomFlavsta");
  writeValueMap(iEvent, myMuons, prodRhos, "prodRhosta");
  writeValueMap(iEvent, myMuons, prodZs, "prodZsta");
  writeValueMap(iEvent, myMuons, momRhos, "momRhosta");
  writeValueMap(iEvent, myMuons, momZs, "momZsta");
  writeValueMap(iEvent, myMuons, tpAssoQualitys, "tpAssoQualitysta");

  if (linkToGenParticles_)
  {
    edm::OrphanHandle < reco::GenParticleCollection > secHandleforSta = iEvent.put(secondariesforSta, "secondariesforSta");

    edm::RefProd < reco::GenParticleCollection > priRPforSta(genParticlesforSta);
    edm::RefProd < reco::GenParticleCollection > secRPforSta(secHandleforSta);
    std::auto_ptr<edm::Association<reco::GenParticleCollection> > outPriforSta(new edm::Association<reco::GenParticleCollection>(priRPforSta)); //
    std::auto_ptr<edm::Association<reco::GenParticleCollection> > outSecforSta(new edm::Association<reco::GenParticleCollection>(secRPforSta)); //
    edm::Association<reco::GenParticleCollection>::Filler fillPriforSta(*outPriforSta), fillSecforSta(*outSecforSta);

    fillPriforSta.insert(myMuons, muToPrimaryforSta.begin(), muToPrimaryforSta.end());
    fillSecforSta.insert(myMuons, muToSecondaryforSta.begin(), muToSecondaryforSta.end());
    fillPriforSta.fill();
    fillSecforSta.fill();
//    iEvent.put(outPriforSta, "toPrimariesSta");
//    iEvent.put(outSecforSta, "toSecondariesSta");
    iEvent.put(outPriforSta, "toPrimaries");
    iEvent.put(outSecforSta, "toSecondaries");



  }
}//end producer

void MuonMCClassifier::myClassification(size_t myNmu, edm::Handle<edm::View<reco::Muon> > myMuons_handle,
    edm::RefToBaseVector<reco::Muon> selMuons, MuonAssociatorByHits::MuonToSimCollection myRecSimColl,
    MuonAssociatorByHits::SimToMuonCollection mySimRecColl, std::vector<int> &myClassif, std::vector<int> &myExt,
    std::vector<int> &myHitsPdgId, std::vector<int> &myMomPdgId, std::vector<int> &myGmomPdgId,
    std::vector<int> &myMomStatus, std::vector<int> &myFlav, std::vector<int> &myMomFlav, std::vector<int> &myGmomFlav,
    std::vector<int> &myHmomFlav, std::vector<int> &myTpId, std::vector<float> &myProdRho, std::vector<float> &myProdZ,
    std::vector<float> &myMomRho, std::vector<float> &myMomZ, std::vector<float> &myTpAssoQuality,
    std::auto_ptr<reco::GenParticleCollection> &secondaries, std::vector<int> &muToPrimary,
    std::vector<int> &muToSecondary, edm::Handle < reco::GenParticleCollection > &genParticles
    )
{

/*
  //initialise the std::vectors

  myExt.clear();
  myClassif.clear();
  myHitsPdgId.clear();
  myMomPdgId.clear();
  myGmomPdgId.clear();
  myMomStatus.clear();
  myFlav.clear();
  myMomFlav.clear();
  myGmomFlav.clear();
  myHmomFlav.clear();
  myTpId.clear();
  myProdRho.clear();
  myProdZ.clear();
  myMomRho.clear();
  myMomZ.clear();
  myTpAssoQuality.clear();
  muToPrimary.clear();
  muToSecondary.clear();
  for (unsigned int k = 0; k < myNmu; k++)
  {

    myClassif[k] = 0;
    myExt[k] = 0;
    myHitsPdgId[k] = 0;
    myMomPdgId[k] = 0;
    myGmomPdgId[k] = 0;
    myMomStatus[k] = 0;
    myMomStatus[k] = 0;
    myMomFlav[k] = 0;
    myGmomFlav[k] = 0;
    myHmomFlav[k] = 0;
    myTpId[k] = -1;
    myProdRho[k] = 0.;
    myProdZ[k] = 0.;
    myMomRho[k] = 0.;
    myMomZ[k] = 0.;
    myTpAssoQuality[k] = -1.;
    muToPrimary[k] = -1;
    muToSecondary[k] = -1;
  }
*/
  std::map<TrackingParticleRef, int> tpToSecondaries;     // map from tp to (index+1) in output collection

  typedef MuonAssociatorByHits::MuonToSimCollection::const_iterator r2s_it;
  typedef MuonAssociatorByHits::SimToMuonCollection::const_iterator s2r_it;

  for (size_t i = 0; i < myNmu; ++i)     //loop over reco::muons
  {
    edm::LogVerbatim("MuonMCClassifier") << "\n reco::Muons # " << i;
    edm::RefToBase < reco::Muon > mu = myMuons_handle->refAt(i);
    //   if (hasMuonCut_ && (std::find(selMuons.begin(), selMuons.end(), mu) == selMuons.end())) //if no match with the selected muons and cuts - fail and go further
    if (std::find(selMuons.begin(), selMuons.end(), mu) == selMuons.end()) //if no match with the selected muons and cuts - fail and go further
    {
      edm::LogVerbatim("MuonMCClassifier") << "\t muon didn't pass the selection. classified as -99 and skipped";
      if (debug)
        std::cout << "\t echo from myClassification" << std::endl;
      myClassif[i] = -99;
      continue;
    }

    TrackingParticleRef tp;
    edm::RefToBase < reco::Muon > muMatchBack;
    r2s_it match = myRecSimColl.find(mu);   //iterator reco to sim collection,
    s2r_it matchback;                     //iterator sim to reco collection
    if (match != myRecSimColl.end())        //tyrsi syvpadenie na reco 4asticata w sim kolekciata
    {
      edm::LogVerbatim("MuonMCClassifier") << "\t Reco to Sim matched Ok...";
      if (debug)
        std::cout << "\t Reco to Sim matched Ok..." << std::endl;
      // match->second is vector, front is first element, first is the ref (second would be the quality)
      tp = match->second.front().first;   //t.e. vizma referencia kym mu
      myTpId[i] = tp.isNonnull() ? tp.key() : -1; // we check, even if null refs should not appear here at all; ako ref e nenuleva wzima Id na 4asticata
      myTpAssoQuality[i] = match->second.front().second;    //vzima quality
//            s2r_it matchback = mySimRecColl.find(tp);
      matchback = mySimRecColl.find(tp);    //sim to reco match; t.e. obratno - tyrsi machnata tp v reco colekciata
      if (debug)
        std::cout << "tp Id\t" << myTpId[i] << std::endl;
      if (matchback != mySimRecColl.end())  //ako ima syvpadenie
      {
        muMatchBack = matchback->second.front().first;    //ako ima syvpadenie vzima referencia kym reco 4asticata
      }
      else
      {
        edm::LogWarning("MuonMCClassifier") << "\n***WARNING:  This I do NOT understand: why no match back? *** \n";
        if (debug)
          std::cout << "\n***WARNING:  This I do NOT understand: why no match back? *** \n" << std::endl;
      }
    }    //end reco to sim match OK

    if (tp.isNonnull())
    {
      bool isGhost = muMatchBack != mu;
      if (isGhost)
        edm::LogVerbatim("MuonMCClassifier") << "\t This seems a GHOST ! myClassif[i] will be < 0";

      myHitsPdgId[i] = tp->pdgId();
      myProdRho[i] = tp->vertex().Rho();
      myProdZ[i] = tp->vertex().Z();
      edm::LogVerbatim("MuonMCClassifier") << "\t TP pdgId = " << myHitsPdgId[i] << ", vertex rho = " << myProdRho[i]
          << ", z = " << myProdZ[i];

      // Try to extract mother and grand mother of this muon.
      // Unfortunately, SIM and GEN histories require diffent code :-(
      if (!tp->genParticles().empty())
      { // Muon is in GEN
        reco::GenParticleRef genp = tp->genParticles()[0];
        reco::GenParticleRef genMom = genp->numberOfMothers() > 0 ? genp->motherRef() : reco::GenParticleRef();
        if (genMom.isNonnull())
        {
          myMomPdgId[i] = genMom->pdgId();
          myMomStatus[i] = genMom->status();
          myMomRho[i] = genMom->vertex().Rho();
          myMomZ[i] = genMom->vz();
          edm::LogVerbatim("MuonMCClassifier") << "\t Particle pdgId = " << myHitsPdgId[i] << " produced at rho = "
              << myProdRho[i] << ", z = " << myProdZ[i] << ", has GEN mother pdgId = " << myMomPdgId[i];
          reco::GenParticleRef genGMom = genMom->numberOfMothers() > 0 ? genMom->motherRef() : reco::GenParticleRef();
          if (genGMom.isNonnull())
          {
            myGmomPdgId[i] = genGMom->pdgId();
            edm::LogVerbatim("MuonMCClassifier") << "\t\t mother prod. vertex rho = " << myMomRho[i] << ", z = "
                << myMomZ[i] << ", grand-mom pdgId = " << myGmomPdgId[i];
          }
          // in this case, we might want to know the heaviest mom flavour
          for (reco::GenParticleRef nMom = genMom; nMom.isNonnull() && abs(nMom->pdgId()) >= 100; // stop when we're no longer looking at hadrons or mesons
              nMom = nMom->numberOfMothers() > 0 ? nMom->motherRef() : reco::GenParticleRef())
          {
            int myFlav = flavour(nMom->pdgId());
            if (myHmomFlav[i] < myFlav)
              myHmomFlav[i] = myFlav;
            edm::LogVerbatim("MuonMCClassifier") << "\t\t backtracking flavour: mom pdgId = " << nMom->pdgId()
                << ", flavour = " << myFlav << ", heaviest so far = " << myHmomFlav[i];

std::cout << "from Gen: backtracking flavour: mom pdgId = " << nMom->pdgId() << ", flavour = " << myFlav << ", heaviest so far = " << myHmomFlav[i] << std::endl;
          }
        }
      }

      else
      { // Muon is in SIM Only
        TrackingParticleRef simMom = getTpMother(tp);
        if (simMom.isNonnull())
        {
          myMomPdgId[i] = simMom->pdgId();
          myMomRho[i] = simMom->vertex().Rho();
          myMomZ[i] = simMom->vertex().Z();
          edm::LogVerbatim("MuonMCClassifier") << "\t Particle pdgId = " << myHitsPdgId[i] << " produced at rho = "
              << myProdRho[i] << ", z = " << myProdZ[i] << ", has SIM mother pdgId = " << myMomPdgId[i]
              << " produced at rho = " << simMom->vertex().Rho() << ", z = " << simMom->vertex().Z();
          if (!simMom->genParticles().empty())
          {
            myMomStatus[i] = simMom->genParticles()[0]->status();
            reco::GenParticleRef genGMom = (
                simMom->genParticles()[0]->numberOfMothers() > 0 ? simMom->genParticles()[0]->motherRef() :
                    reco::GenParticleRef());
            if (genGMom.isNonnull())
              myGmomPdgId[i] = genGMom->pdgId();
            edm::LogVerbatim("MuonMCClassifier") << "\t\t SIM mother is in GEN (status " << myMomStatus[i]
                << "), grand-mom id = " << myGmomPdgId[i];
          }
          else
          {
            myMomStatus[i] = -1;
            TrackingParticleRef simGMom = getTpMother(simMom);
            if (simGMom.isNonnull())
              myGmomPdgId[i] = simGMom->pdgId();
            edm::LogVerbatim("MuonMCClassifier") << "\t\t SIM mother is in SIM only, grand-mom id = " << myGmomPdgId[i];
          }
        }
        else
        {
          edm::LogVerbatim("MuonMCClassifier") << "\t Particle pdgId = " << myHitsPdgId[i] << " produced at rho = "
              << myProdRho[i] << ", z = " << myProdZ[i] << ", has NO mother!";
        }
      }
/////////////////////////////////////////////////////end old version of  Muon is in SIM

      myMomFlav[i] = flavour(myMomPdgId[i]);
      myGmomFlav[i] = flavour(myGmomPdgId[i]);

      // Check first IF this is a muon at all
      if (abs(tp->pdgId()) != 13)
      {
        myClassif[i] = isGhost ? -1 : 1;
        myExt[i] = isGhost ? -1 : 1;
        edm::LogVerbatim("MuonMCClassifier") << "\t This is not a muon. Sorry. myClassif[i] = " << myClassif[i];
        if (debug)
          std::cout << "\t This is not a muon. Sorry ! myClassif[i] = " << myClassif[i] << std::endl;
        continue;
      }

      // Is this SIM muon also a GEN muon, with a mother?
      if (!tp->genParticles().empty() && (myMomPdgId[i] != 0))
      {
        if (abs(myMomPdgId[i]) < 100 && (abs(myMomPdgId[i]) != 15))
        {
          myClassif[i] = isGhost ? -4 : 4;
          myFlav[i] = (abs(myMomPdgId[i]) == 15 ? 15 : 13);
          edm::LogVerbatim("MuonMCClassifier") << "\t This seems PRIMARY MUON ! myClassif[i] = " << myClassif[i];
          if (debug)
            std::cout << "\t This seems PRIMARY MUON ! myClassif[i] = " << myClassif[i] << std::endl;
          myExt[i] = 10;
        }
        else if (myMomFlav[i] == 4 || myMomFlav[i] == 5 || myMomFlav[i] == 15)
        {
          myClassif[i] = isGhost ? -3 : 3;
          myFlav[i] = myMomFlav[i];
          if (myMomFlav[i] == 15)
            myExt[i] = 9; // tau->mu
          else if (myMomFlav[i] == 5)
            myExt[i] = 8; // b->mu
          else if (myHmomFlav[i] == 5)
            myExt[i] = 7; // b->c->mu
          else
            myExt[i] = 6; // c->mu
          edm::LogVerbatim("MuonMCClassifier") << "\t This seems HEAVY FLAVOUR ! myClassif[i] = " << myClassif[i];
          if (debug)
            std::cout << "\t This seems HEAVY FLAVOUR ! myClassif[i] = " << myClassif[i] << std::endl;
        }
        else
        {
          myClassif[i] = isGhost ? -2 : 2;
          myFlav[i] = myMomFlav[i];
          edm::LogVerbatim("MuonMCClassifier") << "\t This seems LIGHT FLAVOUR ! myClassif[i] = " << myClassif[i];
          if (debug)
            std::cout << "\t This seems LIGHT FLAVOUR ! myClassif[i] = " << myClassif[i] << std::endl;
        }
      }
      else
      {
        myClassif[i] = isGhost ? -2 : 2;
        myFlav[i] = myMomFlav[i];
        edm::LogVerbatim("MuonMCClassifier") << "\t This seems LIGHT FLAVOUR ! myClassif[i] = " << myClassif[i];
        if (debug)
          std::cout << "\t This seems LIGHT FLAVOUR ! myClassif[i] = " << myClassif[i] << std::endl;
      }
      // extended classification
      if (myMomPdgId[i] == 0)
        myExt[i] = 2; // if it has no mom, it's not a primary particle so it won't be in ppMuX
      else if (abs(myMomPdgId[i]) < 100)
        myExt[i] = (myMomFlav[i] == 15 ? 9 : 10); // primary mu, or tau->mu
      else if (myMomFlav[i] == 5)
        myExt[i] = 8; // b->mu
      else if (myMomFlav[i] == 4)
        myExt[i] = (myHmomFlav[i] == 5 ? 7 : 6); // b->c->mu and c->mu
      else if (myMomStatus[i] != -1)
      { // primary light particle
        int id = abs(myMomPdgId[i]);
        if (id != /*pi+*/211 && id != /*K+*/321 && id != 130 /*K0L*/)
          myExt[i] = 5; // other light particle, possibly short-lived
        else if (myProdRho[i] < decayRho_ && abs(myProdZ[i]) < decayAbsZ_)
          myExt[i] = 4; // decay a la ppMuX (primary pi/K within a cylinder)
        else
          myExt[i] = 3; // late decay that wouldn't be in ppMuX
      }
      else
        myExt[i] = 2; // decay of non-primary particle, would not be in ppMuX
      if (isGhost)
        myExt[i] = -myExt[i];

      if (linkToGenParticles_ && abs(myExt[i]) >= 2)
      {
        // Link to the genParticle if possible, but not decays in flight (in ppMuX they're in GEN block, but they have wrong parameters)
        if (!tp->genParticles().empty() && abs(myExt[i]) >= 5)
        {
          if (genParticles.id() != tp->genParticles().id())
          {
            throw cms::Exception("Configuration") << "Product ID mismatch between the genParticle collection ("
                << genParticles_ << ", id " << genParticles.id() << ") and the references in the TrackingParticles (id "
                << tp->genParticles().id() << ")\n";
          }
          muToPrimary[i] = tp->genParticles()[0].key();
        }
        else
        {
          // Don't put the same trackingParticle twice!
          int &indexPlus1 = tpToSecondaries[tp]; // will create a 0 if the tp is not in the list already
          if (indexPlus1 == 0)
            indexPlus1 = convertAndPush(*tp, *secondaries, getTpMother(tp), genParticles) + 1;
          muToSecondary[i] = indexPlus1 - 1;
        }
      }
      edm::LogVerbatim("MuonMCClassifier") << "\t Extended clasmyClassifion code = " << myExt[i];
//std::cout << "Extended clasmyClassifion code = " << myExt[i] << std::endl;
    } //ent tp is nonNull
  } //ent loop over myMuons_handle size
}


template<typename T>
void MuonMCClassifier::writeValueMap(edm::Event &iEvent, const edm::Handle<edm::View<reco::Muon> > & handle,
    const std::vector<T> & values, const std::string & label) const
{
  using namespace edm;
  using namespace std;
  auto_ptr<ValueMap<T> > valMap(new ValueMap<T>());
  typename edm::ValueMap<T>::Filler filler(*valMap);
  filler.insert(handle, values.begin(), values.end());
  filler.fill();
  iEvent.put(valMap, label);
}

int MuonMCClassifier::flavour(int pdgId) const
{
  int flav = abs(pdgId);
  // for quarks, leptons and bosons except gluons, take their pdgId
  // muons and taus have themselves as flavour
  if (flav <= 37 && flav != 21)
    return flav;
  // look for barions
  int bflav = ((flav / 1000) % 10);
  if (bflav != 0)
    return bflav;
  // look for mesons
  int mflav = ((flav / 100) % 10);
  if (mflav != 0)
    return mflav;
  return 0;
}

// push secondary in collection.
// if it has a primary mother link to it.
int MuonMCClassifier::convertAndPush(const TrackingParticle &tp, reco::GenParticleCollection &out,
    const TrackingParticleRef & simMom, const edm::Handle<reco::GenParticleCollection> & genParticles) const
{
  out.push_back(reco::GenParticle(tp.charge(), tp.p4(), tp.vertex(), tp.pdgId(), tp.status(), true));
  if (simMom.isNonnull() && !simMom->genParticles().empty())
  {
    if (genParticles.id() != simMom->genParticles().id())
    {
      throw cms::Exception("Configuration") << "Product ID mismatch between the genParticle collection ("
          << genParticles_ << ", id " << genParticles.id() << ") and the references in the TrackingParticles (id "
          << simMom->genParticles().id() << ")\n";
    }
    out.back().addMother(simMom->genParticles()[0]);
  }
  return out.size() - 1;
}

//define this as a plug-in
DEFINE_FWK_MODULE (MuonMCClassifier);

