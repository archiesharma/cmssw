// -*- C++ -*-
//
// Package: MuonMCClassifAndAna
// Class: MuonMCClassifAndAna
//
/**\class MuonMCClassifAndAna MuonMCClassifAndAna.cc MuonAnalysis/MuonAssociators/src/MuonMCClassifAndAna.cc
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
// Original Author: Nov 16 16:12 (lxplus231.cern.ch)
// Created: Sun Nov 16 16:14:09 CET 2008
// $Id: MuonMCClassifAndAna.cc,v 1.9 2013/06/24 12:53:19 speer Exp $
//
//
// system include files
#include <memory>
#include <set>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
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
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//#include <cmath>
#include "TH2F.h"
#include "TH1F.h"
//
// class decleration
class MuonMCClassifAndAna: public edm::EDAnalyzer
{
public:
  explicit MuonMCClassifAndAna(const edm::ParameterSet&);
  bool isTight(const edm::Event& iEvent, reco::Muon myMuon);
  bool isSoft(reco::Muon myMuon, reco::VertexCollection::const_iterator myPVit,
      const reco::VertexCollection myVertices);
  ~MuonMCClassifAndAna();
  void myClassification(size_t myNmu, edm::Handle<edm::View<reco::Muon> > myMuons_handle,
      edm::RefToBaseVector<reco::Muon> selMuons, MuonAssociatorByHits::MuonToSimCollection myRecSimColl,
      MuonAssociatorByHits::SimToMuonCollection mySimRecColl,
      std::vector<int> &myClassif,
      std::vector<int> &myExt, std::vector<int> &myHitsPdgId, std::vector<int> &myMomPdgId,
      std::vector<int> &myGmomPdgId, std::vector<int> &myMomStatus, std::vector<int> &myFlav,
      std::vector<int> &myMomFlav, std::vector<int> &myGmomFlav, std::vector<int> &myHmomFlav, std::vector<int> &myTpId,
      std::vector<float> &myProdRho, std::vector<float> &myProdZ, std::vector<float> &myMomRho,
      std::vector<float> &myMomZ, std::vector<float> &myTpAssoQuality,
      std::auto_ptr<reco::GenParticleCollection> &secondaries, std::vector<int> &muToPrimary,
      std::vector<int> &muToSecondary, edm::Handle<reco::GenParticleCollection> &genParticles, int &binp, int &binl,
      int &binh, int &binpr, int &binn, int &bingst, TH1F *histoPt_p, TH1F *histoPt_l, TH1F *histoPt_h,
      TH1F *histoPt_pr, TH1F *histoPt_n, TH1F *histoPt_gst, TH1F *histoEta_p, TH1F *histoEta_l, TH1F *histoEta_h,
      TH1F *histoEta_pr, TH1F *histoEta_n, TH1F *histoEta_gst, TH1F *histoPhi_p, TH1F *histoPhi_l, TH1F *histoPhi_h,
      TH1F *histoPhi_pr, TH1F *histoPhi_n, TH1F *histoPhi_gst, TH1F *histomatchedstation_p, TH1F *histomatchedstation_l,
      TH1F *histomatchedstation_h, TH1F *histomatchedstation_pr, TH1F *histomatchedstation_n, TH1F *histomatchedstation_gst,
      TH1F *histoMuonsOOT_p, TH1F *histoMuonsOOT_l, TH1F *histoMuonsOOT_h, TH1F *histoMuonsOOT_pr, TH1F *histoMuonsOOT_n,
      TH1F *histoMuonsOOT_gst, TH1F *histoMuonsPtOOT_p, TH1F *histoMuonsPtOOT_l, TH1F *histoMuonsPtOOT_h, TH1F *histoMuonsPtOOT_pr,
      TH1F *histoMuonsPtOOT_n, TH1F *histoMuonsPtOOT_gst, TH1F *histoMuonsEtaOOT_p, TH1F *histoMuonsEtaOOT_l, TH1F *histoMuonsEtaOOT_h, 
      TH1F *histoMuonsEtaOOT_pr, TH1F *histoMuonsEtaOOT_n, TH1F *histoMuonsEtaOOT_gst, TH1F *histoMuonsValidHitsOOT_p,
      TH1F *histoMuonsValidHitsOOT_l, TH1F *histoMuonsValidHitsOOT_h, TH1F *histoMuonsValidHitsOOT_pr, TH1F *histoMuonsValidHitsOOT_n,
      TH1F *histoMuonsValidHitsOOT_gst, TH1F *histoMuonsValidGemHitsOOT_p, TH1F *histoMuonsValidGemHitsOOT_l,
      TH1F *histoMuonsValidGemHitsOOT_h, TH1F *histoMuonsValidGemHitsOOT_pr, TH1F *histoMuonsValidGemHitsOOT_n,
      TH1F *histoMuonsValidGemHitsOOT_gst,TH1F *histoMuonsInTime_p, TH1F *histoMuonsInTime_l, TH1F *histoMuonsInTime_h,
      TH1F *histoMuonsInTime_pr, TH1F *histoMuonsInTime_n, TH1F *histoMuonsInTime_gst, TH1F *histoMuonsPtInTime_p,
      TH1F *histoMuonsPtInTime_l, TH1F *histoMuonsPtInTime_h, TH1F *histoMuonsPtInTime_pr, TH1F *histoMuonsPtInTime_n, 
      TH1F *histoMuonsPtInTime_gst, TH1F *histoMuonsEtaInTime_p, TH1F *histoMuonsEtaInTime_l, TH1F *histoMuonsEtaInTime_h, 
      TH1F *histoMuonsEtaInTime_pr, TH1F *histoMuonsEtaInTime_n, TH1F *histoMuonsEtaInTime_gst, TH1F *histoMuonsValidHitsInTime_p,
      TH1F *histoMuonsValidHitsInTime_l, TH1F *histoMuonsValidHitsInTime_h, TH1F *histoMuonsValidHitsInTime_pr,
      TH1F *histoMuonsValidHitsInTime_n, TH1F *histoMuonsValidHitsInTime_gst, TH1F *histoMuonsValidGemHitsInTime_p,
      TH1F *histoMuonsValidGemHitsInTime_l, TH1F *histoMuonsValidGemHitsInTime_h, TH1F *histoMuonsValidGemHitsInTime_pr,
      TH1F *histoMuonsValidGemHitsInTime_n, TH1F *histoMuonsValidGemHitsInTime_gst, int myCase);
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
// virtual void produce(edm::Event&, const edm::EventSetup&) override;
  bool debug;
//from ana
// int evtsMu = 0;
  int glb_pch;
  int glb_l;
  int glb_h;
  int glb_prm;
  int glb_n;
  int glb_gst;
  int trk_pch;
  int trk_l;
  int trk_h;
  int trk_prm;
  int trk_n;
  int trk_gst;
  int sta_pch;
  int sta_l;
  int sta_h;
  int sta_prm;
  int sta_n;
  int sta_gst;
  int good_pch;
  int good_l;
  int good_h;
  int good_prm;
  int good_n;
  int good_gst;
  int tight_pch;
  int tight_l;
  int tight_h;
  int tight_prm;
  int tight_n;
  int tight_gst;
  int soft_pch;
  int soft_l;
  int soft_h;
  int soft_prm;
  int soft_n;
  int soft_gst;
  int loose_pch;
  int loose_l;
  int loose_h;
  int loose_prm;
  int loose_n;
  int loose_gst;
//test reco
  TH1F *allrecoMuonsEta;
  TH1F *allrecoMuonsPt;
  TH1F *allrecoMuonsPhi;
  TH1F *histoValidMuHits;
  TH1F *histoValidMuGeMHits;
  TH1F *histoMuStaWithValidHits;
  TH1F *histoGemStaWithValidHits;
  TH1F *histoValidMuHitsOut;
  TH1F *histoValidMuGeMHitsOut;
  TH1F *histoMuStaWithValidHitsOut;
  TH1F *histoGemStaWithValidHitsOut;
  TH1F *histoValidMuHitsTight;
  TH1F *histoValidMuGeMHitsTight;
  TH1F *histoMuStaWithValidHitsTight;
  TH1F *histoGemStaWithValidHitsTight;
  TH1F *histoValidMuHitsOutTight;
  TH1F *histoValidMuGeMHitsOutTight;
  TH1F *histoMuStaWithValidHitsOutTight;
  TH1F *histoGemStaWithValidHitsOutTight;
  TH1F *histoValidMuHitsSoft;
  TH1F *histoValidMuGeMHitsSoft;
  TH1F *histoMuStaWithValidHitsSoft;
  TH1F *histoGemStaWithValidHitsSoft;
  TH1F *histoValidMuHitsOutSoft;
  TH1F *histoValidMuGeMHitsOutSoft;
  TH1F *histoMuStaWithValidHitsOutSoft;
  TH1F *histoGemStaWithValidHitsOutSoft;
  TH1F *histoValidMuHitsLoose;
  TH1F *histoValidMuGeMHitsLoose;
  TH1F *histoMuStaWithValidHitsLoose;
  TH1F *histoGemStaWithValidHitsLoose;
  TH1F *histoValidMuHitsOutLoose;
  TH1F *histoValidMuGeMHitsOutLoose;
  TH1F *histoMuStaWithValidHitsOutLoose;
  TH1F *histoGemStaWithValidHitsOutLoose;

//////////// plots for matched stations///////////////////////////////////
  TH1F *GlobalMuons_matchedstation;
  TH1F *StaMuons_matchedstation;  
  TH1F *TightMuons_matchedstation; 
  TH1F *LooseMuons_matchedstation;

/////// Plots for inTime and OOT/////////////////////////////////////////
  TH1F *StaMuonsOOT;
  TH1F *StaMuonsInTime; 
  TH1F *StaMuonsPt_OOT;
  TH1F *StaMuonsPt_InTime;
  TH1F *StaMuonsEta_OOT;
  TH1F *StaMuonsEta_InTime;
  TH1F *StaMuonsValidHits_OOT;
  TH1F *StaMuonsValidHits_InTime;
  TH1F *StaMuonsValidGemHits_OOT;
  TH1F *StaMuonsValidGemHits_InTime;

//test reco end
  TH1F *N_Muons;
  TH1F *N_Muons_norm;
  TH1F *totalevents_withMuons;
  TH1F *numb_MuPerEvent;
  TH1F *numb_PV;
  TH1F *numb_SimPV;
  TH1F *Pt_TrakingParticle;
  TH1F *PtMuons_TrakingParticle;
  TH1F *EtaMuons_TrakingParticle;
  TH1F *allMuonsPt;
  TH1F *allMuonsEta;
  TH1F *allMuonsPhi;
  TH1F *GlobalMuonsPt;
  TH1F *GlobalMuonsEta;
  TH1F *GlobalMuonsPhi;
  TH1F *TrackerMuonsPt;
  TH1F *TrackerMuonsEta;
  TH1F *TrackerMuonsPhi;
  TH1F *StaMuonsPt;
  TH1F *StaMuonsEta;
  TH1F *StaMuonsPhi;
  TH1F *GoodMuonsPt;
  TH1F *GoodMuonsEta;
  TH1F *GoodMuonsPhi;
  TH1F *TightMuonsPt;
  TH1F *TightMuonsEta;
  TH1F *TightMuonsPhi;
  TH1F *SoftMuonsPt;
  TH1F *SoftMuonsEta;
  TH1F *SoftMuonsPhi;
  TH1F *LooseMuonsPt;
  TH1F *LooseMuonsEta;
  TH1F *LooseMuonsPhi;
//pt
//global
  TH1F *muPt_Glb_punch;
  TH1F *muPt_Glb_lightflv;
  TH1F *muPt_Glb_heavyflv;
  TH1F *muPt_Glb_primary;
  TH1F *muPt_Glb_noise;
  TH1F *muPt_Glb_ghost;
//trk
  TH1F *muPt_Trk_punch;
  TH1F *muPt_Trk_lightflv;
  TH1F *muPt_Trk_heavyflv;
  TH1F *muPt_Trk_primary;
  TH1F *muPt_Trk_noise;
  TH1F *muPt_Trk_ghost;
//sta
  TH1F *muPt_Sta_punch;
  TH1F *muPt_Sta_lightflv;
  TH1F *muPt_Sta_heavyflv;
  TH1F *muPt_Sta_primary;
  TH1F *muPt_Sta_noise;
  TH1F *muPt_Sta_ghost;
//Good
  TH1F *muPt_Good_punch;
  TH1F *muPt_Good_lightflv;
  TH1F *muPt_Good_heavyflv;
  TH1F *muPt_Good_primary;
  TH1F *muPt_Good_noise;
  TH1F *muPt_Good_ghost;
//tight
  TH1F *muPt_Tight_punch;
  TH1F *muPt_Tight_lightflv;
  TH1F *muPt_Tight_heavyflv;
  TH1F *muPt_Tight_primary;
  TH1F *muPt_Tight_noise;
  TH1F *muPt_Tight_ghost;
//soft
  TH1F *muPt_Soft_punch;
  TH1F *muPt_Soft_lightflv;
  TH1F *muPt_Soft_heavyflv;
  TH1F *muPt_Soft_primary;
  TH1F *muPt_Soft_noise;
  TH1F *muPt_Soft_ghost;
//loose
  TH1F *muPt_Loose_punch;
  TH1F *muPt_Loose_lightflv;
  TH1F *muPt_Loose_heavyflv;
  TH1F *muPt_Loose_primary;
  TH1F *muPt_Loose_noise;
  TH1F *muPt_Loose_ghost;
///////////////////////////////////////////////////////
//eta
//global
  TH1F *muEta_Glb_punch;
  TH1F *muEta_Glb_lightflv;
  TH1F *muEta_Glb_heavyflv;
  TH1F *muEta_Glb_primary;
  TH1F *muEta_Glb_noise;
  TH1F *muEta_Glb_ghost;
//trk
  TH1F *muEta_Trk_punch;
  TH1F *muEta_Trk_lightflv;
  TH1F *muEta_Trk_heavyflv;
  TH1F *muEta_Trk_primary;
  TH1F *muEta_Trk_noise;
  TH1F *muEta_Trk_ghost;
//sta
  TH1F *muEta_Sta_punch;
  TH1F *muEta_Sta_lightflv;
  TH1F *muEta_Sta_heavyflv;
  TH1F *muEta_Sta_primary;
  TH1F *muEta_Sta_noise;
  TH1F *muEta_Sta_ghost;
//good
  TH1F *muEta_Good_punch;
  TH1F *muEta_Good_lightflv;
  TH1F *muEta_Good_heavyflv;
  TH1F *muEta_Good_primary;
  TH1F *muEta_Good_noise;
  TH1F *muEta_Good_ghost;
//tight
  TH1F *muEta_Tight_punch;
  TH1F *muEta_Tight_lightflv;
  TH1F *muEta_Tight_heavyflv;
  TH1F *muEta_Tight_primary;
  TH1F *muEta_Tight_noise;
  TH1F *muEta_Tight_ghost;
//soft
  TH1F *muEta_Soft_punch;
  TH1F *muEta_Soft_lightflv;
  TH1F *muEta_Soft_heavyflv;
  TH1F *muEta_Soft_primary;
  TH1F *muEta_Soft_noise;
  TH1F *muEta_Soft_ghost;
//loose
  TH1F *muEta_Loose_punch;
  TH1F *muEta_Loose_lightflv;
  TH1F *muEta_Loose_heavyflv;
  TH1F *muEta_Loose_primary;
  TH1F *muEta_Loose_noise;
  TH1F *muEta_Loose_ghost;
//phi
//global
  TH1F *muPhi_Glb_punch;
  TH1F *muPhi_Glb_lightflv;
  TH1F *muPhi_Glb_heavyflv;
  TH1F *muPhi_Glb_primary;
  TH1F *muPhi_Glb_noise;
  TH1F *muPhi_Glb_ghost;
//tracker
  TH1F *muPhi_Trk_punch;
  TH1F *muPhi_Trk_lightflv;
  TH1F *muPhi_Trk_heavyflv;
  TH1F *muPhi_Trk_primary;
  TH1F *muPhi_Trk_noise;
  TH1F *muPhi_Trk_ghost;
//sta
  TH1F *muPhi_Sta_punch;
  TH1F *muPhi_Sta_lightflv;
  TH1F *muPhi_Sta_heavyflv;
  TH1F *muPhi_Sta_primary;
  TH1F *muPhi_Sta_noise;
  TH1F *muPhi_Sta_ghost;
//good
  TH1F *muPhi_Good_punch;
  TH1F *muPhi_Good_lightflv;
  TH1F *muPhi_Good_heavyflv;
  TH1F *muPhi_Good_primary;
  TH1F *muPhi_Good_noise;
  TH1F *muPhi_Good_ghost;
//tight
  TH1F *muPhi_Tight_punch;
  TH1F *muPhi_Tight_lightflv;
  TH1F *muPhi_Tight_heavyflv;
  TH1F *muPhi_Tight_primary;
  TH1F *muPhi_Tight_noise;
  TH1F *muPhi_Tight_ghost;
//soft
  TH1F *muPhi_Soft_punch;
  TH1F *muPhi_Soft_lightflv;
  TH1F *muPhi_Soft_heavyflv;
  TH1F *muPhi_Soft_primary;
  TH1F *muPhi_Soft_noise;
  TH1F *muPhi_Soft_ghost;
//loose
  TH1F *muPhi_Loose_punch;
  TH1F *muPhi_Loose_lightflv;
  TH1F *muPhi_Loose_heavyflv;
  TH1F *muPhi_Loose_primary;
  TH1F *muPhi_Loose_noise;
  TH1F *muPhi_Loose_ghost;
// Matched stations////
  TH1F *muMStation_Glb_punch;
  TH1F *muMStation_Glb_lightflv;
  TH1F *muMStation_Glb_heavyflv;
  TH1F *muMStation_Glb_primary;
  TH1F *muMStation_Glb_noise;
  TH1F *muMStation_Glb_ghost;
//
  TH1F *muMStation_Sta_punch;
  TH1F *muMStation_Sta_lightflv;
  TH1F *muMStation_Sta_heavyflv;
  TH1F *muMStation_Sta_primary;
  TH1F *muMStation_Sta_noise;
  TH1F *muMStation_Sta_ghost;
//
  TH1F *muMStation_Tight_punch;
  TH1F *muMStation_Tight_lightflv;
  TH1F *muMStation_Tight_heavyflv;
  TH1F *muMStation_Tight_primary;
  TH1F *muMStation_Tight_noise;
  TH1F *muMStation_Tight_ghost;
//
  TH1F *muMStation_Loose_punch;
  TH1F *muMStation_Loose_lightflv;
  TH1F *muMStation_Loose_heavyflv;
  TH1F *muMStation_Loose_primary;
  TH1F *muMStation_Loose_noise;
  TH1F *muMStation_Loose_ghost;
//
  TH1F *muOOT_Glb_punch;
  TH1F *muOOT_Glb_lightflv;
  TH1F *muOOT_Glb_heavyflv;
  TH1F *muOOT_Glb_primary;
  TH1F *muOOT_Glb_noise;
  TH1F *muOOT_Glb_ghost;

  TH1F *muPtOOT_Glb_punch;
  TH1F *muPtOOT_Glb_lightflv;
  TH1F *muPtOOT_Glb_heavyflv;
  TH1F *muPtOOT_Glb_primary;
  TH1F *muPtOOT_Glb_noise;
  TH1F *muPtOOT_Glb_ghost;

  TH1F *muEtaOOT_Glb_punch;
  TH1F *muEtaOOT_Glb_lightflv;
  TH1F *muEtaOOT_Glb_heavyflv;
  TH1F *muEtaOOT_Glb_primary;
  TH1F *muEtaOOT_Glb_noise;
  TH1F *muEtaOOT_Glb_ghost;
 
  TH1F *muValidHitsOOT_Glb_punch;
  TH1F *muValidHitsOOT_Glb_lightflv;
  TH1F *muValidHitsOOT_Glb_heavyflv;
  TH1F *muValidHitsOOT_Glb_primary;
  TH1F *muValidHitsOOT_Glb_noise;
  TH1F *muValidHitsOOT_Glb_ghost;

  TH1F *muValidGemHitsOOT_Glb_punch;
  TH1F *muValidGemHitsOOT_Glb_lightflv;
  TH1F *muValidGemHitsOOT_Glb_heavyflv;
  TH1F *muValidGemHitsOOT_Glb_primary;
  TH1F *muValidGemHitsOOT_Glb_noise;
  TH1F *muValidGemHitsOOT_Glb_ghost;  

  TH1F *muInTime_Glb_punch;
  TH1F *muInTime_Glb_lightflv;
  TH1F *muInTime_Glb_heavyflv;
  TH1F *muInTime_Glb_primary;
  TH1F *muInTime_Glb_noise;
  TH1F *muInTime_Glb_ghost;

  TH1F *muPtInTime_Glb_punch;
  TH1F *muPtInTime_Glb_lightflv;
  TH1F *muPtInTime_Glb_heavyflv;
  TH1F *muPtInTime_Glb_primary;
  TH1F *muPtInTime_Glb_noise;
  TH1F *muPtInTime_Glb_ghost;

  TH1F *muEtaInTime_Glb_punch;
  TH1F *muEtaInTime_Glb_lightflv;
  TH1F *muEtaInTime_Glb_heavyflv;
  TH1F *muEtaInTime_Glb_primary;
  TH1F *muEtaInTime_Glb_noise;
  TH1F *muEtaInTime_Glb_ghost;

  TH1F *muValidHitsInTime_Glb_punch;
  TH1F *muValidHitsInTime_Glb_lightflv;
  TH1F *muValidHitsInTime_Glb_heavyflv;
  TH1F *muValidHitsInTime_Glb_primary;
  TH1F *muValidHitsInTime_Glb_noise;
  TH1F *muValidHitsInTime_Glb_ghost;

  TH1F *muValidGemHitsInTime_Glb_punch;
  TH1F *muValidGemHitsInTime_Glb_lightflv;
  TH1F *muValidGemHitsInTime_Glb_heavyflv;
  TH1F *muValidGemHitsInTime_Glb_primary;
  TH1F *muValidGemHitsInTime_Glb_noise;
  TH1F *muValidGemHitsInTime_Glb_ghost;  

//
  TH1F *muOOT_Sta_punch;
  TH1F *muOOT_Sta_lightflv;
  TH1F *muOOT_Sta_heavyflv;
  TH1F *muOOT_Sta_primary;
  TH1F *muOOT_Sta_noise;
  TH1F *muOOT_Sta_ghost;

  TH1F *muPtOOT_Sta_punch;
  TH1F *muPtOOT_Sta_lightflv;
  TH1F *muPtOOT_Sta_heavyflv;
  TH1F *muPtOOT_Sta_primary;
  TH1F *muPtOOT_Sta_noise;
  TH1F *muPtOOT_Sta_ghost;

  TH1F *muEtaOOT_Sta_punch;
  TH1F *muEtaOOT_Sta_lightflv;
  TH1F *muEtaOOT_Sta_heavyflv;
  TH1F *muEtaOOT_Sta_primary;
  TH1F *muEtaOOT_Sta_noise;
  TH1F *muEtaOOT_Sta_ghost;

  TH1F *muValidHitsOOT_Sta_punch;
  TH1F *muValidHitsOOT_Sta_lightflv;
  TH1F *muValidHitsOOT_Sta_heavyflv;
  TH1F *muValidHitsOOT_Sta_primary;
  TH1F *muValidHitsOOT_Sta_noise;
  TH1F *muValidHitsOOT_Sta_ghost;

  TH1F *muValidGemHitsOOT_Sta_punch;
  TH1F *muValidGemHitsOOT_Sta_lightflv;
  TH1F *muValidGemHitsOOT_Sta_heavyflv;
  TH1F *muValidGemHitsOOT_Sta_primary;
  TH1F *muValidGemHitsOOT_Sta_noise;
  TH1F *muValidGemHitsOOT_Sta_ghost;  

  TH1F *muInTime_Sta_punch;
  TH1F *muInTime_Sta_lightflv;
  TH1F *muInTime_Sta_heavyflv;
  TH1F *muInTime_Sta_primary;
  TH1F *muInTime_Sta_noise;
  TH1F *muInTime_Sta_ghost;

  TH1F *muPtInTime_Sta_punch;
  TH1F *muPtInTime_Sta_lightflv;
  TH1F *muPtInTime_Sta_heavyflv;
  TH1F *muPtInTime_Sta_primary;
  TH1F *muPtInTime_Sta_noise;
  TH1F *muPtInTime_Sta_ghost;

  TH1F *muEtaInTime_Sta_punch;
  TH1F *muEtaInTime_Sta_lightflv;
  TH1F *muEtaInTime_Sta_heavyflv;
  TH1F *muEtaInTime_Sta_primary;
  TH1F *muEtaInTime_Sta_noise;
  TH1F *muEtaInTime_Sta_ghost;

  TH1F *muValidHitsInTime_Sta_punch;
  TH1F *muValidHitsInTime_Sta_lightflv;
  TH1F *muValidHitsInTime_Sta_heavyflv;
  TH1F *muValidHitsInTime_Sta_primary;
  TH1F *muValidHitsInTime_Sta_noise;
  TH1F *muValidHitsInTime_Sta_ghost;

  TH1F *muValidGemHitsInTime_Sta_punch;
  TH1F *muValidGemHitsInTime_Sta_lightflv;
  TH1F *muValidGemHitsInTime_Sta_heavyflv;
  TH1F *muValidGemHitsInTime_Sta_primary;
  TH1F *muValidGemHitsInTime_Sta_noise;
  TH1F *muValidGemHitsInTime_Sta_ghost;  

  
  //
  TH1F *muOOT_Tight_punch;
  TH1F *muOOT_Tight_lightflv;
  TH1F *muOOT_Tight_heavyflv;
  TH1F *muOOT_Tight_primary;
  TH1F *muOOT_Tight_noise;
  TH1F *muOOT_Tight_ghost;
  
  TH1F *muPtOOT_Tight_punch;
  TH1F *muPtOOT_Tight_lightflv;
  TH1F *muPtOOT_Tight_heavyflv;
  TH1F *muPtOOT_Tight_primary;
  TH1F *muPtOOT_Tight_noise;
  TH1F *muPtOOT_Tight_ghost;

  TH1F *muEtaOOT_Tight_punch;
  TH1F *muEtaOOT_Tight_lightflv;
  TH1F *muEtaOOT_Tight_heavyflv;
  TH1F *muEtaOOT_Tight_primary;
  TH1F *muEtaOOT_Tight_noise;
  TH1F *muEtaOOT_Tight_ghost;
 
  TH1F *muValidHitsOOT_Tight_punch;
  TH1F *muValidHitsOOT_Tight_lightflv;
  TH1F *muValidHitsOOT_Tight_heavyflv;
  TH1F *muValidHitsOOT_Tight_primary;
  TH1F *muValidHitsOOT_Tight_noise;
  TH1F *muValidHitsOOT_Tight_ghost;

  TH1F *muValidGemHitsOOT_Tight_punch;
  TH1F *muValidGemHitsOOT_Tight_lightflv;
  TH1F *muValidGemHitsOOT_Tight_heavyflv;
  TH1F *muValidGemHitsOOT_Tight_primary;
  TH1F *muValidGemHitsOOT_Tight_noise;
  TH1F *muValidGemHitsOOT_Tight_ghost;  

  TH1F *muInTime_Tight_punch;
  TH1F *muInTime_Tight_lightflv;
  TH1F *muInTime_Tight_heavyflv;
  TH1F *muInTime_Tight_primary;
  TH1F *muInTime_Tight_noise;
  TH1F *muInTime_Tight_ghost;

  TH1F *muPtInTime_Tight_punch;
  TH1F *muPtInTime_Tight_lightflv;
  TH1F *muPtInTime_Tight_heavyflv;
  TH1F *muPtInTime_Tight_primary;
  TH1F *muPtInTime_Tight_noise;
  TH1F *muPtInTime_Tight_ghost;

  TH1F *muEtaInTime_Tight_punch;
  TH1F *muEtaInTime_Tight_lightflv;
  TH1F *muEtaInTime_Tight_heavyflv;
  TH1F *muEtaInTime_Tight_primary;
  TH1F *muEtaInTime_Tight_noise;
  TH1F *muEtaInTime_Tight_ghost;

  TH1F *muValidHitsInTime_Tight_punch;
  TH1F *muValidHitsInTime_Tight_lightflv;
  TH1F *muValidHitsInTime_Tight_heavyflv;
  TH1F *muValidHitsInTime_Tight_primary;
  TH1F *muValidHitsInTime_Tight_noise;
  TH1F *muValidHitsInTime_Tight_ghost;

  TH1F *muValidGemHitsInTime_Tight_punch;
  TH1F *muValidGemHitsInTime_Tight_lightflv;
  TH1F *muValidGemHitsInTime_Tight_heavyflv;
  TH1F *muValidGemHitsInTime_Tight_primary;
  TH1F *muValidGemHitsInTime_Tight_noise;
  TH1F *muValidGemHitsInTime_Tight_ghost;  

//

  TH1F *muOOT_Loose_punch;
  TH1F *muOOT_Loose_lightflv;
  TH1F *muOOT_Loose_heavyflv;
  TH1F *muOOT_Loose_primary;
  TH1F *muOOT_Loose_noise;
  TH1F *muOOT_Loose_ghost;

  TH1F *muPtOOT_Loose_punch;
  TH1F *muPtOOT_Loose_lightflv;
  TH1F *muPtOOT_Loose_heavyflv;
  TH1F *muPtOOT_Loose_primary;
  TH1F *muPtOOT_Loose_noise;
  TH1F *muPtOOT_Loose_ghost;

  TH1F *muEtaOOT_Loose_punch;
  TH1F *muEtaOOT_Loose_lightflv;
  TH1F *muEtaOOT_Loose_heavyflv;
  TH1F *muEtaOOT_Loose_primary;
  TH1F *muEtaOOT_Loose_noise;
  TH1F *muEtaOOT_Loose_ghost;

  TH1F *muValidHitsOOT_Loose_punch;
  TH1F *muValidHitsOOT_Loose_lightflv;
  TH1F *muValidHitsOOT_Loose_heavyflv;
  TH1F *muValidHitsOOT_Loose_primary;
  TH1F *muValidHitsOOT_Loose_noise;
  TH1F *muValidHitsOOT_Loose_ghost;

  TH1F *muValidGemHitsOOT_Loose_punch;
  TH1F *muValidGemHitsOOT_Loose_lightflv;
  TH1F *muValidGemHitsOOT_Loose_heavyflv;
  TH1F *muValidGemHitsOOT_Loose_primary;
  TH1F *muValidGemHitsOOT_Loose_noise;
  TH1F *muValidGemHitsOOT_Loose_ghost;  

  TH1F *muInTime_Loose_punch;
  TH1F *muInTime_Loose_lightflv;
  TH1F *muInTime_Loose_heavyflv;
  TH1F *muInTime_Loose_primary;
  TH1F *muInTime_Loose_noise;
  TH1F *muInTime_Loose_ghost;

  TH1F *muPtInTime_Loose_punch;
  TH1F *muPtInTime_Loose_lightflv;
  TH1F *muPtInTime_Loose_heavyflv;
  TH1F *muPtInTime_Loose_primary;
  TH1F *muPtInTime_Loose_noise;
  TH1F *muPtInTime_Loose_ghost;

  TH1F *muEtaInTime_Loose_punch;
  TH1F *muEtaInTime_Loose_lightflv;
  TH1F *muEtaInTime_Loose_heavyflv;
  TH1F *muEtaInTime_Loose_primary;
  TH1F *muEtaInTime_Loose_noise;
  TH1F *muEtaInTime_Loose_ghost;

  TH1F *muValidHitsInTime_Loose_punch;
  TH1F *muValidHitsInTime_Loose_lightflv;
  TH1F *muValidHitsInTime_Loose_heavyflv;
  TH1F *muValidHitsInTime_Loose_primary;
  TH1F *muValidHitsInTime_Loose_noise;
  TH1F *muValidHitsInTime_Loose_ghost;

  TH1F *muValidGemHitsInTime_Loose_punch;
  TH1F *muValidGemHitsInTime_Loose_lightflv;
  TH1F *muValidGemHitsInTime_Loose_heavyflv;
  TH1F *muValidGemHitsInTime_Loose_primary;
  TH1F *muValidGemHitsInTime_Loose_noise;
  TH1F *muValidGemHitsInTime_Loose_ghost;  

//
  TH1F *StaMuons_validhits;
  TH1F *StaMuons_validhits_atZerotime;
  TH1F *StaMuonsPt_NoValidhits;
  TH1F *StaMuonsPt_NoValidhits_atZeroTime;
  
// The RECO objects
  edm::InputTag trackingParticles_;
  edm::InputTag muons_;
  bool linkToGenParticles_;
  edm::InputTag genParticles_;
  std::string associatorLabel_;
  double decayRho_, decayAbsZ_;
/// Track to use
  int flavour(int pdgId) const;
  double muonPtCut_;
  double muonMinEtaCut_;
  double muonMaxEtaCut_;
  double isoCut_;
  edm::InputTag vxtTag;
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
bool MuonMCClassifAndAna::isTight(const edm::Event& iEvent, reco::Muon myMuon)
{
  bool results = false;
  if (myMuon.isGlobalMuon() && myMuon.isPFMuon())
  {
    if (myMuon.isIsolationValid() && (myMuon.isolationR03().sumPt) / myMuon.muonBestTrack()->pt() < isoCut_)
// if (myMuon.isIsolationValid() && (myMuon.isolationR03().sumPt) / myMuon.pt() < isoCut_)
    {
      reco::TrackRef trackRefGl = myMuon.globalTrack();
      reco::TrackRef trackRefIn = myMuon.innerTrack();
      reco::TrackRef trackRefBest = myMuon.muonBestTrack();
      if (trackRefGl.isNonnull() && trackRefIn.isNonnull() && trackRefBest.isNonnull())
      {
        edm::Handle < reco::VertexCollection > vertexHandle;
        iEvent.getByLabel(vxtTag, vertexHandle);
        const reco::VertexCollection* vertices = vertexHandle.product();
        if (myMuon.numberOfMatchedStations() > 1)
        {
          if (trackRefGl->hitPattern().numberOfValidMuonHits() > 0)
          {
            if (trackRefIn->hitPattern().pixelLayersWithMeasurement(3, 2) > 0)
            {
              if (trackRefIn->hitPattern().trackerLayersWithMeasurement() > 5)
              {
                if (myMuon.globalTrack()->normalizedChi2() < 10.)
// if (myMuon.track()->normalizedChi2() < 10.)
                {
                  if (vertices->size() != 0 && fabs(myMuon.muonBestTrack()->dxy((*vertices)[0].position())) < 0.2
                      && fabs(myMuon.muonBestTrack()->dz((*vertices)[0].position())) < 0.5)
                    results = true;
                }
              }
            }
          }
        }
      }
    }
  }
  return results;
}
bool MuonMCClassifAndAna::isSoft(reco::Muon myMuon, reco::VertexCollection::const_iterator myPVit,
    const reco::VertexCollection myVertices)
{
  bool results = false;
// if (myMuon.muonBestTrack().isNonnull() && myMuon.innerTrack().isNonnull())
  if (myMuon.innerTrack().isNonnull())
  {
    bool isGood = muon::isGoodMuon(myMuon, muon::TMOneStationTight);
    if (isGood)
    {
      bool trkLayMeas = myMuon.muonBestTrack()->hitPattern().trackerLayersWithMeasurement() > 5;
// bool trkLayMeas = myMuon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5;
      bool pxlLayMeas = myMuon.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0;
      bool quality = myMuon.innerTrack()->quality(reco::Track::highPurity);
      if (trkLayMeas && pxlLayMeas && quality)
      {
        reco::Vertex pv = myVertices.front();
        if ((pv.isValid()) && !(pv.isFake())
            && (fabs(pv.position().z()) <= 24 && fabs(pv.position().rho()) <= 2 && pv.ndof() > 4))
        {
          if (fabs(myMuon.innerTrack()->dxy(pv.position())) < 0.3 && fabs(myMuon.innerTrack()->dz(pv.position())) < 20)
            results = true;
        }
      }
    }
  }
  return results;
}
MuonMCClassifAndAna::MuonMCClassifAndAna(const edm::ParameterSet &iConfig) :
    trackingParticles_(iConfig.getParameter < edm::InputTag > ("simLabel")) //
        , muons_(iConfig.getParameter < edm::InputTag > ("muonLabel")) //
        , linkToGenParticles_(iConfig.getParameter<bool>("linkToGenParticles")) //
        , genParticles_(
        linkToGenParticles_ ? iConfig.getParameter < edm::InputTag > ("genParticles") : edm::InputTag("NONE")) //
            , associatorLabel_(iConfig.getParameter < std::string > ("muAssocLabel")) //
        , decayRho_(iConfig.getParameter<double>("decayRho")) //
        , decayAbsZ_(iConfig.getParameter<double>("decayAbsZ")) //
        , muonPtCut_(iConfig.getParameter<double>("muonPtCut")) //
        , muonMinEtaCut_(iConfig.getParameter<double>("muonMinEtaCut")) //
        , muonMaxEtaCut_(iConfig.getParameter<double>("muonMaxEtaCut")) //
        , isoCut_(iConfig.getParameter<double>("isoCut")) //
        , vxtTag(iConfig.getParameter < edm::InputTag > ("vxtTag")) //
{
  debug = false;
  glb_pch = 4;
  glb_l = 5;
  glb_h = 6;
  glb_prm = 7;
  glb_n = 8;
  glb_gst = 9;
  trk_pch = 12;
  trk_l = 13;
  trk_h = 14;
  trk_prm = 15;
  trk_n = 16;
  trk_gst = 17;
  sta_pch = 20;
  sta_l = 21;
  sta_h = 22;
  sta_prm = 23;
  sta_n = 24;
  sta_gst = 25;
  good_pch = 28;
  good_l = 29;
  good_h = 30;
  good_prm = 31;
  good_n = 32;
  good_gst = 33;
  tight_pch = 36;
  tight_l = 37;
  tight_h = 38;
  tight_prm = 39;
  tight_n = 40;
  tight_gst = 41;
  soft_pch = 44;
  soft_l = 45;
  soft_h = 46;
  soft_prm = 47;
  soft_n = 48;
  soft_gst = 49;
  loose_pch = 52;
  loose_l = 53;
  loose_h = 54;
  loose_prm = 55;
  loose_n = 56;
  loose_gst = 57;
}
MuonMCClassifAndAna::~MuonMCClassifAndAna()
{
}
//void MuonMCClassifAndAna::analyze(edm::Event& iEvent, const edm::EventSetup& iSetup)
void MuonMCClassifAndAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout << "New event start !" << std::endl;
   

  using namespace edm;
  using namespace std;
  using namespace reco;
  if (debug)
  {
    std::cout << "in MuonMCClassifAndAna !" << std::endl;
  }
//getting primary vertices
  edm::Handle < reco::VertexCollection > pvHandle;
  iEvent.getByLabel("offlinePrimaryVertices", pvHandle);
  int numberOfVertices = pvHandle->size();
  numb_PV->Fill(numberOfVertices);
  ///const reco::VertexCollection & vertices = *pvHandle.product();
  ///reco::VertexCollection::const_iterator myPV;
//////////Sim Vertex////////////////////
  /*
  edm::Handle < std::vector<SimVertex> > simVertices;
  iEvent.getByLabel("g4SimHits", simVertices);
  int numberOfSimVertices = simVertices->size();
  numb_SimPV->Fill(numberOfSimVertices);
  std::cout << "Simulated vertices !" << numberOfSimVertices << std::endl;
  const std::vector<SimVertex> & simver = *simVertices.product();
  std::vector<SimVertex>::const_iterator isimvtx;
  for (isimvtx = simver.begin(); isimvtx != simver.end(); isimvtx++)
  {
    float x_simPV = isimvtx->position().x();
    float y_simPV = isimvtx->position().y();
    float z_simPV = isimvtx->position().z();
    std::cout << "x position of Simulated vertices = " << x_simPV << std::endl;
    std::cout << "y position of Simulated vertices = " << y_simPV << std::endl;
    std::cout << "z position of Simulated vertices = " << z_simPV << std::endl;
  }
  */
///////////////////////////////////////////////////////////////////////
//getting muon collection
  edm::Handle < edm::View<reco::Muon> > myMuons;
  iEvent.getByLabel(muons_, myMuons);
  size_t nmu = myMuons->size();
  numb_MuPerEvent->Fill(nmu);
  if (nmu > 0)
    N_Muons->Fill(60);
  N_Muons_norm->Fill(60);
// numb_PV->Fill(numberOfVertices);
  N_Muons->Fill(62, numberOfVertices);
  N_Muons_norm->Fill(62, numberOfVertices);
//has valid primary vertices?
// bool hasValidPV = false;
  /*
   for (myPV = vertices.begin(); myPV != vertices.end(); myPV++)
   {
   reco::Vertex pv = (*myPV);
   if (pv.isValid() && (!(pv.isFake())))
   {
   // hasValidPV = true;
   //fill these only if we have at leas one valid PV
   numb_MuPerEvent->Fill(nmu);
   break;
   }
   }
   */
  edm::Handle < TrackingParticleCollection > trackingParticles;
  iEvent.getByLabel(trackingParticles_, trackingParticles);
/////////////////////////////////////////////////////////////////7
  const TrackingParticleCollection& trackingParticleCollection = *(trackingParticles.product());
  for (TrackingParticleCollection::size_type i = 0; i < trackingParticleCollection.size(); ++i)
  {
    TrackingParticleRef tp(trackingParticles, i);
//float Pt_tp = tp->pt();
//std::cout << "pt of tracking particles= " <<Pt_tp<< std::endl;
    Pt_TrakingParticle->Fill(tp->pt());
// std::cout << "pdg Id of tracking particles= " <<tp->pdgId()<< std::endl;
    if (abs(tp->pdgId()) == 13)
    {
      PtMuons_TrakingParticle->Fill(tp->pt());
      EtaMuons_TrakingParticle->Fill(tp->eta());
    }
  }
//////////////////////////////////////////////////////////////////
  edm::Handle < reco::GenParticleCollection > genParticlesforGlb;
  edm::Handle < reco::GenParticleCollection > genParticlesforSta;
  edm::Handle < reco::GenParticleCollection > genParticlesforTight;
  edm::Handle < reco::GenParticleCollection > genParticlesforLoose;
  if (linkToGenParticles_)
  {
    iEvent.getByLabel(genParticles_, genParticlesforGlb);
    iEvent.getByLabel(genParticles_, genParticlesforSta);
    iEvent.getByLabel(genParticles_, genParticlesforTight);
    iEvent.getByLabel(genParticles_, genParticlesforLoose);
  }
  edm::ESHandle < TrackAssociatorBase > associatorBase;
  iSetup.get<TrackAssociatorRecord>().get(associatorLabel_, associatorBase);
  const MuonAssociatorByHits * assoByHits = dynamic_cast<const MuonAssociatorByHits *>(associatorBase.product());
  if (assoByHits == 0)
    throw cms::Exception("Configuration") << "The Track Associator with label '" << associatorLabel_
        << "' is not a MuonAssociatorByHits.\n";
  MuonAssociatorByHits::MuonToSimCollection recSimColl;
  MuonAssociatorByHits::SimToMuonCollection simRecColl;
  edm::RefToBaseVector < reco::Muon > selAllMuons;
  edm::RefToBaseVector < reco::Muon > selGlbMuons;
  edm::RefToBaseVector < reco::Muon > selStaMuons;
  edm::RefToBaseVector < reco::Muon > selTightMuons;
  edm::RefToBaseVector < reco::Muon > selLooseMuons;
// int allMu_sel = 0;
// int glbMu_sel = 0;
// int trkMu_sel = 0;
// int staMu_sel = 0;
// int goodMu_sel = 0;
// int tightMu_sel = 0;
// int softMu_sel = 0;
// int looseMu_sel = 0;
///////////////test reco
//getting muon collection
  edm::Handle < vector<reco::Muon> > myrecoMuons;
  iEvent.getByLabel("muons", myrecoMuons);
// unsigned murecoSize = myMuons->size();
  vector<reco::Muon>::const_iterator amuon1;
  for (amuon1 = myrecoMuons->begin(); amuon1 != myrecoMuons->end(); amuon1++)
  {
    if (amuon1->pt() > muonPtCut_ && fabs(amuon1->eta()) >= muonMinEtaCut_ && fabs(amuon1->eta()) <= muonMaxEtaCut_)
    {
      allrecoMuonsEta->Fill(amuon1->eta());
      allrecoMuonsPt->Fill(amuon1->pt());
      allrecoMuonsPhi->Fill(amuon1->phi());
      reco::TrackRef trackRefGl = amuon1->globalTrack();
      if (trackRefGl.isNonnull())
      {
        histoValidMuHits->Fill(trackRefGl->hitPattern().numberOfValidMuonHits());
        histoValidMuGeMHits->Fill(trackRefGl->hitPattern().numberOfValidMuonGEMHits());
        histoMuStaWithValidHits->Fill(trackRefGl->hitPattern().muonStationsWithValidHits());
        histoGemStaWithValidHits->Fill(trackRefGl->hitPattern().gemStationsWithValidHits());
      }
      reco::TrackRef trackRefOut = amuon1->outerTrack();
      if (trackRefOut.isNonnull())
      {
        histoValidMuHitsOut->Fill(trackRefOut->hitPattern().numberOfValidMuonHits());
        histoValidMuGeMHitsOut->Fill(trackRefOut->hitPattern().numberOfValidMuonGEMHits());
        histoMuStaWithValidHitsOut->Fill(trackRefOut->hitPattern().muonStationsWithValidHits());
        histoGemStaWithValidHitsOut->Fill(trackRefOut->hitPattern().gemStationsWithValidHits());
      }
    }
  }
//////////////test reco end
  for (size_t i = 0, n = myMuons->size(); i < n; ++i)
  {
// evtsMu++;
    N_Muons->Fill(64); //all muons
    N_Muons_norm->Fill(64); //all muons
   
    edm::RefToBase < reco::Muon > rmu = myMuons->refAt(i);
       
    if ((*myMuons)[i].pt() > muonPtCut_ && fabs((*myMuons)[i].eta()) >= muonMinEtaCut_
        && fabs((*myMuons)[i].eta()) <= muonMaxEtaCut_)
    {
// allMu_sel++;
      N_Muons->Fill(1); //all selected muons
      N_Muons_norm->Fill(1); //all selected muons
      allMuonsEta->Fill((*myMuons)[i].eta());
      allMuonsPt->Fill((*myMuons)[i].pt());
      allMuonsPhi->Fill((*myMuons)[i].phi());
      selAllMuons.push_back(rmu);
      if (debug)
        std::cout << "event: run\t" << (iEvent.id()).run() << "\tevent #\t" << (iEvent.id()).event() << std::endl; //
      if (debug)
        std::cout << "I am Selected muon!" << std::endl;
     
    }        

      if ((*myMuons)[i].isGlobalMuon())
       {
// glbMu_sel++;
     reco::TrackRef glbtrackref;
      if ((*myMuons)[i].globalTrack().isNonnull()) glbtrackref = (*myMuons)[i].globalTrack();
     
      if (glbtrackref->pt() > muonPtCut_ && fabs(glbtrackref->eta()) >= muonMinEtaCut_
        && fabs(glbtrackref->eta()) <= muonMaxEtaCut_)
       {
   
        selGlbMuons.push_back(rmu);
        N_Muons->Fill(3);
        N_Muons_norm->Fill(3);
        GlobalMuonsEta->Fill(glbtrackref->eta());
        GlobalMuonsPt->Fill(glbtrackref->pt());
        GlobalMuonsPhi->Fill(glbtrackref->phi());
        GlobalMuons_matchedstation->Fill((*myMuons)[i].numberOfMatchedStations());
        if (debug)
          std::cout << "I am Global muon!" << std::endl;
      }
    }
  
      if ((*myMuons)[i].isStandAloneMuon())
      {
      reco::TrackRef statrackref;
      if ((*myMuons)[i].outerTrack().isNonnull()) statrackref = (*myMuons)[i].outerTrack();

      if (statrackref->pt() > muonPtCut_ && fabs(statrackref->eta()) >= muonMinEtaCut_
        && fabs(statrackref->eta()) <= muonMaxEtaCut_)
       {

        selStaMuons.push_back(rmu);
        N_Muons->Fill(19);
        N_Muons_norm->Fill(19);
        StaMuonsEta->Fill(statrackref->eta());
        StaMuonsPt->Fill(statrackref->pt());
        StaMuonsPhi->Fill(statrackref->phi());
        StaMuons_matchedstation->Fill((*myMuons)[i].numberOfMatchedStations());

        StaMuons_validhits->Fill(statrackref->hitPattern().numberOfValidMuonHits());
        if(statrackref->hitPattern().numberOfValidMuonHits() == 0)StaMuonsPt_NoValidhits->Fill(statrackref->pt());
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if(fabs((*myMuons)[i].time().timeAtIpInOut) > 25+fabs((*myMuons)[i].time().timeAtIpInOutErr)) {
        StaMuonsOOT->Fill((*myMuons)[i].time().timeAtIpInOut);
        StaMuonsPt_OOT->Fill(statrackref->pt());
        StaMuonsEta_OOT->Fill(statrackref->eta());
        StaMuonsValidHits_OOT->Fill(statrackref->hitPattern().numberOfValidMuonHits());
        StaMuonsValidGemHits_OOT->Fill(statrackref->hitPattern().numberOfValidMuonGEMHits());
        }
     
        else{
        StaMuonsInTime->Fill((*myMuons)[i].time().timeAtIpInOut);
        StaMuonsPt_InTime->Fill(statrackref->pt());
        StaMuonsEta_InTime->Fill(statrackref->eta());
        StaMuonsValidHits_InTime->Fill(statrackref->hitPattern().numberOfValidMuonHits());
        StaMuonsValidGemHits_InTime->Fill(statrackref->hitPattern().numberOfValidMuonGEMHits());
        } 

       if((*myMuons)[i].time().timeAtIpInOut == 0)  {
        StaMuons_validhits_atZerotime->Fill(statrackref->hitPattern().numberOfValidMuonHits());
        if(statrackref->hitPattern().numberOfValidMuonHits() == 0)StaMuonsPt_NoValidhits_atZeroTime->Fill(statrackref->pt());
        } 

        if (debug)
          std::cout << "I am Standalone muon!" << std::endl;
      }
     }
    

      bool tight = isTight(iEvent, (*myMuons)[i]);
      if (tight)
      {

     if ((*myMuons)[i].pt() > muonPtCut_ && fabs((*myMuons)[i].eta()) >= muonMinEtaCut_
        && fabs((*myMuons)[i].eta()) <= muonMaxEtaCut_)
       {
        N_Muons->Fill(35);
        N_Muons_norm->Fill(35);
        selTightMuons.push_back(rmu);
        TightMuonsEta->Fill((*myMuons)[i].eta());
        TightMuonsPt->Fill((*myMuons)[i].pt());
        TightMuonsPhi->Fill((*myMuons)[i].phi());
        TightMuons_matchedstation->Fill((*myMuons)[i].numberOfMatchedStations());
/////test reco
        reco::TrackRef trackRefGl = (*myMuons)[i].globalTrack();
        if (trackRefGl.isNonnull())
        {
          histoValidMuHitsTight->Fill(trackRefGl->hitPattern().numberOfValidMuonHits());
          histoValidMuGeMHitsTight->Fill(trackRefGl->hitPattern().numberOfValidMuonGEMHits());
          histoMuStaWithValidHitsTight->Fill(trackRefGl->hitPattern().muonStationsWithValidHits());
          histoGemStaWithValidHitsTight->Fill(trackRefGl->hitPattern().gemStationsWithValidHits());
        }
        reco::TrackRef trackRefOut = (*myMuons)[i].outerTrack();
        if (trackRefOut.isNonnull())
        {
          histoValidMuHitsOutTight->Fill(trackRefOut->hitPattern().numberOfValidMuonHits());
          histoValidMuGeMHitsOutTight->Fill(trackRefOut->hitPattern().numberOfValidMuonGEMHits());
          histoMuStaWithValidHitsOutTight->Fill(trackRefOut->hitPattern().muonStationsWithValidHits());
          histoGemStaWithValidHitsOutTight->Fill(trackRefOut->hitPattern().gemStationsWithValidHits());
        }
///test reco end
        if (debug)
          std::cout << "I am Tight muon!" << std::endl;
      }
     }
 

      if (muon::isLooseMuon ((*myMuons)[i]))
      {

      if ((*myMuons)[i].pt() > muonPtCut_ && fabs((*myMuons)[i].eta()) >= muonMinEtaCut_
        && fabs((*myMuons)[i].eta()) <= muonMaxEtaCut_)
       {
        N_Muons->Fill(51);
        N_Muons_norm->Fill(51);
        selLooseMuons.push_back(rmu);
        LooseMuonsEta->Fill((*myMuons)[i].eta());
        LooseMuonsPt->Fill((*myMuons)[i].pt());
        LooseMuonsPhi->Fill((*myMuons)[i].phi());
        LooseMuons_matchedstation->Fill((*myMuons)[i].numberOfMatchedStations());
/////////////////////////////////////////////////////////////////////////////////////
        reco::TrackRef trackRefGl = (*myMuons)[i].globalTrack();
        if (trackRefGl.isNonnull())
        {
          histoValidMuHitsLoose->Fill(trackRefGl->hitPattern().numberOfValidMuonHits());
          histoValidMuGeMHitsLoose->Fill(trackRefGl->hitPattern().numberOfValidMuonGEMHits());
          histoMuStaWithValidHitsLoose->Fill(trackRefGl->hitPattern().muonStationsWithValidHits());
          histoGemStaWithValidHitsLoose->Fill(trackRefGl->hitPattern().gemStationsWithValidHits());
        }
        reco::TrackRef trackRefOut = (*myMuons)[i].outerTrack();
        if (trackRefOut.isNonnull())
        {
          histoValidMuHitsOutLoose->Fill(trackRefOut->hitPattern().numberOfValidMuonHits());
          histoValidMuGeMHitsOutLoose->Fill(trackRefOut->hitPattern().numberOfValidMuonGEMHits());
          histoMuStaWithValidHitsOutLoose->Fill(trackRefOut->hitPattern().muonStationsWithValidHits());
          histoGemStaWithValidHitsOutLoose->Fill(trackRefOut->hitPattern().gemStationsWithValidHits());
        }
//////////////////////////////////////////////////////////////////////////////////////////////
        if (debug)
          std::cout << "I am Loose muon!" << std::endl;
      }
     
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
//global muon asso by global track
  MuonAssociatorByHits::MuonToSimCollection recSimColl_glb;
  MuonAssociatorByHits::SimToMuonCollection simRecColl_glb;
  edm::LogVerbatim("MuonMCClassifAndAna") << "\n Global Muon association by global track"; //
  assoByHits->associateMuons(recSimColl_glb, simRecColl_glb, selGlbMuons, MuonAssociatorByHits::OuterTk, allTPs,
      &iEvent, &iSetup); //
//standalone muon
  MuonAssociatorByHits::MuonToSimCollection recSimColl_sta;
  MuonAssociatorByHits::SimToMuonCollection simRecColl_sta;
  edm::LogVerbatim("MuonMCClassifAndAna") << "\n Standalone Muon association by outer track "; //
  assoByHits->associateMuons(recSimColl_sta, simRecColl_sta, selStaMuons, MuonAssociatorByHits::OuterTk, allTPs,
      &iEvent, &iSetup); //
//tight muon asso by global track
  MuonAssociatorByHits::MuonToSimCollection recSimColl_tight;
  MuonAssociatorByHits::SimToMuonCollection simRecColl_tight;
  edm::LogVerbatim("MuonMCClassifAndAna") << "\n Tight Muon association by global track "; //
  assoByHits->associateMuons(recSimColl_tight, simRecColl_tight, selTightMuons, MuonAssociatorByHits::OuterTk, allTPs,
      &iEvent, &iSetup); //
//Loose
  MuonAssociatorByHits::MuonToSimCollection recSimColl_loose;
  MuonAssociatorByHits::SimToMuonCollection simRecColl_loose;
  edm::LogVerbatim("MuonMCClassifAndAna") << "\n Loose Muon association by inner track "; //
  assoByHits->associateMuons(recSimColl_loose, simRecColl_loose, selLooseMuons, MuonAssociatorByHits::OuterTk, allTPs,
      &iEvent, &iSetup); //
  
if (debug)
    std::cout << "there are " << nmu << " reco::Muons" << std::endl;
  if (debug)
    std::cout << "=================================" << std::endl;
  edm::LogVerbatim("MuonMCClassifAndAna") << "\n There are " << nmu << " reco::Muons.";
  std::auto_ptr<reco::GenParticleCollection> secondariesforGlb;
  std::auto_ptr<reco::GenParticleCollection> secondariesforSta;
  std::auto_ptr<reco::GenParticleCollection> secondariesforTight;
  std::auto_ptr<reco::GenParticleCollection> secondariesforLoose;
  std::vector<int> muToPrimaryforGlb(nmu, -1), muToSecondaryforGlb(nmu, -1);
  std::vector<int> muToPrimaryforSta(nmu, -1), muToSecondaryforSta(nmu, -1);
  std::vector<int> muToPrimaryforTight(nmu, -1), muToSecondaryforTight(nmu, -1);
  std::vector<int> muToPrimaryforLoose(nmu, -1), muToSecondaryforLoose(nmu, -1);
  if (linkToGenParticles_)
  {
    secondariesforGlb.reset(new reco::GenParticleCollection());
    secondariesforSta.reset(new reco::GenParticleCollection());
    secondariesforTight.reset(new reco::GenParticleCollection());
    secondariesforLoose.reset(new reco::GenParticleCollection());
  }
//for Global MUONS
  std::vector<int> classifgl(nmu, 0), extgl(nmu, 0);
  std::vector<int> hitsPdgIdgl(nmu, 0), momPdgIdgl(nmu, 0), gmomPdgIdgl(nmu, 0), momStatusgl(nmu, 0);
  std::vector<int> flavgl(nmu, 0), momFlavgl(nmu, 0), gmomFlavgl(nmu, 0), hmomFlavgl(nmu, 0);
  std::vector<int> tpIdgl(nmu, -1);
  std::vector<float> prodRhogl(nmu, 0.0), prodZgl(nmu, 0.0), momRhogl(nmu, 0.0), momZgl(nmu, 0.0);
  std::vector<float> tpAssoQualitygl(nmu, -1);
  myClassification(nmu, myMuons, selGlbMuons, recSimColl_glb, simRecColl_glb, 
      classifgl, extgl, hitsPdgIdgl, momPdgIdgl, gmomPdgIdgl, momStatusgl, flavgl, momFlavgl,
      gmomFlavgl, hmomFlavgl, tpIdgl, prodRhogl, prodZgl, momRhogl, momZgl, tpAssoQualitygl, secondariesforGlb,
      muToPrimaryforGlb, muToSecondaryforGlb, genParticlesforGlb, glb_pch, glb_l, glb_h, glb_prm, glb_n, glb_gst,
      muPt_Glb_punch, muPt_Glb_lightflv, muPt_Glb_heavyflv, muPt_Glb_primary, muPt_Glb_noise, muPt_Glb_ghost,
      muEta_Glb_punch, muEta_Glb_lightflv, muEta_Glb_heavyflv, muEta_Glb_primary, muEta_Glb_noise, muEta_Glb_ghost,
      muPhi_Glb_punch, muPhi_Glb_lightflv, muPhi_Glb_heavyflv, muPhi_Glb_primary, muPhi_Glb_noise, muPhi_Glb_ghost,
      muMStation_Glb_punch, muMStation_Glb_lightflv, muMStation_Glb_heavyflv, muMStation_Glb_primary,
      muMStation_Glb_noise, muMStation_Glb_ghost, 
      muOOT_Glb_punch, muOOT_Glb_lightflv, muOOT_Glb_heavyflv, muOOT_Glb_primary, muOOT_Glb_noise,
      muOOT_Glb_ghost, muPtOOT_Glb_punch, muPtOOT_Glb_lightflv, muPtOOT_Glb_heavyflv, muPtOOT_Glb_primary, muPtOOT_Glb_noise,
      muPtOOT_Glb_ghost, muEtaOOT_Glb_punch, muEtaOOT_Glb_lightflv, muEtaOOT_Glb_heavyflv, muEtaOOT_Glb_primary,
      muEtaOOT_Glb_noise, muEtaOOT_Glb_ghost, muValidHitsOOT_Glb_punch, muValidHitsOOT_Glb_lightflv, muValidHitsOOT_Glb_heavyflv,
      muValidHitsOOT_Glb_primary, muValidHitsOOT_Glb_noise, muValidHitsOOT_Glb_ghost, muValidGemHitsOOT_Glb_punch,
      muValidGemHitsOOT_Glb_lightflv, muValidGemHitsOOT_Glb_heavyflv, muValidGemHitsOOT_Glb_primary, muValidGemHitsOOT_Glb_noise,
      muValidGemHitsOOT_Glb_ghost, muInTime_Glb_punch, muInTime_Glb_lightflv, muInTime_Glb_heavyflv,                         
      muInTime_Glb_primary, muInTime_Glb_noise, muInTime_Glb_ghost, muPtInTime_Glb_punch, muPtInTime_Glb_lightflv,                         
      muPtInTime_Glb_heavyflv, muPtInTime_Glb_primary, muPtInTime_Glb_noise, muPtInTime_Glb_ghost, muEtaInTime_Glb_punch,                       
      muEtaInTime_Glb_lightflv, muEtaInTime_Glb_heavyflv, muEtaInTime_Glb_primary, muEtaInTime_Glb_noise,              
      muEtaInTime_Glb_ghost, muValidHitsInTime_Glb_punch, muValidHitsInTime_Glb_lightflv, muValidHitsInTime_Glb_heavyflv,
      muValidHitsInTime_Glb_primary, muValidHitsInTime_Glb_noise, muValidHitsInTime_Glb_ghost, muValidGemHitsInTime_Glb_punch,
      muValidGemHitsInTime_Glb_lightflv, muValidGemHitsInTime_Glb_heavyflv, muValidGemHitsInTime_Glb_primary,
      muValidGemHitsInTime_Glb_noise, muValidGemHitsInTime_Glb_ghost, 1);
//for STANDALONE MUONS
  std::vector<int> classifs(nmu, 0), exts(nmu, 0);
  std::vector<int> hitsPdgIds(nmu, 0), momPdgIds(nmu, 0), gmomPdgIds(nmu, 0), momStatuss(nmu, 0);
  std::vector<int> flavs(nmu, 0), momFlavs(nmu, 0), gmomFlavs(nmu, 0), hmomFlavs(nmu, 0);
  std::vector<int> tpIds(nmu, -1);
  std::vector<float> prodRhos(nmu, 0.0), prodZs(nmu, 0.0), momRhos(nmu, 0.0), momZs(nmu, 0.0);
  std::vector<float> tpAssoQualitys(nmu, -1);
  myClassification(nmu, myMuons, selStaMuons, recSimColl_sta, simRecColl_sta, 
      classifs, exts, hitsPdgIds, momPdgIds, gmomPdgIds, momStatuss, flavs, momFlavs, gmomFlavs,
      hmomFlavs, tpIds, prodRhos, prodZs, momRhos, momZs, tpAssoQualitys, secondariesforSta, muToPrimaryforSta,
      muToSecondaryforSta, genParticlesforSta, sta_pch, sta_l, sta_h, sta_prm, sta_n, sta_gst, muPt_Sta_punch,
      muPt_Sta_lightflv, muPt_Sta_heavyflv, muPt_Sta_primary, muPt_Sta_noise, muPt_Sta_ghost, muEta_Sta_punch,
      muEta_Sta_lightflv, muEta_Sta_heavyflv, muEta_Sta_primary, muEta_Sta_noise, muEta_Sta_ghost, muPhi_Sta_punch,
      muPhi_Sta_lightflv, muPhi_Sta_heavyflv, muPhi_Sta_primary, muPhi_Sta_noise, muPhi_Sta_ghost, muMStation_Sta_punch, 
      muMStation_Sta_lightflv, muMStation_Sta_heavyflv, muMStation_Sta_primary, muMStation_Sta_noise, muMStation_Sta_ghost, 
      muOOT_Sta_punch, muOOT_Sta_lightflv, muOOT_Sta_heavyflv, muOOT_Sta_primary, muOOT_Sta_noise,
      muOOT_Sta_ghost, muPtOOT_Sta_punch, muPtOOT_Sta_lightflv, muPtOOT_Sta_heavyflv, muPtOOT_Sta_primary, muPtOOT_Sta_noise,
      muPtOOT_Sta_ghost, muEtaOOT_Sta_punch, muEtaOOT_Sta_lightflv, muEtaOOT_Sta_heavyflv, muEtaOOT_Sta_primary,
      muEtaOOT_Sta_noise, muEtaOOT_Sta_ghost, muValidHitsOOT_Sta_punch, muValidHitsOOT_Sta_lightflv, muValidHitsOOT_Sta_heavyflv, 
      muValidHitsOOT_Sta_primary, muValidHitsOOT_Sta_noise, muValidHitsOOT_Sta_ghost, muValidGemHitsOOT_Sta_punch,
      muValidGemHitsOOT_Sta_lightflv, muValidGemHitsOOT_Sta_heavyflv, muValidGemHitsOOT_Sta_primary, muValidGemHitsOOT_Sta_noise,
      muValidGemHitsOOT_Sta_ghost, muInTime_Sta_punch, muInTime_Sta_lightflv, muInTime_Sta_heavyflv,                         
      muInTime_Sta_primary, muInTime_Sta_noise, muInTime_Sta_ghost, muPtInTime_Sta_punch, muPtInTime_Sta_lightflv,                         
      muPtInTime_Sta_heavyflv, muPtInTime_Sta_primary, muPtInTime_Sta_noise, muPtInTime_Sta_ghost, muEtaInTime_Sta_punch,                       
      muEtaInTime_Sta_lightflv, muEtaInTime_Sta_heavyflv, muEtaInTime_Sta_primary, muEtaInTime_Sta_noise,   
      muEtaInTime_Sta_ghost, muValidHitsInTime_Sta_punch, muValidHitsInTime_Sta_lightflv, muValidHitsInTime_Sta_heavyflv,
      muValidHitsInTime_Sta_primary, muValidHitsInTime_Sta_noise, muValidHitsInTime_Sta_ghost, muValidGemHitsInTime_Sta_punch,
      muValidGemHitsInTime_Sta_lightflv, muValidGemHitsInTime_Sta_heavyflv, muValidGemHitsInTime_Sta_primary, 
      muValidGemHitsInTime_Sta_noise, muValidGemHitsInTime_Sta_ghost, 2);
//for TIGHT MUONS
  std::vector<int> classift(nmu, 0), extt(nmu, 0);
  std::vector<int> hitsPdgIdt(nmu, 0), momPdgIdt(nmu, 0), gmomPdgIdt(nmu, 0), momStatust(nmu, 0);
  std::vector<int> flavt(nmu, 0), momFlavt(nmu, 0), gmomFlavt(nmu, 0), hmomFlavt(nmu, 0);
  std::vector<int> tpIdt(nmu, -1);
  std::vector<float> prodRhot(nmu, 0.0), prodZt(nmu, 0.0), momRhot(nmu, 0.0), momZt(nmu, 0.0);
  std::vector<float> tpAssoQualityt(nmu, -1);
  myClassification(nmu, myMuons, selTightMuons, recSimColl_tight, simRecColl_tight, 
      classift, extt, hitsPdgIdt, momPdgIdt, gmomPdgIdt, momStatust, flavt, momFlavt,
      gmomFlavt, hmomFlavt, tpIdt, prodRhot, prodZt, momRhot, momZt, tpAssoQualityt, secondariesforTight,
      muToPrimaryforTight, muToSecondaryforTight, genParticlesforTight, tight_pch, tight_l, tight_h, tight_prm, tight_n,
      tight_gst, muPt_Tight_punch, muPt_Tight_lightflv, muPt_Tight_heavyflv, muPt_Tight_primary, muPt_Tight_noise,
      muPt_Tight_ghost, muEta_Tight_punch, muEta_Tight_lightflv, muEta_Tight_heavyflv, muEta_Tight_primary,
      muEta_Tight_noise, muEta_Tight_ghost, muPhi_Tight_punch, muPhi_Tight_lightflv, muPhi_Tight_heavyflv,
      muPhi_Tight_primary, muPhi_Tight_noise, muPhi_Tight_ghost, muMStation_Tight_punch, muMStation_Tight_lightflv, 
      muMStation_Tight_heavyflv, muMStation_Tight_primary, muMStation_Tight_noise, muMStation_Tight_ghost,
      muOOT_Tight_punch, muOOT_Tight_lightflv, muOOT_Tight_heavyflv, muOOT_Tight_primary, muOOT_Tight_noise, muOOT_Tight_ghost,
      muPtOOT_Tight_punch, muPtOOT_Tight_lightflv, muPtOOT_Tight_heavyflv, muPtOOT_Tight_primary, muPtOOT_Tight_noise, muPtOOT_Tight_ghost, 
      muEtaOOT_Tight_punch, muEtaOOT_Tight_lightflv, muEtaOOT_Tight_heavyflv, muEtaOOT_Tight_primary, muEtaOOT_Tight_noise,
      muEtaOOT_Tight_ghost, muValidHitsOOT_Tight_punch, muValidHitsOOT_Tight_lightflv, muValidHitsOOT_Tight_heavyflv, 
      muValidHitsOOT_Tight_primary, muValidHitsOOT_Tight_noise, muValidHitsOOT_Tight_ghost, muValidGemHitsOOT_Tight_punch, 
      muValidGemHitsOOT_Tight_lightflv, muValidGemHitsOOT_Tight_heavyflv, muValidGemHitsOOT_Tight_primary, muValidGemHitsOOT_Tight_noise, 
      muValidGemHitsOOT_Tight_ghost, muInTime_Tight_punch, muInTime_Tight_lightflv, muInTime_Tight_heavyflv, muInTime_Tight_primary, 
      muInTime_Tight_noise, muInTime_Tight_ghost, muPtInTime_Tight_punch, muPtInTime_Tight_lightflv, muPtInTime_Tight_heavyflv,
      muPtInTime_Tight_primary, muPtInTime_Tight_noise, muPtInTime_Tight_ghost, muEtaInTime_Tight_punch, muEtaInTime_Tight_lightflv, 
      muEtaInTime_Tight_heavyflv, muEtaInTime_Tight_primary, muEtaInTime_Tight_noise, muEtaInTime_Tight_ghost, muValidHitsInTime_Tight_punch,
      muValidHitsInTime_Tight_lightflv, muValidHitsInTime_Tight_heavyflv, muValidHitsInTime_Tight_primary, muValidHitsInTime_Tight_noise, 
      muValidHitsInTime_Tight_ghost, muValidGemHitsInTime_Tight_punch, muValidGemHitsInTime_Tight_lightflv,
      muValidGemHitsInTime_Tight_heavyflv, muValidGemHitsInTime_Tight_primary, muValidGemHitsInTime_Tight_noise,
      muValidGemHitsInTime_Tight_ghost, 3);
//for LOOSe MUONS
  std::vector<int> classifls(nmu, 0), extls(nmu, 0);
  std::vector<int> hitsPdgIdls(nmu, 0), momPdgIdls(nmu, 0), gmomPdgIdls(nmu, 0), momStatusls(nmu, 0);
  std::vector<int> flavls(nmu, 0), momFlavls(nmu, 0), gmomFlavls(nmu, 0), hmomFlavls(nmu, 0);
  std::vector<int> tpIdls(nmu, -1);
  std::vector<float> prodRhosl(nmu, 0.0), prodZls(nmu, 0.0), momRhols(nmu, 0.0), momZls(nmu, 0.0);
  std::vector<float> tpAssoQualityls(nmu, -1);
  myClassification(nmu, myMuons, selLooseMuons, recSimColl_loose, simRecColl_loose,
      classifls, extls, hitsPdgIdls, momPdgIdls, gmomPdgIdls, momStatusls, flavs, momFlavls,
      gmomFlavls, hmomFlavls, tpIdls, prodRhosl, prodZls, momRhols, momZls, tpAssoQualityls, secondariesforLoose,
      muToPrimaryforLoose, muToSecondaryforLoose, genParticlesforLoose, loose_pch, loose_l, loose_h, loose_prm, loose_n,
      loose_gst, muPt_Loose_punch, muPt_Loose_lightflv, muPt_Loose_heavyflv, muPt_Loose_primary, muPt_Loose_noise,
      muPt_Loose_ghost, muEta_Loose_punch, muEta_Loose_lightflv, muEta_Loose_heavyflv, muEta_Loose_primary,
      muEta_Loose_noise, muEta_Loose_ghost, muPhi_Loose_punch, muPhi_Loose_lightflv, muPhi_Loose_heavyflv,
      muPhi_Loose_primary, muPhi_Loose_noise, muPhi_Loose_ghost, muMStation_Loose_punch, muMStation_Loose_lightflv,
      muMStation_Loose_heavyflv, muMStation_Loose_primary, muMStation_Loose_noise, muMStation_Loose_ghost,
      muOOT_Loose_punch, muOOT_Loose_lightflv, muOOT_Loose_heavyflv, muOOT_Loose_primary, muOOT_Loose_noise,
      muOOT_Loose_ghost, muPtOOT_Loose_punch, muPtOOT_Loose_lightflv, muPtOOT_Loose_heavyflv, muPtOOT_Loose_primary, muPtOOT_Loose_noise,
      muPtOOT_Loose_ghost, muEtaOOT_Loose_punch, muEtaOOT_Loose_lightflv, muEtaOOT_Loose_heavyflv, muEtaOOT_Loose_primary,
      muEtaOOT_Loose_noise, muEtaOOT_Loose_ghost, muValidHitsOOT_Loose_punch, muValidHitsOOT_Loose_lightflv, muValidHitsOOT_Loose_heavyflv,
      muValidHitsOOT_Loose_primary, muValidHitsOOT_Loose_noise, muValidHitsOOT_Loose_ghost, muValidGemHitsOOT_Loose_punch, 
      muValidGemHitsOOT_Loose_lightflv, muValidGemHitsOOT_Loose_heavyflv, muValidGemHitsOOT_Loose_primary, muValidGemHitsOOT_Loose_noise,
      muValidGemHitsOOT_Loose_ghost, muInTime_Loose_punch, muInTime_Loose_lightflv, muInTime_Loose_heavyflv, 
      muInTime_Loose_primary, muInTime_Loose_noise, muInTime_Loose_ghost, muPtInTime_Loose_punch, muPtInTime_Loose_lightflv, 
      muPtInTime_Loose_heavyflv, muPtInTime_Loose_primary, muPtInTime_Loose_noise, muPtInTime_Loose_ghost, muEtaInTime_Loose_punch,
      muEtaInTime_Loose_lightflv, muEtaInTime_Loose_heavyflv, muEtaInTime_Loose_primary, muEtaInTime_Loose_noise,
      muEtaInTime_Loose_ghost, muValidHitsInTime_Loose_punch, muValidHitsInTime_Loose_lightflv, muValidHitsInTime_Loose_heavyflv, 
      muValidHitsInTime_Loose_primary, muValidHitsInTime_Loose_noise, muValidHitsInTime_Loose_ghost, muValidGemHitsInTime_Loose_punch,
      muValidGemHitsInTime_Loose_lightflv, muValidGemHitsInTime_Loose_heavyflv, muValidGemHitsInTime_Loose_primary,
      muValidGemHitsInTime_Loose_noise, muValidGemHitsInTime_Loose_ghost, 4);
} //end Analyzer

void MuonMCClassifAndAna::myClassification(size_t myNmu, edm::Handle<edm::View<reco::Muon> > myMuons_handle,
    edm::RefToBaseVector<reco::Muon> selMuons, MuonAssociatorByHits::MuonToSimCollection myRecSimColl,
    MuonAssociatorByHits::SimToMuonCollection mySimRecColl,
    std::vector<int> &myClassif, std::vector<int> &myExt,
    std::vector<int> &myHitsPdgId, std::vector<int> &myMomPdgId, std::vector<int> &myGmomPdgId,
    std::vector<int> &myMomStatus, std::vector<int> &myFlav, std::vector<int> &myMomFlav, std::vector<int> &myGmomFlav,
    std::vector<int> &myHmomFlav, std::vector<int> &myTpId, std::vector<float> &myProdRho, std::vector<float> &myProdZ,
    std::vector<float> &myMomRho, std::vector<float> &myMomZ, std::vector<float> &myTpAssoQuality,
    std::auto_ptr<reco::GenParticleCollection> &secondaries, std::vector<int> &muToPrimary,
    std::vector<int> &muToSecondary, edm::Handle<reco::GenParticleCollection> &genParticles, int &binp, int &binl,
    int &binh, int &binpr, int &binn, int &bingst, TH1F *histoPt_p, TH1F *histoPt_l, TH1F *histoPt_h, TH1F *histoPt_pr,
    TH1F *histoPt_n, TH1F *histoPt_gst, TH1F *histoEta_p, TH1F *histoEta_l, TH1F *histoEta_h, TH1F *histoEta_pr,
    TH1F *histoEta_n, TH1F *histoEta_gst, TH1F *histoPhi_p, TH1F *histoPhi_l, TH1F *histoPhi_h, TH1F *histoPhi_pr,
    TH1F *histoPhi_n, TH1F *histoPhi_gst, TH1F *histomatchedstation_p, TH1F *histomatchedstation_l, TH1F *histomatchedstation_h,
    TH1F *histomatchedstation_pr, TH1F *histomatchedstation_n, TH1F *histomatchedstation_gst,
    TH1F *histoMuonsOOT_p, TH1F *histoMuonsOOT_l, TH1F *histoMuonsOOT_h, TH1F *histoMuonsOOT_pr, TH1F *histoMuonsOOT_n, 
    TH1F *histoMuonsOOT_gst, TH1F *histoMuonsPtOOT_p, TH1F *histoMuonsPtOOT_l, TH1F *histoMuonsPtOOT_h, TH1F *histoMuonsPtOOT_pr, 
    TH1F *histoMuonsPtOOT_n, TH1F *histoMuonsPtOOT_gst, TH1F *histoMuonsEtaOOT_p, TH1F *histoMuonsEtaOOT_l, TH1F *histoMuonsEtaOOT_h, 
    TH1F *histoMuonsEtaOOT_pr, TH1F *histoMuonsEtaOOT_n, TH1F *histoMuonsEtaOOT_gst, TH1F *histoMuonsValidHitsOOT_p,
    TH1F *histoMuonsValidHitsOOT_l, TH1F *histoMuonsValidHitsOOT_h, TH1F *histoMuonsValidHitsOOT_pr, TH1F *histoMuonsValidHitsOOT_n,
    TH1F *histoMuonsValidHitsOOT_gst, TH1F *histoMuonsValidGemHitsOOT_p, TH1F *histoMuonsValidGemHitsOOT_l, 
    TH1F *histoMuonsValidGemHitsOOT_h, TH1F *histoMuonsValidGemHitsOOT_pr, TH1F *histoMuonsValidGemHitsOOT_n, 
    TH1F *histoMuonsValidGemHitsOOT_gst, TH1F *histoMuonsInTime_p, TH1F *histoMuonsInTime_l,
    TH1F *histoMuonsInTime_h, TH1F *histoMuonsInTime_pr, TH1F *histoMuonsInTime_n, TH1F *histoMuonsInTime_gst, TH1F *histoMuonsPtInTime_p, 
    TH1F *histoMuonsPtInTime_l, TH1F *histoMuonsPtInTime_h, TH1F *histoMuonsPtInTime_pr, TH1F *histoMuonsPtInTime_n,
    TH1F *histoMuonsPtInTime_gst, TH1F *histoMuonsEtaInTime_p, TH1F *histoMuonsEtaInTime_l, TH1F *histoMuonsEtaInTime_h,
    TH1F *histoMuonsEtaInTime_pr, TH1F *histoMuonsEtaInTime_n, TH1F *histoMuonsEtaInTime_gst, TH1F *histoMuonsValidHitsInTime_p,
    TH1F *histoMuonsValidHitsInTime_l, TH1F *histoMuonsValidHitsInTime_h, TH1F *histoMuonsValidHitsInTime_pr, 
    TH1F *histoMuonsValidHitsInTime_n, TH1F *histoMuonsValidHitsInTime_gst, TH1F *histoMuonsValidGemHitsInTime_p,
    TH1F *histoMuonsValidGemHitsInTime_l, TH1F *histoMuonsValidGemHitsInTime_h, TH1F *histoMuonsValidGemHitsInTime_pr, 
    TH1F *histoMuonsValidGemHitsInTime_n, TH1F *histoMuonsValidGemHitsInTime_gst, int myCase)

{
  std::map<TrackingParticleRef, int> tpToSecondaries; // map from tp to (index+1) in output collection
  typedef MuonAssociatorByHits::MuonToSimCollection::const_iterator r2s_it;
  typedef MuonAssociatorByHits::SimToMuonCollection::const_iterator s2r_it;
  for (size_t i = 0; i < myNmu; ++i) //loop over reco::muons
  {

    edm::LogVerbatim("MuonMCClassifAndAna") << "\n reco::Muons # " << i;
    edm::RefToBase < reco::Muon > mu = myMuons_handle->refAt(i);
    if (std::find(selMuons.begin(), selMuons.end(), mu) == selMuons.end()) //if muon is not from a given type or ID - fail and go further, e.g. it might be tracker but not global and if we check for globals the value will be -99
    {
      edm::LogVerbatim("MuonMCClassifAndAna") << "\t muon didn't pass the selection. classified as -99 and skipped";
      if (debug)
        std::cout << "\t echo from myClassification" << std::endl;
      myClassif[i] = -99;
      continue;
    }
    TrackingParticleRef tp;
    edm::RefToBase < reco::Muon > muMatchBack;
    r2s_it match = myRecSimColl.find(mu);
    s2r_it matchback;

   if (match != myRecSimColl.end()) //look for match with sim particle
     {
        edm::LogVerbatim("MuonMCClassifAndAna") << "\t Reco to Sim matched Ok...";
   // match->second is vector, front is first element, first is the ref (second would be the quality)
    tp = match->second.front().first;
    myTpId[i] = tp.isNonnull() ? tp.key() : -1; // we check, even if null refs should not appear here at all; ako ref e nenuleva wzima Id na 4asticata
    myTpAssoQuality[i] = match->second.front().second; //quality
  
    //s2r_it matchback = mySimRecColl.find(tp);
      matchback = mySimRecColl.find(tp); //sim to reco match; t.e. obratno - tyrsi machnata tp v reco colekciata
           if (matchback != mySimRecColl.end())
               {
                  muMatchBack = matchback->second.front().first; //ako ima syvpadenie vzima referencia kym reco 4asticata
               }
          else
            {
              edm::LogWarning("MuonMCClassifAndAna") << "\n***WARNING: This I do NOT understand: why no match back? *** \n";
            }
         } //end reco to sim match OK

//new part:
//new part
// else if ((trackType_ == MuonAssociatorByHits::GlobalTk) && mu->isGlobalMuon())
// {
// perform a second attempt, matching with the standalone muon
// r2s_it matchSta = UpdSTA_recSimColl.find(mu);
// }
//new part end
    if (tp.isNonnull())
    {
      bool isGhost = muMatchBack != mu;
      if (isGhost)
        edm::LogVerbatim("MuonMCClassifAndAna") << "\t This seems a GHOST ! myClassif[i] will be < 0";
      myHitsPdgId[i] = tp->pdgId();
      myProdRho[i] = tp->vertex().Rho();
      myProdZ[i] = tp->vertex().Z();
      edm::LogVerbatim("MuonMCClassifAndAna") << "\t TP pdgId = " << myHitsPdgId[i] << ", vertex rho = " << myProdRho[i]
          << ", z = " << myProdZ[i];
// Try to extract mother and grand mother of this muon. Unfortunately, SIM and GEN histories require diffent code :-(
      if (!tp->genParticles().empty())
      {
// Muon is in GEN
        reco::GenParticleRef genp = tp->genParticles()[0];
        reco::GenParticleRef genMom = genp->numberOfMothers() > 0 ? genp->motherRef() : reco::GenParticleRef();
        if (genMom.isNonnull())
        {
          myMomPdgId[i] = genMom->pdgId();
          std::cout << "If muon is in GEN then pdg id of genMOM= " << genMom->pdgId() << std::endl;
          myMomStatus[i] = genMom->status();
          myMomRho[i] = genMom->vertex().Rho();
          myMomZ[i] = genMom->vz();
          edm::LogVerbatim("MuonMCClassifAndAna") << "\t Particle pdgId = " << myHitsPdgId[i] << " produced at rho = "
              << myProdRho[i] << ", z = " << myProdZ[i] << ", has GEN mother pdgId = " << myMomPdgId[i];
          reco::GenParticleRef genGMom = genMom->numberOfMothers() > 0 ? genMom->motherRef() : reco::GenParticleRef();
          if (genGMom.isNonnull())
          {
            myGmomPdgId[i] = genGMom->pdgId();
            edm::LogVerbatim("MuonMCClassifAndAna") << "\t\t mother prod. vertex rho = " << myMomRho[i] << ", z = "
                << myMomZ[i] << ", grand-mom pdgId = " << myGmomPdgId[i];
          }
// in this case, we might want to know the heaviest mom flavour
//loop stops when we're no longer looking at hadrons or mesons
          for (reco::GenParticleRef nMom = genMom; nMom.isNonnull() && abs(nMom->pdgId()) >= 100;
              nMom = nMom->numberOfMothers() > 0 ? nMom->motherRef() : reco::GenParticleRef())
          {
            int myFlav = flavour(nMom->pdgId());
            if (myHmomFlav[i] < myFlav)
              myHmomFlav[i] = myFlav;
            edm::LogVerbatim("MuonMCClassifAndAna") << "\t\t backtracking flavour: mom pdgId = " << nMom->pdgId()
                << ", flavour = " << myFlav << ", heaviest so far = " << myHmomFlav[i];
            std::cout << "from Gen: backtracking flavour: mom pdgId = " << nMom->pdgId() << ", flavour = " << myFlav
                << ", heaviest so far = " << myHmomFlav[i] << std::endl;
          }
        }
      }
      else
      {
// Muon is in SIM Only
        TrackingParticleRef simMom = getTpMother(tp);
        if (simMom.isNonnull())
        {
          myMomPdgId[i] = simMom->pdgId();
          std::cout << "If muon is SIM only then pdg id of genMOM= " << simMom->pdgId() << std::endl;
          myMomRho[i] = simMom->vertex().Rho();
          myMomZ[i] = simMom->vertex().Z();
          edm::LogVerbatim("MuonMCClassifAndAna") << "\t Particle pdgId = " << myHitsPdgId[i] << " produced at rho = "
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
            edm::LogVerbatim("MuonMCClassifAndAna") << "\t\t SIM mother is in GEN (status " << myMomStatus[i]
                << "), grand-mom id = " << myGmomPdgId[i];
          }
          else
          {
            myMomStatus[i] = -1;
            TrackingParticleRef simGMom = getTpMother(simMom);
            if (simGMom.isNonnull())
              myGmomPdgId[i] = simGMom->pdgId();
            edm::LogVerbatim("MuonMCClassifAndAna") << "\t\t SIM mother is in SIM only, grand-mom id = "
                << myGmomPdgId[i];
          }
        }
        else
        {
          edm::LogVerbatim("MuonMCClassifAndAna") << "\t Particle pdgId = " << myHitsPdgId[i] << " produced at rho = "
              << myProdRho[i] << ", z = " << myProdZ[i] << ", has NO mother!";
        }
      } // Muon is in SIM Only
      myMomFlav[i] = flavour(myMomPdgId[i]);
      myGmomFlav[i] = flavour(myGmomPdgId[i]);
// Check first IF this is a muon at all
      if (abs(tp->pdgId()) != 13)
      {
        myClassif[i] = isGhost ? -1 : 1;
        myExt[i] = isGhost ? -1 : 1;
        edm::LogVerbatim("MuonMCClassifAndAna") << "\t This is not a muon. Sorry. myClassif[i] = " << myClassif[i];
      }
      else
      {
// Is this SIM muon also a GEN muon, with a mother?
        if (!tp->genParticles().empty() && (myMomPdgId[i] != 0))
        {
// original if (abs(myMomPdgId[i]) < 100 && (abs(myMomPdgId[i]) != 15)) //15 is tau, it will be considered later with b and c mothers
//if ((abs(myMomPdgId[i]) < 100 || abs(myMomPdgId[i]) == 443) && (abs(myMomPdgId[i]) != 15)) //add J/psi to prompt mothers
          if ((abs(myMomPdgId[i]) < 100) && (abs(myMomPdgId[i]) != 15)) // original
          {
            myClassif[i] = isGhost ? -4 : 4;
            myFlav[i] = (abs(myMomPdgId[i]) == 15 ? 15 : 13);
            edm::LogVerbatim("MuonMCClassifAndAna") << "\t This seems PRIMARY MUON ! myClassif[i] = " << myClassif[i];
            if (debug)
              std::cout << "\t This seems PRIMARY MUON ! myClassif[i] = " << myClassif[i] << std::endl;
            myExt[i] = 10;
          }
          else if (myMomFlav[i] == 4 || myMomFlav[i] == 5 || myMomFlav[i] == 15) //if b, c or tau decays. e.g. heavy decaus
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
            edm::LogVerbatim("MuonMCClassifAndAna") << "\t This seems HEAVY FLAVOUR ! myClassif[i] = " << myClassif[i];
          }
          else
          {
            myClassif[i] = isGhost ? -2 : 2;
            myFlav[i] = myMomFlav[i];
            edm::LogVerbatim("MuonMCClassifAndAna") << "\t This seems LIGHT FLAVOUR ! myClassif[i] = " << myClassif[i];
          }
        } // Is this SIM muon also a GEN muon, with a mother?
        else
        {
          myClassif[i] = isGhost ? -2 : 2;
          myFlav[i] = myMomFlav[i];
          edm::LogVerbatim("MuonMCClassifAndAna") << "\t This seems LIGHT FLAVOUR ! myClassif[i] = " << myClassif[i];
          if (debug)
            std::cout << "\t This seems LIGHT FLAVOUR ! myClassif[i] = " << myClassif[i] << std::endl;
        }
// extended classification
        if (myMomPdgId[i] == 0)
          myExt[i] = 2; // if it has no mom, it's not a primary particle so it won't be in ppMuX
        else if (abs(myMomPdgId[i]) < 100)
          myExt[i] = (myMomFlav[i] == 15 ? 9 : 10); // primary mu (10), or tau->mu(9)
        else if (myMomFlav[i] == 5)
          myExt[i] = 8; // b->mu
        else if (myMomFlav[i] == 4)
          myExt[i] = (myHmomFlav[i] == 5 ? 7 : 6); // b->c->mu(b->c->...) and c->mu(6)
        else if (myMomStatus[i] != -1)
        {
// primary light particle
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
        edm::LogVerbatim("MuonMCClassifAndAna") << "\t Extended clasmyClassifion code = " << myExt[i];
//std::cout << "Extended clasmyClassifion code = " << myExt[i] << std::endl;
      } //else statement in place of continue
    } //ent tp is nonNull

    //case == 1 global; case == 2 tracker; case == 3 standalone; case == 4 good; case == 5 tight; case == 6 soft; case == 7 loose
    reco::TrackRef trackref;

    if (!(myCase == 1 || myCase == 2))  //all cases that are not global or standalone, e.g. tracker, good, tight, soft, loose
    {
      if ((*myMuons_handle)[i].muonBestTrack().isNonnull())
        trackref = (*myMuons_handle)[i].muonBestTrack();
      else
        continue;
    }
    else if(myCase == 2)  //standalone case
    {
      if ((*myMuons_handle)[i].outerTrack().isNonnull())
        trackref = (*myMuons_handle)[i].outerTrack();
      else
        continue;
    }
    else if(myCase == 1)  //global case
    {
      if ((*myMuons_handle)[i].globalTrack().isNonnull())
        trackref = (*myMuons_handle)[i].globalTrack();
      else
        continue;
    }


// if(myExt[i] == 1 || myExt[i] == 3)
    if (myClassif[i] == 1)
    {
      N_Muons->Fill(binp);
      N_Muons_norm->Fill(binp);
      histoPt_p->Fill(trackref->pt());
      histoEta_p->Fill(trackref->eta());
      histoPhi_p->Fill(trackref->phi());
      histomatchedstation_p->Fill((*myMuons_handle)[i].numberOfMatchedStations());

     if(fabs((*myMuons_handle)[i].time().timeAtIpInOut) > 25+fabs((*myMuons_handle)[i].time().timeAtIpInOutErr)) {
        histoMuonsOOT_p->Fill((*myMuons_handle)[i].time().timeAtIpInOut);
        histoMuonsPtOOT_p->Fill(trackref->pt());
        histoMuonsEtaOOT_p->Fill(trackref->eta());
        histoMuonsValidHitsOOT_p->Fill(trackref->hitPattern().numberOfValidMuonHits());
        histoMuonsValidGemHitsOOT_p->Fill(trackref->hitPattern().numberOfValidMuonGEMHits()); 
        }

    else{
        histoMuonsInTime_p->Fill((*myMuons_handle)[i].time().timeAtIpInOut);
        histoMuonsPtInTime_p->Fill(trackref->pt());
        histoMuonsEtaInTime_p->Fill(trackref->eta());
        histoMuonsValidHitsInTime_p->Fill(trackref->hitPattern().numberOfValidMuonHits());
        histoMuonsValidGemHitsInTime_p->Fill(trackref->hitPattern().numberOfValidMuonGEMHits());
        }

    }
// if(myExt[i] == 4 || myExt[i] == 5)
    if (myClassif[i] == 2)
    {
      N_Muons->Fill(binl);
      N_Muons_norm->Fill(binl);
      histoPt_l->Fill(trackref->pt());
      histoEta_l->Fill(trackref->eta());
      histoPhi_l->Fill(trackref->phi());
      histomatchedstation_l->Fill((*myMuons_handle)[i].numberOfMatchedStations());

      if(fabs((*myMuons_handle)[i].time().timeAtIpInOut) > 25+fabs((*myMuons_handle)[i].time().timeAtIpInOutErr)) {
        histoMuonsOOT_l->Fill((*myMuons_handle)[i].time().timeAtIpInOut);
        histoMuonsPtOOT_l->Fill(trackref->pt());
        histoMuonsEtaOOT_l->Fill(trackref->eta());
        histoMuonsValidHitsOOT_l->Fill(trackref->hitPattern().numberOfValidMuonHits());
        histoMuonsValidGemHitsOOT_l->Fill(trackref->hitPattern().numberOfValidMuonGEMHits());
        }

    else{
        histoMuonsInTime_l->Fill((*myMuons_handle)[i].time().timeAtIpInOut);
        histoMuonsPtInTime_l->Fill(trackref->pt());
        histoMuonsEtaInTime_l->Fill(trackref->eta());
        histoMuonsValidHitsInTime_l->Fill(trackref->hitPattern().numberOfValidMuonHits());
        histoMuonsValidGemHitsInTime_l->Fill(trackref->hitPattern().numberOfValidMuonGEMHits());

        }

    }
//if(myExt[i] >=6 && myExt[i] <=9)
    if (myClassif[i] == 3)
    {
      N_Muons->Fill(binh);
      N_Muons_norm->Fill(binh);
      histoPt_h->Fill(trackref->pt());
      histoEta_h->Fill(trackref->eta());
      histoPhi_h->Fill(trackref->phi());
      histomatchedstation_h->Fill((*myMuons_handle)[i].numberOfMatchedStations());

     if(fabs((*myMuons_handle)[i].time().timeAtIpInOut) > 25+fabs((*myMuons_handle)[i].time().timeAtIpInOutErr)) {
        histoMuonsOOT_h->Fill((*myMuons_handle)[i].time().timeAtIpInOut);
        histoMuonsPtOOT_h->Fill(trackref->pt());
        histoMuonsEtaOOT_h->Fill(trackref->eta());
        histoMuonsValidHitsOOT_h->Fill(trackref->hitPattern().numberOfValidMuonHits());
        histoMuonsValidGemHitsOOT_h->Fill(trackref->hitPattern().numberOfValidMuonGEMHits());

        }

    else{
        histoMuonsInTime_h->Fill((*myMuons_handle)[i].time().timeAtIpInOut);
        histoMuonsPtInTime_h->Fill(trackref->pt());
        histoMuonsEtaInTime_h->Fill(trackref->eta());
        histoMuonsValidHitsInTime_h->Fill(trackref->hitPattern().numberOfValidMuonHits());
        histoMuonsValidGemHitsInTime_h->Fill(trackref->hitPattern().numberOfValidMuonGEMHits());

        }

    }
// if(myExt[i] == 10)
    if (myClassif[i] == 4)
    {
      N_Muons->Fill(binpr);
      N_Muons_norm->Fill(binpr);
      histoPt_pr->Fill(trackref->pt());
      histoEta_pr->Fill(trackref->eta());
      histoPhi_pr->Fill(trackref->phi());
      histomatchedstation_pr->Fill((*myMuons_handle)[i].numberOfMatchedStations());

     if(fabs((*myMuons_handle)[i].time().timeAtIpInOut) > 25+fabs((*myMuons_handle)[i].time().timeAtIpInOutErr)) {
        histoMuonsOOT_pr->Fill((*myMuons_handle)[i].time().timeAtIpInOut);
        histoMuonsPtOOT_pr->Fill(trackref->pt());
        histoMuonsEtaOOT_pr->Fill(trackref->eta());
        histoMuonsValidHitsOOT_pr->Fill(trackref->hitPattern().numberOfValidMuonHits());
        histoMuonsValidGemHitsOOT_pr->Fill(trackref->hitPattern().numberOfValidMuonGEMHits());

        }

    else{
        histoMuonsInTime_pr->Fill((*myMuons_handle)[i].time().timeAtIpInOut);
        histoMuonsPtInTime_pr->Fill(trackref->pt());
        histoMuonsEtaInTime_pr->Fill(trackref->eta());
        histoMuonsValidHitsInTime_pr->Fill(trackref->hitPattern().numberOfValidMuonHits());
        histoMuonsValidGemHitsInTime_pr->Fill(trackref->hitPattern().numberOfValidMuonGEMHits());

        }

    }
//if (myExt[i] == 2)
    if (myClassif[i] == 0)
    {
      N_Muons->Fill(binn);
      N_Muons_norm->Fill(binn);
      histoPt_n->Fill(trackref->pt());
      histoEta_n->Fill(trackref->eta());
      histoPhi_n->Fill(trackref->phi());
      histomatchedstation_n->Fill((*myMuons_handle)[i].numberOfMatchedStations());
   
    if(fabs((*myMuons_handle)[i].time().timeAtIpInOut) > 25+fabs((*myMuons_handle)[i].time().timeAtIpInOutErr)) {
        histoMuonsOOT_n->Fill((*myMuons_handle)[i].time().timeAtIpInOut);
        histoMuonsPtOOT_n->Fill(trackref->pt());
        histoMuonsEtaOOT_n->Fill(trackref->eta());
        histoMuonsValidHitsOOT_n->Fill(trackref->hitPattern().numberOfValidMuonHits());
        histoMuonsValidGemHitsOOT_n->Fill(trackref->hitPattern().numberOfValidMuonGEMHits());

        }

    else{
        histoMuonsInTime_n->Fill((*myMuons_handle)[i].time().timeAtIpInOut);
        histoMuonsPtInTime_n->Fill(trackref->pt());
        histoMuonsEtaInTime_n->Fill(trackref->eta());
        histoMuonsValidHitsInTime_n->Fill(trackref->hitPattern().numberOfValidMuonHits());
        histoMuonsValidGemHitsInTime_n->Fill(trackref->hitPattern().numberOfValidMuonGEMHits());

        }


    }
//if(myExt[i] < 0 && fabs(myExt[i])<11)
    if (myClassif[i] < 0 && fabs(myClassif[i]) < 5)
    {
      N_Muons->Fill(bingst);
      N_Muons_norm->Fill(bingst);
      histoPt_gst->Fill(trackref->pt());
      histoEta_gst->Fill(trackref->eta());
      histoPhi_gst->Fill(trackref->phi());
      histomatchedstation_gst->Fill((*myMuons_handle)[i].numberOfMatchedStations());
    
     if(fabs((*myMuons_handle)[i].time().timeAtIpInOut) > 25+fabs((*myMuons_handle)[i].time().timeAtIpInOutErr)) {
        histoMuonsOOT_gst->Fill((*myMuons_handle)[i].time().timeAtIpInOut);
        histoMuonsPtOOT_gst->Fill(trackref->pt());
        histoMuonsEtaOOT_gst->Fill(trackref->eta());
        histoMuonsValidHitsOOT_gst->Fill(trackref->hitPattern().numberOfValidMuonHits());
        histoMuonsValidGemHitsOOT_gst->Fill(trackref->hitPattern().numberOfValidMuonGEMHits());

        }

    else{
        histoMuonsInTime_gst->Fill((*myMuons_handle)[i].time().timeAtIpInOut);
        histoMuonsPtInTime_gst->Fill(trackref->pt());
        histoMuonsEtaInTime_gst->Fill(trackref->eta());
        histoMuonsValidHitsInTime_gst->Fill(trackref->hitPattern().numberOfValidMuonHits());
        histoMuonsValidGemHitsInTime_gst->Fill(trackref->hitPattern().numberOfValidMuonGEMHits());

        }



    }
  } //ent loop over myMuons_handle size
 std::cout << "Event end !" << std::endl;
}
int MuonMCClassifAndAna::flavour(int pdgId) const
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
int MuonMCClassifAndAna::convertAndPush(const TrackingParticle &tp, reco::GenParticleCollection &out,
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
void MuonMCClassifAndAna::beginJob()
{
  edm::Service < TFileService > fs;
  N_Muons = fs->make < TH1F > ("N_muons", "Muon collections", 70, 0, 70);
  N_Muons->GetXaxis()->SetBinLabel(2, "All sel muons");
  N_Muons->GetXaxis()->SetBinLabel(4, "Glb muons");
  N_Muons->GetXaxis()->SetBinLabel(5, "Glb mu from punch");
  N_Muons->GetXaxis()->SetBinLabel(6, "Glb mu from light/decay");
  N_Muons->GetXaxis()->SetBinLabel(7, "Glb mu from heavy flv");
  N_Muons->GetXaxis()->SetBinLabel(8, "Glb mu from primary");
  N_Muons->GetXaxis()->SetBinLabel(9, "Glb mu from noise");
  N_Muons->GetXaxis()->SetBinLabel(10, "Glb mu from ghost");
  N_Muons->GetXaxis()->SetBinLabel(12, "Trk muons");
  N_Muons->GetXaxis()->SetBinLabel(13, "Trk mu from punch");
  N_Muons->GetXaxis()->SetBinLabel(14, "Trk mu from light/decay");
  N_Muons->GetXaxis()->SetBinLabel(15, "Trk mu from heavy flv");
  N_Muons->GetXaxis()->SetBinLabel(16, "Trk mu from primary");
  N_Muons->GetXaxis()->SetBinLabel(17, "Trk mu from noise");
  N_Muons->GetXaxis()->SetBinLabel(18, "Trk mu from ghost");
  N_Muons->GetXaxis()->SetBinLabel(20, "Sta muons");
  N_Muons->GetXaxis()->SetBinLabel(21, "Sta mu from punch");
  N_Muons->GetXaxis()->SetBinLabel(22, "Sta mu from light/decay");
  N_Muons->GetXaxis()->SetBinLabel(23, "Sta mu from heavy flv");
  N_Muons->GetXaxis()->SetBinLabel(24, "Sta mu from primary");
  N_Muons->GetXaxis()->SetBinLabel(25, "Sta mu from noise");
  N_Muons->GetXaxis()->SetBinLabel(26, "Sta mu from ghost");
  N_Muons->GetXaxis()->SetBinLabel(28, "Good muons");
  N_Muons->GetXaxis()->SetBinLabel(29, "Good mu from punch");
  N_Muons->GetXaxis()->SetBinLabel(30, "Good mu from light/decay");
  N_Muons->GetXaxis()->SetBinLabel(31, "Good mu from heavy flv");
  N_Muons->GetXaxis()->SetBinLabel(32, "Good mu from primary");
  N_Muons->GetXaxis()->SetBinLabel(33, "Good mu from noise");
  N_Muons->GetXaxis()->SetBinLabel(34, "Good mu from ghost");
  N_Muons->GetXaxis()->SetBinLabel(36, "Tight muons");
  N_Muons->GetXaxis()->SetBinLabel(37, "Tight mu from punch");
  N_Muons->GetXaxis()->SetBinLabel(38, "Tight mu from light/decay");
  N_Muons->GetXaxis()->SetBinLabel(39, "Tight mu from heavy flv");
  N_Muons->GetXaxis()->SetBinLabel(40, "Tight mu from primary");
  N_Muons->GetXaxis()->SetBinLabel(41, "Tight mu from noise");
  N_Muons->GetXaxis()->SetBinLabel(42, "Tight mu from ghost");
  N_Muons->GetXaxis()->SetBinLabel(44, "Soft muons");
  N_Muons->GetXaxis()->SetBinLabel(45, "Soft mu from punch");
  N_Muons->GetXaxis()->SetBinLabel(46, "Soft mu from light/decay");
  N_Muons->GetXaxis()->SetBinLabel(47, "Soft mu from heavy flv");
  N_Muons->GetXaxis()->SetBinLabel(48, "Soft mu from primary");
  N_Muons->GetXaxis()->SetBinLabel(49, "Soft mu from noise");
  N_Muons->GetXaxis()->SetBinLabel(50, "Soft mu from ghost");
  N_Muons->GetXaxis()->SetBinLabel(52, "Loose muons");
  N_Muons->GetXaxis()->SetBinLabel(53, "Loose mu from punch");
  N_Muons->GetXaxis()->SetBinLabel(54, "Loose mu from light/decay");
  N_Muons->GetXaxis()->SetBinLabel(55, "Loose mu from heavy flv");
  N_Muons->GetXaxis()->SetBinLabel(56, "Loose mu from primary");
  N_Muons->GetXaxis()->SetBinLabel(57, "Loose mu from noise");
  N_Muons->GetXaxis()->SetBinLabel(58, "Loose mu from ghost");
  N_Muons->GetXaxis()->SetBinLabel(61, "total muon evts");
  N_Muons->GetXaxis()->SetBinLabel(63, "number of PV");
  N_Muons->GetXaxis()->SetBinLabel(65, "total muons");
////////////////////
  N_Muons_norm = fs->make < TH1F > ("N_Muons_norm", "Muon collections Normalized to the #muons", 70, 0, 70);
  N_Muons_norm->GetXaxis()->SetBinLabel(2, "All sel muons");
  N_Muons_norm->GetXaxis()->SetBinLabel(4, "Glb muons");
  N_Muons_norm->GetXaxis()->SetBinLabel(5, "Glb mu from punch");
  N_Muons_norm->GetXaxis()->SetBinLabel(6, "Glb mu from light/decay");
  N_Muons_norm->GetXaxis()->SetBinLabel(7, "Glb mu from heavy flv");
  N_Muons_norm->GetXaxis()->SetBinLabel(8, "Glb mu from primary");
  N_Muons_norm->GetXaxis()->SetBinLabel(9, "Glb mu from noise");
  N_Muons_norm->GetXaxis()->SetBinLabel(10, "Glb mu from ghost");
  N_Muons_norm->GetXaxis()->SetBinLabel(12, "Trk muons");
  N_Muons_norm->GetXaxis()->SetBinLabel(13, "Trk mu from punch");
  N_Muons_norm->GetXaxis()->SetBinLabel(14, "Trk mu from light/decay");
  N_Muons_norm->GetXaxis()->SetBinLabel(15, "Trk mu from heavy flv");
  N_Muons_norm->GetXaxis()->SetBinLabel(16, "Trk mu from primary");
  N_Muons_norm->GetXaxis()->SetBinLabel(17, "Trk mu from noise");
  N_Muons_norm->GetXaxis()->SetBinLabel(18, "Trk mu from ghost");
  N_Muons_norm->GetXaxis()->SetBinLabel(20, "Sta muons");
  N_Muons_norm->GetXaxis()->SetBinLabel(21, "Sta mu from punch");
  N_Muons_norm->GetXaxis()->SetBinLabel(22, "Sta mu from light/decay");
  N_Muons_norm->GetXaxis()->SetBinLabel(23, "Sta mu from heavy flv");
  N_Muons_norm->GetXaxis()->SetBinLabel(24, "Sta mu from primary");
  N_Muons_norm->GetXaxis()->SetBinLabel(25, "Sta mu from noise");
  N_Muons_norm->GetXaxis()->SetBinLabel(26, "Sta mu from ghost");
  N_Muons_norm->GetXaxis()->SetBinLabel(28, "Good muons");
  N_Muons_norm->GetXaxis()->SetBinLabel(29, "Good mu from punch");
  N_Muons_norm->GetXaxis()->SetBinLabel(30, "Good mu from light/decay");
  N_Muons_norm->GetXaxis()->SetBinLabel(31, "Good mu from heavy flv");
  N_Muons_norm->GetXaxis()->SetBinLabel(32, "Good mu from primary");
  N_Muons_norm->GetXaxis()->SetBinLabel(33, "Good mu from noise");
  N_Muons_norm->GetXaxis()->SetBinLabel(34, "Good mu from ghost");
  N_Muons_norm->GetXaxis()->SetBinLabel(36, "Tight muons");
  N_Muons_norm->GetXaxis()->SetBinLabel(37, "Tight mu from punch");
  N_Muons_norm->GetXaxis()->SetBinLabel(38, "Tight mu from light/decay");
  N_Muons_norm->GetXaxis()->SetBinLabel(39, "Tight mu from heavy flv");
  N_Muons_norm->GetXaxis()->SetBinLabel(40, "Tight mu from primary");
  N_Muons_norm->GetXaxis()->SetBinLabel(41, "Tight mu from noise");
  N_Muons_norm->GetXaxis()->SetBinLabel(42, "Tight mu from ghost");
  N_Muons_norm->GetXaxis()->SetBinLabel(44, "Soft muons");
  N_Muons_norm->GetXaxis()->SetBinLabel(45, "Soft mu from punch");
  N_Muons_norm->GetXaxis()->SetBinLabel(46, "Soft mu from light/decay");
  N_Muons_norm->GetXaxis()->SetBinLabel(47, "Soft mu from heavy flv");
  N_Muons_norm->GetXaxis()->SetBinLabel(48, "Soft mu from primary");
  N_Muons_norm->GetXaxis()->SetBinLabel(49, "Soft mu from noise");
  N_Muons_norm->GetXaxis()->SetBinLabel(50, "Soft mu from ghost");
  N_Muons_norm->GetXaxis()->SetBinLabel(52, "Loose muons");
  N_Muons_norm->GetXaxis()->SetBinLabel(53, "Loose mu from punch");
  N_Muons_norm->GetXaxis()->SetBinLabel(54, "Loose mu from light/decay");
  N_Muons_norm->GetXaxis()->SetBinLabel(55, "Loose mu from heavy flv");
  N_Muons_norm->GetXaxis()->SetBinLabel(56, "Loose mu from primary");
  N_Muons_norm->GetXaxis()->SetBinLabel(57, "Loose mu from noise");
  N_Muons_norm->GetXaxis()->SetBinLabel(58, "Loose mu from ghost");
  N_Muons_norm->GetXaxis()->SetBinLabel(61, "total muon evts");
  N_Muons_norm->GetXaxis()->SetBinLabel(63, "number of PV");
  N_Muons_norm->GetXaxis()->SetBinLabel(65, "total muons");
  totalevents_withMuons = fs->make < TH1F > ("totalevents_withMuons", "Number events with Muons", 100, 0, 100); //
//all reco muons
  numb_MuPerEvent = fs->make < TH1F > ("numb_MuPerEvent", "Number of Reco Muons per Event", 500, 0, 500); //
  numb_PV = fs->make < TH1F > ("numb_PV", "Number of Reconstructed Primary Vertices", 100, 0, 100); //
  numb_SimPV = fs->make < TH1F > ("numb_SimPV", "Number of Simulated Primary Vertices", 100, 0, 100); //
  Pt_TrakingParticle = fs->make < TH1F > ("Pt_TrakingParticle", "p_{T} of Tracking Particles", 1000, 0., 500.);
  PtMuons_TrakingParticle = fs->make < TH1F
      > ("PtMuons_TrakingParticle", "p_{T} of muons from Tracking Particle collection", 1000, 0., 500.);
  EtaMuons_TrakingParticle = fs->make < TH1F
      > ("EtaMuons_TrakingParticle", "#eta of muons from Tracking Particle collection", 70, -3.5, 3.5);

//test reco
  allrecoMuonsEta = fs->make < TH1F > ("allRecoMuonsEta", "#eta RecoMuons", 70, -3.5, 3.5);
  allrecoMuonsPt = fs->make < TH1F > ("allRecoMuonsPt", "p_{T} RecoMuons", 1000, 0., 500.);
  allrecoMuonsPhi = fs->make < TH1F > ("allRecoMuonsPhi", "#varphi RecoMuons", 70, -3.5, 3.5);
  histoValidMuHits = fs->make < TH1F > ("histoValidMuHits", "histoValidMuHits", 50, 0, 50);
  histoValidMuGeMHits = fs->make < TH1F > ("histoValidMuGeMHits", "histoValidMuGeMHits", 50, 0, 50);
  histoMuStaWithValidHits = fs->make < TH1F > ("histoMuStaWithValidHits", "histoMuStaWithValidHits", 50, 0, 50);
  histoGemStaWithValidHits = fs->make < TH1F > ("histoGemStaWithValidHits", "histoGemStaWithValidHits", 50, 0, 50);
  histoValidMuHitsOut = fs->make < TH1F > ("histoValidMuHitsOut", "histoValidMuHitsOut", 50, 0, 50);
  histoValidMuGeMHitsOut = fs->make < TH1F > ("histoValidMuGeMHitsOut", "histoValidMuGeMHitsOut", 50, 0, 50);
  histoMuStaWithValidHitsOut = fs->make < TH1F
      > ("histoMuStaWithValidHitsOut", "histoMuStaWithValidHitsOut", 50, 0, 50);
  histoGemStaWithValidHitsOut = fs->make < TH1F
      > ("histoGemStaWithValidHitsOut", "histoGemStaWithValidHitsOut", 50, 0, 50);
  histoValidMuHitsTight = fs->make < TH1F > ("histoValidMuHitsTight", "histoValidMuHitsTight", 50, 0, 50);
  histoValidMuGeMHitsTight = fs->make < TH1F > ("histoValidMuGeMHitsTight", "histoValidMuGeMHitsTight", 50, 0, 50);
  histoMuStaWithValidHitsTight = fs->make < TH1F
      > ("histoMuStaWithValidHitsTight", "histoMuStaWithValidHitsTight", 50, 0, 50);
  histoGemStaWithValidHitsTight = fs->make < TH1F
      > ("histoGemStaWithValidHitsTight", "histoGemStaWithValidHitsTight", 50, 0, 50);
  histoValidMuHitsOutTight = fs->make < TH1F > ("histoValidMuHitsOutTight", "histoValidMuHitsOutTight", 50, 0, 50);
  histoValidMuGeMHitsOutTight = fs->make < TH1F
      > ("histoValidMuGeMHitsOutTight", "histoValidMuGeMHitsOutTight", 50, 0, 50);
  histoMuStaWithValidHitsOutTight = fs->make < TH1F
      > ("histoMuStaWithValidHitsOutTight", "histoMuStaWithValidHitsOutTight", 50, 0, 50);
  histoGemStaWithValidHitsOutTight = fs->make < TH1F
      > ("histoGemStaWithValidHitsOutTight", "histoGemStaWithValidHitsOutTight", 50, 0, 50);
  histoValidMuHitsSoft = fs->make < TH1F > ("histoValidMuHitsSoft", "histoValidMuHitsSoft", 50, 0, 50);
  histoValidMuGeMHitsSoft = fs->make < TH1F > ("histoValidMuGeMHitsSoft", "histoValidMuGeMHitsSoft", 50, 0, 50);
  histoMuStaWithValidHitsSoft = fs->make < TH1F
      > ("histoMuStaWithValidHitsSoft", "histoMuStaWithValidHitsSoft", 50, 0, 50);
  histoGemStaWithValidHitsSoft = fs->make < TH1F
      > ("histoGemStaWithValidHitsSoft", "histoGemStaWithValidHitsSoft", 50, 0, 50);
  histoValidMuHitsOutSoft = fs->make < TH1F > ("histoValidMuHitsOutSoft", "histoValidMuHitsOutSoft", 50, 0, 50);
  histoValidMuGeMHitsOutSoft = fs->make < TH1F
      > ("histoValidMuGeMHitsOutSoft", "histoValidMuGeMHitsOutSoft", 50, 0, 50);
  histoMuStaWithValidHitsOutSoft = fs->make < TH1F
      > ("histoMuStaWithValidHitsOutSoft", "histoMuStaWithValidHitsOutSoft", 50, 0, 50);
  histoGemStaWithValidHitsOutSoft = fs->make < TH1F
      > ("histoGemStaWithValidHitsOutSoft", "histoGemStaWithValidHitsOutSoft", 50, 0, 50);
  histoValidMuHitsLoose = fs->make < TH1F > ("histoValidMuHitsLoose", "histoValidMuHitsLoose", 50, 0, 50);
  histoValidMuGeMHitsLoose = fs->make < TH1F > ("histoValidMuGeMHitsLoose", "histoValidMuGeMHitsLoose", 50, 0, 50);
  histoMuStaWithValidHitsLoose = fs->make < TH1F
      > ("histoMuStaWithValidHitsLoose", "histoMuStaWithValidHitsLoose", 50, 0, 50);
  histoGemStaWithValidHitsLoose = fs->make < TH1F
      > ("histoGemStaWithValidHitsLoose", "histoGemStaWithValidHitsLoose", 50, 0, 50);
  histoValidMuHitsOutLoose = fs->make < TH1F > ("histoValidMuHitsOutLoose", "histoValidMuHitsOutLoose", 50, 0, 50);
  histoValidMuGeMHitsOutLoose = fs->make < TH1F
      > ("histoValidMuGeMHitsOutLoose", "histoValidMuGeMHitsOutLoose", 50, 0, 50);
  histoMuStaWithValidHitsOutLoose = fs->make < TH1F
      > ("histoMuStaWithValidHitsOutLoose", "histoMuStaWithValidHitsOutLoose", 50, 0, 50);
  histoGemStaWithValidHitsOutLoose = fs->make < TH1F
      > ("histoGemStaWithValidHitsOutLoose", "histoGemStaWithValidHitsOutLoose", 50, 0, 50);

///////////// for matched stations/////////////////////////////

  GlobalMuons_matchedstation = fs->make < TH1F > ("GlobalMuons_matchedstation", "Matched stations for global muon collection", 20, 0., 20.);
  StaMuons_matchedstation = fs->make < TH1F > ("StaMuons_matchedstation", "Matched stations for STA muon collection", 20, 0., 20.);
  TightMuons_matchedstation = fs->make < TH1F > ("TightMuons_matchedstation", "Matched stations for tight muon collection", 20, 0., 20.);
  LooseMuons_matchedstation = fs->make < TH1F > ("LooseMuons_matchedstation", "Matched stations for Loose muon collection", 20, 0., 20.);

/////////////// Plots for InTime and OOT /////////////////////

  StaMuonsOOT = fs->make < TH1F > ("StaMuonsOOT", "OOT for Standalone muons", 400, -200.,200);
  StaMuonsInTime = fs->make < TH1F > ("StaMuonsInTime", "In Time for Standalone muons", 400, -200.,200);
  StaMuonsPt_OOT = fs->make < TH1F > ("StaMuonsPt_OOT", "Pt for OOT Standalone muons", 1000, 0.,500);
  StaMuonsPt_InTime = fs->make < TH1F > ("StaMuonsPt_InTime", "Pt for in time Standalone muons", 1000, 0.,500);
  StaMuonsEta_OOT = fs->make < TH1F > ("StaMuonsEta_OOT", "Eta for OOT Standalone muons", 70, -3.5,3.5);
  StaMuonsEta_InTime = fs->make < TH1F > ("StaMuonsEta_InTime", "Eta for in time Standalone muons", 70, -3.5,3.5);

 StaMuonsValidHits_OOT = fs->make < TH1F > ("StaMuonsValidHits_OOT", "Valid Hits for OOT Standalone muons", 50, 0, 50);
 StaMuonsValidHits_InTime = fs->make < TH1F > ("StaMuonsValidHits_InTime", "Valid Hits for In Time Standalone muons", 50, 0, 50); 
 StaMuonsValidGemHits_OOT = fs->make < TH1F > ("StaMuonsValidGemHits_OOT", "Valid Gem Hits for OOT Standalone muons", 50, 0, 50);
 StaMuonsValidGemHits_InTime = fs->make < TH1F > ("StaMuonsValidGemHits_InTime", "Valid Gem Hits for In Time Standalone muons", 50, 0, 50);

///////// Plots for valid hits and Invalid hits /////////////////////////////

  StaMuons_validhits = fs->make < TH1F > ("StaMuons_validhits", "valid hits of STA muons", 50, 0, 50);
  StaMuons_validhits_atZerotime = fs->make < TH1F > ("StaMuons_validhits_atZerotime", "valid hits of STA muons at zero time", 50, 0, 50);
  StaMuonsPt_NoValidhits = fs->make < TH1F > ("StaMuonsPt_NoValidhits", "pt of STA muons with zero valid hits", 1000, 0., 500.);
  StaMuonsPt_NoValidhits_atZeroTime = fs->make < TH1F > ("StaMuonsPt_NoValidhits_atZeroTime", "pt of STA muons with zero valid hits at zero time", 1000, 0., 500.);

//test reco end
////////
  allMuonsPt = fs->make < TH1F > ("allMuonsPt", "p_{T} for All Selected Reco Muons", 1000, 0., 500.);
  allMuonsPt->SetXTitle("p_{T} [GeV]");
  allMuonsPt->SetYTitle("# of muons");
  allMuonsEta = fs->make < TH1F > ("allMuonsEta", "#eta for All Selected Reco Muons", 70, -3.5, 3.5);
  allMuonsEta->SetXTitle("#eta");
  allMuonsEta->SetYTitle("# of muons");
  allMuonsPhi = fs->make < TH1F > ("allMuonsPhi", "#varphi for All Selected Reco Muons", 70, -3.5, 3.5);
  allMuonsPhi->SetXTitle("#varphi");
  GlobalMuonsPt = fs->make < TH1F > ("GlobalMuonsPt", "p_{T} for Global Reco Muons", 1000, 0., 500.);
  GlobalMuonsEta = fs->make < TH1F > ("GlobalMuonsEta", "#eta for Global Reco Muons", 70, -3.5, 3.5);
  GlobalMuonsPhi = fs->make < TH1F > ("GlobalMuonsPhi", "#varphi for Global Reco Muons", 70, -3.5, 3.5);
  TrackerMuonsPt = fs->make < TH1F > ("TrackerMuonsPt", "p_{T} for Tracker Reco Muons", 1000, 0., 500.);
  TrackerMuonsEta = fs->make < TH1F > ("TrackerMuonsEta", "#eta for Tracker Reco Muons", 70, -3.5, 3.5);
  TrackerMuonsPhi = fs->make < TH1F > ("TrackerMuonsPhi", "#varphi for Tracker Reco Muons", 70, -3.5, 3.5);
  StaMuonsPt = fs->make < TH1F > ("StaMuonsPt", "p_{T} for Stand alone Reco Muons", 1000, 0., 500.);
  StaMuonsEta = fs->make < TH1F > ("StaMuonsEta", "#eta for Stand alone Reco Muons", 70, -3.5, 3.5);
  StaMuonsPhi = fs->make < TH1F > ("StaMuonsPhi", "#varphi for Stand alone Reco Muons", 70, -3.5, 3.5);
  GoodMuonsPt = fs->make < TH1F > ("GoodMuonsPt", "p_{T} for Good Reco Muons", 1000, 0., 500.);
  GoodMuonsEta = fs->make < TH1F > ("GoodMuonsEta", "#eta for Good Reco Muons", 70, -3.5, 3.5);
  GoodMuonsPhi = fs->make < TH1F > ("GoodMuonsPhi", "#varphi for Good Reco Muons", 70, -3.5, 3.5);
  TightMuonsPt = fs->make < TH1F > ("TightMuonsPt", "p_{T} for Tight Reco Muons", 1000, 0., 500.);
  TightMuonsEta = fs->make < TH1F > ("TightMuonsEta", "#eta for Tight Reco Muons", 70, -3.5, 3.5);
  TightMuonsPhi = fs->make < TH1F > ("TightMuonsPhi", "#varphi for Tight Reco Muons", 70, -3.5, 3.5);
  SoftMuonsPt = fs->make < TH1F > ("SoftMuonsPt", "p_{T} for Soft Reco Muons", 1000, 0., 500.);
  SoftMuonsEta = fs->make < TH1F > ("SoftMuonsEta", "#eta for Soft Reco Muons", 70, -3.5, 3.5);
  SoftMuonsPhi = fs->make < TH1F > ("SoftMuonsPhi", "#varphi for Soft Reco Muons", 70, -3.5, 3.5);
  LooseMuonsPt = fs->make < TH1F > ("LooseMuonsPt", "p_{T} for Loose Reco Muons", 1000, 0., 500.);
  LooseMuonsEta = fs->make < TH1F > ("LooseMuonsEta", "#eta for Loose Reco Muons", 70, -3.5, 3.5);
  LooseMuonsPhi = fs->make < TH1F > ("LooseMuonsPhi", "#varphi for Loose Reco Muons", 70, -3.5, 3.5);
  muPt_Glb_punch = fs->make < TH1F > ("muPt_Glb_punch", "p_{T} for Global punch thorugh Muons", 1000, 0., 500.);
  muPt_Glb_lightflv = fs->make < TH1F > ("muPt_Glb_lightflv", "p_{T} for Global light flavour Muons", 1000, 0., 500.);
  muPt_Glb_heavyflv = fs->make < TH1F > ("muPt_Glb_heavyflv", "p_{T} for Global heavy flavour Muons", 1000, 0., 500.);
  muPt_Glb_primary = fs->make < TH1F > ("muPt_Glb_primary", "p_{T} for Global primary Muons", 1000, 0., 500.);
  muPt_Glb_noise = fs->make < TH1F > ("muPt_Glb_noise", "p_{T} for Global noise", 1000, 0., 500.);
  muPt_Glb_ghost = fs->make < TH1F > ("muPt_Glb_ghost", "p_{T} for Global ghosts", 1000, 0., 500.);
  muPt_Trk_punch = fs->make < TH1F > ("muPt_Trk_punch", "p_{T} for Tracker punch thorugh Muons", 1000, 0., 500.);
  muPt_Trk_lightflv = fs->make < TH1F > ("muPt_Trk_lightflv", "p_{T} for Tracker light flavour Muons", 1000, 0., 500.);
  muPt_Trk_heavyflv = fs->make < TH1F > ("muPt_Trk_heavyflv", "p_{T} for Tracker heavy flavour Muons", 1000, 0., 500.);
  muPt_Trk_primary = fs->make < TH1F > ("muPt_Trk_primary", "p_{T} for Tracker primary Muons", 1000, 0., 500.);
  muPt_Trk_noise = fs->make < TH1F > ("muPt_Trk_noise", "p_{T} for Tracker noise", 1000, 0., 500.);
  muPt_Trk_ghost = fs->make < TH1F > ("muPt_Trk_ghost", "p_{T} for Tracker ghosts", 1000, 0., 500.);
  muPt_Sta_punch = fs->make < TH1F > ("muPt_Sta_punch", "p_{T} for Sta punch thorugh Muons", 1000, 0., 500.);
  muPt_Sta_lightflv = fs->make < TH1F > ("muPt_Sta_lightflv", "p_{T} for Sta light flavour Muons", 1000, 0., 500.);
  muPt_Sta_heavyflv = fs->make < TH1F > ("muPt_Sta_heavyflv", "p_{T} for Sta heavy flavour Muons", 1000, 0., 500.);
  muPt_Sta_primary = fs->make < TH1F > ("muPt_Sta_primary", "p_{T} for Sta primary Muons", 1000, 0., 500.);
  muPt_Sta_noise = fs->make < TH1F > ("muPt_Sta_noise", "p_{T} for Sta noise", 1000, 0., 500.);
  muPt_Sta_ghost = fs->make < TH1F > ("muPt_Sta_ghost", "p_{T} for Sta ghost", 1000, 0., 500.);
  muPt_Good_punch = fs->make < TH1F > ("muPt_Good_punch", "p_{T} for Good punch thorugh Muons", 1000, 0., 500.);
  muPt_Good_lightflv = fs->make < TH1F > ("muPt_Good_lightflv", "p_{T} for Good light flavour Muons", 1000, 0., 500.);
  muPt_Good_heavyflv = fs->make < TH1F > ("muPt_Good_heavyflv", "p_{T} for Good heavy flavour Muons", 1000, 0., 500.);
  muPt_Good_primary = fs->make < TH1F > ("muPt_Good_primary", "p_{T} for Good primary Muons", 1000, 0., 500.);
  muPt_Good_noise = fs->make < TH1F > ("muPt_Good_noise", "p_{T} for Good noise", 1000, 0., 500.);
  muPt_Good_ghost = fs->make < TH1F > ("muPt_Good_ghost", "p_{T} for Good ghosts", 1000, 0., 500.);
////
  muPt_Tight_punch = fs->make < TH1F > ("muPt_Tight_punch", "p_{T} for Tight punch thorugh Muons", 1000, 0., 500.);
  muPt_Tight_lightflv = fs->make < TH1F
      > ("muPt_Tight_lightflv", "p_{T} for Tight light flavour Muons", 1000, 0., 500.);
  muPt_Tight_heavyflv = fs->make < TH1F
      > ("muPt_Tight_heavyflv", "p_{T} for Tight heavy flavour Muons", 1000, 0., 500.);
  muPt_Tight_primary = fs->make < TH1F > ("muPt_Tight_primary", "p_{T} for Tight primary Muons", 1000, 0., 500.);
  muPt_Tight_noise = fs->make < TH1F > ("muPt_Tight_noise", "p_{T} for Tight noise", 1000, 0., 500.);
  muPt_Tight_ghost = fs->make < TH1F > ("muPt_Tight_ghost", "p_{T} for Tight ghosts", 1000, 0., 500.);
  muPt_Soft_punch = fs->make < TH1F > ("muPt_Soft_punch", "p_{T} for Soft punch thorugh Muons", 1000, 0., 500.);
  muPt_Soft_lightflv = fs->make < TH1F > ("muPt_Soft_lightflv", "p_{T} for Soft light flavour Muons", 1000, 0., 500.);
  muPt_Soft_heavyflv = fs->make < TH1F > ("muPt_Soft_heavyflv", "p_{T} for Soft heavy flavour Muons", 1000, 0., 500.);
  muPt_Soft_primary = fs->make < TH1F > ("muPt_Soft_primary", "p_{T} for Soft primary Muons", 1000, 0., 500.);
  muPt_Soft_noise = fs->make < TH1F > ("muPt_Soft_noise", "p_{T} for Soft noise", 1000, 0., 500.);
  muPt_Soft_ghost = fs->make < TH1F > ("muPt_Soft_ghost", "p_{T} for Soft ghosts", 1000, 0., 500.);
  muPt_Loose_punch = fs->make < TH1F > ("muPt_Loose_punch", "p_{T} for Loose punch thorugh Muons", 1000, 0., 500.);
  muPt_Loose_lightflv = fs->make < TH1F
      > ("muPt_Loose_lightflv", "p_{T} for Loose light flavour Muons", 1000, 0., 500.);
  muPt_Loose_heavyflv = fs->make < TH1F
      > ("muPt_Loose_heavyflv", "p_{T} for Loose heavy flavour Muons", 1000, 0., 500.);
  muPt_Loose_primary = fs->make < TH1F > ("muPt_Loose_primary", "p_{T} for Loose primary Muons", 1000, 0., 500.);
  muPt_Loose_noise = fs->make < TH1F > ("muPt_Loose_noise", "p_{T} for Loose noise", 1000, 0., 500.);
  muPt_Loose_ghost = fs->make < TH1F > ("muPt_Loose_ghost", "p_{T} for Loose ghosts", 1000, 0., 500.);
//eta
  muEta_Glb_punch = fs->make < TH1F > ("muEta_Glb_punch", "#eta for Global punch through Muons", 70, -3.5, 3.5);
  muEta_Glb_lightflv = fs->make < TH1F > ("muEta_Glb_lightflv", "#eta for Global light flavour Muons", 70, -3.5, 3.5);
  muEta_Glb_heavyflv = fs->make < TH1F > ("muEta_Glb_heavyflv", "#eta for Global heavy flavour Muons", 70, -3.5, 3.5);
  muEta_Glb_primary = fs->make < TH1F > ("muEta_Glb_primary", "#eta for Global primary Muons", 70, -3.5, 3.5);
  muEta_Glb_noise = fs->make < TH1F > ("muEta_Glb_noise", "#eta for Global noise", 70, -3.5, 3.5);
  muEta_Glb_ghost = fs->make < TH1F > ("muEta_Glb_ghost", "#eta for Global ghosts", 70, -3.5, 3.5);
  muEta_Trk_punch = fs->make < TH1F > ("muEta_Trk_punch", "#eta for Tracker punch through Muons", 70, -3.5, 3.5);
  muEta_Trk_lightflv = fs->make < TH1F > ("muEta_Trk_lightflv", "#eta for Tracker light flavour Muons", 70, -3.5, 3.5);
  muEta_Trk_heavyflv = fs->make < TH1F > ("muEta_Trk_heavyflv", "#eta for Tracker heavy flavour Muons", 70, -3.5, 3.5);
  muEta_Trk_primary = fs->make < TH1F > ("muEta_Trk_primary", "#eta for Tracker primary Muons", 70, -3.5, 3.5);
  muEta_Trk_noise = fs->make < TH1F > ("muEta_Trk_noise", "#eta for Tracker noise", 70, -3.5, 3.5);
  muEta_Trk_ghost = fs->make < TH1F > ("muEta_Trk_ghost", "#eta for Tracker ghosts", 70, -3.5, 3.5);
  muEta_Sta_punch = fs->make < TH1F > ("muEta_Sta_punch", "#eta for Sta punch through Muons", 70, -3.5, 3.5);
  muEta_Sta_lightflv = fs->make < TH1F > ("muEta_Sta_lightflv", "#eta for Sta light flavour Muons", 70, -3.5, 3.5);
  muEta_Sta_heavyflv = fs->make < TH1F > ("muEta_Sta_heavyflv", "#eta for Sta heavy flavour Muons", 70, -3.5, 3.5);
  muEta_Sta_primary = fs->make < TH1F > ("muEta_Sta_primary", "#eta for Sta primary Muons", 70, -3.5, 3.5);
  muEta_Sta_noise = fs->make < TH1F > ("muEta_Sta_noise", "#eta for Sta noise", 70, -3.5, 3.5);
  muEta_Sta_ghost = fs->make < TH1F > ("muEta_Sta_ghost", "#eta for Sta ghosts", 70, -3.5, 3.5);
  muEta_Good_punch = fs->make < TH1F > ("muEta_Good_punch", "#eta for Good punch through Muons", 70, -3.5, 3.5);
  muEta_Good_lightflv = fs->make < TH1F > ("muEta_Good_lightflv", "#eta for Good light flavour Muons", 70, -3.5, 3.5);
  muEta_Good_heavyflv = fs->make < TH1F > ("muEta_Good_heavyflv", "#eta for Good heavy flavour Muons", 70, -3.5, 3.5);
  muEta_Good_primary = fs->make < TH1F > ("muEta_Good_primary", "#eta for Good primary Muons", 70, -3.5, 3.5);
  muEta_Good_noise = fs->make < TH1F > ("muEta_Good_noise", "#eta for Good noise", 70, -3.5, 3.5);
  muEta_Good_ghost = fs->make < TH1F > ("muEta_Good_ghost", "#eta for Good ghosts", 70, -3.5, 3.5);
  muEta_Tight_punch = fs->make < TH1F > ("muEta_Tight_punch", "#eta for Tight punch through Muons", 70, -3.5, 3.5);
  muEta_Tight_lightflv = fs->make < TH1F
      > ("muEta_Tight_lightflv", "#eta for Tight light flavour Muons", 70, -3.5, 3.5);
  muEta_Tight_heavyflv = fs->make < TH1F
      > ("muEta_Tight_heavyflv", "#eta for Tight heavy flavour Muons", 70, -3.5, 3.5);
  muEta_Tight_primary = fs->make < TH1F > ("muEta_Tight_primary", "#eta for Tight primary Muons", 70, -3.5, 3.5);
  muEta_Tight_noise = fs->make < TH1F > ("muEta_Tight_noise", "#eta for Tight noise", 70, -3.5, 3.5);
  muEta_Tight_ghost = fs->make < TH1F > ("muEta_Tight_ghost", "#eta for Tight ghosts", 70, -3.5, 3.5);
  muEta_Soft_punch = fs->make < TH1F > ("muEta_Soft_punch", "#eta for Soft punch through Muons", 70, -3.5, 3.5);
  muEta_Soft_lightflv = fs->make < TH1F > ("muEta_Soft_lightflv", "#eta for Soft light flavour Muons", 70, -3.5, 3.5);
  muEta_Soft_heavyflv = fs->make < TH1F > ("muEta_Soft_heavyflv", "#eta for Soft heavy flavour Muons", 70, -3.5, 3.5);
  muEta_Soft_primary = fs->make < TH1F > ("muEta_Soft_primary", "#eta for Soft primary Muons", 70, -3.5, 3.5);
  muEta_Soft_noise = fs->make < TH1F > ("muEta_Soft_noise", "#eta for Soft noise", 70, -3.5, 3.5);
  muEta_Soft_ghost = fs->make < TH1F > ("muEta_Soft_ghost", "#eta for Soft ghosts", 70, -3.5, 3.5);
  muEta_Loose_punch = fs->make < TH1F > ("muEta_Loose_punch", "#eta for Loose punch through Muons", 70, -3.5, 3.5);
  muEta_Loose_lightflv = fs->make < TH1F
      > ("muEta_Loose_lightflv", "#eta for Loose light flavour Muons", 70, -3.5, 3.5);
  muEta_Loose_heavyflv = fs->make < TH1F
      > ("muEta_Loose_heavyflv", "#eta for Loose heavy flavour Muons", 70, -3.5, 3.5);
  muEta_Loose_primary = fs->make < TH1F > ("muEta_Loose_primary", "#eta for Loose primary Muons", 70, -3.5, 3.5);
  muEta_Loose_noise = fs->make < TH1F > ("muEta_Loose_noise", "#eta for Loose noise", 70, -3.5, 3.5);
  muEta_Loose_ghost = fs->make < TH1F > ("muEta_Loose_ghost", "#eta for Loose ghosts", 70, -3.5, 3.5);
//phi
  muPhi_Glb_punch = fs->make < TH1F > ("muPhi_Glb_punch", "#varphi for Global punch through Muons", 70, -3.5, 3.5);
  muPhi_Glb_lightflv = fs->make < TH1F
      > ("muPhi_Glb_lightflv", "#varphi for Global light flavour Muons", 70, -3.5, 3.5); //
  muPhi_Glb_heavyflv = fs->make < TH1F
      > ("muPhi_Glb_heavyflv", "#varphi for Global heavy flavour Muons", 70, -3.5, 3.5); //
  muPhi_Glb_primary = fs->make < TH1F > ("muPhi_Glb_primary", "#varphi for Global primary Muons", 70, -3.5, 3.5);
  muPhi_Glb_noise = fs->make < TH1F > ("muPhi_Glb_noise", "#varphi for Global noise", 70, -3.5, 3.5);
  muPhi_Glb_ghost = fs->make < TH1F > ("muPhi_Glb_ghost", "#varphi for Global ghosts", 70, -3.5, 3.5);
  muPhi_Trk_punch = fs->make < TH1F > ("muPhi_Trk_punch", "#varphi for Tracker punch through Muons", 70, -3.5, 3.5);
  muPhi_Trk_lightflv = fs->make < TH1F
      > ("muPhi_Trk_lightflv", "#varphi for Tracker light flavour Muons", 70, -3.5, 3.5); //
  muPhi_Trk_heavyflv = fs->make < TH1F
      > ("muPhi_Trk_heavyflv", "#varphi for Tracker heavy flavour Muons", 70, -3.5, 3.5); //
  muPhi_Trk_primary = fs->make < TH1F > ("muPhi_Trk_primary", "#varphi for Tracker primary Muons", 70, -3.5, 3.5);
  muPhi_Trk_noise = fs->make < TH1F > ("muPhi_Trk_noise", "#varphi for Tracker noise", 70, -3.5, 3.5);
  muPhi_Trk_ghost = fs->make < TH1F > ("muPhi_Trk_ghost", "#varphi for Tracker ghosts", 70, -3.5, 3.5);
  muPhi_Sta_punch = fs->make < TH1F > ("muPhi_Sta_punch", "#varphi for Sta punch through Muons", 70, -3.5, 3.5);
  muPhi_Sta_lightflv = fs->make < TH1F > ("muPhi_Sta_lightflv", "#varphi for Sta light flavour Muons", 70, -3.5, 3.5);
  muPhi_Sta_heavyflv = fs->make < TH1F > ("muPhi_Sta_heavyflv", "#varphi for Sta heavy flavour Muons", 70, -3.5, 3.5);
  muPhi_Sta_primary = fs->make < TH1F > ("muPhi_Sta_primary", "#varphi for Sta primary Muons", 70, -3.5, 3.5);
  muPhi_Sta_noise = fs->make < TH1F > ("muPhi_Sta_noise", "#varphi for Sta noise", 70, -3.5, 3.5);
  muPhi_Sta_ghost = fs->make < TH1F > ("muPhi_Sta_ghost", "#varphi for Sta ghosts", 70, -3.5, 3.5);
  muPhi_Good_punch = fs->make < TH1F > ("muPhi_Good_punch", "#varphi for Good punch through Muons", 70, -3.5, 3.5);
  muPhi_Good_lightflv = fs->make < TH1F
      > ("muPhi_Good_lightflv", "#varphi for Good light flavour Muons", 70, -3.5, 3.5); //
  muPhi_Good_heavyflv = fs->make < TH1F
      > ("muPhi_Good_heavyflv", "#varphi for Good heavy flavour Muons", 70, -3.5, 3.5); //
  muPhi_Good_primary = fs->make < TH1F > ("muPhi_Good_primary", "#varphi for Good primary Muons", 70, -3.5, 3.5);
  muPhi_Good_noise = fs->make < TH1F > ("muPhi_Good_noise", "#varphi for Good noise", 70, -3.5, 3.5);
  muPhi_Good_ghost = fs->make < TH1F > ("muPhi_Good_ghost", "#varphi for Good ghosts", 70, -3.5, 3.5);
////
  muPhi_Tight_punch = fs->make < TH1F > ("muPhi_Tight_punch", "#varphi for Tight punch through Muons", 70, -3.5, 3.5);
  muPhi_Tight_lightflv = fs->make < TH1F
      > ("muPhi_Tight_lightflv", "#varphi for Tight light flavour Muons", 70, -3.5, 3.5); //
  muPhi_Tight_heavyflv = fs->make < TH1F
      > ("muPhi_Tight_heavyflv", "#varphi for Tight heavy flavour Muons", 70, -3.5, 3.5); //
  muPhi_Tight_primary = fs->make < TH1F > ("muPhi_Tight_primary", "#varphi for Tight primary Muons", 70, -3.5, 3.5);
  muPhi_Tight_noise = fs->make < TH1F > ("muPhi_Tight_noise", "#varphi for Tight noise", 70, -3.5, 3.5);
  muPhi_Tight_ghost = fs->make < TH1F > ("muPhi_Tight_ghost", "#varphi for Tight ghosts", 70, -3.5, 3.5);
  muPhi_Soft_punch = fs->make < TH1F > ("muPhi_Soft_punch", "#varphi for Soft punch through Muons", 70, -3.5, 3.5);
  muPhi_Soft_lightflv = fs->make < TH1F
      > ("muPhi_Soft_lightflv", "#varphi for Soft light flavour Muons", 70, -3.5, 3.5); //
  muPhi_Soft_heavyflv = fs->make < TH1F
      > ("muPhi_Soft_heavyflv", "#varphi for Soft heavy flavour Muons", 70, -3.5, 3.5); //
  muPhi_Soft_primary = fs->make < TH1F > ("muPhi_Soft_primary", "#varphi for Soft primary Muons", 70, -3.5, 3.5);
  muPhi_Soft_noise = fs->make < TH1F > ("muPhi_Soft_noise", "#varphi for Soft noise", 70, -3.5, 3.5);
  muPhi_Soft_ghost = fs->make < TH1F > ("muPhi_Soft_ghost", "#varphi for Soft ghosts", 70, -3.5, 3.5);
  muPhi_Loose_punch = fs->make < TH1F > ("muPhi_Loose_punch", "#varphi for Loose punch through Muons", 70, -3.5, 3.5);
  muPhi_Loose_lightflv = fs->make < TH1F
      > ("muPhi_Loose_lightflv", "#varphi for Loose light flavour Muons", 70, -3.5, 3.5); //
  muPhi_Loose_heavyflv = fs->make < TH1F
      > ("muPhi_Loose_heavyflv", "#varphi for Loose heavy flavour Muons", 70, -3.5, 3.5); //
  muPhi_Loose_primary = fs->make < TH1F > ("muPhi_Loose_primary", "#varphi for Loose primary Muons", 70, -3.5, 3.5);
  muPhi_Loose_noise = fs->make < TH1F > ("muPhi_Loose_noise", "#varphi for Loose noise", 70, -3.5, 3.5);
  muPhi_Loose_ghost = fs->make < TH1F > ("muPhi_Loose_ghost", "#varphi for Loose ghosts", 70, -3.5, 3.5);
////matched stations
  muMStation_Glb_punch = fs->make < TH1F > ("muMStation_Glb_punch", "Matched stations for Global punch through Muons", 20, 0., 20.);
  muMStation_Glb_lightflv = fs->make < TH1F
      > ("muMStation_Glb_lightflv", "Matched stations for Global light flavour Muons", 20, 0., 20.); //
  muMStation_Glb_heavyflv = fs->make < TH1F
      > ("muMStation_Glb_heavyflv", "Matched stations for Global heavy flavour Muons", 20, 0., 20.); //
  muMStation_Glb_primary = fs->make < TH1F > ("muMStation_Glb_primary", "Matched stations for Global primary Muons", 20, 0., 20.);
  muMStation_Glb_noise = fs->make < TH1F > ("muMStation_Glb_noise", "Matched stations for Global noise", 20, 0., 20.);
  muMStation_Glb_ghost = fs->make < TH1F > ("muMStation_Glb_ghost", "Matched stations for Global ghosts", 20, 0., 20.);
  muMStation_Sta_punch = fs->make < TH1F > ("muMStation_Sta_punch", "Matched stations for Sta punch through Muons", 20, 0., 20.);
  muMStation_Sta_lightflv = fs->make < TH1F
      > ("muMStation_Sta_lightflv", "Matched stations for Sta light flavour Muons", 20, 0., 20.); //
  muMStation_Sta_heavyflv = fs->make < TH1F
      > ("muMStation_Sta_heavyflv", "Matched stations for Sta heavy flavour Muons", 20, 0., 20.); //
  muMStation_Sta_primary = fs->make < TH1F > ("muMStation_Sta_primary", "Matched stations for Sta primary Muons", 20, 0., 20.);
  muMStation_Sta_noise = fs->make < TH1F > ("muMStation_Sta_noise", "Matched stations for Sta noise", 20, 0., 20.);
  muMStation_Sta_ghost = fs->make < TH1F > ("muMStation_Sta_ghost", "Matched stations for Sta ghosts", 20, 0., 20.);
  muMStation_Tight_punch = fs->make < TH1F > ("muMStation_Tight_punch", "Matched stations for Tight punch through Muons", 20, 0., 20.);
  muMStation_Tight_lightflv = fs->make < TH1F
      > ("muMStation_Tight_lightflv", "Matched stations for Tight light flavour Muons", 20, 0., 20.); //
  muMStation_Tight_heavyflv = fs->make < TH1F
      > ("muMStation_Tight_heavyflv", "Matched stations for Tight heavy flavour Muons", 20, 0., 20.); //
  muMStation_Tight_primary = fs->make < TH1F > ("muMStation_Tight_primary", "Matched stations for Tight primary Muons", 20, 0., 20.);
  muMStation_Tight_noise = fs->make < TH1F > ("muMStation_Tight_noise", "Matched stations for Tight noise", 20, 0., 20.);
  muMStation_Tight_ghost = fs->make < TH1F > ("muMStation_Tight_ghost", "Matched stations for Tight ghosts", 20, 0., 20.);
  muMStation_Loose_punch = fs->make < TH1F > ("muMStation_Loose_punch", "Matched stations for Loose punch through Muons", 20, 0., 20.);
  muMStation_Loose_lightflv = fs->make < TH1F
      > ("muMStation_Loose_lightflv", "Matched stations for Loose light flavour Muons", 20, 0., 20.); //
  muMStation_Loose_heavyflv = fs->make < TH1F
      > ("muMStation_Loose_heavyflv", "Matched stations for Loose heavy flavour Muons", 20, 0., 20.); //
  muMStation_Loose_primary = fs->make < TH1F > ("muMStation_Loose_primary", "Matched stations for Loose primary Muons", 20, 0., 20.);
  muMStation_Loose_noise = fs->make < TH1F > ("muMStation_Loose_noise", "Matched stations for Loose noise", 20, 0., 20.);
  muMStation_Loose_ghost = fs->make < TH1F > ("muMStation_Loose_ghost", "Matched stations for Loose ghosts", 20, 0., 20.);

  muOOT_Glb_punch = fs->make < TH1F > ("muOOT_Glb_punch", "OOT for global muons from punch through", 400, -200.,200);
  muOOT_Glb_lightflv = fs->make < TH1F > ("muOOT_Glb_lightflv", "OOT for global muons from light flavour muons", 400, -200.,200);
  muOOT_Glb_heavyflv = fs->make < TH1F > ("muOOT_Glb_heavyflv", "OOT for global muons from heavy flavour muons", 400, -200.,200);
  muOOT_Glb_primary = fs->make < TH1F > ("muOOT_Glb_primary", "OOT for global muons from primary", 400, -200.,200);
  muOOT_Glb_noise = fs->make < TH1F > ("muOOT_Glb_noise", "OOT for global muons from noise", 400, -200.,200);
  muOOT_Glb_ghost = fs->make < TH1F > ("muOOT_Glb_ghost", "OOT for global muons from ghosts", 400, -200.,200);

  muPtOOT_Glb_punch = fs->make < TH1F > ("muPtOOT_Glb_punch", "Pt of OOT for global muons from punch through", 1000, 0., 500.);
  muPtOOT_Glb_lightflv = fs->make < TH1F > ("muPtOOT_Glb_lightflv", "Pt of OOT for global muons from light flavour muons", 1000, 0., 500.);
  muPtOOT_Glb_heavyflv = fs->make < TH1F > ("muPtOOT_Glb_heavyflv", "Pt of OOT for global muons from heavy flavour muons", 1000, 0., 500.);
  muPtOOT_Glb_primary = fs->make < TH1F > ("muPtOOT_Glb_primary", "Pt of OOT for global muons from primary", 1000, 0., 500.);
  muPtOOT_Glb_noise = fs->make < TH1F > ("muPtOOT_Glb_noise", "Pt of OOT for global muons from noise", 1000, 0., 500.);
  muPtOOT_Glb_ghost = fs->make < TH1F > ("muPtOOT_Glb_ghost", "Pt of OOT for global muons from ghosts", 1000, 0., 500.);

  muEtaOOT_Glb_punch = fs->make < TH1F > ("muEtaOOT_Glb_punch", "Eta of OOT for global muons from punch through", 70, -3.5, 3.5);
  muEtaOOT_Glb_lightflv = fs->make < TH1F > ("muEtaOOT_Glb_lightflv", "Eta of OOT for global muons from light flavour muons", 70, -3.5, 3.5);
  muEtaOOT_Glb_heavyflv = fs->make < TH1F > ("muEtaOOT_Glb_heavyflv", "Eta of OOT for global muons from heavy flavour muons", 70, -3.5, 3.5);
  muEtaOOT_Glb_primary = fs->make < TH1F > ("muEtaOOT_Glb_primary", "Eta of OOT for global muons from primary", 70, -3.5, 3.5);
  muEtaOOT_Glb_noise = fs->make < TH1F > ("muEtaOOT_Glb_noise", "Eta of OOT for global muons from noise", 70, -3.5, 3.5);
  muEtaOOT_Glb_ghost = fs->make < TH1F > ("muEtaOOT_Glb_ghost", "Eta of OOT for global muons from ghosts", 70, -3.5, 3.5);

  muValidHitsOOT_Glb_punch = fs->make < TH1F > ("muValidHitsOOT_Glb_punch", "Valid Hits of OOT for Glb  muons from punch through", 50, 0, 50);
  muValidHitsOOT_Glb_lightflv = fs->make < TH1F > ("muValidHitsOOT_Glb_lightflv", "Valid Hits of OOT for Glb  muons from light flavour muons", 50, 0, 50 );
  muValidHitsOOT_Glb_heavyflv = fs->make < TH1F > ("muValidHitsOOT_Glb_heavyflv", "Valid Hits of OOT for Glb  muons from heavy flavour muons", 50, 0, 50);
  muValidHitsOOT_Glb_primary = fs->make < TH1F > ("muValidHitsOOT_Glb_primary", "Valid Hits of OOT for Glb  muons from primary", 50, 0, 50);
  muValidHitsOOT_Glb_noise = fs->make < TH1F > ("muValidHitsOOT_Glb_noise", "Valid Hits of OOT for Glb  muons from noise", 50, 0, 50);
  muValidHitsOOT_Glb_ghost = fs->make < TH1F > ("muValidHitsOOT_Glb_ghost", "Valid Hits of OOT for Glb  muons from ghosts", 50, 0, 50);

  muValidGemHitsOOT_Glb_punch = fs->make < TH1F > ("muValidGemHitsOOT_Glb_punch", "Valid Gem Hits of OOT for Glb  muons from punch through", 50, 0, 50);
  muValidGemHitsOOT_Glb_lightflv = fs->make < TH1F > ("muValidGemHitsOOT_Glb_lightflv", "Valid Gem Hits of OOT for Glb  muons from light flavour muons", 50, 0, 50 );
  muValidGemHitsOOT_Glb_heavyflv = fs->make < TH1F > ("muValidGemHitsOOT_Glb_heavyflv", "Valid Gem Hits of OOT for Glb  muons from heavy flavour muons", 50, 0, 50);
  muValidGemHitsOOT_Glb_primary = fs->make < TH1F > ("muValidGemHitsOOT_Glb_primary", "Valid Gem Hits of OOT for Glb  muons from primary", 50, 0, 50);
  muValidGemHitsOOT_Glb_noise = fs->make < TH1F > ("muValidGemHitsOOT_Glb_noise", "Valid Gem Hits of OOT for Glb  muons from noise", 50, 0, 50);
  muValidGemHitsOOT_Glb_ghost = fs->make < TH1F > ("muValidGemHitsOOT_Glb_ghost", "Valid Gem Hits of OOT for Glb  muons from ghosts", 50, 0, 50);
 
  muInTime_Glb_punch = fs->make < TH1F > ("muInTime_Glb_punch", "In Time for global muons from punch through", 400, -200.,200);
  muInTime_Glb_lightflv = fs->make < TH1F > ("muInTime_Glb_lightflv", "In Time for global muons from light flavour muons", 400, -200.,200);
  muInTime_Glb_heavyflv = fs->make < TH1F > ("muInTime_Glb_heavyflv", "In Time for global muons from heavy flavour muons", 400, -200.,200);
  muInTime_Glb_primary = fs->make < TH1F > ("muInTime_Glb_primary", "In Time for global muons from primary", 400, -200.,200);
  muInTime_Glb_noise = fs->make < TH1F > ("muInTime_Glb_noise", "In Time for global muons from noise", 400, -200.,200);
  muInTime_Glb_ghost = fs->make < TH1F > ("muInTime_Glb_ghost", "In Time for global muons from ghosts", 400, -200.,200);

  muPtInTime_Glb_punch = fs->make < TH1F > ("muPtInTime_Glb_punch", "Pt of In Time for global muons from punch through", 1000, 0., 500.);
  muPtInTime_Glb_lightflv = fs->make < TH1F > ("muPtInTime_Glb_lightflv", "Pt of In Time for global muons from light flavour muons", 1000, 0., 500.);
  muPtInTime_Glb_heavyflv = fs->make < TH1F > ("muPtInTime_Glb_heavyflv", "Pt of In Time for global muons from heavy flavour muons", 1000, 0., 500.);
  muPtInTime_Glb_primary = fs->make < TH1F > ("muPtInTime_Glb_primary", "Pt of In Time for global muons from primary", 1000, 0., 500.);
  muPtInTime_Glb_noise = fs->make < TH1F > ("muPtInTime_Glb_noise", "Pt of In Time for global muons from noise", 1000, 0., 500.);
  muPtInTime_Glb_ghost = fs->make < TH1F > ("muPtInTime_Glb_ghost", "Pt of In Time for global muons from ghosts", 1000, 0., 500.);

  muEtaInTime_Glb_punch = fs->make < TH1F > ("muEtaInTime_Glb_punch", "Eta of In Time for global muons from punch through", 70, -3.5, 3.5);
  muEtaInTime_Glb_lightflv = fs->make < TH1F > ("muEtaInTime_Glb_lightflv", "Eta of In Time for global muons from light flavour muons", 70, -3.5, 3.5);
  muEtaInTime_Glb_heavyflv = fs->make < TH1F > ("muEtaInTime_Glb_heavyflv", "Eta of In Time for global muons from heavy flavour muons", 70, -3.5, 3.5);
  muEtaInTime_Glb_primary = fs->make < TH1F > ("muEtaInTme_Glb_primary", "Eta of In Time for global muons from primary", 70, -3.5, 3.5);
  muEtaInTime_Glb_noise = fs->make < TH1F > ("muEtaInTime_Glb_noise", "Eta of In Time for global muons from noise", 70, -3.5, 3.5);
  muEtaInTime_Glb_ghost = fs->make < TH1F > ("muEtaInTime_Glb_ghost", "Eta of In Time for global muons from ghosts", 70, -3.5, 3.5);

  muValidHitsInTime_Glb_punch = fs->make < TH1F > ("muValidHitsInTime_Glb_punch", "Valid Hits of In Time for Glb  muons from punch through", 50, 0, 50);
  muValidHitsInTime_Glb_lightflv = fs->make < TH1F > ("muValidHitsInTime_Glb_lightflv", "Valid Hits of In Time for Glb  muons from light flavour muons", 50, 0, 50 );
  muValidHitsInTime_Glb_heavyflv = fs->make < TH1F > ("muValidHitsInTime_Glb_heavyflv", "Valid Hits of In Time for Glb  muons from heavy flavour muons", 50, 0, 50);
  muValidHitsInTime_Glb_primary = fs->make < TH1F > ("muValidHitsInTime_Glb_primary", "Valid Hits of In Time for Glb  muons from primary", 50, 0, 50);
  muValidHitsInTime_Glb_noise = fs->make < TH1F > ("muValidHitsInTime_Glb_noise", "Valid Hits of In Time for Glb  muons from noise", 50, 0, 50);
  muValidHitsInTime_Glb_ghost = fs->make < TH1F > ("muValidHitsInTime_Glb_ghost", "Valid Hits of In Time for Glb  muons from ghosts", 50, 0, 50);

  muValidGemHitsInTime_Glb_punch = fs->make < TH1F > ("muValidGemHitsInTime_Glb_punch", "Valid Gem Hits of In Time for Glb  muons from punch through", 50, 0, 50);
  muValidGemHitsInTime_Glb_lightflv = fs->make < TH1F > ("muValidGemHitsInTime_Glb_lightflv", "Valid Gem Hits of In Time for Glb  muons from light flavour muons", 50, 0, 50 );
  muValidGemHitsInTime_Glb_heavyflv = fs->make < TH1F > ("muValidGemHitsInTime_Glb_heavyflv", "Valid Gem Hits of In Time for Glb  muons from heavy flavour muons", 50, 0, 50);
  muValidGemHitsInTime_Glb_primary = fs->make < TH1F > ("muValidGemHitsInTime_Glb_primary", "Valid Gem Hits of In Time for Glb  muons from primary", 50, 0, 50);
  muValidGemHitsInTime_Glb_noise = fs->make < TH1F > ("muValidGemHitsInTime_Glb_noise", "Valid Gem Hits of In Time for Glb  muons from noise", 50, 0, 50);
  muValidGemHitsInTime_Glb_ghost = fs->make < TH1F > ("muValidGemHitsInTime_Glb_ghost", "Valid Gem Hits of In Time for Glb  muons from ghosts", 50, 0, 50);

////////////////////

  muOOT_Sta_punch = fs->make < TH1F > ("muOOT_Sta_punch", "OOT for Sta  muons from punch through", 400, -200.,200);
  muOOT_Sta_lightflv = fs->make < TH1F > ("muOOT_Sta_lightflv", "OOT for Sta  muons from light flavour muons", 400, -200.,200);
  muOOT_Sta_heavyflv = fs->make < TH1F > ("muOOT_Sta_heavyflv", "OOT for Sta  muons from heavy flavour muons", 400, -200.,200);
  muOOT_Sta_primary = fs->make < TH1F > ("muOOT_Sta_primary", "OOT for Sta  muons from primary", 400, -200.,200);
  muOOT_Sta_noise = fs->make < TH1F > ("muOOT_Sta_noise", "OOT for Sta  muons from noise", 400, -200.,200);
  muOOT_Sta_ghost = fs->make < TH1F > ("muOOT_Sta_ghost", "OOT for Sta  muons from ghosts", 400, -200.,200);

  muPtOOT_Sta_punch = fs->make < TH1F > ("muPtOOT_Sta_punch", "Pt of OOT for Sta  muons from punch through", 1000, 0., 500.);
  muPtOOT_Sta_lightflv = fs->make < TH1F > ("muPtOOT_Sta_lightflv", "Pt of OOT for Sta  muons from light flavour muons", 1000, 0., 500.);
  muPtOOT_Sta_heavyflv = fs->make < TH1F > ("muPtOOT_Sta_heavyflv", "Pt of OOT for Sta  muons from heavy flavour muons", 1000, 0., 500.);
  muPtOOT_Sta_primary = fs->make < TH1F > ("muPtOOT_Sta_primary", "Pt of OOT for Sta  muons from primary", 1000, 0., 500.);
  muPtOOT_Sta_noise = fs->make < TH1F > ("muPtOOT_Sta_noise", "Pt of OOT for Sta  muons from noise", 1000, 0., 500.);
  muPtOOT_Sta_ghost = fs->make < TH1F > ("muPtOOT_Sta_ghost", "Pt of OOT for Sta  muons from ghosts", 1000, 0., 500.);

  muEtaOOT_Sta_punch = fs->make < TH1F > ("muEtaOOT_Sta_punch", "Eta of OOT for Sta  muons from punch through", 70, -3.5, 3.5);
  muEtaOOT_Sta_lightflv = fs->make < TH1F > ("muEtaOOT_Sta_lightflv", "Eta of OOT for Sta  muons from light flavour muons", 70, -3.5, 3.5);
  muEtaOOT_Sta_heavyflv = fs->make < TH1F > ("muEtaOOT_Sta_heavyflv", "Eta of OOT for Sta  muons from heavy flavour muons", 70, -3.5, 3.5);
  muEtaOOT_Sta_primary = fs->make < TH1F > ("muEtaOOT_Sta_primary", "Eta of OOT for Sta  muons from primary", 70, -3.5, 3.5);
  muEtaOOT_Sta_noise = fs->make < TH1F > ("muEtaOOT_Sta_noise", "Eta of OOT for Sta  muons from noise", 70, -3.5, 3.5);
  muEtaOOT_Sta_ghost = fs->make < TH1F > ("muEtaOOT_Sta_ghost", "Eta of OOT for Sta  muons from ghosts", 70, -3.5, 3.5);

  muValidHitsOOT_Sta_punch = fs->make < TH1F > ("muValidHitsOOT_Sta_punch", "Valid Hits of OOT for Sta  muons from punch through", 50, 0, 50);
  muValidHitsOOT_Sta_lightflv = fs->make < TH1F > ("muValidHitsOOT_Sta_lightflv", "Valid Hits of OOT for Sta  muons from light flavour muons", 50, 0, 50 );
  muValidHitsOOT_Sta_heavyflv = fs->make < TH1F > ("muValidHitsOOT_Sta_heavyflv", "Valid Hits of OOT for Sta  muons from heavy flavour muons", 50, 0, 50);
  muValidHitsOOT_Sta_primary = fs->make < TH1F > ("muValidHitsOOT_Sta_primary", "Valid Hits of OOT for Sta  muons from primary", 50, 0, 50);
  muValidHitsOOT_Sta_noise = fs->make < TH1F > ("muValidHitsOOT_Sta_noise", "Valid Hits of OOT for Sta  muons from noise", 50, 0, 50);
  muValidHitsOOT_Sta_ghost = fs->make < TH1F > ("muValidHitsOOT_Sta_ghost", "Valid Hits of OOT for Sta  muons from ghosts", 50, 0, 50);

  muValidGemHitsOOT_Sta_punch = fs->make < TH1F > ("muValidGemHitsOOT_Sta_punch", "Valid Gem Hits of OOT for Sta  muons from punch through", 50, 0, 50);
  muValidGemHitsOOT_Sta_lightflv = fs->make < TH1F > ("muValidGemHitsOOT_Sta_lightflv", "Valid Gem Hits of OOT for Sta  muons from light flavour muons", 50, 0, 50 );
  muValidGemHitsOOT_Sta_heavyflv = fs->make < TH1F > ("muValidGemHitsOOT_Sta_heavyflv", "Valid Gem Hits of OOT for Sta  muons from heavy flavour muons", 50, 0, 50);
  muValidGemHitsOOT_Sta_primary = fs->make < TH1F > ("muValidGemHitsOOT_Sta_primary", "Valid Gem Hits of OOT for Sta  muons from primary", 50, 0, 50);
  muValidGemHitsOOT_Sta_noise = fs->make < TH1F > ("muValidGemHitsOOT_Sta_noise", "Valid Gem Hits of OOT for Sta  muons from noise", 50, 0, 50);
  muValidGemHitsOOT_Sta_ghost = fs->make < TH1F > ("muValidGemHitsOOT_Sta_ghost", "Valid Gem Hits of OOT for Sta  muons from ghosts", 50, 0, 50);
 
  muInTime_Sta_punch = fs->make < TH1F > ("muInTime_Sta_punch", "In Time for Sta  muons from punch through", 400, -200.,200);
  muInTime_Sta_lightflv = fs->make < TH1F > ("muInTime_Sta_lightflv", "In Time for Sta  muons from light flavour muons", 400, -200.,200);
  muInTime_Sta_heavyflv = fs->make < TH1F > ("muInTime_Sta_heavyflv", "In Time for Sta  muons from heavy flavour muons", 400, -200.,200);
  muInTime_Sta_primary = fs->make < TH1F > ("muInTime_Sta_primary", "In Time for Sta  muons from primary", 400, -200.,200);
  muInTime_Sta_noise = fs->make < TH1F > ("muInTime_Sta_noise", "In Time for Sta  muons from noise", 400, -200.,200);
  muInTime_Sta_ghost = fs->make < TH1F > ("muInTime_Sta_ghost", "In Time for Sta  muons from ghosts", 400, -200.,200);

  muPtInTime_Sta_punch = fs->make < TH1F > ("muPtInTime_Sta_punch", "Pt of In Time for Sta  muons from punch through", 1000, 0., 500.);
  muPtInTime_Sta_lightflv = fs->make < TH1F > ("muPtInTime_Sta_lightflv", "Pt of In Time for Sta  muons from light flavour muons", 1000, 0., 500.);
  muPtInTime_Sta_heavyflv = fs->make < TH1F > ("muPtInTime_Sta_heavyflv", "Pt of In Time for Sta  muons from heavy flavour muons", 1000, 0., 500.);
  muPtInTime_Sta_primary = fs->make < TH1F > ("muPtInTime_Sta_primary", "Pt of In Time for Sta  muons from primary", 1000, 0., 500.);
  muPtInTime_Sta_noise = fs->make < TH1F > ("muPtInTime_Sta_noise", "Pt of In Time for Sta  muons from noise", 1000, 0., 500.);
  muPtInTime_Sta_ghost = fs->make < TH1F > ("muPtInTime_Sta_ghost", "Pt of In Time for Sta muons from ghosts", 1000, 0., 500.);

  muEtaInTime_Sta_punch = fs->make < TH1F > ("muEtaInTime_Sta_punch", "Eta of In Time for Sta  muons from punch through", 70, -3.5, 3.5);
  muEtaInTime_Sta_lightflv = fs->make < TH1F > ("muEtaInTime_Sta_lightflv", "Eta of In Time for Sta  muons from light flavour muons", 70, -3.5, 3.5);
  muEtaInTime_Sta_heavyflv = fs->make < TH1F > ("muEtaInTime_Sta_heavyflv", "Eta of In Time for Sta  muons from heavy flavour muons", 70, -3.5, 3.5);
  muEtaInTime_Sta_primary = fs->make < TH1F > ("muEtaInTme_Sta_primary", "Eta of In Time for Sta  muons from primary", 70, -3.5, 3.5);
  muEtaInTime_Sta_noise = fs->make < TH1F > ("muEtaInTime_Sta_noise", "Eta of In Time for Sta  muons from noise", 70, -3.5, 3.5);
  muEtaInTime_Sta_ghost = fs->make < TH1F > ("muEtaInTime_Sta_ghost", "Eta of In Time for Sta  muons from ghosts", 70, -3.5, 3.5);

  muValidHitsInTime_Sta_punch = fs->make < TH1F > ("muValidHitsInTime_Sta_punch", "Valid Hits of In Time for Sta  muons from punch through", 50, 0, 50);
  muValidHitsInTime_Sta_lightflv = fs->make < TH1F > ("muValidHitsInTime_Sta_lightflv", "Valid Hits of In Time for Sta  muons from light flavour muons", 50, 0, 50 );
  muValidHitsInTime_Sta_heavyflv = fs->make < TH1F > ("muValidHitsInTime_Sta_heavyflv", "Valid Hits of In Time for Sta  muons from heavy flavour muons", 50, 0, 50);
  muValidHitsInTime_Sta_primary = fs->make < TH1F > ("muValidHitsInTime_Sta_primary", "Valid Hits of In Time for Sta  muons from primary", 50, 0, 50);
  muValidHitsInTime_Sta_noise = fs->make < TH1F > ("muValidHitsInTime_Sta_noise", "Valid Hits of In Time for Sta  muons from noise", 50, 0, 50);
  muValidHitsInTime_Sta_ghost = fs->make < TH1F > ("muValidHitsInTime_Sta_ghost", "Valid Hits of In Time for Sta  muons from ghosts", 50, 0, 50);

  muValidGemHitsInTime_Sta_punch = fs->make < TH1F > ("muValidGemHitsInTime_Sta_punch", "Valid Gem Hits of In Time for Sta  muons from punch through", 50, 0, 50);
  muValidGemHitsInTime_Sta_lightflv = fs->make < TH1F > ("muValidGemHitsInTime_Sta_lightflv", "Valid Gem Hits of In Time for Sta  muons from light flavour muons", 50, 0, 50 );
  muValidGemHitsInTime_Sta_heavyflv = fs->make < TH1F > ("muValidGemHitsInTime_Sta_heavyflv", "Valid Gem Hits of In Time for Sta  muons from heavy flavour muons", 50, 0, 50);
  muValidGemHitsInTime_Sta_primary = fs->make < TH1F > ("muValidGemHitsInTime_Sta_primary", "Valid Gem Hits of In Time for Sta  muons from primary", 50, 0, 50);
  muValidGemHitsInTime_Sta_noise = fs->make < TH1F > ("muValidGemHitsInTime_Sta_noise", "Valid Gem Hits of In Time for Sta  muons from noise", 50, 0, 50);
  muValidGemHitsInTime_Sta_ghost = fs->make < TH1F > ("muValidGemHitsInTime_Sta_ghost", "Valid Gem Hits of In Time for Sta  muons from ghosts", 50, 0, 50);

///////////////////////////////////////////////////////

  muOOT_Tight_punch = fs->make < TH1F > ("muOOT_Tight_punch", "OOT for tight  muons from punch through", 400, -200.,200);
  muOOT_Tight_lightflv = fs->make < TH1F > ("muOOT_Tight_lightflv", "OOT for tight  muons from light flavour muons", 400, -200.,200);
  muOOT_Tight_heavyflv = fs->make < TH1F > ("muOOT_Tight_heavyflv", "OOT for tight  muons from heavy flavour muons", 400, -200.,200);
  muOOT_Tight_primary = fs->make < TH1F > ("muOOT_Tight_primary", "OOT for tight  muons from primary", 400, -200.,200);
  muOOT_Tight_noise = fs->make < TH1F > ("muOOT_Tight_noise", "OOT for tight  muons from noise", 400, -200.,200);
  muOOT_Tight_ghost = fs->make < TH1F > ("muOOT_Tight_ghost", "OOT for tight  muons from ghosts", 400, -200.,200);

  muPtOOT_Tight_punch = fs->make < TH1F > ("muPtOOT_Tight_punch", "Pt of OOT for tight  muons from punch through", 1000, 0., 500.);
  muPtOOT_Tight_lightflv = fs->make < TH1F > ("muPtOOT_Tight_lightflv", "Pt of OOT for tight  muons from light flavour muons", 1000, 0., 500.);
  muPtOOT_Tight_heavyflv = fs->make < TH1F > ("muPtOOT_Tight_heavyflv", "Pt of OOT for tight  muons from heavy flavour muons", 1000, 0., 500.);
  muPtOOT_Tight_primary = fs->make < TH1F > ("muPtOOT_Tight_primary", "Pt of OOT for tight  muons from primary", 1000, 0., 500.);
  muPtOOT_Tight_noise = fs->make < TH1F > ("muPtOOT_Tight_noise", "Pt of OOT for tight  muons from noise", 1000, 0., 500.);
  muPtOOT_Tight_ghost = fs->make < TH1F > ("muPtOOT_Tight_ghost", "Pt of OOT for tight  muons from ghosts", 1000, 0., 500.);

  muEtaOOT_Tight_punch = fs->make < TH1F > ("muEtaOOT_Tight_punch", "Eta of OOT for tight  muons from punch through", 70, -3.5, 3.5);
  muEtaOOT_Tight_lightflv = fs->make < TH1F > ("muEtaOOT_Tight_lightflv", "Eta of OOT for tight  muons from light flavour muons", 70, -3.5, 3.5);
  muEtaOOT_Tight_heavyflv = fs->make < TH1F > ("muEtaOOT_Tight_heavyflv", "Eta of OOT for tight  muons from heavy flavour muons", 70, -3.5, 3.5);
  muEtaOOT_Tight_primary = fs->make < TH1F > ("muEtaOOT_Tight_primary", "Eta of OOT for tight  muons from primary", 70, -3.5, 3.5);
  muEtaOOT_Tight_noise = fs->make < TH1F > ("muEtaOOT_Tight_noise", "Eta of OOT for tight  muons from noise", 70, -3.5, 3.5);
  muEtaOOT_Tight_ghost = fs->make < TH1F > ("muEtaOOT_Tight_ghost", "Eta of OOT for tight  muons from ghosts", 70, -3.5, 3.5);

  muValidHitsOOT_Tight_punch = fs->make < TH1F > ("muValidHitsOOT_Tight_punch", "Valid Hits of OOT for tight  muons from punch through", 50, 0, 50);
  muValidHitsOOT_Tight_lightflv = fs->make < TH1F > ("muValidHitsOOT_Tight_lightflv", "Valid Hits of OOT for tight  muons from light flavour muons", 50, 0, 50 );
  muValidHitsOOT_Tight_heavyflv = fs->make < TH1F > ("muValidHitsOOT_Tight_heavyflv", "Valid Hits of OOT for tight  muons from heavy flavour muons", 50, 0, 50);
  muValidHitsOOT_Tight_primary = fs->make < TH1F > ("muValidHitsOOT_Tight_primary", "Valid Hits of OOT for tight  muons from primary", 50, 0, 50);
  muValidHitsOOT_Tight_noise = fs->make < TH1F > ("muValidHitsOOT_Tight_noise", "Valid Hits of OOT for tight  muons from noise", 50, 0, 50);
  muValidHitsOOT_Tight_ghost = fs->make < TH1F > ("muValidHitsOOT_Tight_ghost", "Valid Hits of OOT for tight  muons from ghosts", 50, 0, 50);

  muValidGemHitsOOT_Tight_punch = fs->make < TH1F > ("muValidGemHitsOOT_Tight_punch", "Valid Gem Hits of OOT for tight  muons from punch through", 50, 0, 50);
  muValidGemHitsOOT_Tight_lightflv = fs->make < TH1F > ("muValidGemHitsOOT_Tight_lightflv", "Valid Gem Hits of OOT for tight  muons from light flavour muons", 50, 0, 50 );
  muValidGemHitsOOT_Tight_heavyflv = fs->make < TH1F > ("muValidGemHitsOOT_Tight_heavyflv", "Valid Gem Hits of OOT for tight  muons from heavy flavour muons", 50, 0, 50);
  muValidGemHitsOOT_Tight_primary = fs->make < TH1F > ("muValidGemHitsOOT_Tight_primary", "Valid Gem Hits of OOT for tight  muons from primary", 50, 0, 50);
  muValidGemHitsOOT_Tight_noise = fs->make < TH1F > ("muValidGemHitsOOT_Tight_noise", "Valid Gem Hits of OOT for tight  muons from noise", 50, 0, 50);
  muValidGemHitsOOT_Tight_ghost = fs->make < TH1F > ("muValidGemHitsOOT_Tight_ghost", "Valid Gem Hits of OOT for tight  muons from ghosts", 50, 0, 50);

  muInTime_Tight_punch = fs->make < TH1F > ("muInTime_Tight_punch", "In Timefor tight  muons from punch through", 400, -200.,200);
  muInTime_Tight_lightflv = fs->make < TH1F > ("muInTime_Tight_lightflv", "In Timefor tight  muons from light flavour muons", 400, -200.,200);
  muInTime_Tight_heavyflv = fs->make < TH1F > ("muInTime_Tight_heavyflv", "In Timefor tight  muons from heavy flavour muons", 400, -200.,200);
  muInTime_Tight_primary = fs->make < TH1F > ("muInTime_Tight_primary", "In Timefor tight  muons from primary", 400, -200.,200);
  muInTime_Tight_noise = fs->make < TH1F > ("muInTime_Tight_noise", "In Timefor tight  muons from noise", 400, -200.,200);
  muInTime_Tight_ghost = fs->make < TH1F > ("muInTime_Tight_ghost", "In Timefor tight  muons from ghosts", 400, -200.,200);

  muPtInTime_Tight_punch = fs->make < TH1F > ("muPtInTime_Tight_punch", "Pt of In Timefor tight  muons from punch through", 1000, 0., 500.);
  muPtInTime_Tight_lightflv = fs->make < TH1F > ("muPtInTime_Tight_lightflv", "Pt of In Timefor tight  muons from light flavour muons", 1000, 0., 500.);
  muPtInTime_Tight_heavyflv = fs->make < TH1F > ("muPtInTime_Tight_heavyflv", "Pt of In Timefor tight  muons from heavy flavour muons", 1000, 0., 500.);
  muPtInTime_Tight_primary = fs->make < TH1F > ("muPtInTime_Tight_primary", "Pt of In Timefor tight  muons from primary", 1000, 0., 500.);
  muPtInTime_Tight_noise = fs->make < TH1F > ("muPtInTime_Tight_noise", "Pt of In Timefor tight  muons from noise", 1000, 0., 500.);
  muPtInTime_Tight_ghost = fs->make < TH1F > ("muPtInTime_Tight_ghost", "Pt of In Timefor tight muons from ghosts", 1000, 0., 500.);

  muEtaInTime_Tight_punch = fs->make < TH1F > ("muEtaInTime_Tight_punch", "Eta of In Timefor tight  muons from punch through", 70, -3.5, 3.5);
  muEtaInTime_Tight_lightflv = fs->make < TH1F > ("muEtaInTime_Tight_lightflv", "Eta of In Timefor tight  muons from light flavour muons", 70, -3.5, 3.5);
  muEtaInTime_Tight_heavyflv = fs->make < TH1F > ("muEtaInTime_Tight_heavyflv", "Eta of In Timefor tight  muons from heavy flavour muons", 70, -3.5, 3.5);
  muEtaInTime_Tight_primary = fs->make < TH1F > ("muEtaInTme_Tight_primary", "Eta of In Timefor tight  muons from primary", 70, -3.5, 3.5);
  muEtaInTime_Tight_noise = fs->make < TH1F > ("muEtaInTime_Tight_noise", "Eta of In Timefor tight  muons from noise", 70, -3.5, 3.5);
  muEtaInTime_Tight_ghost = fs->make < TH1F > ("muEtaInTime_Tight_ghost", "Eta of In Timefor tight  muons from ghosts", 70, -3.5, 3.5);

  muValidHitsInTime_Tight_punch = fs->make < TH1F > ("muValidHitsInTime_Tight_punch", "Valid Hits of In Time for tight  muons from punch through", 50, 0, 50);
  muValidHitsInTime_Tight_lightflv = fs->make < TH1F > ("muValidHitsInTime_Tight_lightflv", "Valid Hits of In Time for tight  muons from light flavour muons", 50, 0, 50 );
  muValidHitsInTime_Tight_heavyflv = fs->make < TH1F > ("muValidHitsInTime_Tight_heavyflv", "Valid Hits of In Time for tight  muons from heavy flavour muons", 50, 0, 50);
  muValidHitsInTime_Tight_primary = fs->make < TH1F > ("muValidHitsInTime_Tight_primary", "Valid Hits of In Time for tight  muons from primary", 50, 0, 50);
  muValidHitsInTime_Tight_noise = fs->make < TH1F > ("muValidHitsInTime_Tight_noise", "Valid Hits of In Time for tight  muons from noise", 50, 0, 50);
  muValidHitsInTime_Tight_ghost = fs->make < TH1F > ("muValidHitsInTime_Tight_ghost", "Valid Hits of In Time for tight  muons from ghosts", 50, 0, 50);

  muValidGemHitsInTime_Tight_punch = fs->make < TH1F > ("muValidGemHitsInTime_Tight_punch", "Valid Gem Hits of In Time for tight  muons from punch through", 50, 0, 50);
  muValidGemHitsInTime_Tight_lightflv = fs->make < TH1F > ("muValidGemHitsInTime_Tight_lightflv", "Valid Gem Hits of In Time for tight  muons from light flavour muons", 50, 0, 50 );
  muValidGemHitsInTime_Tight_heavyflv = fs->make < TH1F > ("muValidGemHitsInTime_Tight_heavyflv", "Valid Gem Hits of In Time for tight  muons from heavy flavour muons", 50, 0, 50);
  muValidGemHitsInTime_Tight_primary = fs->make < TH1F > ("muValidGemHitsInTime_Tight_primary", "Valid Gem Hits of In Time for tight  muons from primary", 50, 0, 50);
  muValidGemHitsInTime_Tight_noise = fs->make < TH1F > ("muValidGemHitsInTime_Tight_noise", "Valid Gem Hits of In Time for tight  muons from noise", 50, 0, 50);
  muValidGemHitsInTime_Tight_ghost = fs->make < TH1F > ("muValidGemHitsInTime_Tight_ghost", "Valid Gem Hits of In Time for tight  muons from ghosts", 50, 0, 50);

///////////////////////////////////////////////////////////////////77

 
  muOOT_Loose_punch = fs->make < TH1F > ("muOOT_Loose_punch", "OOT for loose muons from punch through", 400, -200.,200);
  muOOT_Loose_lightflv = fs->make < TH1F > ("muOOT_Loose_lightflv", "OOT for loose muons from light flavour muons", 400, -200.,200);
  muOOT_Loose_heavyflv = fs->make < TH1F > ("muOOT_Loose_heavyflv", "OOT for loose  muons from heavy flavour muons", 400, -200.,200);
  muOOT_Loose_primary = fs->make < TH1F > ("muOOT_Loose_primary", "OOT for loose muons from primary", 400, -200.,200);
  muOOT_Loose_noise = fs->make < TH1F > ("muOOT_Loose_noise", "OOT for loose muons from noise", 400, -200.,200);
  muOOT_Loose_ghost = fs->make < TH1F > ("muOOT_Loose_ghost", "OOT for loose muons from ghosts", 400, -200.,200);

  muPtOOT_Loose_punch = fs->make < TH1F > ("muPtOOT_Loose_punch", "Pt of OOT for loose  muons from punch through", 1000, 0., 500.);
  muPtOOT_Loose_lightflv = fs->make < TH1F > ("muPtOOT_Loose_lightflv", "Pt of OOT for loose  muons from light flavour muons", 1000, 0., 500.);
  muPtOOT_Loose_heavyflv = fs->make < TH1F > ("muPtOOT_Loose_heavyflv", "Pt of OOT for loose  muons from heavy flavour muons", 1000, 0., 500.);
  muPtOOT_Loose_primary = fs->make < TH1F > ("muPtOOT_Loose_primary", "Pt of OOT for loose  muons from primary", 1000, 0., 500.);
  muPtOOT_Loose_noise = fs->make < TH1F > ("muPtOOT_Loose_noise", "Pt of OOT for loose  muons from noise", 1000, 0., 500.);
  muPtOOT_Loose_ghost = fs->make < TH1F > ("muPtOOT_Loose_ghost", "Pt of OOT for loose  muons from ghosts", 1000, 0., 500.);

  muEtaOOT_Loose_punch = fs->make < TH1F > ("muEtaOOT_Loose_punch", "Eta of OOT for loose  muons from punch through", 70, -3.5, 3.5);
  muEtaOOT_Loose_lightflv = fs->make < TH1F > ("muEtaOOT_Loose_lightflv", "Eta of OOT for loose  muons from light flavour muons", 70, -3.5, 3.5);
  muEtaOOT_Loose_heavyflv = fs->make < TH1F > ("muEtaOOT_Loose_heavyflv", "Eta of OOT for loose  muons from heavy flavour muons", 70, -3.5, 3.5);
  muEtaOOT_Loose_primary = fs->make < TH1F > ("muEtaOOT_Loose_primary", "Eta of OOT for loose  muons from primary", 70, -3.5, 3.5);
  muEtaOOT_Loose_noise = fs->make < TH1F > ("muEtaOOT_Loose_noise", "Eta of OOT for loose  muons from noise", 70, -3.5, 3.5);
  muEtaOOT_Loose_ghost = fs->make < TH1F > ("muEtaOOT_Loose_ghost", "Eta of OOT for loose  muons from ghosts", 70, -3.5, 3.5);

  muValidHitsOOT_Loose_punch = fs->make < TH1F > ("muValidHitsOOT_Loose_punch", "Valid Hits of OOT for loose  muons from punch through", 50, 0, 50);
  muValidHitsOOT_Loose_lightflv = fs->make < TH1F > ("muValidHitsOOT_Loose_lightflv", "Valid Hits of OOT for loose  muons from light flavour muons", 50, 0, 50 );
  muValidHitsOOT_Loose_heavyflv = fs->make < TH1F > ("muValidHitsOOT_Loose_heavyflv", "Valid Hits of OOT for loose  muons from heavy flavour muons", 50, 0, 50);
  muValidHitsOOT_Loose_primary = fs->make < TH1F > ("muValidHitsOOT_Loose_primary", "Valid Hits of OOT for loose  muons from primary", 50, 0, 50);
  muValidHitsOOT_Loose_noise = fs->make < TH1F > ("muValidHitsOOT_Loose_noise", "Valid Hits of OOT for loose  muons from noise", 50, 0, 50);
  muValidHitsOOT_Loose_ghost = fs->make < TH1F > ("muValidHitsOOT_Loose_ghost", "Valid Hits of OOT for loose  muons from ghosts", 50, 0, 50);

  muValidGemHitsOOT_Loose_punch = fs->make < TH1F > ("muValidGemHitsOOT_Loose_punch", "Valid Gem Hits of OOT for loose  muons from punch through", 50, 0, 50);
  muValidGemHitsOOT_Loose_lightflv = fs->make < TH1F > ("muValidGemHitsOOT_Loose_lightflv", "Valid Gem Hits of OOT for loose  muons from light flavour muons", 50, 0, 50 );
  muValidGemHitsOOT_Loose_heavyflv = fs->make < TH1F > ("muValidGemHitsOOT_Loose_heavyflv", "Valid Gem Hits of OOT for loose  muons from heavy flavour muons", 50, 0, 50);
  muValidGemHitsOOT_Loose_primary = fs->make < TH1F > ("muValidGemHitsOOT_Loose_primary", "Valid Gem Hits of OOT for loose  muons from primary", 50, 0, 50);
  muValidGemHitsOOT_Loose_noise = fs->make < TH1F > ("muValidGemHitsOOT_Loose_noise", "Valid Gem Hits of OOT for loose  muons from noise", 50, 0, 50);
  muValidGemHitsOOT_Loose_ghost = fs->make < TH1F > ("muValidGemHitsOOT_Loose_ghost", "Valid Gem Hits of OOT for loose  muons from ghosts", 50, 0, 50);

  muInTime_Loose_punch = fs->make < TH1F > ("muInTime_Loose_punch", "In Time for loose  muons from punch through", 400, -200.,200);
  muInTime_Loose_lightflv = fs->make < TH1F > ("muInTime_Loose_lightflv", "In Time for loose  muons from light flavour muons", 400, -200.,200);
  muInTime_Loose_heavyflv = fs->make < TH1F > ("muInTime_Loose_heavyflv", "In Time for loose  muons from heavy flavour muons", 400, -200.,200);
  muInTime_Loose_primary = fs->make < TH1F > ("muInTime_Loose_primary", "In Time for loose  muons from primary", 400, -200.,200);
  muInTime_Loose_noise = fs->make < TH1F > ("muInTime_Loose_noise", "In Time for loose  muons from noise", 400, -200.,200);
  muInTime_Loose_ghost = fs->make < TH1F > ("muInTime_Loose_ghost", "In Time for loose  muons from ghosts", 400, -200.,200);

  muPtInTime_Loose_punch = fs->make < TH1F > ("muPtInTime_Loose_punch", "Pt of In Time for loose  muons from punch through", 1000, 0., 500.);
  muPtInTime_Loose_lightflv = fs->make < TH1F > ("muPtInTime_Loose_lightflv", "Pt of In Time for loose  muons from light flavour muons", 1000, 0., 500.);
  muPtInTime_Loose_heavyflv = fs->make < TH1F > ("muPtInTime_Loose_heavyflv", "Pt of In Time for loose  muons from heavy flavour muons", 1000, 0., 500.);
  muPtInTime_Loose_primary = fs->make < TH1F > ("muPtInTime_Loose_primary", "Pt of In Time for loose  muons from primary", 1000, 0., 500.);
  muPtInTime_Loose_noise = fs->make < TH1F > ("muPtInTime_Loose_noise", "Pt of In Time for loose  muons from noise", 1000, 0., 500.);
  muPtInTime_Loose_ghost = fs->make < TH1F > ("muPtInTime_Loose_ghost", "Pt of In Time for loose muons from ghosts", 1000, 0., 500.);

  muEtaInTime_Loose_punch = fs->make < TH1F > ("muEtaInTime_Loose_punch", "Eta of In Time for loose  muons from punch through", 70, -3.5, 3.5);
  muEtaInTime_Loose_lightflv = fs->make < TH1F > ("muEtaInTime_Loose_lightflv", "Eta of In Time for loose  muons from light flavour muons", 70, -3.5, 3.5);
  muEtaInTime_Loose_heavyflv = fs->make < TH1F > ("muEtaInTime_Loose_heavyflv", "Eta of In Time for loose  muons from heavy flavour muons", 70, -3.5, 3.5);
  muEtaInTime_Loose_primary = fs->make < TH1F > ("muEtaInTme_Loose_primary", "Eta of In Time for loose  muons from primary", 70, -3.5, 3.5);
  muEtaInTime_Loose_noise = fs->make < TH1F > ("muEtaInTime_Loose_noise", "Eta of In Time for loose  muons from noise", 70, -3.5, 3.5);
  muEtaInTime_Loose_ghost = fs->make < TH1F > ("muEtaInTime_Loose_ghost", "Eta of In Time for loose  muons from ghosts", 70, -3.5, 3.5);

  muValidHitsInTime_Loose_punch = fs->make < TH1F > ("muValidHitsInTime_Loose_punch", "Valid Hits of In Time for loose  muons from punch through", 50, 0, 50);
  muValidHitsInTime_Loose_lightflv = fs->make < TH1F > ("muValidHitsInTime_Loose_lightflv", "Valid Hits of In Time for loose  muons from light flavour muons", 50, 0, 50 );
  muValidHitsInTime_Loose_heavyflv = fs->make < TH1F > ("muValidHitsInTime_Loose_heavyflv", "Valid Hits of In Time for loose  muons from heavy flavour muons", 50, 0, 50);
  muValidHitsInTime_Loose_primary = fs->make < TH1F > ("muValidHitsInTime_Loose_primary", "Valid Hits of In Time for loose  muons from primary", 50, 0, 50);
  muValidHitsInTime_Loose_noise = fs->make < TH1F > ("muValidHitsInTime_Loose_noise", "Valid Hits of In Time for loose  muons from noise", 50, 0, 50);
  muValidHitsInTime_Loose_ghost = fs->make < TH1F > ("muValidHitsInTime_Loose_ghost", "Valid Hits of In Time for loose  muons from ghosts", 50, 0, 50);

  muValidGemHitsInTime_Loose_punch = fs->make < TH1F > ("muValidGemHitsInTime_Loose_punch", "Valid Gem Hits of In Time for loose  muons from punch through", 50, 0, 50);
  muValidGemHitsInTime_Loose_lightflv = fs->make < TH1F > ("muValidGemHitsInTime_Loose_lightflv", "Valid Gem Hits of In Time for loose  muons from light flavour muons", 50, 0, 50 );
  muValidGemHitsInTime_Loose_heavyflv = fs->make < TH1F > ("muValidGemHitsInTime_Loose_heavyflv", "Valid Gem Hits of In Time for loose  muons from heavy flavour muons", 50, 0, 50);
  muValidGemHitsInTime_Loose_primary = fs->make < TH1F > ("muValidGemHitsInTime_Loose_primary", "Valid Gem Hits of In Time for loose  muons from primary", 50, 0, 50);
  muValidGemHitsInTime_Loose_noise = fs->make < TH1F > ("muValidGemHitsInTime_Loose_noise", "Valid Gem Hits of In Time for loose  muons from noise", 50, 0, 50);
  muValidGemHitsInTime_Loose_ghost = fs->make < TH1F > ("muValidGemHitsInTime_Loose_ghost", "Valid Gem Hits of In Time for loose  muons from ghosts", 50, 0, 50);
}
void MuonMCClassifAndAna::endJob()
{
  double myScale = N_Muons->GetBinContent(65);
  N_Muons_norm->Scale(1. / myScale);
}
void MuonMCClassifAndAna::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
//The following says we do not know what parameters are allowed so do no validation
// Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE (MuonMCClassifAndAna);

