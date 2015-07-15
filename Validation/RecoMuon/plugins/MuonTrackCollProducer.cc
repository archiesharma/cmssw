
#include "Validation/RecoMuon/plugins/MuonTrackCollProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h" 
#include <sstream>

std::vector<double> MuonTrackCollProducer::findSimVtx(edm::Event& iEvent){

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);
  std::vector<double> vtxCoord;
  vtxCoord.push_back(0);
  vtxCoord.push_back(0);
  vtxCoord.push_back(0);
  vtxCoord.push_back(0);
  vtxCoord.push_back(0);
  vtxCoord.push_back(0);
  vtxCoord.push_back(0);

  if(genParticles.isValid()){

  	for(reco::GenParticleCollection::const_iterator itg = genParticles->begin(); itg != genParticles->end(); ++itg ){

		int id = itg->pdgId();
		int status = itg->status();
		//int nDaughters = itg->numberOfDaughters();
		int nMothers = itg->numberOfMothers();
		//double phiGen = itg->phi();
		//double etaGen = itg->eta();
		//std::cout<<"id "<<id<<" "<<status<<" "<<nMothers<<" "<<phiGen<<" "<<etaGen<<std::endl;

		if((abs(id) == 23 || abs(id) == 22) && status == 3){

	 		vtxCoord[0] = 1;

			vtxCoord[4] = (double)(itg->vx()); 
			vtxCoord[5] = (double)(itg->vy());
			vtxCoord[6] = (double)(itg->vz());

		}

		//if(fabs(id) == 13) std::cout<<"ID "<<id<<" Status "<<status<<std::endl;

		else if(abs(id) == 13 && status == 1 && nMothers == 0){

	 		vtxCoord[0] = 2;

			vtxCoord[1] = (double)(itg->vx()); 
			vtxCoord[2] = (double)(itg->vy());
			vtxCoord[3] = (double)(itg->vz());

		}

		else if(abs(id) == 12 && status == 1 && nMothers == 0){

	 		vtxCoord[0] = 3;

			vtxCoord[1] = (double)(itg->vx()); 
			vtxCoord[2] = (double)(itg->vy());
			vtxCoord[3] = (double)(itg->vz());

		}

	}

  }

  //std::cout<<vtxCoord.size()<<" "<<vtxCoord[0]<<std::endl;
  return vtxCoord;

}

bool MuonTrackCollProducer::isLoose(edm::Event& iEvent, reco::MuonCollection::const_iterator muon)
{
  bool isPF = muon->isPFMuon();
  bool isGLB = muon->isGlobalMuon();
  bool isTrk = muon->isTrackerMuon();

  return (isPF && (isGLB || isTrk) );
}

bool MuonTrackCollProducer::isSoft(edm::Event& iEvent, reco::MuonCollection::const_iterator muon, bool useIPxy, bool useIPz)
{
  bool result = false;

  if (muon->muonBestTrack().isNonnull() && muon->innerTrack().isNonnull()){

 	edm::Handle<reco::VertexCollection> vertexHandle;
  	iEvent.getByLabel(vxtTag,vertexHandle);
  	const reco::VertexCollection* vertices = vertexHandle.product();

	bool isGood = muon::isGoodMuon((*muon), muon::TMOneStationTight); 
	bool trkLayMeas = muon->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5; 
	bool pxlLayMeas = muon->innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0; 
	bool quality = muon->innerTrack()->quality(reco::Track::highPurity);
	bool ipxy = false;
	bool ipz = false;
	if(vertices->size()!=0 && useIPxy) ipxy = fabs(muon->muonBestTrack()->dxy((*vertices)[0].position())) < 0.2;
	else ipxy = true;
 	if(vertices->size()!=0 && useIPz) ipz = fabs(muon->muonBestTrack()->dz((*vertices)[0].position())) < 0.5;
	else ipz = true;
	if(isGood && trkLayMeas && pxlLayMeas && quality && ipxy && ipz) result = true;

  }

  return result;
}

bool MuonTrackCollProducer::isTight(edm::Event& iEvent, reco::MuonCollection::const_iterator muon, bool useIPxy, bool useIPz)
{
  bool result = false;

  if (muon->muonBestTrack().isNonnull() && muon->innerTrack().isNonnull() && muon->globalTrack().isNonnull()){

	std::vector<double> vtxCoord = findSimVtx(iEvent);
        GlobalPoint point(vtxCoord[1],vtxCoord[2],vtxCoord[3]);
        GlobalPoint pointDY(vtxCoord[4],vtxCoord[5],vtxCoord[6]);

	//double muonX = muon->vx();
	//double muonY = muon->vy();
	//double muonZ = muon->vz();

	double muonZ = pointDY.z();
		
 	edm::Handle<reco::VertexCollection> vertexHandle;
  	iEvent.getByLabel(vxtTag,vertexHandle);
  	const reco::VertexCollection* vertices = vertexHandle.product();

	double distInit = 24;
	int indexFinal = 0;
	for(int i = 0; i < (int)vertices->size(); i++){

		//double vtxX = (*vertices)[i].x();
		//double vtxY = (*vertices)[i].y();
		double vtxZ = (*vertices)[i].z();

		double dist = fabs(muonZ - vtxZ);
		//std::cout<<"dist "<<dist<<std::endl;
		if(dist < distInit){
				
			distInit = dist;
			indexFinal = i;

		}

	}
	//std::cout<<distInit<<" "<<indexFinal<<std::endl;

	double ipxySim = 999;
	double ipzSim = 999;
	
	if(vtxCoord[0] > 1.5 && vtxCoord[0] < 3.5){

        	ipxySim = fabs(muon->muonBestTrack()->dxy(math::XYZPoint(point.x(),point.y(),point.z())));
        	ipzSim = fabs(muon->muonBestTrack()->dz(math::XYZPoint(point.x(),point.y(),point.z())));
	}
	else if(vtxCoord[0] > 0.5 && vtxCoord[0] < 1.5){

		ipxySim = fabs(muon->muonBestTrack()->dxy(math::XYZPoint(pointDY.x(),pointDY.y(),pointDY.z())));
        	ipzSim = fabs(muon->muonBestTrack()->dz(math::XYZPoint(pointDY.x(),pointDY.y(),pointDY.z())));

	}
	bool ipxySimBool = ipxySim < 0.2;
	bool ipzSimBool = ipzSim < 0.5;
        //std::cout<<"vx: "<<point.x()<<" vy: "<<point.y()<<" vz: "<<point.z()<<" |Dxy|: "<<ipxySim<<" "<<ipxySimBool<<" |Dz|: "<<ipzSim<<" "<<ipzSimBool<<std::endl;
        //std::cout<<"vx: "<<pointDY.x()<<" vy: "<<pointDY.y()<<" vz: "<<pointDY.z()<<" |Dxy|: "<<ipxySim<<" "<<ipxySimBool<<" |Dz|: "<<ipzSim<<" "<<ipzSimBool<<std::endl;

	bool trkLayMeas = muon->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5; 
	bool isGlb = muon->isGlobalMuon(); 
	bool isPF = muon->isPFMuon(); 
	bool chi2 = muon->globalTrack()->normalizedChi2() < 10.; 
	bool validHits = muon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0; 
	bool matchedSt = muon->numberOfMatchedStations() > 1;

	bool ipxy = false;
	bool ipz = false;
	if(vertices->size() !=0 && useIPxy == true){
	//if(useIPxy){

		if(vtxCoord[0] > 0.5 && vtxCoord[0] < 1.5) ipxy = fabs(muon->muonBestTrack()->dxy((*vertices)[indexFinal].position())) < 0.2;
		else if(vtxCoord[0] > 2.5 && vtxCoord[0] < 3.5) ipxy = fabs(muon->muonBestTrack()->dxy((*vertices)[0].position())) < 0.2;
		else ipxy = ipxySimBool;
		//ipxy = ipxySimBool;

	}
	else if(vertices->size() == 0 && useIPxy == true) ipxy = false;
	else if(useIPxy == false) ipxy = true;

 	if(vertices->size() !=0 && useIPz == true){
 	//if(useIPz){

		if(vtxCoord[0] > 0.5 && vtxCoord[0] < 1.5) ipz = fabs(muon->muonBestTrack()->dz((*vertices)[indexFinal].position())) < 0.5;
		else if(vtxCoord[0] > 2.5 && vtxCoord[0] < 3.5) ipz = fabs(muon->muonBestTrack()->dz((*vertices)[0].position())) < 0.5;	
		else ipz = ipzSimBool;
		//ipz = ipzSimBool;

        	//std::cout<<"vx: "<<pointDY.x()<<" vy: "<<pointDY.y()<<" vz: "<<pointDY.z()<<" |Dz|: "<<fabs(muon->muonBestTrack()->dz((*vertices)[indexFinal].position()))<<" "<<ipz<<std::endl;

	}
	else if(vertices->size() == 0 && useIPz == true) ipz = false;
	else if(useIPz == false) ipz = true;
	//bool validPxlHit = muon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0;
	bool validPxlHit = muon->innerTrack()->hitPattern().pixelLayersWithMeasurement(3,2) > 0;
	//bool validPxlHit = muon->innerTrack()->hitPattern().pixelLayersWithMeasurement(4,3) > 0;

	if(trkLayMeas && isGlb && isPF && chi2 && validHits && matchedSt && ipxy && ipz && validPxlHit) result = true;

  }

  return result;
}

bool MuonTrackCollProducer::isLoose2(edm::Event& iEvent, reco::MuonCollection::const_iterator muon, bool useIPxy, bool useIPz)
{
  bool result = false;

  if (muon->muonBestTrack().isNonnull() && muon->innerTrack().isNonnull() && muon->globalTrack().isNonnull()){

	std::vector<double> vtxCoord = findSimVtx(iEvent);
        GlobalPoint point(vtxCoord[1],vtxCoord[2],vtxCoord[3]);
        GlobalPoint pointDY(vtxCoord[4],vtxCoord[5],vtxCoord[6]);

	//double muonX = muon->vx();
	//double muonY = muon->vy();
	//double muonZ = muon->vz();

	double muonZ = pointDY.z();
		
 	edm::Handle<reco::VertexCollection> vertexHandle;
  	iEvent.getByLabel(vxtTag,vertexHandle);
  	const reco::VertexCollection* vertices = vertexHandle.product();

	double distInit = 24;
	int indexFinal = 0;
	for(int i = 0; i < (int)vertices->size(); i++){

		//double vtxX = (*vertices)[i].x();
		//double vtxY = (*vertices)[i].y();
		double vtxZ = (*vertices)[i].z();

		double dist = fabs(muonZ - vtxZ);
		//std::cout<<"dist "<<dist<<std::endl;
		if(dist < distInit){
				
			distInit = dist;
			indexFinal = i;

		}

	}
	//std::cout<<distInit<<" "<<indexFinal<<std::endl;

	double ipxySim = 999;
	double ipzSim = 999;
	
	if(vtxCoord[0] > 1.5 && vtxCoord[0] < 3.5){

        	ipxySim = fabs(muon->muonBestTrack()->dxy(math::XYZPoint(point.x(),point.y(),point.z())));
        	ipzSim = fabs(muon->muonBestTrack()->dz(math::XYZPoint(point.x(),point.y(),point.z())));
	}
	else if(vtxCoord[0] > 0.5 && vtxCoord[0] < 1.5){

		ipxySim = fabs(muon->muonBestTrack()->dxy(math::XYZPoint(pointDY.x(),pointDY.y(),pointDY.z())));
        	ipzSim = fabs(muon->muonBestTrack()->dz(math::XYZPoint(pointDY.x(),pointDY.y(),pointDY.z())));

	}
	bool ipxySimBool = ipxySim < 0.2;
	bool ipzSimBool = ipzSim < 0.5;
        //std::cout<<"vx: "<<point.x()<<" vy: "<<point.y()<<" vz: "<<point.z()<<" |Dxy|: "<<ipxySim<<" "<<ipxySimBool<<" |Dz|: "<<ipzSim<<" "<<ipzSimBool<<std::endl;
        //std::cout<<"vx: "<<pointDY.x()<<" vy: "<<pointDY.y()<<" vz: "<<pointDY.z()<<" |Dxy|: "<<ipxySim<<" "<<ipxySimBool<<" |Dz|: "<<ipzSim<<" "<<ipzSimBool<<std::endl;

	bool isGlb = muon->isGlobalMuon(); 
	bool isPF = muon->isPFMuon(); 
  	//bool isTrk = muon->isTrackerMuon();

	bool ipxy = false;
	bool ipz = false;
	if(vertices->size() !=0 && useIPxy == true){
	//if(useIPxy){

		if(vtxCoord[0] > 0.5 && vtxCoord[0] < 1.5) ipxy = fabs(muon->muonBestTrack()->dxy((*vertices)[indexFinal].position())) < 0.2;
		else if(vtxCoord[0] > 2.5 && vtxCoord[0] < 3.5) ipxy = fabs(muon->muonBestTrack()->dxy((*vertices)[0].position())) < 0.2;
		else ipxy = ipxySimBool;
		//ipxy = ipxySimBool;

	}
	else if(vertices->size() == 0 && useIPxy == true) ipxy = false;
	else if(useIPxy == false) ipxy = true;

 	if(vertices->size() !=0 && useIPz == true){
 	//if(useIPz){

		if(vtxCoord[0] > 0.5 && vtxCoord[0] < 1.5) ipz = fabs(muon->muonBestTrack()->dz((*vertices)[indexFinal].position())) < 0.5;
		else if(vtxCoord[0] > 2.5 && vtxCoord[0] < 3.5) ipz = fabs(muon->muonBestTrack()->dz((*vertices)[0].position())) < 0.5;	
		else ipz = ipzSimBool;
		//ipz = ipzSimBool;

        	//std::cout<<"vx: "<<pointDY.x()<<" vy: "<<pointDY.y()<<" vz: "<<pointDY.z()<<" |Dz|: "<<fabs(muon->muonBestTrack()->dz((*vertices)[indexFinal].position()))<<" "<<ipz<<std::endl;

	}
	else if(vertices->size() == 0 && useIPz == true) ipz = false;
	else if(useIPz == false) ipz = true;

	//if((isGlb || isTrk) && isPF && ipxy && ipz) result = true;
	if(isGlb && isPF && ipxy && ipz) result = true;

  }

  return result;
}

MuonTrackCollProducer::MuonTrackCollProducer(const edm::ParameterSet& parset) :
  muonsTag(parset.getParameter< edm::InputTag >("muonsTag")),
  vxtTag(parset.getParameter< edm::InputTag >("vxtTag")),
  useIPxy(parset.getUntrackedParameter< bool >("useIPxy", true)),
  useIPz(parset.getUntrackedParameter< bool >("useIPz", true)),
  selectionTags(parset.getParameter< std::vector<std::string> >("selectionTags")),
  trackType(parset.getParameter< std::string >("trackType")),
  parset_(parset)
{
  produces<reco::TrackCollection>();
  //produces<reco::TrackExtraCollection>();
  //produces<TrackingRecHitCollection>();
}

MuonTrackCollProducer::~MuonTrackCollProducer() {
}

void MuonTrackCollProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  iEvent.getByLabel(muonsTag,muonCollectionH);
  
  std::auto_ptr<reco::TrackCollection> selectedTracks(new reco::TrackCollection);
  //std::auto_ptr<reco::TrackExtraCollection> selectedTrackExtras( new reco::TrackExtraCollection() );
  //std::auto_ptr<TrackingRecHitCollection> selectedTrackHits( new TrackingRecHitCollection() );

  reco::TrackRefProd rTracks = iEvent.getRefBeforePut<reco::TrackCollection>();
  //reco::TrackExtraRefProd rTrackExtras = iEvent.getRefBeforePut<reco::TrackExtraCollection>();
  //TrackingRecHitRefProd rHits = iEvent.getRefBeforePut<TrackingRecHitCollection>();

  //edm::Ref<reco::TrackExtraCollection>::key_type idx = 0;
  //edm::Ref<reco::TrackExtraCollection>::key_type hidx = 0;

  edm::LogVerbatim("MuonTrackProducer") <<"\nThere are "<< muonCollectionH->size() <<" reco::Muons.";
  unsigned int muon_index = 0;
  for(reco::MuonCollection::const_iterator muon = muonCollectionH->begin();
       muon != muonCollectionH->end(); ++muon, muon_index++) {
    edm::LogVerbatim("MuonTrackProducer") <<"\n******* muon index : "<<muon_index;
    
    std::vector<bool> isGood;
    for(unsigned int index=0; index<selectionTags.size(); ++index) {
      isGood.push_back(false);

      muon::SelectionType muonType = muon::selectionTypeFromString(selectionTags[index]);
      isGood[index] = muon::isGoodMuon(*muon, muonType);
    }

    bool isGoodResult=true;
    for(unsigned int index=0; index<isGood.size(); ++index) {
      edm::LogVerbatim("MuonTrackProducer") << "selectionTag = "<<selectionTags[index]<< ": "<<isGood[index]<<"\n";
      isGoodResult *= isGood[index];
    }

    if (isGoodResult) {
      // new copy of Track
      reco::TrackRef trackref;

      bool loose = isLoose(iEvent, muon);
      bool loose2 = isLoose2(iEvent, muon, useIPxy, useIPz);
      bool soft = isSoft(iEvent, muon, useIPxy, useIPz);
      bool tight = isTight(iEvent, muon, useIPxy, useIPz);
      if (trackType == "innerTrack") {
        if (muon->innerTrack().isNonnull()) trackref = muon->innerTrack();
        else continue;
      }
      else if (trackType == "outerTrack") {
        if (muon->outerTrack().isNonnull()) trackref = muon->outerTrack();
        else continue;
      }
      else if (trackType == "globalTrack") {
        if (muon->globalTrack().isNonnull()) trackref = muon->globalTrack();
        else continue;
      }
      else if (trackType == "innerTrackPlusSegments") {
	if (muon->innerTrack().isNonnull()) trackref = muon->innerTrack();
	else continue;
      }
      else if (trackType == "bestMuon") {
	if (muon->muonBestTrack().isNonnull()) trackref = muon->muonBestTrack();
	else continue;
      }
      else if (trackType == "bestMuonLoose") {
	if (muon->muonBestTrack().isNonnull() && loose) trackref = muon->muonBestTrack();
	else continue;
      }
      else if (trackType == "bestMuonLoose2") {
	if (muon->muonBestTrack().isNonnull() && loose2) trackref = muon->muonBestTrack();
	else continue;
      }
      else if (trackType == "bestMuonSoft") {
	if (muon->muonBestTrack().isNonnull() && soft) trackref = muon->muonBestTrack();
	else continue;
      }
      else if (trackType == "bestMuonTight") {
	if (muon->muonBestTrack().isNonnull() && tight) trackref = muon->muonBestTrack();
	else continue;
      }
      else if (trackType == "bestMuonTuneP") {
	if (muon->tunePMuonBestTrack().isNonnull()) trackref = muon->tunePMuonBestTrack();
	else continue;
      }
      else if (trackType == "cocktailBest") {

  	//None, InnerTrack, OuterTrack, CombinedTrack,TPFMS, Picky, DYT
	reco::Muon::MuonTrackType type = muon->muonBestTrackType();

	switch((int) type){

		case 1: 
			if (muon->innerTrack().isNonnull()) trackref = muon->innerTrack();
			else continue;
			break;
		case 3:
			if (muon->globalTrack().isNonnull()) trackref = muon->globalTrack();
			else continue;
			break;
		case 4:
			if (muon->tpfmsTrack().isNonnull()) trackref = muon->tpfmsTrack();
			else continue;
			break;
		case 5:
			if (muon->pickyTrack().isNonnull()) trackref = muon->pickyTrack();
			else continue;
			break;
		default: continue;

	}

      }
      else if (trackType == "cocktailTuneP") {

  	//None, InnerTrack, OuterTrack, CombinedTrack,TPFMS, Picky, DYT
	reco::Muon::MuonTrackType type = muon->tunePMuonBestTrackType();
	
	switch((int) type){

		case 1: 
			if (muon->innerTrack().isNonnull()) trackref = muon->innerTrack();
			else continue;
			break;
		case 3:
			if (muon->globalTrack().isNonnull()) trackref = muon->globalTrack();
			else continue;
			break;
		case 4:
			if (muon->tpfmsTrack().isNonnull()) trackref = muon->tpfmsTrack();
			else continue;
			break;
		case 5:
			if (muon->pickyTrack().isNonnull()) trackref = muon->pickyTrack();
			else continue;
			break;
		default: continue;

	}

      }

      const reco::Track* trk = &(*trackref);
      // pointer to old track:
      reco::Track* newTrk = new reco::Track(*trk);

      //newTrk->setExtra( reco::TrackExtraRef( rTrackExtras, idx++ ) );
      //PropagationDirection seedDir = trk->seedDirection();
      // new copy of track Extras
      /*reco::TrackExtra * newExtra = new reco::TrackExtra( trk->outerPosition(), trk->outerMomentum(),
                                        trk->outerOk(), trk->innerPosition(),
                                        trk->innerMomentum(), trk->innerOk(),
                                        trk->outerStateCovariance(), trk->outerDetId(),
                                        trk->innerStateCovariance(), trk->innerDetId() , seedDir ) ;*/

      // new copy of the silicon hits; add hit refs to Extra and hits to hit collection
      //unsigned int index_hit = 0;
      
      // edm::LogVerbatim("MuonTrackProducer")<<"\n printing initial hit_pattern";
      // trk->hitPattern().print();

      /*for (trackingRecHit_iterator iHit = trk->recHitsBegin(); iHit != trk->recHitsEnd(); iHit++) {
        TrackingRecHit* hit = (*iHit)->clone();
        index_hit++;
        selectedTrackHits->push_back( hit );

        newExtra->add( TrackingRecHitRef( rHits, hidx++ ) );
      }*/

      // edm::LogVerbatim("MuonTrackProducer")<<"\n printing final hit_pattern";
      // newTrk->hitPattern().print();
      
      selectedTracks->push_back( *newTrk );
      //selectedTrackExtras->push_back( *newExtra );

    } // if (isGoodResult)
  } // loop on reco::MuonCollection
  
  iEvent.put(selectedTracks);
  //iEvent.put(selectedTrackExtras);
  //iEvent.put(selectedTrackHits);
}
