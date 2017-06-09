// -*- C++ -*-
//
// Package:    Validation/GEMSegmentStudies
// Class:      GEMSegmentStudies
// 
/**\class GEMSegmentStudies GEMSegmentStudies.cc Validation/GEMSegmentStudies/plugins/GEMSegmentStudies.cc
Description: [one line class summary]
Implementation:
[Notes on implementation]
*/
//
// Original Author:  kaur amandeepkalsi
//         Created:  Wed, 05 Apr 2017 11:06:48 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <memory>
#include <fstream>
#include <sys/time.h>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>


// root include files
#include "TFile.h"
#include "TH1F.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include <Geometry/GEMGeometry/interface/GEMEtaPartition.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <DataFormats/MuonDetId/interface/GEMDetId.h>
#include <DataFormats/GEMRecHit/interface/GEMSegmentCollection.h>
#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include <DataFormats/CSCRecHit/interface/CSCSegmentCollection.h>
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLink.h"
#include "SimDataFormats/GEMDigiSimLink/interface/GEMDigiSimLink.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/GEMDigiSimLink/interface/GEMDigiSimLink.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include <TMath.h>
#include <map>
#include <set>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class GEMSegmentStudies : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
	public:
		//explicit GEMSegmentStudies(const edm::ParameterSet& iConfig);
		explicit GEMSegmentStudies(const edm::ParameterSet&);
                ~GEMSegmentStudies();
                
                GEMSegmentStudies(const edm::Event& iEvent, const edm::EventSetup& iSetup);                               

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

                std::set<GEMDigiSimLink> findGEMDigiSimLink(uint32_t gemDetId, int strip, int bx) const;

	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
                virtual void initEvent(const edm::Event&, const edm::EventSetup&);
		virtual void endJob() override;

		// ----------member data ---------------------------
		edm::ESHandle<GEMGeometry> gemGeom;
                edm::Handle < edm::DetSetVector<GEMDigiSimLink> > theSimlinkDigis;

		edm::EDGetTokenT<GEMSegmentCollection>        GEMSegment_Token;
		edm::EDGetTokenT<edm::SimTrackContainer>      SIMTrack_Token;
		edm::EDGetTokenT<edm::PSimHitContainer> SIMHitsGEM;
		edm::EDGetTokenT<GEMDigiCollection> gemDigiToken_;
		edm::EDGetTokenT<edm::DetSetVector<GEMDigiSimLink> > gemDigiSimLinkToken_;
		bool isSimMatched(edm::SimTrackContainer::const_iterator simTrack, edm::PSimHitContainer::const_iterator itHit) ;
		
                typedef std::map<edm::SimTrackContainer::const_iterator,edm::PSimHitContainer> MapTypeSim;
                typedef std::map<GEMSegmentCollection::const_iterator,std::vector<GEMDigiSimLink> > MapTypeSeg;


                //std::vector<GEMSegmentStudies::SimHitIdpr> GEMSegmentStudies::associateRecHit(const TrackingRecHit & hit);
	        //typedef std::map<GEMSegmentCollection::const_iterator,std::vector<GEMRecHit>> MapTypeSim;
		//typedef std::map<edm::DetSetVector<GEMDigiSimLink>::const_iterator,edm::DetSet<GEMDigiSimLink> > MapTypeSeg;
		int countST;
		int segInNum; 

                TH1F *gem_matchedsimsegment_eta, *gem_matchedsimsegment_pt, *gem_matchedsimsegment_phi;
                TH1F *gem_matchedsegment_eta, *gem_matchedsegment_pt, *gem_matchedsegment_phi;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GEMSegmentStudies::GEMSegmentStudies(const edm::ParameterSet& iConfig)

        /*
	  GEMSegment_Token(consumes<GEMSegmentCollection>(pset.getParameter < edm::InputTag > ("gemSegments")))
	, SIMTrack_Token(consumes<edm::SimTrackContainer>(pset.getParameter < edm::InputTag > ("g4SimHits")))
	, SIMHitsGEM(consumes<edm::PSimHitContainer>(edm::InputTag("g4SimHits","MuonGEMHits")))
	, gemDigiToken_(consumes<GEMDigiCollection>(pset.getParameter < edm::InputTag > ("simMuonGEMDigis")))
	, gemDigiSimLinkToken_(consumes<edm::DetSetVector<GEMDigiSimLink> >(edm::InputTag("simMuonGEMDigis","GEM")))
        */

{       

        GEMSegment_Token  = consumes<GEMSegmentCollection>(edm::InputTag("gemSegments"));
	SIMTrack_Token    = consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits"));
	SIMHitsGEM  = consumes<edm::PSimHitContainer>(edm::InputTag("g4SimHits","MuonGEMHits"));
	gemDigiToken_ = consumes<GEMDigiCollection>(edm::InputTag("simMuonGEMDigis"));
        gemDigiSimLinkToken_ = consumes<edm::DetSetVector<GEMDigiSimLink> >(edm::InputTag("simMuonGEMDigis","GEM")); 

        // usesResource("TFileService");
        edm::Service<TFileService> fs;
        if(!fs){
           throw edm::Exception( edm::errors::Configuration,
                             "TFileService is not registered in cfg file" );
        }

        gem_matchedsimsegment_eta = fs->make<TH1F>("gem_matchedsimsegment_eta", "Matched SimSegment Eta Distribution; #eta; entries", 10, 1.5, 2.5);
        gem_matchedsimsegment_pt = fs->make<TH1F>("gem_matchedsimsegment_pt", "Matched SimSegment pT Distribution; p_{T}; entries", 100, 0.0, 100.0);
        gem_matchedsimsegment_phi = fs->make<TH1F>("gem_matchedsimsegment_phi", "Matched SimSegment Phi Distribution; #phi; entries", 18,-M_PI,+M_PI);

        gem_matchedsegment_eta = fs->make<TH1F>("gem_matchedsegment_eta", "Matched Segment Eta Distribution; #eta; entries", 10, 1.5, 2.5);
        gem_matchedsegment_pt = fs->make<TH1F>("gem_matchedsegment_pt", "Matched Segment pT Distribution; p_{T}; entries", 100, 0.0, 100.0);
        gem_matchedsegment_phi = fs->make<TH1F>("gem_matchedsegment_phi", "Matched Segment Phi Distribution; #phi; entries", 18,-M_PI,+M_PI);


        // edm::Handle < edm::DetSetVector<GEMDigiSimLink> > theSimlinkDigis;
        // iEvent.getByToken(gemDigiSimLinkToken_, theSimlinkDigis);
         

}

GEMSegmentStudies::GEMSegmentStudies(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  initEvent(iEvent,iSetup);
}

GEMSegmentStudies::~GEMSegmentStudies()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}

void GEMSegmentStudies::initEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup)
  {

      edm::Handle < edm::DetSetVector<GEMDigiSimLink> > theSimlinkDigis;
      iEvent.getByToken(gemDigiSimLinkToken_, theSimlinkDigis);  
   
 }

//
// member functions
//

// ------------ method called for each event  ------------
	void
GEMSegmentStudies::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	iSetup.get<MuonGeometryRecord>().get(gemGeom);

	edm::Handle<GEMSegmentCollection> gemSegmentCollection;
	iEvent.getByToken(GEMSegment_Token, gemSegmentCollection);
	std::cout <<"Number of GEM Segments in this event: "<<gemSegmentCollection->size()<<"\n"<<std::endl;

	// vector<PSimHit>                       "g4SimHits"                 "MuonGEMHits"     "SIM"                                           

	edm::Handle < GEMDigiCollection > digis;
	iEvent.getByToken(gemDigiToken_, digis);

	edm::Handle < edm::DetSetVector<GEMDigiSimLink> > theSimlinkDigis;
	iEvent.getByToken(gemDigiSimLinkToken_, theSimlinkDigis);

	edm::Handle<edm::PSimHitContainer> theSimHits;
	iEvent.getByToken(SIMHitsGEM,theSimHits);
 
	edm::Handle<edm::SimTrackContainer> SimTk;
	iEvent.getByToken(SIMTrack_Token,SimTk);

        MapTypeSim myMap;
        MapTypeSeg myMapSeg;

	edm::SimTrackContainer::const_iterator simTrack;
	for (simTrack = SimTk->begin(); simTrack != SimTk->end(); ++simTrack){

	//	edm::PSimHitContainer GEMSimHits;
                edm::PSimHitContainer selectedGEMSimHits;

		if (std::abs(simTrack->type()) != 13) continue;
		int count = 0;

		for (edm::PSimHitContainer::const_iterator itHit = theSimHits->begin(); itHit != theSimHits->end(); ++itHit){

			int particleType_sh = itHit->particleType();
			// int evtId_sh = itHit->eventId().event();
			int bx_sh = itHit->eventId().bunchCrossing();
			// int procType_sh = itHit->processType();
			if(!(abs(particleType_sh) == 13 && bx_sh == 0)) continue;

			if(isSimMatched(simTrack, itHit)){

				++count;
				selectedGEMSimHits.push_back(*itHit);;

			}


		}//End loop SHs
		if(selectedGEMSimHits.size() >= 2) {
                   countST++;
                   //selectedGEMSimHits.push_back(GEMSimHits);
                   myMap.insert(MapTypeSim::value_type(simTrack,selectedGEMSimHits));
               }
	}
	// Loop over GEM Segments
	//   // ======================
	
	
	for (auto gems = gemSegmentCollection->begin(); gems != gemSegmentCollection->end(); ++gems) {
		auto gemrhs = gems->specificRecHits();
		std::cout<<"gemrhs::"<<gemrhs.size()<<std::endl;
		std::vector<GEMDigiSimLink> theGEMRecHits;
		for (auto rh = gemrhs.begin(); rh!= gemrhs.end(); rh++){
			// GEMDetId gemId((*rh).geographicalId());

                        GEMDetId gemDetId = rh->gemId();
             		int bx = rh->BunchX();
			int fstrip = rh->firstClusterStrip();
			int cls = rh->clusterSize();
			//auto rhLP = rh->localPosition();
                        if(!(bx == 0)) continue;

			std::cout<<gemDetId<<"::"<<bx<<"::"<<fstrip<<"::"<<cls<<std::endl;

                        for(int i = fstrip; i < fstrip+cls; ++i) {      /////// loop over the fired strips 

                        std::cout<<"in the loop over digi strips ::"<<std::endl;

                        std::set<GEMDigiSimLink> links = findGEMDigiSimLink(gemDetId.rawId(),i,bx);

                        if (links.empty()) edm::LogInfo("GEMSegmentStudies")
                        <<"*** WARNING in GEMSegmentStudies::associateRecHit, GEMRecHit "<<*rh<<", strip "<<i<<" has no associated GEMDigiSimLink !"<<std::endl;


                        std::cout<<"matched digi found !! "<<std::endl;
                        for(std::set<GEMDigiSimLink>::iterator itlink = links.begin(); itlink != links.end(); ++itlink) {

                                   int bx_matched = itlink->getBx();
                                   int particletype_matched = itlink->getParticleType();
                                   if(!(abs(particletype_matched) == 13 && bx_matched == 0 )) continue;
			           theGEMRecHits.push_back(*itlink); 
			           myMapSeg.insert(MapTypeSeg::value_type(gems,theGEMRecHits));
          }
	}
   
   }

}
            bool isThereOneSegmentMatched = false;
            float quality = 0;
            int num_sh = 0;
            int matched_sh = 0;
            float simsegEta = 0.0;
            float simsegPt = 0.0;
            float simsegPhi = 0.0;
            float segEta = 0.0;
            float segPt = 0.0;
            float segPhi = 0.0;

             for(auto const& st : myMap) {
              num_sh = st.second.size();
              simsegEta = (*(st.first)).momentum().eta();
              simsegPt = (*(st.first)).momentum().pt();
              simsegPhi = (*(st.first)).momentum().phi();
             } 

             for(auto const& seg : myMap) {
              matched_sh = seg.second.size();
              segEta = (*(seg.first)).momentum().eta();
              segPt = (*(seg.first)).momentum().pt(); 
              segPhi = (*(seg.first)).momentum().phi();
             }           

            //if(selectedGEMSimHits.size() >= 2) int num_sh = selectedGEMSimHits.size();
            //int matched_sh = theGEMRecHits.size();
            if(num_sh != 0) quality = matched_sh/(1.0*num_sh);
            if(quality > 0) isThereOneSegmentMatched = true;

           
               if(isThereOneSegmentMatched){
               gem_matchedsimsegment_eta->Fill(std::abs(simsegEta));
               gem_matchedsimsegment_pt->Fill(simsegPt);
               gem_matchedsimsegment_phi->Fill(simsegPhi);
                    
         
               gem_matchedsegment_eta->Fill(std::abs(segEta));
               gem_matchedsegment_pt->Fill(segPt);
               gem_matchedsegment_phi->Fill(segPhi);
             

           }


/*	// gem digi sim links
	for (edm::DetSetVector<GEMDigiSimLink>::const_iterator itsimlink = theSimlinkDigis->begin();
			itsimlink != theSimlinkDigis->end(); itsimlink++)
	{
		DetSet<GEMDigiSimLink> linkHits;                                                          

		int detid = itsimlink->detId();
		for (edm::DetSet<GEMDigiSimLink>::const_iterator link_iter = itsimlink->data.begin();
				link_iter != itsimlink->data.end(); ++link_iter)
		{

			int strip = link_iter->getStrip();
			int processtype = link_iter->getProcessType();
			int particletype = link_iter->getParticleType();
			int bx_sim = link_iter->getBx();
			double partTof = link_iter->getTimeOfFlight();
			double myEnergyLoss = link_iter->getEnergyLoss();                                                             
			auto locMuonEntry = link_iter->getEntryPoint();

			// consider only muons
			if(abs(particletype) == 13 && bx_sim == 0 ) {
				std::cout<<"ENTRY!!!!!"<<std::endl;
				std::cout<<strip<<processtype<<particletype<<bx_sim<<partTof<<myEnergyLoss<<locMuonEntry<<std::endl;
				//	theGEMRecHits.push_back(*rh);			
				linkHits.push_back(*link_iter);
			}
		}
		myMapSeg.insert(MapTypeSeg::value_type(itsimlink,linkHits));
	} // digi-sim links

int num_sh_matched;
	for(auto const& seg : myMapSeg) { 
		int num_sh = seg.second.size();
		bool isThereOneSegmentMatched = false;
		num_sh_matched=0;
		for(auto const& st : myMap) {


			for(auto const& rh : seg.second) { // digi sim Hits
				for(auto const& sh : st.second) {  // rec Hits
					auto digisimStrip = rh.getStrip();
					auto firstStrp = sh.firstClusterStrip();
					auto clsSize = sh.clusterSize();
					for(int i = firstStrp; i < (firstStrp+clsSize); ++i) {
						if(int(digisimStrip) != i) continue;
						std::cout<<"Matched Stuff:"<<std::endl;
						//
						num_sh_matched++;
					}

				}
			}
		}
		if(num_sh_matched > 0 ) segInNum++;
	}

*/

#ifdef THIS_IS_AN_EVENT_EXAMPLE
	Handle<ExampleData> pIn;
	iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	ESHandle<SetupData> pSetup;
	iSetup.get<SetupRecord>().get(pSetup);
#endif
}

std::set<GEMDigiSimLink>  GEMSegmentStudies::findGEMDigiSimLink(uint32_t gemDetId, int strip, int bx) const {

   std::cout<<"in the loop of gem segment ::"<<std::endl; 
   std::set<GEMDigiSimLink> links;

   std::cout<<"digi sim link Iter ::"<<theSimlinkDigis->size()<<std::endl;

     for (edm::DetSetVector<GEMDigiSimLink>::const_iterator itlink = theSimlinkDigis->begin(); itlink != theSimlinkDigis->end(); itlink++){

       std::cout<<"in the loop of digi sim link  ::"<<std::endl;
       for(edm::DetSet<GEMDigiSimLink>::const_iterator digi_iter=itlink->data.begin();digi_iter != itlink->data.end();++digi_iter){

       std::cout<<"in the loop of the data of digi sim link  ::"<<std::endl;

       uint32_t detid = digi_iter->getDetUnitId();
       int str = digi_iter->getStrip();
       int bunchx = digi_iter->getBx();
       if(detid == gemDetId && str == strip && bunchx == bx){
         links.insert(*digi_iter);
         }
       }
    }	
 
    return links;
}	


bool GEMSegmentStudies::isSimMatched(edm::SimTrackContainer::const_iterator simTrack, edm::PSimHitContainer::const_iterator itHit)
{

	bool result = false;
	int trackId = simTrack->trackId();
	int trackId_sim = itHit->trackId();
	if(trackId == trackId_sim) result = true;
	return result;

}
/*
   std::vector<GEMSegmentStudies::SimHitIdpr> GEMSegmentStudies::associateRecHit(const TrackingRecHit & hit) const {
   std::vector<SimHitIdpr> matched;
   const TrackingRecHit * hitp = &hit;
   const GEMRecHit * gemrechit = dynamic_cast<const GEMRecHit *>(hitp);
   if (gemrechit) {
   GEMDetId gemDetId = gemrechit->gemId();
   int fstrip = gemrechit->firstClusterStrip();
   int cls = gemrechit->clusterSize();
   int bx = gemrechit->BunchX();
   DigiSimLinks::const_iterator layerLinks = theDigiSimLinks->find(gemDetId);
   if (layerLinks != theDigiSimLinks->end()) {
   for(int i = fstrip; i < (fstrip+cls); ++i) {
   for(LayerLinks::const_iterator itlink = layerLinks->begin(); itlink != layerLinks->end(); ++itlink) {
   int ch = static_cast<int>(itlink->channel());
   if(ch != i) continue;
   SimHitIdpr currentId(itlink->SimTrackId(), itlink->eventId());
   if(find(matched.begin(),matched.end(),currentId ) == matched.end())
   matched.push_back(currentId);
   }
   }
   }else edm::LogWarning("GEMSegmentStudies")
   <<"*** WARNING in GEMSegmentStudies: GEM layer "<<gemDetId<<" has no DigiSimLinks !"<<std::endl;
   } else edm::LogWarning("GEMSegmentStudies")<<"*** WARNING in GEMSegmentStudies::associateRecHit, null dynamic_cast !";
   }
   return  matched;
   }
   */
// ------------ method called once each job just before starting event loop  ------------
	void 
GEMSegmentStudies::beginJob()
{
	countST = 0;
	segInNum = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
GEMSegmentStudies::endJob() 
{
	std::cout<<"Total Segments in denominator:"<<countST<<std::endl;
	std::cout<<"Total Segments in numerator:"<<segInNum<<std::endl;
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GEMSegmentStudies::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GEMSegmentStudies);
