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
		explicit GEMSegmentStudies(const edm::ParameterSet&);
		~GEMSegmentStudies();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

		// ----------member data ---------------------------
		edm::ESHandle<GEMGeometry> gemGeom;

		edm::EDGetTokenT<GEMSegmentCollection>        GEMSegment_Token;
		edm::EDGetTokenT<edm::SimTrackContainer>      SIMTrack_Token;
		edm::EDGetTokenT<edm::PSimHitContainer> SIMHitsGEM;
		edm::EDGetTokenT<GEMDigiCollection> gemDigiToken_;
		edm::EDGetTokenT<edm::DetSetVector<GEMDigiSimLink> > gemDigiSimLinkToken_;
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

{
	//now do what ever initialization is needed
	usesResource("TFileService");
	GEMSegment_Token  = consumes<GEMSegmentCollection>(edm::InputTag("gemSegments"));
	SIMTrack_Token    = consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits"));
	SIMHitsGEM  = consumes<edm::PSimHitContainer>(edm::InputTag("g4SimHits","MuonGEMHits"));
	gemDigiToken_ = consumes<GEMDigiCollection>(edm::InputTag("simMuonGEMDigis"));
	gemDigiSimLinkToken_ = consumes<edm::DetSetVector<GEMDigiSimLink> >(edm::InputTag("simMuonGEMDigis","GEM")); 
}


GEMSegmentStudies::~GEMSegmentStudies()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

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
	// Loop over GEM Segments
	//   // ======================
	for (auto gems = gemSegmentCollection->begin(); gems != gemSegmentCollection->end(); ++gems) {
		auto gemrhs = gems->specificRecHits();
		std::cout<<"gemrhs::"<<gemrhs.size()<<std::endl;

	}
#ifdef THIS_IS_AN_EVENT_EXAMPLE
	Handle<ExampleData> pIn;
	iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	ESHandle<SetupData> pSetup;
	iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
	void 
GEMSegmentStudies::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
GEMSegmentStudies::endJob() 
{
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
