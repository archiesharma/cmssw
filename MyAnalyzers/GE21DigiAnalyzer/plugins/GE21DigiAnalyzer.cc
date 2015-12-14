// -*- C++ -*-
//
// Package:    GE21DigiAnalyzer
// Class:      GE21DigiAnalyzer
// 
/**\class GE21DigiAnalyzer GE21DigiAnalyzer.cc MyAnalyzers/GE21DigiAnalyzer/plugins/GE21DigiAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Archie Sharma
//         Created:  Wed, 25 Nov 2015 09:51:20 GMT
// $Id$
//
//


// system include files
#include <memory>
#include <fstream>
#include <sys/time.h>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

// root include files
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"

//user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonTime.h"
#include "DataFormats/GEMDigi/interface/GEMDigiPreRecoCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"


//
// class declaration
//

using namespace reco;
using namespace std;
class GE21DigiAnalyzer : public edm::EDAnalyzer {
   public:
      explicit GE21DigiAnalyzer(const edm::ParameterSet&);
      ~GE21DigiAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

      edm::Service<TFileService> fs;

     TH1F *hFillTimingDigi, *hFillrPosDigi;

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
GE21DigiAnalyzer::GE21DigiAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


GE21DigiAnalyzer::~GE21DigiAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GE21DigiAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{	
 using namespace edm; 
// Access ME0Digis
// -----------------------
edm::Handle<GEMDigiPreRecoCollection> gemdigis;
iEvent.getByLabel("simMuonGEMDigiswithTiming", gemdigis);

edm::ESHandle<GEMGeometry> gemGeom;
iSetup.get<MuonGeometryRecord>().get(gemGeom);

GEMDigiPreRecoCollection::DigiRangeIterator detUnitIt;
// Loop over the DetUnits

for (detUnitIt = gemdigis->begin(); detUnitIt != gemdigis->end(); ++detUnitIt) {
    const GEMDetId& id = (*detUnitIt).first;
    const GEMEtaPartition* roll = gemGeom->etaPartition(id);

// Loop over the digis of this DetUnit

const GEMDigiPreRecoCollection::Range& range = (*detUnitIt).second;
    for (GEMDigiPreRecoCollection::const_iterator digiIt = range.first; digiIt!=range.second; ++digiIt) {
      double     hitTime = digiIt->tof();
      int          hitId = digiIt->pdgid();
      std::cout <<"pdg Id of the digis" <<  hitId << std::endl;
      LocalPoint      lp(digiIt->x(),digiIt->y(),0);
      GlobalPoint     gp = roll->toGlobal(lp);
      double        rpos = sqrt(pow(gp.x(),2)+pow(gp.y(),2));
      if(fabs(hitId)==13)        { hFillrPosDigi->Fill(rpos);   hFillTimingDigi->Fill(hitTime); }  
      
    }
   } 

} //main analyzer loop


// ------------ method called once each job just before starting event loop  ------------
void 
GE21DigiAnalyzer::beginJob()
{

hFillTimingDigi = fs->make<TH1F>("hFillTimingDigi","hFillTimingDigi",5000,-350,150);
hFillrPosDigi = fs->make<TH1F>("hFillrPosDigi","hFillrPosDigi",300,50,200); 

}

// ------------ method called once each job just after ending the event loop  ------------
void 
GE21DigiAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
GE21DigiAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
GE21DigiAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
GE21DigiAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
GE21DigiAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GE21DigiAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GE21DigiAnalyzer);
