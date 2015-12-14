/** \file
 *
 *  $Date: 2013/04/24 17:16:35 $
 *  $Revision: 1.1 $
 *  \author M. Maggi -- INFN Bari
*/

#include "GEMRecHitProducerPreReco.h"


#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMDigiPreRecoCollection.h"

#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHit.h"

#include "RecoLocalMuon/GEMRecHit/interface/GEMRecHitBaseAlgoPreReco.h"
#include "RecoLocalMuon/GEMRecHit/interface/GEMRecHitAlgoFactoryPreReco.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"

// #include "CondFormats/GEMObjects/interface/GEMMaskedStrips.h"
// #include "CondFormats/DataRecord/interface/GEMMaskedStripsRcd.h"
// #include "CondFormats/GEMObjects/interface/GEMDeadStrips.h"
// #include "CondFormats/DataRecord/interface/GEMDeadStripsRcd.h"

#include <string>


using namespace edm;
using namespace std;


GEMRecHitProducerPreReco::GEMRecHitProducerPreReco(const ParameterSet& config){

  // Set verbose output

  produces<GEMRecHitCollection>();

  thePreGEMDigiLabel = config.getParameter<InputTag>("gemPreDigiLabel");

  // Get the concrete reconstruction algo from the factory

  string theAlgoName = config.getParameter<string>("prerecAlgo");
  theAlgo = GEMRecHitAlgoFactoryPreReco::get()->create(theAlgoName,
						config.getParameter<ParameterSet>("recAlgoConfig"));

}


GEMRecHitProducerPreReco::~GEMRecHitProducerPreReco(){

  delete theAlgo;
  // delete GEMMaskedStripsObj;
  // delete GEMDeadStripsObj;

}



void GEMRecHitProducerPreReco::beginRun(const edm::Run& r, const edm::EventSetup& setup){

}



void GEMRecHitProducerPreReco::produce(Event& event, const EventSetup& setup) {

  // Get the GEM Geometry

  ESHandle<GEMGeometry> gemGeom;
  setup.get<MuonGeometryRecord>().get(gemGeom);

  // Get the digis from the event

  Handle<GEMDigiPreRecoCollection> predigis;
  event.getByLabel(thePreGEMDigiLabel,predigis);

  // Pass the EventSetup to the algo

  theAlgo->setES(setup);

  // Create the pointer to the collection which will store the rechits

  auto_ptr<GEMRecHitCollection> recHitCollection(new GEMRecHitCollection());

  // Iterate through all digi collections ordered by LayerId   

  GEMDigiPreRecoCollection::DigiRangeIterator gemdgIt;
  for (gemdgIt = predigis->begin(); gemdgIt != predigis->end();
       ++gemdgIt){
       
    // The layerId
    const GEMDetId& gemId = (*gemdgIt).first;

   // Get the iterators over the digis associated with this LayerId
    const GEMDigiPreRecoCollection::Range& range = (*gemdgIt).second;


    // Call the reconstruction algorithm    

    OwnVector<GEMRecHit> recHits =
      theAlgo->reconstruct(gemId, range);
    
    if(recHits.size() > 0) //FIXME: is it really needed?
      recHitCollection->put(gemId, recHits.begin(), recHits.end());
  }

  event.put(recHitCollection);

}

