#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Common/interface/Handle.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"

#include "SimMuon/GEMDigitizer/interface/GEMDigiPreRecoProducer.h"
#include "SimMuon/GEMDigitizer/interface/GEMDigiPreRecoModelFactory.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"

#include "CLHEP/Random/RandomEngine.h"

#include <sstream>
#include <string>
#include <map>
#include <vector>


GEMDigiPreRecoProducer::GEMDigiPreRecoProducer(const edm::ParameterSet& ps)
  : collectionXF_(ps.getParameter<std::string>("inputCollection"))
  , digiPreRecoModelString_(ps.getParameter<std::string>("digiPreRecoModelString"))
{
  produces<GEMDigiPreRecoCollection>();

  edm::Service<edm::RandomNumberGenerator> rng;
  if (!rng.isAvailable()){
    throw cms::Exception("Configuration")
      << "GEMDigiPreRecoProducer::GEMPreRecoDigiProducer() - RandomNumberGeneratorService is not present in configuration file.\n"
      << "Add the service in the configuration file or remove the modules that require it.";
  }
  CLHEP::HepRandomEngine& engine = rng->getEngine();
  gemDigiPreRecoModel_ = GEMDigiPreRecoModelFactory::get()->create("GEM" + digiPreRecoModelString_ + "Model", ps);
  LogDebug("GEMDigiPreRecoProducer") << "Using GEM" + digiPreRecoModelString_ + "Model";
  gemDigiPreRecoModel_->setRandomEngine(engine);
}


GEMDigiPreRecoProducer::~GEMDigiPreRecoProducer()
{
  delete gemDigiPreRecoModel_;
}


void GEMDigiPreRecoProducer::produce(edm::Event& e, const edm::EventSetup& eventSetup)
{
  // set geometry
  edm::ESHandle<GEMGeometry> hGeom;
  eventSetup.get<MuonGeometryRecord>().get(hGeom);
  gemDigiPreRecoModel_->setGeometry(&*hGeom);
  gemDigiPreRecoModel_->setup();

  edm::Handle<CrossingFrame<PSimHit> > cf;
  e.getByLabel("mix", collectionXF_, cf);

  std::auto_ptr<MixCollection<PSimHit> > hits( new MixCollection<PSimHit>(cf.product()) );

  // Create empty output
  std::auto_ptr<GEMDigiPreRecoCollection> digis(new GEMDigiPreRecoCollection());

  // arrange the hits by eta partition
  std::map<uint32_t, edm::PSimHitContainer> hitMap;
  for(auto &hit: *hits){
    hitMap[hit.detUnitId()].push_back(hit);
  }
  
  // simulate signal and noise for each eta partition
  const auto & etaPartitions(gemDigiPreRecoModel_->getGeometry()->etaPartitions());
  
  for(auto &roll: etaPartitions){
    const GEMDetId detId(roll->id());
    const uint32_t rawId(detId.rawId());
    const auto & simHits(hitMap[rawId]);
    
    LogDebug("GEMDigiPreRecoProducer") 
      << "GEMDigiPreRecoProducer: found " << simHits.size() << " hit(s) in eta partition" << rawId;
    
    gemDigiPreRecoModel_->simulateSignal(roll, simHits);
    //gemDigiPreRecoModel_->simulateNoise(roll); //to be added later
    gemDigiPreRecoModel_->fillDigis(rawId, *digis);
  }
  
  // store them in the event
  e.put(digis);
}

