#ifndef GEMDigitizer_GEMDigiPreRecoProducer_h
#define GEMDigitizer_GEMDigiPreRecoProducer_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimMuon/GEMDigitizer/interface/GEMDigiPreRecoModel.h"
#include "DataFormats/Common/interface/DetSetVector.h"

#include "string"

class GEMGeometry;

class GEMDigiPreRecoProducer : public edm::EDProducer
{
public:

  explicit GEMDigiPreRecoProducer(const edm::ParameterSet& ps);

  virtual ~GEMDigiPreRecoProducer();

  virtual void beginRun( edm::Run&, const edm::EventSetup& ) {};

  virtual void produce(edm::Event& e, const edm::EventSetup& c);

  virtual void endRun( edm::Run&, const edm::EventSetup& ) {}

private:

  //Name of Collection used for create the XF 
  std::string collectionXF_;
  std::string digiPreRecoModelString_;
  
  GEMDigiPreRecoModel* gemDigiPreRecoModel_;

};

#endif

