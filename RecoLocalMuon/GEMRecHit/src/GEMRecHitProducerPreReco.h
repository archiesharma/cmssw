#ifndef RecoLocalMuon_GEMRecHitProducerPreReco_h
#define RecoLocalMuon_GEMRecHitProducerPreReco_h

/** \class GEMRecHitProducerPreReco
 *  Module for GEMRecHit production. 
 *  
 *  $Date: 2013/06/10 15:53:23 $
 *  $Revision: 1.2 $
 *  \author M. Maggim -- INFN Bari
 */


#include <memory>
#include <fstream>
#include <iostream>
#include <stdint.h>
#include <cstdlib>
#include <bitset>
#include <map>

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"

// #include "CondFormats/GEMObjects/interface/GEMMaskedStrips.h"
// #include "CondFormats/DataRecord/interface/GEMMaskedStripsRcd.h"
// #include "CondFormats/GEMObjects/interface/GEMDeadStrips.h"
// #include "CondFormats/DataRecord/interface/GEMDeadStripsRcd.h"



namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

class GEMRecHitBaseAlgoPreReco;

class GEMRecHitProducerPreReco : public edm::EDProducer {

public:
  /// Constructor
  GEMRecHitProducerPreReco(const edm::ParameterSet& config);

  /// Destructor
  virtual ~GEMRecHitProducerPreReco();

  // Method that access the EventSetup for each run
  virtual void beginRun(const edm::Run&, const edm::EventSetup& ) override;

  /// The method which produces the rechits
  virtual void produce(edm::Event& event, const edm::EventSetup& setup) override;

private:

  // The label to be used to retrieve GEM digis from the event
  edm::InputTag thePreGEMDigiLabel;
  // The reconstruction algorithm
  GEMRecHitBaseAlgoPreReco *theAlgo;
  //   static std::string theAlgoName;


};

#endif

