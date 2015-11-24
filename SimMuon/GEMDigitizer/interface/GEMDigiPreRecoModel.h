#ifndef GEMDigitizer_GEMDigiPreRecoModel_h
#define GEMDigitizer_GEMDigiPreRecoModel_h

/** 
 *  \class GEMDigiPreRecoModel
 *
 *  Base Class for the GEM strip response simulation 
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/GEMDigi/interface/GEMDigiPreRecoCollection.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "CLHEP/Random/RandomEngine.h"

#include <map>
#include <set>

class GEMEtaPartition;
class GEMGeometry;
class PSimHit;

class GEMDigiPreRecoModel
{
public:
  
  virtual ~GEMDigiPreRecoModel() {}

  void setGeometry(const GEMGeometry *geom) {geometry_ = geom;}

  const GEMGeometry* getGeometry() {return geometry_;}

  virtual void simulateSignal(const GEMEtaPartition*, const edm::PSimHitContainer&) = 0;

//  virtual void simulateNoise(const GEMEtaPartition*) = 0; // to be added later
  
  //  virtual std::vector<std::pair<int,int> > 
  //   simulateClustering(const ME0EtaPartition*, const PSimHit*, const int) = 0;

  virtual void setRandomEngine(CLHEP::HepRandomEngine&) = 0;

  void fillDigis(int rollDetId, GEMDigiPreRecoCollection&);

  virtual void setup() = 0;

protected:
  std::set< GEMDigiPreReco> digi_;

  GEMDigiPreRecoModel(const edm::ParameterSet&) {}

  const GEMGeometry * geometry_;
  
  //  DetectorHitMap detectorHitMap_;
};
#endif
