#include "SimMuon/GEMDigitizer/interface/GEMPreRecoNoSmearModel.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "DataFormats/GEMDigi/interface/GEMDigiPreReco.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <cmath>
#include <utility>
#include <map>


GEMPreRecoNoSmearModel::GEMPreRecoNoSmearModel(const edm::ParameterSet& config) 
  : GEMDigiPreRecoModel(config)
{
}

void 
GEMPreRecoNoSmearModel::simulateSignal(const GEMEtaPartition* roll,
				const edm::PSimHitContainer& simHits)
{

  for (const auto & hit: simHits)
  {
    if (std::abs(hit.particleType()) != 13) continue;
    auto entry = hit.entryPoint();
    float x=entry.x();
    float y=entry.y(); 
    float ex=0.001;
    float ey=0.001;
    float corr=0.;
    float t = hit.timeOfFlight();
    int pdgid=hit.particleType();
     // please keep hit time always 0 for this model
    GEMDigiPreReco digi(x,y,ex,ey,corr,t,pdgid);
    digi_.insert(digi);
  }
}

/*
void 
GEMPreRecoNoSmearModel::simulateNoise(const GEMEtaPartition* roll)
{
}
*/


