/*
 *  See header file for a description of this class.
 *
 *  $Date: 2013/04/24 17:16:35 $
 *  $Revision: 1.1 $
 *  \author M. Maggi -- INFN Bari
 */



#include "RecoLocalMuon/GEMRecHit/interface/GEMRecHitBaseAlgoPreReco.h"
#include "RecoLocalMuon/GEMRecHit/src/GEMClusterContainer.h"
#include "RecoLocalMuon/GEMRecHit/src/GEMCluster.h"
#include "RecoLocalMuon/GEMRecHit/src/GEMClusterizer.h"
#include "RecoLocalMuon/GEMRecHit/src/GEMMaskReClusterizer.h"
#include "DataFormats/GEMDigi/interface/GEMDigiPreRecoCollection.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


GEMRecHitBaseAlgoPreReco::GEMRecHitBaseAlgoPreReco(const edm::ParameterSet& config) {
  //  theSync = GEMTTrigSyncFactory::get()->create(config.getParameter<string>("tTrigMode"),
  //config.getParameter<ParameterSet>("tTrigModeConfig"));
  recOnlyMuons = config.getParameter<bool>("recOnlyMuons");
}

GEMRecHitBaseAlgoPreReco::~GEMRecHitBaseAlgoPreReco(){}


// Build all hits in the range associated to the layerId, at the 1st step.
edm::OwnVector<GEMRecHit> GEMRecHitBaseAlgoPreReco::reconstruct(const GEMDetId& gemId,
                   						      const GEMDigiPreRecoCollection::Range& digiRange) {
  edm::OwnVector<GEMRecHit> result; 



  for (GEMDigiPreRecoCollection::const_iterator digi = digiRange.first;
       digi != digiRange.second;digi++) {

   std::cout << "inside the digi loop " << std::endl;

    LocalError tmpErr;
    LocalPoint point;
    bool OK = this->compute(*digi,point, tmpErr);
    std::cout << "compute the digi" << std::endl;
    if (!OK) continue;  

    if (recOnlyMuons && std::abs(digi->pdgid()) != 13)  continue;
     std::cout << "pgd id and recoOnlyMuons" << std::endl;

    if (recOnlyPrompt && !digi->prompt()) continue;
   std::cout << "digi prompt and recoOnlyPrompt" << std::endl;

    GEMRecHit* recHit = new GEMRecHit(gemId,digi->tof());
    result.push_back(recHit);
    std::cout << "rec hit result " << std::endl;
  }

  return result;
}
