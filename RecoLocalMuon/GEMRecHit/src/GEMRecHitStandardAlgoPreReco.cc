/*
 *  See header file for a description of this class.
 *
 *  $Date: 2013/04/24 17:16:35 $
 *  $Revision: 1.1 $
 *  \author M. Maggi -- INFN
 */

#include "GEMCluster.h"
#include "RecoLocalMuon/GEMRecHit/src/GEMRecHitStandardAlgoPreReco.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/Exception.h"


GEMRecHitStandardAlgoPreReco::GEMRecHitStandardAlgoPreReco(const edm::ParameterSet& config) :
  GEMRecHitBaseAlgoPreReco(config) 
{
}

GEMRecHitStandardAlgoPreReco::~GEMRecHitStandardAlgoPreReco() 
{
}

void GEMRecHitStandardAlgoPreReco::setES(const edm::EventSetup& setup) 
{
}

// First Step

bool GEMRecHitStandardAlgoPreReco::compute(const GEMDigiPreReco& digi,
				    LocalPoint& Point,
				    LocalError& error)  const
{
  LocalPoint loctemp2(digi.x(),digi.y(),0.);
  Point = loctemp2;
  LocalError loerr2(digi.ex()*digi.ex(),digi.corr()*digi.ex()*digi.ey(),digi.ey()*digi.ey());
  error = loerr2;
  return true;
}

bool GEMRecHitStandardAlgoPreReco::compute(const GEMDigiPreReco& digi,
                                           const float& angle,
				           const GlobalPoint& globPos, 
				           LocalPoint& Point,
				           LocalError& error ) const
{
std::cout << "compute working fine " << std::endl;
  if (globPos.z()<0){ } // Fake use to avoid warnings
  if (angle<0.){ }      // Fake use to avoid warnings
  this->compute(digi,Point,error);
  return true;
}
