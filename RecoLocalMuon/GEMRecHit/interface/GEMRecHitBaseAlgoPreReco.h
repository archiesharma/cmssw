#ifndef RecoLocalMuon_GEMRecHitBaseAlgoPreReco_H
#define RecoLocalMuon_GEMRecHitBaseAlgoPreReco_H

/** \class GEMRecHitBaseAlgoPreReco
 *  Abstract algorithmic class to compute Rec Hit
 *  form a GEM digi
 *
 *  $Date: 2013/04/24 17:16:32 $
 *  $Revision: 1.1 $
 *  \author M. Maggi -- INFN Bari
 */


#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometrySurface/interface/LocalError.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMDigiPreRecoCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHit.h"
#include "DataFormats/Common/interface/OwnVector.h"

#include "RecoLocalMuon/GEMRecHit/src/GEMEtaPartitionMask.h"
#include "RecoLocalMuon/GEMRecHit/src/GEMMaskReClusterizer.h"

class GEMDetId;

namespace edm {
  class ParameterSet;
  class EventSetup;
}


class GEMRecHitBaseAlgoPreReco {

 public:
  
  /// Constructor
  GEMRecHitBaseAlgoPreReco(const edm::ParameterSet& config);

  /// Destructor
  virtual ~GEMRecHitBaseAlgoPreReco();  

  /// Pass the Event Setup to the algo at each event
  virtual void setES(const edm::EventSetup& setup) = 0;

  /// Build all hits in the range associated to the gemId, at the 1st step.
  virtual edm::OwnVector<GEMRecHit> reconstruct(const GEMDetId& gemId,
						const GEMDigiPreRecoCollection::Range& digiRange );


  /// standard local recHit computation

  virtual bool compute(const GEMDigiPreReco& digi,
                       LocalPoint& Point,
                       LocalError& error) const = 0; 
  virtual bool compute(const GEMDigiPreReco& digi,
                       const float& angle,
                       const GlobalPoint& globPos, 
                       LocalPoint& Point,
                       LocalError& error) const = 0;

  bool recOnlyMuons;
  bool recOnlyPrompt;

};
#endif
