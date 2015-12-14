#ifndef RecoLocalMuon_GEMRecHitStandardAlgoPreReco_H
#define RecoLocalMuon_GEMRecHitStandardAlgoPreReco_H

/** \class GEMRecHitStandardAlgoPreReco
 *  Concrete implementation of GEMRecHitBaseAlgo.
 *
 *  $Date: 2013/04/24 17:16:36 $
 *  $Revision: 1.1 $
 *  \author M. Maggi -- INFN Bari
 */

#include "RecoLocalMuon/GEMRecHit/interface/GEMRecHitBaseAlgoPreReco.h"



class GEMRecHitStandardAlgoPreReco : public GEMRecHitBaseAlgoPreReco {
 public:
  /// Constructor
  GEMRecHitStandardAlgoPreReco(const edm::ParameterSet& config);

  /// Destructor
  virtual ~GEMRecHitStandardAlgoPreReco();

  // Operations

  /// Pass the Event Setup to the algo at each event
  virtual void setES(const edm::EventSetup& setup);


  virtual bool compute(const GEMDigiPreReco& digi,
                       LocalPoint& point,
                       LocalError& error) const;

  virtual bool compute(const GEMDigiPreReco& digi,
                       const float& angle,
                       const GlobalPoint& globPos, 
                       LocalPoint& point,
                       LocalError& error) const;
};
#endif


