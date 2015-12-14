#ifndef RecoLocalMuon_GEMRecHitAlgoFactoryPreReco_H
#define RecoLocalMuon_GEMRecHitAlgoFactoryPreReco_H

/** \class GEMRecHitAlgoFactoryPreReco
 *  Factory of seal plugins for 1D RecHit reconstruction algorithms.
 *  The plugins are concrete implementations of GEMRecHitBaseAlgo base class.
 *
 *  $Date: 2013/04/24 17:16:32 $
 *  $Revision: 1.1 $
 *  \author G. Cerminara - INFN Torino
 */
#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "RecoLocalMuon/GEMRecHit/interface/GEMRecHitBaseAlgoPreReco.h"

typedef edmplugin::PluginFactory<GEMRecHitBaseAlgoPreReco *(const edm::ParameterSet &)> GEMRecHitAlgoFactoryPreReco;
#endif




