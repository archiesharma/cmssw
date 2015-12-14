#ifndef GEMDigitizer_GEMDigiPreRecoModelFactory_h
#define GEMDigitizer_GEMDigiPreRecoModelFactory_h

/** \class GEMDigiPreRecoModelFactory
 *
 * Factory of seal plugins for GEMPreRecoDigitizer
 */
#include "FWCore/PluginManager/interface/PluginFactory.h"

namespace edm{
  class ParameterSet;
}

class GEMDigiPreRecoModel;

typedef edmplugin::PluginFactory<GEMDigiPreRecoModel *(const edm::ParameterSet &)> GEMDigiPreRecoModelFactory;

#endif
