#ifndef GEMDigitizer_GEMPreRecoNoSmearModel_h
#define GEMDigitizer_GEMPreRecoNoSmearModel_h

/**
 * \class GEMPreRecoNoSmearModel
 *
 * Class for the GEM NoSmear response simulation as pre-reco step 
 */

#include "SimMuon/GEMDigitizer/interface/GEMDigiPreRecoModel.h"

class GEMGeometry;

class GEMPreRecoNoSmearModel: public GEMDigiPreRecoModel
{
public:

  GEMPreRecoNoSmearModel(const edm::ParameterSet&);

  ~GEMPreRecoNoSmearModel() {}

  void simulateSignal(const GEMEtaPartition*, const edm::PSimHitContainer&);

  void setRandomEngine(CLHEP::HepRandomEngine&) {}

//  void simulateNoise(const GEMEtaPartition*);   // to be added later

  void setup() {}

private:
};
#endif
