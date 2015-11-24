#ifndef GEMDigitizer_GEMPreRecoGaussianModel_h
#define GEMDigitizer_GEMPreRecoGaussianModel_h
/**
 * \class GEMPreRecoGaussianModel
 *
 * Class for the GEM Gaussian response simulation as pre-reco step
 */
#include "SimMuon/GEMDigitizer/interface/GEMDigiPreRecoModel.h"
class GEMGeometry;
namespace CLHEP
{
  class HepRandomEngine;
  class RandFlat;
  class RandGaussQ;
  class RandPoissonQ;
}
class GEMPreRecoGaussianModel: public GEMDigiPreRecoModel
{
public:
  GEMPreRecoGaussianModel(const edm::ParameterSet&);
  ~GEMPreRecoGaussianModel();
  void simulateSignal(const GEMEtaPartition*, const edm::PSimHitContainer&);
  void setRandomEngine(CLHEP::HepRandomEngine&);
//  void simulateNoise(const GEMEtaPartition*);
  void setup()
  {
  }
private:
  double sigma_t;
  double sigma_u;
  double sigma_v;
  bool corr;
  bool etaproj;
  bool digitizeOnlyMuons_;
  double averageEfficiency_;
  // bool simulateIntrinsicNoise_; // not implemented
  // double averageNoiseRate_;     // not implemented
  bool simulateElectronBkg_;   // to be added later
  bool simulateNeutralBkg_;    // to be added later
  int minBunch_;
  int maxBunch_;
  CLHEP::RandGaussQ* gauss_;
  CLHEP::RandFlat* flat1_;
  CLHEP::RandFlat* flat2_;
  CLHEP::RandFlat* poisson_;
  // params for the simple pol6 model of neutral bkg for GEM:
  std::vector<double> neuBkg, eleBkg; // to be added later
};
#endif

