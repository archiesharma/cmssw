#include "SimMuon/GEMDigitizer/interface/GEMPreRecoGaussianModel.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "DataFormats/GEMDigi/interface/GEMDigiPreReco.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include <cmath>
#include <utility>
#include <map>


const double cspeed = 29.9792458; // [cm/ns]
const int bxwidth   = 25;         // [ns]

GEMPreRecoGaussianModel::GEMPreRecoGaussianModel(const edm::ParameterSet& config) :
    GEMDigiPreRecoModel(config), 
    sigma_t(config.getParameter<double>("timeResolution")), 
    sigma_u(config.getParameter<double>("phiResolution")), 
    sigma_v(config.getParameter<double>("etaResolution")), 
    corr(config.getParameter<bool>("useCorrelation")), 
    etaproj(config.getParameter<bool>("useEtaProjectiveGEO")), 
    digitizeOnlyMuons_(config.getParameter<bool>("digitizeOnlyMuons")), 
    averageEfficiency_(config.getParameter<double>("averageEfficiency")), 
    // simulateIntrinsicNoise_(config.getParameter<bool>("simulateIntrinsicNoise")),
    // averageNoiseRate_(config.getParameter<double>("averageNoiseRate")), 
    simulateElectronBkg_(config.getParameter<bool>("simulateElectronBkg")), 
    simulateNeutralBkg_(config.getParameter<bool>("simulateNeutralBkg")), 
    minBunch_(config.getParameter<int>("minBunch")), 
    maxBunch_(config.getParameter<int>("maxBunch"))
{
  // polynomial parametrisation of neutral (n+g) and electron background
  // neuBkg.push_back(899644.0);     neuBkg.push_back(-30841.0);     neuBkg.push_back(441.28); 
  // neuBkg.push_back(-3.3405);      neuBkg.push_back(0.0140588);    neuBkg.push_back(-3.11473e-05); neuBkg.push_back(2.83736e-08);
  // eleBkg.push_back(4.68590e+05);  eleBkg.push_back(-1.63834e+04); eleBkg.push_back(2.35700e+02);
  // eleBkg.push_back(-1.77706e+00); eleBkg.push_back(7.39960e-03);  eleBkg.push_back(-1.61448e-05); eleBkg.push_back(1.44368e-08);
  neuBkg.push_back(5.69e+06);     neuBkg.push_back(-293334);     neuBkg.push_back(6279.6);
  neuBkg.push_back(-71.2928);     neuBkg.push_back(0.452244);    neuBkg.push_back(-0.0015191);  neuBkg.push_back(2.1106e-06);
  eleBkg.push_back(3.77712e+06);  eleBkg.push_back(-199280);     eleBkg.push_back(4340.69);
  eleBkg.push_back(-49.922);      eleBkg.push_back(0.319699);    eleBkg.push_back(-0.00108113); eleBkg.push_back(1.50889e-06);

}
GEMPreRecoGaussianModel::~GEMPreRecoGaussianModel()
{
  if (flat1_)   delete flat1_;
  if (flat2_)   delete flat2_;
  if (gauss_)   delete gauss_;
  if (poisson_) delete poisson_;
}
void GEMPreRecoGaussianModel::setRandomEngine(CLHEP::HepRandomEngine& eng)
{
  flat1_ = new CLHEP::RandFlat(eng);
  flat2_ = new CLHEP::RandFlat(eng);
  gauss_ = new CLHEP::RandGaussQ(eng);
  poisson_ = new CLHEP::RandFlat(eng);
}
void GEMPreRecoGaussianModel::simulateSignal(const GEMEtaPartition* roll, const edm::PSimHitContainer& simHits)
{
  for (const auto & hit : simHits)
  {
    // Digitize only Muons?
    if (std::abs(hit.particleType()) != 13 && digitizeOnlyMuons_) continue;
    // is GEM efficient?
    if (flat1_->fire(1) > averageEfficiency_) continue;
    // create digi
    auto entry = hit.entryPoint();
    float x = gauss_->fire(entry.x(), sigma_u);
    float y = gauss_->fire(entry.y(), sigma_v);
    float ex = sigma_u;
    float ey = sigma_v;
    float corr = 0.;
    float tof = gauss_->fire(hit.timeOfFlight(), sigma_t);
    int pdgid = hit.particleType();
    GEMDigiPreReco digi(x, y, ex, ey, corr, tof, pdgid);
    digi_.insert(digi);
  }
}

/*
void GEMPreRecoGaussianModel::simulateNoise(const GEMEtaPartition* roll)
{
  double trArea(0.0);

  // Extract detailed information from the Strip Topology:
  // base_bottom, base_top, height, strips, pads 
  // note that (0,0) is in the middle of the roll ==> all param are at all half length
  const TrapezoidalStripTopology* top_(dynamic_cast<const TrapezoidalStripTopology*>(&(roll->topology())));

  auto& parameters(roll->specs()->parameters());
  float bottomLength(parameters[0]); bottomLength = 2*bottomLength;
  float topLength(parameters[1]);    topLength    = 2*topLength;
  float height(parameters[2]);       height       = 2*height;
  float myTanPhi    = (topLength - bottomLength) / (height * 2);
  double rollRadius = top_->radius();
  trArea = height * (topLength + bottomLength) / 2.0;

  // simulate intrinsic noise and background hits in all BX that are being read out
  for(int bx=minBunch_; bx<maxBunch_+1; ++bx) {

    // 1) Intrinsic Noise ... Not implemented right now
    // ------------------------------------------------
    // if (simulateIntrinsicNoise_)
    // {
    // }

    // 2) Background Noise 
    // ----------------------------
 
    // 2a) electron background
    // -----------------------
    if (simulateElectronBkg_) {

      float myRandY = flat2_->fire(0., 1.);
      float yy_rand = height * (myRandY - 0.5); // random Y coord in Local Coords
      double yy_glob = rollRadius + yy_rand;    // random Y coord in Global Coords

      // Extract / Calculate the Average Electron Rate 
      // for the given global Y coord from Parametrization
      double averageElectronRatePerRoll = 0.0;
      for(int j=0; j<7; ++j) { averageElectronRatePerRoll += eleBkg[j]*pow(yy_glob,j); }

      // Rate [Hz/cm^2] * 25*10^-9 [s] * Area [cm] = # hits in this roll 
      const double averageElecRate(averageElectronRatePerRoll * (bxwidth*1.0e-9) * trArea); 
      int n_elechits(poisson_->fire(averageElecRate));

      // max length in x for given y coordinate (cfr trapezoidal eta partition)
      double xMax = topLength/2.0 - (height/2.0 - yy_rand) * myTanPhi;

      // loop over amount of electron hits in this roll
      for (int i = 0; i < n_elechits; ++i) {
	//calculate xx_rand at a given yy_rand
	float myRandX = flat1_->fire(0., 1.);
	float xx_rand = 2 * xMax * (myRandX - 0.5);
	float ex = sigma_u;
	float ey = sigma_v;
	float corr = 0.;
	// extract random time in this BX
	float myrandT = flat1_->fire(0., 1.);
	float minBXtime = (bx-0.5)*bxwidth;      // float maxBXtime = (bx+0.5)*bxwidth;
	float time = myrandT*bxwidth+minBXtime;
	float myrandP = flat1_->fire(0., 1.);
	int pdgid = 0;
	if (myrandP <= 0.5) pdgid = -11; // electron
	else 	            pdgid = 11;  // positron
	GEMDigiPreReco digi(xx_rand, yy_rand, ex, ey, corr, time, pdgid);
	digi_.insert(digi);
      }
    }

    // 2b) neutral (n+g) background
    // ----------------------------
    if (simulateNeutralBkg_) {

      float myRandY = flat2_->fire(0., 1.);
      float yy_rand = height * (myRandY - 0.5); // random Y coord in Local Coords
      double yy_glob = rollRadius + yy_rand;    // random Y coord in Global Coords

      // Extract / Calculate the Average Electron Rate 
      // for the given global Y coord from Parametrization
      double averageNeutralRatePerRoll = 0.0;
      for(int j=0; j<7; ++j) { averageNeutralRatePerRoll += neuBkg[j]*pow(yy_glob,j); }

      // Rate [Hz/cm^2] * 25*10^-9 [s] * Area [cm] = # hits in this roll
      const double averageNeutrRate(averageNeutralRatePerRoll * (bxwidth*1.0e-9) * trArea);
      int n_hits(poisson_->fire(averageNeutrRate));

      // max length in x for given y coordinate (cfr trapezoidal eta partition)
      double xMax = topLength/2.0 - (height/2.0 - yy_rand) * myTanPhi;

      // loop over amount of neutral hits in this roll
      for (int i = 0; i < n_hits; ++i) {
	//calculate xx_rand at a given yy_rand
	float myRandX = flat1_->fire(0., 1.);
	float xx_rand = 2 * xMax * (myRandX - 0.5);
	float ex = sigma_u;
	float ey = sigma_v;
	float corr = 0.;
	// extract random time in this BX
        float myrandT = flat1_->fire(0., 1.);
        float minBXtime = (bx-0.5)*bxwidth;
	float time = myrandT*bxwidth+minBXtime;
	int pdgid = 0;
        float myrandP = flat1_->fire(0., 1.);
        if (myrandP <= 0.08) pdgid = 2112; // neutrons: GEM sensitivity for neutrons: 0.08%
        else                 pdgid = 22;   // photons:  GEM sensitivity for photons:  1.04% ==> neutron fraction = (0.08 / 1.04) = 0.077 = 0.08
        GEMDigiPreReco digi(xx_rand, yy_rand, ex, ey, corr, time, pdgid);
        digi_.insert(digi);
      }
    }

  } // end loop over bx
}
*/
