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
#include "TMath.h"

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
    gaussianSmearing_(config.getParameter<bool>("gaussianSmearing")),
    averageEfficiency_(config.getParameter<double>("averageEfficiency")), 
    // simulateIntrinsicNoise_(config.getParameter<bool>("simulateIntrinsicNoise")),
    // averageNoiseRate_(config.getParameter<double>("averageNoiseRate")), 
    simulateElectronBkg_(config.getParameter<bool>("simulateElectronBkg")), 
    simulateNeutralBkg_(config.getParameter<bool>("simulateNeutralBkg")), 
    doBkgNoise_(config.getParameter<bool>("doBkgNoise")),
    minBunch_(config.getParameter<int>("minBunch")), 
    maxBunch_(config.getParameter<int>("maxBunch"))
{

  //initialise parameters from the fit:
  //params for pol3 model of electron bkg for GE1/1:
    GE11ElecBkgParam0 = 3402.63;
    GE11ElecBkgParam1 = -42.9979;
    GE11ElecBkgParam2 = 0.188475;
    GE11ElecBkgParam3 = -0.0002822;
  //params for expo model of electron bkg for GE2/1:
    constElecGE21 = 9.74156e+02;
    slopeElecGE21 = -1.18398e-02;

  //Neutral Bkg
  //High Rate model L=5x10^{34}cm^{-2}s^{-1}
  //params for expo model of neutral bkg for GE1/1:
    constNeuGE11_highRate = 1.02603e+04;
    slopeNeuGE11_highRate = -1.62806e-02;
  //params for pol5 model of neutral bkg for GE2/1:
    GE21ModNeuBkgParam0 = 21583.2;
    GE21ModNeuBkgParam1 = -476.59;
    GE21ModNeuBkgParam2 = 4.24037;
    GE21ModNeuBkgParam3 = -0.0185558;
    GE21ModNeuBkgParam4 = 3.97809e-05;
    GE21ModNeuBkgParam5 = -3.34575e-08;
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
    // Digitize only in [minBunch,maxBunch] window
    // window is: [(2n-1)*bxw/2, (2n+1)*bxw/2], n = [minBunch, maxBunch]
    if(hit.timeOfFlight() < (2*minBunch_-1)*bxwidth*1.0/2 || hit.timeOfFlight() > (2*maxBunch_+1)*bxwidth*1.0/2) continue;
    // is GEM efficient?
    if (flat1_->fire(1) > averageEfficiency_) continue;
    // create digi
    auto entry = hit.entryPoint();
    double x=0.0, y=0.0;
    if(gaussianSmearing_) { // Gaussian Smearing
      x=gauss_->fire(entry.x(), sigma_u);
      y=gauss_->fire(entry.x(), sigma_u);
    }
    else { // Uniform Smearing ... use the sigmas as boundaries
      x=entry.x()+(flat1_->fire(0., 1.)-0.5)*sigma_u;
      y=entry.y()+(flat1_->fire(0., 1.)-0.5)*sigma_v;
    }
    float ex = sigma_u;
    float ey = sigma_v;
    float corr = 0.;
    float tof = gauss_->fire(hit.timeOfFlight(), sigma_t);
    int pdgid = hit.particleType();
    GEMDigiPreReco digi(x, y, ex, ey, corr, tof, pdgid, 1);
    digi_.insert(digi);
  }
}


void GEMPreRecoGaussianModel::simulateNoise(const GEMEtaPartition* roll)
{
 if (!doBkgNoise_)
    return;
 const GEMDetId gemId(roll->id()); 
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
     double averageElectronRatePerRoll = 0.0;

       if (gemId.station() == 1) {
          averageElectronRatePerRoll = GE11ElecBkgParam0 + GE11ElecBkgParam1 * rollRadius
          + GE11ElecBkgParam2 * rollRadius * rollRadius + GE11ElecBkgParam3 * rollRadius * rollRadius * rollRadius;
       }
       if (gemId.station() == 2 || gemId.station() == 3){
       averageElectronRatePerRoll = constElecGE21 * TMath::Exp(slopeElecGE21 * rollRadius);
       }

      // Rate [Hz/cm^2] * 25*10^-9 [s] * Area [cm] = # hits in this roll 
      const double averageElecRate(averageElectronRatePerRoll * (bxwidth*1.0e-9) * trArea);
      int n_elechits(poisson_->fire(averageElecRate));

      float myRandY = flat2_->fire(0., 1.);
      float yy_rand = height * (myRandY - 0.5); // random Y coord in Local Coords
      //double yy_glob = rollRadius + yy_rand;    // random Y coord in Global Coords

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
	GEMDigiPreReco digi(xx_rand, yy_rand, ex, ey, corr, time, pdgid, 0);
	digi_.insert(digi);

        //edm::LogVerbatim("GEMPreRecoGaussianModelNoise") << "[GEMPreRecoDigi :: simulateNoise :: ele bkg] :: electron hit in "<<roll->id()<<" pdgid = "<<pdgid<<" bx = "<<bx<<" ==> digitized at loc x = "<<xx_rand<<" loc y = "<<yy_rand<<" time = "<<time<<" [ns]";

      }  ///////for (int i = 0; i < n_elechits; ++i)
    }  ////if (simulateElectronBkg_)


/*
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
*/
  } // end loop over bx
}  // end loop on void class
