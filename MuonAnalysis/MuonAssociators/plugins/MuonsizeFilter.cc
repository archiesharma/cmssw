// -*- C++ -*-
//
// Package:    MuonsizeFilter
// Class:      MuonsizeFilter
// 
/**\class MuonsizeFilter MuonsizeFilter.cc MuonAnalysis/MuonAssociators/plugins/MuonsizeFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Roumyana Mileva Hadjiiska
//         Created:  Wed, 08 Oct 2014 08:50:53 GMT
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/MuonReco/interface/Muon.h"

//
// class declaration
//

class MuonsizeFilter : public edm::EDFilter {
   public:
      explicit MuonsizeFilter(const edm::ParameterSet&);
      ~MuonsizeFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;


      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonsizeFilter::MuonsizeFilter(const edm::ParameterSet& iConfig)
{

}


MuonsizeFilter::~MuonsizeFilter()
{
}

bool
MuonsizeFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  bool result = false;   
  using namespace edm;
  using namespace std;
  edm::Handle<std::vector<reco::Muon> > myMuons;
  iEvent.getByLabel("muons", myMuons);
  unsigned int muSize = myMuons->size();
  if (muSize > 0) result = true;
  return result;
}
void 
MuonsizeFilter::beginJob()
{
}
void 
MuonsizeFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
MuonsizeFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
MuonsizeFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
MuonsizeFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
MuonsizeFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonsizeFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(MuonsizeFilter);
