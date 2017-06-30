// -*- C++ -*-
//
// Package:    MuonNewGEMDigis
// Class:      MuonNewGEMDigis
// 
/**\class MuonNewGEMDigis MuonNewGEMDigis.cc Validation/MuonNewGEMDigis/plugins/MuonNewGEMDigis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Archie Sharma
//         Created:  Wed, 02 Mar 2016 18:50:43 GMT
// $Id$
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DQMServices/Core/interface/DQMStore.h"
//#include "DQMServices/Core/interface/MonitorElement.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/MonitorElement.h"


//Data Format
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/GeometrySurface/interface/LocalError.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
//#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

///Geometry
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"

///Log messages
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphErrors.h"

//
// class declaration
//

class MuonNewGEMDigis : public edm::EDAnalyzer {
public:
  explicit MuonNewGEMDigis(const edm::ParameterSet&);
  ~MuonNewGEMDigis();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
  virtual void beginJob() override;
  virtual void analyze(const edm::Event& e, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  // bool isSimTrackGood(const SimTrack &);
  //  bool isGEMDigiMatched(MyGEMDigi gem_dg, MyGEMSimHit gem_sh);
  
  //  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------


  

  bool debug_;
  
  //edm::Handle<edm::PSimHitContainer> GEMHits;
  //edm::Handle<GEMDigiCollection> gem_digis;
  //edm::Handle<edm::SimTrackContainer> sim_tracks;
  //edm::Handle<edm::SimVertexContainer> sim_vertices;
  edm::ESHandle<GEMGeometry> gem_geom;
  
  edm::ParameterSet cfg_;


private:
  
  edm::EDGetTokenT<edm::PSimHitContainer> gemSimHitInput_;
  edm::EDGetTokenT<edm::SimTrackContainer> simTrackInput_; 
  edm::EDGetTokenT<GEMDigiCollection> gemDigiInput_;
 
  //double simTrackMinPt_;
  //double simTrackMaxPt_;
  //double simTrackMinEta_;
  //double simTrackMaxEta_;
  //double simTrackOnlyMuon_;
  
  const GEMGeometry* gem_geometry_;
  
  //MyGEMDigi gem_digi_;
  //MyGEMSimHit gem_sh;
  //MySimTrack track_;
  
  //bool hasGEMGeometry_;

  TH1F *hProces; 
  
  //  TH2F *hNstripEtaParts;
  TH1F *hNstripEtaPart1;
  TH1F *hNstripEtaPart2;
  TH1F *hNstripEtaPart3;
  TH1F *hNstripEtaPart4;
  TH1F *hNstripEtaPart5;
  TH1F *hNstripEtaPart6;
  TH1F *hNstripEtaPart7;
  TH1F *hNstripEtaPart8;
  TH1F *hBx;
  TH2F *hRadiusEtaPartVsNdigi;
  TH2F *hRadiusEtaPartVsNdigiOvTrArea;
  TH1F *hRadiusEtaPart;
  TH1F *hdeltaXEntryPointVsCentreStrip;
  TH1F *hResidualsSimPhi;
  TH1F *hResidualsDigiPhi;
  TH1F *hResidualsSimVsDigiPhi;
  TH1F *hglobalR_Sta1_long;
  TH1F *hglobalR_Sta1_short;
  TH1F *hglobalR_Sta2_lay1;
  TH1F *hglobalR_Sta2_lay2;
  TGraphErrors *grRatePerRoll_sta1;
  TGraphErrors *grRatePerRoll_sta2; 
 
  int numbEvents;
  int ndigi1_sta1, ndigi2_sta1, ndigi3_sta1, ndigi4_sta1, ndigi5_sta1, ndigi6_sta1, ndigi7_sta1, ndigi8_sta1;
  double ndigiVsArea1_sta1, ndigiVsArea2_sta1, ndigiVsArea3_sta1, ndigiVsArea4_sta1, ndigiVsArea5_sta1, ndigiVsArea6_sta1, ndigiVsArea7_sta1, ndigiVsArea8_sta1;
  double rollRadius1_sta1, rollRadius2_sta1, rollRadius3_sta1, rollRadius4_sta1, rollRadius5_sta1, rollRadius6_sta1, rollRadius7_sta1, rollRadius8_sta1;
  
  int ndigi1_sta2, ndigi2_sta2, ndigi3_sta2, ndigi4_sta2, ndigi5_sta2, ndigi6_sta2, ndigi7_sta2, ndigi8_sta2, ndigi9_sta2, ndigi10_sta2, ndigi11_sta2, ndigi12_sta2;
  double ndigiVsArea1_sta2, ndigiVsArea2_sta2, ndigiVsArea3_sta2, ndigiVsArea4_sta2, ndigiVsArea5_sta2, ndigiVsArea6_sta2, ndigiVsArea7_sta2, ndigiVsArea8_sta2, ndigiVsArea9_sta2, ndigiVsArea10_sta2, ndigiVsArea11_sta2, ndigiVsArea12_sta2;
  double rollRadius1_sta2, rollRadius2_sta2, rollRadius3_sta2, rollRadius4_sta2, rollRadius5_sta2, rollRadius6_sta2, rollRadius7_sta2, rollRadius8_sta2, rollRadius9_sta2, rollRadius10_sta2, rollRadius11_sta2, rollRadius12_sta2;

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
MuonNewGEMDigis::MuonNewGEMDigis(const edm::ParameterSet& cfg): 
   debug_(cfg.getParameter<bool>("debug")),   
  //now do what ever initialization is needed
  
  //now do what ever initialization is needed
   //dbe = edm::Service<DQMStore>().operator->();
  
  //simTrackMinPt_ = cfg.getParameter<double>("minPt");
  //  simTrackMaxPt_ = cfg.getParameter<double>("maxPt");
  //  simTrackMinEta_ = cfg.getParameter<double>("minEta");
  //  simTrackMaxEta_ = cfg.getParameter<double>("maxEta");
  //    simTrackOnlyMuon_ = cfg.getParameter<bool>("onlyMuon");
 
	
gemSimHitInput_(consumes<edm::PSimHitContainer>(cfg.getParameter<edm::InputTag>("simInputLabel"))),
simTrackInput_(consumes<edm::SimTrackContainer>(cfg.getParameter<edm::InputTag>("simTrackLabel"))),
gemDigiInput_(consumes<GEMDigiCollection>(cfg.getParameter<edm::InputTag>("digiInputLabel")))

//hasGEMGeometry_=false;
{

   edm::Service < TFileService > fs; 

  hProces = fs->make < TH1F > ("hProces", "Process type for all the simHits", 20, 0, 20);
//  hNstripEtaParts = fs->make <TH2F> ("NstripEtaParts", "Nstrips in each EtaPartition ", 40, 0.5, 10.5, 770, 1, 770);
  hNstripEtaPart1 = fs->make <TH1F> ("NstripEtaPart1", "Nstrips in EtaPartition 1", 770, 1, 770);
  hNstripEtaPart2 = fs->make <TH1F> ("NstripEtaPart2", "Nstrips in EtaPartition 2", 770, 1, 770);
  hNstripEtaPart3 = fs->make <TH1F> ("NstripEtaPart3", "Nstrips in EtaPartition 3", 770, 1, 770);
  hNstripEtaPart4 = fs->make <TH1F> ("NstripEtaPart4", "Nstrips in EtaPartition 4", 770, 1, 770);
  hNstripEtaPart5 = fs->make <TH1F> ("NstripEtaPart5", "Nstrips in EtaPartition 5", 770, 1, 770);
  hNstripEtaPart6 = fs->make <TH1F> ("NstripEtaPart6", "Nstrips in EtaPartition 6", 770, 1, 770);
  hNstripEtaPart7 = fs->make <TH1F> ("NstripEtaPart7", "Nstrips in EtaPartition 7", 770, 1, 770);
  hNstripEtaPart8 = fs->make <TH1F> ("NstripEtaPart8", "Nstrips in EtaPartition 8", 770, 1, 770);
  hBx = fs->make <TH1F> ("hBx", "bx from digi - for all #eta partiotions", 9, -5.5, 3.5 );
  hRadiusEtaPartVsNdigi = fs->make <TH2F> ("hRadiusEtaPartVsNdigi", "Radius Eta Partition vs Ndigi", 2500, 0., 250., 200, 0., 20. );//MM 
  hRadiusEtaPartVsNdigiOvTrArea = fs->make <TH2F> ("hRadiusEtaPartVsNdigiOvTrArea", "Ndigi/TrArea vs Radius Eta Partition", 2500, 0., 250., 1000, 0., 0.1 );
  hRadiusEtaPart = fs->make <TH1F> ("hRadiusEtaPart", "Radius Eta Partition", 200, 0., 200. );
  hdeltaXEntryPointVsCentreStrip = fs->make <TH1F> ("deltaX", "delta X Residuals", 200, -10., 10. );
  hResidualsSimPhi= fs->make <TH1F> ("ResidualsSimPhi", "Global SimMuon Phi", 200, -10., 10. );
  hResidualsDigiPhi= fs->make <TH1F> ("ResidualsDigiPhi", "Global DigiMuon Phi", 200, -10., 10. );
  hResidualsSimVsDigiPhi= fs->make <TH1F> ("ResidualsSimVsDigiPhi", "Residuals (SimMuon-Digi) Phi", 50000, -0.5, 0.5 );

  grRatePerRoll_sta1 = fs->make<TGraphErrors> (8);
  grRatePerRoll_sta1->SetName("grRatePerRoll_sta1");
  grRatePerRoll_sta1->SetTitle("GEM Rate vs Roll Radius for Station1 - BKG model");

  grRatePerRoll_sta2 = fs->make<TGraphErrors> (12);
  grRatePerRoll_sta2->SetName("grRatePerRoll_sta2");
  grRatePerRoll_sta2->SetTitle("GEM Rate vs Roll Radius for Station2 - BKG model");
  

  float bins1[] = {130.2, 141.5, 152.8, 166.3, 179.7, 195.8, 211.9, 231.3, 250.8};
  int binnum1 = sizeof(bins1)/sizeof(float) - 1;

  hglobalR_Sta1_long = fs->make <TH1F> ("hglobalR_Sta1_long", "Rate from long chambers in station 1; Radius [cm]; Entries", binnum1, bins1);

  float bins2[] = {130.2, 140.4, 150.6, 162.6, 174.5, 188.6, 202.6, 219.2, 235.8};
  int binnum2 = sizeof(bins2)/sizeof(float) - 1;

  hglobalR_Sta1_short = fs->make <TH1F> ("hglobalR_Sta1_short", "Rate from short chambers in station 1; Radius [cm]; Entries", binnum2, bins2);
 

  float bins3[] = {136.5, 158.1, 183.2, 204.8, 229.9, 251.4, 276.6, 298.1, 319.7};
  int binnum3 = sizeof(bins3)/sizeof(float) - 1;

  hglobalR_Sta2_lay1 = fs->make <TH1F> ("hglobalR_Sta2_lay1", "Rate from layer 1 chambers in station 2; Radius [cm]; Entries", binnum3, bins3);

  float bins4[] = {136.5, 156.1, 179.1, 198.7, 221.8, 245.4, 272.5, 296.1, 319.7};
  int binnum4 = sizeof(bins4)/sizeof(float) - 1;

  hglobalR_Sta2_lay2 = fs->make <TH1F> ("hglobalR_Sta2_lay2", "Rate from layer 2 chambers in station 2; Radius [cm]; Entries", binnum4, bins4);


 // hglobalR_Sta2

  numbEvents = 0;

  ndigi1_sta1 = 0; ndigi2_sta1 = 0; ndigi3_sta1 = 0; ndigi4_sta1 = 0; ndigi5_sta1 = 0; ndigi6_sta1 = 0; ndigi7_sta1 = 0; ndigi8_sta1 = 0;
  ndigiVsArea1_sta1 = 0.; ndigiVsArea2_sta1 = 0.; ndigiVsArea3_sta1 = 0.; ndigiVsArea4_sta1 = 0.; ndigiVsArea5_sta1 = 0.; ndigiVsArea6_sta1 = 0.; ndigiVsArea7_sta1 = 0.; ndigiVsArea8_sta1 = 0.;
  rollRadius1_sta1 = 0.; rollRadius2_sta1 = 0.; rollRadius3_sta1 = 0.; rollRadius4_sta1 = 0.; rollRadius5_sta1 = 0.; rollRadius6_sta1 = 0.; rollRadius7_sta1 = 0.; rollRadius8_sta1 = 0.; 

 
  ndigi1_sta2 = 0; ndigi2_sta2 = 0; ndigi3_sta2 = 0; ndigi4_sta2 = 0; ndigi5_sta2 = 0; ndigi6_sta2 = 0; ndigi7_sta2 = 0; ndigi8_sta2 = 0; ndigi9_sta2 = 0; ndigi10_sta2 = 0; ndigi11_sta2 = 0; ndigi12_sta2 = 0;
  ndigiVsArea1_sta2 = 0.; ndigiVsArea2_sta2 = 0.; ndigiVsArea3_sta2 = 0.; ndigiVsArea4_sta2 = 0.; ndigiVsArea5_sta2 = 0.; ndigiVsArea6_sta2 = 0.; ndigiVsArea7_sta2 = 0.; ndigiVsArea8_sta2 = 0.; ndigiVsArea9_sta2 = 0.; ndigiVsArea10_sta2 = 0.; ndigiVsArea11_sta2 = 0.; ndigiVsArea12_sta2 = 0.;
  rollRadius1_sta2 = 0.; rollRadius2_sta2 = 0.; rollRadius3_sta2 = 0.; rollRadius4_sta2 = 0.; rollRadius5_sta2 = 0.; rollRadius6_sta2 = 0.; rollRadius7_sta2 = 0.; rollRadius8_sta2 = 0.; rollRadius9_sta2 = 0.; rollRadius10_sta2 = 0.; rollRadius11_sta2 = 0.; rollRadius12_sta2 = 0.;

 
}


MuonNewGEMDigis::~MuonNewGEMDigis()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

// ------------ method called once each job just before starting event loop  ------------
void
MuonNewGEMDigis::beginJob()
{
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuonNewGEMDigis::analyze(const edm::Event& e, const edm::EventSetup& iSetup)
{
    using namespace edm;
   

   // e.getByToken(simTrackInput_, sim_tracks);
    //e.getByToken(simTrackInput_, sim_vertices);
    //e.getByToken(gemSimHitInput_, GEMHits);
    //e.getByToken(gemDigiInput_, gem_digis);

    edm::ESHandle<GEMGeometry> pDD;
    iSetup.get<MuonGeometryRecord>().get( pDD ); 

    edm::Handle<GEMDigiCollection> gem_digis;
    e.getByToken(gemDigiInput_, gem_digis);

    edm::Handle<edm::PSimHitContainer> GEMHits;
    e.getByToken(gemSimHitInput_, GEMHits);

    edm::Handle<edm::SimTrackContainer> sim_tracks;
    e.getByToken(simTrackInput_, sim_tracks);


    std::vector<int> trackIds;
    std::vector<int> trackType;
    //const edm::SimTrackContainer & sim_trks = *sim_tracks.product();

    
   
/* 
   for (auto& t: sim_trks)
    {
        if (!isSimTrackGood(t)) continue;
        trackType.push_back(t.type());
        trackIds.push_back(t.trackId());
    }

*/

    int countRoll1 = 0;
    int countRoll2 = 0;
    int countRoll3 = 0;
    int countRoll4 = 0;
    int countRoll5 = 0;
    int countRoll6 = 0;
    int countRoll7 = 0;
    int countRoll8 = 0;
    
    for(GEMDigiCollection::DigiRangeIterator detUnitIt = gem_digis->begin(); detUnitIt != gem_digis->end(); ++detUnitIt){

        const GEMDetId& id = (*detUnitIt).first;

        const GeomDet* gdet = pDD->idToDet(id);
           if ( gdet == nullptr) { 
                std::cout<<"Getting DetId failed. Discard this gem strip hit.Maybe it comes from unmatched geometry."<<std::endl;
            continue; 
           }

        
        const BoundPlane & surface = gdet->surface();
        const GEMEtaPartition* roll = pDD->etaPartition(id);

        //const GeomDet* gdet = gem_geom->idToDet(id);
        //const BoundPlane & surface = gdet->surface();

        //int region = (Short_t) id.region();
        std::cout <<" region " << (Short_t) id.region() <<std::endl; 
        std::cout <<" station " << (Short_t) id.station() <<std::endl;  
        std::cout <<" Roll " << (Short_t) id.roll() <<std::endl;
 
	int ndigi = 0;
	double trArea(0.0);
	double trStripArea(0.0);
	Local3DPoint locMuonEntry(0., 0., 0.);
	GlobalPoint globMuonEntry(0., 0., 0.);
	//double simMuPhi = -99.;
	//double deltaPhi = -99.;
	Local3DPoint locDigi(0., 0., 0.);
	GlobalPoint pointDigiHit;


        const TrapezoidalStripTopology* top_(dynamic_cast<const TrapezoidalStripTopology*> (&(roll->topology())));

        const float rollRadius = top_->radius();
	const float striplength(top_->stripLength());
	const int nstrips = roll->nstrips();
        float rollPitch = roll->pitch();     
 
        std::cout <<" number of strips " << nstrips <<std::endl;
        std::cout <<" strip length " << top_->stripLength() <<std::endl;
        std::cout <<" pitch " << roll->pitch() <<std::endl;

	trStripArea = (roll->pitch()) * striplength;
	trArea = trStripArea * nstrips;
        std::cout <<" trArea " << trArea <<std::endl;
        std::cout <<" roll radius " << rollRadius <<std::endl;

	if(id.roll() == 1) { countRoll1++;}
	if(id.roll() == 2) { countRoll2++;}
	if(id.roll() == 3) { countRoll3++;}
	if(id.roll() == 4) { countRoll4++;}
	if(id.roll() == 5) { countRoll5++;}
	if(id.roll() == 6) { countRoll6++;}
	if(id.roll() == 7) { countRoll7++;}
	if(id.roll() == 8) { countRoll8++;}

      
ndigi1_sta1 = 0; ndigi2_sta1 = 0; ndigi3_sta1 = 0; ndigi4_sta1 = 0; ndigi5_sta1 = 0; ndigi6_sta1 = 0; ndigi7_sta1 = 0; ndigi8_sta1 = 0;
ndigi1_sta2 = 0; ndigi2_sta2 = 0; ndigi3_sta2 = 0; ndigi4_sta2 = 0; ndigi5_sta2 = 0; ndigi6_sta2 = 0; ndigi7_sta2 = 0; ndigi8_sta2 = 0; ndigi9_sta2 = 0; ndigi10_sta2 = 0; ndigi11_sta2 = 0; ndigi12_sta2 = 0;  
	const GEMDigiCollection::Range& range = (*detUnitIt).second;
	for (GEMDigiCollection::const_iterator digiIt = range.first; digiIt!=range.second; ++digiIt)
	  {


            LocalPoint lp = roll->centreOfStrip(digiIt->strip());
            GlobalPoint gp = surface.toGlobal(lp);
            Float_t g_r = (Float_t) gp.perp();
             	   
            
            if((Short_t) id.station() == 1){ 

            //if(roll->pitch() == 0.111622 || roll->pitch() == 0.102626 || roll->pitch() == 0.0943929 || roll->pitch() == 0.0869219 || roll->pitch() == 0.0800638 || roll->pitch() == 0.0738186 || roll->pitch() == 0.0680703 || roll->pitch() == 0.0628188){
            std::cout <<" roll pitch " << rollPitch <<std::endl;
            //float rollpitch1 = 0.111622; 


           if( fabs( rollPitch- 0.111622) < 0.00001 || fabs( rollPitch - 0.102626) < 0.00001 || fabs( rollPitch - 0.0943929) < 0.00001 || fabs( rollPitch - 0.0869219) < 0.00001 || fabs( rollPitch - 0.0800638) < 0.00001 || fabs( rollPitch - 0.0738186) < 0.00001 || fabs( rollPitch - 0.0680703) < 0.00001 || fabs( rollPitch - 0.0628188) < 0.00001){
             std::cout <<" global r "<< g_r << " station = " << id.station() << " pitch = " <<roll->pitch() << std::endl;
             hglobalR_Sta1_long->Fill(fabs(g_r));
            }
            //hglobalR_Sta1->Fill(fabs(g_r));
            //

            if( fabs( rollPitch- 0.105415) < 0.00001 || fabs( rollPitch - 0.0976903) < 0.00001 || fabs( rollPitch - 0.0905665) < 0.00001 || fabs( rollPitch - 0.0840435) < 0.00001 || fabs( rollPitch - 0.0780059) < 0.00001 || fabs( rollPitch - 0.0724538) < 0.00001 || fabs( rollPitch - 0.0673066) < 0.00001 || fabs( rollPitch - 0.0625643) < 0.00001){
             std::cout <<" global r "<< g_r << " station = " << id.station() << " pitch = " <<roll->pitch() << std::endl;
             hglobalR_Sta1_short->Fill(fabs(g_r));
            }



	    if(id.roll() == 1) {hNstripEtaPart1->Fill(digiIt->strip()); ndigi1_sta1++;}
	    if(id.roll() == 2) {hNstripEtaPart2->Fill(digiIt->strip()); ndigi2_sta1++;}
	    if(id.roll() == 3) {hNstripEtaPart3->Fill(digiIt->strip()); ndigi3_sta1++;}
	    if(id.roll() == 4) {hNstripEtaPart4->Fill(digiIt->strip()); ndigi4_sta1++;}
	    if(id.roll() == 5) {hNstripEtaPart5->Fill(digiIt->strip()); ndigi5_sta1++;}
	    if(id.roll() == 6) {hNstripEtaPart6->Fill(digiIt->strip()); ndigi6_sta1++;}
	    if(id.roll() == 7) {hNstripEtaPart7->Fill(digiIt->strip()); ndigi7_sta1++;}
	    if(id.roll() == 8) {hNstripEtaPart8->Fill(digiIt->strip()); ndigi8_sta1++;}




	    }

 
            if((Short_t) id.station() == 2){
 
     
            if((Short_t) id.layer() == 1){
            std::cout <<" global r "<< g_r << " station = " << id.station() << " layer = " << id.layer() << std::endl;
            hglobalR_Sta2_lay1->Fill(fabs(g_r));

            }

            if((Short_t) id.layer() == 2){
            std::cout <<" global r "<< g_r << " station = " << id.station() << " layer = " << id.layer() << std::endl;
            hglobalR_Sta2_lay2->Fill(fabs(g_r));

            }

            if(id.roll() == 1) {ndigi1_sta2++;}
            if(id.roll() == 2) {ndigi2_sta2++;}
            if(id.roll() == 3) {ndigi3_sta2++;}
            if(id.roll() == 4) {ndigi4_sta2++;}
            if(id.roll() == 5) {ndigi5_sta2++;}
            if(id.roll() == 6) {ndigi6_sta2++;}
            if(id.roll() == 7) {ndigi7_sta2++;}
            if(id.roll() == 8) {ndigi8_sta2++;}
            if(id.roll() == 9) {ndigi9_sta2++;}
            if(id.roll() == 10) {ndigi10_sta2++;}
            if(id.roll() == 11) {ndigi11_sta2++;}
            if(id.roll() == 12) {ndigi12_sta2++;}

            }



	    
	    hBx->Fill(digiIt->bx()); 
	    ndigi++;


	    
	    if (digiIt->strip() < 1 || digiIt->strip() > roll->nstrips() )
	     		 {
		std::cout <<" XXXXXXXXXXXXX Problemt with "<<id<<"  a digi has strip# = "<<digiIt->strip()<<std::endl;
	      }

	  }
	
	hRadiusEtaPartVsNdigi->Fill(rollRadius, ndigi);
	hRadiusEtaPartVsNdigiOvTrArea->Fill(rollRadius,ndigi/trArea);
	hRadiusEtaPart->Fill(rollRadius);

        std::cout <<" digi in roll 1 in station 1 "<< ndigi1_sta1 << " tr area " << trArea << std::endl;
        std::cout <<" digi in roll 2 in station 1 "<< ndigi2_sta1 << " tr area " << trArea << std::endl;
	std::cout <<" digi in roll 3 in station 1 "<< ndigi3_sta1 << " tr area " << trArea << std::endl;
        std::cout <<" digi in roll 4 in station 1 "<< ndigi4_sta1 << " tr area " << trArea << std::endl;
        std::cout <<" digi in roll 5 in station 1 "<< ndigi5_sta1 << " tr area " << trArea << std::endl;
        std::cout <<" digi in roll 6 in station 1 "<< ndigi6_sta1 << " tr area " << trArea << std::endl;
        std::cout <<" digi in roll 7 in station 1 "<< ndigi7_sta1 << " tr area " << trArea << std::endl;
        std::cout <<" digi in roll 8 in station 1 "<< ndigi8_sta1 << " tr area " << trArea << std::endl;

        if((Short_t) id.station() == 1){
        if(id.roll() == 1) {ndigiVsArea1_sta1 = ndigiVsArea1_sta1 + ndigi1_sta1*1./trArea; rollRadius1_sta1 = rollRadius;}
	if(id.roll() == 2) {ndigiVsArea2_sta1 = ndigiVsArea2_sta1 + ndigi2_sta1*1./trArea; rollRadius2_sta1 = rollRadius;}
	if(id.roll() == 3) {ndigiVsArea3_sta1 = ndigiVsArea3_sta1 + ndigi3_sta1*1./trArea; rollRadius3_sta1 = rollRadius;}
	if(id.roll() == 4) {ndigiVsArea4_sta1 = ndigiVsArea4_sta1 + ndigi4_sta1*1./trArea; rollRadius4_sta1 = rollRadius;}
	if(id.roll() == 5) {ndigiVsArea5_sta1 = ndigiVsArea5_sta1 + ndigi5_sta1*1./trArea; rollRadius5_sta1 = rollRadius;}
	if(id.roll() == 6) {ndigiVsArea6_sta1 = ndigiVsArea6_sta1 + ndigi6_sta1*1./trArea; rollRadius6_sta1 = rollRadius;}
	if(id.roll() == 7) {ndigiVsArea7_sta1 = ndigiVsArea7_sta1 + ndigi7_sta1*1./trArea; rollRadius7_sta1 = rollRadius;}
	if(id.roll() == 8) {ndigiVsArea8_sta1 = ndigiVsArea8_sta1 + ndigi8_sta1*1./trArea; rollRadius8_sta1 = rollRadius;}


        //std::cout <<" digi in roll 1 "<< ndigi1 << "in station " << (Short_t) id.station() << std::endl;
        //std::cout <<" digi in roll 2 "<< ndigi2 << "in station " << (Short_t) id.station() << std::endl;
        //std::cout <<" digi in roll 3 "<< ndigi3 << "in station " << (Short_t) id.station() << std::endl;
        //std::cout <<" digi in roll 4 "<< ndigi4 << "in station " << (Short_t) id.station() << std::endl;
        //std::cout <<" digi in roll 5 "<< ndigi5 << "in station " << (Short_t) id.station() << std::endl;
        //std::cout <<" digi in roll 6 "<< ndigi6 << "in station " << (Short_t) id.station() << std::endl;
        //std::cout <<" digi in roll 7 "<< ndigi7 << "in station " << (Short_t) id.station() << std::endl;
        //std::cout <<" digi in roll 8 "<< ndigi8 << "in station " << (Short_t) id.station() << std::endl;



        }


       if((Short_t) id.station() == 2){
        if(id.roll() == 1) {ndigiVsArea1_sta2 = ndigiVsArea1_sta2 + ndigi1_sta2*1./trArea; rollRadius1_sta2 = rollRadius;}
        if(id.roll() == 2) {ndigiVsArea2_sta2 = ndigiVsArea2_sta2 + ndigi2_sta2*1./trArea; rollRadius2_sta2 = rollRadius;}
        if(id.roll() == 3) {ndigiVsArea3_sta2 = ndigiVsArea3_sta2 + ndigi3_sta2*1./trArea; rollRadius3_sta2 = rollRadius;}
        if(id.roll() == 4) {ndigiVsArea4_sta2 = ndigiVsArea4_sta2 + ndigi4_sta2*1./trArea; rollRadius4_sta2 = rollRadius;}
        if(id.roll() == 5) {ndigiVsArea5_sta2 = ndigiVsArea5_sta2 + ndigi5_sta2*1./trArea; rollRadius5_sta2 = rollRadius;}
        if(id.roll() == 6) {ndigiVsArea6_sta2 = ndigiVsArea6_sta2 + ndigi6_sta2*1./trArea; rollRadius6_sta2 = rollRadius;}
        if(id.roll() == 7) {ndigiVsArea7_sta2 = ndigiVsArea7_sta2 + ndigi7_sta2*1./trArea; rollRadius7_sta2 = rollRadius;}
        if(id.roll() == 8) {ndigiVsArea8_sta2 = ndigiVsArea8_sta2 + ndigi8_sta2*1./trArea; rollRadius8_sta2 = rollRadius;}
        if(id.roll() == 9) {ndigiVsArea9_sta2 = ndigiVsArea9_sta2 + ndigi9_sta2*1./trArea; rollRadius9_sta2 = rollRadius;}
        if(id.roll() == 10) {ndigiVsArea10_sta2 = ndigiVsArea10_sta2 + ndigi10_sta2*1./trArea; rollRadius10_sta2 = rollRadius;}
        if(id.roll() == 11) {ndigiVsArea11_sta2 = ndigiVsArea11_sta2 + ndigi11_sta2*1./trArea; rollRadius11_sta2 = rollRadius;}
        if(id.roll() == 12) {ndigiVsArea12_sta2 = ndigiVsArea12_sta2 + ndigi12_sta2*1./trArea; rollRadius12_sta2 = rollRadius;}

        }  
  
	
        }// for eta partitions (rolls)
    

    //std::cout << "roll 1 numbers = " <<  countRoll1 << "\tndigi = " << ndigiVsArea1 << std::endl;
    //std::cout << "roll 2 numbers = " <<  countRoll2 << "\tndigi = " << ndigiVsArea2 << std::endl;
    //std::cout << "roll 3 numbers = " <<  countRoll3 << "\tndigi = " << ndigiVsArea3 << std::endl;
    //std::cout << "roll 4 numbers = " <<  countRoll4 << "\tndigi = " << ndigiVsArea4 << std::endl;
    //std::cout << "roll 5 numbers = " <<  countRoll5 << "\tndigi = " << ndigiVsArea5 << std::endl;
    //std::cout << "roll 6 numbers = " <<  countRoll6 << "\tndigi = " << ndigiVsArea6 << std::endl;
    //std::cout << "roll 7 numbers = " <<  countRoll7 << "\tndigi = " << ndigiVsArea7 << std::endl;
    //std::cout << "roll 8 numbers = " <<  countRoll8 << "\tndigi = " << ndigiVsArea8 << std::endl;
    
    numbEvents++;
}



// ------------ method called once each job just after ending the event loop  ------------
void
MuonNewGEMDigis::endJob()
{
  
  std::cout << "number of events = " << numbEvents << std::endl;
  std::cout << "--------------" << std::endl;
  

  std::vector<double> myRadii_sta1, myRates_sta1;
 
  
  std::cout << "ndigiVsArea1_sta1 = " << ndigiVsArea1_sta1;
  ndigiVsArea1_sta1 = ndigiVsArea1_sta1/(numbEvents * 3 *25 * 1.0e-9 * 2. * 72. );
  std::cout << "\tRate [Hz/cm2] = " << ndigiVsArea1_sta1 << std::endl;


  myRadii_sta1.push_back(rollRadius1_sta1); myRates_sta1.push_back(ndigiVsArea1_sta1);
  
  std::cout << "ndigiVsArea2_sta1 = " << ndigiVsArea2_sta1;
  ndigiVsArea2_sta1 = ndigiVsArea2_sta1/(numbEvents * 3 *25 * 1.0e-9 * 2. * 72. );
  std::cout << "\tRate [Hz/cm2] = " << ndigiVsArea2_sta1 << std::endl;

  myRadii_sta1.push_back(rollRadius2_sta1); myRates_sta1.push_back(ndigiVsArea2_sta1);

  std::cout << "ndigiVsArea3_sta1 = " << ndigiVsArea3_sta1;
  ndigiVsArea3_sta1 = ndigiVsArea3_sta1/(numbEvents * 3 *25 * 1.0e-9 * 2. * 72. );
  std::cout << "\tRate [Hz/cm2] = " << ndigiVsArea3_sta1 << std::endl;

  myRadii_sta1.push_back(rollRadius3_sta1); myRates_sta1.push_back(ndigiVsArea3_sta1);

  std::cout << "ndigiVsArea4_sta1 = " << ndigiVsArea4_sta1;
  ndigiVsArea4_sta1 = ndigiVsArea4_sta1/(numbEvents * 3 *25 * 1.0e-9 * 2. * 72. );
  std::cout << "\tRate [Hz/cm2] = " << ndigiVsArea4_sta1 << std::endl;
  
  myRadii_sta1.push_back(rollRadius4_sta1); myRates_sta1.push_back(ndigiVsArea4_sta1);

  std::cout << "ndigiVsArea5_sta1 = " << ndigiVsArea5_sta1;
  ndigiVsArea5_sta1 = ndigiVsArea5_sta1/(numbEvents * 3 *25 * 1.0e-9 * 2. * 72. );
  std::cout << "\tRate [Hz/cm2] = " << ndigiVsArea5_sta1 << std::endl;

  myRadii_sta1.push_back(rollRadius5_sta1); myRates_sta1.push_back(ndigiVsArea5_sta1);

  std::cout << "ndigiVsArea6_sta1 = " << ndigiVsArea6_sta1;
  ndigiVsArea6_sta1 = ndigiVsArea6_sta1/(numbEvents * 3 *25 * 1.0e-9 * 2. * 72. );
  std::cout << "\tRate [Hz/cm2] = " << ndigiVsArea6_sta1 << std::endl;

  myRadii_sta1.push_back(rollRadius6_sta1); myRates_sta1.push_back(ndigiVsArea6_sta1);

  std::cout << "ndigiVsArea7_sta1 = " << ndigiVsArea7_sta1;
  ndigiVsArea7_sta1 = ndigiVsArea7_sta1/(numbEvents * 3 *25 * 1.0e-9 * 2. * 72. );
  std::cout << "\tRate [Hz/cm2] = " << ndigiVsArea7_sta1 << std::endl;

 
  myRadii_sta1.push_back(rollRadius7_sta1); myRates_sta1.push_back(ndigiVsArea7_sta1);
  
  std::cout << "ndigiVsArea8_sta1 = " << ndigiVsArea8_sta1;
  ndigiVsArea8_sta1 = ndigiVsArea8_sta1/(numbEvents * 3 *25 * 1.0e-9 * 2. * 72. );
  std::cout << "\tRate [Hz/cm2] = " << ndigiVsArea8_sta1 << std::endl;

  myRadii_sta1.push_back(rollRadius8_sta1); myRates_sta1.push_back(ndigiVsArea8_sta1);

  std::cout << "rollRadius_sta1[cm]\tRate_sta1[Hz/cm2]" << std::endl;
  std::cout << rollRadius1_sta1 << "\t" << ndigiVsArea1_sta1 << std::endl;
  std::cout << rollRadius2_sta1 << "\t" << ndigiVsArea2_sta1 << std::endl;
  std::cout << rollRadius3_sta1 << "\t" << ndigiVsArea3_sta1 << std::endl;
  std::cout << rollRadius4_sta1 << "\t" << ndigiVsArea4_sta1 << std::endl;
  std::cout << rollRadius5_sta1 << "\t" << ndigiVsArea5_sta1 << std::endl;
  std::cout << rollRadius6_sta1 << "\t" << ndigiVsArea6_sta1 << std::endl;
  std::cout << rollRadius7_sta1 << "\t" << ndigiVsArea7_sta1 << std::endl;
  std::cout << rollRadius8_sta1 << "\t" << ndigiVsArea8_sta1 << std::endl;
 
  for (unsigned int i = 0; i < myRadii_sta1.size(); i++)
  {
    std::cout << "radius_sta1 = " << myRadii_sta1[i] << "\tRate_sta1 = " << myRates_sta1[i] << std::endl;
    grRatePerRoll_sta1->SetPoint(i, myRadii_sta1[i], myRates_sta1[i]);

  }


  //++++++++++++++++++++++++++++ for Station 2 +++++++++++++++++++++++++++++++++++++++++++++++++//


   std::vector<double> myRadii_sta2, myRates_sta2;

   std::cout << "ndigiVsArea1_sta2 = " << ndigiVsArea1_sta2;
   ndigiVsArea1_sta2 = ndigiVsArea1_sta2/(numbEvents * 3 *25 * 1.0e-9 * 2. * 36. );
   myRadii_sta2.push_back(rollRadius1_sta2); myRates_sta2.push_back(ndigiVsArea1_sta2);
   std::cout << "\tRate [Hz/cm2] = " << ndigiVsArea1_sta2 << std::endl;

   std::cout << "ndigiVsArea2_sta2 = " << ndigiVsArea2_sta2;
   ndigiVsArea2_sta2 = ndigiVsArea2_sta2/(numbEvents * 9 *25 * 1.0e-9 * 2. * 36. );
   myRadii_sta2.push_back(rollRadius2_sta2); myRates_sta2.push_back(ndigiVsArea2_sta2);
   std::cout << "\tRate [Hz/cm2] = " << ndigiVsArea2_sta2 << std::endl;   


   std::cout << "ndigiVsArea3_sta2 = " << ndigiVsArea3_sta2;
   ndigiVsArea3_sta2 = ndigiVsArea3_sta2/(numbEvents * 9 *25 * 1.0e-9 * 2. * 36. );
   myRadii_sta2.push_back(rollRadius3_sta2); myRates_sta2.push_back(ndigiVsArea3_sta2);
   std::cout << "\tRate [Hz/cm2] = " << ndigiVsArea3_sta2 << std::endl;


   std::cout << "ndigiVsArea4_sta2 = " << ndigiVsArea4_sta2;
   ndigiVsArea4_sta2 = ndigiVsArea4_sta2/(numbEvents * 9 *25 * 1.0e-9 * 2. * 36. );
   myRadii_sta2.push_back(rollRadius4_sta2); myRates_sta2.push_back(ndigiVsArea4_sta2);
   std::cout << "\tRate [Hz/cm2] = " << ndigiVsArea4_sta2 << std::endl;

   std::cout << "ndigiVsArea5_sta2 = " << ndigiVsArea5_sta2;
   ndigiVsArea5_sta2 = ndigiVsArea5_sta2/(numbEvents * 9 *25 * 1.0e-9 * 2. * 36. );
   myRadii_sta2.push_back(rollRadius5_sta2); myRates_sta2.push_back(ndigiVsArea5_sta2);
   std::cout << "\tRate [Hz/cm2] = " << ndigiVsArea5_sta2 << std::endl;

   std::cout << "ndigiVsArea6_sta2 = " << ndigiVsArea6_sta2;
   ndigiVsArea6_sta2 = ndigiVsArea6_sta2/(numbEvents * 9 *25 * 1.0e-9 * 2. * 36. );
   myRadii_sta2.push_back(rollRadius6_sta2); myRates_sta2.push_back(ndigiVsArea6_sta2);
   std::cout << "\tRate [Hz/cm2] = " << ndigiVsArea6_sta2 << std::endl;

   std::cout << "ndigiVsArea7_sta2 = " << ndigiVsArea7_sta2;
   ndigiVsArea7_sta2 = ndigiVsArea7_sta2/(numbEvents * 9 *25 * 1.0e-9 * 2. * 36. );
   myRadii_sta2.push_back(rollRadius7_sta2); myRates_sta2.push_back(ndigiVsArea7_sta2);
   std::cout << "\tRate [Hz/cm2] = " << ndigiVsArea7_sta2 << std::endl;


   std::cout << "ndigiVsArea8_sta2 = " << ndigiVsArea8_sta2;
   ndigiVsArea8_sta2 = ndigiVsArea8_sta2/(numbEvents * 9 *25 * 1.0e-9 * 2. * 36. );
   myRadii_sta2.push_back(rollRadius8_sta2); myRates_sta2.push_back(ndigiVsArea8_sta2);
   std::cout << "\tRate [Hz/cm2] = " << ndigiVsArea8_sta2 << std::endl;

   std::cout << "ndigiVsArea9_sta2 = " << ndigiVsArea9_sta2;
   ndigiVsArea9_sta2 = ndigiVsArea9_sta2/(numbEvents * 9 *25 * 1.0e-9 * 2. * 36. );
   myRadii_sta2.push_back(rollRadius9_sta2); myRates_sta2.push_back(ndigiVsArea9_sta2);
   std::cout << "\tRate [Hz/cm2] = " << ndigiVsArea9_sta2 << std::endl;

   std::cout << "ndigiVsArea10_sta2 = " << ndigiVsArea10_sta2;
   ndigiVsArea10_sta2 = ndigiVsArea10_sta2/(numbEvents * 9 *25 * 1.0e-9 * 2. * 36. );
   myRadii_sta2.push_back(rollRadius10_sta2); myRates_sta2.push_back(ndigiVsArea10_sta2);
   std::cout << "\tRate [Hz/cm2] = " << ndigiVsArea10_sta2 << std::endl;

   std::cout << "ndigiVsArea11_sta2 = " << ndigiVsArea11_sta2;
   ndigiVsArea11_sta2 = ndigiVsArea11_sta2/(numbEvents * 9 *25 * 1.0e-9 * 2. * 36. );
   myRadii_sta2.push_back(rollRadius11_sta2); myRates_sta2.push_back(ndigiVsArea11_sta2);
   std::cout << "\tRate [Hz/cm2] = " << ndigiVsArea11_sta2 << std::endl;

   std::cout << "ndigiVsArea12_sta2 = " << ndigiVsArea12_sta2;
   ndigiVsArea12_sta2 = ndigiVsArea12_sta2/(numbEvents * 9 *25 * 1.0e-9 * 2. * 36. );
   myRadii_sta2.push_back(rollRadius12_sta2); myRates_sta2.push_back(ndigiVsArea12_sta2);
   std::cout << "\tRate [Hz/cm2] = " << ndigiVsArea12_sta2 << std::endl;

  std::cout << "rollRadius_sta2[cm]\tRate_sta2[Hz/cm2]" << std::endl;
  std::cout << rollRadius1_sta2 << "\t" << ndigiVsArea1_sta2 << std::endl;
  std::cout << rollRadius2_sta2 << "\t" << ndigiVsArea2_sta2 << std::endl;
  std::cout << rollRadius3_sta2 << "\t" << ndigiVsArea3_sta2 << std::endl;
  std::cout << rollRadius4_sta2 << "\t" << ndigiVsArea4_sta2 << std::endl;
  std::cout << rollRadius5_sta2 << "\t" << ndigiVsArea5_sta2 << std::endl;
  std::cout << rollRadius6_sta2 << "\t" << ndigiVsArea6_sta2 << std::endl;
  std::cout << rollRadius7_sta2 << "\t" << ndigiVsArea7_sta2 << std::endl;
  std::cout << rollRadius8_sta2 << "\t" << ndigiVsArea8_sta2 << std::endl;
  std::cout << rollRadius9_sta2 << "\t" << ndigiVsArea9_sta2 << std::endl;
  std::cout << rollRadius10_sta2 << "\t" << ndigiVsArea10_sta2 << std::endl;
  std::cout << rollRadius11_sta2 << "\t" << ndigiVsArea11_sta2 << std::endl;
  std::cout << rollRadius12_sta2 << "\t" << ndigiVsArea12_sta2 << std::endl;



 
  for (unsigned int i = 0; i < myRadii_sta2.size(); i++)
  {
    std::cout << "radius_sta2 = " << myRadii_sta2[i] << "\tRate_sta2 = " << myRates_sta2[i] << std::endl;
    grRatePerRoll_sta2->SetPoint(i, myRadii_sta2[i], myRates_sta2[i]);

  }


}



// ------------ method called when starting to processes a luminosity block  ------------
/*
void
MuonNewGEMDigis::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
MuonNewGEMDigis::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
 MuonNewGEMDigis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
     edm::ParameterSetDescription desc;
     desc.setUnknown();
     descriptions.addDefault(desc);
     }
           

//define this as a plug-in
DEFINE_FWK_MODULE(MuonNewGEMDigis);
