// -*- C++ -*-
//
// Package:    ME0TimingAnalysis
// Class:      ME0TimingAnalysis
// 
/**\class ME0TimingAnalysis ME0TimingAnalysis.cc Analyzer/ME0TimingAnalysis/plugins/ME0TimingAnalysis.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author: Archie Sharma 
//         Created:  Tue, 22 Sep 2015 14:16:56 GMT
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>
#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "RecoMuon/MuonIdentification/interface/ME0MuonSelector.h"
#include <DataFormats/MuonReco/interface/ME0Muon.h>
#include <DataFormats/MuonReco/interface/ME0MuonCollection.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartition.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "DataFormats/MuonDetId/interface/ME0DetId.h"
#include "RecoMuon/MuonIdentification/plugins/ME0MuonSelector.cc"
//
// class declaration
//
using namespace reco;
using namespace std;
class ME0TimingAnalysis : public edm::EDAnalyzer {
	public:
		explicit ME0TimingAnalysis(const edm::ParameterSet&);
		~ME0TimingAnalysis();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;
		// ----------member data ---------------------------
		int tmpindex;
                int matchindex; 
		bool MatchedMuon(vector<int> me0muons, int recomuon) ;
		edm::Service<TFileService> fs;
    
        std::string wp_;
        std::string fD_;
        std::string TD_;
        double etaMin_;
        double etaMax_;
        double dr_;
        double dr1_;
        double ptMin_;
        double ptMinGen_;
        double timeMin_;
        double timeMax_;

        int evts = 0;
        TH1F *hFillEvents;
		TH1F *hFillSignalMuontime,*hFillPUMuontime;
		TH1F *hFillSignalMuontimeErr,*hFillPUMuontimeErr;
		TH1F *hFillZMass , *hFillZGenMass, *hFillRecoEta, *hFillRecoPt, *hFillAllRecoEta, *hFillAllRecoPt;
		TH1F *SignalMuonTime, *BGMuonTime;
		TH1F *hFillZMassInWindow,*hFillZMassOutWindow;
		TH1F *hFillPUMuonPt, *hFillGenMuEta, *hFillGenMuPt, *hFillPUMuonPtInWindow, *hFillPUMuonPtOutWindow;
        TH1F *hFillPUMuonEta, *hFillGenMuinME0Eta, *hFillGenMuinME0Pt, *hFillPUMuonEtaInWindow, *hFillPUMuonEtaOutWindow;
        TH1F *hFillPUMuonPt_rebin, *hFillPUMuonPtInWindow_rebin, *hFillPUMuonPtOutWindow_rebin;
        TH1F *hFillGenMuEta_allEta, *hFillGenMuPt_allEta, *hFillGenMuEta_me0Eta, *hFillGenMuPt_me0Eta;
        TH1F *hFillRecoMuonTimeErr, *hFillRecoMuonTime;
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
ME0TimingAnalysis::ME0TimingAnalysis(const edm::ParameterSet& iConfig):
  wp_(iConfig.getParameter<std::string>("wp")),
  fD_(iConfig.getParameter<std::string>("fD")), 
  TD_(iConfig.getParameter<std::string>("TD")),
  etaMin_(iConfig.getParameter<double>("etaMin")),
  etaMax_(iConfig.getParameter<double>("etaMax")),
  dr_(iConfig.getParameter<double>("dr")),
  dr1_(iConfig.getParameter<double>("dr1")),
  ptMin_(iConfig.getParameter<double>("ptMin")),
  ptMinGen_(iConfig.getParameter<double>("ptMinGen")),
  timeMin_(iConfig.getParameter<double>("timeMin")),
  timeMax_(iConfig.getParameter<double>("timeMax"))
{
	//now do what ever initialization is needed

}


ME0TimingAnalysis::~ME0TimingAnalysis()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
	void
ME0TimingAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   std::cout << "<===== New event started =======>" << std::endl;
    evts++;
    hFillEvents->Fill(evts);
	using namespace edm;
    
	/*
	   vector<reco::ME0Muon>                 "me0SegmentMatching"        ""                "RECO"       
    */
    
	edm::Handle <std::vector<ME0Muon> > OurMuons;
	iEvent.getByLabel <std::vector<ME0Muon> > ("me0SegmentMatching", OurMuons);

	edm::Handle <reco::GenParticleCollection> genparticles;
	iEvent.getByLabel("genParticles",genparticles);

    edm::ESHandle<ME0Geometry> me0geom;
    iSetup.get<MuonGeometryRecord>().get(me0geom);
    
    edm::Handle <std::vector<reco::Muon> > muons;
    iEvent.getByLabel("muons", muons );
 
 	TLorentzVector genmuon1, genmuon2;                                                             

	vector<unsigned int> indexmu;
    indexmu.clear();

	vector<unsigned int> indexme0mu;
	indexme0mu.clear();
    
    vector<unsigned int> indexme0muOut;
    indexme0muOut.clear();
    
	//============ gen particle collection  =============//

	for(unsigned int i = 0; i < genparticles->size();i++) {
        
        //std::cout<<"pdgID: "<<genparticles->at(i).pdgId()<<" Status: "<<genparticles->at(i).status()<<" Eta: "<<genparticles->at(i).eta()<<" Phi: "<<genparticles->at(i).phi()<<std::endl;
        
        //if(!(genparticles->at(i).pt() > ptMin_)) continue;
        
        if(!(abs(genparticles->at(i).pdgId()) == 13 && genparticles->at(i).status() == 1 && genparticles->at(i).numberOfMothers() > 0)) continue;
        
        if(fabs(genparticles->at(i).mother()->pdgId()) == 23 || fabs(genparticles->at(i).mother()->pdgId()) == 22) indexmu.push_back(i);
        else if(abs(genparticles->at(i).pdgId()) == abs(genparticles->at(i).mother()->pdgId())) {
            
            if(genparticles->at(i).mother()->numberOfMothers() > 0) {

                if(abs(genparticles->at(i).mother()->mother()->pdgId()) == 23 || fabs(genparticles->at(i).mother()->pdgId()) == 22) indexmu.push_back(i);

            }
            
            if(genparticles->at(i).mother()->mother()->numberOfMothers() > 0) {
		  
                if(abs(genparticles->at(i).mother()->mother()->mother()->pdgId()) == 23 || fabs(genparticles->at(i).mother()->pdgId()) == 22) indexmu.push_back(i);
                
            }
            
        }
		
	}//close for gen particle loop

    std::cout << "size of all muon vector " << indexmu.size() << std::endl;
        
	for(unsigned int j = 0 ; j < indexmu.size(); j++){
        
        hFillGenMuEta_allEta->Fill(fabs(genparticles->at(indexmu.at(j)).eta()));
        hFillGenMuPt_allEta->Fill(genparticles->at(indexmu.at(j)).pt());
        std::cout << "eta of selected  gen muon " << genparticles->at(indexmu.at(j)).eta() << std::endl;
        std::cout << "pt of selected  gen muon " << genparticles->at(indexmu.at(j)).pt() << std::endl;       
     //   std::cout << "status of selected  gen muon " << genparticles->at(indexmu.at(j)).status() << " decay in to particle with pgdId " << genparticles->at(indexmu.at(j)).daughter()->pdgId() << " with pt = "<< genparticles->at(indexmu.at(j)).daughter()->pt() << std::endl;
 
        if(fabs(genparticles->at(indexmu.at(j)).eta()) > etaMin_ && fabs(genparticles->at(indexmu.at(j)).eta()) < etaMax_ && genparticles->at(indexmu.at(j)).pt()> ptMin_) {

            std::cout << "me0 muon eta " << genparticles->at(indexmu.at(j)).eta() << "at j " << j << std::endl;
            indexme0mu.push_back(indexmu.at(j));
            
        }
        else indexme0muOut.push_back(indexmu.at(j));

    }
    std::cout << "size of ME0 muon vector " << indexme0mu.size() << std::endl;
    
    
    //check if there is at least one mu in me0 and the other is elsewhere
    bool isOneMuinME0=false;
    bool isTwoMuinME0=false;
    if(indexme0mu.size() == 1 && indexme0muOut.size() == 1) isOneMuinME0 = true;
    if(indexme0mu.size() == 2 && indexme0muOut.size() == 0) isTwoMuinME0 = true;
    
    //match the elsewhere muon with the reco::muon
    uint kk = 0;
    std::vector<int> recomuout;
    
    if(isOneMuinME0 == true){
        
        std::cout<<" <---------------------- case one mu in me0, one mu out ----------------------> "<<std::endl;
	    double DRmintmp = dr1_;
        
	    for ( std::vector<reco::Muon>::const_iterator mu = muons->begin(); mu != muons->end(); ++mu, ++kk){
            
            int idxprovareco = -99;
            double drReco =  reco::deltaR(genparticles->at(indexme0muOut[0]), *mu);
            if(drReco < DRmintmp){
                
                DRmintmp = drReco;
                idxprovareco = kk;
                //std::cout<<" drreco "<<DRmintmp<<" idx reco "<<idxprovareco<<"  idxgen "<<genmuout[0]<<std::endl;
                recomuout.push_back(idxprovareco);
                
            }
	    
	    //std::cout<<" DR MIN FUORI ME0 cov. "<<DRmintmp<<" idx "<<idxprovareco<<std::endl;	   
        }
	    if(recomuout.size()) std::cout<<"OUT::recomuidx="<<recomuout[0]<<" genmu="<<indexme0muOut[0]<<" dr="<< reco::deltaR(genparticles->at(indexme0muOut[0]),muons->at(recomuout[0]))<<" eta:"<<muons->at(recomuout[0]).eta()<<std::endl;

    }//elsewhere muon event
    

    //if(indexme0mu.size() == 0) return;
    if(indexmu.size() != 2) return;
    
    for(unsigned int n = 0 ; n < indexme0mu.size(); n++){
        
        if((genparticles->at(indexme0mu.at(n)).pt() > ptMin_)) hFillGenMuEta_me0Eta->Fill(fabs(genparticles->at(indexme0mu.at(n)).eta()));
        hFillGenMuPt_me0Eta->Fill(genparticles->at(indexme0mu.at(n)).pt());
        std::cout << "Signal muon eta ==> " << genparticles->at(indexme0mu.at(n)).eta() << "  Signal muon pt ==> " << genparticles->at(indexme0mu.at(n)).pt() << "  Signal muon phi ==> " << genparticles->at(indexme0mu.at(n)).phi() << std::endl;

    }

	//double DRtmp = dr_;
	
	std::vector<bool> IsMatched;
	std::vector<int> me0muons;
        std::vector<int> me0muonsall;
        std::vector<int> assGenMuons;
        me0muonsall.clear();
	me0muons.clear(); 
        assGenMuons.clear();
       
    
	for( unsigned int j = 0;  j < indexme0mu.size(); j++) {

         std::cout << "inside gen particle loop " << std::endl;

        std::cout << "gen ME0 muon eta " << genparticles->at(indexme0mu.at(j)).eta() << "  gen ME0 muon pt " << genparticles->at(indexme0mu.at(j)).pt() << "  gen ME0 muon phi " << genparticles->at(indexme0mu.at(j)).phi() <<  std::endl;

           tmpindex = -1;
           matchindex = -1;
         bool allme0mu = false;
         double dr = 0.0;
         double pTRes = 0.0;
         double DRtmp = dr_;
         double tmpPTRes = 10;
 	for(unsigned int t = 0; t < OurMuons->size(); t++) {
     
               
           allme0mu = true;
           bool isintime = false;
          
           if(OurMuons->at(t).me0segment().time() >= timeMin_  && OurMuons->at(t).me0segment().time() <= timeMax_) isintime = true;              
           bool decisionTime = (TD_ == "true") ? isintime :  allme0mu;
                    
           if(!decisionTime)continue;
           if(!(OurMuons->at(t).pt() > ptMin_))continue;
            
           
            bool isLoose = muon::isGoodMuon(me0geom, OurMuons->at(t), muon::Loose);
            bool isTight = muon::isGoodMuon(me0geom, OurMuons->at(t), muon::Tight);
            bool decision = (wp_ == "loose") ? isLoose : isTight;
            bool decisionN = (!decision);
            bool decisionF = (fD_ == "ID") ? decision : decisionN;        
 
            if (!decisionF) continue;
            //if(OurMuons->at(t).pt() < ptMin_) continue;
        
            if(int(t) == tmpindex) continue;
                      
          //  std::cout << "reco ME0 muon eta " << OurMuons->at(t).eta() << "  reco ME0 muon pt " << OurMuons->at(t).pt() << "  reco ME0 muon phi " << OurMuons->at(t).phi() <<  std::endl;
            
		 dr = reco::deltaR(genparticles->at(indexme0mu.at(j)).eta(), genparticles->at(indexme0mu.at(j)).phi(),OurMuons->at(t).eta(),OurMuons->at(t).phi());
           
                      std::cout << " dR for all gen reco pairs" << dr << std::endl;           
 
			if(dr < DRtmp) {
                
                            // std::cout << " dr is less then DRtemp" << dr << std::endl;
                            // DRtmp = dr;
		             tmpindex = t;
                             if(tmpindex != -1) me0muonsall.push_back(tmpindex);                                                       
                
			}

                 
             //std::cout << "minimum dR " << DRtmp << std::endl;
            //std::cout << "eta of matched me0 muon inside dr loop " << OurMuons->at(t).eta() << std::endl;
            
		}   /////////////////////////////////////////
           
           std::cout << "size of matched me0 muon vector " << me0muonsall.size() << std::endl;
           for( unsigned int d = 0;  d < me0muonsall.size(); d++) {

              if(int(d) == matchindex) continue;
              double ptrec_meoMuons = OurMuons->at(me0muonsall.at(d)).pt();               
              double ptsim_me0Muons = genparticles->at(indexme0mu.at(j)).pt();   
              pTRes = (ptsim_me0Muons - ptrec_meoMuons)/ptsim_me0Muons ;
             // std::cout << "pt Resolution " << pTRes << " sim pt for pt Resolution " << ptsim_me0Muons << " rec pt for pt Resolution " << ptrec_meoMuons << std::endl;
              if (pTRes < tmpPTRes){
                 
                  tmpPTRes = pTRes;
                  matchindex = d; 
                  }

              std::cout << "Min. pt Resolution " << pTRes << "  for me0 index " << matchindex <<std::endl;
             }
          
           //std::cout << "final pt Resolution " << pTRes << "  for me0 index " << matchindex << std::endl;
           if(matchindex != -1) {
            
            me0muons.push_back(me0muonsall.at(matchindex));
            assGenMuons.push_back(indexme0mu.at(j));
         //   std::cout << "eta of me0muons vector " << OurMuons->at(me0muons.at(matchindex)).eta() << " pt of me0muons vector " << OurMuons->at(me0muons.at(matchindex)).pt() << " phi of me0muons vector " << OurMuons->at(me0muons.at(matchindex)).phi() << " index of  me0muons vector " << matchindex <<  std::endl;
        
        }
        
	}
    
    for( unsigned int k = 0;  k < assGenMuons.size(); k++) {
 
       std::cout << "Matched Signal muon eta ==> " << genparticles->at(assGenMuons.at(k)).eta() << "  Matched Signal muon pt ==> " << genparticles->at(assGenMuons.at(k)).pt() << "  Matched Signal muon phi ==> " << genparticles->at(assGenMuons.at(k)).phi() << std::endl;

        if((genparticles->at(assGenMuons.at(k)).pt() > ptMin_)) hFillGenMuEta->Fill(fabs(genparticles->at(assGenMuons.at(k)).eta()));
        hFillGenMuPt->Fill(genparticles->at(assGenMuons.at(k)).pt());

    }

    //std::cout << "Size of reco ME0 muon vector " << me0muons.size() << std::endl;
    //std::cout << "Size of ass gen muon vector " << assGenMuons.size() << std::endl;
    
	/// again filling

    if(me0muons.size() > 0 && me0muons.at(0) != -1){
      
        std::cout << "eta of first matched me0 reco muon " << OurMuons->at(me0muons.at(0)).eta() << " pt of first matched me0 reco muon " << OurMuons->at(me0muons.at(0)).pt() <<" phi of first matched me0 reco muon " << OurMuons->at(me0muons.at(0)).phi() <<std::endl;
        hFillSignalMuontime->Fill(OurMuons->at(me0muons.at(0)).me0segment().time());
        hFillSignalMuontimeErr->Fill(OurMuons->at(me0muons.at(0)).me0segment().timeErr());
        SignalMuonTime->Fill(OurMuons->at(me0muons.at(0)).me0segment().time(), OurMuons->at(me0muons.at(0)).me0segment().timeErr());
  
    }
    if(me0muons.size() > 1 && me0muons.at(1) != -1){
      
        std::cout << "eta of second matched me0 reco muon " << OurMuons->at(me0muons.at(1)).eta() << std::endl;
        hFillSignalMuontime->Fill(OurMuons->at(me0muons.at(1)).me0segment().time());
        hFillSignalMuontimeErr->Fill(OurMuons->at(me0muons.at(1)).me0segment().timeErr());
    
    }
    
// SignalMuonTime->Fill(OurMuons->at(me0muons.at(0)).me0segment().time(), OurMuons->at(me0muons.at(0)).me0segment().timeErr()); }
//	if(me0muons.size() > 2) {hFillSignalMuontime->Fill(OurMuons->at(me0muons.at(2)).me0segment().time()); hFillSignalMuontimeErr->Fill(OurMuons->at(me0muons.at(2)).me0segment().timeErr());SignalMuonTime->Fill(OurMuons->at(me0muons.at(0)).me0segment().time(), OurMuons->at(me0muons.at(0)).me0segment().timeErr());}

	//now fill the Z mass plots
	TLorentzVector recomu1, recomu2, ZCandReco, genmu1, genmu2, ZCandGen;
	if(isOneMuinME0 && recomuout.size() > 0 && me0muons.size() > 0 ){
        
       // std::cout<<"genmu out eta="<<genparticles->at(indexme0muOut[0]).eta()<<" recomu out eta="<< muons->at(recomuout[0]).eta()<<std::endl;
       // std::cout<<"genmu me0 eta="<<genparticles->at(indexme0mu[0]).eta()<<" recomu me0 eta="<<OurMuons->at(me0muons[0]).eta()<<std::endl;

        recomu1.SetPtEtaPhiM(muons->at(recomuout[0]).pt(), muons->at(recomuout[0]).eta(),muons->at(recomuout[0]).phi(),0.105);
        recomu2.SetPtEtaPhiM(OurMuons->at(me0muons[0]).pt(), OurMuons->at(me0muons[0]).eta(),OurMuons->at(me0muons[0]).phi(), 0.105);
        genmu1.SetPtEtaPhiM(genparticles->at(indexme0muOut[0]).pt(),genparticles->at(indexme0muOut[0]).eta(), genparticles->at(indexme0muOut[0]).phi(), 0.105);
        genmu2.SetPtEtaPhiM(genparticles->at(indexme0mu[0]).pt(),genparticles->at(indexme0mu[0]).eta(), genparticles->at(indexme0mu[0]).phi(), 0.105);
        ZCandReco = recomu1+recomu2;
        ZCandGen = genmu1+genmu2;
        std::cout<<"TipoA: zcand reco  "<< ZCandReco.M()<<" ZCand gen"<<ZCandGen.M()<<endl;

        hFillZGenMass->Fill((ZCandGen).M());
        hFillZMass->Fill((ZCandReco).M());

        if((OurMuons->at(me0muons.at(0)).me0segment().time() >= timeMin_  && OurMuons->at(me0muons.at(0)).me0segment().time() <= timeMax_)){
          
            hFillZMassInWindow->Fill((ZCandReco).M());
          
        } else {
          
            hFillZMassOutWindow->Fill((ZCandReco).M());
          
        }
	}// EVENTO A
    
    if (isTwoMuinME0 && me0muons.size() > 1){
        
        recomu1.SetPtEtaPhiM(OurMuons->at(me0muons[0]).pt(),OurMuons->at(me0muons[0]).eta(),OurMuons->at(me0muons[0]).phi(), 0.105);
        recomu2.SetPtEtaPhiM(OurMuons->at(me0muons[1]).pt(),OurMuons->at(me0muons[1]).eta(),OurMuons->at(me0muons[1]).phi(), 0.105);
        genmu1.SetPtEtaPhiM(genparticles->at(indexme0mu[0]).pt(),genparticles->at(indexme0mu[0]).eta(), genparticles->at(indexme0mu[0]).phi(), 0.105   );
        genmu2.SetPtEtaPhiM(genparticles->at(indexme0mu[1]).pt(),genparticles->at(indexme0mu[1]).eta(), genparticles->at(indexme0mu[1]).phi(), 0.105   );
        ZCandReco = recomu1+recomu2;
        ZCandGen = genmu1+genmu2;
        std::cout<<"TipoB: zcand reco  "<< ZCandReco.M()<<" ZCand gen"<<ZCandGen.M()<<endl;
        hFillZGenMass->Fill((ZCandGen).M());
        hFillZMass->Fill((ZCandReco).M());

        if((OurMuons->at(me0muons.at(0)).me0segment().time() >= timeMin_  && OurMuons->at(me0muons.at(0)).me0segment().time() <= timeMax_) && ( OurMuons->at(me0muons.at(1)).me0segment().time() >= timeMin_  && OurMuons->at(me0muons.at(1)).me0segment().time() <= timeMax_ )  ){
            
            hFillZMassInWindow->Fill((ZCandReco).M());
            
        } else {
            
            hFillZMassOutWindow->Fill((ZCandReco).M());
            
        }
	}//EVENTO B

        bool allme0mu = false;
	for (unsigned int t = 0; t < OurMuons->size() ; t++){
        
        
        allme0mu = true;      
  
        hFillAllRecoEta->Fill(OurMuons->at(t).eta());
        hFillAllRecoPt->Fill(OurMuons->at(t).pt());

         bool isintime = false;
         if(OurMuons->at(t).me0segment().time() >= timeMin_  && OurMuons->at(t).me0segment().time() <= timeMax_) isintime = true;
         bool decisionTime = (TD_ == "true") ? isintime :  allme0mu;

         if(!decisionTime)continue;
        //if(!decisionTime) std::cout << "reco muon in time " << std::endl;
         

        bool isLoose = muon::isGoodMuon(me0geom, OurMuons->at(t), muon::Loose);
        bool isTight = muon::isGoodMuon(me0geom, OurMuons->at(t), muon::Tight);
        bool decision = (wp_ == "loose") ? isLoose : isTight;
        bool decisionN = (!decision);
        bool decisionF = (fD_ == "ID") ? decision : decisionN;

   
        if (!decisionF) continue;
        if(!(OurMuons->at(t).pt() > ptMin_))continue;
        //if (OurMuons->at(t).pt() < ptMinGen_) continue;
             
        hFillRecoEta->Fill(OurMuons->at(t).eta());
        hFillRecoPt->Fill(OurMuons->at(t).pt());
        hFillRecoMuonTime->Fill(OurMuons->at(t).me0segment().time());
        hFillRecoMuonTimeErr->Fill(OurMuons->at(t).me0segment().timeErr());     

        if(MatchedMuon(me0muons, int(t))) {
           std::cout << "Matched reco muon eta ==> " << OurMuons->at(t).eta() << "  Matched Signal muon pt ==> " << OurMuons->at(t).pt() << "  Matched Signal muon phi ==> " << OurMuons->at(t).phi() << std::endl;
          }
   
        if(!MatchedMuon(me0muons, int(t))) {
            
            hFillPUMuontime->Fill(OurMuons->at(t).me0segment().time());
            hFillPUMuonPt->Fill(OurMuons->at(t).pt());
            hFillPUMuonEta->Fill(fabs(OurMuons->at(t).eta()));
            hFillPUMuonPt_rebin->Fill(OurMuons->at(t).pt());
            hFillPUMuontimeErr->Fill(OurMuons->at(t).me0segment().timeErr());
            BGMuonTime->Fill(OurMuons->at(t).me0segment().time(),OurMuons->at(t).me0segment().timeErr());
            
            if(OurMuons->at(t).me0segment().time() >= timeMin_ && OurMuons->at(t).me0segment().time() <= timeMax_){
            
                hFillPUMuonPtInWindow->Fill(OurMuons->at(t).pt());
                hFillPUMuonPtInWindow_rebin->Fill(OurMuons->at(t).pt());
                hFillPUMuonEtaInWindow->Fill(fabs(OurMuons->at(t).eta()));
                
            }
            else {
                
                hFillPUMuonPtOutWindow->Fill(OurMuons->at(t).pt());
                hFillPUMuonPtOutWindow_rebin->Fill(OurMuons->at(t).pt());
                hFillPUMuonEtaOutWindow->Fill(fabs(OurMuons->at(t).eta()));
                
            }
            
		}
        
    }

}


// ------------ method called once each job just before starting event loop  ------------
	void 
ME0TimingAnalysis::beginJob()
{
    
    hFillEvents = fs->make<TH1F>("hFillEvents","Number of events",10,0,10);
	hFillSignalMuontime = fs->make<TH1F>("hFillSignalMuontime","hFillSignalMuontime",1000,-300,300);  
	hFillPUMuontime = fs->make<TH1F>("hFillPUMuontime","hFillPUMuontime",1000,-300,300);
	hFillSignalMuontimeErr = fs->make<TH1F>("hFillSignalMuontimeErr","hFillSignalMuontimeErr",1000,0,10);  
	hFillPUMuontimeErr = fs->make<TH1F>("hFillPUMuontimeErr","hFillPUMuontimeErr",1000,0,10);

	SignalMuonTime  = fs->make<TH1F>("SignalMuonTime","SignalMuonTime",1000,-300,300);  
	BGMuonTime  = fs->make<TH1F>("BGMuonTime","BGMuonTime",1000,-300,300);  
	hFillZMass = fs->make<TH1F>("hFillZMass","hFillZMass",500,0,250); 
	hFillRecoEta = fs->make<TH1F>("hFillRecoEta","hFillRecoEta",500,-5,5);
    hFillRecoPt = fs->make<TH1F>("hFillRecoPt","hFillRecoPt",500,0,300);
    hFillAllRecoEta = fs->make<TH1F>("hFillAllRecoEta","hFillRecoEta",500,-5,5);
    hFillAllRecoPt = fs->make<TH1F>("hFillAllRecoPt","hFillRecoPt",500,0,300);
    hFillGenMuEta = fs->make<TH1F>("hFillGenMuEta","hFillGenMuEta",250,0,5);
    hFillGenMuPt = fs->make<TH1F>("hFillGenMuPt","hFillGenMuPt",500,0,300);
    hFillZGenMass = fs->make<TH1F>("hFillZGenMass","hFillZGenMass",500,0,250);    

   hFillRecoMuonTime = fs->make<TH1F>("hFillRecoMuonTime","hFillRecoMuonTime",1000,-300,300);
   hFillRecoMuonTimeErr = fs->make<TH1F>("hFillRecoMuonTimeErr","hFillRecoMuonTimeErr",1000,0,10);

    hFillGenMuEta_allEta = fs->make<TH1F>("hFillGenMuEta_allEta","hFillGenMuEta_allEta",250,0,5);
    hFillGenMuPt_allEta = fs->make<TH1F>("hFillGenMuPt_allEta","hFillGenMuPt_allEta",500,0,300);
 
    hFillGenMuEta_me0Eta = fs->make<TH1F>("hFillGenMuEta_me0Eta","hFillGenMuEta_meoEta",250,0,5);
    hFillGenMuPt_me0Eta = fs->make<TH1F>("hFillGenMuPt_me0Eta","hFillGenMuPt_me0Eta",500,0,300);

    hFillGenMuinME0Eta = fs->make<TH1F>("hFillGenMuinME0Eta","hFillGenMuinME0Eta",250,0,5);
    hFillGenMuinME0Pt = fs->make<TH1F>("hFillGenMuinME0Pt","hFillGenMuinME0Pt",500,0,300);

	hFillZMassInWindow = fs->make<TH1F>("hFillZMassInWindow","hFillZMassInWindow",500,0,250);
	hFillZMassOutWindow = fs->make<TH1F>("hFillZMassOutWindow","hFillZMassOutWindow",500,0,250);
	hFillPUMuonPt = fs->make<TH1F>("hFillPUMuonPt","hFillPUMuonPt",500,0,300);
    hFillPUMuonPtInWindow = fs->make<TH1F>("hFillPUMuonPtInWindow","hFillPUMuonPtInWindow",500,0,300);
    hFillPUMuonPtOutWindow = fs->make<TH1F>("hFillPUMuonPtOutWindow","hFillPUMuonPtOutWindow",500,0,300);

    hFillPUMuonEta = fs->make<TH1F>("hFillPUMuonEta","hFillPUMuonEta",250,0,5);
    hFillPUMuonEtaInWindow = fs->make<TH1F>("hFillPUMuonEtaInWindow","hFillPUMuonEtaInWindow",250,0,5);
    hFillPUMuonEtaOutWindow = fs->make<TH1F>("hFillPUMuonEtaOutWindow","hFillPUMuonEtaOutWindow",250,0,5);

    //Float_t bins[] = {0,1,3,5,10,20,40,60};
    Float_t bins[] = {0,5,10,20,40,70};
    hFillPUMuonPt_rebin = fs->make<TH1F>("hFillPUMuonPt_rebin","hFillPUMuonPt_rebin",5,bins);
    hFillPUMuonPtInWindow_rebin = fs->make<TH1F>("hFillPUMuonPtInWindow_rebin","hFillPUMuonPtInWindow_rebin",5,bins);
    hFillPUMuonPtOutWindow_rebin = fs->make<TH1F>("hFillPUMuonPtOutWindow_rebin","hFillPUMuonPtOutWindow_rebin",5,bins);

}


// ------------ method called once each job just after ending the event loop  ------------
	void 
ME0TimingAnalysis::endJob() 
{
}

void
ME0TimingAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

bool ME0TimingAnalysis::MatchedMuon(vector<int> me0muons, int recomuon) {
	bool ismatched = false;
	for(unsigned int sd = 0; sd < me0muons.size(); sd++) {
		if(int(recomuon) == int(me0muons.at(sd))) {
			ismatched = true;
			break;
		}
	}
	return ismatched;
}
//define this as a plug-in
DEFINE_FWK_MODULE(ME0TimingAnalysis);
