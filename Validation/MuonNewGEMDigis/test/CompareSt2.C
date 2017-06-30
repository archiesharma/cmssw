#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TList.h"
#include "TGaxis.h"
#include "TMath.h"
#include <stdio.h>
#include <TROOT.h>
#include "TSystem.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>


void savePlot(TCanvas*, TString);
void prepareHisto (TH1F* hist, Int_t rebin, Double_t weight, Color_t color);


void CompareSt2()
{

  TH1::SetDefaultSumw2();
  cout<<"ciao"<<endl;
  TFile *file = TFile::Open("testGEMDigiReader_NB.root");
  file->cd("/dumper");
  TH1F *BkgNB = (TH1F*)gDirectory->Get("hglobalR_Sta2_lay2");
  
  TFile *file1 = TFile::Open("testGEMDigiReader_NoNB.root");
  file1->cd("/dumper");
  TH1F *BkgNoNB = (TH1F*)gDirectory->Get("hglobalR_Sta2_lay2");

 
  
  Double_t weight = 8900*3*25e-09;
  Double_t weight2 = 8450*3*25e-09;

  prepareHisto(BkgNB,1,weight,kRed+1);
  prepareHisto(BkgNoNB,1,weight2,kBlue+1);

  BkgNB->Add(BkgNoNB,-1);
  BkgNB->Sumw2(); 

  cout<<"ciao"<<endl;

 
  
  TCanvas *c = new TCanvas("BkgHitRate", "BkgHitRate", 600, 600);
  c->cd();
  c->SetTicks(1,1);
  c->SetGrid(1,1);
  c->SetLogy();
  c->SetTitle("GEM St1 Neutron Background -long chambers");
  // c->SetTitle("");

    
  TF1 *f1 = new TF1("totalFunc", "(97.0505+1440.44)-(0.452612+7.48607)*x+(0.000550599+0.0103078)*x*x", 130., 320.);
   
  //f1->GetYaxis()->SetRangeUser(1.01,1.2E5);
  f1->GetYaxis()->SetTitleOffset(1.35);
  f1->GetYaxis()->SetTitle("Rate (Hz/cm^{2})");
  f1->GetXaxis()->SetTitle("R (cm)");
  c->Update();
  //savePlot(c,"BkgParamSLHC26p3");

  TCanvas* c_all = new TCanvas("Superposition","Superposition",600,600);
  c_all->cd();
  c_all->SetTicks(1,1);
  c_all->SetGrid(1,1);
  c_all->SetLogy();
  c_all->SetTitle("GEM St2 Neutron Background Parameterization");
  f1->Draw("FC");
  BkgNB->Draw("9p1e0same");
  
   TLegend *leg2 = new TLegend(0.60,0.60,0.85,0.85);
  leg2->SetFillColor(0);
  leg2->SetHeader("GEM St 2 Neutron Bkg Param and Hits - Layer1");
  leg2->AddEntry("totalFunc","n/#gamma Bkg Model","l");
  leg2->AddEntry(BkgNB,"From DIGI Model in CMSSW","pl");
  leg2->Draw();
  c_all->Update();


}


void prepareHisto (TH1F* hist, Int_t rebin, Double_t weight, Color_t color){

  hist->Rebin(rebin);
  //std::cout<<" BinWidth dopo Rebin: "<<hist->GetBinContent(30)<<std::endl;
  //hist->Scale(1/weight);
  hist->SetStats(kFALSE);
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(kFullTriangleDown);
  hist->SetMarkerSize(0.7);
  hist->GetXaxis()->SetRangeUser(130,320);
  std::cout << " Number of bins : "<<hist->GetNbinsX() << std::endl;
  std::cout << " Histogram :  "<<hist->GetName() << " -------------  ";
  for (Int_t bin = 1; bin<=hist->GetNbinsX();++bin){
    //  if (bin == 1) minBase = bottomLenght;
    //  else minBase = maxBase;
    //  if (bin==hist->GetNbinsX()) maxBase = topLenght;

      Double_t R_min = hist->GetBinCenter(bin) - 0.5 * hist->GetBinWidth(bin);
      Double_t R_max = hist->GetBinCenter(bin) + 0.5 * hist->GetBinWidth(bin);
    //  Double_t DeltaR = R_max - R_min;
    //  maxBase = 2*DeltaR*0.1763269 + minBase ;

    //  Double_t trapArea = (0.5*DeltaR*(minBase + maxBase))*18*6*2; //18*6*2 fattori geometrici per tenere conto dei 6 layers, 18 camere e le 2 regioni
           
      Double_t Area     = (TMath::Pi() * ( TMath::Power(R_max,2) - TMath::Power(R_min,2) ))*4;
     
      
      std::cout << "\n #bin = "<<bin
		<<"\n\tArea = "<<Area <<"\n\tWeight = "<<weight<<"  Weight = "<< Area*weight
                <<"\n\tBinContent = "<<hist->GetBinContent(bin)<<"  SetBinContent = "<<(hist->GetBinContent(bin))/(Area*weight)<<std::endl;
      hist->SetBinContent(bin,(hist->GetBinContent(bin))/(Area*weight));
  }

}
