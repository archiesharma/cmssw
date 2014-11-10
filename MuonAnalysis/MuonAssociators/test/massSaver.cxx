/**
* filename is the name of the root file
* histopath is the path within the root document to the directory, containing
*            all the histograms
*
* Execution example:
* [] .x massSaver.cxx("test.root", "someHistoDir") DQMData/Muons/MuonSegEff/Residuals/run2014_RE12res_cls123
//run2014_RE12res_Angle
run2014_RE12res_Chamber
run2014_RE12res_Sector
run2014_RE12res_cls123
run2014_RE4res
run2014_xyViews	

*/
void massSaver(const char* fileName, const char* histoPath)
{
	TFile *file = TFile::Open(fileName);
	gDirectory->cd(histoPath);
	TIter iter(gDirectory->GetListOfKeys());
	TKey *key;
	TClass *cl;
	TCanvas *c = new TCanvas("c", "c",1200, 600);
//        c->SetGridx();
//        c->SetGridy();
	c->cd();
	while((key = (TKey*) iter.Next()))
	{
		cl = gROOT->GetClass(key->GetClassName());
		if (cl->InheritsFrom("TH1"))
		{
			TH1 *h = (TH1*)key->ReadObj();
//			std::string outName = std::string(key->GetName()) + ".C";
//			h->SaveAs(outName.c_str());
			h->Draw();
                        h->SetDrawOption("colztext");
			std::string outName = std::string(key->GetName()) + ".png";
			c->SaveAs(outName.c_str());
			delete h;		
		}
	}
	file->Close();
	delete c;
	delete cl;
	delete key;
	delete file;
}
