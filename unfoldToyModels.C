#include "include/ToyModel.hh"
#include "include/Unfolding.hh"

using namespace GlobalVariables;

void unfoldToyModels(int binningType,bool closure)
{
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	//gROOT->SetBatch(true);
	
	TString saveFile = "data/toyModelRecoBin";
	saveFile += binningType;
	saveFile += ".root";

	cout << "***********************************" << endl;
	cout << "Opening file: " << saveFile << endl;
	cout << "***********************************" << endl;

	TFile*file = new TFile(saveFile,"update");
	TH1F*hTrue;
	TH1F*hReco;
	if(closure){
		hTrue = (TH1F*)file->Get("hTrue");
		hReco = (TH1F*)file->Get("hReco");
	}
	else{
		hTrue = (TH1F*)file->Get("hTrueRandom");
		hReco = (TH1F*)file->Get("hRecoRandom");
	}
	TH2F*hMatrix = (TH2F*)file->Get("hMatrix");

	Unfold*unfold = new Unfold();
	Unfold::RegType regType = unfold->NO_REG;

	TH1F*hUnfolded;
	if(binningType<100 || binningType>199){
		hUnfolded = unfold->unfoldTUnfold(regType,hReco,hTrue,hMatrix);
		hUnfolded->SetMarkerStyle(25);
		hUnfolded->SetMarkerColor(kBlue+2);
		hUnfolded->SetLineColor(kBlue+2);
	}
	else {
		hUnfolded = unfold->unfoldInversion(hReco,hTrue,hMatrix);
		hUnfolded->SetMarkerStyle(25);
		hUnfolded->SetMarkerColor(kBlue+2);
		hUnfolded->SetLineColor(kBlue+2);
	}

	bool logPlot = false;
	TCanvas*c3 = unfold->plotUnfolded("c3","Unfolded results",hReco,hTrue,hUnfolded,logPlot);
	TString saveName = "plots/unfoldedRecoBin";
	saveName += binningType;
	if(closure) saveName+= "_Closure";
	saveName += ".png";
	c3->SaveAs(saveName);
	
	file->cd();
	hUnfolded->Write();
	file->Close();	
}
