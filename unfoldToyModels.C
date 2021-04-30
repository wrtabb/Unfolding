#include "include/ToyModel.hh"
#include "include/Unfolding.hh"


void unfoldToyModels(int binningType)
{
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	TString saveFile = "./data/toyModelRecoBin";
	saveFile += binningType;
	saveFile += ".root";

	TFile*file = new TFile(saveFile);
	TH1F*hTrue = (TH1F*)file->Get("hTrueRandom");
	TH1F*hReco = (TH1F*)file->Get("hRecoRandom");
	TH2F*hMatrix = (TH2F*)file->Get("hMatrix");

	Unfold*unfold = new Unfold();
	Unfold::RegType regType = unfold->NO_REG;

	TH1F*hUnfolded = unfold->unfoldTUnfold(regType,hReco,hTrue,hMatrix);
	hUnfolded->SetMarkerStyle(25);
	hUnfolded->SetMarkerColor(kBlue+2);
	hUnfolded->SetLineColor(kBlue+2);

	TCanvas*c3 = unfold->plotUnfolded("c3","Unfolded results",hReco,hTrue,hUnfolded);
	TString saveName = "plots/unfoldedRecoBin";
	saveName += binningType;
	saveName += ".png";
	c3->SaveAs(saveName);
}
