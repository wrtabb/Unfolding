#include "include/ToyModel.hh"
#include "include/Unfolding.hh"

const TString loadName = "data/toyModel.root";
TCanvas*PlotProjections(TString canvasName,TH2F*hMatrix,TH1F*hTrue,TH1F*hReco);
TCanvas*PlotMatrix(TString canvasName,TString plotTitle,TH2F*hist);

void unfoldToyModels(int nBins,double resSigma)
{
	TH1::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	gROOT->SetBatch(true);
	
	TString saveNameTag = "_nBins_";
	saveNameTag += nBins;
	saveNameTag += "_ResSmear_";
	saveNameTag += resSigma;
	saveNameTag += ".png";
	
	TFile*file = new TFile(loadName);
	TH1F*hTrue = (TH1F*)file->Get("hTrue");;
	TH1F*hReco = (TH1F*)file->Get("hReco");;
	TH1F*hTrueRandom = (TH1F*)file->Get("hTrueRandom");
	TH1F*hRecoRandom = (TH1F*)file->Get("hRecoRandom");
	TH2F*hMatrix = (TH2F*)file->Get("hMatrix");

	//Create unfolding object
	Unfold*unfold = new Unfold();
	//Set regularization type
	Unfold::RegType regType = unfold->NO_REG;
	//Carry out unfolding closure test and place in a histogram
	TH1F*hUnfoldedClosure = unfold->unfoldTUnfold(regType,hReco,hTrue,hMatrix);
	TH1F*hUnfolded = unfold->unfoldTUnfold(regType,hRecoRandom,hTrueRandom,hMatrix);

	//plot unfolded results
	TCanvas*c5 = unfold->plotUnfolded("c5","TUnfold unfolded (closure)",hReco,hTrue,
					  hUnfoldedClosure);
	c5->SaveAs("plots/unfoldedClosure"+saveNameTag);
	TCanvas*c6 = unfold->plotUnfolded("c6","TUnfold unfolded",hRecoRandom,hTrueRandom,
					  hUnfolded);
	c6->SaveAs("plots/unfolded"+saveNameTag);

	//Now with inversion method
	TH1F*hUnfoldedInversion = unfold->unfoldInversion(hRecoRandom,hTrueRandom,hMatrix);
	TH1F*hUnfoldedInversionClosure = unfold->unfoldInversion(hReco,hTrue,hMatrix);

	//plot unfolded results
	TCanvas*c7 = unfold->plotUnfolded("c7","Inversion unfolded (closure)",hReco,hTrue,
					  hUnfoldedInversionClosure);
	c7->SaveAs("plots/unfoldedInversionClosure"+saveNameTag);
	TCanvas*c8 = unfold->plotUnfolded("c8","Inversion unfolded",hRecoRandom,hTrueRandom,
					  hUnfoldedInversion);
	c8->SaveAs("plots/unfoldedInversion"+saveNameTag);
}

