#include "include/ToyModel.hh"
#include "include/Unfolding.hh"

const TString loadMatrix = "data/toyModel.root";
const TString loadDists  = "data/toyModelRandomSeeds.root";
const TString saveLocation = "data/ratioPlotsForRandomSeeds.root";

void unfoldToyModelsRandomSeeds(int nSeeds)
{
	TH1::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	gROOT->SetBatch(true);
	
	TString saveNameTag = "_nBins_50.png";
	
	TFile*fileMatrix = new TFile(loadMatrix);
	TH2F*hMatrix = (TH2F*)fileMatrix->Get("hMatrix");

	TFile*fileDists = new TFile(loadDists);
	TH1F*hTrue[nSeeds];
	TH1F*hReco[nSeeds];
	TString hTrueLoad;
	TString hRecoLoad;
	for(int i=0;i<nSeeds;i++){
		hTrueLoad = "histTrue";
		hTrueLoad += i;
		hRecoLoad = "histReco";
		hRecoLoad += i;

		hTrue[i]=(TH1F*)fileDists->Get(hTrueLoad);
		hReco[i]=(TH1F*)fileDists->Get(hRecoLoad);
	}

	//Create unfolding object
	Unfold*unfold = new Unfold();
	//Set regularization type
	Unfold::RegType regType = unfold->NO_REG;

	//Unfold the distributions
	TH1F*hUnfolded[nSeeds];
	TH1F*hRatio[nSeeds];
	TCanvas*canvas[nSeeds];
	TFile*saveFile = new TFile(saveLocation,"recreate");
	for(int i=0;i<nSeeds;i++){
		TString cName = "canvas";
		cName += i;
		TString ratioName = "ratio";
		ratioName += i;
		TString saveName = "plots/randomSeedUnfold";
		saveName += i;
		saveName += saveNameTag;
		hUnfolded[i] = unfold->unfoldTUnfold(regType,hReco[i],hTrue[i],hMatrix);
		hRatio[i] = (TH1F*)hUnfolded[i]->Clone(ratioName);
		hRatio[i]->Divide(hTrue[i]);
		canvas[i] = unfold->plotUnfolded(cName,"",hReco[i],hTrue[i],hUnfolded[i]);
		canvas[i]->SaveAs(saveName);
		saveFile->cd();
		hRatio[i]->Write();
	}
	saveFile->Close();
}

