#include "include/Unfolding.hh"

TString fileName = "data/toyModelDistributions.root";
TString saveFileName = "data/unfoldedOutput.root";

void calcUnfold()
{
	//Initialize file and get histograms needed for unfolding from it
	TFile*file = new TFile(fileName);
	TH1F*hTrue =   (TH1F*)file->Get("hTrue");	
	TH1F*hRecoInversion =   (TH1F*)file->Get("hRecoInversion");	
	TH1F*hRecoTUnfold =   (TH1F*)file->Get("hRecoTUnfold");	
	TH2F*hMatrix = (TH2F*)file->Get("migrationHist");	

	//Initialize Unfold object
	Unfold*unfold = new Unfold();
	//Set regularization to "no regularization". This should be the default
	//We will later include regularization if we determine we need to
	Unfold::RegType regType = unfold->NO_REG;
	//unfoldType is set for the plotting function only
	//I am rethinking the way this is implemented, so it may change in the future
	Unfold::UnfoldType unfoldType = unfold->TUNFOLD;

	//Do unfolding using TUnfold and place results in a histogram
	TH1F*hUnfoldedTUnfold = unfold->unfoldTUnfold(regType,hRecoTUnfold,hTrue,hMatrix);
	//Do unfolding using inversion method and place results in a histogram
	TH1F*hUnfoldedInversion = unfold->unfoldInversion(hRecoInversion,hTrue,hMatrix);
//	unfold->plotUnfolded(hReco,hTrue,hUnfolded,unfoldType,true);


	//Save all files together in one place
	TFile*saveFile = new TFile(saveFileName,"recreate");
	hTrue->Write();
	hRecoInversion->Write();
	hRecoTUnfold->Write();
	hMatrix->Write();
	hUnfoldedTUnfold->Write();
	hUnfoldedInversion->Write();
	saveFile->Close();
}
