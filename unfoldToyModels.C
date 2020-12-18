#include "include/Unfolding.hh"

void unfoldToyModels()
{
	const TString fileLocation = "./data/toyModelDistributions.root";
	TFile*file = new TFile(fileLocation);

	TH1D*hReco =     (TH1D*)file->Get("hReco");
	TH1D*hClosure =  (TH1D*)file->Get("hClosure");
	TH1D*hTrue =     (TH1D*)file->Get("hTrue");
	TH2D*hMatrix =   (TH2D*)file->Get("hMatrix");
	TH2D*hResponse = (TH2D*)file->Get("hResponse");

	Unfold*toyUnfold = new Unfold();
	Unfold::RegType regType = toyUnfold->NO_REG;
	Unfold::UnfoldType unfoldType = toyUnfold->TUNFOLD;
	//TH1F*hTUnfoldClosure = toyUnfold->unfoldTUnfold(regType,hClosure,hTrue,hMatrix);
	TH1F*hInversionClosure = toyUnfold->unfoldInversion(hClosure,hTrue,hResponse);
	//toyUnfold->plotUnfolded(hClosure,hTrue,hTUnfoldClosure,unfoldType,true);	
	toyUnfold->plotUnfolded(hClosure,hTrue,hInversionClosure,unfoldType,true);	
}
