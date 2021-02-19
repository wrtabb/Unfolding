#include "include/ToyModel.hh"
const TString saveName = "data/toyModelDistributions.root";

void makeToyModels()
{
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	
	//number of true bins
	int nBins = 50;

	double norm = 100;
	double expDecay = 20;
	double peakNormRel = 10;
	double mean = 60;
	double sigma = 5;
	double resSigma = 3;
	double xmin = 0;
	double xmax = 100;

	ToyModel*model = new ToyModel(norm,expDecay,peakNormRel,mean,sigma,resSigma,xmin,xmax,
				      nBins);
	TH1F*hTrue = model->GetTrueHist();
	TH1F*hRecoInversion = model->GetRecoHist(false);
	TH1F*hRecoTUnfold   = model->GetRecoHist(true);
	TH2F*hMatrixInversion = model->GetMigrationMatrix(false);
	TH2F*hMatrixTUnfold = model->GetMigrationMatrix(true);

	TFile*saveFile = new TFile(saveName,"recreate");
	hTrue->Write();
	hRecoInversion->Write();
	hRecoTUnfold->Write();
	hMatrixInversion->Write();
	hMatrixTUnfold->Write();
	saveFile->Close();
}
