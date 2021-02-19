#include "include/ToyModel.hh"
const TString saveName = "data/toyModelDistributions.root";

void makeToyModels()
{
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	int nBins = 100;

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
	TH2F*hMatrix = model->GetMigrationMatrix(false);

	TFile*saveFile = new TFile(saveName,"recreate");
	hTrue->Write();
	hRecoInversion->Write();
	hRecoTUnfold->Write();
	hMatrix->Write();
}
