#include "include/ToyModel.hh"
#include "include/Unfolding.hh"

const TString saveName = "data/toyModelDistributions.root";
TCanvas*PlotProjections(TString canvasName,TH2F*hMatrix,TH1F*hTrue,TH1F*hReco);
TCanvas*PlotMatrix(TString canvasName,TH2F*hist);

void makeToyModels()
{
	TH1::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	
	//number of true bins
	int nBins = 30;

	double norm = 1000;
	double power = 1;
	double peakNormRel = 0.3;
	double mean = 91;
	double sigma = 15;
	double resSigma = 2;
	double xmin = 0;
	double xmax = 180;

	ToyModel*model = new ToyModel(norm,power,peakNormRel,mean,sigma,resSigma,xmin,xmax,
				      nBins);
	TH1F*hTrue = model->GetTrueHist();
	TH1F*hRecoInversion = model->GetRecoHist(false);
	TH1F*hRecoTUnfold   = model->GetRecoHist(true);
	TH2F*hMatrixInversion = model->GetMigrationMatrix(false);
	TH2F*hMatrixTUnfold = model->GetMigrationMatrix(true);
	TH2F*hResponseInversion = model->GetResponseMatrix(hMatrixInversion);
	TH2F*hResponseTUnfold = model->GetResponseMatrix(hMatrixTUnfold);

	TH2F*hMatrixCompare = (TH2F*)hMatrixInversion->Clone();
	hMatrixCompare->Add(hResponseInversion,-1);

	TCanvas*c1 = PlotMatrix("c1",hResponseTUnfold);
	TCanvas*c2 = PlotMatrix("c2",hMatrixTUnfold);
	TCanvas*c4 = PlotProjections("c4",hMatrixTUnfold,hTrue,hRecoTUnfold);

	Unfold*unfold = new Unfold();
	Unfold::RegType regType = unfold->NO_REG;
	TH1F*hUnfoldedTUnfold = unfold->unfoldTUnfold(regType,hRecoTUnfold,hTrue,
						      hMatrixTUnfold);

	TFile*saveFile = new TFile(saveName,"recreate");
	hTrue->Write();
	hRecoInversion->Write();
	hRecoTUnfold->Write();
	hMatrixInversion->Write();
	hMatrixTUnfold->Write();
	hUnfoldedTUnfold->Write();
	saveFile->Close();

	TCanvas*c5 = unfold->plotUnfolded("c5",hRecoTUnfold,hTrue,hUnfoldedTUnfold);
}

TCanvas*PlotMatrix(TString canvasName,TH2F*hist)
{
	TCanvas*c1 = new TCanvas(canvasName,"",0,0,1000,1000);
	c1->SetGrid();
	hist->Draw("colz");	
	return c1;
}

TCanvas*PlotProjections(TString canvasName,TH2F*hMatrix,TH1F*hTrue,TH1F*hReco)
{
	int nBinsTrue = hTrue->GetNbinsX();
	int nBinsReco = hReco->GetNbinsX();
	TH1F*hRecoRebin = (TH1F*)hReco->Clone("hRecoRebin");
	TH2F*hMatrixRebin = (TH2F*)hMatrix->Clone("hMatrixRebin");
	if(nBinsReco!=nBinsTrue){ 
		//this loop assumes that nBinsReco = 2*nBinsTrue
		//This should be made more general with the use of my
		//custom rebinning function
		//Need to implement it in this package
		hRecoRebin->Rebin(2);
		hMatrixRebin->RebinX(2);
	}
	TCanvas*canvas = new TCanvas(canvasName,"",0,0,1000,1000);
	canvas->SetGrid();
	TH1F*projX = (TH1F*)hMatrixRebin->ProjectionX();
	TH1F*projY = (TH1F*)hMatrixRebin->ProjectionY();
	projX->Scale(hReco->Integral()/projX->Integral());
	projX->SetMarkerStyle(20);
	projY->SetMarkerStyle(20);
	projX->SetMarkerColor(kBlue);
	projY->SetMarkerColor(kRed);
	projX->SetLineColor(kBlue);
	projY->SetLineColor(kRed);
	hTrue->SetLineColor(kRed);
	hRecoRebin->SetLineColor(kBlue);
	hTrue->SetFillColor(kWhite);
	hRecoRebin->SetFillColor(kWhite);

	hTrue->GetYaxis()->SetRangeUser(0,150);
	hTrue->Draw("hist");
	hRecoRebin->Draw("hist,same");
	projX->Draw("pe,same");
	projY->Draw("pe,same");
	return canvas;
}

