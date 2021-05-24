#include "include/ToyModel.hh"
#include "include/Unfolding.hh"

using namespace GlobalVariables;
using namespace Utilities;

void peakComparisonUnfolding()
{
	TH1::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	//gROOT->SetBatch(true);

	TString loadName = "data/toyModelRecoBin0.root";
	TString loadNameDY = "data/toyModelRecoBin5.root";
	TString saveName = "data/peakComparisonInversion.root";

	TFile*file = new TFile(loadName);
	TH1F*hTrue = (TH1F*)file->Get("hTrueRandom");
	TH1F*hReco = (TH1F*)file->Get("hRecoRandom");

	TFile*fileDY = new TFile(loadNameDY);
	TH1F*hTrueDY = (TH1F*)fileDY->Get("hTrueRandom");
	TH1F*hRecoDY = (TH1F*)fileDY->Get("hRecoRandom");

	vector<double> binTrueDY = _peakMassTrueDY;
	//vector<double> binRecoDY = _peakMassRecoDY;
	vector<double> binRecoDY = binTrueDY;
	vector<double> binTrue = _peakMassTrue;
	//vector<double> binReco = _peakMassReco;
	vector<double> binReco = binTrue;

	const int sizeTrue = binTrue.size();
	double binTrueArray[sizeTrue];
	const int sizeReco = binReco.size();
	double binRecoArray[sizeReco];

	for(int i=0;i<sizeTrue;i++){
		binTrueArray[i] = binTrue.at(i);
	}

	for(int i=0;i<sizeReco;i++){
		binRecoArray[i] = binReco.at(i);
	}

	const int sizeTrueDY = binTrueDY.size();
	double binTrueDYArray[sizeTrueDY];
	const int sizeRecoDY = binRecoDY.size();
	double binRecoDYArray[sizeRecoDY];

	for(int i=0;i<sizeTrueDY;i++){
		binTrueDYArray[i] = binTrueDY.at(i);
	}

	for(int i=0;i<sizeRecoDY;i++){
		binRecoDYArray[i] = binRecoDY.at(i);
	}

	int nBinsTrue = sizeTrue - 1;
	int nBinsReco = sizeReco - 1;
	int nBinsTrueDY = sizeTrueDY - 1;
	int nBinsRecoDY = sizeRecoDY - 1;

	double norm = 1000000;//overall scaling factor
	double peakNormRel = 0.00006;//scaling factor for the peak only
	double mean = 91;//mean of the peak
	double sigma = 10;//spread of the peak
	double shift = 6*sigma;//shift of the asymptote of the power function
	double resSigma = 2.0;

	ToyModel*model   = new ToyModel(norm,shift,peakNormRel,mean,sigma,resSigma,binTrue,
					binReco);	
	ToyModel*modelDY = new ToyModel(norm,shift,peakNormRel,mean,sigma,resSigma,binTrueDY,
					binRecoDY);	

	TH1F*hTrueRandom = new TH1F("hTrueRandom","",nBinsTrue,binTrueArray);
	hTrueRandom->SetFillColor(kRed+2);
	hTrueRandom->SetLineColor(kRed+2);
	TH1F*hRecoRandom = new TH1F("hRecoRandom","",nBinsReco,binRecoArray);
	hRecoRandom->SetMarkerColor(kBlack);
	hRecoRandom->SetMarkerStyle(20);

	TH1F*hTrueRandomDY = new TH1F("hTrueRandomDY","",nBinsTrueDY,binTrueDYArray);
	hTrueRandomDY->SetFillColor(kRed+2);
	hTrueRandomDY->SetLineColor(kRed+2);
	TH1F*hRecoRandomDY = new TH1F("hRecoRandomDY","",nBinsRecoDY,binRecoDYArray);
	hRecoRandomDY->SetMarkerColor(kBlack);
	hRecoRandomDY->SetMarkerStyle(20);

	TH2F*hMatrix = model->GetMigrationMatrix("hMatrix");
	TH2F*hResponse = makeResponseMatrix(hMatrix);
	TH2F*hResponseRebin = RebinTH2(hResponse,"hResponseRebin",hTrueRandom);

	TH2F*hMatrixDY = modelDY->GetMigrationMatrix("hMatrixDY");
	TH2F*hResponseDY = makeResponseMatrix(hMatrixDY);
	TH2F*hResponseDYRebin = RebinTH2(hResponseDY,"hResponseDYRebin",hTrueRandomDY);

	//Let's look at the response matrices
	TCanvas*c1 = PlotMatrix("c1","Response Matrix",hResponseRebin,false);
	c1->SaveAs("plots/responseMatrix.png");
	TCanvas*c2 = PlotMatrix("c2","DY Response Matrix",hResponseDYRebin,false);
	c2->SaveAs("plots/responseMatrixDY.png");

	double contentTrue;
	double contentReco;
	double errorTrue;
	double errorReco;
	int binShiftTrue = 12; 
	for(int i=0;i<=nBinsTrue;i++){
		contentTrue = hTrue->GetBinContent(i+binShiftTrue);
		errorTrue = hTrue->GetBinError(i+binShiftTrue);
		hTrueRandom->SetBinContent(i,contentTrue);
		hTrueRandom->SetBinError(i,errorTrue);
		contentReco = hReco->GetBinContent(i+binShiftTrue);
		errorReco = hReco->GetBinError(i+binShiftTrue);
		hRecoRandom->SetBinContent(i,contentReco);
		hRecoRandom->SetBinError(i,errorReco);
	}
/*
	int binShiftReco = 24; 
	for(int i=0;i<=nBinsReco;i++){
		content = hReco->GetBinContent(i+binShiftReco);
		error = hReco->GetBinError(i+binShiftReco);
		hRecoRandom->SetBinContent(i,content);
		hRecoRandom->SetBinError(i,error);
	}
*/

	int binShiftTrueDY = 9; 
	for(int i=0;i<=nBinsTrueDY;i++){
		contentTrue = hTrueDY->GetBinContent(i+binShiftTrueDY);
		errorTrue = hTrueDY->GetBinError(i+binShiftTrueDY);
		hTrueRandomDY->SetBinContent(i,contentTrue);
		hTrueRandomDY->SetBinError(i,errorTrue);
		contentReco = hRecoDY->GetBinContent(i+binShiftTrueDY);
		errorReco = hRecoDY->GetBinError(i+binShiftTrueDY);
		hRecoRandomDY->SetBinContent(i,contentReco);
		hRecoRandomDY->SetBinError(i,errorReco);
	}
/*
	int binShiftRecoDY = 18; 
	for(int i=0;i<=nBinsRecoDY;i++){
		content = hRecoDY->GetBinContent(i+binShiftRecoDY);
		error = hRecoDY->GetBinError(i+binShiftRecoDY);
		hRecoRandomDY->SetBinContent(i,content);
		hRecoRandomDY->SetBinError(i,error);
	}
*/
	Unfold*unfold = new Unfold();
	Unfold::RegType regType = unfold->NO_REG;
	Unfold*unfoldDY = new Unfold();
	Unfold::RegType regTypeDY = unfold->NO_REG;
	cout << "Condition number: " << unfold->GetConditionNumber(hResponseRebin) << endl;
	cout << "Condition number: " << unfoldDY->GetConditionNumber(hResponseDYRebin) << endl;

	//TH1F*hUnfolded = unfold->unfoldTUnfold(regType,hRecoRandom,hTrueRandom,hMatrix);
	//TH1F*hUnfoldedDY = unfoldDY->unfoldTUnfold(regTypeDY,hRecoRandomDY,hTrueRandomDY,
	TH1F*hUnfolded = unfold->unfoldInversion(hRecoRandom,hTrueRandom,hMatrix);
	TH1F*hUnfoldedDY = unfoldDY->unfoldInversion(hRecoRandomDY,hTrueRandomDY,
						   hMatrixDY);

	TCanvas*c3 = unfold->plotUnfolded("c3","Unfolded results",hRecoRandom,hTrueRandom,
					  hUnfolded,false);
	TCanvas*c4 = unfoldDY->plotUnfolded("c4","DY Unfolded results",hRecoRandomDY,
					  hTrueRandomDY,hUnfoldedDY,false);
	c3->SaveAs("plots/unfoldedResultsInversion.png");
	c4->SaveAs("plots/unfoldedResultsDYInversion.png");
	
	hUnfolded->SetName("hUnfolded");
	hUnfoldedDY->SetName("hUnfoldedDY");
	TFile*saveFile = new TFile(saveName,"recreate");
	hTrueRandom->Write();
	hRecoRandom->Write();
	hMatrix->Write();
	hTrueRandomDY->Write();
	hRecoRandomDY->Write();
	hMatrixDY->Write();
	hUnfolded->Write();
	hUnfoldedDY->Write();
	saveFile->Close();
}


