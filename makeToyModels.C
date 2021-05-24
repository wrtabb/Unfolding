#include "include/ToyModel.hh"
#include "include/Unfolding.hh"
#include "include/BinningModels.h"

using namespace Utilities;

TCanvas*PlotProjections(TString canvasName,TH2F*hMatrix,TH1F*hTrue,TH1F*hReco,bool logPlot);

// for binning types, see include/GlobalVariables.h
void makeToyModels(int binningType)
{
	TH1::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	//gROOT->SetBatch(true);

	TString saveName = "data/toyModelRecoBin";
	saveName += binningType;
	saveName += ".root";
	
	vector<double> binTrue;
	vector<double> binReco;

	// Range 15-3000
	if(binningType==0){
		binReco = binningmodels::_massbinningReco0;
		binTrue = binningmodels::_massbinningTrue;
	}
	else if(binningType==1){
		binReco = binningmodels::_massbinningReco1;
		binTrue = binningmodels::_massbinningTrue;
	}
	else if(binningType==2){
		binReco = binningmodels::_massbinningReco2;
		binTrue = binningmodels::_massbinningTrue;
	}
	else if(binningType==3){
		binReco = binningmodels::_massbinningReco3;
		binTrue = binningmodels::_massbinningTrue;
	}
	// Range 15-3000, inversion method
	else if(binningType==100){
		binReco = binningmodels::_massbinningTrue;
		binTrue = binningmodels::_massbinningTrue;
	}
	// Range 60-120
	else if(binningType==200){
		binReco = peakmodels::_peakbinningRecoDY;
		binTrue = peakmodels::_peakbinningTrueDY;
	}
	else if(binningType==201){
		binReco = peakmodels::_peakbinningRecoDY_SplitLow;
		binTrue = peakmodels::_peakbinningTrueDY;
	}
	else if(binningType==202){
		binReco = peakmodels::_peakbinningRecoDY_SplitHigh;
		binTrue = peakmodels::_peakbinningTrueDY;
	}
	else if(binningType==203){
		binReco = peakmodels::_peakbinningRecoDY_SplitLowAndHigh;
		binTrue = peakmodels::_peakbinningTrueDY;
	}
	else if(binningType==204){
		binReco = peakmodels::_peakbinningReco5GeV;
		binTrue = peakmodels::_peakbinningTrue5GeV;
	}
	else if(binningType==205){
		binReco = peakmodels::_peakbinningReco4GeV;
		binTrue = peakmodels::_peakbinningTrue4GeV;
	}

	else {
		cout << "binningType = " << binningType << " does not exist" << endl;
		cout << "Please see include/GlobalVariables.h for a list of binning types" << endl;
		return;
	}

	const int sizeTrue = binTrue.size();
	cout << "True size: " << sizeTrue << endl;
	double binTrueArray[sizeTrue];

	const int sizeReco = binReco.size();
	cout << "Reco size = " << sizeReco << endl;
	double binRecoArray[sizeReco];

	for(int i=0;i<sizeTrue;i++){
		binTrueArray[i] = binTrue.at(i);
	}

	for(int i=0;i<sizeReco;i++){
		binRecoArray[i] = binReco.at(i);
	}

	int nBinsTrue = sizeTrue - 1;
	int nBinsReco = sizeReco - 1;

	double norm = 1000000;//overall scaling factor
	double peakNormRel = 0.00001;//scaling factor for the peak only
	double mean = 91;//mean of the peak
	double sigma = 3.8;//spread of the peak
	double resSigma = 2.0;
	double shift = 10*resSigma;//shift of the asymptote of the power function

	//Create model using parameters above
	ToyModel*model = new ToyModel(norm,shift,peakNormRel,mean,sigma,resSigma,binTrue,
				      binReco);
	
	//Define all needed histogams from model
	//At minimum, true,reco, and matrix are required for unfolding
	TH1F*hTrue = model->GetTrueHist("hTrue");
	TH1F*hReco = model->GetRecoHist("hReco");
	TH1F*hRecoRebin = RebinTH1(hReco,"hRecoRebin",hTrue);
	TH2F*hMatrix = model->GetMigrationMatrix("hMatrix");
	TH2F*hResponse = makeResponseMatrix(hMatrix);
	TH2F*hResponseT = makeResponseMatrixT(hMatrix);

	TF1*trueFunction = model->GetTrueFunction();
	TF1*recoFunction = model->GetRecoFunction();

	//By default, the reco binning is finer than true binning
	//This is because TUnfold has this requirement
	//Here we rebin them for plotting
	TH2F*hResponseRebin = RebinTH2(hResponse,"hResponseRebin",hTrue);
	TH2F*hResponseTRebin = RebinTH2(hResponseT,"hResponseTRebin",hTrue);
	TH2F*hMatrixRebin = RebinTH2(hMatrix,"hMatrixRebin",hTrue);

	TString saveNameTag = "_binning";
	saveNameTag += binningType;
	saveNameTag += ".png";

	bool logPlot;
	logPlot = true;
	//Plot matrices and migration matrix projections alongside true and reco
	TCanvas*c1 = PlotMatrix("c1","response matrix",hResponseRebin,logPlot);
	c1->SaveAs("plots/responseMatrix"+saveNameTag);
	//TCanvas*c2 = PlotMatrix("c2","response matrix norm inverted",hResponseTRebin,logPlot);
	//c2->SaveAs("plots/responseMatrixT"+saveNameTag);

	TCanvas*c3 = PlotProjections("c3",hMatrix,hTrue,hReco,true);
	//Now produce randomly filled distributions from the model
	//To test unfolding on distributions which are not identical
	//to those used to make the matrix
	Long64_t nEvents = 1e7;
	double trueContent;
	double recoContent;
	TH1F*hTrueRandom = new TH1F("hTrueRandom","",nBinsTrue,binTrueArray);
	hTrueRandom->SetFillColor(kRed+2);
	hTrueRandom->SetLineColor(kRed+2);
	TH1F*hRecoRandom = new TH1F("hRecoRandom","",nBinsReco,binRecoArray);
	hRecoRandom->SetMarkerColor(kBlack);
	hRecoRandom->SetMarkerStyle(20);

	for(Long64_t i=0;i<nEvents;i++){
		trueContent = trueFunction->GetRandom();
		recoContent = recoFunction->GetRandom();
		hTrueRandom->Fill(trueContent);
		hRecoRandom->Fill(recoContent);	
	}
	
	TFile*saveFile = new TFile(saveName,"recreate");
	hTrue->Write();
	hReco->Write();
	hTrueRandom->Write();
	hRecoRandom->Write();
	hMatrix->Write();
	saveFile->Close();
}

TCanvas*PlotProjections(TString canvasName,TH2F*hMatrix,TH1F*hTrue,TH1F*hReco,bool logPlot)
{
	int nBinsTrue = hTrue->GetNbinsX();
	int nBinsReco = hReco->GetNbinsX();

	TCanvas*canvas = new TCanvas(canvasName,"",0,0,1000,1000);
	canvas->SetGrid();
	if(logPlot){
		canvas->SetLogx();
		canvas->SetLogy();
	}
	TH1F*projX = (TH1F*)hMatrix->ProjectionX();
	TH1F*projY = (TH1F*)hMatrix->ProjectionY();
	projX->Scale(hReco->Integral()/projX->Integral());
	projX->SetMarkerStyle(20);
	projY->SetMarkerStyle(20);
	projX->SetMarkerColor(kBlue);
	projY->SetMarkerColor(kRed);
	projX->SetLineColor(kBlue);
	projY->SetLineColor(kRed);
	TH1F*hTrueClone = (TH1F*)hTrue->Clone();
	hTrueClone->SetLineColor(kRed);
	hReco->SetLineColor(kBlue);
	hTrueClone->SetFillColor(kWhite);
	hReco->SetFillColor(kWhite);

	float maxY = 0.0;
	float binContent;
	for(int i=1;i<=nBinsTrue;i++){
		binContent = hTrueClone->GetBinContent(i);
		if(binContent > maxY) maxY = binContent;
	}
	float maxYRange = maxY*1.2;
	hTrueClone->GetXaxis()->SetTitle("mass [GeV]");
	hTrueClone->GetYaxis()->SetRangeUser(0.01,maxYRange);
	hTrueClone->SetTitle("matrix projections with distributions");
	hTrueClone->Draw("hist");
	hReco->Draw("hist,same");
	projX->Draw("pe,same");
	projY->Draw("pe,same");
	TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
	legend->SetTextSize(0.02);
	legend->AddEntry(projX,"observed");
	legend->AddEntry(projY,"true");
	legend->Draw("same");
	return canvas;
}

