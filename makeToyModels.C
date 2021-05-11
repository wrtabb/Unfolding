#include "include/ToyModel.hh"
#include "include/Unfolding.hh"

using namespace Utilities;
using namespace GlobalVariables;

TCanvas*PlotProjections(TString canvasName,TH2F*hMatrix,TH1F*hTrue,TH1F*hReco,bool logPlot);
TCanvas*PlotMatrix(TString canvasName,TString plotTitle,TH2F*hist,bool logPlot);

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

	if(binningType>4) binTrue = _massbinningTrue;
	else binTrue = _massbinningTrue0;

	if(binningType==0) binReco = _massbinningReco0;
	else if(binningType==1) binReco = _massbinningReco1;
	else if(binningType==2) binReco = _massbinningReco2;
	else if(binningType==3) binReco = _massbinningReco3;
	else if(binningType==4) binReco = _massbinningReco4;
	else if(binningType==5) binReco = _massbinningReco;
	else if(binningType==6) binReco = binTrue;
	else if(binningType==7) binReco = _massbinningReco5;
	else if(binningType==8) binReco = _massbinningReco6;
	else if(binningType==9) binReco = _massbinningReco7;
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
	double peakNormRel = 0.00006;//scaling factor for the peak only
	double mean = 91;//mean of the peak
	double sigma = 10;//spread of the peak
	double shift = 6*sigma;//shift of the asymptote of the power function
	double resSigma = 2.0;

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

	TF1*trueFunction = model->GetTrueFunction();
	TF1*recoFunction = model->GetRecoFunction();

	//By default, the reco binning is finer than true binning
	//This is because TUnfold has this requirement
	//Here we rebin them for plotting
	TH2F*hResponseRebin = RebinTH2(hResponse,"hResponseRebin",hTrue);
	TH2F*hMatrixRebin = RebinTH2(hMatrix,"hMatrixRebin",hTrue);

	TString saveNameTag = "_binning";
	saveNameTag += binningType;
	saveNameTag += ".png";

	bool logPlot;
	if(binningType>4) logPlot = true;	
	else logPlot = false;
	//Plot matrices and migration matrix projections alongside true and reco
	TCanvas*c1 = PlotMatrix("c1","response matrix",hResponseRebin,logPlot);
	c1->SaveAs("plots/responseMatrix"+saveNameTag);
	TCanvas*c2 = PlotMatrix("c2","migration matrix",hMatrixRebin,logPlot);
	c2->SaveAs("plots/migrationMatrix"+saveNameTag);

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

TCanvas*PlotMatrix(TString canvasName,TString plotTitle,TH2F*hist,bool logPlot)
{
	TCanvas*c1 = new TCanvas(canvasName,"",0,0,1000,1000);
	if(logPlot){
		c1->SetLogx();
		c1->SetLogy();
		hist->GetYaxis()->SetNoExponent();
		hist->GetXaxis()->SetNoExponent();
		hist->GetYaxis()->SetMoreLogLabels();
		hist->GetXaxis()->SetMoreLogLabels();
	}
	c1->SetGrid();
	c1->SetRightMargin(0.13);
	c1->SetLeftMargin(0.13);
	hist->SetTitle(plotTitle);
	hist->GetYaxis()->SetTitle("mass^{true} [GeV]");
	hist->GetXaxis()->SetTitle("mass^{obs} [GeV]");
	hist->Draw("colz");	
	return c1;
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

