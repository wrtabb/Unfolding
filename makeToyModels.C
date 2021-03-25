#include "include/ToyModel.hh"
#include "include/Unfolding.hh"

const TString saveName = "data/toyModel.root";
TCanvas*PlotProjections(TString canvasName,TH2F*hMatrix,TH1F*hTrue,TH1F*hReco);
TCanvas*PlotMatrix(TString canvasName,TString plotTitle,TH2F*hist);

void makeToyModels()
{
	TH1::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	gROOT->SetBatch(true);
	
	//number of true bins
	//For reco distribution for TUnfold, nBinsReco = 2*nBinsTrue
	int nBins = 50;
	int nBinsReco = 2*nBins;

	//Set model parameters
	double norm = 10000000;//overall scaling factor
	double peakNormRel = 0.00003;//scaling factor for the peak only
	double mean = 91;//mean of the peak
	double sigma = 15;//spread of the peak
	double resSigma = 1;//resolution of the smearing from true to reco
	double shift = 6*sigma;//shift of the asymptote of the power function
	double xmin = 0;//minimum x axis value
	double xmax = 180;//maximum x axis value

	TString saveNameTag = "_nBins_";
	saveNameTag += nBins;
	saveNameTag += "_ResSmear_";
	saveNameTag += resSigma;
	saveNameTag += ".png";

	//Create model from parameters above
	ToyModel*model = new ToyModel(norm,shift,peakNormRel,mean,sigma,resSigma,xmin,xmax,
				      nBins);
	
	//Define all needed histogams from model
	//At minimum, true,reco, and matrix are required for unfolding
	TH1F*hTrue = model->GetTrueHist();
	TH1F*hReco = model->GetRecoHist();
	TH2F*hMatrix = model->GetMigrationMatrix();

	TH2F*hResponse = model->GetResponseMatrix(hMatrix);

	//By default, the reco binning is finer than true binning
	//This is because TUnfold has this requirement
	//Here we rebin them for plotting
	TH2F*hResponseRebin = (TH2F*)hResponse->Clone("hResponseRebin");
	hResponseRebin->RebinX(2);
	TH2F*hMatrixRebin = (TH2F*)hMatrix->Clone("hMatrixRebin");
	hMatrixRebin->RebinX(2);

	//Plot matrices and migration matrix projections alongside true and reco
	TCanvas*c1 = PlotMatrix("c1","response matrix",hResponseRebin);
	c1->SaveAs("plots/responseMatrix"+saveNameTag);
	TCanvas*c2 = PlotMatrix("c2","migration matrix",hMatrixRebin);
	c2->SaveAs("plots/migrationMatrix"+saveNameTag);
	TCanvas*c4 = PlotProjections("c4",hMatrix,hTrue,hReco);
	c4->SaveAs("plots/projectionsVsDistributions"+saveNameTag);

	//Now produce randomly filled distributions from the model
	//To test unfolding on distributions which are not identical
	//to those used to make the matrix
	Long64_t nEvents = 1e7;
	double trueContent;
	double recoContent;
	TH1F*hTrueRandom = new TH1F("hTrueRandom","",nBins,xmin,xmax);
	TH1F*hRecoRandom = new TH1F("hRecoRandom","",nBinsReco,xmin,xmax);
	for(Long64_t i=0;i<nEvents;i++){
		trueContent = hTrue->GetRandom();
		recoContent = hReco->GetRandom();
		hTrueRandom->Fill(trueContent);
		hRecoRandom->Fill(recoContent);	
	}
	
	//Save the results to a root file
	TFile*saveFile = new TFile(saveName,"recreate");
	hTrue->Write();
	hReco->Write();
	hMatrix->Write();
	hResponse->Write();
	hTrueRandom->Write();
	hRecoRandom->Write();
	saveFile->Close();
}

TCanvas*PlotMatrix(TString canvasName,TString plotTitle,TH2F*hist)
{
	TCanvas*c1 = new TCanvas(canvasName,"",0,0,1000,1000);
	c1->SetGrid();
	c1->SetRightMargin(0.13);
	c1->SetLeftMargin(0.13);
	hist->SetTitle(plotTitle);
	hist->GetYaxis()->SetTitle("mass^{true} [GeV]");
	hist->GetXaxis()->SetTitle("mass^{obs} [GeV]");
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
	TH1F*hTrueClone = (TH1F*)hTrue->Clone();
	hTrueClone->SetLineColor(kRed);
	hRecoRebin->SetLineColor(kBlue);
	hTrueClone->SetFillColor(kWhite);
	hRecoRebin->SetFillColor(kWhite);

	float maxY = 0.0;
	float binContent;
	for(int i=1;i<=nBinsTrue;i++){
		binContent = hTrueClone->GetBinContent(i);
		if(binContent > maxY) maxY = binContent;
	}
	float maxYRange = maxY*1.2;
	hTrueClone->GetXaxis()->SetTitle("mass [GeV]");
	hTrueClone->GetYaxis()->SetRangeUser(0,maxYRange);
	hTrueClone->SetTitle("matrix projections with distributions");
	hTrueClone->Draw("hist");
	hRecoRebin->Draw("hist,same");
	projX->Draw("pe,same");
	projY->Draw("pe,same");
	return canvas;
}

