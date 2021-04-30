#include "include/ToyModel.hh"
#include "include/Unfolding.hh"


TCanvas*PlotProjections(TString canvasName,TH2F*hMatrix,TH1F*hTrue,TH1F*hReco);
TCanvas*PlotMatrix(TString canvasName,TString plotTitle,TH2F*hist);

// for binning types, see include/GlobalVariables.h
void makeToyModels(int binningType)
{
	TH1::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	gROOT->SetBatch(true);

	TString saveName = "data/toyModelRecoBin";
	saveName += binningType;
	saveName += ".root";
	
	vector<double> binTrue;
	vector<double> binReco;

	binTrue = _massbinningTrue0;

	if(binningType==0) binReco = _massbinningReco0;
	else if(binningType==1) binReco = _massbinningReco1;
	else if(binningType==2) binReco = _massbinningReco2;
	else if(binningType==3) binReco = _massbinningReco3;
	else if(binningType==4) binReco = _massbinningReco4;
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

	int nBinsTrue = sizeTrue - 1;;
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

	//By default, the reco binning is finer than true binning
	//This is because TUnfold has this requirement
	//Here we rebin them for plotting
	TH2F*hResponseRebin = RebinTH2(hResponse,"hResponseRebin",hTrue);
	TH2F*hMatrixRebin = RebinTH2(hMatrix,"hMatrixRebin",hTrue);

	TString saveNameTag = "_binning";
	saveNameTag += binningType;
	saveNameTag += ".png";

	//Plot matrices and migration matrix projections alongside true and reco
	TCanvas*c1 = PlotMatrix("c1","response matrix",hResponseRebin);
	c1->SaveAs("plots/responseMatrix"+saveNameTag);
	TCanvas*c2 = PlotMatrix("c2","migration matrix",hMatrixRebin);
	c2->SaveAs("plots/migrationMatrix"+saveNameTag);
	TCanvas*c4 = PlotProjections("c4",hMatrixRebin,hTrue,hRecoRebin);
	c4->SaveAs("plots/projectionsVsDistributions"+saveNameTag);

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
		trueContent = hTrue->GetRandom();
		recoContent = hReco->GetRandom();
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

TCanvas*PlotMatrix(TString canvasName,TString plotTitle,TH2F*hist)
{
	TCanvas*c1 = new TCanvas(canvasName,"",0,0,1000,1000);
//	c1->SetLogx();
//	c1->SetLogy();
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
	TH1F*hRecoRebin = RebinTH1(hReco,"hRecoRebin",hTrue);
	TH2F*hMatrixRebin = RebinTH2(hMatrix,"hMatrixRebin",hTrue);
	TCanvas*canvas = new TCanvas(canvasName,"",0,0,1000,1000);
	canvas->SetGrid();
//	canvas->SetLogx();
//	canvas->SetLogy();
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
	TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
	legend->SetTextSize(0.02);
	legend->AddEntry(projX,"observed");
	legend->AddEntry(projY,"true");
	legend->Draw("same");
	return canvas;
}

