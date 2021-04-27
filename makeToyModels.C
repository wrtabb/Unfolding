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
	//gROOT->SetBatch(true);
	
	vector<double> binTrue;
	vector<double> binReco;

	double trueBinWidth = 5.0;
	double recoBinWidth = 2.5;
	double trueBinEdge = 0.0;
	cout << "*************************" << endl;
	cout << "* Getting true binning: *" << endl;
	cout << "*************************" << endl;
	for(int i=0;i<1000000;i++){
		if(trueBinEdge > 165) break;
		cout << trueBinEdge << endl;
		binTrue.push_back(trueBinEdge);
		trueBinEdge += trueBinWidth;
	}

	const int sizeTrue = binTrue.size();
	cout << "Array size = " << sizeTrue << endl;
	double binTrueArray[sizeTrue];

	double recoBinEdge = 0.0;
	cout << endl;
	cout << "*************************" << endl;
	cout << "* Getting reco binning: *" << endl;
	cout << "*************************" << endl;
	for(int i=0;i<1000000;i++){
		if(recoBinEdge > 165) break;
		cout << recoBinEdge << endl;
		binReco.push_back(recoBinEdge);
		recoBinEdge += recoBinWidth;
	}
	const int sizeReco = binReco.size();
	cout << "Array size = " << sizeReco << endl;
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
	double peakNormRel = 0.00003;//scaling factor for the peak only
	double mean = 91;//mean of the peak
	double sigma = 10;//spread of the peak
	double shift = 6*sigma;//shift of the asymptote of the power function
	double resSigma = 2.0;

	//Create model from parameters above
	ToyModel*model = new ToyModel(norm,shift,peakNormRel,mean,sigma,resSigma,binTrue,
				      binReco);
	
	//Define all needed histogams from model
	//At minimum, true,reco, and matrix are required for unfolding
	TH1F*hTrue = model->GetTrueHist("hTrue");
	TH1F*hReco = model->GetRecoHist("hReco");
	TH1F*hRecoRebin = RebinTH1(hReco,"hRecoRebin",hTrue);
	TH2F*hMatrix = model->GetMigrationMatrix("hMatrix");

	TH2F*hResponse = model->GetResponseMatrix(hMatrix);

	//By default, the reco binning is finer than true binning
	//This is because TUnfold has this requirement
	//Here we rebin them for plotting
	TH2F*hResponseRebin = RebinTH2(hResponse,"hResponseRebin",hTrue);
	TH2F*hMatrixRebin = RebinTH2(hMatrix,"hMatrixRebin",hTrue);

	TString saveNameTag = "testNewBinningModel.png";
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
	TH1F*hRecoRebin = (TH1F*)hReco->Clone("hRecoRebin");
	TH2F*hMatrixRebin = (TH2F*)hMatrix->Clone("hMatrixRebin");
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

