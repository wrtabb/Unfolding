#include "include/ToyModel.hh"
#include "include/Unfolding.hh"

TCanvas*PlotProjections(TString canvasName,TH2F*hMatrix,TH1F*hTrue,TH1F*hReco,bool logPlot);
vector<double> CreateConstantWidthBinning(float lowEdge,float highEdge,float binWidth);

void testToyModelCreation()
{
	TH1::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);	

	float lowEdge = 15;
	float highEdge = 3000;
	float trueBinWidth = 5;
	float recoBinWidth = 2.5;
	vector<double> binTrue = CreateConstantWidthBinning(lowEdge,highEdge,trueBinWidth);
	vector<double> binReco = CreateConstantWidthBinning(lowEdge,highEdge,recoBinWidth);
	
	binTrue = peakmodels::_peakbinningTrueDY;
	binReco = peakmodels::_peakbinningRecoDY;
	//binTrue = binningmodels::_massbinningTrue; 
	//binReco = binningmodels::_massbinningReco0;
	int nBinsTrue = binTrue.size()-1;
	int nBinsReco = binReco.size()-1;

	cout << endl;
	cout << "******************************" << endl;
	cout << "True binning: " << endl;
	for(int i=0;i<=nBinsTrue;i++){
		cout << binTrue.at(i) << endl;
	}
	cout << endl;
	cout << "******************************" << endl;
	cout << "Reco binning: " << endl;
	for(int i=0;i<=nBinsReco;i++){
		cout << binReco.at(i) << endl;
	}

	double norm = 1e6;//overall scaling factor
        double peakNormRel = 0.00002;//scaling factor for the peak only
        double mean = 91;//mean of the peak
        double sigma = 3.8;//spread of the peak
        double resSigma = 2;
        double shift = 15;//shift of the asymptote of the power function

	cout << "Asymptote shift = " << shift << endl;
	ToyModel*model = new ToyModel(norm,shift,peakNormRel,mean,sigma,resSigma,binTrue,
                                      binReco);
	Long64_t nEntries = 1e7;
	TH1F*hTrue = model->GetTrueHistRandom("hTrue",nEntries);
        TH1F*hReco = model->GetRecoHistRandom("hReco",nEntries);
	TH2F*hMatrix = model->GetMigrationMatrix("hMatrix");
        TH2F*hResponse = makeResponseMatrix(hMatrix);
	TH2F*hResponseRebin = Utilities::RebinTH2(hResponse,"hResponseRebin",hTrue);

	TH1F*hTrueClosed = model->GetTrueHist("hTrueClosed");
        TH1F*hRecoClosed = model->GetRecoHist("hRecoClosed");
	
	bool logPlot = false;
	TCanvas*c1 = PlotProjections("c1",hMatrix,hTrueClosed,hRecoClosed,logPlot);
	TCanvas*c2 = Utilities::PlotMatrix("c2","response matrix",hResponseRebin,logPlot);

	Unfold*unfClosure = new Unfold();
	Unfold::RegType regType = unfClosure->NO_REG;
	TH1F*hUnfoldedClosed = unfClosure->unfoldTUnfold(regType,hRecoClosed,hTrueClosed,
							  hMatrix);
	TCanvas*c3 = unfClosure->plotUnfolded("c3","unfolded closure test",hRecoClosed,hTrueClosed,hUnfoldedClosed,logPlot);
	Unfold*unfold = new Unfold();
	TH1F*hUnfolded = unfold->unfoldTUnfold(regType,hReco,hTrue,hMatrix);
	TCanvas*c4 = unfold->plotUnfolded("c4","unfolded",hReco,hTrue,hUnfolded,logPlot);

	TString baseSave = "plots/test";
	TString tagSave = "_PeakRegion";
	TString projSave = "Projections";
	TString respSave = "ResponseMatrix";
	TString unfClosureSave = "Unfold_Closure";
	TString unfoldSave = "Unfold";
	
	TString saveName = baseSave+projSave+tagSave+".png";
	c1->SaveAs(saveName);
	saveName = baseSave+respSave+tagSave+".png";
	c2->SaveAs("plots/testResponseMatrix.png");
	saveName = baseSave+unfClosureSave+tagSave+".png";
	c3->SaveAs("plots/testUnfold_Closure.png");
	saveName = baseSave+unfoldSave+tagSave+".png";
	c4->SaveAs("plots/testUnfold.png");
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
        projY->SetMarkerStyle(25);
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
	cout << endl;
	cout << "Plotting projections" << endl;
	cout << "Max histogram value: " << maxY << endl;
        float maxYRange = maxY*1.2;
	cout << "Setting top of canvas to " << maxYRange;
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

	int recoOver = nBinsReco;
	int trueOver = nBinsTrue;
	cout << "*********************************************" << endl;
	cout << "Now check underflow and overflow comparisons:" << endl;
	cout << "*********************************************" << endl;
	cout << "hTrue underflow = " << hTrueClone->GetBinContent(0) << endl;
	cout << "projY underflow = " << projY->GetBinContent(0) << endl;
	cout << endl;
	cout << "hReco underflow = " << hReco->GetBinContent(0) << endl;
	cout << "projX underflow = " << projX->GetBinContent(0) << endl;
	cout << endl;
	cout << "hTrue overflow = " << hTrueClone->GetBinContent(trueOver) << endl;
	cout << "projY overflow = " << projY->GetBinContent(trueOver) << endl;
	cout << endl;
	cout << "hReco overflow = " << hReco->GetBinContent(recoOver) << endl;
	cout << "projX overflow = " << projX->GetBinContent(recoOver) << endl;
	cout << endl;
        return canvas;

}

vector<double> CreateConstantWidthBinning(float lowEdge,float highEdge,float binWidth)
{
	vector<double> binning;
	float binEdge = lowEdge;
	for(int i=0;i<1e6;i++){
		if(binEdge>highEdge){
			binEdge = lowEdge;
			break;
		}
		binning.push_back(binEdge);
		binEdge += binWidth;
	}
	return binning;
}
