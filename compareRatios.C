const TString fileNamePeak   = "data/peakComparison.root";
const TString fileNameFull = "data/toyModelRecoBin5.root";
TCanvas*plotWithRatio(TH1F*hist1,TH1F*hist2);

void compareRatios()
{
	TFile*filePeak = new TFile(fileNamePeak);
	TFile*fileFull = new TFile(fileNameFull);

	TH1F*hUnfoldedPeak = (TH1F*)filePeak->Get("hUnfoldedDY");
	hUnfoldedPeak->SetMarkerStyle(20);
	hUnfoldedPeak->SetMarkerColor(kBlack);
	hUnfoldedPeak->SetLineColor(kBlack);
	TH1F*hUnfoldedFull = (TH1F*)fileFull->Get("hUnfoldedTUnfold");
	hUnfoldedFull->SetFillColor(kRed+2);
	hUnfoldedFull->SetLineColor(kRed+2);
	hUnfoldedFull->SetMarkerColor(kRed+2);
	hUnfoldedFull->GetXaxis()->SetRangeUser(60,120);
	hUnfoldedFull->GetXaxis()->SetNoExponent();
	hUnfoldedFull->GetXaxis()->SetMoreLogLabels();
	hUnfoldedFull->SetTitle("unfolded peak vs. unfolded full");
	
	TCanvas*c1=new TCanvas("c1","",0,0,1200,1000);
	c1->SetGrid();
	TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
        legend->SetTextSize(0.02);
	legend->AddEntry(hUnfoldedFull,"unfolded full");
	legend->AddEntry(hUnfoldedPeak,"unfolded peak");

	hUnfoldedFull->Draw("hist");
	hUnfoldedPeak->Draw("pe,same");
	legend->Draw("same");

	c1->SaveAs("plots/comparePeaks.png");
}

