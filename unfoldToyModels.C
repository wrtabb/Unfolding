#include "include/ToyModel.hh"
#include "include/Unfolding.hh"

const TString saveFile = "./data/toyModel.root";

void unfoldToyModels()
{
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);

	TFile*file = new TFile(saveFile);
	TH1F*hTrue = (TH1F*)file->Get("hTrue");
	TH1F*hReco = (TH1F*)file->Get("hReco");
	TH2F*hMatrix = (TH2F*)file->Get("hMatrix");

	Unfold*unfold = new Unfold();
	Unfold::RegType regType = unfold->NO_REG;

	TH1F*hUnfolded = unfold->unfoldTUnfold(regType,hReco,hTrue,hMatrix);
	//TH1F*hUnfolded = unfold->unfoldInversion(hReco,hTrue,hMatrix);
	hUnfolded->SetMarkerStyle(25);
	hUnfolded->SetMarkerColor(kBlue+2);
	hUnfolded->SetLineColor(kBlue+2);

	//Here, the RebinTH1 function is called and does not throw an error
	TH1F*hRecoRebin = RebinTH1(hReco,"hRecoRebin",hTrue);
	hRecoRebin->SetMarkerStyle(20);
	hRecoRebin->SetMarkerColor(kBlack);
	//Here, the RebinTH2 function is called and does not throw an error
	TH2F*hMatrixRebin = RebinTH2(hMatrix,"hMatrixRebin",hTrue);

	//TCanvas*c3 = unfold->plotUnfolded("c3","Unfolded results",hReco,hTrue,hUnfolded);
/*
	TCanvas*c1 = new TCanvas("c1","",0,0,1000,1000);
	c1->SetGrid();
	hTrue->Draw("hist");
	hRecoRebin->Draw("pe,same");
	hUnfolded->Draw("pe,same");

	TH1F*ratioPlot = (TH1F*)hUnfolded->Clone("ratio");
	ratioPlot->Divide(hTrue);
	ratioPlot->SetMarkerStyle(20);
	ratioPlot->SetMinimum(0.8);
	ratioPlot->SetMaximum(1.2);
	ratioPlot->GetYaxis()->SetTitle("unfolded/true");
	ratioPlot->SetTitle("");
	ratioPlot->GetXaxis()->SetTitle("mass [GeV]");
	double x1 = hTrue->GetBinLowEdge(1);
	double x2 = hTrue->GetBinLowEdge(hTrue->GetNbinsX()+1);
	double y1 = 1.0;
	double y2 = y1;
	TLine*line = new TLine(x1,y1,x2,y2);
	line->SetLineColor(kRed);

	TCanvas*c2 = new TCanvas("c2","",500,0,1000,1000);
	c2->SetGrid();
	ratioPlot->Draw("pe");	
	line->Draw("same");
*/
}
