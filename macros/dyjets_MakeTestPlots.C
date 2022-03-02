
void HistOptions(TH1D*hist,int kColor);
void PlotProjections(TH1D*hProjX,TH1D*hProjY,TH1D*hGen,TH1D*hRec,TString saveName);
TCanvas*MakeCanvas(TString cName);
TString file_name = "dyjets-DYJets.root";

void dyjets_MakeTestPlots()
{
    TFile*load_file = new TFile(file_name);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);

    // Get gen and reco plots
    TH1D*hGen = (TH1D*)load_file->Get("gen_mass_wide_range_exc0jet");
    HistOptions(hGen,kRed);
    TH1D*hRec = (TH1D*)load_file->Get("reco_mass_wide_range_exc0jet");
    HistOptions(hRec,kBlue);

    // Get migration matrix plot
    TH2D*hMatrix = (TH2D*)load_file->Get("reco_mass_wide_range_exc0jet-matrix");
    hMatrix->RebinY(2);
    hMatrix->GetXaxis()->SetTitle("m_{reco} [GeV]");
    hMatrix->GetYaxis()->SetTitle("m_{gen} [GeV]");
    hMatrix->SetTitle("muon migration matrix");
    hMatrix->GetXaxis()->SetNoExponent();
    hMatrix->GetXaxis()->SetMoreLogLabels();
    hMatrix->GetYaxis()->SetNoExponent();
    hMatrix->GetYaxis()->SetMoreLogLabels();

    // Draw migration matrix
    TCanvas*c2 = MakeCanvas("c2");
    hMatrix->Draw("colz");
    c2->SaveAs("testMassMatrix.png");

    // Make x and y projections of migration matrix
    TH1D*hProjX = (TH1D*)hMatrix->ProjectionX();
    hProjX->SetLineColor(kBlue);
    TH1D*hProjY = (TH1D*)hMatrix->ProjectionY();
    hProjY->SetLineColor(kRed);

    TString saveProjectionsName = "testMassPlot_MatrixProjections.png";
    PlotProjections(hProjX,hProjY,hGen,hRec,saveProjectionsName);
}

void PlotProjections(TH1D*hProjX,TH1D*hProjY,TH1D*hGen,TH1D*hRec,TString saveName)
{
    TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
    legend->SetTextSize(0.02);
    legend->AddEntry(hGen,"gen");
    legend->AddEntry(hRec,"reco");
    legend->AddEntry(hProjX,"x-projection");
    legend->AddEntry(hProjY,"y-projection");

    TLine*line = new TLine(15,1,3000,1);
    line->SetLineColor(kBlack);

    double ratioRange = 0.5;
    double upperBound = 1.0-ratioRange;
    double lowerBound = 1.0+ratioRange;
    TH1F*hRatioTrue = (TH1F*)hGen->Clone("trueRatio");
    hRatioTrue->Divide(hProjY);
    hRatioTrue->SetMarkerStyle(20);
    hRatioTrue->SetMarkerColor(kRed);
    hRatioTrue->SetMinimum(1.0-ratioRange);
    hRatioTrue->SetMaximum(1.0+ratioRange);

    TH1F*hRatioReco = (TH1F*)hRec->Clone("recoRatio");
    hRatioReco->Divide(hProjX);
    hRatioReco->SetMarkerStyle(20);
    hRatioReco->SetMarkerColor(kBlue);
    hRatioReco->SetMinimum(1.0-ratioRange);
    hRatioReco->SetMaximum(1.0+ratioRange);

    TCanvas*canvas = new TCanvas("canvas","",0,0,1000,1000);
    const float padmargins = 0.03;
    const float yAxisMinimum = 0.1;
    const float yAxisMaximum = 1e7;
    TPad*pad1 = new TPad("","",0,0.3,1.0,1.0);
    pad1->SetLogx();
    pad1->SetLogy();
    pad1->SetBottomMargin(padmargins);
    pad1->SetGrid();
    pad1->SetTicks(1,1);
    pad1->Draw();
    pad1->cd();

    hProjX->SetTitle("matrix projections vs 1D distributions");
    hProjX->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
    hProjX->GetXaxis()->SetNoExponent();
    hProjX->GetXaxis()->SetMoreLogLabels();
    hProjX->Draw("hist");
    hProjY->Draw("hist,same");
    hGen->Draw("pe,same");
    hRec->Draw("pe,same");
    legend->Draw("same");

    double ratioSplit = 0.18;
    canvas->cd();
    TPad*pad2 = new TPad("","",0,ratioSplit,1,0.3);
    pad2->SetLogx();
    pad2->SetTopMargin(padmargins);
    pad2->SetBottomMargin(0.2);
    pad2->SetGrid();
    pad2->SetTicks(1,1);
    pad2->Draw();
    pad2->cd();
    hRatioTrue->GetYaxis()->SetLabelSize(0.06);
    hRatioTrue->GetYaxis()->SetTitleSize(0.08);
    hRatioTrue->GetYaxis()->SetTitleOffset(0.3);
    hRatioTrue->GetYaxis()->SetTitle("gen/matrix projection");
    hRatioTrue->GetXaxis()->SetLabelSize(0);
    hRatioTrue->GetXaxis()->SetTitleSize(0);
    hRatioTrue->Draw("pe");
    line->Draw("same");

    canvas->cd();
    TPad*pad3 = new TPad("","",0,0.05,1,ratioSplit);
    pad3->SetLogx();
    pad3->SetTopMargin(padmargins);
    pad3->SetBottomMargin(0.2);
    pad3->SetGrid();
    pad3->SetTicks(1,1);
    pad3->Draw();
    pad3->cd();
    hRatioReco->GetYaxis()->SetLabelSize(0.06);
    hRatioReco->GetYaxis()->SetTitleSize(0.08);
    hRatioReco->GetYaxis()->SetTitleOffset(0.3);
    hRatioReco->GetYaxis()->SetTitle("reco/matrix projection");
    hRatioReco->GetXaxis()->SetLabelSize(0.1);
    hRatioReco->GetXaxis()->SetTitleSize(0.1);
    hRatioReco->GetXaxis()->SetNoExponent();
    hRatioReco->GetXaxis()->SetMoreLogLabels();
    hRatioReco->GetXaxis()->SetTitle("mass [GeV]");
    hRatioReco->Draw("pe");
    line->Draw("same");

    canvas->SaveAs(saveName);
}
void HistOptions(TH1D*hist,int kColor)
{
    hist->SetTitle("dimuon invariant mass");

    hist->SetLineColor(kColor);
    hist->SetMarkerColor(kColor);
    hist->SetMarkerStyle(20);
}

TCanvas*MakeCanvas(TString cName)
{
    TCanvas*canvas = new TCanvas(cName,"",0,0,1000,1000);
    canvas->SetGrid();
    canvas->SetLogy();
    canvas->SetLogx();
    canvas->SetLogz();

    return canvas;
}
