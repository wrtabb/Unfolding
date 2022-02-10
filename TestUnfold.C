#include "include/Unfolding.hh"
#include "include/Utilities.h"

void TestUnfold()
{
    gROOT->SetBatch(true);
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);

    TString savePrefix = "plots/testUnfold";
    TString saveSuffix = ".png";

    TString data_name = "data.root";
    TString dyll_name = "dyll.root";
    TFile*dataFile = new TFile(data_name);
    TFile*dyllFile = new TFile(dyll_name);

    TH1F*hReco   = (TH1F*)dyllFile->Get("histInvMassReco"); 
    TH1F*hTrue   = (TH1F*)dyllFile->Get("histInvMassDressed");
    TH2F*hMatrix = (TH2F*)dyllFile->Get("histMatrixInvMassDressed");

    Unfold*unf = new Unfold(hReco,hTrue,hMatrix);

    Unfold::UnfoldType unfTUnfold = Unfold::TUNFOLD;
    unf->EngageUnfolding(unfTUnfold);
    TCanvas*c2 = unf->plotUnfolded("c2","TUnfold Test",true);
    TString tunfoldSave = savePrefix;
    tunfoldSave += "_TUnfold_Unfolded";
    tunfoldSave += saveSuffix;
    c2->SaveAs(tunfoldSave);

    TH2F*hResponse = unf->ReturnSquareResponseMatrix();
    TString responseSave = savePrefix;
    responseSave += "_ResponseMatrix";
    responseSave += saveSuffix;
    unf->plotMatrix(hResponse,responseSave,true);
}
