#include "include/Unfolding.hh"
#include "include/Utilities.h"

using namespace Utilities;

void testUnf2()
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

    Unfold*unf = new Unfold();
    unf->SetMatrix(hMatrix);
//    cout << "trueVert = " << unf->ReturnTrueVert() << endl;
//    TH2F*hResponse = unf->ReturnSquareResponseMatrix();

cout << "That's it!" << endl;
}
