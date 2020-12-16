#include "include/Unfolding.hh"

void makeToyModels()
{
 gStyle->SetOptStat(0);
 gROOT->SetBatch(true);
 double mean = 25.0;
 double xlow = 0.0;
 double xhigh = 50.0;
 int nBinsReco = 50;
 Unfold*models = new Unfold();
 models->makeToyModels(mean,xlow,xhigh,nBinsReco);
}
