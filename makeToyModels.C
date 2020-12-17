#include "include/Unfolding.hh"

void makeToyModels()
{
 gStyle->SetOptStat(0);
 gROOT->SetBatch(true);
 double mean = 30.0;
 double xlow = 0.0;
 double xhigh = 50.0;
 int nBinsReco = 60;
 Unfold*models = new Unfold();
 models->makeToyModels(mean,xlow,xhigh,nBinsReco);
}
