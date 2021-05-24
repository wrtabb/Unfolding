#include "include/ToyModel.hh"
#include "include/Unfolding.hh"

using namespace Utilities;
using namespace GlobalVariables;

void inversionTUnfoldComparison()
{
	TH1::SetDefaultSumw2();
        gStyle->SetOptStat(0);
        gStyle->SetPalette(1);
//        gROOT->SetBatch(true);

	vector<double> binTrue = _massbinningTrue;
        vector<double> binReco = _massbinningReco;
	
	int sizeTrue = binTrue.size();
	int sizeReco = binReco.size();
	int nBinsTrue = sizeTrue - 1;
	int nBinsReco = sizeReco - 1;

	double truebinning[sizeTrue];
	double recobinning[sizeReco];

	for(int i=0;i<sizeTrue;i++){
		truebinning[i] = binTrue.at(i);
	}
	for(int i=0;i<sizeReco;i++){
		recobinning[i] = binReco.at(i);
	}

	TH1F*hTrue = new TH1F("hTrue","",nBinsTrue,truebinning);
	TH1F*hReco = new TH1F("hReco","",nBinsReco,recobinning);
	TH1F*hRecoInv = RebinTH1(hReco,"hRecoInv",hTrue);

	double norm = 1000000;//overall scaling factor
        double peakNormRel = 0.00001;//scaling factor for the peak only
        double mean = 91;//mean of the peak
        double sigma = 3.8;//spread of the peak
        double shift = 60;//shift of the asymptote of the power function
        double resSigma = 2.0;//width of the smearing function
	ToyModel*model = new ToyModel(norm,shift,peakNormRel,mean,sigma,resSigma,binTrue,
                                      binReco);
	ToyModel*modelInv = new ToyModel(norm,shift,peakNormRel,mean,sigma,resSigma,binTrue,
                                         binTrue);
	
	TH2F*hMatrix = model->GetMigrationMatrix("hMatrix");
	TH2F*hMatrixInv = modelInv->GetMigrationMatrix("hMatrixInv");
	TF1*trueDist = model->GetTrueFunction();
	TF1*recoDist = model->GetRecoFunction();

	Long64_t nEntries = 1e7;
	double trueContent;
	double recoContent;
	for(Long64_t i=0;i<nEntries;i++){
		trueContent = trueDist->GetRandom();
		recoContent = recoDist->GetRandom();
		hTrue->Fill(trueContent);
		hReco->Fill(recoContent);
		hRecoInv->Fill(recoContent);		
	}

	Unfold*unfold = new Unfold();
	Unfold*unfInv = new Unfold();
	Unfold::RegType regType = unfold->VAR_REG_LCURVE;
	TH1F*hUnfolded = unfold->unfoldTUnfold(regType,hReco,hTrue,hMatrix);
	TH1F*hUnfoldedInv = unfInv->unfoldInversion(hReco,hTrue,hMatrix);

	TCanvas*c1 = unfold->plotUnfolded("c1","unfolded results",hReco,hTrue,hUnfolded,true);
	TCanvas*c2 = unfold->plotUnfolded("c2","unfolded results inversion",hRecoInv,hTrue,hUnfoldedInv,true);
}
