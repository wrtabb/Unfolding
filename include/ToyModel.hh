#ifndef ToyModel_HH
#define ToyModel_HH

#include <TMath.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <Math/PdfFuncMathCore.h>

#include "FunctionalModels.h"
#include "GlobalVariables.h"
#include "Utilities.h"
#include "BinningModels.h"

class ToyModel
{
	public:
		ToyModel();
		ToyModel(double distNorm,double power,double peakNormRel,
			 double distMean,double distSigma,double resSigma,
			 vector<double>binningTrue,vector<double>binningReco);
		void SetModelParameters(double distNorm,double power,double peakNormRel,
					double distMean,double distSigma,double resSigma,
					vector<double>binningTrue,vector<double>binningReco);
		TH1F*GetTrueHist(TString trueName);
		TH1F*GetTrueHistRandom(TString trueName,Long64_t nEntries);
		TH1F*GetRecoHist(TString recoName);
		TH1F*GetRecoHistRandom(TString recoName,Long64_t nEntries);
		TH2F*GetMigrationMatrix(TString matrixName);
		TH2F*GetResponseMatrix(TH2F*hist);
		TH2F*GetResponseMatrixT(TH2F*hist);
		TF1*GetTrueFunction();
		TF1*GetRecoFunction();
		int GetNSigma();
	private:
		int _nSigma = 7; //number of resSigma away to integrate over
		double _distNorm; //Normalization factor
		double _power; //power of the power law
		double _shift; //how far to shift asymptote of power function
		double _peakNormRel; //normalization for peak relative to power function
		double _distMean; //mean of the peak
		double _distSigma; //width of the peak
		double _resSigma; //width of the smearing function
		double _xmin;
		double _xmax;
		double _xMin;
		double _xMax;
		int _nBinsTrue;
		int _nBinsReco;
		vector<double>_binningTrue;	
		vector<double>_binningReco;	
		TF1*_trueFunc;
		TF1*_recoFunc;
		TF2*_matrixFunc;
		void SetModelFunctions();
		TF1*Get1DIntegrand();
		TF2*Get2DIntegrand();
};

double trueDistribution(double*var,double*par);
double resolutionFunction(double*var,double*par);
double recoDistribution(double*var,double*par);
double integrand1D(double*var,double*par);
double integrand2D(double*var,double*par);










#endif
