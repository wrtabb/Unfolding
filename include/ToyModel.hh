#ifndef ToyModel_HH
#define ToyModel_HH
#include <TMath.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <Math/PdfFuncMathCore.h>
#include "FunctionalModels.h"

class ToyModel
{
	public:
		ToyModel(double distNorm,double power,double peakNormRel,
			 double distMean,double distSigma,double resSigma,
			 double xmin,double xmax,int nBins);
		TH1F*GetTrueHist();
		TH1F*GetRecoHist(bool useTUnfold = false);
		TH2F*GetMigrationMatrix(bool useTUnfold = false);
		TH2F*GetResponseMatrix(TH2F*hist);
	private:
		double _distNorm;
		double _power;
		double _peakNormRel;
		double _distMean;
		double _distSigma;
		double _resSigma;
		double _xmin;
		double _xmax;
		int _nBins;
		
		TF1*_trueFunc;
		TF1*_recoFunc;
		TF2*_matrixFunc;
		void SetModelParameters(double distNorm,double power,double peakNormRel,
					double distMean,double distSigma,double resSigma,
					double xmin,double xmax,int nBins);
		void SetModelFunctions(TF1*trueFunc,TF1*recoFunc,TF2*matrixFunc);
		TF1*Get1DIntegrand();
		TF2*Get2DIntegrand();
};

double trueDistribution(double*var,double*par);
double resolutionFunction(double*var,double*par);
double recoDistribution(double*var,double*par);
double integrand1D(double*var,double*par);
double integrand2D(double*var,double*par);










#endif
