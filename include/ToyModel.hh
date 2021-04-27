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

class ToyModel
{
	public:
		ToyModel(double distNorm,double power,double peakNormRel,
			 double distMean,double distSigma,double resSigma,
			 vector<double>binningTrue,vector<double>binningReco);
		TH1F*GetTrueHist(TString trueName);
		TH1F*GetRecoHist(TString recoName);
		TH2F*GetMigrationMatrix(TString matrixName);
		TH2F*GetResponseMatrix(TH2F*hist);
		TH2F*GetResponseMatrixT(TH2F*hist);
	private:
		double _distNorm;
		double _power;
		double _peakNormRel;
		double _distMean;
		double _distSigma;
		double _resSigma;
		double _xmin;
		double _xmax;
		int _nBinsTrue;
		int _nBinsReco;
		vector<double>_binningTrue;	
		vector<double>_binningReco;	
		TF1*_trueFunc;
		TF1*_recoFunc;
		TF2*_matrixFunc;
		void SetModelParameters(double distNorm,double power,double peakNormRel,
					double distMean,double distSigma,double resSigma,
					vector<double>binningTrue,vector<double>binningReco);
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
