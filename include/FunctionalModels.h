#ifndef FunctionalModels_H
#define FunctionalModels_H
#include "TMath.h"
#include "TF1.h"


double trueDistribution(double*var,double*par)
{
	double norm = par[0];
	double expDecay = par[1];
	double peakNormRel = par[2];
	double peakMean = par[3];
	double peakWidth = par[4];

        double x = var[0];
        bool isNormalized = true;

        double value = norm*exp(-x/expDecay)+norm*peakNormRel*
		       ROOT::Math::breitwigner_pdf(x,peakWidth,peakMean); 

        return value;
}

double resolutionFunction(double*var,double*par)
{

        double mean = par[0];
        double sigma = par[1];

        double x = var[0];

        bool isNormalized = true;

        double value = TMath::Gaus(x,mean,sigma,isNormalized);
        return value;
}

double integrand1D(double*var,double*par)
{
	double x = par[6];

	double norm = par[0];
	double expDecay = par[1];
	double peakNormRel = par[2];
	double peakMean = par[3];
	double peakWidth = par[4];

	double y = var[0];
	double parTrue[5] = {norm, expDecay, peakNormRel, peakMean, peakWidth};
	double varTrue[1] = {y};

	const double resMean = 0.0;
	double resWidth = par[5];
	double parRes[2] = {resMean, resWidth};
	double varRes[1] = {x-y};

	double value = trueDistribution(varTrue,parTrue)*resolutionFunction(varRes,parRes);
	return value;
}

double recoDistribution(double*var,double*par)
{
	double norm = par[0];
	double expDecay = par[1];
	double peakNormRel = par[2];
	double peakMean = par[3];
	double peakWidth = par[4];
	double resWidth  = par[5];
	double x = var[0];

	TF1 *integrand1DFunc = new TF1("integrand1DFunc",integrand1D,-10,150,7);
	integrand1DFunc->SetParameters(norm,expDecay,peakNormRel,peakMean,peakWidth,resWidth, 
				       x);

	double value = integrand1DFunc->Integral(x-7*resWidth, x+7*resWidth);
	return value;
}

Double_t integrand2D(Double_t *var, Double_t *par)
{
	double norm = par[0];
	double expDecay = par[1];
	double peakNormRel = par[2];
	double peakMean = par[3];
	double peakWidth = par[4];

	double x = var[0];
	double y = var[1];

	double parTrue[5] = {norm,expDecay,peakNormRel,peakMean,peakWidth};
	double varTrue[1] = {y};

	const double resMean = 0.0;
	double resWidth = par[5];
	double parRes[2] = {resMean,resWidth};
	double varRes[1] = {x-y};

	double varue = trueDistribution(varTrue,parTrue)*resolutionFunction(varRes,parRes);
	return varue;
}

#endif
