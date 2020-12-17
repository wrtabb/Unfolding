#ifndef Unfold_HH
#define Unfold_HH

#include <TH1.h>
#include <TH2.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TLatex.h>
#include <TF1.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <iostream>
#include <TUnfold.h>
#include <TUnfoldDensity.h>
#include <TSpline.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <TDecompSVD.h>
#include <TString.h>
#include <TStyle.h>

class Unfold
{
	public:
		//-----Variables-----//
		double mean;
		double condition;
		int nBinsReco;
		int nBinsTrue;
		
		//-------Enums-------//
		enum RegType {      //Strength of regularization
  		NO_REG,           //No regularization
  		CONST_REG,        //User defined regularization
  		VAR_REG_LCURVE,   //TUnfoldDensity determines choice of regularization strength
  		VAR_REG_SCANSURE, //TUnfoldDensity determines choice of regularization strength
  		VAR_REG_SCANTAU   //TUnfoldDensity determines choice of regularization strength
		};
		enum UnfoldType {
  		TUNFOLD,
  		INVERSION,
		DNN
		};
		//-----Functions-----//
		Unfold();
		void makeToyModels(double mean,double xlow,double xhigh,int nBinsReco);
		void plotUnfolded(TH1D*hReco,TH1D*hTrue,TH1F*hUnfoldedE,UnfoldType unfoldType,
				  bool closure);
		TH1F*unfoldTUnfold(RegType regType,TH1D*hReco,TH1D*hTrue,TH2D*hMatrix);
		TH1F*unfoldInversion(TH1D*hReco,TH1D*hTrue,TH2D*hResponse);
		void plotMatrix(TH2D*hMatrix,TString saveName,bool printCondition);
		double GetConditionNumber(TH2D*hResponse);
		TMatrixD makeMatrixFromHist(TH2D*hist);
		TVectorD makeVectorFromHist(TH1D*hist);
		TH2D*makeResponseMatrix(TH2D*hist);
		TH1F*makeHistFromVector(TVectorD vec,TH1D*hist);
		TH2*RebinTH2(TH2*hist,TString histName,std::vector<double> binning);
		TH2*RebinTH2(TH2*hist,TString histName,TH2*hBinning);
		TH1*RebinTH1(TH1*hist,TString histName,std::vector<double> binning);
		TH1*RebinTH1(TH1*hist,TString histName,TH1*hBinning);

};//end class Unfold

#endif
