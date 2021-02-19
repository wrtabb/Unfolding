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
		void plotUnfolded(TH1F*hReco,TH1F*hTrue,TH1F*hUnfoldedE,UnfoldType unfoldType,
				  bool closure);
		TH1F*unfoldTUnfold(RegType regType,TH1F*hReco,TH1F*hTrue,TH2F*hMatrix);
		TH1F*unfoldInversion(TH1F*hReco,TH1F*hTrue,TH2F*hResponse);
		void plotMatrix(TH2F*hMatrix,TString saveName,bool printCondition);
		double GetConditionNumber(TH2F*hResponse);
		TMatrixD makeMatrixFromHist(TH2F*hist);
		TVectorD makeVectorFromHist(TH1F*hist);
		TH2F*makeResponseMatrix(TH2F*hist);
		TH1F*makeHistFromVector(TVectorD vec,TH1F*hist);
		TH2F*RebinTH2(TH2F*hist,TString histName,std::vector<double> binning);
		TH2F*RebinTH2(TH2F*hist,TString histName,TH2F*hBinning);
		TH1F*RebinTH1(TH1F*hist,TString histName,std::vector<double> binning);
		TH1F*RebinTH1(TH1F*hist,TString histName,TH1F*hBinning);

};//end class Unfold

#endif
