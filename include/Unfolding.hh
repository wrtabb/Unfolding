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
#include <TString.h>

#include "GlobalVariables.h"
#include "Utilities.h"//custom utilities used in analysis

class Unfold
{
	public:
		//-------Enums-------//
		enum RegType {        //Strength of regularization
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

        /// Constructed for Unfold
		Unfold();
		Unfold(TH1F*hReco,TH1F*hTrue,TH2F*hMatrix);
        /** 
         * \\Performs the unfolding using TUnfold
         * RegType regType is the type of regularization to be used
         * TH1F*hReco is the reconstructed distribution and must have more bins than hTrue
         * TH1F*hTrue is the true or gen-level distributions
         * TH2F*hMatrix is the matrix of migrations for unfolding. 
         * This does not need to be normalized
         */ 
		void unfoldTUnfold(RegType regType);
        /**
         * \\ Produces plot of unfolded distribution with true and reco distributions
             * TString canvasName is the name given to the TCanvas
             * TString titleName is the title printed at the top of the plot
         * TH1F*hReco is the reconstructed distribution
         * TH1F*hTrue is the true distribution
         * TH1F*hUnfolded is the unfolded distribution
         * bool logPlot determines if the x-axis is defined on a log scale or not
         */
		TCanvas*plotUnfolded(TString canvasName,TString titleName,bool logPlot);
        /**
         * \\Calculates the condition number of the response matrix
         * TH2F*hResponse is the response matrix
         * This function obtains the condition number of the response matrix
         * This is important because it tells us about how diagonal the matrix is
         * which gives an idea of how well unfolding will work
         * This helps us figure out if we need to use regularization or not
         * Small enough condition numbers mean no regularization is needed
         */
		double GetConditionNumber(TH2F*hResponse);

        /**
         * \\Normalizes the migration matrix to create the response matrix
         * TH2F*hist is the migration matrix
         * bool trueVert is true if the true distribution is on the y-axis
         * bool trueVert is false if the true distribution is on the x-axis
         * Currently this only normalizes rows, but may add functionality to do both
         * For now place true distribution along y-axis and reco along x-axis
         */ 
		TH2F*makeResponseMatrix(TH2F*hist);

        /**
         * \\Make a matrix from a given 2D histogram
         * TH2F*hist is the histogram to be made into a matrix
         * This is used because for the inversion method of unfolding
         * TMatrixD objects are used for the calculations
         */
		TMatrixD makeMatrixFromHist(TH2F*hist);

        /**
         * \\Make a vector from a given 1D histogram
         * TH1F*hist is the histogram to be made into a vector
         * This is used because for the inversion method of unfolding
         * TVectorD objects are used for the calculations
         */
		TVectorD makeVectorFromHist(TH1F*hist);

        /**
         * \\Makes a 1D histogram from a vector
         * TVectorD vec is the vector to be made into a histogram
         * Takes a TVectorD object and returns a histogram
         * TH1F*hist is used to get the desired binning correct
         */
		TH1F*makeHistFromVector(TVectorD vec,TH1F*hist);

        /**
         * \\Carries out unfolding using the inversion method
         * TH1F*hReco is reconstructed histogram
         * TH1F*hTrue is the true histogram
         * TH2F*hMatrix is the migration matrix
         */
		TH1F*unfoldInversion(TH1F*hReco,TH1F*hTrue,TH2F*hResponse);

        /**
         * \\Makes plot of a 2D histogram
         * TH2F*hmatrix is the matrix to be plotted
         * TString saveName is the save name for the plot
         * bool printCondition tells the function whether to print the condition number
         * on the plot or not
         */
		void plotMatrix(TH2F*hMatrix,TString saveName,bool printCondition);

        void SetReco(TH1F*hist);
        void SetTrue(TH1F*hist);
        void SetBackground(TH1F*hist);
        void SetMatrix(TH2F*hist);
        double ReturnCondition();

    private:
		//-----Variables-----//
		double _mean;
		double _condition;
		int _nBinsReco;
		int _nBinsTrue;
        bool _trueVert;
        bool _backgroundSubtraction = false;

        //-----Histograms-----//
        TH1F*_hReco;
        TH1F*_hTrue;
        TH1F*_hBack;
        TH1F*_hBlank;
        TH1F*_hUnfolded;
        TH2F*_hMatrix;
        TH2F*_hResponse;
		
};//end class Unfold

#endif
