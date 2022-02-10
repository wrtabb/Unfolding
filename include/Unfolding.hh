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
		enum UnfoldType{    // Type of unfolding
            TUNFOLD,        // TUnfold (the primary method for this study)
            INVERSION,      // Matrix inversion method
            DNN             // Deep neural network (not implemented for this study)
		};

		//-----Functions-----//

        /**
        * \\Blank constructor for unfold
        * All necessary parameters can be set later
        * Necessary parameters:
        * reconstructed or observed histogram
        * true histogram
        * migration matrix
        * RegType regType defined the regularization to be used
        */

		Unfold();
        /**
        * \\Constructor for unfold
        * hReco is observed or reconstructed histogram
        * hTrue is true histogram
        * hMatrix is matrix of migrations
        * regType is the regularization type to be used
        * */
		Unfold(TH1F*hReco,TH1F*hTrue,TH2F*hMatrix,RegType regType);

        /**
        * \\ Carries out the unfolding
        * unfoldType is the type of unfolding to be done
        * unfoldType must be: TUNFOLD, INVERSION, or DNN 
        * DNN not implemented
        */
        void EngageUnfolding(UnfoldType unfoldType);

        /**
         * \\ Produces plot of unfolded distribution with true and reco distributions
         * TString canvasName is the name given to the TCanvas
         * TString titleName is the title printed at the top of the plot
         * bool logPlot determines if the x-axis is defined on a log scale or not
         */
		TCanvas*plotUnfolded(TString canvasName,TString titleName,bool logPlot);

        /**
         * \\Calculates the condition number of the response matrix
         * This function obtains the condition number of the response matrix
         * This is important because it tells us about how diagonal the matrix is
         * This helps us figure out if we need to use regularization or not
         * Small enough condition numbers mean no regularization is needed
         */
		void SetConditionNumber();

        /**
         * \\Normalizes the migration matrix to create the response matrix
         * TH2F*hist is the migration matrix
         * The axis to be normalized along is determined automatically
         */ 
		void makeResponseMatrix(TH2F*hist);

        /**
         * \\Make a matrix from a given 2D histogram
         * TH2F*hist is the histogram to be made into a matrix
         */
		TMatrixD makeMatrixFromHist(TH2F*hist,bool underflow,bool overflow);

        /**
         * \\Make a vector from a given 1D histogram
         * TH1F*hist is the histogram to be made into a vector
         */
		TVectorD makeVectorFromHist(TH1F*hist);

        /**
         * \\Makes a 1D histogram from a vector
         * Takes a TVectorD object and returns a histogram
         * TVectorD vec is the vector to be made into a histogram
         * TH1F*hist is used to get the desired binning correct
         */
		TH1F*makeHistFromVector(TVectorD vec,TH1F*hist);

        /**
         * \\Makes plot of a 2D histogram
         * TH2F*hmatrix is the matrix to be plotted
         * TString saveName is the save name for the plot
         * bool printCondition tells the function whether to print the condition number
         * on the plot or not
         */
		void plotMatrix(TH2F*hMatrix,TString saveName,bool printCondition);

        /*
         * \\ Sets the reco or observed distribution for unfolding
         * hist is the reco distribution
         * This function is only needed if it wasn't set by the constructor
         */ 
        void SetReco(TH1F*hist);

        /*
         * \\ Sets the true distribution for unfolding
         * hist is the true distribution
         * This function is only needed if it wasn't set by the constructor
         */ 
        void SetTrue(TH1F*hist);

        /*
         * \\ Sets the type of regularization to be used
         * regTypes are: 
         * NO_REG (default)  //No regularization used
         * CONST_REG,        //User defined regularization
         * VAR_REG_LCURVE,   //TUnfoldDensity determines choice of regularization strength
         * VAR_REG_SCANSURE, //TUnfoldDensity determines choice of regularization strength
         * VAR_REG_SCANTAU   //TUnfoldDensity determines choice of regularization strength
         */
        void SetRegularizationType(RegType regType);

        /*
         * \\ Sets the background distribution
         * hist is a histogram of all backgrounds summed together
         */ 
        void SetBackground(TH1F*hist);

        /*
         * \\ Sets the migration matrix for unfolding
         * hist is the 2D histogram of the migration matrix
         * This function is only needed if it wasn't set by the constructor
         */ 
        void SetMatrix(TH2F*hist);

        /*
         * \\ Sets the attributes for the response matrix 
         * hist is the 2D histogram of the repsonse matrix
         */ 
        void SetMatrixPlotAttributes(TH2F*hist);

        /*
         * \\Returns the condition number of the response matrix
         */   
        double ReturnCondition();

        /*
         * \\Returns the determinant of the response matrix
         */   
        double ReturnDeterminant();

        /*
         * \\Returns the unfolded distribution 
         */   
        TH1F*ReturnUnfolded();

        /*
         * \\Returns the response matrix which is the normalized migration matrix
         */   
        TH2F*ReturnResponseMatrix();

        /*
         * \\Returns the response matrix rebinned to be a square matrix 
         */   
        TH2F*ReturnSquareResponseMatrix();

        /*
         * \\Returns the reco distribution 
         */   
        TH1F*ReturnReco();

        /*
         * \\Returns the true distribution 
         */   
        TH1F*ReturnTrue();

        /*
         * \\Returns the background distribution 
         */   
        TH1F*ReturnBackground();

        /*
         * \\Returns the migration matrix 
         */   
        TH2F*ReturnMatrix();

        bool ReturnTrueVert();
    private:
		//-----Variables-----//
		double _condition   = -1000;
        double _determinant = -1000;
		int _nBinsReco      = -1000;
		int _nBinsTrue      = -1000;

        bool _trueVert = false;
        bool _backgroundSubtraction = false;
        RegType _regType = NO_REG;

        //-----Histograms-----//
        TH1F*_hReco             = NULL;
        TH1F*_hTrue             = NULL;
        TH1F*_hBack             = NULL;
        TH1F*_hBlank            = NULL;
        TH1F*_hUnfolded         = NULL;
        TH2F*_hMatrix           = NULL;
        TH2F*_hResponse         = NULL;
        TH2F*_hResponseSquare   = NULL;

	    //-----Functions-----//	
        /** 
         * \\Performs the unfolding using TUnfold
         */ 
		void unfoldTUnfold();

        /**
         * \\Carries out unfolding using the inversion method
         */
		void unfoldInversion();
};//end class Unfold

#endif
