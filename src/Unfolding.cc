#include "../include/Unfolding.hh"
using namespace Utilities;

Unfold::Unfold()
{

}

Unfold::Unfold(TH1F*hReco,TH1F*hTrue,TH2F*hMatrix,RegType regType)
{
    // calclate global parameters
    _hReco = hReco; // reconstructed distribution
    _hTrue = hTrue; // true distribution
    _hMatrix = hMatrix;// matrix of migrations
    makeResponseMatrix(_hMatrix); // normalized migration matrix
    _hResponseSquare = RebinTH2(_hResponse,"hResponseSquare",_hTrue);
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << "response matrix x-axis bins: " << _hResponse->GetNbinsY() << endl;
    cout << "square matrix x-axis bins: " << _hResponseSquare->GetNbinsY() << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    _condition = GetConditionNumber(); // condition number
    _nBinsReco = _hReco->GetNbinsX(); // number of reco bins
    _nBinsTrue = _hTrue->GetNbinsX(); // number of true bins
    _regType = regType; // regularization type

    // Determine if true distribution is on the vertical or horizontal axis
    int nBinsX = _hMatrix->GetNbinsX();
    int nBinsY = _hMatrix->GetNbinsY();
    if(nBinsY == _nBinsTrue) _trueVert = true;
    else if(nBinsY == _nBinsReco) _trueVert = false;

    if(_nBinsReco==_nBinsTrue){
        cout << "For TUnfold, the observed histogram must have more bins than the true histogram" << endl;
        cout << "Input bins: " << _nBinsReco << endl;
        cout << "Output bins: " << _nBinsTrue << endl;
    }
}

void Unfold::unfoldTUnfold()
{
    cout << endl;
    cout << "*****************************************" << endl;
    cout << "Beginning unfolding process using TUnfold" << endl;
    cout << "*****************************************" << endl;
    cout << endl;

    if(_hReco == NULL || _hTrue == NULL || _hMatrix == NULL){
        cout << "The reco and true distributions and the matrix of migrations must be specified to carry out unfolding" << endl;
        cout << "Please check that these are properly defined and try again" << endl;
        return;
    }

    TH1F*hReco = _hReco;
    TH1F*hTrue = _hTrue;
    TH2F*hMatrix = _hMatrix;

    ////////////////////////////
    //  Regularization Modes  //
    ////////////////////////////
    //TUnfold::ERegMode regMode = TUnfold::kRegModeNone; //breaks 
    //TUnfold::ERegMode regMode = TUnfold::kRegModeSize;
    //TUnfold::ERegMode regMode = TUnfold::kRegModeDerivative;
    TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
    //TUnfold::ERegMode regMode = TUnfold::kRegModeMixed; //breaks

    ///////////////////////////
    //  Types of Constraint  //
    ///////////////////////////
    TUnfold::EConstraint constraintMode = TUnfold::kEConstraintNone;
    //TUnfold::EConstraint constraintMode = TUnfold::kEConstraintArea;

    /////////////////////
    //  Density Modes  //
    /////////////////////
    //TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeNone;
    TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;
    //TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeUser;
    //TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidthAndUser;
    /////////////////////////////////////
    //  Horizontal vs Vertical Output  //
    /////////////////////////////////////
    TUnfold::EHistMap outputMap;
    if(_trueVert) outputMap = TUnfold::kHistMapOutputVert;
    else outputMap = TUnfold::kHistMapOutputHoriz;


    //////////////////////////////////////
    //  Constructor for TUnfoldDensity  //
    //////////////////////////////////////
    TUnfoldDensity unfold(hMatrix,outputMap,regMode,constraintMode,densityFlags);
    unfold.SetInput(hReco);//the measured distribution
    double backScale = 1.0;
    double backScaleError = 0.0;//scale error for background
    if(_backgroundSubtraction) 
        unfold.SubtractBackground(_hBack,"background",backScale,backScaleError);

    ///////////////////////
    //  Begin Unfolding  //
    ///////////////////////
    Int_t iBest;
    TSpline *logTauX,*logTauY;
    TGraph *lCurve;
    TGraph*bestLcurve;
    TGraph*bestLogTauLogChi2;
    if(_regType == VAR_REG_LCURVE){
        Int_t nScan=30;//This number chosen only because it was given in the tutorial
        Double_t tauMin = 0.0;//If tauMin=tauMax, TUnfold automatically chooses a range
        Double_t tauMax = 0.0;//Not certain how TUnfold chooses the range
        iBest=unfold.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
        cout<< "tau=" << unfold.GetTau() << endl;
        cout<<"chi**2="<<unfold.GetChi2A()<<"+"<<unfold.GetChi2L()<<" / "<<unfold.GetNdf()<<"\n";
        Double_t t[1],x[1],y[1];
        logTauX->GetKnot(iBest,t[0],x[0]);
        logTauY->GetKnot(iBest,t[0],y[0]);
        bestLcurve=new TGraph(1,x,y);
        bestLogTauLogChi2=new TGraph(1,t,x);
    }
    else if(_regType == VAR_REG_SCANSURE){
        cout << "Option VAR_REG_SCANSURE not yet implemented" << endl;
        return;
    }
    else if(_regType == VAR_REG_SCANTAU){
        cout << "Option VAR_REG_SCANTAU not yet implemented" << endl;
        return;
    }
    else if(_regType == NO_REG){
        double tau = 0;
        unfold.DoUnfold(tau,hReco);
    }
    else{//user defined
        //choose a tau
        //this isn't implemented more elegantly because I don't use it
        double tau = 1e-8;//larger tau introduces more MC bias
        unfold.DoUnfold(tau,hReco);
    }

    //Create unfolded histogram
    TH1*hUnfolded = unfold.GetOutput("hUnfolded");

    //Create error matrices
    TH2*histEmatStat=unfold.GetEmatrixInput("hEmatrixInput");
    TH2*histEmatTotal=unfold.GetEmatrixTotal("hEmatrixTotal");

    //Create unfolding histogram with errors
    TH1F*hUnfoldedE = (TH1F*)hUnfolded->Clone("hUnfoldedTUnfold");

    //loop over unfolded histogram bins and assign errors to each one
    for(int i=0;i<=_nBinsTrue;i++){
        double binError = TMath::Sqrt(histEmatTotal->GetBinContent(i+1,i+1));
        hUnfoldedE->SetBinError(i+1,binError);
    }

    _hUnfolded = hUnfoldedE;
}//end unfoldTUnfold

TCanvas*Unfold::plotUnfolded(TString canvasName,TString titleName,bool logPlot)
{
    gStyle->SetOptStat(0);

    //Get parameters for histograms
    //These are used to define the canvases and pads
    //as well as the locations of drawn objects
    int binLow = _hTrue->GetBinLowEdge(1);
    int binHigh = _hTrue->GetBinLowEdge(_nBinsTrue)+_hTrue->GetBinWidth(_nBinsTrue);
    float peakMax = 0;	
    float binContent;
    for(int i=1;i<=_nBinsTrue;i++){
        binContent = _hTrue->GetBinContent(i);
        if(binContent > peakMax) peakMax = binContent;
    }

    //Rebin reco to plot alongside true for easy comparison
    TH1F*hRecoRebin2 = RebinTH1(_hReco,"hRecoRebin2",_hTrue);

    //set histogram drawing options
    _hTrue->SetFillColor(kRed+2);
    _hTrue->SetLineColor(kRed+2);
    hRecoRebin2->SetMarkerStyle(20);
    hRecoRebin2->SetMarkerColor(kBlack);
    hRecoRebin2->SetLineColor(kBlack);
    _hUnfolded->SetMarkerStyle(25);
    _hUnfolded->SetMarkerColor(kBlue+2);
    _hUnfolded->SetLineColor(kBlue+2);
    _hUnfolded->SetFillColor(kWhite);

    //define the ratio plot
    TH1F*ratio = (TH1F*)_hUnfolded->Clone("ratio");
    ratio->Divide(_hTrue);
    ratio->SetTitle("");

    //define the legend
    TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
    legend->SetTextSize(0.02);
    legend->AddEntry(_hTrue,"True Distribution");
    legend->AddEntry(hRecoRebin2,"Observed Distribution");
    legend->AddEntry(_hUnfolded,"Unfolded Distribution");

    //Create a label that shows the chi^2 value to print on graph
    float xMax = _hTrue->GetXaxis()->GetXmax();
    float xMin = _hTrue->GetXaxis()->GetXmin();
    float xRange = xMax-xMin;
    float yMax = 1.1*peakMax;
    double xChiLabel;
    double yChiLabel;

    if(logPlot){
        xChiLabel = xRange*0.2+xMin;
        yChiLabel = yMax*0.15;
    }
    else {
        xChiLabel = xRange*0.66+xMin;
        yChiLabel = yMax*0.7;
    }
    double x[_nBinsTrue],res[_nBinsTrue];
    double chi = _hUnfolded->Chi2Test(_hTrue,"CHI2/NDF",res);//chi2/ndf to print on plot
    double pValues = _hUnfolded->Chi2Test(_hTrue,"P",res);//outputs chi2,prob,ndf,igood
    TLatex*chiLabel = new TLatex(xChiLabel,yChiLabel,Form("#chi^{2}/ndf = %lg", chi));

    //Draw canvas and pads to make plot
    TCanvas*canvas = new TCanvas(canvasName,"",0,0,1000,1000);
    const float padmargins = 0.03;
    const float yAxisMinimum = 100;
    const float yAxisMaximum = 1e7;
    TPad*pad1 = new TPad("","",0,0.3,1.0,1.0);
    if(logPlot){
        pad1->SetLogx();
        pad1->SetLogy();
    }

    pad1->SetBottomMargin(padmargins);
    pad1->SetGrid();
    pad1->SetTicks(1,1);
    pad1->Draw();
    pad1->cd();
    _hTrue->SetLabelSize(0);
    _hTrue->SetTitleSize(0);
    _hTrue->SetMinimum(yAxisMinimum);
    _hTrue->SetMaximum(yAxisMaximum);
    _hTrue->SetTitle(titleName);
    _hTrue->Draw("hist");
    hRecoRebin2->Draw("pe,same");
    _hUnfolded->Draw("pe,same");	
    legend->Draw("same");
    chiLabel->Draw("same");

    canvas->cd();
    TPad*pad2 = new TPad("","",0,0.05,1,0.3);
    if(logPlot) pad2->SetLogx();
    pad2->SetTopMargin(padmargins);
    pad2->SetBottomMargin(0.2);
    pad2->SetGrid();
    pad2->SetTicks(1,1);
    pad2->Draw();
    pad2->cd();
    ratio->SetMinimum(0.7);
    ratio->SetMaximum(1.3);
    ratio->GetYaxis()->SetLabelSize(0.06);
    ratio->GetYaxis()->SetTitleSize(0.08);
    ratio->GetYaxis()->SetTitleOffset(0.3);
    ratio->GetYaxis()->SetTitle("Unfolded/Truth");
    ratio->GetXaxis()->SetLabelSize(0.1);
    ratio->GetXaxis()->SetTitleSize(0.1);
    ratio->GetXaxis()->SetNoExponent();
    ratio->GetXaxis()->SetMoreLogLabels();
    ratio->GetXaxis()->SetTitle("mass [GeV]");
    ratio->SetMarkerStyle(20);
    ratio->SetMarkerColor(kBlack);
    ratio->SetLineColor(kBlack);
    ratio->Draw("PE");
    TLine*line = new TLine(binLow,1,binHigh,1);
    line->SetLineColor(kRed);
    line->Draw("same");

    return canvas;
}//end plotUnfolded

double Unfold::GetConditionNumber()
{
	TMatrixD matrix = makeMatrixFromHist(_hResponseSquare);

	TDecompSVD decomp(matrix);
	_condition = decomp.Condition();
	cout << "The condition number: " << _condition << endl;

	double determinant;
	TMatrixD mInverse = matrix.Invert(&determinant);
	cout << "The determinant: " << determinant << endl;
	return _condition;
}//end GetConditionNumber

void Unfold::makeResponseMatrix(TH2F*hist)
{
    bool trueVert = _trueVert;
    TH2F*hResponse = (TH2F*)hist->Clone("hResponse");
    int nBinsX = hResponse->GetNbinsX();
    int nBinsY = hResponse->GetNbinsY();
    double nEntriesX,nEntriesY;
    double binContent;
    double scaledContent;

    // true distribution on y-axis
    // reco distribution on x-axis
    if(trueVert){
        //Loop over all true bins (the y-axis)
        for(int j=0;j<=nBinsY+1;j++){
            nEntriesX = 0.0;
            //for each true bin, sum up the number of events across all reco bins
            for(int i=0;i<=nBinsX+1;i++){
                binContent = hist->GetBinContent(i,j);
                nEntriesX += binContent;
            }//end first loop over reco bins
            //For each true bin scale the bin content by the number of entries in all
            //reco bins, then place this content into the new matrix
            for(int i=0;i<=nBinsX+1;i++){
                scaledContent = hist->GetBinContent(i,j)/nEntriesX;
                hResponse->SetBinContent(i,j,scaledContent);
            }//end second loop over reco bins
        }//end loop over true bins
    }// end if trueVert

    // true distribution on x-axis
    // reco distribution on y-axis
    else{
        //Loop over all true bins (the x-axis)
        for(int i=0;i<=nBinsX+1;i++){
            nEntriesY = 0.0;
            //for each true bin, sum up the number of events across all reco bins
            for(int j=0;j<=nBinsY+1;j++){
                binContent = hist->GetBinContent(i,j);
                nEntriesY += binContent;
            }//end first loop over reco bins
            //For each true bin scale the bin content by the number of entries in all
            //reco bins, then place this content into the new matrix
            for(int j=0;j<=nBinsY+1;j++){
                scaledContent = hist->GetBinContent(i,j)/nEntriesY;
                hResponse->SetBinContent(i,j,scaledContent);
            }//end second loop over reco bins
        }//end loop over true bins
    }// end if !trueVert
    _hResponse = hResponse;
}//end makeResponseMatrix

TMatrixD Unfold::makeMatrixFromHist(TH2F*hist)
{
	int nBinsX = hist->GetNbinsX();
	int nBinsY = hist->GetNbinsY();
    if(nBinsX!=nBinsY){
        cout << endl;
        cout << "x and y bins do not match" << endl;
        cout << "in Unfold::makeMatrixFromHist()" << endl;
        cout << endl;
    }

	TMatrixD matrix(nBinsY+2,nBinsX+2);
	for(int i=0;i<=nBinsX+1;i++){
		for(int j=0;j<=nBinsY+1;j++){
			matrix(j,i) = hist->GetBinContent(i,j);
		}
	}
	return matrix;
}//end makeMatrixFromHist

TVectorD Unfold::makeVectorFromHist(TH1F*hist)
{
	int nBins = hist->GetNbinsX();
	TVectorD vec(nBins+2);
	for(int i=0;i<=nBins+1;i++){
		vec(i) = hist->GetBinContent(i);
	}
	return vec;
}//end makeVectorFromHist

TH1F*Unfold::makeHistFromVector(TVectorD vec,TH1F*hist)
{
	TH1F*hReturn = (TH1F*)hist->Clone("hUnfolded");
	int nBins = vec.GetNrows();
	for(int i=0;i<nBins+1;i++){
		hReturn->SetBinContent(i,vec(i));
	}
	return hReturn;
}//end makeHistFromVector

void Unfold::unfoldInversion()
{
	std::cout << endl;
	std::cout << "**************************************************" << endl;
	std::cout << "Beginning unfolding process using Matrix inversion" << endl;
	std::cout << "**************************************************" << endl;
	std::cout << endl;

	TH1F*hist1_oldBinning = (TH1F*)_hReco->Clone();
	TH1F*hist2 = (TH1F*)_hTrue->Clone();
	TH2F*hist3_oldBinning = (TH2F*)_hResponse->Clone();

	//For inversion method, we need the matrix to be square
	//So here the reco and matrix are rebinned
	TH1F*hist1 = RebinTH1(hist1_oldBinning,"hRecoRebin",_hTrue);
	TH2F*hist3 = RebinTH2(hist3_oldBinning,"hMatrixRebin",_hTrue);
	int nBins = hist2->GetNbinsX();

	//Make Vy, input covariance assuming for now that it is diagonal
	//This is used to calculate the output covariance later to get errors
	TMatrixD Vy(nBins+2,nBins+2);
	double entry;
	for(int i=0;i<=nBins;i++){
		for(int j=0;j<=nBins;j++){
			if(j==i) entry = (hist1->GetBinError(i))*(hist1->GetBinError(i));
			else entry = 0;
			Vy(i,j) = entry;	
		}
	}
	
	int nBinsTrue = hist2->GetNbinsX();
	int nBinsReco = hist1->GetNbinsX();
	double conditionNumber = _condition;
	cout << "Condition number: " << conditionNumber << endl;
	cout << "Number of output bins: " << nBinsTrue << endl;
	cout << "Number of input bins: " << nBinsReco << endl;

	//Turn histograms into matrices and vectors
	TMatrixD responseM = makeMatrixFromHist(hist3);
	TVectorD trueV = makeVectorFromHist(hist2);
	TVectorD recoV = makeVectorFromHist(hist1);

	//transpose response matrix to mutliply by the reco vector
	responseM.T();

	//Invert the response matrix
	TMatrixD invertedM = responseM.Invert();

	//Multiply inverse Matrix to reconstructed vector to get unfolded vector
	TVectorD unfoldedV = invertedM*recoV;

	TMatrixD invertedMT = invertedM.T();
	TMatrixD Vx_temp = Vy*invertedMT;
	TMatrixD Vx = invertedM*Vx_temp;

	TCanvas*c10 = new TCanvas("c10","",0,0,1000,1000);
	c10->SetGrid();
	Vx.Draw("colz");
	c10->SaveAs("plots/inversionCovarianceTest.png");
	TH1F*hUnfolded = makeHistFromVector(unfoldedV,_hTrue);
	TH1F*hUnfoldedE = (TH1F*)hUnfolded->Clone("hUnfoldedInversion");

	//Get error bars from diagonal of covariance matrix for unfolded histogram
	double binError;
	for(int i=0;i<=nBinsTrue;i++){
		binError = TMath::Sqrt(Vx(i,i));
		hUnfoldedE->SetBinError(i,binError);
		cout << "bin: " << i << ", % error: " << 100*binError/hUnfoldedE->GetBinContent(i) << endl;
	}
    _hUnfolded = hUnfoldedE;
}//end unfoldInversion

void Unfold::plotMatrix(TH2F*hMatrix,TString saveName,bool printCondition)
{
	gStyle->SetPalette(1);
	double xPosition,yPosition;
	TLatex*conditionLabel;
	TCanvas*canvas = new TCanvas("canvas","",0,0,1000,1000);
	canvas->SetGrid();
	canvas->SetLogy();
	canvas->SetLogx();
	canvas->SetRightMargin(0.15);
	canvas->SetLeftMargin(0.15);
	hMatrix->Draw("colz");
	
	if(printCondition){ _condition = GetConditionNumber();
		xPosition = 15;
		yPosition = 5;
		conditionLabel = new TLatex(xPosition,yPosition,
                                	    Form("condition number = %lg",_condition));
		conditionLabel->Draw("same");
	}

	TString saveDirectory = "./plots/";
	TString save = saveDirectory+saveName+".png";
	canvas->SaveAs(save);	
	delete canvas;
}

void Unfold::SetBackground(TH1F*hist)
{
    _hBack = hist;
    _backgroundSubtraction = true;
}

void Unfold::SetMatrix(TH2F*hist)
{
    _hMatrix = hist;
    makeResponseMatrix(_hMatrix); // normalized migration matrix
    _hResponseSquare = RebinTH2(_hResponse,"hResponseSquare",_hTrue);
    int nBinsX = _hMatrix->GetNbinsX();
    int nBinsY = _hMatrix->GetNbinsY();

    if(nBinsY == _nBinsTrue) _trueVert = true;
    else if(nBinsY == _nBinsReco) _trueVert = false;

    if(_trueVert){
        _nBinsReco = hist->GetNbinsX();
        _nBinsTrue = hist->GetNbinsY();
    }
    else{
        _nBinsReco = hist->GetNbinsY();
        _nBinsTrue = hist->GetNbinsX();
    }
}

void Unfold::SetReco(TH1F*hist)
{
    _hReco = hist;
	_nBinsReco = hist->GetNbinsX();
}

void Unfold::SetTrue(TH1F*hist)
{
    _hTrue = hist;
	_nBinsTrue = hist->GetNbinsY();
}

void Unfold::SetRegularizationType(RegType regType)
{
    _regType = regType;
}

double Unfold::ReturnCondition()
{
    return _condition;
}

TH1F*Unfold::ReturnUnfolded()
{
    return _hUnfolded;
}

TH2F*Unfold::ReturnResponseMatrix()
{
    return _hResponse;
}

TH1F*Unfold::ReturnReco()
{
    return _hReco;
}

TH1F*Unfold::ReturnTrue()
{
    return _hTrue;
}

TH1F*Unfold::ReturnBackground()
{
    if(_hBack == NULL)
        cout << "Error in ReturnBackground(): A background was not specified for unfolding" << endl;
    return _hBack;
}

TH2F*Unfold::ReturnMatrix()
{
    return _hMatrix;
}

void Unfold::DoUnfold(UnfoldType unfType)
{
    if(unfType == TUNFOLD) unfoldTUnfold();
    else if(unfType == INVERSION) unfoldInversion();
}
