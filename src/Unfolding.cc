#include "../include/Unfolding.hh"

Unfold::Unfold()
{

}

TH1F*Unfold::unfoldTUnfold(RegType regType,TH1F*hReco,TH1F*hTrue,TH2F*hMatrix)
{
	cout << endl;
	cout << "*****************************************" << endl;
	cout << "Beginning unfolding process using TUnfold" << endl;
	cout << "*****************************************" << endl;
	cout << endl;

	TH1F*hBlank = new TH1F("hBlank","",1,0,1);
	int nBinsReco = hReco->GetNbinsX();
	int nBinsTrue = hTrue->GetNbinsX();
	
	if(nBinsReco<=nBinsTrue){
		cout << "For TUnfold, the observed histogram must have more bins than the true histogram" << endl;
		cout << "Input bins: " << nBinsReco << endl;
		cout << "Output bins: " << nBinsTrue << endl;
		cout << "Closing unfolding script" << endl;
		return hBlank;
	}

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
	TUnfold::EHistMap outputMap = TUnfold::kHistMapOutputVert;
	//TUnfold::EHistMap outputMap = TUnfold::kHistMapOutputHoriz;
	
	//////////////////////////////////////
	//  Constructor for TUnfoldDensity  //
	//////////////////////////////////////
	TUnfoldDensity unfold(hMatrix,outputMap,regMode,constraintMode,densityFlags);
	unfold.SetInput(hReco);//the measured distribution

	///////////////////////
	//  Begin Unfolding  //
	///////////////////////
	Int_t iBest;
	TSpline *logTauX,*logTauY;
	TGraph *lCurve;
	TGraph*bestLcurve;
	TGraph*bestLogTauLogChi2;
	if(regType == VAR_REG_LCURVE){
		Int_t nScan=30;//This number chosen only because it was given in the tutorial
		Double_t tauMin = 0.0;//If tauMin=tauMax, TUnfold automatically chooses a range
		Double_t tauMax = 0.0;//Not certain how they choose the range
		iBest=unfold.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
		cout<< "tau=" << unfold.GetTau() << endl;
		cout<<"chi**2="<<unfold.GetChi2A()<<"+"<<unfold.GetChi2L()<<" / "<<unfold.GetNdf()<<"\n";
		Double_t t[1],x[1],y[1];
		logTauX->GetKnot(iBest,t[0],x[0]);
		logTauY->GetKnot(iBest,t[0],y[0]);
		bestLcurve=new TGraph(1,x,y);
		bestLogTauLogChi2=new TGraph(1,t,x);
	}
	else if(regType == VAR_REG_SCANSURE){
		cout << "Option VAR_REG_SCANSURE not yet implemented" << endl;
		return hBlank;
	}
	else if(regType == VAR_REG_SCANTAU){
		cout << "Option VAR_REG_SCANTAU not yet implemented" << endl;
		return hBlank;
	}
	else if(regType == NO_REG){
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
	for(int i=0;i<=nBinsTrue;i++){
		double binError = TMath::Sqrt(histEmatTotal->GetBinContent(i+1,i+1));
  		hUnfoldedE->SetBinError(i+1,binError);
 	}

	TFile*errorFile = new TFile("data/errorMatrixTUnfold.root","recreate");
	hUnfoldedE->Write();
	histEmatStat->Write();
	histEmatTotal->Write();
	errorFile->Close();

	return hUnfoldedE;
}//end unfoldTUnfold

TCanvas*Unfold::plotUnfolded(TString canvasName,TString titleName,TH1F*hReco,TH1F*hTrue,
			     TH1F*hUnfolded)
{
	gStyle->SetOptStat(0);
	//Get parameters for histograms
	//These are used to define the canvases and pads
	//as well as the locations of drawn objects
	int nBinsTrue = hTrue->GetNbinsX();
	int nBinsReco = hReco->GetNbinsX();
	int binLow = hTrue->GetBinLowEdge(1);
	int binHigh = hTrue->GetBinLowEdge(nBinsTrue)+hTrue->GetBinWidth(nBinsTrue);
	float peakMax = 0;	
	float binContent;
	for(int i=1;i<=nBinsTrue;i++){
		binContent = hTrue->GetBinContent(i);
		if(binContent > peakMax) peakMax = binContent;
	}

	//For TUnfold, nBinsReco = 2*nBinsTrue but we want to plot them all together
	//So to look nice, we rebin nBinsReco to match the binning of the true distribution
	TH1F*hRecoRebin2;
	hRecoRebin2 = RebinTH1(hReco,"hRecoRebin2",hTrue);

	//set histogram drawing options
	hTrue->SetFillColor(kRed+2);
	hTrue->SetLineColor(kRed+2);
	hRecoRebin2->SetMarkerStyle(20);
	hRecoRebin2->SetMarkerColor(kBlack);
	hRecoRebin2->SetLineColor(kBlack);
	hUnfolded->SetMarkerStyle(25);
	hUnfolded->SetMarkerColor(kBlue+2);
	hUnfolded->SetLineColor(kBlue+2);
	hUnfolded->SetFillColor(kWhite);

	//define the ratio plot to easily see how well the unfolded result matches the true one
	TH1F*ratio = (TH1F*)hUnfolded->Clone("ratio");
	ratio->Divide(hTrue);
	ratio->SetTitle("");


	TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
	legend->SetTextSize(0.02);
	legend->AddEntry(hTrue,"True Distribution");
	legend->AddEntry(hRecoRebin2,"Observed Distribution");
	legend->AddEntry(hUnfolded,"Unfolded Distribution");

	//Create a label that shows the chi^2 value to print on graph
	float xMax = hTrue->GetXaxis()->GetXmax();
	float yMax = 1.1*peakMax;
	double xChiLabel = xMax*0.70;
	double yChiLabel = yMax*0.75;
	double x[nBinsTrue],res[nBinsTrue];
	double chi = hUnfolded->Chi2Test(hTrue,"CHI2/NDF",res);//chi2/ndf to print on plot
	double pValues = hUnfolded->Chi2Test(hTrue,"P",res);//outputs chi2,prob,ndf,igood
	TLatex*chiLabel = new TLatex(xChiLabel,yChiLabel,Form("#chi^{2}/ndf = %lg", chi));

	//Draw canvas and pads to make plot
	TCanvas*canvas = new TCanvas(canvasName,"",0,0,1000,1000);
	const float padmargins = 0.03;
	const float yAxisMinimum = 0.1;
	const float yAxisMaximum = peakMax*1.1;
	TPad*pad1 = new TPad("","",0,0.3,1.0,1.0);
	pad1->SetBottomMargin(padmargins);
	pad1->SetGrid();
	pad1->SetTicks(1,1);
	pad1->Draw();
	pad1->cd();
	hTrue->SetLabelSize(0);
	hTrue->SetTitleSize(0);
	hTrue->SetMinimum(yAxisMinimum);
	hTrue->SetMaximum(yAxisMaximum);
	hTrue->SetTitle(titleName);
	hTrue->Draw("hist");
	hRecoRebin2->Draw("pe,same");
	hUnfolded->Draw("pe,same");	
	legend->Draw("same");
	chiLabel->Draw("same");

	canvas->cd();
	TPad*pad2 = new TPad("","",0,0.05,1,0.3);
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

double Unfold::GetConditionNumber(TH2F*hResponse)
{
	//This function obtains the condition number of the response matrix
	//This is important because it tells us about how diagonal the matrix is
	//which gives an idea of how well unfolding will work
	//This helps us figure out if we need to use regularization or not
	//Small enough condition numbers mean no regularization is needed
	
	TString histName = hResponse->GetName();
	int nBinsX = hResponse->GetNbinsX();
	int nBinsY = hResponse->GetNbinsY();
	TMatrixD matrix = makeMatrixFromHist(hResponse);

	TDecompSVD decomp(matrix);
	double condition = decomp.Condition();
	cout << "The condition number for " << histName << ": " << condition << endl;

	double determinant;
	TMatrixD mInverse = matrix.Invert(&determinant);
	cout << "The determinant of " << histName << " is " << determinant << endl;
	return condition;
}//end GetConditionNumber

TH2F*Unfold::makeResponseMatrix(TH2F*hist)
{
	TH2F*hResponse = (TH2F*)hist->Clone("hResponse");
	int nBinsX = hResponse->GetNbinsX();
	int nBinsY = hResponse->GetNbinsY();
	double nEntriesX;
	double binContent;
	double scaledContent;
	
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

	return hResponse;
}//end makeResponseMatrix

TMatrixD Unfold::makeMatrixFromHist(TH2F*hist)
{
	int nBinsX = hist->GetNbinsX();
	int nBinsY = hist->GetNbinsY();
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

TH1F*Unfold::unfoldInversion(TH1F*hReco,TH1F*hTrue,TH2F*hMatrix)
{
	std::cout << endl;
	std::cout << "**************************************************" << endl;
	std::cout << "Beginning unfolding process using Matrix inversion" << endl;
	std::cout << "**************************************************" << endl;
	std::cout << endl;

	TH1F*hist1_oldBinning = (TH1F*)hReco->Clone();
	TH1F*hist2 = (TH1F*)hTrue->Clone();
	TH2F*hist3_oldBinning = (TH2F*)hMatrix->Clone();


	//For inversion method, we need the matrix to be square
	//So here they are rebinned
	TH1F*hist1 = RebinTH1(hist1_oldBinning,"hRecoRebin",hTrue);
	TH2F*hist3 = RebinTH2(hist3_oldBinning,"hMatrixRebin",hTrue);

	TH2F*hNorm = makeResponseMatrix(hist3);
	int nBinsTrue = hist2->GetNbinsX();
	int nBinsReco = hist1->GetNbinsX();
	double conditionNumber = GetConditionNumber(hNorm);
	cout << "Condition number: " << conditionNumber << endl;
	cout << "Number of output bins: " << nBinsTrue << endl;
	cout << "Number of input bins: " << nBinsReco << endl;

	//Turn histograms into matrices and vectors
	TMatrixD responseM = makeMatrixFromHist(hNorm);
	TVectorD trueV = makeVectorFromHist(hist2);
	TVectorD recoV = makeVectorFromHist(hist1);

	//transpose response matrix to mutliply by the reco vector
	responseM.T();

	//Invert the response matrix
	TMatrixD invertedM = responseM.Invert();

	//Multiply inverse Matrix to reconstructed vector to get unfolded vector
	TVectorD unfoldedV = invertedM*recoV;

	TMatrixD invertedMT = invertedM.T();
	TMatrixD Vx = invertedM*invertedMT;

	TH1F*hUnfolded = makeHistFromVector(unfoldedV,hTrue);
	TH1F*hUnfoldedE = (TH1F*)hUnfolded->Clone("hUnfoldedInversion");

	//Get error bars from diagonal of covariance matrix for unfolded histogram
	TFile*loadError = new TFile("data/errorMatrixTUnfold.root");
	TH2D*hError = (TH2D*)loadError->Get("hEmatrixTotal");
	double binError;
	for(int i=1;i<=nBinsTrue;i++){
		binError = TMath::Sqrt(hError->GetBinContent(i,i));
		hUnfoldedE->SetBinError(i,binError);
	}
	return hUnfoldedE;
}//end unfoldInversion

void Unfold::plotMatrix(TH2F*hMatrix,TString saveName,bool printCondition)
{
	gStyle->SetPalette(1);
	double xPosition,yPosition;
	TLatex*conditionLabel;
	TCanvas*canvas = new TCanvas("canvas","",0,0,1000,1000);
	canvas->SetGrid();
	canvas->SetRightMargin(0.15);
	canvas->SetLeftMargin(0.15);
	hMatrix->Draw("colz");
	

	if(printCondition){ condition = GetConditionNumber(hMatrix);
		xPosition = 15;
		yPosition = 5;
		conditionLabel = new TLatex(xPosition,yPosition,
                                	    Form("condition number = %lg",condition));
		conditionLabel->Draw("same");
	}

	TString saveDirectory = "./plots/";
	TString save = saveDirectory+saveName+".png";
	canvas->SaveAs(save);	
	delete canvas;
}

TH2F*Unfold::makeResponseMatrixT(TH2F*hist)
{
	//Here we normalize for each reco bin instead of each true bin
	//This is to investigate how events migrate into a given reco bin
	TH2F*hResponseT = (TH2F*)hist->Clone("hResponseT");
	int nBinsX = hResponseT->GetNbinsX();
	int nBinsY = hResponseT->GetNbinsY();
	double nEntriesY;
	double binContent;
	double scaledContent;
	
	//Loop over all true bins (the y-axis)
	for(int i=1;i<=nBinsX;i++){
		nEntriesY = 0.0;
		//for each reco bin, sum up the number of events across all reco bins
		for(int j=1;j<=nBinsY;j++){
			binContent = hist->GetBinContent(i,j);
			nEntriesY += binContent;
		}//end first loop over reco bins
		//For each true bin scale the bin content by the number of entries in all
		//reco bins, then place this content into the new matrix
		for(int j=1;j<=nBinsY;j++){
			scaledContent = hist->GetBinContent(i,j)/nEntriesY;
			hResponseT->SetBinContent(i,j,scaledContent);
		}//end second loop over reco bins
	}//end loop over true bins

	return hResponseT;
}//end makeResponseMatrix
