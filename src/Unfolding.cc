#include "../include/Unfolding.hh"

Unfold::Unfold()
{

}

TH1F* Unfold::unfoldTUnfold(RegType regType,TH1F*hReco,TH1F*hTrue,TH2F*hMatrix)
{
	cout << endl;
	cout << "*****************************************" << endl;
	cout << "Beginning unfolding process using TUnfold" << endl;
	cout << "*****************************************" << endl;
	cout << endl;

	TH1F*hBlank = new TH1F("hBlank","",1,0,1);
	int nBinsReco = hReco->GetNbinsX();
	int nBinsTrue = hTrue->GetNbinsX();

	if(nBinsReco==nBinsTrue){
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
	TH1*hUnfolded = unfold.GetOutput("Unfolded");

	//Create error matrices
	TH2*histEmatStat=unfold.GetEmatrixInput("unfolding stat error matrix");
	TH2*histEmatTotal=unfold.GetEmatrixTotal("unfolding total error matrix");

	TFile*errorFile = new TFile("data/errorMatrixTUnfold.root","recreate");
	histEmatStat->Write();
	histEmatTotal->Write();
	errorFile->Close();

	//Create unfolding histogram with errors
	TH1F*hUnfoldedE = (TH1F*)hUnfolded->Clone("hUnfoldedTUnfold");

	//loop over unfolded histogram bins and assign errors to each one
	for(int i=0;i<=nBinsTrue;i++){
		double binError = TMath::Sqrt(histEmatTotal->GetBinContent(i+1,i+1));
  		hUnfoldedE->SetBinError(i+1,binError);
 	}
	return hUnfoldedE;
}//end unfoldTUnfold

TCanvas*Unfold::plotUnfolded(TString canvasName,TH1F*hReco,TH1F*hTrue,TH1F*hUnfolded)
{
	gStyle->SetOptStat(0);
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

	TH1F*hRecoRebin2;
	if(nBinsReco!=nBinsTrue) hRecoRebin2 = RebinTH1(hReco,"hRecoRebin2",hTrue);
	else hRecoRebin2 = (TH1F*)hReco->Clone();

	hTrue->SetFillColor(kRed+2);
	hTrue->SetLineColor(kRed+2);
	hRecoRebin2->SetMarkerStyle(20);
	hRecoRebin2->SetMarkerColor(kBlack);
	hRecoRebin2->SetLineColor(kBlack);
	hUnfolded->SetMarkerStyle(25);
	hUnfolded->SetMarkerColor(kBlue+2);
	hUnfolded->SetLineColor(kBlue+2);

	TH1F*ratio = (TH1F*)hUnfolded->Clone("ratio");
	ratio->Divide(hTrue);

	TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
	legend->SetTextSize(0.02);
	legend->AddEntry(hTrue,"True Distribution");
	legend->AddEntry(hReco,"Observed Distribution");
	legend->AddEntry(hUnfolded,"Unfolded Distribution");

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
	hTrue->SetTitle("unfolded results");
	hTrue->Draw("hist");
	hRecoRebin2->Draw("pe,same");
	hUnfolded->Draw("pe,same");	
	legend->Draw("same");

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
	for(int i=1;i<=nBinsX;i++){
		double nEntriesX = 0;
		for(int j=1;j<=nBinsY;j++){
			nEntriesX += hResponse->GetBinContent(i,j);
		}
		double sum = 0;
		for(int j=1;j<=nBinsY;j++){
			double scaledContent = hResponse->GetBinContent(i,j)/nEntriesX;
			hResponse->SetBinContent(i,j,scaledContent);
			sum += scaledContent;
		}
	}
	return hResponse;
}//end makeResponseMatrix

TMatrixD Unfold::makeMatrixFromHist(TH2F*hist)
{
	int nBinsX = hist->GetNbinsX();
	int nBinsY = hist->GetNbinsY();
	TMatrixD matrix(nBinsY,nBinsX);
	for(int i=1;i<=nBinsX;i++){
		for(int j=1;j<=nBinsY;j++){
			matrix(j-1,i-1) = hist->GetBinContent(i,j);
		}
	}
	return matrix;
}//end makeMatrixFromHist

TVectorD Unfold::makeVectorFromHist(TH1F*hist)
{
	int nBins = hist->GetNbinsX();
	TVectorD vec(nBins);
	for(int i=1;i<=nBins;i++){
		vec(i-1) = hist->GetBinContent(i);
	}
	return vec;
}//end makeVectorFromHist

TH1F*Unfold::makeHistFromVector(TVectorD vec,TH1F*hist)
{
	TH1F*hReturn = (TH1F*)hist->Clone("hUnfolded");
	int nBins = vec.GetNrows();
	for(int i=0;i<nBins;i++){
		hReturn->SetBinContent(i+1,vec(i));
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

	TString unfoldType = "Inversion";
	TH2F*hResponse = makeResponseMatrix(hMatrix);
	int nBinsTrue = hTrue->GetNbinsX();
	int nBinsReco = hReco->GetNbinsX();
	double conditionNumber = GetConditionNumber(hResponse);
	cout << "Condition number: " << conditionNumber << endl;
	cout << "Number of output bins: " << nBinsTrue << endl;
	cout << "Number of input bins: " << nBinsReco << endl;

	//Turn histograms into matrices and vectors
	TMatrixD responseM = makeMatrixFromHist(hResponse);
	TVectorD trueV = makeVectorFromHist(hTrue);
	TVectorD recoV = makeVectorFromHist(hReco);

	//Invert
	TMatrixD invertedM = responseM.Invert();
	TVectorD unfoldedV = invertedM*recoV;

	TMatrixD invertedMT = invertedM.T();
	TMatrixD Vx = invertedM*invertedMT;

	TH1F*hUnfolded = makeHistFromVector(unfoldedV,hTrue);
	TH1F*hUnfoldedE = (TH1F*)hUnfolded->Clone("hUnfoldedInversion");

	//Get error bars from diagonal of covariance matrix for unfolded histogram
	for(int i=0;i<nBinsTrue;i++){
		hUnfoldedE->SetBinError(i+1,TMath::Sqrt(Vx(i,i)));
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


TH1F*Unfold::RebinTH1(TH1F*hist,TString histName,TH1F*hBinning)
{
	int nBinsOld = hist->GetNbinsX();
	int nBinsNew = hBinning->GetNbinsX();
	if(nBinsNew > nBinsOld){
		cout << "*********************************************************" << endl;
		cout << "ERROR: new binning must have fewer bins than old binning!" << endl;
		cout << "*********************************************************" << endl;
		return hist;
	}
	double newbinning[nBinsNew];
	for(int i=0;i<=nBinsNew;i++){
		if(i==0) newbinning[i] = hBinning->GetBinLowEdge(i+1);
		else newbinning[i] = newbinning[i-1]+hBinning->GetBinWidth(i);
	}
	TH1F*hRebin = new TH1F(histName,"",nBinsNew,newbinning);
	double y,x;
	double nEntries;
	double histErrors[nBinsOld+2];
	for(int j=1;j<=nBinsOld;j++){
		histErrors[j] = hist->GetBinError(j);
		x = hist->GetXaxis()->GetBinCenter(j);
		nEntries = hist->GetBinContent(j);
		hRebin->Fill(x,nEntries);
		hRebin->GetBin(x);
	} //end x bin loop

	double histBinWidth;
	double newBinWidth;
	double nBinsOldInNew[nBinsNew];
	double newBinUpperEdge,oldBinUpperEdge;
	int k = 1;
	for(int i=1;i<=nBinsNew;i++){
		newBinUpperEdge = hRebin->GetBinLowEdge(i)+hRebin->GetBinWidth(i);
		double newError2 = 0;
		for(int j=k;j<=nBinsOld;j++){
			oldBinUpperEdge = hist->GetBinLowEdge(j)+hist->GetBinWidth(j);
			if(oldBinUpperEdge > newBinUpperEdge) continue;
			k = j;
			newError2 += hist->GetBinError(j)*hist->GetBinError(j);
		}
		hRebin->SetBinError(i,sqrt(newError2));
	}
	return hRebin;

}

