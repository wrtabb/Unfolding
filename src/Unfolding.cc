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
	//TUnfold::EConstraint constraintMode = TUnfold::kEConstraintNone;
	TUnfold::EConstraint constraintMode = TUnfold::kEConstraintArea;

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

	//Create unfolding histogram with errors
	TH1F*hUnfoldedE = (TH1F*)hUnfolded->Clone("hUnfoldedTUnfold");

	//loop over unfolded histogram bins and assign errors to each one
	for(int i=0;i<=nBinsTrue;i++){
		double binError = TMath::Sqrt(histEmatTotal->GetBinContent(i+1,i+1));
  		hUnfoldedE->SetBinError(i+1,binError);
 	}
	return hUnfoldedE;
}//end unfoldTUnfold

void Unfold::plotUnfolded(TH1F*hReco,TH1F*hTrue,TH1F*hUnfoldedE,UnfoldType unfoldType,
			  bool closure)
{
	gStyle->SetOptStat(0);
	gROOT->SetBatch(true);
	//if hReco is from TUnfold, it has more bins that hTrue and needs to be 
	//rebinned for plotting
	TH1F*hRecoRebin = (TH1F*)hReco->Clone("hRecoRebin");
	//need to update custom rebinning function to calculate error and use it here
	if(unfoldType==TUNFOLD) hRecoRebin->Rebin(2);; 

	hUnfoldedE->SetMarkerStyle(25);
	hUnfoldedE->SetMarkerColor(kBlue+2);
	hUnfoldedE->SetMarkerSize(1);
	hUnfoldedE->SetFillColor(kWhite);
	hRecoRebin->SetMarkerStyle(20);
	hRecoRebin->SetMarkerColor(kBlack);
	hRecoRebin->SetLineColor(kBlack);
	hRecoRebin->SetFillColor(kWhite);
	hTrue->SetFillColor(kRed+2);
	hTrue->SetLineColor(kRed+2);
	hTrue->SetTitle("");
	hTrue->SetLabelSize(0);
	hTrue->SetTitleSize(0);

	//Ratio of unfolded histogram over true histogram
	TH1F*ratio = (TH1F*)hUnfoldedE->Clone("ratio");
	ratio->Divide(hTrue);

	double xChiLabel = 35;
	double yChiLabel = 5e5;
	double x[nBinsTrue],res[nBinsTrue];
	double chi = hUnfoldedE->Chi2Test(hTrue,"CHI2/NDF",res);//chi2/ndf to print on plot
	TLatex*chiLabel = new TLatex(xChiLabel,yChiLabel,Form("#chi^{2}/ndf = %lg", chi));

	const float padmargins = 0.03;
	TCanvas*canvas1 = new TCanvas("canvas1","",10,10,1200,1000);
	TPad*pad1 = new TPad("","",0,0.3,1.0,1.0);
	pad1->SetBottomMargin(padmargins);
	pad1->SetGrid();
	pad1->SetTicks(1,1);
	pad1->Draw();
	pad1->cd();
	double binLow = hUnfoldedE->GetXaxis()->GetBinLowEdge(1);
	double binHigh = hUnfoldedE->GetXaxis()->GetBinUpEdge(hUnfoldedE->GetNbinsX());
	TLine*line = new TLine(binLow,1,binHigh,1);
	line->SetLineColor(kRed);
	TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
	legend->SetTextSize(0.02);
	legend->AddEntry(hTrue,"True Distribution");
	legend->AddEntry(hRecoRebin,"Measured Distribution");
	legend->AddEntry(hUnfoldedE,"Unfolded Distribution");
	hTrue->Draw("hist");
	hRecoRebin->Draw("PE,same");
	hUnfoldedE->Draw("PE,same");
	legend->Draw("same");
	chiLabel->Draw("same");

	canvas1->cd();
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
	ratio->SetMarkerStyle(20);
	ratio->SetMarkerColor(kBlack);
	line->Draw();
	ratio->Draw("PE,same");
	TString saveName = "plots/unfolded";
	if(unfoldType==TUNFOLD) saveName += "TUnfold";
	else if(unfoldType==INVERSION) saveName += "Inversion";
	if(closure) saveName += "ClosureTest";
	saveName += "_NoRegularization_Mean";
	saveName += (int)mean;
	saveName += "_RecoBins";
	saveName += (int)nBinsReco;
	saveName += ".png";
	canvas1->SaveAs(saveName);
	delete canvas1;
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

TH2F*Unfold::RebinTH2(TH2F*hist,TString histName,TH2F*hBinning)
{
	int nBinsHist = hist->GetNbinsY();
	int nBinsReco = hBinning->GetNbinsY();
	double newbinning[nBinsReco];
	double binningTrue[nBinsTrue];
	for(int i=0;i<=nBinsTrue;i++){
		if(i==nBinsTrue) binningTrue[i] = hist->GetXaxis()->GetBinUpEdge(i);
		else binningTrue[i] = hist->GetXaxis()->GetBinLowEdge(i+1);
	}

	TH1D*hBinningProjection = hBinning->ProjectionY();
	for(int i=0;i<=nBinsReco;i++){
		if(i==0) newbinning[i] = hBinningProjection->GetBinLowEdge(i+1);
		else newbinning[i] = newbinning[i-1]+hBinningProjection->GetBinWidth(i);
	}
	TH2F*hRebin = new TH2F(histName,"",nBinsTrue,binningTrue,nBinsReco,newbinning);

	double y,x;
	for(int i=1;i<=nBinsTrue;i++){
		for(int j=1;j<=nBinsHist;j++){
			x = hist->GetXaxis()->GetBinCenter(i);
			y = hist->GetYaxis()->GetBinCenter(j);
			int nEntries = hist->GetBinContent(i,j);
			for(int k=0;k<nEntries;k++){
				hRebin->Fill(x,y);
			}//end filling hRebin
		}//end y bin loop
	} //end x bin loop
	return hRebin;
}

TH2F*Unfold::RebinTH2(TH2F*hist,TString histName,std::vector<double> binning)
{
 int nBinsHist = hist->GetNbinsY();
 int nBinsReco = binning.size()-1;
 double newbinning[nBinsReco];
 double binningTrue[nBinsTrue];
 for(int i=0;i<=nBinsTrue;i++){
  if(i==nBinsTrue) binningTrue[i] = hist->GetXaxis()->GetBinUpEdge(i);
  else binningTrue[i] = hist->GetXaxis()->GetBinLowEdge(i+1);
 }
 for(int i=0;i<=nBinsReco;i++){
  newbinning[i] = binning.at(i);
 }
 TH2F*hRebin = new TH2F(histName,"",nBinsTrue,binningTrue,nBinsReco,newbinning);

 double y,x;
 for(int i=1;i<=nBinsTrue;i++){
  for(int j=1;j<=nBinsHist;j++){
   x = hist->GetXaxis()->GetBinCenter(i);
   y = hist->GetYaxis()->GetBinCenter(j);
   int nEntries = hist->GetBinContent(i,j);
   for(int k=0;k<nEntries;k++){
    hRebin->Fill(x,y);
   }//end filling hRebin
  }//end y bin loop
 } //end x bin loop
 return hRebin;
}

TH1F*Unfold::RebinTH1(TH1F*hist,TString histName,std::vector<double> binning)
{
 int nBinsHist = hist->GetNbinsX();
 int nBinsReco = binning.size()-1;
 double newbinning[nBinsReco];
 for(int i=0;i<=nBinsReco;i++){
  newbinning[i] = binning.at(i);
 }
 TH1F*hRebin = new TH1F(histName,"",nBinsReco,newbinning);

 int bin;
 double y,x;
 for(int j=1;j<=nBinsHist;j++){
  x = hist->GetXaxis()->GetBinCenter(j);
  int nEntries = hist->GetBinContent(j);
  for(int k=0;k<nEntries;k++){
   hRebin->Fill(x);
  }//end filling hRebin
 } //end x bin loop
 return hRebin;
}

TH1F*Unfold::RebinTH1(TH1F*hist,TString histName,TH1F*hBinning)
{
 //This function takes a histogram,hist, and defines its binning according to a given 
 //histogram, hBinning
 //Currently, it does not keep the errors, so this still needs to be added.
 int nBinsHist = hist->GetNbinsX();
 int nBinsNew = hBinning->GetNbinsX();
 double newbinning[nBinsNew];
 for(int i=0;i<=nBinsNew;i++){
  if(i==0) newbinning[i] = hBinning->GetBinLowEdge(i+1);
  else newbinning[i] = newbinning[i-1]+hBinning->GetBinWidth(i);
 }
 TH1F*hRebin = new TH1F(histName,"",nBinsNew,newbinning);
 
 double y,x;
 double nEntries;
 for(int j=1;j<=nBinsHist;j++){
  x = hist->GetXaxis()->GetBinCenter(j);
  nEntries = hist->GetBinContent(j);
  hRebin->Fill(x,nEntries);
 } //end x bin loop
 return hRebin;
}

