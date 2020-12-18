#include "../include/Unfolding.hh"

Unfold::Unfold()
{

}

//Makes model with a Gaussian peak on top of a 1/x function
//Choose desired mean, domain, and number of reco bins
//Function assumes reco bins are twice the number of true bins
void Unfold::makeToyModels(double mean,double xlow,double xhigh,int nBinsReco)
{
	cout << endl;
	cout << "*************************************" << endl;
	cout << "Creating toy models to test unfolding" << endl;
	cout << "*************************************" << endl;
 	cout << endl;
	
	int nBinsTrue = nBinsReco/2;
	//measured distribution
	TH1D*hReco =        new TH1D("hReco","",nBinsReco,xlow,xhigh);
	//measured distribution using same seed as migration matrix
 	TH1D*hRecoClosure = new TH1D("hClosure","",nBinsReco,xlow,xhigh);	
	//true distribution
	TH1D*hTrue =        new TH1D("hTrue","",nBinsTrue,xlow,xhigh);
	//2D matrix of true versus reco distributions
	TH2D*hMatrix =      new TH2D("hMatrix","",nBinsTrue,xlow,xhigh,nBinsReco,xlow,xhigh);

	//-----Parameters for toy model-----//
	const double sigma = 5;//Sigma of gaussian part of distribution
	const double mean_smeared = 0;//mean of smearing used to make reco
	const double sigma_smeared = 1.0;//sigma of smearing used to make reco
	const Long64_t nEntries = 1e7;//number of entries	

	//-----Function model-----//
	TF1*func = new TF1("func","1/(x+1)+gaus(0)",0,50);
	func->SetParameters(1.0,mean,sigma);

	TRandom3 gen1;
	TRandom3 gen2;
	gen1.SetSeed(82);
	gen2.SetSeed(1981);

	//-----Save for later retrieval-----//
	TFile*saveFile = new TFile("data/toyModelDistributions.root","recreate");	

	//-----Fill distributions with random numbers-----//
	double peak,peakReco,peak_smeared,peakReco_smeared;
	for(int i=0;i<nEntries;i++){
 		peak = func->GetRandom();
 		peakReco = func->GetRandom();
 		peak_smeared = peak+gen1.Gaus(mean_smeared,sigma_smeared);
 		peakReco_smeared = peakReco+gen2.Gaus(mean_smeared,sigma_smeared);

 		hRecoClosure->Fill(peak_smeared);
 		hReco->Fill(peakReco_smeared);
 		hTrue->Fill(peak);
 		hMatrix->Fill(peak,peak_smeared);
	}//end loop over entries

	//
	TH2D*hResponse = makeResponseMatrix(hMatrix);
	hMatrix->GetYaxis()->SetTitle("reco");
	hMatrix->GetXaxis()->SetTitle("true");
	hResponse->GetYaxis()->SetTitle("reco");
	hResponse->GetXaxis()->SetTitle("true");
	int meanName = mean;
	int binName = nBinsReco;
	TString matrixSaveName = "migrationMatrix_Mean";
	matrixSaveName += meanName;
	matrixSaveName += "_RecoBins";
	matrixSaveName += nBinsReco;
	plotMatrix(hMatrix,matrixSaveName,false);

	TString responseSaveName = "responseMatrix_Mean";
	responseSaveName += meanName;
	responseSaveName += "_RecoBins";
	responseSaveName += nBinsReco;
	plotMatrix(hResponse,responseSaveName,true);

	//-----Save the distributions-----//
	saveFile->cd();
	hReco->Write();
	hRecoClosure->Write();
	hTrue->Write();
	hMatrix->Write();
	hResponse->Write();
	saveFile->Close();
}//end makeToyModels

TH1F* Unfold::unfoldTUnfold(RegType regType,TH1D*hReco,TH1D*hTrue,TH2D*hMatrix)
{
	cout << endl;
	cout << "*****************************************" << endl;
	cout << "Beginning unfolding process using TUnfold" << endl;
	cout << "*****************************************" << endl;
	cout << endl;

	TH1F*hBlank = new TH1F("hBlank","",1,0,1);

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
	//TUnfold::EHistMap outputMap = TUnfold::kHistMapOutputVert;
	TUnfold::EHistMap outputMap = TUnfold::kHistMapOutputHoriz;
	
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
	TH1F*hUnfoldedE = (TH1F*)hUnfolded->Clone("Unfolded with errors");

	//loop over unfolded histogram bins and assign errors to each one
	for(int i=0;i<=nBinsTrue;i++){
		double binError = TMath::Sqrt(histEmatTotal->GetBinContent(i+1,i+1));
  		hUnfoldedE->SetBinError(i+1,binError);
 	}
	return hUnfoldedE;
}//end unfoldTUnfold

void Unfold::plotUnfolded(TH1D*hReco,TH1D*hTrue,TH1F*hUnfoldedE,UnfoldType unfoldType,
			  bool closure)
{
	gStyle->SetOptStat(0);
	gROOT->SetBatch(true);
	//if hReco is from TUnfold, it has more bins that hTrue and needs to be 
	//rebinned for plotting
	TH1D*hRecoRebin = (TH1D*)hReco->Clone("hRecoRebin");
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
	TH1D*ratio = (TH1D*)hUnfoldedE->Clone("ratio");
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
	//ratio->SetMinimum(0.7);
	//ratio->SetMaximum(1.3);
	//ratio->GetYaxis()->SetLabelSize(0.06);
	//ratio->GetYaxis()->SetTitleSize(0.08);
	//ratio->GetYaxis()->SetTitleOffset(0.3);
	//ratio->GetYaxis()->SetTitle("Unfolded/Truth");
	//ratio->GetXaxis()->SetLabelSize(0.1);
	//ratio->GetXaxis()->SetTitleSize(0.1);
	//ratio->SetMarkerStyle(20);
	//ratio->SetMarkerColor(kBlack);
	line->Draw();
	//ratio->Draw("PE,same");
	TString saveName = "plots/unfolded";
	if(unfoldType==TUNFOLD) saveName += "TUnfold";
	else if(unfoldType==INVERSION) saveName += "Inversion";
	if(closure) saveName += "ClosureTest";
	saveName += "_NoRegularization_Mean";
	//saveName += (int)mean;
	saveName += "_RecoBins";
	//saveName += (int)nBinsReco;
	saveName += ".png";
	canvas1->SaveAs(saveName);
	delete canvas1;
}//end plotUnfolded

double Unfold::GetConditionNumber(TH2D*hResponse)
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

TH2D*Unfold::makeResponseMatrix(TH2D*hist)
{
	TH2D*hResponse = (TH2D*)hist->Clone("hResponse");
	hResponse->RebinY(2);
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

TMatrixD Unfold::makeMatrixFromHist(TH2D*hist)
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

TVectorD Unfold::makeVectorFromHist(TH1D*hist)
{
	int nBins = hist->GetNbinsX();
	TVectorD vec(nBins);
	for(int i=1;i<=nBins;i++){
		vec(i-1) = hist->GetBinContent(i);
	}
	return vec;
}//end makeVectorFromHist

TH1F*Unfold::makeHistFromVector(TVectorD vec,TH1D*hist)
{
	TH1F*hReturn = (TH1F*)hist->Clone("hUnfolded");
	int nBins = vec.GetNrows();
	for(int i=0;i<nBins;i++){
		hReturn->SetBinContent(i+1,vec(i));
	}
	return hReturn;
}//end makeHistFromVector

TH1F*Unfold::unfoldInversion(TH1D*hReco,TH1D*hTrue,TH2D*hResponse)
{
	std::cout << endl;
	std::cout << "**************************************************" << endl;
	std::cout << "Beginning unfolding process using Matrix inversion" << endl;
	std::cout << "**************************************************" << endl;
	std::cout << endl;

	TString unfoldType = "Inversion";

	//Make response matrix histogram from input histogram matrix
	//Meaning we normalize each column (true bins) 
	//Here I am assuming the reco has twice as many bins as tru
	//This is because this is what I'm using for this specific task
	//This would have to be treated differently if this wasn't the case
	TH1D*hRecoRebin = (TH1D*)hReco->Clone("hRecoRebin");
	hRecoRebin->Rebin(2);//need to update this to use custom rebin function

	//Turn histograms into matrices and vectors
	TMatrixD responseM = makeMatrixFromHist(hResponse);
	TVectorD trueV = makeVectorFromHist(hTrue);
	TVectorD recoV = makeVectorFromHist(hRecoRebin);

	//Invert
	TMatrixD invertedM = responseM.Invert();
	TVectorD unfoldedV = invertedM*recoV;

	//Get covariance (assuming Vy = identity)
	TMatrixD invertedMT = invertedM.T();
	TMatrixD Vx = invertedM*invertedMT;

	TH1F*hUnfolded = makeHistFromVector(unfoldedV,hTrue);
	TH1F*hUnfoldedE = (TH1F*)hUnfolded->Clone("unfolded with errors");
	std::cout << "Number of bins = " << nBinsTrue << endl;

	//Get error bars from diagonal of covariance matrix for unfolded histogram
	for(int i=0;i<nBinsTrue;i++){
		hUnfoldedE->SetBinError(i+1,TMath::Sqrt(Vx(i,i)));
	}
	return hUnfoldedE;
}//end unfoldInversion

void Unfold::plotMatrix(TH2D*hMatrix,TString saveName,bool printCondition)
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

TH2*Unfold::RebinTH2(TH2*hist,TString histName,TH2*hBinning)
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
	TH2D*hRebin = new TH2D(histName,"",nBinsTrue,binningTrue,nBinsReco,newbinning);

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

TH2*Unfold::RebinTH2(TH2*hist,TString histName,std::vector<double> binning)
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
 TH2D*hRebin = new TH2D(histName,"",nBinsTrue,binningTrue,nBinsReco,newbinning);

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

TH1*Unfold::RebinTH1(TH1*hist,TString histName,std::vector<double> binning)
{
 int nBinsHist = hist->GetNbinsX();
 int nBinsReco = binning.size()-1;
 double newbinning[nBinsReco];
 for(int i=0;i<=nBinsReco;i++){
  newbinning[i] = binning.at(i);
 }
 TH1D*hRebin = new TH1D(histName,"",nBinsReco,newbinning);

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

TH1*Unfold::RebinTH1(TH1*hist,TString histName,TH1*hBinning)
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
 TH1D*hRebin = new TH1D(histName,"",nBinsNew,newbinning);
 
 double y,x;
 double nEntries;
 for(int j=1;j<=nBinsHist;j++){
  x = hist->GetXaxis()->GetBinCenter(j);
  nEntries = hist->GetBinContent(j);
  hRebin->Fill(x,nEntries);
 } //end x bin loop
 return hRebin;
}

