#include "../include/ToyModel.hh"

ToyModel::ToyModel(double distNorm,double power,double peakNormRel,double distMean,
		   double distSigma,double resSigma,double xmin,double xmax,int nBins)
{
	TH1::SetDefaultSumw2();
	SetModelParameters(distNorm,power,peakNormRel,distMean,distSigma,resSigma,xmin,xmax,
			   nBins);

	TF1*trueFunc = new TF1("trueFunction",trueDistribution,-10,250,5);	
	trueFunc->SetParameters(distNorm,power,peakNormRel,distMean,distSigma);
	
	TF1*recoFunc = new TF1("recoFunction",recoDistribution,-10,250,6);
	recoFunc->SetParameters(distNorm,power,peakNormRel,distMean,distSigma,resSigma);

	TF2*matrixFunc = new TF2("integrand2DFunc", integrand2D,-10,250,-10,250,6);
	matrixFunc->SetParameters(distNorm,power,peakNormRel,distMean,distSigma,resSigma);

	SetModelFunctions(trueFunc,recoFunc,matrixFunc);
}

void ToyModel::SetModelFunctions(TF1*trueFunc,TF1*recoFunc,TF2*matrixFunc)
{
	_trueFunc = trueFunc;
	_recoFunc = recoFunc;
	_matrixFunc = matrixFunc;
}

void ToyModel::SetModelParameters(double distNorm,double power,double peakNormRel,
				  double distMean,double distSigma,double resSigma,double xmin, 				 double xmax,int nBins)
{
	_distNorm = distNorm;
	_power = power;
	_peakNormRel = peakNormRel;
	_distMean = distMean;
	_distSigma = distSigma;
	_resSigma = resSigma;
	_xmin = xmin;
	_xmax = xmax;
	_nBins = nBins;
}

TH1F*ToyModel::GetTrueHist()
{
	TF1*func = _trueFunc;
	func->SetParameters(_distNorm,_power,_peakNormRel,_distMean,_distSigma,_resSigma);	
	TH1F*hist = new TH1F("hTrue","",_nBins,_xmin,_xmax);
	hist->SetFillColor(kRed+2);
	hist->SetLineColor(kRed+2);

	int nBinsTrue = _nBins;
	double nEntries;
	double binMin;
	double binMax;
	double binCenter;
	double xMin = _xmin-7*_resSigma;
	double xMax = _xmax+7*_resSigma;
	for(int i=0;i<=nBinsTrue;i++){
		if(i==0) binMin = xMin;
		else binMin = hist->GetXaxis()->GetBinLowEdge(i);
		if(i==nBinsTrue+1) binMax = xMax;
		else binMax = hist->GetXaxis()->GetBinUpEdge(i);
		binCenter = (binMax+binMin)/2.0;
		nEntries = func->Integral(binMin,binMax);
		hist->SetBinContent(i,nEntries);
		if(nEntries!=0) hist->SetBinError(i,TMath::Sqrt(nEntries));
	}

	return hist;
}


TH1F*ToyModel::GetRecoHist()
{
	int nBinsTrue = _nBins;
	int nBinsReco = 2*nBinsTrue;

	TF1*func = new TF1("func",recoDistribution,_xmin,_xmax,6);
	func->SetParameters(_distNorm,_power,_peakNormRel,_distMean,_distSigma,_resSigma);	
	TString histName = "hReco";
	TH1F*hist = new TH1F(histName,"",nBinsReco,_xmin,_xmax);
	hist->SetLineColor(kBlack);
	hist->SetMarkerStyle(20);
	hist->SetMarkerColor(kBlack);

	double nEntries;
	double binMin;
	double binMax;
	double binCenter;
	double xMin = _xmin-7*_resSigma;
	double xMax = _xmax+7*_resSigma;

	for(int i=0;i<=nBinsReco+1;i++){
		if(i==0) binMin = xMin;
		else binMin = hist->GetXaxis()->GetBinLowEdge(i);
		if(i==nBinsReco+1) binMax = xMax;
		else binMax = hist->GetXaxis()->GetBinUpEdge(i);
		nEntries = func->Integral(binMin,binMax);
		binCenter = (binMax+binMin)/2.0;
		hist->SetBinContent(i,nEntries);
		if(nEntries!=0) hist->SetBinError(i,TMath::Sqrt(nEntries));
	}
	return hist;
}

TH2F*ToyModel::GetMigrationMatrix()
{
	int nBinsTrue = _nBins;
	int nBinsReco = 2*nBinsTrue;

	double xMin = _xmin-7*_resSigma;
	double xMax = _xmax+7*_resSigma;
	TF2*integrand2DFunc = new TF2("integrand2DFunc",integrand2D,xMin,xMax,xMin,xMax,6);
	integrand2DFunc->SetParameters(_distNorm,_power,_peakNormRel,_distMean,_distSigma,
				       _resSigma);

	TString histName = "hMatrix";
	TH2F*migrationHist = new TH2F(histName, "",nBinsReco,_xmin,_xmax,
			              nBinsTrue,_xmin,_xmax);
	double xlow,xhi,ylow,yhi,yield;
	double xBinCenter,yBinCenter;
        for(int iBin = 0; iBin <= nBinsReco+1; iBin++){//loop over columns
                for(int jBin = 0; jBin <= nBinsTrue+1; jBin++){//looop over rows
                        if(iBin==0) xlow = xMin;
                        else xlow = migrationHist->GetXaxis()->GetBinLowEdge(iBin);

                        if(iBin==nBinsReco+1) xhi = xMax;
                        else xhi  = migrationHist->GetXaxis()->GetBinUpEdge (iBin);

                        if(jBin==0) ylow = xMin;
                        else ylow = migrationHist->GetYaxis()->GetBinLowEdge(jBin);

                        if(jBin==nBinsTrue+1) yhi = xMax;
                        else  yhi  = migrationHist->GetYaxis()->GetBinUpEdge (jBin);

			xBinCenter = (xhi+xlow)/2.0;
			yBinCenter = (yhi+ylow)/2.0;
                        double nEntries = integrand2DFunc->Integral(xlow,xhi,ylow,yhi);
			migrationHist->SetBinContent(iBin,jBin,nEntries);
                }//end loop over rows
        }//end loop over columns
	return migrationHist;
}

TH2F*ToyModel::GetResponseMatrix(TH2F*hist)
{
	int nBinsTrue = hist->GetNbinsY(); //hard coded such tat true is on y-axis
	int nBinsReco = hist->GetNbinsX(); //and reco is on x-axis	
	TH2F*hResponse = (TH2F*)hist->Clone("hResponse");
	double sum;
	double binContent;
	double newContent;

	for(int j=1;j<=nBinsTrue;j++){
		sum = 0.0;
		for(int i=1;i<=nBinsReco;i++){
			sum += hist->GetBinContent(i,j);
		}//end i loop
		if(sum!=0){
			for(int i=1;i<=nBinsReco;i++){
				binContent = hist->GetBinContent(i,j);
				newContent = binContent/sum;
				hResponse->SetBinContent(i,j,newContent); 
			}//end i loop
		}//end if sum!=0
		else{
			cout << "Migration matrix not normalizable" << endl;
			return hist;
		}
	}//end j loop
	return hResponse;
}//end function
