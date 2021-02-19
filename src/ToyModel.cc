#include "../include/ToyModel.hh"

ToyModel::ToyModel(double distNorm,double expDecay,double peakNormRel,double distMean,
		   double distSigma,double resSigma,double xmin,double xmax,int nBins)
{
	SetModelParameters(distNorm,expDecay,peakNormRel,distMean,distSigma,resSigma,xmin,xmax,
			   nBins);

	TF1*trueFunc = new TF1("trueFunction",trueDistribution,0,100,5);	
	trueFunc->SetParameters(distNorm,expDecay,peakNormRel,distMean,distSigma);
	
	TF1*recoFunc = new TF1("recoFunction",recoDistribution,-10,150,6);
	recoFunc->SetParameters(distNorm,expDecay,peakNormRel,distMean,distSigma,resSigma);

	TF2*matrixFunc = new TF2("integrand2DFunc", integrand2D,0,100,-10,150,6);
	matrixFunc->SetParameters(distNorm,expDecay,peakNormRel,distMean,distSigma,resSigma);

	SetModelFunctions(trueFunc,recoFunc,matrixFunc);
}

void ToyModel::SetModelFunctions(TF1*trueFunc,TF1*recoFunc,TF2*matrixFunc)
{
	_trueFunc = trueFunc;
	_recoFunc = recoFunc;
	_matrixFunc = matrixFunc;
}

void ToyModel::SetModelParameters(double distNorm,double expDecay,double peakNormRel,
				  double distMean,double distSigma,double resSigma,double xmin, 				 double xmax,int nBins)
{
	_distNorm = distNorm;
	_expDecay = expDecay;
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
	func->SetParameters(_distNorm,_expDecay,_peakNormRel,_distMean,_distSigma,_resSigma);	
	TH1F*hist = new TH1F("hTrue","",_nBins,0,100);
	hist->SetFillColor(kRed+2);
	hist->SetLineColor(kRed+2);

	double content;
	double binMin;
	double binMax;
	for(int i=1;i<=_nBins;i++){
		binMin = hist->GetXaxis()->GetBinLowEdge(i);
		binMax = hist->GetXaxis()->GetBinUpEdge(i);
		content = func->Integral(binMin,binMax);
		hist->SetBinContent(i,content);
		hist->SetBinError(i,TMath::Sqrt(content));
	}

	return hist;
}


TH1F*ToyModel::GetRecoHist(bool useTUnfold)
{
	int nBins;
	if(!useTUnfold) nBins = _nBins; 
	else nBins = 2*_nBins;//TUnfold requires more reco bins than true
	TF1*func = new TF1("func",recoDistribution,_xmin,_xmax,6);
	func->SetParameters(_distNorm,_expDecay,_peakNormRel,_distMean,_distSigma,_resSigma);	
	TString histName;
	if(useTUnfold) histName = "hRecoTUnfold";
	else histName = "hRecoInversion";
	TH1F*hist = new TH1F(histName,"",nBins,_xmin,_xmax);
	hist->SetLineColor(kBlack);
	hist->SetMarkerStyle(20);
	hist->SetMarkerColor(kBlack);

	double content;
	double binMin;
	double binMax;
	for(int i=1;i<=nBins;i++){
		binMin = hist->GetXaxis()->GetBinLowEdge(i);
		binMax = hist->GetXaxis()->GetBinUpEdge(i);
		content = func->Integral(binMin,binMax);
		hist->SetBinContent(i,content);
		hist->SetBinError(i,TMath::Sqrt(content));
	}
	return hist;
}

TH2F*ToyModel::GetMigrationMatrix(bool useTUnfold)
{
	int nBinsTrue = _nBins;
	int nBinsReco;
	if(useTUnfold) nBinsReco = 2*nBinsTrue;
	else nBinsReco = nBinsTrue;
	double xMin = _xmin-7*_resSigma;
	double xMax = _xmax+7*_resSigma;
	TF2*integrand2DFunc = new TF2("integrand2DFunc",integrand2D,xMin,xMax,xMin,xMax,6);
	integrand2DFunc->SetParameters(_distNorm,_expDecay,_peakNormRel,_distMean,_distSigma,
				       _resSigma);

	TString histName;
	if(useTUnfold) histName = "hMatrixTUnfold";
	else histName = "hMatrixInversion";
	TH2F*migrationHist = new TH2F(histName, "",nBinsReco,_xmin,_xmax,
			              nBinsTrue,_xmin,_xmax);
	double xlow,xhi,ylow,yhi,yield;
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

                        double content = integrand2DFunc->Integral(xlow, xhi, ylow, yhi);
                        migrationHist->SetBinContent(iBin,jBin,content);
			migrationHist->SetBinError(iBin,jBin,TMath::Sqrt(content));
                }//end loop over rows
        }//end loop over columns
	return migrationHist;
}
