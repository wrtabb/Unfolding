#include "../include/ToyModel.hh"

ToyModel::ToyModel(double distNorm,double power,double peakNormRel,double distMean,
		   double distSigma,double resSigma,vector<double>binningTrue,
		   vector<double>binningReco)
{
	TH1::SetDefaultSumw2();
	SetModelParameters(distNorm,power,peakNormRel,distMean,distSigma,resSigma,binningTrue,
			   binningReco);
}

void ToyModel::SetModelFunctions(TF1*trueFunc,TF1*recoFunc,TF2*matrixFunc)
{
	_trueFunc = trueFunc;
	_recoFunc = recoFunc;
	_matrixFunc = matrixFunc;
}

void ToyModel::SetModelParameters(double distNorm,double power,double peakNormRel,
				  double distMean,double distSigma,double resSigma,
				  vector<double>binningTrue,vector<double>binningReco)
{
	int lastBinTrue = binningTrue.size();
	int lastBinReco = binningReco.size();
	int nBinsTrue = lastBinTrue - 1;
	int nBinsReco = lastBinReco - 1;

	_distNorm = distNorm;
	_power = power;
	_peakNormRel = peakNormRel;
	_distMean = distMean;
	_distSigma = distSigma;
	_resSigma = resSigma;
	_xmin = binningTrue.at(0);
	_xmax = binningReco.at(nBinsReco);
	_nBinsTrue = nBinsTrue;
	_nBinsReco = nBinsReco;
	_binningTrue = binningTrue;
	_binningReco = binningReco;
}

TH1F*ToyModel::GetTrueHist(TString trueName)
{
	double xMin = _xmin-7*_resSigma;
	double xMax = _xmax+7*_resSigma;
	TF1*func = new TF1("func",trueDistribution,xMin,xMax,5);
	func->SetParameters(_distNorm,_power,_peakNormRel,_distMean,_distSigma);	

	int nBinsTrue = _nBinsTrue;
	double binning[nBinsTrue+1];
	for(int i=0;i<=nBinsTrue;i++){
		binning[i] = _binningTrue.at(i);
	}

	TH1F*hist = new TH1F(trueName,"",nBinsTrue,binning);
	hist->SetFillColor(kRed+2);
	hist->SetLineColor(kRed+2);

	double nEntries;
	double binMin;
	double binMax;
	double binCenter;

	for(int i=0;i<=nBinsTrue+1;i++){
		if(i==0) binMin = xMin;
		else binMin = hist->GetXaxis()->GetBinLowEdge(i);
		if(i==nBinsTrue+1) binMax = xMax;
		else binMax = hist->GetXaxis()->GetBinUpEdge(i);
		binCenter = (binMax+binMin)/2.0;
		nEntries = func->Integral(binMin,binMax);
		hist->SetBinContent(i,nEntries);
		hist->SetBinError(i,TMath::Sqrt(nEntries));
	}

	return hist;
}


TH1F*ToyModel::GetRecoHist(TString recoName)
{
	double xMin = _xmin-7*_resSigma;
	double xMax = _xmax+7*_resSigma;
	TF1*func = new TF1("func",recoDistribution,xMin,xMax,6);
	func->SetParameters(_distNorm,_power,_peakNormRel,_distMean,_distSigma,_resSigma);

	int nBinsReco = _nBinsReco;
	double binning[nBinsReco+1];
	for(int i=0;i<=nBinsReco;i++){
		binning[i] = _binningReco.at(i);
	}

	TH1F*hist = new TH1F(recoName,"",nBinsReco,binning);
	hist->SetLineColor(kBlack);
	hist->SetMarkerStyle(20);
	hist->SetMarkerColor(kBlack);

	double nEntries;
	double binMin;
	double binMax;
	double binCenter;

	for(int i=0;i<=nBinsReco+1;i++){
		if(i==0) binMin = xMin;
		else binMin = hist->GetXaxis()->GetBinLowEdge(i);
		if(i==nBinsReco+1) binMax = xMax;
		else binMax = hist->GetXaxis()->GetBinUpEdge(i);
		nEntries = func->Integral(binMin,binMax);
		binCenter = (binMax+binMin)/2.0;
		hist->SetBinContent(i,nEntries);
		hist->SetBinError(i,TMath::Sqrt(nEntries));
	}
	return hist;
}

TH2F*ToyModel::GetMigrationMatrix(TString matrixName)
{
	double xMin = _xmin-7*_resSigma;
	double xMax = _xmax+7*_resSigma;
	TF2*integrand2DFunc = new TF2("integrand2DFunc",integrand2D,xMin,xMax,xMin,xMax,6);
	integrand2DFunc->SetParameters(_distNorm,_power,_peakNormRel,_distMean,_distSigma,
				       _resSigma);
	int nBinsTrue = _nBinsTrue;
	int nBinsReco = _nBinsReco;
	double binningTrue[nBinsTrue+1];
	double binningReco[nBinsReco+1];

	for(int i=0;i<=nBinsTrue;i++){
		binningTrue[i] = _binningTrue.at(i);
	}

	for(int i=0;i<=nBinsReco;i++){
		binningReco[i] = _binningReco.at(i);
	}

	TH2F*migrationHist = new TH2F(matrixName, "",nBinsReco,binningReco,nBinsTrue,
				      binningTrue);
	double xlow,xhi,ylow,yhi,yield;
        for(int i=0;i<=nBinsReco+1;i++){//loop over columns
                for(int j=0;j<=nBinsTrue+1;j++){//looop over rows
                        if(i==0) xlow = xMin;
                        else xlow = migrationHist->GetXaxis()->GetBinLowEdge(i);

                        if(i==nBinsReco+1) xhi = xMax;
                        else xhi  = migrationHist->GetXaxis()->GetBinUpEdge (i);

                        if(j==0) ylow = xMin;
                        else ylow = migrationHist->GetYaxis()->GetBinLowEdge(j);

                        if(j==nBinsTrue+1) yhi = xMax;
                        else  yhi  = migrationHist->GetYaxis()->GetBinUpEdge (j);

                        double nEntries = integrand2DFunc->Integral(xlow,xhi,ylow,yhi);
			migrationHist->SetBinContent(i,j,nEntries);
                }//end loop over rows
        }//end loop over columns
	makeResponseMatrix(migrationHist);
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

	for(int j=0;j<=nBinsTrue+1;j++){
		sum = 0.0;
		for(int i=0;i<=nBinsReco+1;i++){
			sum += hist->GetBinContent(i,j);
		}//end i loop
		if(sum!=0){
			for(int i=0;i<=nBinsReco+1;i++){
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

TH2F*ToyModel::GetResponseMatrixT(TH2F*hist)
{
	int nBinsTrue = hist->GetNbinsY(); //hard coded such tat true is on y-axis
	int nBinsReco = hist->GetNbinsX(); //and reco is on x-axis	
	TH2F*hResponseT = (TH2F*)hist->Clone("hResponseT");
	double sum;
	double binContent;
	double newContent;

	for(int i=0;i<=nBinsReco+1;i++){
		sum = 0.0;
		for(int j=0;j<=nBinsTrue+1;j++){
			sum += hist->GetBinContent(i,j);
		}//end i loop
		if(sum!=0){
			for(int j=0;j<=nBinsTrue+1;j++){
				binContent = hist->GetBinContent(i,j);
				newContent = binContent/sum;
				hResponseT->SetBinContent(i,j,newContent); 
			}//end i loop
		}//end if sum!=0
		else{
			cout << "Migration matrix not normalizable" << endl;
			return hist;
		}
	}//end j loop
	return hResponseT;
}//end function
