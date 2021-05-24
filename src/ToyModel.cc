#include "../include/ToyModel.hh"
using namespace Utilities;
using namespace GlobalVariables;

ToyModel::ToyModel()
{
	TH1::SetDefaultSumw2();
	SetModelFunctions();
};

ToyModel::ToyModel(double distNorm,double shift,double peakNormRel,double distMean,
		   double distSigma,double resSigma,vector<double>binningTrue,
		   vector<double>binningReco)
{
	SetModelParameters(distNorm,shift,peakNormRel,distMean,distSigma,resSigma,binningTrue,
			   binningReco);
	SetModelFunctions();
}


void ToyModel::SetModelParameters(double distNorm,double shift,double peakNormRel,
				  double distMean,double distSigma,double resSigma,
				  vector<double>binningTrue,vector<double>binningReco)
{
	int lastBinTrue = binningTrue.size();
	int lastBinReco = binningReco.size();
	int nBinsTrue = lastBinTrue - 1;
	int nBinsReco = lastBinReco - 1;
	_distNorm = distNorm;
	_shift = shift;
	_peakNormRel = peakNormRel;
	_distMean = distMean;
	_distSigma = distSigma;
	_resSigma = resSigma;

	_xmin = binningTrue.at(0);//lower bound of first bin
	_xmax = binningReco.at(nBinsReco);//upper bound of last bin
	cout << "Distribution mass range: " << _xmin << " to " << _xmax << " GeV" << endl; 
	_xMin = _xmin-_nSigma*_resSigma;//lower bound for integration
	_xMax = _xmax+_nSigma*_resSigma;//upper bound for ingetration
	cout << "Functions defined on range: " << _xMin << " to " << _xMax << " GeV" << endl;
	_nBinsTrue = nBinsTrue;//number of bins in true histogram
	_nBinsReco = nBinsReco;//number of bins in reco histogram
	_binningTrue = binningTrue;//vector of bin boundaries for true histogram
	_binningReco = binningReco;//vector of bin boundaries for reco histogram
}

void ToyModel::SetModelFunctions()
{
	TF1*trueFunc = new TF1("trueFunc",trueDistribution,_xMin,_xMax,5);
	trueFunc->SetParameters(_distNorm,_shift,_peakNormRel,_distMean,_distSigma);	

	TF1*recoFunc = new TF1("recFunc",recoDistribution,_xMin,_xMax,6);
	recoFunc->SetParameters(_distNorm,_shift,_peakNormRel,_distMean,_distSigma,_resSigma);

	TF2*integrand2DFunc = new TF2("integrand2DFunc",integrand2D,_xMin,_xMax,_xMin,_xMax,6);
	integrand2DFunc->SetParameters(_distNorm,_shift,_peakNormRel,_distMean,_distSigma,
				       _resSigma);

	_trueFunc = trueFunc;
	_recoFunc = recoFunc;
	_matrixFunc = integrand2DFunc;
}

TF1*ToyModel::GetTrueFunction()
{
	return _trueFunc;
}

TF1*ToyModel::GetRecoFunction()
{
	return _recoFunc;
}

int ToyModel::GetNSigma()
{
	return _nSigma;
}

TH1F*ToyModel::GetTrueHist(TString trueName)
{
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

	cout << "***** Getting true histogram *****" << endl;
	for(int i=0;i<=nBinsTrue+1;i++){
		if(i==0) binMin = _xMin;
		else binMin = hist->GetXaxis()->GetBinLowEdge(i);
		if(i==nBinsTrue+1) binMax = _xMax;
		else binMax = hist->GetXaxis()->GetBinUpEdge(i);
		nEntries = _trueFunc->Integral(binMin,binMax);
		hist->SetBinContent(i,nEntries);
		hist->SetBinError(i,TMath::Sqrt(nEntries));
		cout << "bin " << i << ": ["<< binMin << ", " << binMax << "]" << endl;
	}

	return hist;
}

TH1F*ToyModel::GetTrueHistRandom(TString trueName,Long64_t nEntries)
{
	int nBins = _nBinsTrue;
	double binning[nBins+1];
	for(int i=0;i<=nBins;i++){
		binning[i] = _binningTrue.at(i);
	}

	TH1F*hist = new TH1F(trueName,"",nBins,binning);
	hist->SetFillColor(kRed+2);
	hist->SetLineColor(kRed+2);

	double content;
	for(Long64_t i=0;i<nEntries;i++){
		content = _trueFunc->GetRandom();
		hist->Fill(content);
	}
	return hist;
}

TH1F*ToyModel::GetRecoHist(TString recoName)
{
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

	cout << "***** Getting reco histogram *****" << endl;
	for(int i=0;i<=nBinsReco+1;i++){
		if(i==0) binMin = _xMin;
		else binMin = hist->GetXaxis()->GetBinLowEdge(i);
		if(i==nBinsReco+1) binMax = _xMax;
		else binMax = hist->GetXaxis()->GetBinUpEdge(i);
		nEntries = _recoFunc->Integral(binMin,binMax);
		hist->SetBinContent(i,nEntries);
		hist->SetBinError(i,TMath::Sqrt(nEntries));
		cout << "bin " << i << ": ["<< binMin << ", " << binMax << "]" << endl;
	}
	return hist;
}

TH1F*ToyModel::GetRecoHistRandom(TString recoName,Long64_t nEntries)
{
	int nBins = _nBinsReco;
	double binning[nBins+1];
	for(int i=0;i<=nBins;i++){
		binning[i] = _binningReco.at(i);
	}

	TH1F*hist = new TH1F(recoName,"",nBins,binning);
	hist->SetLineColor(kBlack);
	hist->SetMarkerStyle(20);
	hist->SetMarkerColor(kBlack);

	double content;
	for(Long64_t i=0;i<nEntries;i++){
		content = _recoFunc->GetRandom();
		hist->Fill(content);
	}
	return hist;
}

TH2F*ToyModel::GetMigrationMatrix(TString matrixName)
{
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
	cout << "***** Getting migration histogram *****" << endl;
        for(int i=0;i<=nBinsReco+1;i++){//loop over columns
                for(int j=0;j<=nBinsTrue+1;j++){//looop over rows
                        if(i==0) xlow = _xMin;
                        else xlow = migrationHist->GetXaxis()->GetBinLowEdge(i);

                        if(i==nBinsReco+1) xhi = _xMax;
                        else xhi  = migrationHist->GetXaxis()->GetBinUpEdge(i);

                        if(j==0) ylow = _xMin;
                        else ylow = migrationHist->GetYaxis()->GetBinLowEdge(j);

                        if(j==nBinsTrue+1) yhi = _xMax;
                        else  yhi  = migrationHist->GetYaxis()->GetBinUpEdge(j);

                        double nEntries = _matrixFunc->Integral(xlow,xhi,ylow,yhi);
			migrationHist->SetBinContent(i,j,nEntries);
		if(i==0) cout << "bin " << j << ": ["<< ylow << ", " << yhi << "]" << endl;
                }//end loop over rows
		cout << "bin " << i << ": ["<< xlow << ", " << xhi << "]" << endl;
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
