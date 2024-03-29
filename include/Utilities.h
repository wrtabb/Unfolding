#ifndef Utilities_H
#define Utilities_H
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>

namespace Utilities{

	TH1F*RebinTH1(TH1F*hist,TString histName,TH1F*hBinning)
	{
		int nBinsOld = hist->GetNbinsX();
		int nBinsNew = hBinning->GetNbinsX();
		if(nBinsNew > nBinsOld){
			cout << "*********************************************************" << endl;
			cout << "ERROR: new binning must have fewer bins than old binning!" << endl;
			cout << "*********************************************************" << endl;
			return hist;
		}
		else if(nBinsNew == nBinsOld){
			cout << "************************************************" << endl;
			cout << "Histogram already has the same binning as source" << endl;
			cout << "Returning histogram" << endl;
			cout << "************************************************" << endl;
			return hist;
		}
		double newbinning[nBinsNew];
		for(int i=0;i<=nBinsNew;i++){
			if(i==0) newbinning[i] = hBinning->GetBinLowEdge(i+1);
			else newbinning[i] = newbinning[i-1]+hBinning->GetBinWidth(i);
		}
		TH1F*hRebin = new TH1F(histName,"",nBinsNew,newbinning);
		double x;
		double nEntries;
		for(int j=0;j<=nBinsOld+1;j++){
			x = hist->GetXaxis()->GetBinCenter(j);
			nEntries = hist->GetBinContent(j);
			hRebin->Fill(x,nEntries);
		} //end x bin loop

		double histBinWidth;
		double newBinWidth;
		double nBinsOldInNew[nBinsNew];
		double newBinUpperEdge,oldBinUpperEdge;
		int k = 1;
		for(int i=0;i<=nBinsNew+1;i++){
			newBinUpperEdge = hRebin->GetBinLowEdge(i)+hRebin->GetBinWidth(i);
			double newError2 = 0;
			for(int j=k;j<=nBinsOld;j++){
				oldBinUpperEdge = hist->GetBinLowEdge(j)+hist->GetBinWidth(j);
				if(oldBinUpperEdge > newBinUpperEdge) continue;
				k = j;
				newError2 += hist->GetBinError(j)*hist->GetBinError(j);
			}// end loop over j bins
			hRebin->SetBinError(i,sqrt(newError2));
		}// end loop over i bins
		return hRebin;

	}// end Rebin1D

	TH2F*RebinTH2X(TH2F*hist,TString histName,TH1F*hBinning)
	{
		int nBinsOld = hist->GetNbinsX();
		int nBinsNew = hBinning->GetNbinsX();
		if(nBinsNew > nBinsOld){
			cout << "*********************************************************" << endl;
			cout << "ERROR: new binning must have fewer bins than old binning!" << endl;
			cout << "*********************************************************" << endl;
			return hist;
		}
		else if(nBinsNew == nBinsOld){
			cout << "************************************************" << endl;
			cout << "Histogram already has the same binning as source" << endl;
			cout << "Returning histogram" << endl;
			cout << "************************************************" << endl;
			return hist;
		}

		double newbinning[nBinsNew];
		for(int i=0;i<=nBinsNew;i++){
			if(i==0) newbinning[i] = hBinning->GetBinLowEdge(i+1);
			else newbinning[i] = newbinning[i-1]+hBinning->GetBinWidth(i);
		}
		TH2F*hRebin = new TH2F(histName,"",nBinsNew,newbinning,nBinsNew,newbinning);
		double y,x;
		double nEntries;
		for(int i=0;i<=nBinsOld+1;i++){//x loop
			x = hist->GetXaxis()->GetBinCenter(i);
			for(int j=0;j<=nBinsNew+1;j++){//y loop
				y = hist->GetYaxis()->GetBinCenter(j);
				nEntries = hist->GetBinContent(i,j);
				hRebin->Fill(x,y,nEntries);
			}//end y bin loop
		} //end x bin loop

		double histBinWidth;
		double newBinWidth;
		double nBinsOldInNew[nBinsNew];
		double newBinUpperEdge,oldBinUpperEdge;
		int k = 1;
		for(int i=1;i<=nBinsNew;i++){
			newBinUpperEdge = hBinning->GetBinLowEdge(i)+hBinning->GetBinWidth(i);
			double newError2 = 0;
			for(int j=k;j<=nBinsOld;j++){
				oldBinUpperEdge = hist->GetXaxis()->GetBinLowEdge(j) +
						  hist->GetXaxis()->GetBinWidth(j);
				if(oldBinUpperEdge > newBinUpperEdge) continue;
				k = j;
				newError2 += hist->GetBinError(i,j)*hist->GetBinError(i,j);
				hRebin->SetBinError(i,j,sqrt(newError2));
			}// end j loop
		}// end i loop
		return hRebin;

	}// end Rebin2DX

	TH2F*RebinTH2Y(TH2F*hist,TString histName,TH1F*hBinning)
	{
		int nBinsOld = hist->GetNbinsY();
		int nBinsNew = hBinning->GetNbinsX();
		if(nBinsNew > nBinsOld){
			cout << "*********************************************************" << endl;
			cout << "ERROR: new binning must have fewer bins than old binning!" << endl;
			cout << "*********************************************************" << endl;
			return hist;
		}
		else if(nBinsNew == nBinsOld){
			cout << "************************************************" << endl;
			cout << "Histogram already has the same binning as source" << endl;
			cout << "Returning histogram" << endl;
			cout << "************************************************" << endl;
			return hist;
		}

		double newbinning[nBinsNew];
		for(int i=0;i<=nBinsNew;i++){
			if(i==0) newbinning[i] = hBinning->GetBinLowEdge(i+1);
			else newbinning[i] = newbinning[i-1]+hBinning->GetBinWidth(i);
		}
		TH2F*hRebin = new TH2F(histName,"",nBinsNew,newbinning,nBinsNew,newbinning);
		double y,x;
		double nEntries;
		for(int i=0;i<=nBinsOld+1;i++){//y loop
			y = hist->GetYaxis()->GetBinCenter(i);
			for(int j=0;j<=nBinsNew+1;j++){//x loop
				x = hist->GetXaxis()->GetBinCenter(j);
				nEntries = hist->GetBinContent(j,i);
				hRebin->Fill(x,y,nEntries);
			}//end y bin loop
		} //end x bin loop

		double histBinWidth;
		double newBinWidth;
		double nBinsOldInNew[nBinsNew];
		double newBinUpperEdge,oldBinUpperEdge;
		int k = 1;
		for(int i=1;i<=nBinsNew;i++){
			newBinUpperEdge = hBinning->GetBinLowEdge(i)+hBinning->GetBinWidth(i);
			double newError2 = 0;
			for(int j=k;j<=nBinsOld;j++){
				oldBinUpperEdge = hist->GetXaxis()->GetBinLowEdge(j) +
						  hist->GetXaxis()->GetBinWidth(j);
				if(oldBinUpperEdge > newBinUpperEdge) continue;
				k = j;
				newError2 += hist->GetBinError(j,i)*hist->GetBinError(j,i);
				hRebin->SetBinError(j,i,sqrt(newError2));
			}// end j loop
		}// end i loop
		return hRebin;

	}// end Rebin2DY

    TH2F*RebinTH2(TH2F*hist,TString histName,TH1F*hBinning,bool trueVert)
    {
        using namespace Utilities;
        if(trueVert) return RebinTH2X(hist,histName,hBinning);
        else return RebinTH2Y(hist,histName,hBinning);
    }

}//end namespace

#endif
