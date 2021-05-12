#ifndef Utilities_H
#define Utilities_H
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>

namespace Utilities{
	// RebinTH1 works in the macros in the main directory
	// But does not work in the class Unfold in src/Unfolding.cc
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
			}
			hRebin->SetBinError(i,sqrt(newError2));
		}
		return hRebin;

	}

	// RebinTH2 seems to work everywhere I use it
	TH2F*RebinTH2(TH2F*hist,TString histName,TH1F*hBinning)
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
			}
		}
		return hRebin;

	}

	TH2F*makeResponseMatrix(TH2F*hist)
	{
		TH2F*hResponse = (TH2F*)hist->Clone("hResponse");
		int nBinsX = hResponse->GetNbinsX();
		int nBinsY = hResponse->GetNbinsY();
		double nEntriesX;
		double binContent;
		double scaledContent;

		for(int j=0;j<=nBinsY+1;j++){
			nEntriesX = 0.0;
			for(int i=0;i<=nBinsX+1;i++){
				hResponse->SetBinContent(i,j,0);
				binContent = hist->GetBinContent(i,j);
				nEntriesX += binContent;
			}//end first loop over reco bins
			for(int i=0;i<=nBinsX+1;i++){
				scaledContent = hist->GetBinContent(i,j)/nEntriesX;
				if(scaledContent > 1e-4) 
					hResponse->SetBinContent(i,j,scaledContent);
				hResponse->SetBinContent(i,j,scaledContent);
			}//end second loop over reco bins
		}//end loop over true bins

		return hResponse;
	}//end makeResponseMatrix

	TH2F*makeResponseMatrixT(TH2F*hist)
	{
		TH2F*hResponse = (TH2F*)hist->Clone("hResponse");
		int nBinsX = hResponse->GetNbinsX();
		int nBinsY = hResponse->GetNbinsY();
		double nEntriesY;
		double binContent;
		double scaledContent;

		for(int i=0;i<=nBinsX+1;i++){
			nEntriesY = 0.0;
			for(int j=0;j<=nBinsY+1;j++){
				hResponse->SetBinContent(i,j,0);
				binContent = hist->GetBinContent(i,j);
				nEntriesY += binContent;
			}//end first loop over reco bins
			for(int j=0;j<=nBinsY+1;j++){
				scaledContent = hist->GetBinContent(i,j)/nEntriesY;
				if(scaledContent > 1e-4) 
					hResponse->SetBinContent(i,j,scaledContent);
			}//end second loop over reco bins
		}//end loop over true bins

		return hResponse;
	}//end makeResponseMatrix

TCanvas*PlotMatrix(TString canvasName,TString plotTitle,TH2F*hist,bool logPlot)
{
        TCanvas*c1 = new TCanvas(canvasName,"",0,0,1000,1000);
        if(logPlot){
                c1->SetLogx();
                c1->SetLogy();
                hist->GetYaxis()->SetNoExponent();
                hist->GetXaxis()->SetNoExponent();
                hist->GetYaxis()->SetMoreLogLabels();
                hist->GetXaxis()->SetMoreLogLabels();
        }
        c1->SetGrid();
        c1->SetRightMargin(0.13);
        c1->SetLeftMargin(0.13);
        hist->SetTitle(plotTitle);
        hist->GetYaxis()->SetTitle("mass^{true} [GeV]");
        hist->GetXaxis()->SetTitle("mass^{obs} [GeV]");
        hist->Draw("colz");
        return c1;
}

}//end namespace

#endif
