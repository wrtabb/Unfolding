
void Counter(Long64_t count,Long64_t totalCount);

void makeToyModelsRandomSeeds(int nSeeds)
{
	TH1::SetDefaultSumw2();
	TFile*file = new TFile("data/toyModel.root");
	TH1F*_hTrue = (TH1F*)file->Get("hTrue");
	TH1F*_hReco = (TH1F*)file->Get("hReco");

	vector<TH1F*> hTrue;
	vector<TH1F*> hReco;
	
	Long64_t nEvents = 1e7;
	double trueContent;
	double recoContent;

	int nBinsTrue = _hTrue->GetNbinsX();
	int nBinsReco = 2*nBinsTrue;
	double xmin = _hTrue->GetBinLowEdge(1);
	double xmax = _hTrue->GetBinLowEdge(nBinsTrue)+_hTrue->GetBinWidth(nBinsTrue);;
	int output = 0;
	Long64_t nTotal = nSeeds*nEvents;
	Long64_t count = 0;

	for(int i=0;i<nSeeds;i++){
		TString trueName = "histTrue";
		trueName += i;
		TString recoName = "histReco";
		recoName += i;

		hTrue.push_back(new TH1F(trueName,"",nBinsTrue,xmin,xmax));
		hReco.push_back (new TH1F(recoName,"",nBinsReco,xmin,xmax));
		gRandom->SetSeed(i);
		for(Long64_t j=0;j<nEvents;j++){
			trueContent = _hTrue->GetRandom();
			recoContent = _hReco->GetRandom();
			hTrue.at(i)->Fill(trueContent);
			hReco.at(i)->Fill(recoContent);

			if(count%(nTotal/100)==0){
				output += 1;
				cout << output << "%" << endl;
			}
			count++;
		}	
	}
	TFile*saveFile = new TFile("data/toyModelRandomSeeds.root","recreate");
	for(int i=0;i<nSeeds;i++){
		hTrue.at(i)->Write();
		hReco.at(i)->Write();
	}
	saveFile->Close();
}

void Counter(Long64_t count,Long64_t totalCount)
{
	//This function does not work for some reason
	//It prints 0% every time
	double output;
	output = 100.0*(count/totalCount);
	if(count%(totalCount/100)==0) cout << output << "%" << endl;
	else return;
}
