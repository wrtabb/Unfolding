
#include "include/GlobalVariables.h"
using namespace GlobalVariables;

void createBinning()
{
	float lowEdge;
	float highEdge;
	float newBin;
	float nBinsNew = 4.0;
	for(int i=0;i<_massbinningTrue.size();i++){
		lowEdge = _massbinningTrue.at(i);
		highEdge = _massbinningTrue.at(i+1);	
		newBin = (highEdge-lowEdge)/nBinsNew;

		if(i==0) cout << lowEdge << ", " << endl;
		for(int j=1;j<nBinsNew;j++){
			cout << lowEdge+newBin << ", " << endl;
			newBin = newBin+newBin;
		}
		cout << highEdge << ", " << endl;
	}
}
