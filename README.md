# Creating a toy model
Toy models are creating using the ToyModel class, defined like this:  
	`ToyModel*model = new ToyModel(distNorm, power, peakNormRel, distMean, distSigma, resSigma, binningTrue, binningReco)`  

Variable definitions:  
distNorm: A normalization factor for the whole distribution  
power: The power used for the power law  
peakNormRel: A normalization just for the peak. This allows one to define the peak height compared to the power law  
distMean: This is the mean of the peak of the distribution  
distSigma: This is the width of the peak of the distribution  
resSigma: The is the width of the smearing function used to create the reco distribution from the true distribution  
binningTrue: A c++ vector containing the bin boundaries of the true distribution  
binningReco: A c++ vector containing the bin boundaries of the reco distribution  

To carry out unfolding, you need a true histogram, a reco (or observed) histogram, and a matrix of migrations. To get each of these, use the following functions:  
	`TH1F*histTrue = model->GetTrueHist("trueName");  
	TH1F*histReco = model->GetRecoHist("recoName");  
	TH2F*histMatrix = model->GetMigrationMatrix("matrixName");`  

Then if you want to see the normalized matrix of migrations, the response matrix, use this:  
	`model->GetResponseMatrix(histMatrix);`  
