# Creating a toy model
Toy models are creating using the ToyModel class, defined like this:  
```
ToyModel*model = new ToyModel(distNorm, power, peakNormRel, distMean, distSigma, resSigma, binningTrue, binningReco);
``` 

### Variable definitions:  
**distNorm:** A normalization factor for the whole distribution  
**power:** The power used for the power law  
**peakNormRel:** A normalization just for the peak. This allows one to define the peak height compared to the power law  
**distMean:** This is the mean of the peak of the distribution  
**distSigma:** This is the width of the peak of the distribution  
**resSigma:** The is the width of the smearing function used to create the reco distribution from the true distribution  
**binningTrue:** A c++ vector containing the bin boundaries of the true distribution  
**binningReco:** A c++ vector containing the bin boundaries of the reco distribution  

To carry out unfolding, you need a true histogram, a reco (or observed) histogram, and a matrix of migrations. To get each of these, use the following functions:  
```
TH1F*histTrue = model->GetTrueHist("trueName");   
TH1F*histReco = model->GetRecoHist("recoName");  
TH2F*histMatrix = model->GetMigrationMatrix("matrixName");
```  

Then if you want to see the normalized matrix of migrations, the response matrix, use this:  
```
model->GetResponseMatrix(histMatrix);
```
The response matrix is not needed for unfolding as it will be calculated automatically by the unfolding functions, but it is useful to have.  

# Unfolding
This package has two methods of unfolding: TUnfold(https://root.cern.ch/doc/master/classTUnfold.html) and a simple matrix inversion method.  
First, create an Unfold class object:  
```
Unfold*unfold = new Unfold();
```
You need to define a regularization type. At the moment, three types are implemented: NO_REG (no regularization), CONST_REG (user defined value of the regularization parameter. this is currently hard-coded because I never use it, but it should be changed), and VAR_REG_LCURVE (used the L-curve scan method to find the best value of the regularization parameter). 
