
{
	gROOT->SetBatch("true");
	gROOT->ProcessLine(".L src/ToyModel.cc+");
	gROOT->ProcessLine(".L src/Unfolding.cc+");
}
