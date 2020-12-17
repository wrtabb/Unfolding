

void integral()
{
 double mean = 20;
 double sigma = 4;
 TF1 func("func","1/(x+1)+gaus(0)",0,50);
 func.SetParameters(1.0,mean,sigma);
 
 ROOT::Math::WrappedTF1 wf1(func);

 //ROOT::Math::GaussIntegrator ig;
 ROOT::Math::GaussLegendreIntegrator ig;
 ig.SetFunction(wf1);
 ig.SetRelTolerance(0.000000001);
 
 cout << ig.Integral(0,50) << endl;
}
