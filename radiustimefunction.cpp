#include "radiustimefunction.h"

//pass in the time spectrum as a vector 
TF1* radius_for_time(vector<float> times, float t0)
{
//make shifted histogram
TH1F* shiftgraph;
shiftgraph = new TH1F("", "...", 230, 0, 3600);
	
for( int k = 0 ; k < times.size() ; k++)
	shiftgraph->Fill( times.at(k) - t0);

//integrate shifted histogram
TH1* h_cumulativeshift;
h_cumulativeshift =  new TH1F("", "...", 230, 0, 3600); 

shiftgraph->GetXaxis()->SetRange(0, shiftgraph->GetNbinsX()+1);
h_cumulativeshift = shiftgraph->GetCumulative();

//rescale integrated histogram
double inte = h_cumulativeshift->GetMaximum();
h_cumulativeshift->Scale( ( 1.46 / inte ) );
h_cumulativeshift->Sumw2(0);

//interpolate 
TF1* f ;
f = new TF1("f",
[=](double *x, double */*p*/){return h_cumulativeshift->Interpolate(x[0]);},
	h_cumulativeshift->GetXaxis()->GetXmin(), h_cumulativeshift->GetXaxis()->GetXmax(), 0);
	
return f;
}