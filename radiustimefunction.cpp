#include "radiustimefunction.h"

// pass in the time spectrum as a vector
TF1 *radius_for_time(vector<float> times, float t0)
{
	// make shifted histogram
	TH1F *shiftgraph;
	shiftgraph = new TH1F("", "...", 230, 0, 3600);

	for (int k = 0; k < times.size(); k++)
		shiftgraph->Fill(times.at(k) - t0);

	// integrate shifted histogram
	TH1 *h_cumulativeshift;
	h_cumulativeshift = new TH1F("", "...", 230, 0, 3600);

	shiftgraph->GetXaxis()->SetRange(0, shiftgraph->GetNbinsX() + 1);
	h_cumulativeshift = shiftgraph->GetCumulative();

	// rescale integrated histogram
	double inte = h_cumulativeshift->GetMaximum();
	h_cumulativeshift->Scale((1.46 / inte));
	h_cumulativeshift->Sumw2(0);

	// interpolate
	TF1 *f;
	f = new TF1(
		"Radius as a function of time",
		[=](double *x, double * /*p*/)
		{ return h_cumulativeshift->Interpolate(x[0]); },
		h_cumulativeshift->GetXaxis()->GetXmin(), h_cumulativeshift->GetXaxis()->GetXmax(), 0);

	return f;
}

// overload for making graph
TF1 *radius_for_time(vector<float> times, float t0, int chambernum)
{
	// make shifted histogram
	TH1F *shiftgraph;
	shiftgraph = new TH1F("", "...", 230, 0, 3600);

	for (int k = 0; k < times.size(); k++)
		shiftgraph->Fill(times.at(k) - t0);

	// integrate shifted histogram
	TH1 *h_cumulativeshift;
	h_cumulativeshift = new TH1F("", "...", 230, 0, 3600);

	shiftgraph->GetXaxis()->SetRange(0, shiftgraph->GetNbinsX() + 1);
	h_cumulativeshift = shiftgraph->GetCumulative();

	// rescale integrated histogram
	double inte = h_cumulativeshift->GetMaximum();
	h_cumulativeshift->Scale((1.46 / inte));
	h_cumulativeshift->Sumw2(0);

	// interpolate
	TF1 *f;
	f = new TF1(
		"f",
		[=](double *x, double * /*p*/)
		{ return h_cumulativeshift->Interpolate(x[0]); },
		h_cumulativeshift->GetXaxis()->GetXmin(), h_cumulativeshift->GetXaxis()->GetXmax(), 0);

	TCanvas *c1 = new TCanvas();

	f->Draw("");
	c1->Modified();
	c1->Update();

	h_cumulativeshift->Draw("same");

	c1->Modified();
	c1->Update();

	string title = "r(t) for Chamber" + std::to_string(chambernum) + ".png";
	string titleFullPath = "./Output/" + title;
	string grtitle = "r(t) for Chamber" + std::to_string(chambernum);

	c1->SetTitle(grtitle.c_str());
	f->GetYaxis()->SetTitle("Distance (cm)");
	f->GetXaxis()->SetTitle("Time (ns)");

	c1->Modified();
	c1->Update();
	// c1->cd(0);
	c1->SaveAs(titleFullPath.c_str());

	delete c1;
	//////////////////

	TCanvas *c2 = new TCanvas();

	string grtitle1 = "Drift Velocity for Time of Chamber " + std::to_string(chambernum);
	c2->SetTitle(grtitle1.c_str());

	TGraph *derivf = (TGraph *)f->DrawDerivative();

	f->SetTitle(grtitle1.c_str());
	derivf->SetTitle(grtitle1.c_str());
	derivf->GetYaxis()->SetTitle("Velocity (cm/ns)");
	derivf->GetXaxis()->SetTitle("Time (ns)");
	derivf->Draw("");
	c2->Modified();
	c2->Update();

	// h_cumulativeshift->Draw("same");

	c2->Modified();
	c2->Update();

	string title1 = "Drift Velocity for Time of Chamber " + std::to_string(chambernum) + ".png";
	string title1FullPath = "./Output/" + title1;

	c2->Modified();
	c2->Update();
	// c1->cd(0);
	c2->SaveAs(title1FullPath.c_str());

	return f;
}