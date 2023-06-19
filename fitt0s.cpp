// given a time distribution h, fit_t0 returns a struct with t0 and t0s error and the fitted graph -- populate t0s will populate everything
// input not a pointer as to not directly change the graph -- will return the fit function (TF1*) at the end in the struct
#include "fitt0s.h"
#include "TCanvas.h"

// pass in the vector of the time spectrum
Fittedt0s fit_t0(vector<float> times, int chambnum)
{

	TH1F *h;
	h = new TH1F("", "...", 230, 0, 3600);

	// TODO: DO WE NEED ERROR BARS ON THESE HISTOGRAMS

	for (int k = 0; k < times.size(); k++)
		h->Fill(times.at(k));

	// set the fitting function
	TF1 *fitfn;
	fitfn = new TF1("fit", "([2] / ( 1 + TMath::Exp( -[0]*(x-[1] ) ) ) )");

	// do the fits
	int binmax = h->GetMaximumBin();
	double t0_guess = h->GetXaxis()->GetBinUpEdge(binmax); // this initial t0 is a guess

	// set guess parameters
	fitfn->SetParameters(1, t0_guess, 4000);

	// here we fit
	//h->SetBit(TH1::kNoStats);
	h->Fit(fitfn, "0", "S", t0_guess - 300., t0_guess + 10);
	//h->ResetBit(TH1::kNoStats);

	double xmin;
	double xmax;

	h->GetFunction("fit")->GetRange(xmin, xmax);

	// plot and save the time

	double xvals[h->FindBin(xmax) - h->FindBin(xmin)];
	double yvals[h->FindBin(xmax) - h->FindBin(xmin)];

	for (int ibin = h->FindBin(xmin); ibin <= h->FindBin(xmax); ++ibin)
	{
		double deez = h->GetBinContent(ibin);

		double xm = h->GetBinCenter(ibin);

		double res = (deez - h->GetFunction("fit")->Eval(h->GetBinCenter(ibin)));

		xvals[ibin - h->FindBin(xmin)] = xm;
		yvals[ibin - h->FindBin(xmin)] = res;
	}

	TGraph *graph = new TGraph(h->FindBin(xmax) - h->FindBin(xmin), xvals, yvals);

	TCanvas *c1 = new TCanvas();
	c1->Divide(1, 2);
	c1->cd(1);

	//h->SetBit(TH1::kNoStats);

	h->Draw("");
	c1->Modified();
	c1->Update();

	h->GetFunction("fit")->Draw("same");
	c1->Modified();
	c1->Update();

	string title = "Fit for Chamber " + std::to_string(chambnum + 1) + ".png";
	string titleFullPath = "./Output/" + title;

	string graphtitle = "t0 Fit for Chamber " + std::to_string(chambnum + 1);

	h->SetTitle(graphtitle.c_str());
	h->GetXaxis()->SetTitle("Time (ns)");
	h->GetYaxis()->SetTitle("Counts");

	c1->Modified();
	c1->Update();

	c1->cd(2);

	graph->SetMarkerStyle(8);
	graph->Draw("");
	c1->Modified();
	c1->Update();

	graph->SetTitle("Fit Residuals");
	graph->GetXaxis()->SetTitle("Time (ns)");
	graph->GetYaxis()->SetTitle("Fit Residuals");

	c1->Modified();
	c1->Update();

	c1->cd(0);

	

	c1->SaveAs(titleFullPath.c_str());

	//h->ResetBit(TH1::kNoStats);

	// get parameters from fits
	double t0;
	double t0err;

	// set up the return thing
	t0 = fitfn->GetParameter(1);
	t0err = fitfn->GetParError(1);

	// return object
	Fittedt0s fittedt0s;

	fittedt0s.t0 = t0;
	fittedt0s.t0err = t0err;
	fittedt0s.fitfn = fitfn;
	fittedt0s.h = h;

	return fittedt0s;
}

vector<vector<float>> get_times(vector<NewEvent> events)
{
	vector<vector<float>> times;
	for (int i = 0; i < 144; i++)
	{
		vector<float> tubetimes;
		times.push_back(tubetimes);
	}

	for (int i = 0; i < events.size(); i++)
	{
		NewEvent currentevent = events.at(i);

		for (int j = 0; j < currentevent.t.size(); j++)
		{
			int tubenum = currentevent.chamber.at(j) * 48 + currentevent.layer.at(j) * 16 + currentevent.tube.at(j);
			times.at(tubenum).push_back(currentevent.t.at(j));
		}
	}
	return times;
}

// pass in the events vector, and names of what you populate the branches with as well as the branches themselves (everything in order)
vector<TF1> populate_t0s(vector<NewEvent> events, double &t0, double &t0err, TF1 &t0fit, TBranch *branch1, TBranch *branch2, TBranch *branch3)
{
	auto times = get_times(events);

	vector<vector<float>> newtimes;
	for (int i = 0; i < 3; ++i)
	{
		vector<float> guh;
		newtimes.push_back(guh);
	}

	for (int i = 0; i < times.size(); ++i)
	{
		for (int j = 0; j < times.at(i).size(); j++)
		{
			if (i < 48)
			{

				newtimes.at(0).push_back(times.at(i).at(j));
			}

			if (48 <= i && i < 48 + 48)
			{
				newtimes.at(1).push_back(times.at(i).at(j));
			}

			if (48 + 48 <= i && i < 144)
			{
				newtimes.at(2).push_back(times.at(i).at(j));
			}
		}
	}

	times = newtimes;

	vector<TF1> rfuncs;
	vector<TF1> rfunccs;

	vector<double> t0sc;

	for (int i = 0; i < times.size(); i++) // times.size() = number of chamebrs
	{
		// need to also return t0 error -> then will need error fits too(?)
		auto t0_ = fit_t0(times.at(i), i);

		t0 = t0_.t0;
		t0err = t0_.t0err;
		t0fit = *t0_.fitfn;
		TH1F *h = t0_.h;

		TF1 *function = radius_for_time(times.at(i), t0_.t0, i + 1);
		TF1 r_function = *function;

		rfunccs.push_back(*function);
		t0sc.push_back(t0);

		branch1->Fill();
		branch2->Fill();
		branch3->Fill();
	}

	for (int i = 0; i < 144; i++)
	{
		if (i < 48)
		{
			rfuncs.push_back(rfunccs.at(0));
			extern_t0s.push_back(t0sc.at(0));
		}

		if (48 <= i && i < 48 + 48)
		{
			rfuncs.push_back(rfunccs.at(1));
			extern_t0s.push_back(t0sc.at(1));
		}

		if (48 + 48 <= i && i < 144)
		{
			rfuncs.push_back(rfunccs.at(2));
			extern_t0s.push_back(t0sc.at(2));
		}
	}

	return rfuncs;
}
