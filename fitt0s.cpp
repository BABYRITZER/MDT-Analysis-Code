// given a time distribution h, fit_t0 returns a struct with t0 and t0s error and the fitted graph -- populate t0s will populate everything
// input not a pointer as to not directly change the graph -- will return the fit function (TF1*) at the end in the struct
#include "fitt0s.h"

// pass in the vector of the time spectrum
Fittedt0s fit_t0(vector<float> times)
{

	TH1F *h;
	h = new TH1F("", "...", 230, 0, 3600);

	for (int k = 0; k < times.size(); k++)
		h->Fill(times.at(k));

	// set the fitting function
	TF1 *fitfn;
	fitfn = new TF1("fit", "([2] / ( 1 + TMath::Exp( -[0]*(x-[1] ) ) ) )");

	// do the fits
	int binmax = h->GetMaximumBin();
	double t0_guess = h->GetXaxis()->GetBinCenter(binmax); // this initial t0 is a guess

	// set guess parameters
	fitfn->SetParameters(1, t0_guess, binmax);

	// here we fit
	h->SetBit(TH1::kNoStats);
	h->Fit(fitfn, "0", "S", 0, t0_guess + 15);
	h->ResetBit(TH1::kNoStats);

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

	vector<TF1> rfuncs;

	for (int i = 0; i < times.size(); i++) // times.size() = #tubes -- 144
	{
		// need to also return t0 error -> then will need error fits too(?)
		auto t0_ = fit_t0(times.at(i));

		t0 = t0_.t0;
		t0err = t0_.t0err;
		t0fit = *t0_.fitfn;

		TF1 *function = radius_for_time(times.at(i), t0_.t0);
		TF1 r_function = *function;

		rfuncs.push_back(*function);
		extern_t0s.push_back(t0);

		branch1->Fill();
		branch2->Fill();
		branch3->Fill();

	}

	return rfuncs;
}
