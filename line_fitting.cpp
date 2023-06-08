#include "line_fitting.h"

static double d0(float x0, float y0, double *par)
{
	// params are in form ax + by + c = 0 ------> (-1/b)*(ax + c) = y

	//	double value = abs(par[0] * x0 + par[1] * y0 + 1) / sqrt(par[0] * par[0] + par[1] * par[1]);
	double value = abs(par[0] * x0 + y0 + par[1]) / sqrt(par[0] * par[0] + 1);

	return value;
}

static double d0_noslope_fit(float x0, float y0, vector<float> a, double *par)
{
	// params are in form ax + y + c = 0 ------> (-ax - c) = y

	//	double value = abs(par[0] * x0 + par[1] * y0 + 1) / sqrt(par[0] * par[0] + par[1] * par[1]);
	double value = abs(a.at(what_entrynum) * x0 + y0 + par[0]) / sqrt(a.at(what_entrynum) * a.at(what_entrynum) + 1);

	return value;
}

static void fcn2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
	const Int_t nbins = xvs.size();
	Int_t i;

	// for res divide 1.46 / time spec width -> 200 / that -> convert to cm
	// 200 because expect 200 micron resolution

	// v_d * 10
	double sigmas = 1; // this is 1mm resolution

	// calculate chisquare
	Double_t chisq = 0;
	Double_t delta;
	for (i = 0; i < nbins; i++)
	{
		delta = (rvs.at(i) - d0_noslope_fit(xvs.at(i), yvs.at(i), as, par)) / sigmas;
		chisq += delta * delta;
	}
	f = chisq;
}

static void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
	const Int_t nbins = xvs.size();
	Int_t i;

	// for res divide 1.46 / time spec width -> 200 / that -> convert to cm
	// 200 because expect 200 micron resolution

	// v_d * 10
	double sigmas = 1; // this is 1mm resolution

	// calculate chisquare
	Double_t chisq = 0;
	Double_t delta;
	for (i = 0; i < nbins; i++)
	{
		delta = (rvs.at(i) - d0(xvs.at(i), yvs.at(i), par)) / sigmas;
		chisq += delta * delta;
	}
	f = chisq;
}

LineParts justfitlines(int setfn)
{
	TMinuit *gMinuit = new TMinuit(2);
	gMinuit->SetPrintLevel(-1);
	if (setfn == 0)
		gMinuit->SetFCN(fcn);
	else
		gMinuit->SetFCN(fcn2);

	Double_t arglist[10];
	Int_t ierflg = 0;

	arglist[0] = 2;
	gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
	/// gMinuit->SetErrorDef(1);

	// Set starting values and step sizes for parameters
	static Double_t vstart[2] = {100, 500};
	static Double_t step[2] = {0.01, 0.01};
	if (setfn == 0)
	{
		gMinuit->mnparm(0, "a", vstart[0], step[0], 0, 0, ierflg);
		gMinuit->mnparm(1, "b", vstart[1], step[1], 0, 0, ierflg);
	}
	else
	{
		gMinuit->mnparm(0, "b", vstart[1], step[1], 0, 0, ierflg);
	}

	// Now ready for minimization step
	arglist[0] = 1000000;
	arglist[1] = 1.;
	gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

	LineParts loine;

	if (setfn == 0)
	{
		double a;
		double b;

		double erra;
		double errb;

		double chisq;

		gMinuit->GetParameter(0, a, erra);
		gMinuit->GetParameter(1, b, errb);

		chisq = gMinuit->fAmin;

		loine.a = a;
		loine.b = b;

		loine.aerr = erra;
		loine.berr = errb;

		loine.chisq = chisq;
	}
	else
	{
		double b;

		double errb;

		double chisq;

		gMinuit->GetParameter(0, b, errb);

		chisq = gMinuit->fAmin;

		loine.b = b;

		loine.berr = errb;

		loine.chisq = chisq;
	}

	delete gMinuit;

	return loine;
}

vector<float> getTubeCoords(int chamber, int layer, int tube)
{
	vector<float> coord;
	coord.push_back(layer == 1 ? 1.5 + tube * 3 : 3 + tube * 3);
	coord.push_back(69.6 - (1.5 + 2.6 * layer + (22.5 + 8.2) * chamber));

	return coord;
}
// THis is the function that fits the lines for each chamber. Does all at once.
void fit_chamber(vector<NewEvent> events, vector<TF1> rfuncs, LineParts &lineparams, TBranch *branch_a, TBranch *branch_aerr, TBranch *branch_b, TBranch *branch_berr, TBranch *branch_chisq)
{
	std::cout << "fitting " << events.size() << " lines for the whole apparatus " << std::endl;
	for (int i = 0; i < events.size(); i++)
	{

		if (i % 5000 == 0)
			std::cout << "you are at fit for event " << i << " of " << events.size() << std::endl;

		int hitsnum = events.at(i).t.size();

		for (int j = 0; j < hitsnum; j++)
		{

			if (events.at(i).is_inlier.at(j) == 1)
			{

				float time = events.at(i).t.at(j);

				vector<float> xy = getTubeCoords(events.at(i).chamber.at(j), events.at(i).layer.at(j), events.at(i).tube.at(j));

				int tubenum = events.at(i).chamber.at(j) * 48 + events.at(i).layer.at(j) * 16 + events.at(i).tube.at(j);

				xvs.push_back(xy.at(0));
				yvs.push_back(xy.at(1));

				rvs.push_back(rfuncs.at(tubenum).Eval(time));
			}
		}

		// fit all chambers
		lineparams = justfitlines();

		branch_a->Fill();
		branch_aerr->Fill();
		branch_b->Fill();
		branch_berr->Fill();
		branch_chisq->Fill();

		xvs.clear();
		yvs.clear();
		rvs.clear();
	}
}

vector<LineParts> fit_single_chamber(int chambernumber, int setfn, vector<NewEvent> events, vector<TF1> rfuncs, LineParts &lineparamsc1,
									 TBranch *branch_ac1, TBranch *branch_aerrc1, TBranch *branch_bc1, TBranch *branch_berrc1, TBranch *branch_chisqc1,
									 TBranch *branch_evnumc1)
{

	std::cout << "fitting " << events.size() << " lines for chamber " << chambernumber << std::endl;

	vector<LineParts> c1fits;

	for (int i = 0; i < events.size(); i++)
	{
		LineParts wow;
		c1fits.push_back(wow);
	}

	LineParts guh;

	for (int i = 0; i < events.size(); i++)
		c1fits.push_back(guh);

	for (int i = 0; i < events.size(); i++)
	{
		if (i % 5000 == 0)
			std::cout << "you are at fit for event " << i << " of " << events.size() << std::endl;

		int numhits = events.at(i).t.size();

		for (int j = 0; j < numhits; j++)
		{
			if (events.at(i).chamber.at(j) == chambernumber && events.at(i).is_inlier.at(j) == 1)
			{
				float time = events.at(i).t.at(j);

				vector<float> xy = getTubeCoords(events.at(i).chamber.at(j), events.at(i).layer.at(j), events.at(i).tube.at(j));

				int tubenum = events.at(i).chamber.at(j) * 48 + events.at(i).layer.at(j) * 16 + events.at(i).tube.at(j);

				xvs.push_back(xy.at(0));
				yvs.push_back(xy.at(1));
				rvs.push_back(rfuncs.at(tubenum).Eval(time));
			}
		}

		lineparamsc1 = justfitlines(setfn);
		lineparamsc1.eventNum = i;

		branch_ac1->Fill();
		branch_aerrc1->Fill();
		branch_bc1->Fill();
		branch_berrc1->Fill();
		branch_chisqc1->Fill();
		branch_evnumc1->Fill();

		c1fits.at(i) = (lineparamsc1);

		xvs.clear();
		yvs.clear();
		rvs.clear();
	}

	return c1fits;
}

vector<LineParts> fit_single_chamber(int chambernumber, int setfn, vector<NewEvent> events, vector<TF1> rfuncs, LineParts &lineparamsc1, TBranch *branch_bc1, TBranch *branch_berrc1)
{

	std::cout << "fitting " << events.size() << " lines for chamber " << chambernumber << std::endl;

	vector<LineParts> c1fits;

	for (int i = 0; i < events.size(); i++)
	{
		LineParts wow;
		c1fits.push_back(wow);
	}

	LineParts guh;

	for (int i = 0; i < events.size(); i++)
		c1fits.push_back(guh);

	for (int i = 0; i < events.size(); i++)
	{

		if (i % 5000 == 0)
			std::cout << "you are at fit for event " << i << " of " << events.size() << std::endl;

		int numhits = events.at(i).t.size();

		if (as.at(i) == 0)
		{
			what_entrynum += 1;
			continue;
		}

		for (int j = 0; j < numhits; j++)
		{
			if (events.at(i).chamber.at(j) == chambernumber && events.at(i).is_inlier.at(j) == 1)
			{
				float time = events.at(i).t.at(j);

				vector<float> xy = getTubeCoords(events.at(i).chamber.at(j), events.at(i).layer.at(j), events.at(i).tube.at(j));

				int tubenum = events.at(i).chamber.at(j) * 48 + events.at(i).layer.at(j) * 16 + events.at(i).tube.at(j);

				xvs.push_back(xy.at(0));
				yvs.push_back(xy.at(1));
				rvs.push_back(rfuncs.at(tubenum).Eval(time));
			}
		}

		lineparamsc1 = justfitlines(setfn);
		lineparamsc1.eventNum = i;

		branch_bc1->Fill();
		branch_berrc1->Fill();

		c1fits.at(i) = (lineparamsc1);
		xvs.clear();
		yvs.clear();
		rvs.clear();

		what_entrynum += 1;
	}

	return c1fits;
}