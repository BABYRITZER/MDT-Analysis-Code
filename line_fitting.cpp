#include "line_fitting.h"

static double d0(float x0, float y0, double *par)
{
	// params are in form ax + by + c = 0 ------> (-1/b)*(ax + c) = y

	//	double value = abs(par[0] * x0 + par[1] * y0 + 1) / sqrt(par[0] * par[0] + par[1] * par[1]);
	double value = abs(x0 + par[0] * y0 + par[1]) / sqrt(par[0] * par[0] + 1);

	return value;
}

static double d0_noslope_fit(float x0, float y0, vector<float> a, double *par)
{
	// params are in form ax + y + c = 0 ------> (-ax - c) = y

	//	double value = abs(par[0] * x0 + par[1] * y0 + 1) / sqrt(par[0] * par[0] + par[1] * par[1]);
	double value = abs(x0 + a.at(what_entrynum) * y0 + par[0]) / sqrt(a.at(what_entrynum) * a.at(what_entrynum) + 1);

	return value;
}

static void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
	const Int_t nbins = xvs.size();
	Int_t i;

	// for res do v_d * 10

	// calculate chisquare
	Double_t chisq = 0;
	Double_t delta;
	for (i = 0; i < nbins; i++)
	{
		float sigma = sigmas.at(i);

		// TODO: CHANGED SIGMA ON CHISQUARE

		delta = (rvs.at(i) - d0(xvs.at(i), yvs.at(i), par)) / (double)sigma;
		// delta = (rvs.at(i) - d0(xvs.at(i), yvs.at(i), par)) / 0.1;

		chisq += delta * delta;

		// sstd::cout << rvs.at(i) << " " << d0(xvs.at(i), yvs.at(i), par) << " " << xvs.at(i) << " " << yvs.at(i) << " " << par[0] << " " << par[1] << std::endl;
	}
	// std::cout << "chisq is " << chisq << std::endl;
	f = chisq;
}

LineParts justfitlines(int setfn)
{
	TMinuit *gMinuit = new TMinuit(2);
	// TODO: changed print level of minuit
	gMinuit->SetPrintLevel(-1);
	gMinuit->SetFCN(fcn);

	Double_t arglist[10];
	Int_t ierflg = 0;

	arglist[0] = 2;
	gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
	// TODO: TRYING THIS, WHATEVER IT IS
	gMinuit->SetErrorDef(1);

	// Set starting values and step sizes for parameters
	// static Double_t vstart[2] = {0, 25};
	static Double_t vstart[2] = {0, 0};
	// TODO: TESTING IF STARTING AT 0 0 GETS RID OF THE WEIRD ASS BUMP
	static Double_t step[2] = {0.001, 0.001};
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
	// TODO: CHANGED MINIMIZATION TOLERANCE -- before it was like 1
	arglist[1] = 0.00001;

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

// TODO: i think that maybe the fitting is whats making everything kind of go bad but i am not sure i need to subtract or rotate first

// This is the function that fits the lines for each chamber. Does all at once.
// vector<LineParts> fit_chamber(vector<NewEvent> events, vector<TF1> rfuncs, LineParts &lineparams, TBranch *branch_a, TBranch *branch_aerr, TBranch *branch_b, TBranch *branch_berr, TBranch *branch_chisq, float meanc1_angle, float meanc3_angle, float meanc1_b_diffs, float meanc3_b_diffs)
vector<LineParts> fit_chamber(vector<NewEvent> events, vector<TF1> rfuncs, LineParts &lineparams, TBranch *branch_a, TBranch *branch_aerr, TBranch *branch_b, TBranch *branch_berr, TBranch *branch_chisq, float c1_angle, float c3_angle, float meanc1_b_diffs, float meanc3_b_diffs)
{
	vector<LineParts> lines;
	std::cout << "fitting " << events.size() << " lines for the whole apparatus " << std::endl;
	for (int i = 0; i < events.size(); i++)
	{

		if (i % 5000 == 0)
			std::cout << "you are at fit for event " << i << " of " << events.size() << std::endl;

		int hitsnum = events.at(i).t.size();

		double rotationmatc1[2][2] = {{cos(c1_angle), -sin(c1_angle)}, {sin(c1_angle), cos(c1_angle)}};
		double rotationmatc3[2][2] = {{cos(c3_angle), -sin(c3_angle)}, {sin(c3_angle), cos(c3_angle)}};

		for (int j = 0; j < hitsnum; j++)
		{

			if (events.at(i).is_inlier.at(j) == 1)
			{

				float time = events.at(i).t.at(j);
				// TODO: DO WE ROTATE OR TRANSLATE FIRST
				vector<float> xy = getTubeCoords(events.at(i).chamber.at(j), events.at(i).layer.at(j), events.at(i).tube.at(j));

				// rotation first
				/*
								if (events.at(i).chamber.at(j) == 0)
								{
									xy.at(0) = (rotationmatc1[0][0] * xy.at(0) + rotationmatc1[0][1] * xy.at(1)) - meanc1_b_diffs;
									xy.at(1) = rotationmatc1[1][0] * xy.at(0) + rotationmatc1[1][1] * xy.at(1);
								}
								if (events.at(i).chamber.at(j) == 2)
								{
									xy.at(0) = (rotationmatc3[0][0] * xy.at(0) + rotationmatc3[0][1] * xy.at(1)) - meanc3_b_diffs;
									xy.at(1) = rotationmatc3[1][0] * xy.at(0) + rotationmatc3[1][1] * xy.at(1);
								}
				*/

				// translation first
				if (events.at(i).chamber.at(j) == 0)
				{
					xy.at(0) = xy.at(0) - meanc1_b_diffs;
					xy.at(0) = (rotationmatc1[0][0] * xy.at(0) + rotationmatc1[0][1] * xy.at(1));
					xy.at(1) = rotationmatc1[1][0] * xy.at(0) + rotationmatc1[1][1] * xy.at(1);
				}
				if (events.at(i).chamber.at(j) == 2)
				{
					xy.at(0) = xy.at(0) - meanc1_b_diffs;
					xy.at(0) = (rotationmatc3[0][0] * xy.at(0) + rotationmatc3[0][1] * xy.at(1)) - meanc3_b_diffs;
					xy.at(1) = rotationmatc3[1][0] * xy.at(0) + rotationmatc3[1][1] * xy.at(1);
				}

				int tubenum = events.at(i).chamber.at(j) * 48 + events.at(i).layer.at(j) * 16 + events.at(i).tube.at(j);

				xvs.push_back(xy.at(0));
				yvs.push_back(xy.at(1));

				rvs.push_back(rfuncs.at(tubenum).Eval(time - extern_t0s.at(tubenum)));
				sigmas.push_back(rfuncs.at(tubenum).Derivative(time - extern_t0s.at(tubenum)) * 25. / sqrt(12.));
			}
		}

		// fit all chambers
		lineparams = justfitlines();

		lines.push_back(lineparams); // save the line parameters

		branch_a->Fill();
		branch_aerr->Fill();
		branch_b->Fill();
		branch_berr->Fill();
		branch_chisq->Fill();

		xvs.clear();
		yvs.clear();
		rvs.clear();
		sigmas.clear();
	}

	return lines;
}

// SINGLE CHAMBER FIT, RETURNS EVERY PARAMETER INTO A BRANCH
// ONLY FITS IF THE CHAMBER HAS ALL THREE HITS
vector<LineParts> fit_single_chamber(int chambernumber, int setfn, vector<NewEvent> events, vector<TF1> rfuncs, LineParts &lineparamsc1,
									 TBranch *branch_ac1, TBranch *branch_aerrc1, TBranch *branch_bc1, TBranch *branch_berrc1, TBranch *branch_chisqc1,
									 TBranch *branch_evnumc1)
{

	std::cout << "fitting " << events.size() << " lines for chamber " << chambernumber << std::endl;
	// TODO
	vector<LineParts> c1fits;

	for (int i = 0; i < events.size(); i++)
	{
		LineParts wow;
		c1fits.push_back(wow);
	}

	LineParts guh;

	for (int i = 0; i < events.size(); i++)
	{

		if (i % 5000 == 0)
			std::cout << "you are at fit for event " << i << " of " << events.size() << std::endl;

		LineParts wow;

		int numhits = events.at(i).t.size();

		int chambhits = 0;
		for (int z = 0; z < numhits; ++z)
		{
			if (events.at(i).chamber.at(z) == chambernumber)
				chambhits = chambhits + 1;
		}

		if (chambhits == 3)
		{
			for (int j = 0; j < numhits; j++)
			{
				if (events.at(i).chamber.at(j) == chambernumber && events.at(i).is_inlier.at(j) == 1)
				{

					float time = events.at(i).t.at(j);

					vector<float> xy = getTubeCoords(events.at(i).chamber.at(j), events.at(i).layer.at(j), events.at(i).tube.at(j));

					int tubenum = events.at(i).chamber.at(j) * 48 + events.at(i).layer.at(j) * 16 + events.at(i).tube.at(j);

					xvs.push_back(xy.at(0));
					yvs.push_back(xy.at(1));
					rvs.push_back(rfuncs.at(tubenum).Eval(time - extern_t0s.at(tubenum)));
					sigmas.push_back(rfuncs.at(tubenum).Derivative(time - extern_t0s.at(tubenum)) * 25. / sqrt(12.));
				}
			}

			// TODO: HAVE SET IT SO THAT IT ONLY FITS LINES IF THERE ARE MORE THAN OR EQUAL TO TWO POINTS IN THE CHAMBER
			// if (xvs.size() > 2)
			wow = justfitlines(setfn);
		}

		lineparamsc1 = wow;

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
		sigmas.clear();
	}

	return c1fits;
}

// SINGLE CHAMBER FIT, RETURNS ONLY INTERCEPT ---->ROTATION
vector<LineParts> fit_single_chamber(int chambernumber, int setfn, double rotationangle, vector<NewEvent> events, vector<TF1> rfuncs, LineParts &lineparamsc1, TBranch *branch_bc1, TBranch *branch_berrc1)
{

	std::cout << "fitting " << events.size() << " lines for chamber " << chambernumber << std::endl;

	vector<LineParts> c1fits;

	for (int i = 0; i < events.size(); i++)
	{
		LineParts wow;
		c1fits.push_back(wow);
	}

	double rotationmat[2][2] = {{cos(rotationangle), -sin(rotationangle)}, {sin(rotationangle), cos(rotationangle)}};

	for (int i = 0; i < events.size(); i++)
	{
		LineParts wow;

		if (i % 5000 == 0)
			std::cout << "you are at fit for event " << i << " of " << events.size() << std::endl;

		int numhits = events.at(i).t.size();

		int chambhits = 0;
		for (int z = 0; z < numhits; ++z)
		{
			if (events.at(i).chamber.at(z) == chambernumber)
				chambhits = chambhits + 1;
		}

		if (chambhits == 3)
		{
			for (int j = 0; j < numhits; j++)
			{
				if (events.at(i).chamber.at(j) == chambernumber && events.at(i).is_inlier.at(j) == 1)
				{
					float time = events.at(i).t.at(j);

					vector<float> xy = getTubeCoords(events.at(i).chamber.at(j), events.at(i).layer.at(j), events.at(i).tube.at(j));

					xy.at(0) = rotationmat[0][0] * xy.at(0) + rotationmat[0][1] * xy.at(1);
					xy.at(1) = rotationmat[1][0] * xy.at(0) + rotationmat[1][1] * xy.at(1);

					int tubenum = events.at(i).chamber.at(j) * 48 + events.at(i).layer.at(j) * 16 + events.at(i).tube.at(j);

					xvs.push_back(xy.at(0));
					yvs.push_back(xy.at(1));
					rvs.push_back(rfuncs.at(tubenum).Eval(time - extern_t0s.at(tubenum)));
					sigmas.push_back(rfuncs.at(tubenum).Derivative(time - extern_t0s.at(tubenum)) * 25. / sqrt(12.));
				}
			}

			wow = justfitlines(setfn);
		}

		lineparamsc1 = wow;

		lineparamsc1.eventNum = i;

		branch_bc1->Fill();
		branch_berrc1->Fill();

		c1fits.at(i) = (lineparamsc1);
		xvs.clear();
		yvs.clear();
		rvs.clear();
		sigmas.clear();

		what_entrynum += 1;
	}

	return c1fits;
}
