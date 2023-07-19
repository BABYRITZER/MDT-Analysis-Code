#include "line_fitting.h"

LineParts subtractlineparts(LineParts line1, LineParts line2)
{
	LineParts wow;

	wow.a = line1.a - line2.a;
	wow.b = line1.b - line2.b;
	wow.c = line1.c - line2.c;

	wow.aerr = sqrt(line1.aerr * line1.aerr + line2.aerr * line2.aerr);
	wow.berr = sqrt(line1.berr * line1.berr + line2.berr * line2.berr);
	wow.cerr = sqrt(line1.cerr * line1.cerr + line2.cerr * line2.cerr);

	wow.eventNum = line1.eventNum;

	return wow;
}

// Custom right-tailed function
double RightTailedFunction(double *x, double *par)
{
	// Parameters
	double mean = par[0];	   // Mean of the distribution
	double sigma = par[1];	   // Standard deviation of the distribution
	double tailParam = par[2]; // Parameter specific to the right-tail distribution
	double amplitude = par[3];

	// Calculate the value of the function
	double arg = (x[0] - mean) / sigma;
	double exponent = TMath::Exp(-0.5 * arg * arg);
	// double tail = TMath::Power(x[0], tailParam);

	double tail = 0.5 * (1. + TMath::Erf(tailParam * arg / sqrt(2)));

	double asdf = amplitude * exponent * tail;

	return asdf;
}

static double d0(float x0, float y0, double *par, int chamber)
{
	double value;

	// double value = abs(x0 + par[0] * y0 - par[1]) / sqrt(par[0] * par[0] + 1.);

	// c1 shift is par[2]
	// c3 shift is par[3]

	value = abs(par[1] + par[0] * y0 - x0) / sqrt(par[0] * par[0] + 1.);

	return value;
}

static void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
	const Int_t nbins = xvs.size();
	Int_t i;

	// calculate chisquare
	Double_t chisq = 0;
	Double_t delta;

	for (i = 0; i < nbins; i++)
	{
		float sigma = sigmas.at(i);

		// delta = (rvs.at(i) - d0(xvs.at(i), yvs.at(i), par)) / (double)sigma;

		delta = (rvs.at(i) - d0(xvs.at(i), yvs.at(i), par, chlist.at(i))) / 0.1;

		chisq += delta * delta;

		// std::cout << "tube is: " << i << " delta is: " << delta << " sigma is: " << sigma << std::endl;
		// std::cout << "radius is " << rvs.at(i) << " d0 is " << d0(xvs.at(i), yvs.at(i), par) << " x is: " << xvs.at(i) << " y is: " << yvs.at(i) << " par[0] is: " << par[0] << " par[1] is: " << par[1] << std::endl;
	}
	// std::cout << "chisq is " << chisq << std::endl;
	f = (1. / ((double)xvs.size() - 2.)) * chisq;
}

LineParts justfitlines(int setfn, int event, vector<float> gransac_lineparams)
{
	TMinuit *gMinuit = new TMinuit(2);

	gMinuit->SetPrintLevel(-1);
	gMinuit->SetFCN(fcn);

	Double_t arglist[10];
	Int_t ierflg = 0;

	arglist[0] = 1;
	gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

	gMinuit->SetErrorDef(1);

	float x0 = 0;
	float m = 0;

	if (0.5 * (xvs.at(0) + xvs.at(2)) < xvs.at(1))
	{
		x0 = xvs.at(1) - rvs.at(1);
		m = (xvs.at(0) + rvs.at(0) - xvs.at(2) + rvs.at(2)) / (yvs.at(2) - yvs.at(0));
	}

	else
	{
		x0 = xvs.at(1) + rvs.at(1);
		m = ((xvs.at(0) - rvs.at(0)) - (xvs.at(2) - rvs.at(2))) / (yvs.at(2) - yvs.at(0));
	}

	x0 = x0 + m * yvs.at(1);

	static Double_t vstart[2];

	vstart[0] = m;
	vstart[1] = x0;

	static Double_t step[2] = {0.00001, 0.001};

	gMinuit->mnparm(0, "a", vstart[0], step[0], -0.5, 0.5, ierflg);
	gMinuit->mnparm(1, "b", vstart[1], step[1], -10, 50, ierflg);

	// gMinuit->mnparm(2, "ch1loc", 0., step[1], -5., 5., ierflg);
	// gMinuit->mnparm(3, "ch3loc", 0., step[1], -5., 5., ierflg);

	// Now ready for minimization step
	arglist[0] = 1000000;

	// arglist[1] = 0.001;
	arglist[1] = 0.000000000001;

	gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

	LineParts loine;

	if (setfn == 0)
	{
		double a;
		double b;

		double ch1;
		double ch3;

		double erra;
		double errb;

		double errch1;
		double errch3;

		double chisq;

		gMinuit->GetParameter(0, a, erra);
		gMinuit->GetParameter(1, b, errb);

		// gMinuit->GetParameter(2, ch1, errch1);
		// gMinuit->GetParameter(3, ch3, errch3);

		chisq = gMinuit->fAmin;

		loine.a = a;
		loine.b = b;

		// loine.ch1 = ch1;
		// loine.ch3 = errch3;

		loine.aerr = erra;
		loine.berr = errb;

		// loine.ch1err = errch1;
		// loine.ch3err = errch3;

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

	// std::cout << "---------------------------------------------------" << std::endl;
	delete gMinuit;

	return loine;
}

// Gets the coordinates of the tube based on the chamber, layer, and tube number. Chamber = 0, 1, 2. Layer = 0, 1, 2 Tube = 0-15.
vector<float> getTubeCoords(int chamber, int layer, int tube)
{
	vector<float> coord;
	coord.push_back(layer == 1 ? 1.5 + tube * 3 : 3 + tube * 3);
	coord.push_back(69.6 - (1.5 + 2.6 * layer + (22.5 + 8.2) * chamber));

	return coord;
}

// This is the function that fits the lines for each chamber. Does all at once.
vector<LineParts> fit_chamber(vector<NewEvent> events, vector<TF1> rfuncs, LineParts &lineparams, TBranch *branch_a, TBranch *branch_aerr,
							  TBranch *branch_b, TBranch *branch_berr,
							  TBranch *branch_chisq, float c1_angle, float c3_angle, float meanc1_b_diffs,
							  float meanc3_b_diffs, vector<float> gransac_lineparams, TBranch *branch_ch1,
							  TBranch *branch_ch1err, TBranch *branch_ch3, TBranch *branch_ch3err)
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

				vector<float> xy = getTubeCoords(events.at(i).chamber.at(j), events.at(i).layer.at(j), events.at(i).tube.at(j));

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
				chlist.push_back(events.at(i).chamber.at(j));
			}
		}

		// fit all chambers
		if (xvs.size() > 4)
			lineparams = justfitlines(0, i, gransac_lineparams);

		lines.push_back(lineparams); // save the line parameters

		branch_a->Fill();
		branch_aerr->Fill();
		branch_b->Fill();
		branch_berr->Fill();
		branch_chisq->Fill();

		branch_ch1->Fill();
		branch_ch1err->Fill();
		branch_ch3->Fill();
		branch_ch3err->Fill();

		xvs.clear();
		yvs.clear();
		rvs.clear();
		chlist.clear();
		sigmas.clear();
	}

	return lines;
}

// This version of the function ignores the specified layer number.
vector<LineParts> fit_chamber(vector<NewEvent> events, vector<TF1> rfuncs, int layer_to_ignore)
{
	vector<LineParts> lines;
	// std::cout << "fitting " << events.size() << " lines for the whole apparatus " << std::endl;
	for (int i = 0; i < events.size(); i++)
	{

		float c1_angle = 0.;
		float c3_angle = 0;
		float meanc1_b_diffs = 0.;
		float meanc3_b_diffs = 0.;

		// if (i % 5000 == 0)
		// std::cout << "you are at fit for event " << i << " of " << events.size() << std::endl;

		int hitsnum = events.at(i).t.size();

		double rotationmatc1[2][2] = {{cos(c1_angle), -sin(c1_angle)}, {sin(c1_angle), cos(c1_angle)}};
		double rotationmatc3[2][2] = {{cos(c3_angle), -sin(c3_angle)}, {sin(c3_angle), cos(c3_angle)}};

		for (int j = 0; j < hitsnum; j++)
		{

			if (events.at(i).is_inlier.at(j) == 1)
			{
				int layernum = events.at(i).layer.at(j) + events.at(i).chamber.at(j) * 3; // We need to get the absolute layer number 0-8 instead of relative to the chamber 0-2

				if (layernum == layer_to_ignore) // Important: ignore the specified layer number
					continue;

				float time = events.at(i).t.at(j);

				vector<float> xy = getTubeCoords(events.at(i).chamber.at(j), events.at(i).layer.at(j), events.at(i).tube.at(j));

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
				// TODO:
				rvs.push_back(rfuncs.at(tubenum).Eval(time - extern_t0s.at(tubenum)));
				sigmas.push_back(rfuncs.at(tubenum).Derivative(time - extern_t0s.at(tubenum)) * 25. / sqrt(12.));
				chlist.push_back(events.at(i).chamber.at(j));
			}
		}

		LineParts lineparams;

		// fit all chambers
		if (xvs.size() >= 5)
			lineparams = justfitlines(0, i, vector<float>());

		lines.push_back(lineparams); // save the line parameters

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
									 TBranch *branch_evnumc1, vector<float> gransac_lineparams)
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

		/*
				int chambhits = 0;
				for (int z = 0; z < numhits; ++z)
				{
					if (events.at(i).chamber.at(z) == chambernumber)
						chambhits = chambhits + 1;
				}
		*/
		// if (chambhits == 3)
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
			if (xvs.size() == 3)
				wow = justfitlines(setfn, i, gransac_lineparams);
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

		// std::cout << "event is: " << i << std::endl;

		// if (i > 10)
		//{
		//	std::cout << "done" << std::endl;
		//	return c1fits;
		// }
	}

	return c1fits;
}

// SINGLE CHAMBER FIT, RETURNS ONLY INTERCEPT ---->ROTATION
vector<LineParts> fit_single_chamber(int chambernumber, int setfn, double rotationangle, vector<NewEvent> events, vector<TF1> rfuncs, LineParts &lineparamsc1, TBranch *branch_bc1, TBranch *branch_berrc1)
{

	std::cout << "fitting " << events.size() << " lines for chamber " << chambernumber << std::endl;

	vector<LineParts> c1fits;
	// fitguff is to pass into the justfitlines fn
	vector<float> fitguff;

	for (int i = 0; i < events.size(); i++)
	{
		LineParts wow;
		c1fits.push_back(wow);
		fitguff.push_back(0);
		fitguff.push_back(0);
	}

	double rotationmat[2][2] = {{cos(rotationangle), -sin(rotationangle)}, {sin(rotationangle), cos(rotationangle)}};

	for (int i = 0; i < events.size(); i++)
	{
		LineParts wow;

		if (i % 5000 == 0)
			std::cout << "you are at fit for event " << i << " of " << events.size() << std::endl;

		int numhits = events.at(i).t.size();

		// int chambhits = 0;
		// for (int z = 0; z < numhits; ++z)
		//{
		//	if (events.at(i).chamber.at(z) == chambernumber)
		//		chambhits = chambhits + 1;
		// }

		// if (chambhits == 3)
		//{
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

		if (xvs.size() == 3)
			wow = justfitlines(setfn, i, fitguff);
		//}

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

float fittwochambers(int c1, int c2, int c_fit, vector<NewEvent> events, vector<TF1> rfuncs)
{
	vector<float> fitguff;

	for (int i = 0; i < events.size(); i++)
	{
		fitguff.push_back(0);
		fitguff.push_back(0);
	}

	vector<LineParts> all_fits;

	string graftitle = "Distance of chamber " + std::to_string(c_fit) + " points to " + std::to_string(c1) + " and " + std::to_string(c2) + "fit";

	TH1F *chambdistances = new TH1F("", graftitle.c_str(), 200, -15., 15.);

	std::cout << "fitting " << events.size() << " lines for chambers " << c1 << " and " << c2 << std::endl;

	for (int i = 0; i < events.size(); i++)
	{

		LineParts wow;

		if (i % 5000 == 0)
			std::cout << "you are at fit for event " << i << " of " << events.size() << std::endl;

		int numhits = events.at(i).t.size();

		vector<float> offchamberhitsx;
		vector<float> offchamberhitsy;
		vector<float> offchamberhitsrads;

		for (int j = 0; j < numhits; j++)
		{
			if ((events.at(i).chamber.at(j) == c1 || events.at(i).chamber.at(j) == c2) && events.at(i).is_inlier.at(j) == 1)
			{
				float time = events.at(i).t.at(j);

				vector<float> xy = getTubeCoords(events.at(i).chamber.at(j), events.at(i).layer.at(j), events.at(i).tube.at(j));

				int tubenum = events.at(i).chamber.at(j) * 48 + events.at(i).layer.at(j) * 16 + events.at(i).tube.at(j);

				xvs.push_back(xy.at(0));
				yvs.push_back(xy.at(1));
				rvs.push_back(rfuncs.at(tubenum).Eval(time - extern_t0s.at(tubenum)));
				sigmas.push_back(rfuncs.at(tubenum).Derivative(time - extern_t0s.at(tubenum)) * 25. / sqrt(12.));
			}
			else if (events.at(i).chamber.at(j) == c_fit && events.at(i).is_inlier.at(j) == 1)
			{
				float time = events.at(i).t.at(j);

				vector<float> xy = getTubeCoords(events.at(i).chamber.at(j), events.at(i).layer.at(j), events.at(i).tube.at(j));

				int tubenum = events.at(i).chamber.at(j) * 48 + events.at(i).layer.at(j) * 16 + events.at(i).tube.at(j);

				offchamberhitsx.push_back(xy.at(0));
				offchamberhitsy.push_back(xy.at(1));
				offchamberhitsrads.push_back(rfuncs.at(tubenum).Eval(time - extern_t0s.at(tubenum)));
			}
		}

		if (xvs.size() > 4)
			wow = justfitlines(0, i, fitguff);

		all_fits.push_back(wow);

		if (wow.a != 0)
		{
			for (int w = 0; w < offchamberhitsx.size(); ++w)
			{
				// float d1 = abs(-offchamberhitsx.at(w) + wow.a * offchamberhitsy.at(w) + wow.b) / sqrt(wow.a * wow.a + 1.) + offchamberhitsrads.at(w);
				// float d2 = abs(-offchamberhitsx.at(w) + wow.a * offchamberhitsy.at(w) + wow.b) / sqrt(wow.a * wow.a + 1.) - offchamberhitsrads.at(w);

				// chambdistances->Fill(abs(d1) < abs(d2) ? d1 : d2);

				float d1 = (-offchamberhitsx.at(w) + wow.a * offchamberhitsy.at(w) + wow.b) / sqrt(wow.a * wow.a + 1.) - offchamberhitsrads.at(w);

				chambdistances->Fill(d1);
			}
		}

		xvs.clear();
		yvs.clear();
		rvs.clear();
		sigmas.clear();
	}

	/*	 // Define the custom function
		TF1 *rightTailed = new TF1("rightTailed", RightTailedFunction, -1.5, 15., 4);

		rightTailed->SetParameters(-.6, 4, 5.5, 2400);
		rightTailed->SetRange(-15,15);

		// Perform the fitting
		chambdistances->Fit(rightTailed, "R");
	*/

	vector<TF1 *> fits;

	for (int c = 0; c < 50; ++c)
	{
		double par[6];

		double rangemin = -10 + (double)(20. * c) / 50.;

		string total_title = "total_" + std::to_string(c);
		string g1_title = "g1_" + std::to_string(c);
		string g2_title = "g2_" + std::to_string(c);

		TF1 *leftgaus = new TF1(g1_title.c_str(), "gaus", -15., rangemin);
		TF1 *rightgaus = new TF1(g2_title.c_str(), "gaus", rangemin, 15.);

		TF1 *total = new TF1(total_title.c_str(), "gaus(0)+gaus(3)", -15, 15);

		chambdistances->Fit(leftgaus, "R");
		chambdistances->Fit(rightgaus, "R+");

		leftgaus->GetParameters(&par[0]);
		rightgaus->GetParameters(&par[3]);

		total->SetParameters(par);
		total->SetLineColor(7);

		chambdistances->Fit(total, "R+");

		fits.push_back(total);
	}

	double lowest_chisq = 1000000000000.;
	int best_fit = 0;

	for (int c = 0; c < fits.size(); ++c)
	{
		double chi = fits.at(c)->GetChisquare();

		if (chi < lowest_chisq)
		{
			best_fit = c;
			lowest_chisq = chi;
			std::cout << chi << " " << lowest_chisq << std::endl;
		}
	}

	TCanvas *graph = new TCanvas();

	string bestfit_title = "total_" + std::to_string(best_fit);

	chambdistances->GetXaxis()->SetTitle("Distance to fitted line");
	chambdistances->GetYaxis()->SetTitle("Counts");

	fits.at(best_fit)->Draw();
	chambdistances->Draw();

	graph->Modified();
	graph->Update();

	string title = "./Output/chambers_" + std::to_string(c1) + "_" + std::to_string(c2) + "_distances.png";

	graph->SaveAs(title.c_str());

	std::cout << chambdistances->GetFunction(bestfit_title.c_str())->GetParameter(4) << std::endl;

	return 0; //( chambdistances->GetFunction(bestfit_title.c_str())->GetParameter(1) + chambdistances->GetFunction(bestfit_title.c_str())->GetParameter(4) ) / 2.;
}
