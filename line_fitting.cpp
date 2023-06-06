#include "line_fitting.h"

static double d0(float x0, float y0, double *par)
{
	// params are in form ax + by + c = 0 ------> (-1/b)*(ax + c) = y

	//	double value = abs(par[0] * x0 + par[1] * y0 + 1) / sqrt(par[0] * par[0] + par[1] * par[1]);
	double value = abs(par[0] * x0 + y0 + par[1]) / sqrt(par[0] * par[0] + 1);

	return value;
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

LineParts justfitlines()
{
	TMinuit *gMinuit = new TMinuit(2);
	gMinuit->SetPrintLevel(-1);
	gMinuit->SetFCN(fcn);

	Double_t arglist[10];
	Int_t ierflg = 0;

	arglist[0] = 2;
	gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
	/// gMinuit->SetErrorDef(1);

	// Set starting values and step sizes for parameters
	static Double_t vstart[3] = {100, 500};
	static Double_t step[3] = {0.01, 0.01};
	gMinuit->mnparm(0, "a", vstart[0], step[0], 0, 0, ierflg);
	gMinuit->mnparm(1, "b", vstart[1], step[1], 0, 0, ierflg);
	// gMinuit->mnparm(2, "c", vstart[2], step[2], 0, 0, ierflg);

	// Now ready for minimization step
	arglist[0] = 1000000;
	arglist[1] = 1.;
	gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

	double a;
	double b;
	// double c;

	double erra;
	double errb;
	// double errc;

	double chisq;

	gMinuit->GetParameter(0, a, erra);
	gMinuit->GetParameter(1, b, errb);
	// gMinuit->GetParameter(2, c, errc);

	chisq = gMinuit->fAmin;

	LineParts loine;
	loine.a = a;
	loine.b = b;
	// loine.c = c;

	loine.aerr = erra;
	loine.berr = errb;
	// loine.cerr = errc;

	loine.chisq = chisq;

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

vector<vector<LineParts>> fit_single_chambers(vector<NewEvent> events, vector<TF1> rfuncs, LineParts &lineparamsc1, LineParts &lineparamsc2, LineParts &lineparamsc3,
						 TBranch *branch_ac1, TBranch *branch_aerrc1, TBranch *branch_bc1, TBranch *branch_berrc1, TBranch *branch_chisqc1,
						 TBranch *branch_ac2, TBranch *branch_aerrc2, TBranch *branch_bc2, TBranch *branch_berrc2, TBranch *branch_chisqc2,
						 TBranch *branch_ac3, TBranch *branch_aerrc3, TBranch *branch_bc3, TBranch *branch_berrc3, TBranch *branch_chisqc3)
{

	std::cout << "fitting " << events.size() << " lines for each chamber " << std::endl;

	vector<vector<LineParts>> cfits;

	vector<LineParts> c1fits;
	vector<LineParts> c2fits;
	vector<LineParts> c3fits;

	for (int i = 0; i < events.size(); i++)
	{
		if (i % 5000 == 0)
			std::cout << "you are at fit for event " << i << " of " << events.size() << std::endl;

		int lastchambnum = 0;

		int numhits = events.at(i).t.size();

		for (int j = 0; j < numhits; j++)
		{
			// TODO
			if (lastchambnum != events.at(i).chamber.at(j))
			{
				switch (lastchambnum)
				{
				case (0):
					lineparamsc1 = justfitlines();
					branch_ac1->Fill();
					branch_aerrc1->Fill();
					branch_bc1->Fill();
					branch_berrc1->Fill();
					branch_chisqc1->Fill();

					c1fits.push_back(lineparamsc1);

					break;
				case (1):
					lineparamsc2 = justfitlines();
					branch_ac2->Fill();
					branch_aerrc2->Fill();
					branch_bc2->Fill();
					branch_berrc2->Fill();
					branch_chisqc2->Fill();

					c2fits.push_back(lineparamsc2);

					break;
				}

				xvs.clear();
				yvs.clear();
				rvs.clear();

				lastchambnum = events.at(i).chamber.at(j);
			}

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

		lineparamsc3 = justfitlines();
		branch_ac3->Fill();
		branch_aerrc3->Fill();
		branch_bc3->Fill();
		branch_berrc3->Fill();
		branch_chisqc3->Fill();

		c3fits.push_back(lineparamsc3);

		xvs.clear();
		yvs.clear();
		rvs.clear();
	}

	cfits.push_back(c1fits);
	cfits.push_back(c2fits);
	cfits.push_back(c3fits);

	return cfits;

}