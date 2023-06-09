#include "Event.h"
#include "line_fitting.h"
#include "load_from_root_file.h"
#include "fitt0s.h"
#include "gransac_implementations.h"
#include <TApplication.h>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TF1.h>
#include <iostream>
#include <fstream>
#include <vector>

using std::cout;
using std::endl;
using std::string;
using std::vector;

vector<float> xvs;
vector<float> yvs;
vector<float> rvs;
vector<float> as;
int what_entrynum;

// int main(int argc, char **argv)
int main()
{

	//	TApplication app("Root app", &argc, argv);

	std::cout << "type the run number" << std::endl;

	vector<string> filenames;
	string foilenaem;
	std::cin >> foilenaem;
	filenames.push_back("reco_run" + foilenaem);

	std::cout << "Starting program" << std::endl;

	std::cout << filenames.size() << std::endl;

	int whatfile = 0;

	std::cout << "Opening file " << std::endl;
	// open the tree
	string name = filenames.at(whatfile) + "_analysis.root";
	TFile *file = new TFile(name.c_str(), "RECREATE");
	TTree *tree1 = new TTree("myTree", "My Tree Title");

	std::cout << "Loading event info " << std::endl;

	// get the events
	string filename = filenames.at(whatfile) + ".root";
	vector<NewEvent> events = loadTreeFromFile(filename);

	file->cd();

	// populating the t0 for each tube
	vector<TF1> rfuncs;
	vector<TF1 *> rfuncsst;
	double t0;
	double t0err;
	TF1 t0fit;

	TBranch *brancht0 = tree1->Branch("t0s", &t0);
	TBranch *brancht0_err = tree1->Branch("t0_errors", &t0err);
	TBranch *brancht0_fits = tree1->Branch("t0_fits", &t0fit);
	TBranch *branchr_funcs = tree1->Branch("r_functions", &rfuncsst);

	rfuncs = populate_t0s(events, t0, t0err, t0fit, brancht0, brancht0_err, brancht0_fits);

	for (int q = 0; q < 144; q++)
		rfuncsst.push_back(&rfuncs.at(q));

	branchr_funcs->Fill();

	std::cout << "finding outlier points" << std::endl;
	// get number of events and add them to tree
	long int NEvents = events.size();

	std::cout << NEvents << " events have been loaded " << std::endl;
	// pattern recog here
	// populates the tree with new entries with pattern recognized hits

	NewEvent event;

	TBranch *branch_eventnum = tree1->Branch("event_num", &event.eventnum);
	TBranch *branch_t = tree1->Branch("time", &event.t);
	TBranch *branch_chg = tree1->Branch("charge", &event.charge);
	TBranch *branch_chmb = tree1->Branch("chamber", &event.chamber);
	TBranch *branch_layer = tree1->Branch("layer", &event.layer);
	TBranch *branch_tube = tree1->Branch("tube", &event.tube);
	TBranch *branch_is_inliner = tree1->Branch("is_inliner", &event.is_inlier);

	findInliers(events, event, branch_eventnum, branch_t, branch_chg, branch_chmb, branch_layer, branch_tube, branch_is_inliner);

	std::cout << "fitting chambers " << std::endl;

	xvs.clear();
	yvs.clear();
	rvs.clear();

	// fit single chambers
	LineParts lineparamsc1;
	LineParts lineparamsc2;
	LineParts lineparamsc3;

	TBranch *branch_ac1 = tree1->Branch("a_c1", &lineparamsc1.a);
	TBranch *branch_aerrc1 = tree1->Branch("aerr_c1", &lineparamsc1.aerr);
	TBranch *branch_bc1 = tree1->Branch("b_c1", &lineparamsc1.b);
	TBranch *branch_berrc1 = tree1->Branch("berr_c1", &lineparamsc1.berr);
	TBranch *branch_chisqc1 = tree1->Branch("c1_chisq", &lineparamsc1.chisq);
	TBranch *branch_evnumc1 = tree1->Branch("c1_eventnums", &lineparamsc1.eventNum);

	TBranch *branch_ac2 = tree1->Branch("a_c2", &lineparamsc2.a);
	TBranch *branch_aerrc2 = tree1->Branch("aerr_c2", &lineparamsc2.aerr);
	TBranch *branch_bc2 = tree1->Branch("b_c2", &lineparamsc2.b);
	TBranch *branch_berrc2 = tree1->Branch("berr_c2", &lineparamsc2.berr);
	TBranch *branch_chisqc2 = tree1->Branch("c2_chisq", &lineparamsc2.chisq);
	TBranch *branch_evnumc2 = tree1->Branch("c2_eventnums", &lineparamsc2.eventNum);

	TBranch *branch_ac3 = tree1->Branch("a_c3", &lineparamsc3.a);
	TBranch *branch_aerrc3 = tree1->Branch("aerr_c3", &lineparamsc3.aerr);
	TBranch *branch_bc3 = tree1->Branch("b_c3", &lineparamsc3.b);
	TBranch *branch_berrc3 = tree1->Branch("berr_c3", &lineparamsc3.berr);
	TBranch *branch_chisqc3 = tree1->Branch("c3_chisq", &lineparamsc3.chisq);
	TBranch *branch_evnumc3 = tree1->Branch("c3_eventnums", &lineparamsc3.eventNum);

	vector<LineParts> c1_fit_params = fit_single_chamber(0, 0, events, rfuncs, lineparamsc1, branch_ac1, branch_aerrc1, branch_bc1, branch_berrc1, branch_chisqc1, branch_evnumc1);

	xvs.clear();
	yvs.clear();
	rvs.clear();

	vector<LineParts> c2_fit_params = fit_single_chamber(1, 0, events, rfuncs, lineparamsc2, branch_ac2, branch_aerrc2, branch_bc2, branch_berrc2, branch_chisqc2, branch_evnumc2);

	xvs.clear();
	yvs.clear();
	rvs.clear();

	vector<LineParts> c3_fit_params = fit_single_chamber(2, 0, events, rfuncs, lineparamsc3, branch_ac3, branch_aerrc3, branch_bc3, branch_berrc3, branch_chisqc3, branch_evnumc3);

	xvs.clear();
	yvs.clear();
	rvs.clear();

	vector<vector<LineParts>> chamber_fit_params = {c1_fit_params, c2_fit_params, c3_fit_params};

	// calibration stuff
	// not each chamber always has hits so they have diff amount of fits - may need to associate an event num with each fit

	float c2mc1d;
	vector<float> c2mc1ds;
	float c2mc3d;
	vector<float> c2mc3ds;
	TBranch *branch_bdiff_c2_c1 = tree1->Branch("cham2_minus_chamb1_slopes", &c2mc1d);
	TBranch *branch_bdiff_c2_c3 = tree1->Branch("cham2_minus_chamb3_slopes", &c2mc3d);

	/*
		for (int event = 0; event < events.size(); event++)
		{
			if (chamber_fit_params.at(1).at(event).a != 0 && chamber_fit_params.at(0).at(event).a != 0)
			{
				c2mc1d = (chamber_fit_params.at(1).at(event).a) - (chamber_fit_params.at(0).at(event).a);
				c2mc1ds.push_back(c2mc1d);
				branch_bdiff_c2_c1->Fill();
			}

			if ((chamber_fit_params.at(2).at(event).a) != 0 && (chamber_fit_params.at(1).at(event).a) != 0)
			{
				c2mc3d = (chamber_fit_params.at(1).at(event).a) - (chamber_fit_params.at(2).at(event).a);
				c2mc3ds.push_back(c2mc3d);
				branch_bdiff_c2_c3->Fill();
			}
		}
	*/

	for (int event = 0; event < events.size(); event++)
	{
		if (chamber_fit_params.at(1).at(event).a != 0 && chamber_fit_params.at(0).at(event).a != 0)
		{
			c2mc1d = (chamber_fit_params.at(1).at(event).a) - (chamber_fit_params.at(0).at(event).a);

			// Get angular difference between the two lines
			// c2mc1d = atan((chamber_fit_params.at(1).at(event).a) - (chamber_fit_params.at(0).at(event).a) / (1 + (chamber_fit_params.at(1).at(event).a) * (chamber_fit_params.at(0).at(event).a)));
			c2mc1ds.push_back(c2mc1d);
			branch_bdiff_c2_c1->Fill();
		}

		if ((chamber_fit_params.at(2).at(event).a) != 0 && (chamber_fit_params.at(1).at(event).a) != 0)
		{
			c2mc3d = (chamber_fit_params.at(1).at(event).a) - (chamber_fit_params.at(2).at(event).a);

			// Get angular difference between the two lines
			// c2mc3d = atan((chamber_fit_params.at(1).at(event).a) - (chamber_fit_params.at(2).at(event).a) / (1 + (chamber_fit_params.at(1).at(event).a) * (chamber_fit_params.at(2).at(event).a)));

			c2mc3ds.push_back(c2mc3d);
			branch_bdiff_c2_c3->Fill();
		}
	}

	// fit histograms of these diffs to get the middle

	TH1F *twominus1 = new TH1F("", "h1 title", 200, -10.0, 10.0);
	TH1F *twominus3 = new TH1F("", "h1 title", 200, -10.0, 10.0);

	for (int i = 0; i < c2mc1ds.size(); i++)
	{
		if ((chamber_fit_params.at(1).at(i).a) != 0) //&& (c2mc1ds.at(i) < 10. && c2mc1ds.at(i) > -10.))
		{
			twominus1->Fill(c2mc1ds.at(i));
		}
	}
	for (int i = 0; i < c2mc3ds.size(); i++)
	{
		if ((chamber_fit_params.at(1).at(i).a) != 0) //&& (c2mc3ds.at(i) < 10 && c2mc3ds.at(i) > -10))
		{
			twominus3->Fill(c2mc3ds.at(i));
		}
	}

	twominus1->Fit("gaus");
	twominus3->Fit("gaus");

	// Get the parameter 1 of the fit
	float meanc1_a_diffs = twominus1->GetFunction("gaus")->GetParameter(1);
	float meanc3_a_diffs = twominus3->GetFunction("gaus")->GetParameter(1);

	std::cout << meanc1_a_diffs << " " << meanc3_a_diffs << " " << std::endl;

	LineParts c1_f_p;
	LineParts c3_f_p;

	TBranch *branch_bpc1 = tree1->Branch("b'_c1", &c1_f_p.b);
	TBranch *branch_bpec1 = tree1->Branch("b'_err_c1", &c1_f_p.berr);

	TBranch *branch_bpc3 = tree1->Branch("b'_c3", &c3_f_p.b);
	TBranch *branch_bpec3 = tree1->Branch("b'_err_c3", &c3_f_p.berr);

	// TODO: Should we be taking the average of the differences between the slopes and subtracting or finding the
	// average angle between the lines and then rotating the points themselves when fitting

	// c1
	what_entrynum = 0;
	for (int i = 0; i < events.size(); i++)
	{

		as.push_back(chamber_fit_params.at(0).at(i).a - meanc1_a_diffs); // This is a global variable defined in line_fitting.h TODO: remove

		// Get the angle of the line and subtract the mean difference. Convert back to slope
		// as.push_back(tan(atan(chamber_fit_params.at(0).at(i).a) - meanc1_a_diffs)); // This is a global variable defined in line_fitting.h
	}

	cout << "This is the slow bit" << endl;
	vector<LineParts> c1bprime = fit_single_chamber(0, 1, events, rfuncs, c1_f_p, branch_bpc1, branch_bpec1);
	as.clear();

	// c3
	what_entrynum = 0;
	for (int i = 0; i < events.size(); i++)
	{
		as.push_back(chamber_fit_params.at(2).at(i).a - meanc1_a_diffs); // This is a global variable defined in line_fitting.h TODO: remove

		// Get the angle of the line and subtract the mean difference. Convert back to slope
		// as.push_back(tan(atan(chamber_fit_params.at(2).at(i).a) - meanc3_a_diffs)); // This is a global variable defined in line_fitting.h
	}

	vector<LineParts> c3bprime = fit_single_chamber(2, 1, events, rfuncs, c3_f_p, branch_bpc3, branch_bpec3);
	as.clear();

	what_entrynum = 0;

	// make histograms of the b
	TH1F *bprime_chamber1 = new TH1F("", "h1 title", 200, -100., 100.);
	TH1F *bprime_chamber3 = new TH1F("", "h1 title", 200, -100., 100.);

	for (int i = 0; i < events.size(); i++)
	{
		if (c1bprime.at(i).b != 0)
			bprime_chamber1->Fill(c1bprime.at(i).b);
		if (c3bprime.at(i).b != 0)
			bprime_chamber3->Fill(c3bprime.at(i).b);
	}

	bprime_chamber1->Fit("gaus");
	bprime_chamber3->Fit("gaus");

	// Get the parameter 1 of the fit
	float meanc1_b_diffs = twominus1->GetFunction("gaus")->GetParameter(1);
	float meanc3_b_diffs = twominus3->GetFunction("gaus")->GetParameter(1);

	TBranch *branch_bpfc1 = tree1->Branch("fitted_b'_c1", &meanc1_b_diffs);
	TBranch *branch_bpfc3 = tree1->Branch("fitted_b'_c3", &meanc3_b_diffs);

	branch_bpfc1->Fill();
	branch_bpfc3->Fill();

	// fit all chambers together

	std::cout << "fitting all chambers together " << std::endl;

	LineParts lineparams;

	TBranch *branch_a = tree1->Branch("a", &lineparams.a);
	TBranch *branch_aerr = tree1->Branch("aerr", &lineparams.aerr);
	TBranch *branch_b = tree1->Branch("b", &lineparams.b);
	TBranch *branch_berr = tree1->Branch("berr", &lineparams.berr);
	TBranch *branch_chisq = tree1->Branch("chisq", &lineparams.chisq);

	vector<LineParts> fittedlines = fit_chamber(events, rfuncs, lineparams, branch_a, branch_aerr, branch_b, branch_berr, branch_chisq, meanc1_b_diffs, meanc3_b_diffs);

	// Checking all the tubes the line intersects with for efficiency calculations

	vector<int> dist_lessthan_counter;
	vector<int> tube_hit_counter;

	for (int i = 0; i < 144; i++)
	{
		dist_lessthan_counter.push_back(0);
		tube_hit_counter.push_back(0);
	}

	for (int i = 0; i < events.size(); i++)
	{
		LineParts line = fittedlines.at(i);

		vector<int> hittubenums;
		for (int j = 0; j < events.at(i).t.size(); j++)
		{
			int chambnum = events.at(i).chamber.at(j);
			int layernum = events.at(i).layer.at(j);
			int tubenum = events.at(i).tube.at(j);

			int realtubenum = 48 * (chambnum) + 16 * (layernum) + tubenum;

			hittubenums.push_back(realtubenum);
		}

		for (int chamb = 0; chamb < 3; chamb++)
			for (int layer = 0; layer < 3; layer++)
				for (int tube = 0; tube < 16; tube++)
				{
					int tubenum = chamb * 48 + layer * 16 + tube;
					vector<float> xy = getTubeCoords(chamb, layer, tube);
					double d_line_tube = abs(line.a * xy.at(0) + xy.at(1) + line.b) / sqrt(line.a * line.a + 1);
					if (d_line_tube < 1.5) // tube radius
					{
						dist_lessthan_counter.at(tubenum) += 1;

						if (std::find(hittubenums.begin(), hittubenums.end(), tubenum) != hittubenums.end())
							tube_hit_counter.at(tubenum) += 1;
					}
				}
	}

	tree1->SetEntries(events.size());

	// Write the tree to the file and close the file to save.
	tree1->Write();
	file->Close();

	std::cout << "done" << std::endl;

	// app.Run();

	return 0;
}