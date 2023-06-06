#include "Event.h"
#include "line_fitting.h"
#include "load_from_root_file.h"
#include "fitt0s.h"
#include "gransac_implementations.h"

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TF1.h>
#include <iostream>
#include <fstream>
#include <vector>

using std::string;
using std::vector;

vector<float> xvs;
vector<float> yvs;
vector<float> rvs;

int main()
{
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

	// TBranch *branch0 = tree1->Branch("Events", &event);

	// TBranch *branch0 = tree1->Branch("Events", &event, "eventnum/I:t/F:charge/F:chamber/I:layer/I:tube/I");
	// addeventstotree(events, event, branch0);

	std::cout << NEvents << " events have been loaded " << std::endl;
	// pattern recog here
	// populates the tree with new entries with pattern recognized hits
	//	TBranch *branch1 = tree1->Branch("Processed_Events", &processed_events);

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

	TBranch *branch_ac2 = tree1->Branch("a_c2", &lineparamsc2.a);
	TBranch *branch_aerrc2 = tree1->Branch("aerr_c2", &lineparamsc2.aerr);
	TBranch *branch_bc2 = tree1->Branch("b_c2", &lineparamsc2.b);
	TBranch *branch_berrc2 = tree1->Branch("berr_c2", &lineparamsc2.berr);
	TBranch *branch_chisqc2 = tree1->Branch("c2_chisq", &lineparamsc2.chisq);

	TBranch *branch_ac3 = tree1->Branch("a_c3", &lineparamsc1.a);
	TBranch *branch_aerrc3 = tree1->Branch("aerr_c3", &lineparamsc1.aerr);
	TBranch *branch_bc3 = tree1->Branch("b_c3", &lineparamsc1.b);
	TBranch *branch_berrc3 = tree1->Branch("berr_c3", &lineparamsc1.berr);
	TBranch *branch_chisqc3 = tree1->Branch("c3_chisq", &lineparamsc1.chisq);

	vector<vector<LineParts>> chamber_fit_params = fit_single_chambers(events, rfuncs, lineparamsc1, lineparamsc2, lineparamsc3,
																	   branch_ac1, branch_aerrc1, branch_bc1, branch_berrc1, branch_chisqc1,
																	   branch_ac2, branch_aerrc2, branch_bc2, branch_berrc2, branch_chisqc2,
																	   branch_ac3, branch_aerrc3, branch_bc3, branch_berrc3, branch_chisqc3);

	xvs.clear();
	yvs.clear();
	rvs.clear();

//	std::cout << chamber_fit_params.size() << " " << chamber_fit_params.at(0).size() << " " << chamber_fit_params.at(1).size() << " " << chamber_fit_params.at(2).size() << std::endl;
//	return 0;

	// calibration stuff
	//not each chamber always has hits so they have diff amount of fits - may need to associate an event num with each fit

	/*
		float c2mc1d;
		float c2mc3d;
		TBranch *branch_bdiff_c2_c1 = tree1->Branch("cham2_minus_chamb1_slopes", &c2mc1d);
		TBranch *branch_bdiff_c2_c3 = tree1->Branch("cham2_minus_chamb3_slopes", &c2mc3d);

		for (int event = 0; event < events.size(); event++)
		{
			c2mc1d = (-1. / chamber_fit_params.at(1).at(event).a) - (-1. / chamber_fit_params.at(0).at(event).a);
			c2mc3d = (-1. / chamber_fit_params.at(1).at(event).a) - (-1. / chamber_fit_params.at(2).at(event).a);
			branch_bdiff_c2_c1->Fill();
			branch_bdiff_c2_c3->Fill();
		}
	*/

	// after taking these slope differences to get the angle offsets relative to chamber two we then have to refit taking into account the angle differences to get a value for the intercept
	// then fit all chambers together after adjusting

	// fit all chambers together
	std::cout << "fitting all chambers together " << std::endl;

	LineParts lineparams;

	TBranch *branch_a = tree1->Branch("a", &lineparams.a);
	TBranch *branch_aerr = tree1->Branch("aerr", &lineparams.aerr);
	TBranch *branch_b = tree1->Branch("b", &lineparams.b);
	TBranch *branch_berr = tree1->Branch("berr", &lineparams.berr);
	TBranch *branch_chisq = tree1->Branch("chisq", &lineparams.chisq);

	fit_chamber(events, rfuncs, lineparams, branch_a, branch_aerr, branch_b, branch_berr, branch_chisq);

	tree1->SetEntries(events.size());
	tree1->Write();
	file->Close();

	std::cout << "done" << std::endl;
	return 0;
}