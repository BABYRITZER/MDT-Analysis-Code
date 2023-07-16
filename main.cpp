#include "Event.h"
#include "line_fitting.h"
#include "load_from_root_file.h"
#include "fitt0s.h"
#include "gransac_implementations.h"
#include "efficiency_calc.h"
#include <TApplication.h>
#include <TStyle.h>
#include <utility>
#include <tuple>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TF1.h>
#include <iostream>
#include <fstream>
#include <vector>

using std::cout;
using std::endl;
using std::get;
using std::string;
using std::tuple;
using std::vector;

vector<float> xvs;
vector<float> yvs;
vector<float> rvs;
vector<float> sigmas;
vector<float> extern_t0s;
vector<float> as;
int what_entrynum;
vector<int> chlist;

// int main(int argc, char **argv)
int main()
{

	// Create the folder "Output" if it doesn't exist. We are going to put all images and plots in it. Only works on Linux and Mac.
	system("mkdir -p Output");

	gStyle->SetOptFit();

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
	vector<NewEvent> events_plus25 = events;
	vector<NewEvent> events_minus25 = events;

	for (int q = 0; q < events.size(); ++q)
	{
		int hitsize = events.at(q).t.size();
		for (int w = 0; w < hitsize; ++w)
		{
			events_plus25.at(q).t.at(w) = events_plus25.at(q).t.at(w) + (25. / sqrt(12.));
			events_minus25.at(q).t.at(w) = events_minus25.at(q).t.at(w) + (25. / sqrt(12.));
		}
	}

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

	/*
		vector<TF1> rfuncs_p25;
		vector<TF1 *> rfuncs_p25st;
		vector<TF1> rfuncs_m25;
		vector<TF1 *> rfuncs_m25st;

		double cc;
		double dd;
		TF1 ee;

		rfuncs_p25 = populate_t0s(events_plus25, cc, dd, ee, brancht0, brancht0_err, brancht0_fits);
		rfuncs_m25 = populate_t0s(events_minus25, cc, dd, ee, brancht0, brancht0_err, brancht0_fits);
	*/
	for (int q = 0; q < 144; q++)
	{
		rfuncsst.push_back(&rfuncs.at(q));

		// rfuncs_p25st.push_back(&rfuncs_p25.at(q));
		// rfuncs_m25st.push_back(&rfuncs_p25.at(q));
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	branchr_funcs->Fill();

	std::cout << "finding outlier points" << std::endl;
	// get number of events and add them to tree
	long int NEvents = events.size();

	std::cout << NEvents << " events have been loaded " << std::endl;
	// pattern recog here
	// populates the tree with new entries with pattern recognized hits

	NewEvent event;

	int nhits;

	TBranch *branch_eventnum = tree1->Branch("event_num", &event.eventnum);
	TBranch *branch_t = tree1->Branch("time", &event.t);
	TBranch *branch_chg = tree1->Branch("charge", &event.charge);
	TBranch *branch_chmb = tree1->Branch("chamber", &event.chamber);
	TBranch *branch_layer = tree1->Branch("layer", &event.layer);
	TBranch *branch_tube = tree1->Branch("tube", &event.tube);
	TBranch *branch_is_inliner = tree1->Branch("is_inliner", &event.is_inlier);
	TBranch *branch_nhits = tree1->Branch("nhits", &nhits);

	vector<float> gransac_lineparams = findInliers(events, event, branch_eventnum, branch_t, branch_chg, branch_chmb, branch_layer, branch_tube, branch_is_inliner);

	// vector<float> gransac_lineparams ;
	// fill in the nhits branch
	for (int a = 0; a < events.size(); ++a)
	{
		gransac_lineparams.push_back(0);
		gransac_lineparams.push_back(0);

		nhits = events.at(a).t.size();
		branch_nhits->Fill();
	}

	xvs.clear();
	yvs.clear();
	rvs.clear();

	/*

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
	//TODO:


		float mean = fittwochambers(0, 2, 1, events, rfuncs);
		mean = fittwochambers(1, 2, 0, events, rfuncs);
		mean = fittwochambers(0, 1, 2, events, rfuncs);

		return 0;




		vector<LineParts> c1_fit_params = fit_single_chamber(1, 0, events, rfuncs, lineparamsc1, branch_ac1, branch_aerrc1, branch_bc1, branch_berrc1, branch_chisqc1, branch_evnumc1, gransac_lineparams);

		xvs.clear();
		yvs.clear();
		rvs.clear();
		sigmas.clear();

		vector<LineParts> c2_fit_params = fit_single_chamber(1, 0, events, rfuncs, lineparamsc2, branch_ac2, branch_aerrc2, branch_bc2, branch_berrc2, branch_chisqc2, branch_evnumc2, gransac_lineparams);

		xvs.clear();
		yvs.clear();
		rvs.clear();
		sigmas.clear();

		vector<LineParts> c3_fit_params = fit_single_chamber(2, 0, events, rfuncs, lineparamsc3, branch_ac3, branch_aerrc3, branch_bc3, branch_berrc3, branch_chisqc3, branch_evnumc3, gransac_lineparams);

		xvs.clear();
		yvs.clear();
		rvs.clear();
		sigmas.clear();

		vector<vector<LineParts>> chamber_fit_params = {c1_fit_params, c2_fit_params, c3_fit_params};

		// calibration stuff

		float c2mc1d;
		float c2mc3d;

		TBranch *branch_bdiff_c2_c1 = tree1->Branch("cham2_minus_chamb1_slopes", &c2mc1d);
		TBranch *branch_bdiff_c2_c3 = tree1->Branch("cham2_minus_chamb3_slopes", &c2mc3d);

		TH1F *twominus1 = new TH1F("SoC2MC1", "Slope of C2 Minus C1", 200, -2 * 3.14, 2 * 3.14);
		TH1F *twominus3 = new TH1F("SoC2MC3", "Slope of C2 Minus C3", 200, -2 * 3.14, 2 * 3.14);

		TH1F *c1intercepts = new TH1F("c1int", "C2-C1 Unadjusted Fit Intercepts", 200, -50, 50);
		TH1F *c3intercepts = new TH1F("c3int", "C2-C3 Unadjusted Fit Intercepts", 200, -50, 50);

		// TAKING ARCTANS HERE
		for (int i = 0; i < chamber_fit_params.at(1).size(); ++i)
		{
			// to ignore if fit failed or if the fn didnt fit cuz not enough hits
			if (chamber_fit_params.at(1).at(i).a != 0 && chamber_fit_params.at(0).at(i).a != 0)
			{
				c2mc1d = atan(chamber_fit_params.at(1).at(i).a) - atan(chamber_fit_params.at(0).at(i).a);
				twominus1->Fill(c2mc1d);
				branch_bdiff_c2_c1->Fill();
			}
		}

		for (int i = 0; i < chamber_fit_params.at(1).size(); ++i)
		{
			// to ignore if fit failed or if the fn didnt fit cuz not enough hits
			if (chamber_fit_params.at(1).at(i).a != 0 && chamber_fit_params.at(2).at(i).a != 0)
			{
				c2mc3d = atan(chamber_fit_params.at(1).at(i).a) - atan(chamber_fit_params.at(2).at(i).a);
				twominus3->Fill(c2mc3d);
				branch_bdiff_c2_c3->Fill();
			}
		}

		for (int i = 0; i < chamber_fit_params.at(1).size(); ++i)
		{
			if (chamber_fit_params.at(1).at(i).b != 0 && chamber_fit_params.at(0).at(i).b != 0)
				c1intercepts->Fill(chamber_fit_params.at(1).at(i).b - chamber_fit_params.at(0).at(i).b);

			if (chamber_fit_params.at(1).at(i).b != 0 && chamber_fit_params.at(2).at(i).b != 0)
				c3intercepts->Fill(chamber_fit_params.at(1).at(i).b - chamber_fit_params.at(2).at(i).b);
		}

		TCanvas *c1int = new TCanvas();

		c1intercepts->GetXaxis()->SetTitle("Fitted Intercepts of Chamber 2 - Chamber 1 Before Angle Adjustment (cm)");
		c1intercepts->GetYaxis()->SetTitle("Counts");
		c1intercepts->Draw("");
		c1int->Modified();
		c1int->Update();

		c1int->SaveAs("./Output/c1intercepts_noangle.png");
		c1int->Close();

		TCanvas *c3int = new TCanvas();

		c3intercepts->GetXaxis()->SetTitle("Fitted Intercepts of Chamber 2 - Chamber 1 Before Angle Adjustment (cm)");
		c3intercepts->GetYaxis()->SetTitle("Counts");
		c3intercepts->Draw("");
		c3int->Modified();
		c3int->Update();

		c3int->SaveAs("./Output/c3intercepts_noangle.png");
		c3int->Close();

		// fitting the slope diffs
		twominus1->Fit("gaus");
		twominus3->Fit("gaus");

		TCanvas *cc1 = new TCanvas();

		cc1->Divide(1, 2);

		cc1->cd(1);

		//twominus1->SetBit(TH1::kNoStats);

		twominus1->GetFunction("gaus")->Draw();
		twominus1->GetXaxis()->SetTitle("Chamber 2 and 1 slope difference (rad)");
		twominus1->GetYaxis()->SetTitle("Counts");
		twominus1->Draw("");
		cc1->Modified();
		cc1->Update();

		cc1->cd(2);

		double xvals[200];
		double yvals[200];

		for (int ibin = 0; ibin <= 200; ++ibin)
		{
			double deez = twominus1->GetBinContent(ibin);
			// TODO: check that this is correct
			double xm = -(2 * 3.14) + (4 * 3.14 * (double)ibin / 200.);

			double res = (deez - twominus1->GetFunction("gaus")->Eval(twominus1->GetBinCenter(ibin)));

			xvals[ibin] = xm;
			yvals[ibin] = res;
		}

		TGraph *graph = new TGraph(200, xvals, yvals);

		graph->SetTitle("Fit Residuals");
		graph->Draw("A*");
		graph->SetMarkerStyle(3);
		cc1->Modified();
		cc1->Update();

		cc1->cd(0);
		cc1->SaveAs("./Output/slope_diffs_chamber1.png");

		// end c2-c1 stuff

		TCanvas *c2 = new TCanvas();

		c2->Divide(1, 2);

		c2->cd(1);

		// twominus3->SetBit(TH1::kNoStats);

		twominus3->GetFunction("gaus")->Draw();
		twominus3->GetXaxis()->SetTitle("Chamber 2 and 3 slope difference (rad)");
		twominus3->GetYaxis()->SetTitle("Counts");
		twominus3->Draw("");
		c2->Modified();
		c2->Update();

		c2->cd(2);

		// Get residuals
		double xvals3[200];
		double yvals3[200];

		for (int ibin = 0; ibin <= 200; ++ibin)
		{
			double deez = twominus3->GetBinContent(ibin);
			// TODO::check this is ok
			double xm = -(2 * 3.14) + (4 * 3.14 * (double)ibin / 200.);

			double res = (deez - twominus3->GetFunction("gaus")->Eval(twominus3->GetBinCenter(ibin)));

			xvals3[ibin] = xm;
			yvals3[ibin] = res;
		}

		TGraph *graph3 = new TGraph(200, xvals3, yvals3);

		graph3->SetTitle("Fit Residuals");
		graph3->GetYaxis()->SetTitle("Residuals");
		graph3->GetXaxis()->SetTitle("Chamber 2 and 3 slope difference");
		graph3->Draw("A*");
		graph3->SetMarkerStyle(3);

		c2->Modified();
		c2->Update();
		// end c2 - c3 stuff

		c2->cd(0);
		c2->SaveAs("./Output/slope_diffs_chamber3.png");

		// Get the parameter 1 of the fit
		float meanc1_a_diffs = twominus1->GetFunction("gaus")->GetParameter(1);
		float meanc3_a_diffs = twominus3->GetFunction("gaus")->GetParameter(1);

		TBranch *branch_afc1 = tree1->Branch("fitted_a_diff_c1", &meanc1_a_diffs);
		TBranch *branch_afc3 = tree1->Branch("fitted_a_diff_c3", &meanc3_a_diffs);

		branch_afc1->Fill();
		branch_afc3->Fill();

		std::cout << meanc1_a_diffs << " " << meanc3_a_diffs << " " << std::endl;

		// start doing the b' stuff
		LineParts c1_f_p;
		LineParts c3_f_p;

		TBranch *branch_bpc1 = tree1->Branch("b'_c1", &c1_f_p.b);
		TBranch *branch_bpec1 = tree1->Branch("b'_err_c1", &c1_f_p.berr);

		TBranch *branch_bpc3 = tree1->Branch("b'_c3", &c3_f_p.b);
		TBranch *branch_bpec3 = tree1->Branch("b'_err_c3", &c3_f_p.berr);

		// c1
		vector<LineParts> c1bprimefits = fit_single_chamber(0, 1, meanc1_a_diffs, events, rfuncs, c1_f_p, branch_bpc1, branch_bpec1);

		// c3
		vector<LineParts> c3bprimefits = fit_single_chamber(2, 1, meanc3_a_diffs, events, rfuncs, c3_f_p, branch_bpc3, branch_bpec3);

		// make histograms of the b
		TH1F *bprime_chamber1 = new TH1F("", "h1 title", 200, -3, 3);
		TH1F *bprime_chamber3 = new TH1F("", "h1 title", 200, -8, 8);

		// original settings for these (bad fits) were new TH1F("", "h1 title", 100, -60., 60.);
		// after changed range to -2,2

		for (int i = 0; i < events.size(); i++)
		{
			if (c1bprimefits.at(i).b != 0)
				bprime_chamber1->Fill(chamber_fit_params.at(1).at(i).b - c1bprimefits.at(i).b);
			if (c3bprimefits.at(i).b != 0)
				bprime_chamber3->Fill(chamber_fit_params.at(1).at(i).b - c3bprimefits.at(i).b);
		}

		bprime_chamber1->Fit("gaus");
		bprime_chamber3->Fit("gaus");

		TCanvas *c3 = new TCanvas("c3", "c3");
		TCanvas *c4 = new TCanvas("c4", "c4");

		c3->Divide(1, 2);
		c4->Divide(1, 2);

		c3->cd(1);

		//bprime_chamber1->SetBit(TH1::kNoStats);

		bprime_chamber1->SetTitle("Fitted Line Intercepts Difference wrt Chamber 2 after Rotation for Chamber 1");
		bprime_chamber1->GetXaxis()->SetTitle("Line Intercept (cm)");
		bprime_chamber1->GetYaxis()->SetTitle("Counts");
		bprime_chamber1->Draw("");
		bprime_chamber1->GetFunction("gaus")->Draw("same");

		c3->Modified();
		c3->Update();

		c3->cd(2);

		// Get residuals
		double xvals4[200];
		double yvals4[200];

		for (int ibin = 0; ibin <= 200; ++ibin)
		{
			double deez = bprime_chamber1->GetBinContent(ibin);

			double xm = bprime_chamber1->GetBinCenter(ibin);

			double res = (deez - bprime_chamber1->GetFunction("gaus")->Eval(bprime_chamber1->GetBinCenter(ibin)));

			xvals4[ibin] = xm;
			yvals4[ibin] = res;
		}

		TGraph *graph_bprime_chamber_1 = new TGraph(200, xvals4, yvals4);

		graph_bprime_chamber_1->SetTitle("Fit Residuals");
		graph_bprime_chamber_1->GetXaxis()->SetTitle("Line Intercept (cm)");
		graph_bprime_chamber_1->GetYaxis()->SetTitle("Residuals");
		graph_bprime_chamber_1->Draw("A*");

		graph_bprime_chamber_1->SetMarkerStyle(3);

		c3->Modified();
		c3->Update();

		c4->cd(1);

		//bprime_chamber3->SetBit(TH1::kNoStats);

		bprime_chamber3->SetTitle("Fitted Line Intercepts Difference wrt Chamber 2 After Rotation for Chamber 3");
		bprime_chamber3->GetXaxis()->SetTitle("Line Intercept (cm)");
		bprime_chamber3->GetYaxis()->SetTitle("Counts");
		bprime_chamber3->Draw();
		bprime_chamber3->GetFunction("gaus")->Draw("same");

		c4->Modified();

		c4->cd(2);

		// Get residuals

		double xvals5[200];
		double yvals5[200];

		for (int ibin = 0; ibin <= 200; ++ibin)
		{
			double deez = bprime_chamber3->GetBinContent(ibin);

			double xm = bprime_chamber3->GetBinCenter(ibin);

			double res = (deez - bprime_chamber3->GetFunction("gaus")->Eval(bprime_chamber3->GetBinCenter(ibin)));

			xvals5[ibin] = xm;
			yvals5[ibin] = res;
		}

		TGraph *graph_bprime_chamber_3 = new TGraph(200, xvals5, yvals5);

		graph_bprime_chamber_3->SetTitle("Fit Residuals");
		graph_bprime_chamber_3->GetXaxis()->SetTitle("Line Intercept (cm)");
		graph_bprime_chamber_3->GetYaxis()->SetTitle("Residuals");
		graph_bprime_chamber_3->Draw("A*");
		graph_bprime_chamber_3->SetMarkerStyle(3);

		c4->Modified();
		c4->Update();

		c3->SaveAs("./Output/bprime_chamber1.png");
		c4->SaveAs("./Output/bprime_chamber3.png");

		// Get the parameter 1 of the fit
		float meanc1_b_diffs = bprime_chamber1->GetFunction("gaus")->GetParameter(1);
		float meanc3_b_diffs = bprime_chamber3->GetFunction("gaus")->GetParameter(1);

		TBranch *branch_bpfc1 = tree1->Branch("fitted_b'_c1", &meanc1_b_diffs);
		TBranch *branch_bpfc3 = tree1->Branch("fitted_b'_c3", &meanc3_b_diffs);

		branch_bpfc1->Fill();
		branch_bpfc3->Fill();

		// finally, fit all chambers together
	*/

	std::cout << "~~~~~fitting all chambers together~~~~~" << std::endl;

	LineParts lineparams;
	LineParts linepar_dummy;

	TBranch *branch_a = tree1->Branch("a", &lineparams.a);
	TBranch *branch_aerr = tree1->Branch("aerr", &lineparams.aerr);
	TBranch *branch_b = tree1->Branch("b", &lineparams.b);
	TBranch *branch_berr = tree1->Branch("berr", &lineparams.berr);
	TBranch *branch_chisq = tree1->Branch("chisq", &lineparams.chisq);

	TBranch *branch_c1 = tree1->Branch("c1_offset", &lineparams.ch1);
	TBranch *branch_c1err = tree1->Branch("c1_offseterr", &lineparams.ch1err);
	TBranch *branch_c3 = tree1->Branch("c3_offset", &lineparams.ch3);
	TBranch *branch_c3err = tree1->Branch("c3_offseterr", &lineparams.ch3err);

	float meanc1_angle = 0.;
	float meanc3_angle = 0.;

	float meanc1_b_diffs = 0.;
	float meanc3_b_diffs = 0.;

	vector<LineParts> fittedlines = fit_chamber(events, rfuncs, lineparams, branch_a, branch_aerr, branch_b, branch_berr,
												branch_chisq, meanc1_angle, meanc3_angle, meanc1_b_diffs, meanc3_b_diffs, gransac_lineparams,
												branch_c1, branch_c1err, branch_c3, branch_c3err);

	vector<LineParts> fittedlines_plus25 = fit_chamber(events_plus25, rfuncs, linepar_dummy, branch_a, branch_aerr, branch_b, branch_berr,
													   branch_chisq, meanc1_angle, meanc3_angle, meanc1_b_diffs, meanc3_b_diffs, gransac_lineparams, branch_c1, branch_c1err, branch_c3, branch_c3err);

	vector<LineParts> fittedlines_minus25 = fit_chamber(events_minus25, rfuncs, linepar_dummy, branch_a, branch_aerr, branch_b, branch_berr,
														branch_chisq, meanc1_angle, meanc3_angle, meanc1_b_diffs, meanc3_b_diffs, gransac_lineparams,
														branch_c1, branch_c1err, branch_c3, branch_c3err);

	vector<LineParts> fittedlines_systematicerror_plus;
	vector<LineParts> fittedlines_systematicerror_minus;

	for (int e = 0; e < fittedlines.size(); ++e)
	{
		fittedlines_systematicerror_plus.push_back(subtractlineparts(fittedlines_plus25.at(e), fittedlines.at(e)));
		fittedlines_systematicerror_minus.push_back(subtractlineparts(fittedlines.at(e), fittedlines_minus25.at(e)));
	}

	// now calculate per tube efficiency

	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~efficiency calculations~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

	int eff_top;
	int eff_bot;

	int dummy;

	TBranch *branch_top = tree1->Branch("eff_top", &eff_top);
	TBranch *branch_bot = tree1->Branch("eff_bot", &eff_bot);

	auto effs_errors_tuple = eff_calc(events, fittedlines, rfuncs, branch_top, branch_bot, eff_top, eff_bot);
	auto effs_errors_tuple_plus25 = eff_calc(events_plus25, fittedlines, rfuncs, branch_top, branch_bot, dummy, dummy);
	auto effs_errors_tuple_minus25 = eff_calc(events_minus25, fittedlines, rfuncs, branch_top, branch_bot, dummy, dummy);

	vector<float> efficiency_systematicerror_plus;
	vector<float> efficiency_systematicerror_minus;

	// Getting efficiencies from the tuples (looks weird but this is the way to do it)
	vector<float> effs = get<0>(effs_errors_tuple);
	vector<float> effs_plus25 = get<0>(effs_errors_tuple_plus25);
	vector<float> effs_minus25 = get<0>(effs_errors_tuple_minus25);

	// Getting errors from the tuples (1 is +, 2 is -) only for effs
	vector<float> effs_errors_plus = get<1>(effs_errors_tuple);
	vector<float> effs_errors_minus = get<2>(effs_errors_tuple);

	for (int e = 0; e < effs.size(); ++e)
	{
		efficiency_systematicerror_plus.push_back(effs_plus25.at(e) - effs.at(e));
		efficiency_systematicerror_minus.push_back(effs.at(e) - effs_minus25.at(e));
	}

	/* Write everything to branches in the trees:
	 * 1. fitted line parameters systematic errors (a, b)
	 * 2. efficiencies (values, errors, systematic errors)
	 */

	// 1. fitted line parameters systematic errors (a, b)

	double a_sys_plus;
	double b_sys_plus;
	double a_sys_minus;
	double b_sys_minus;

	TBranch *branch_a_sys_plus = tree1->Branch("a_sys_plus", &a_sys_plus);
	TBranch *branch_b_sys_plus = tree1->Branch("b_sys_plus", &b_sys_plus);
	TBranch *branch_a_sys_minus = tree1->Branch("a_sys_minus", &a_sys_minus);
	TBranch *branch_b_sys_minus = tree1->Branch("b_sys_minus", &b_sys_minus);

	for (int g = 0; g < fittedlines.size(); ++g)
	{
		a_sys_plus = fittedlines_systematicerror_plus.at(g).a;
		branch_a_sys_plus->Fill();

		b_sys_plus = fittedlines_systematicerror_plus.at(g).b;
		branch_b_sys_plus->Fill();

		a_sys_minus = fittedlines_systematicerror_minus.at(g).a;
		branch_a_sys_minus->Fill();

		b_sys_minus = fittedlines_systematicerror_minus.at(g).b;
		branch_b_sys_minus->Fill();
	}

	// 2. efficiencies (values, errors, systematic errors)

	double eff;
	double eff_err_plus;
	double eff_err_minus;

	double eff_sys_plus;
	double eff_sys_minus;

	TBranch *branch_eff = tree1->Branch("eff", &eff);
	TBranch *branch_eff_err_plus = tree1->Branch("eff_err_plus", &eff_err_plus);
	TBranch *branch_eff_err_minus = tree1->Branch("eff_err_minus", &eff_err_minus);

	TBranch *branch_eff_sys_plus = tree1->Branch("eff_sys_plus", &eff_sys_plus);
	TBranch *branch_eff_sys_minus = tree1->Branch("eff_sys_minus", &eff_sys_minus);

	for (int g = 0; g < effs.size(); ++g)
	{
		// Values
		eff = effs.at(g);
		branch_eff->Fill();

		// Errors
		eff_err_plus = effs_errors_plus.at(g);
		branch_eff_err_plus->Fill();

		eff_err_minus = effs_errors_minus.at(g);
		branch_eff_err_minus->Fill();

		// Systematic errors
		eff_sys_plus = efficiency_systematicerror_plus.at(g);
		branch_eff_sys_plus->Fill();

		eff_sys_minus = efficiency_systematicerror_minus.at(g);
		branch_eff_sys_minus->Fill();
	}

	auto layereffs = layer_effcalc(events, fittedlines, rfuncs);

	float layereff_top;
	float layereff_bot;

	float layereff_top_sysplus;

	float layereff_top_sysminus;

	TBranch *branch_layereff_top = tree1->Branch("layer_eff_top", &layereff_top);
	TBranch *branch_layereff_bot = tree1->Branch("layer_eff_bot", &layereff_bot);

	// systematic error comes from the MIGRAD error in the fit parameters
	TBranch *branch_layereff_top_sys_plus = tree1->Branch("layer_eff_top_sys_plus", &layereff_top_sysplus);
	TBranch *branch_layereff_top_sys_minus = tree1->Branch("layer_eff_top_sys_minus", &layereff_top_sysminus);

	for (int w = 0; w < 9; ++w)
	{
		//    return std::make_tuple(topnumbers, bottomnumbers, topnumssysplus, bottomnumssysplus, topnumssysminus, bottomnumssysminus);
		// 0 1 2 4
		layereff_top = get<0>(layereffs).at(w);
		layereff_bot = get<1>(layereffs).at(w);

		layereff_top_sysplus = get<2>(layereffs).at(w);
		layereff_top_sysminus = get<4>(layereffs).at(w);

		branch_layereff_top->Fill();
		branch_layereff_bot->Fill();

		branch_layereff_top_sys_plus->Fill();
		branch_layereff_top_sys_minus->Fill();
	}

	tree1->SetEntries(events.size());

	// Write the tree to the file and close the file to save.
	tree1->Write();
	file->Close();

	std::cout << "done" << std::endl;

	// app.Run();

	return 0;
}
