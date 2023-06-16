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
vector<float> sigmas;
vector<float> extern_t0s;
vector<float> as;
int what_entrynum;

// int main(int argc, char **argv)
int main()
{

	// Create the folder "Output" if it doesn't exist. We are going to put all images and plots in it. Only works on Linux and Mac.
	system("mkdir -p Output");

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
	sigmas.clear();

	vector<LineParts> c2_fit_params = fit_single_chamber(1, 0, events, rfuncs, lineparamsc2, branch_ac2, branch_aerrc2, branch_bc2, branch_berrc2, branch_chisqc2, branch_evnumc2);

	xvs.clear();
	yvs.clear();
	rvs.clear();
	sigmas.clear();

	vector<LineParts> c3_fit_params = fit_single_chamber(2, 0, events, rfuncs, lineparamsc3, branch_ac3, branch_aerrc3, branch_bc3, branch_berrc3, branch_chisqc3, branch_evnumc3);

	xvs.clear();
	yvs.clear();
	rvs.clear();
	sigmas.clear();

	vector<vector<LineParts>> chamber_fit_params = {c1_fit_params, c2_fit_params, c3_fit_params};

	// calibration stuff

	// TODO: the fitting still kinda sucks
	//////////////////////////////////////////////////////
	float c2mc1d;
	float c2mc3d;

	TBranch *branch_bdiff_c2_c1 = tree1->Branch("cham2_minus_chamb1_slopes", &c2mc1d);
	TBranch *branch_bdiff_c2_c3 = tree1->Branch("cham2_minus_chamb3_slopes", &c2mc3d);

	TH1F *twominus1 = new TH1F("SoC2MC1", "Slope of C2 Minus C1", 100, -10, 10);
	TH1F *twominus3 = new TH1F("SoC2MC3", "Slope of C2 Minus C3", 100, -10, 10);

	std::cout << chamber_fit_params.at(0).size() << " " << chamber_fit_params.at(1).size() << " " << chamber_fit_params.at(2).size() << std::endl;

	for (int i = 0; i < chamber_fit_params.at(1).size(); ++i)
	{
		c2mc1d = chamber_fit_params.at(1).at(i).a - chamber_fit_params.at(0).at(i).a;
		if (c2mc1d != 0)
			twominus1->Fill(c2mc1d);
		branch_bdiff_c2_c1->Fill();
	}

	for (int i = 0; i < chamber_fit_params.at(1).size(); ++i)
	{
		c2mc3d = chamber_fit_params.at(1).at(i).a - chamber_fit_params.at(2).at(i).a;
		if (c2mc3d != 0)
			twominus3->Fill(c2mc3d);
		branch_bdiff_c2_c3->Fill();
	}

	twominus1->Fit("gaus");
	twominus3->Fit("gaus");

	TCanvas *cc1 = new TCanvas();

	cc1->Divide(1, 2);

	cc1->cd(1);

	twominus1->SetBit(TH1::kNoStats);

	twominus1->GetFunction("gaus")->Draw();
	twominus1->GetXaxis()->SetTitle("Chamber 2 and 1 slope difference");
	twominus1->GetYaxis()->SetTitle("Counts");
	twominus1->Draw("");
	cc1->Modified();
	cc1->Update();

	cc1->cd(2);

	double xvals[100];
	double yvals[100];

	for (int ibin = 0; ibin <= 100; ++ibin)
	{
		double deez = twominus1->GetBinContent(ibin);

		double xm = -10. + (20. * (double)ibin / 100.);

		double res = (deez - twominus1->GetFunction("gaus")->Eval(twominus1->GetBinCenter(ibin)));

		xvals[ibin] = xm;
		yvals[ibin] = res;
	}

	TGraph *graph = new TGraph(100, xvals, yvals);

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

	twominus3->SetBit(TH1::kNoStats);

	twominus3->GetFunction("gaus")->Draw();
	twominus3->GetXaxis()->SetTitle("Chamber 2 and 3 slope difference");
	twominus3->GetYaxis()->SetTitle("Counts");
	twominus3->Draw("");
	c2->Modified();
	c2->Update();

	c2->cd(2);

	// Get residuals
	double xvals3[100];
	double yvals3[100];

	for (int ibin = 0; ibin <= 100; ++ibin)
	{
		double deez = twominus3->GetBinContent(ibin);

		double xm = -10. + (20. * (double)ibin / 100.);

		double res = (deez - twominus3->GetFunction("gaus")->Eval(twominus3->GetBinCenter(ibin)));

		xvals3[ibin] = xm;
		yvals3[ibin] = res;
	}

	TGraph *graph3 = new TGraph(100, xvals3, yvals3);

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
	vector<LineParts> c1bprimefits = fit_single_chamber(0, 1, atan(meanc1_a_diffs), events, rfuncs, c1_f_p, branch_bpc1, branch_bpec1);

	// c3
	vector<LineParts> c3bprimefits = fit_single_chamber(2, 1, atan(meanc3_a_diffs), events, rfuncs, c3_f_p, branch_bpc3, branch_bpec3);

	// make histograms of the b
	TH1F *bprime_chamber1 = new TH1F("", "h1 title", 100, -2., 2.);
	TH1F *bprime_chamber3 = new TH1F("", "h1 title", 100, -2., 2.);

	//original settings for these (bad fits) were new TH1F("", "h1 title", 100, -60., 60.);

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

	bprime_chamber1->SetBit(TH1::kNoStats);

	bprime_chamber1->SetTitle("Fitted Line Intercepts after Rotation for Chamber 1");
	bprime_chamber1->GetXaxis()->SetTitle("Line Intercept (cm)");
	bprime_chamber1->GetYaxis()->SetTitle("Counts");
	bprime_chamber1->Draw("");
	bprime_chamber1->GetFunction("gaus")->Draw("same");

	c3->Modified();
	c3->Update();

	c3->cd(2);

	// Get residuals
	double xvals4[100];
	double yvals4[100];

	for (int ibin = 0; ibin <= 100; ++ibin)
	{
		double deez = bprime_chamber1->GetBinContent(ibin);

		double xm = bprime_chamber1->GetBinCenter(ibin);

		double res = (deez - bprime_chamber1->GetFunction("gaus")->Eval(bprime_chamber1->GetBinCenter(ibin)));

		xvals4[ibin] = xm;
		yvals4[ibin] = res;
	}

	TGraph *graph_bprime_chamber_1 = new TGraph(100, xvals4, yvals4);

	graph_bprime_chamber_1->SetTitle("Fit Residuals");
	graph_bprime_chamber_1->GetXaxis()->SetTitle("Line Intercept (cm)");
	graph_bprime_chamber_1->GetYaxis()->SetTitle("Residuals");
	graph_bprime_chamber_1->Draw("A*");

	graph_bprime_chamber_1->SetMarkerStyle(3);

	c3->Modified();
	c3->Update();

	c4->cd(1);

	bprime_chamber3->SetBit(TH1::kNoStats);

	bprime_chamber3->SetTitle("Fitted Line Intercepts After Rotation for Chamber 3");
	bprime_chamber3->GetXaxis()->SetTitle("Line Intercept (cm)");
	bprime_chamber3->GetYaxis()->SetTitle("Counts");
	bprime_chamber3->Draw();
	bprime_chamber3->GetFunction("gaus")->Draw("same");

	c4->Modified();

	c4->cd(2);

	// Get residuals

	double xvals5[100];
	double yvals5[100];

	for (int ibin = 0; ibin <= 100; ++ibin)
	{
		double deez = bprime_chamber3->GetBinContent(ibin);

		double xm = bprime_chamber3->GetBinCenter(ibin);

		double res = (deez - bprime_chamber3->GetFunction("gaus")->Eval(bprime_chamber3->GetBinCenter(ibin)));

		xvals5[ibin] = xm;
		yvals5[ibin] = res;
	}

	TGraph *graph_bprime_chamber_3 = new TGraph(100, xvals5, yvals5);

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

	std::cout << "fitting all chambers together " << std::endl;

	LineParts lineparams;

	TBranch *branch_a = tree1->Branch("a", &lineparams.a);
	TBranch *branch_aerr = tree1->Branch("aerr", &lineparams.aerr);
	TBranch *branch_b = tree1->Branch("b", &lineparams.b);
	TBranch *branch_berr = tree1->Branch("berr", &lineparams.berr);
	TBranch *branch_chisq = tree1->Branch("chisq", &lineparams.chisq);

	// TODO: should we be taking arctan here
	float meanc1_angle =  atan(meanc1_a_diffs);
	float meanc3_angle = atan(meanc3_a_diffs);

	//meanc1_b_diffs = 0.5;//-0.504;
	//meanc3_b_diffs = 0.06;//1.01;
	//testing hard coded offsets -- didnt really help 

	vector<LineParts> fittedlines = fit_chamber(events, rfuncs, lineparams, branch_a, branch_aerr, branch_b, branch_berr, branch_chisq, meanc1_angle, meanc3_angle, meanc1_b_diffs, meanc3_b_diffs);

	// now calculate per tube efficiency

	vector<int> dist_lessthan_counter;
	vector<int> tube_hit_counter;

	double rotationmatc3[2][2] = {{cos(meanc3_angle), -sin(meanc3_angle)}, {sin(meanc3_angle), cos(meanc3_angle)}};
	double rotationmatc1[2][2] = {{cos(meanc1_angle), -sin(meanc1_angle)}, {sin(meanc1_angle), cos(meanc1_angle)}};

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
					// rotate and shift the xy
					if (chamb == 0)
					{
						xy.at(0) = xy.at(0) - meanc1_b_diffs;
						xy.at(0) = (rotationmatc1[0][0] * xy.at(0) + rotationmatc1[0][1] * xy.at(1)) ;
						xy.at(1) = rotationmatc1[1][0] * xy.at(0) + rotationmatc1[1][1] * xy.at(1);
					}
					
					if (chamb == 2)
					{
						xy.at(0) =  xy.at(0) - meanc3_b_diffs;
						xy.at(0) = (rotationmatc3[0][0] * xy.at(0) + rotationmatc3[0][1] * xy.at(1));
						xy.at(1) = rotationmatc3[1][0] * xy.at(0) + rotationmatc3[1][1] * xy.at(1);
					}

					double d_line_tube = abs(xy.at(0) + line.a * xy.at(1) + line.b) / sqrt(line.a * line.a + 1);
					if (d_line_tube < 1.5) // TODO: tube radius
					{
						dist_lessthan_counter.at(tubenum) += 1;

						if (std::find(hittubenums.begin(), hittubenums.end(), tubenum) != hittubenums.end())
							tube_hit_counter.at(tubenum) += 1;
					}
				}
	}

	double tubeeff;
	TBranch *tube_eff = tree1->Branch("tube_eff", &tubeeff);

	for (int i = 0; i < 144; i++)
	{
		tubeeff = (double)tube_hit_counter.at(i) / (double)dist_lessthan_counter.at(i);
		tube_eff->Fill();
	}

	tree1->SetEntries(events.size());

	// Write the tree to the file and close the file to save.
	tree1->Write();
	file->Close();

	std::cout << "done" << std::endl;

	// app.Run();

	return 0;
}
