#pragma once
#include "TMinuit.h"
#include "TF1.h"
#include <TTree.h>
#include <TClonesArray.h>
#include "TCanvas.h"
#include <vector>
#include "Event.h"
#include "TH1F.h"

using std::string;
using std::vector;

extern vector<float> xvs;
extern vector<float> yvs;
extern vector<float> rvs;
extern vector<float> extern_t0s;
extern vector<float> sigmas;
extern vector<float> as;
extern int what_entrynum;
extern vector<int> chlist;

ClassImp(NewEvent) // this is needed for root to be able to read the class from a tree. Generates a dictionary for the class.
	/**
	 * struct
	 *
	 */
	typedef struct Hit
{
	long int eventnum;
	float t;
	float charge;
	int chamber;
	int layer;
	int tube;
} Hit;

typedef struct Event
{
	vector<Hit> hits;
} Event;

// //this will become the actual events
// typedef struct NewEvent
// {
// 	int eventnum;
// 	vector<float> t;
// 	vector<float> charge;
// 	vector<int> chamber;
// 	vector<int> layer;
// 	vector<int> tube ;
// } NewEvent;

typedef struct LineParts
{
	double a = 0;
	double b = 0;
	double c = 0; 
	double ch1 = 0;
	double ch3 = 0;

	double aerr = 0;
	double berr = 0;
	double cerr = 0;
	double ch1err = 0;
	double ch3err = 0;

	double chisq = 0;

	int eventNum = 0;
} LineParts;

LineParts subtractlineparts(LineParts line1, LineParts line2);

vector<float> getTubeCoords(int chamber, int layer, int tube);

LineParts justfitlines(int setfn, int event, vector<float> gransac_lineparams);

// This is the function that fits the lines for each chamber. Ignores specified layer in layer_to_ignore variable, as per the professor's request. Layer is 0-8, 0 being the top layer.
vector<LineParts> fit_chamber(vector<NewEvent> events, vector<TF1> rfuncs, int layer_to_ignore);

vector<LineParts> fit_chamber(vector<NewEvent> events, vector<TF1> rfuncs, LineParts &lineparams, TBranch *branch_a, TBranch *branch_aerr,
							  TBranch *branch_b, TBranch *branch_berr,
							  TBranch *branch_chisq, float c1_angle, float c3_angle, float meanc1_b_diffs,
							  float meanc3_b_diffs, vector<float> gransac_lineparams, TBranch *branch_ch1,
							  TBranch *branch_ch1err, TBranch *branch_ch3, TBranch *ch3err);

vector<vector<LineParts>> fit_single_chambers(vector<NewEvent> events, vector<TF1> rfuncs, LineParts &lineparamsc1, LineParts &lineparamsc2, LineParts &lineparamsc3,
											  TBranch *branch_ac1, TBranch *branch_aerrc1, TBranch *branch_bc1, TBranch *branch_berrc1, TBranch *branch_chisqc1,
											  TBranch *branch_ac2, TBranch *branch_aerrc2, TBranch *branch_bc2, TBranch *branch_berrc2, TBranch *branch_chisqc2,
											  TBranch *branch_ac3, TBranch *branch_aerrc3, TBranch *branch_bc3, TBranch *branch_berrc3, TBranch *branch_chisqc3,
											  TBranch *branch_evnumc1, TBranch *branch_evnumc2, TBranch *branch_evnumc3);

vector<LineParts> fit_single_chamber(int chambernumber, int setfn, vector<NewEvent> events, vector<TF1> rfuncs, LineParts &lineparamsc1,
									 TBranch *branch_ac1, TBranch *branch_aerrc1, TBranch *branch_bc1, TBranch *branch_berrc1, TBranch *branch_chisqc1,
									 TBranch *branch_evnumc1, vector<float> gransac_lineparams);

vector<LineParts> fit_single_chamber(int chambernumber, int setfn, double rotationangle, vector<NewEvent> events, vector<TF1> rfuncs, LineParts &lineparamsc1, TBranch *branch_bc1, TBranch *branch_berrc1);

float fittwochambers(int c1, int c2, int c_fit, vector<NewEvent> events, vector<TF1> rfuncs);