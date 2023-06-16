#pragma once
#include "TMinuit.h"
#include "TF1.h"
#include <TTree.h>
#include <TClonesArray.h>
#include <vector>
#include "Event.h"

using std::string;
using std::vector;

extern vector<float> xvs;
extern vector<float> yvs;
extern vector<float> rvs;
extern vector<float> extern_t0s;
extern vector<float> sigmas;
extern vector<float> as;
extern int what_entrynum;

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

	double aerr = 0;
	double berr = 0;
	double cerr = 0;

	double chisq = 0;

	int eventNum = 0;
} LineParts;

vector<float> getTubeCoords(int chamber, int layer, int tube);

LineParts justfitlines(int setfn = 0);

vector<LineParts> fit_chamber(vector<NewEvent> events, vector<TF1> rfuncs, LineParts &lineparams, TBranch *branch_a, TBranch *branch_aerr, TBranch *branch_b, TBranch *branch_berr, TBranch *branch_chisq, float c1_angle, float c3_angle, float meanc1_b_diffs, float meanc3_b_diffs);

vector<vector<LineParts>> fit_single_chambers(vector<NewEvent> events, vector<TF1> rfuncs, LineParts &lineparamsc1, LineParts &lineparamsc2, LineParts &lineparamsc3,
											  TBranch *branch_ac1, TBranch *branch_aerrc1, TBranch *branch_bc1, TBranch *branch_berrc1, TBranch *branch_chisqc1,
											  TBranch *branch_ac2, TBranch *branch_aerrc2, TBranch *branch_bc2, TBranch *branch_berrc2, TBranch *branch_chisqc2,
											  TBranch *branch_ac3, TBranch *branch_aerrc3, TBranch *branch_bc3, TBranch *branch_berrc3, TBranch *branch_chisqc3,
											  TBranch *branch_evnumc1, TBranch *branch_evnumc2, TBranch *branch_evnumc3);

vector<LineParts> fit_single_chamber(int chambernumber, int setfn, vector<NewEvent> events, vector<TF1> rfuncs, LineParts &lineparamsc1,
									 TBranch *branch_ac1, TBranch *branch_aerrc1, TBranch *branch_bc1, TBranch *branch_berrc1, TBranch *branch_chisqc1,
									 TBranch *branch_evnumc1);

vector<LineParts> fit_single_chamber(int chambernumber, int setfn, double rotationangle, vector<NewEvent> events, vector<TF1> rfuncs, LineParts &lineparamsc1, TBranch *branch_bc1, TBranch *branch_berrc1);
