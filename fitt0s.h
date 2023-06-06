#pragma once
#include "TH1F.h"
#include "TF1.h"
#include <vector>
#include "line_fitting.h"
#include "radiustimefunction.h"

using std::vector;
using std::string;

typedef struct Fittedt0s
{
	double t0;
	double t0err;
	TF1* fitfn;
} Fittedt0s;

Fittedt0s fit_t0( vector<float> times );
vector<TF1> populate_t0s(vector<NewEvent> events, double &t0, double &t0err, TF1 &t0fit, TBranch *branch1, TBranch *branch2, TBranch *branch3);