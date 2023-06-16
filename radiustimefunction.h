#pragma once
#include "TH1F.h"
#include "TF1.h"
#include <vector>
#include "TCanvas.h"
#include "TPad.h"
#include "TGraph.h"

using std::vector;
using std::string;

TF1* radius_for_time(vector<float> times, float t0);

//if u want to save the grpahs
TF1* radius_for_time(vector<float> times, float t0, int chambernum);
