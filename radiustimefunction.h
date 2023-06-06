#pragma once
#include "TH1F.h"
#include "TF1.h"
#include <vector>

using std::vector;
using std::string;

TF1* radius_for_time(vector<float> times, float t0);