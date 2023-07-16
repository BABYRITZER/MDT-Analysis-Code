#pragma once
#include "Event.h"
#include <vector>
#include "line_fitting.h"
#include <utility>

// Calculates the efficiency of the tubes. Returns tuple of (efficiency, +error, -error) vectors. Also, fills branches with the raw numbers for efficiency.
std::tuple<vector<float>, vector<float>, vector<float>> eff_calc(vector<NewEvent> events, vector<LineParts> fittedlines, vector<TF1> rfuncs, TBranch *nhitsbranch, TBranch *distlessbranch, int &efftop, int &effbot);

// Returns the distance from a point to a line and the error in that distance (in that order)
std::pair<float, float> dist_calc_with_error(float x, float y, LineParts line);

std::tuple<vector<int>, vector<int>, vector<int>, vector<int>, vector<int>, vector<int>>
layer_effcalc(vector<NewEvent> events, vector<LineParts> fittedlines, vector<TF1> rfuncs);
