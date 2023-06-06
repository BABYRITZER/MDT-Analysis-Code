#include <iostream>
#include <vector>
#include "line_fitting.h"
#include "gransac/GRANSAC.hpp"
#include "gransac/LineModel.hpp"
#include "gransac/AbstractModel.hpp"

void findInliers(vector<NewEvent> events, NewEvent &processed_event, TBranch *branch_eventnum, TBranch *branch_t, TBranch *branch_chg, TBranch *branch_chmb, TBranch *branch_layer, TBranch *branch_tube, TBranch *branch_is_inliner);