#include "line_fitting.h"
#include <TFile.h>
#include <TTree.h>

vector<NewEvent> loadTreeFromFile(string filename);
void addeventstotree(vector<NewEvent> events, NewEvent &event, TBranch *branch);