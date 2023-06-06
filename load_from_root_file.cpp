#include "load_from_root_file.h"

// generate a vector of type "event", populated with relevant data
vector<NewEvent> loadTreeFromFile(string filename)
{
	UInt_t eventNumber = 0;
	UInt_t nhits = 0;
	vector<unsigned int> *bcid = 0;
	vector<unsigned int> *tdc = 0;
	vector<unsigned int> *channel = 0;
	vector<unsigned int> *coarse = 0;
	vector<unsigned int> *fine = 0;
	vector<unsigned int> *chamber = 0;
	vector<unsigned int> *layer = 0;
	vector<unsigned int> *tube = 0;
	vector<float> *times = 0;
	vector<float> *charge = 0;
	vector<bool> *leading = 0;

	TBranch *b_eventNumber;
	TBranch *b_nhits;
	TBranch *b_bcid;
	TBranch *b_tdc;
	TBranch *b_channel;
	TBranch *b_coarse;
	TBranch *b_fine;
	TBranch *b_chamber;
	TBranch *b_layer;
	TBranch *b_tube;
	TBranch *b_times;
	TBranch *b_charge;
	TBranch *b_leading;

	TFile *f = new TFile((filename.c_str()));
	TTree *tree = (TTree *)f->Get("mdt");

	tree->SetBranchAddress("chamber", &chamber, &b_chamber);
	tree->SetBranchAddress("layer", &layer, &b_layer);
	tree->SetBranchAddress("tube", &tube, &b_tube);
	tree->SetBranchAddress("time", &times, &b_times);
	tree->SetBranchAddress("charge", &charge, &b_charge);
	tree->SetBranchAddress("leading", &leading, &b_leading);

	long int entries = tree->GetEntriesFast();

	std::cout << "entries is " << entries << std::endl;

	vector<Event> events;

	vector<NewEvent> evs;

	for (int entry = 0; entry < entries; entry++)
	{

		tree->LoadTree(entry);

		tree->GetEntry(entry);

		int numhits = times->size();

		NewEvent newevent;

		for (int i = 0; i < numhits; i++)
		{
			if (leading->at(i) > 0 && charge->at(i) > 40 && charge->at(i) < 600) //&& numhits <= 9)
			{
				int chambnum = chamber->at(i) - 1;
				int layernum = layer->at(i) - 1;
				int tubenum = tube->at(i) - 1;

				newevent.eventnum = entry;
				newevent.t.push_back(times->at(i));
				newevent.charge.push_back(charge->at(i));
				newevent.chamber.push_back(chambnum);
				newevent.layer.push_back(layernum);
				newevent.tube.push_back(tubenum);
				newevent.is_inlier.push_back(1);
			}
		}

		evs.push_back(newevent);
	}

	delete tree;

	f->Close();

	delete f;

	return evs;
}

void addeventstotree(vector<NewEvent> events, NewEvent &event, TBranch *branch)
{
	for (int i = 0; i < events.size(); i++)
	{
		event = events.at(i);
		branch->Fill();
	}
}
