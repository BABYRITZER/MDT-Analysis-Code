#pragma once
#include <vector>
#include <TObject.h>

using std::vector;
/**
 * Event class. Inherits TObject so that it is recognizable in the TTree. Contains the event number, the time, charge, chamber, layer and tube of each hit.
 *  Member information:
    *      - Int_t eventnum: event number
 *       - vector<float> t: time of each hit
 *       - vector<float> charge: charge of each hit
 *       - vector<int> chamber: chamber of each hit
 *       - vector<int> layer: layer of each hit
 *      - vector<int> tube: tube of each hit
 */
class NewEvent : public TObject
{
public:
    Int_t eventnum;
    std::vector<Float_t> t;
    std::vector<Float_t> charge;
    std::vector<Int_t> chamber;
    std::vector<Int_t> layer;
    std::vector<Int_t> tube;
    std::vector<Int_t> is_inlier;

    NewEvent()
    {
        eventnum = 0;
        t.clear();
        charge.clear();
        chamber.clear();
        layer.clear();
        tube.clear();
        is_inlier.clear();
    };

    NewEvent(Int_t eventnum, std::vector<Float_t> t, std::vector<Float_t> charge, std::vector<Int_t> chamber, std::vector<Int_t> layer, std::vector<Int_t> tube, std::vector<Int_t> is_inlier)
    {
        this->eventnum = eventnum;
        this->t = t;
        this->charge = charge;
        this->chamber = chamber;
        this->layer = layer;
        this->tube = tube;
        this->is_inlier = is_inlier;
    };

    ClassDef(NewEvent, 1)
};