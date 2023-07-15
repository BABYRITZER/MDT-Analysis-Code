#include "efficiency_calc.h"
#include <cmath>
#include <utility>
#include "line_fitting.h"

using std::make_pair;
using std::pair;
using std::tuple;

pair<float, float> dist_calc_with_error(float x, float y, LineParts line)
{
    float dist = abs(line.a * y + line.b - x) / sqrt(line.a * line.a + 1);

    float dx = 0.;
    float dy = 0.;
    float da = line.aerr;
    float db = line.berr;

    double error_term1 = line.a * y + line.b - x;
    double error_term2 = sqrt(pow(line.a, 2) + 1);

    double err = abs(error_term1) * sqrt(
                                        -pow(1 / (error_term1 * error_term2), 2) * pow(dx, 2) +
                                        pow(1 / (error_term1 * error_term2), 2) * pow(db, 2) +
                                        pow(line.a, 2) * pow(1 / (error_term1 * error_term2), 2) * pow(dy, 2) +
                                        pow(pow(error_term2, 2) * y / (error_term1 * error_term2) - line.a, 2) * pow(da, 2) / pow(error_term2, 6));

    return make_pair(dist, err);
}

std::tuple<vector<float>, vector<float>, vector<float>> eff_calc(vector<NewEvent> events, vector<LineParts> fittedlines, vector<TF1> rfuncs,
                                                                 TBranch *topbranch, TBranch *botbranch, int &efftop, int &effbot)
{

    vector<int> dist_lessthan_counter;
    vector<int> tube_hit_counter;

    vector<int> dist_lessthan_counter_pluserr;
    vector<int> tube_hit_counter_pluserr;

    vector<int> dist_lessthan_counter_minuserr;
    vector<int> tube_hit_counter_minuserr;

    // we didnt calculate these so we just set them to 0
    double meanc3_angle = 0.;
    double meanc1_angle = 0.;

    double meanc1_b_diffs = 0.;
    double meanc3_b_diffs = 0.;

    double rotationmatc3[2][2] = {{cos(meanc3_angle), -sin(meanc3_angle)}, {sin(meanc3_angle), cos(meanc3_angle)}};
    double rotationmatc1[2][2] = {{cos(meanc1_angle), -sin(meanc1_angle)}, {sin(meanc1_angle), cos(meanc1_angle)}};

    for (int i = 0; i < 144; i++)
    {
        dist_lessthan_counter.push_back(0);
        tube_hit_counter.push_back(0);
        dist_lessthan_counter_pluserr.push_back(0);
        dist_lessthan_counter_minuserr.push_back(0);
        tube_hit_counter_pluserr.push_back(0);
        tube_hit_counter_minuserr.push_back(0);
    }

    for (int i = 0; i < events.size(); i++)
    {

        // FOR EFFICIENCY CUTTING CHISQUARE AT 50 -- using 0.1 resolution
        if (fittedlines.at(i).chisq > 50)
            continue;

        // TODO:
        if (events.at(i).t.size() < 6)
            continue;

        vector<int> hittubenums;

        for (int j = 0; j < events.at(i).t.size(); j++)
        {
            if (events.at(i).is_inlier.at(j) == 1)
            {
                int chambnum = events.at(i).chamber.at(j);
                int layernum = events.at(i).layer.at(j);
                int tubenum = events.at(i).tube.at(j);

                int realtubenum = 48 * (chambnum) + 16 * (layernum) + tubenum;

                hittubenums.push_back(realtubenum);
            }
        }

        for (int chamb = 0; chamb < 3; chamb++)
            for (int layer = 0; layer < 3; layer++)
            {
                int layer_to_ignore = layer;

                vector<NewEvent> ev;
                ev.push_back(events.at(i));

                vector<LineParts> wowowow = fit_chamber(ev, rfuncs, layer_to_ignore); //.at(0);\

                LineParts line = wowowow.at(0);

                ev.clear();

                for (int tube = 0; tube < 16; tube++)
                {
                    int tubenum = chamb * 48 + layer * 16 + tube;
                    vector<float> xy = getTubeCoords(chamb, layer, tube);
                    // rotate and shift the xy
                    if (chamb == 0)
                    {
                        xy.at(0) = xy.at(0) - meanc1_b_diffs;
                        xy.at(0) = (rotationmatc1[0][0] * xy.at(0) + rotationmatc1[0][1] * xy.at(1));
                        xy.at(1) = rotationmatc1[1][0] * xy.at(0) + rotationmatc1[1][1] * xy.at(1);
                    }

                    if (chamb == 2)
                    {
                        xy.at(0) = xy.at(0) - meanc3_b_diffs;
                        xy.at(0) = (rotationmatc3[0][0] * xy.at(0) + rotationmatc3[0][1] * xy.at(1));
                        xy.at(1) = rotationmatc3[1][0] * xy.at(0) + rotationmatc3[1][1] * xy.at(1);
                    }

                    // check if the line passes inside the tube at all
                    // double d_line_tube = abs(xy.at(0) + line.a * xy.at(1) - line.b) / sqrt(line.a * line.a + 1);

                    pair<float, float> dist_err = dist_calc_with_error(xy.at(0), xy.at(1), line);

                    float d_line_tube = dist_err.first;
                    float d_line_tube_err = dist_err.second;

                    if (d_line_tube < 1.46)
                    {
                        dist_lessthan_counter.at(tubenum) += 1;

                        if (std::find(hittubenums.begin(), hittubenums.end(), tubenum) != hittubenums.end())
                            tube_hit_counter.at(tubenum) += 1;
                    }

                    // Error stuff
                    if (abs(d_line_tube + d_line_tube_err) < 1.46)
                    {
                        dist_lessthan_counter_pluserr.at(tubenum) += 1;

                        if (std::find(hittubenums.begin(), hittubenums.end(), tubenum) != hittubenums.end())
                            tube_hit_counter_pluserr.at(tubenum) += 1;
                    }

                    if (abs(d_line_tube - d_line_tube_err) < 1.46)
                    {
                        dist_lessthan_counter_minuserr.at(tubenum) += 1;

                        if (std::find(hittubenums.begin(), hittubenums.end(), tubenum) != hittubenums.end())
                            tube_hit_counter_minuserr.at(tubenum) += 1;
                    }
                }
            }
    }

    vector<float> tube_eff;
    vector<float> tube_eff_pluserr;
    vector<float> tube_eff_minuserr;

    for (int i = 0; i < 144; i++)
        // what proportion are hits divided by what the line actually gives
        tube_eff.push_back((float)tube_hit_counter.at(i) / (float)dist_lessthan_counter.at(i));

    for (int i = 0; i < 144; i++)
    {
        efftop = tube_hit_counter_pluserr.at(i);
        effbot = dist_lessthan_counter_pluserr.at(i);

        topbranch->Fill();
        botbranch->Fill();

        tube_eff_pluserr.push_back((float)tube_hit_counter_pluserr.at(i) / (float)dist_lessthan_counter_pluserr.at(i));
        tube_eff_minuserr.push_back((float)tube_hit_counter_minuserr.at(i) / (float)dist_lessthan_counter_minuserr.at(i));
    }

    return std::make_tuple(tube_eff, tube_eff_pluserr, tube_eff_minuserr);
}

/*vector<float> layer_effcalc(vector<NewEvent> events, vector<LineParts> fittedlines, vector<TF1> rfuncs)
{

    int bottomnumber = 0;
    int topnumber;

    for (int i = 0; i < events.size(); i++)
    {

        // FOR EFFICIENCY CUTTING CHISQUARE AT 50 -- using 0.1 resolution
        if (fittedlines.at(i).chisq > 50)
            continue;

        // here since we need at least 5 hits in the line, we require now that we need 6 -- so that the line still fits with 5 points
        if (events.at(i).t.size() < 6)
            continue;

        vector<int> hittubenums;

        for (int j = 0; j < events.at(i).t.size(); j++)
        {
            if (events.at(i).is_inlier.at(j) == 1)
            {
                int chambnum = events.at(i).chamber.at(j);
                int layernum = events.at(i).layer.at(j);
                int tubenum = events.at(i).tube.at(j);

                int realtubenum = 48 * (chambnum) + 16 * (layernum) + tubenum;

                hittubenums.push_back(realtubenum);
            }
        }

        for (int chamb = 0; chamb < 3; chamb++)
            for (int layer = 0; layer < 3; layer++)
            {
                int layer_to_ignore = layer;

                vector<NewEvent> ev;
                ev.push_back(events.at(i));
                vector<LineParts> wowowow = fit_chamber(ev, rfuncs, layer_to_ignore);
                LineParts line = wowowow.at(0);
                ev.clear();

                bottomnumber += 1;

                int hitscount = 0;

                for (int q = 0; q < hittubenums.size(); ++q)

                    if (hitscount > 0)
                    {
                        topnumber += 1;
                    }
            }
    }
}*/
