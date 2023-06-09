#include <TROOT.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TPaveLabel.h>
#include <algorithm>
#include "TMinuit.h"
#include "TBranch.h"
#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include <vector>
#include "Event.h"
#include <TF1.h>
#include <TChain.h>
#include <TLeaf.h>
#include <TApplication.h>
#include <TGraph.h>

#include "fitt0s.h"
#include "radiustimefunction.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;


int main(int argc, char **argv)
{
    TApplication app("Root app", &argc, argv);

    std::cout << " guh " << std::endl;

    TFile *f = new TFile(("reco_run988_analysis.root"));
    TTree *tree = (TTree *)f->Get("myTree");

    Double_t t0s;
    Double_t t0_errors;
    TF1 *t0_fits;
    vector<TF1 *> *r_functions;
    Int_t event_num;
    vector<float> *time;
    vector<float> *charge;
    vector<int> *chamber;
    vector<int> *layer;
    vector<int> *tube;
    vector<int> *is_inliner;
    Double_t a_c1;
    Double_t aerr_c1;
    Double_t b_c1;
    Double_t berr_c1;
    Double_t c1_chisq;
    Double_t a_c2;
    Double_t aerr_c2;
    Double_t b_c2;
    Double_t berr_c2;
    Double_t c2_chisq;
    Double_t a_c3;
    Double_t aerr_c3;
    Double_t b_c3;
    Double_t berr_c3;
    Double_t c3_chisq;
    Double_t a;
    Double_t aerr;
    Double_t b;
    Double_t berr;
    Double_t chisq;

    TBranch *b_t0s;         //!
    TBranch *b_t0_errors;   //!
    TBranch *b_t0_fits;     //!
    TBranch *b_r_functions; //!
    TBranch *b_event_num;   //!
    TBranch *b_time;        //!
    TBranch *b_charge;      //!
    TBranch *b_chamber;     //!
    TBranch *b_layer;       //!
    TBranch *b_tube;        //!
    TBranch *b_is_inliner;  //!
    TBranch *b_a_c1;        //!
    TBranch *b_aerr_c1;     //!
    TBranch *b_b_c1;        //!
    TBranch *b_berr_c1;     //!
    TBranch *b_c1_chisq;    //!
    TBranch *b_a_c2;        //!
    TBranch *b_aerr_c2;     //!
    TBranch *b_b_c2;        //!
    TBranch *b_berr_c2;     //!
    TBranch *b_c2_chisq;    //!
    TBranch *b_a_c3;        //!
    TBranch *b_aerr_c3;     //!
    TBranch *b_b_c3;        //!
    TBranch *b_berr_c3;     //!
    TBranch *b_c3_chisq;    //!
    TBranch *b_a;           //!
    TBranch *b_aerr;        //!
    TBranch *b_b;           //!
    TBranch *b_berr;        //!
    TBranch *b_chisq;       //!

    // Set branch addresses and branch pointers
    t0_fits = 0;
    r_functions = 0;
    time = 0;
    charge = 0;
    chamber = 0;
    layer = 0;
    tube = 0;
    is_inliner = 0;

    tree->SetBranchAddress("t0s", &t0s, &b_t0s);
    tree->SetBranchAddress("t0_errors", &t0_errors, &b_t0_errors);
    tree->SetBranchAddress("t0_fits", &t0_fits, &b_t0_fits);
    tree->SetBranchAddress("r_functions", &r_functions, &b_r_functions);
    tree->SetBranchAddress("event_num", &event_num, &b_event_num);
    tree->SetBranchAddress("time", &time, &b_time);
    tree->SetBranchAddress("charge", &charge, &b_charge);
    tree->SetBranchAddress("chamber", &chamber, &b_chamber);
    tree->SetBranchAddress("layer", &layer, &b_layer);
    tree->SetBranchAddress("tube", &tube, &b_tube);
    tree->SetBranchAddress("is_inliner", &is_inliner, &b_is_inliner);
    tree->SetBranchAddress("a_c1", &a_c1, &b_a_c1);
    tree->SetBranchAddress("aerr_c1", &aerr_c1, &b_aerr_c1);
    tree->SetBranchAddress("b_c1", &b_c1, &b_b_c1);
    tree->SetBranchAddress("berr_c1", &berr_c1, &b_berr_c1);
    tree->SetBranchAddress("c1_chisq", &c1_chisq, &b_c1_chisq);
    tree->SetBranchAddress("a_c2", &a_c2, &b_a_c2);
    tree->SetBranchAddress("aerr_c2", &aerr_c2, &b_aerr_c2);
    tree->SetBranchAddress("b_c2", &b_c2, &b_b_c2);
    tree->SetBranchAddress("berr_c2", &berr_c2, &b_berr_c2);
    tree->SetBranchAddress("c2_chisq", &c2_chisq, &b_c2_chisq);
    tree->SetBranchAddress("a_c3", &a_c3, &b_a_c3);
    tree->SetBranchAddress("aerr_c3", &aerr_c3, &b_aerr_c3);
    tree->SetBranchAddress("b_c3", &b_c3, &b_b_c3);
    tree->SetBranchAddress("berr_c3", &berr_c3, &b_berr_c3);
    tree->SetBranchAddress("c3_chisq", &c3_chisq, &b_c3_chisq);
    tree->SetBranchAddress("a", &a, &b_a);
    tree->SetBranchAddress("aerr", &aerr, &b_aerr);
    tree->SetBranchAddress("b", &b, &b_b);
    tree->SetBranchAddress("berr", &berr, &b_berr);
    tree->SetBranchAddress("chisq", &chisq, &b_chisq);

    std::cout << " guh " << std::endl;

    long int entries = 100000;

    std::cout << entries << std::endl;

    // get rt fns
    vector<vector<float>> times;
    vector<float> itm;
    for (int i = 0; i < 144; i++)
        times.push_back(itm);

    vector<float> t0_vals;
    for (int i = 0; i < entries; i++)
    {
        tree->GetEntry(i);

        for (int j = 0; j < time->size(); j++)
        {
            int tubenum = chamber->at(j) * 48 + layer->at(j) * 16 + tube->at(j);
            times.at(tubenum).push_back(time->at(j));
            t0_vals.push_back(t0s);
        }
    }

    vector<TF1 *> rfns;
    for (int i = 0; i < 144; i++)
        rfns.push_back(radius_for_time(times.at(i), t0_vals.at(i)));

    // start
    TCanvas *c1 = new TCanvas("c1", "c1", 495, 696);
    c1->Range(0, 0, 49.5, 69.6);

    for (int entry = 0; entry < entries; entry++)
    {
        tree->LoadTree(entry);

        tree->GetEntry(entry);

        std::cout << "this is event " << entry << std::endl;

        vector<int> hittubenums;
        vector<float> hittubetimes;

        for (int i = 0; i < time->size(); i++)
        {
            int chambnum = chamber->at(i) - 1;
            int layernum = layer->at(i) - 1;
            int tubenum = tube->at(i) - 1;
            // 48*(chamber -1)+16*layer-1+tube-1==0 && charge > 40 && leading > 1
            // int finallayer = 3 * chambnum + layernum;
            int realtubenum = 48 * (chambnum) + 16 * (layernum) + tubenum;

            hittubenums.push_back(realtubenum);
            //TODO
            float toime = time->at(i) - 740.;
            hittubetimes.push_back(toime);
            // hitsreconst.at(jentry)->Fill( tubenum, finallayer );
        }

        int tubeNum;
        int layerNum = 0;
        int chambNum;
        float x;
        float y;

        for (chambNum = 0; chambNum < 3; chambNum++)
        {
            for (layerNum = 0; layerNum < 3; layerNum++)
            {
                for (tubeNum = 0; tubeNum < 16; tubeNum++)
                {
                    x = layerNum == 1 ? 1.5 + tubeNum * 3 : 3 + tubeNum * 3;
                    y = 1.5 + 2.6 * layerNum + (22.5 + 8.2) * chambNum;

                    TEllipse *el = new TEllipse(x, 69.6 - y, 1.5, 1.5);

                    int tubenum = 48 * (chambNum - 1) + 16 * (layerNum - 1) + (tubeNum - 1);
                    int tubenum2 = 48 * chambNum + 16 * layerNum + tubeNum;

                    if (std::find(hittubenums.begin(), hittubenums.end(), tubenum) != hittubenums.end())
                    {
                        int tubeloc = std::find(hittubenums.begin(), hittubenums.end(), tubenum) - hittubenums.begin();
                        float tubetime = hittubetimes.at(tubeloc);

                        float radius = (rfns.at(tubeloc)->Eval(tubetime));

                        TEllipse *innercircle = new TEllipse(x, 69.6 - y, radius, radius);

                        el->SetFillStyle(1001);
                        if (is_inliner->at(tubeloc) == 1)
                            el->SetFillColor(6);
                        else
                            el->SetFillColor(2);

                        el->SetLineColor(6);
                        el->Draw();

                        innercircle->SetFillStyle(1001);
                        innercircle->SetFillColor(4);
                        innercircle->SetLineColor(4);
                        innercircle->Draw();
                    }

                    else
                    {
                        el->SetFillStyle(1001);
                        el->SetLineColor(4);
                        el->Draw();
                    }
                }
            }
        }

        int bins = 30;

        float yvalues[bins];
        float xvaluesc1[bins];
        float xvaluesc2[bins];
        float xvaluesc3[bins];
        float xvaluesc[bins];

        for (int i = 0; i < bins; i++)
        {
            y = ((float)i * 69.6) / ((float)bins);
            yvalues[i] = y;

            xvaluesc1[i] = (-1. * (a_c1 * y + b_c1));
            xvaluesc2[i] = (-1. * (a_c2 * y + b_c2));
            xvaluesc3[i] = (-1. * (a_c3 * y + b_c3));
            xvaluesc[i] = (-1. * (a * y + b));
        }


        std::cout << "c1 fit in green, c2 fit in red, c3 fit in blue and fit over all chambers in cyan." << std::endl;

        std::cout << "chisquare for c1 fit is " << c1_chisq << std::endl;
        std::cout << "chisquare for c2 fit is " << c2_chisq << std::endl;
        std::cout << "chisquare for c3 fit is " << c3_chisq << std::endl;
        std::cout << "chisquare for all chambers fit is " << chisq << std::endl;

        TGraph *guh0 = new TGraph(bins, xvaluesc1, yvalues);
        guh0->SetLineColor(kGreen);
        guh0->Draw("same");

        TGraph *guh1 = new TGraph(bins, xvaluesc2, yvalues);
        guh1->SetLineColor(kRed);
        guh1->Draw("same");

        TGraph *guh2 = new TGraph(bins, xvaluesc3, yvalues);
        guh2->SetLineColor(kBlue);
        guh2->Draw("same");

        TGraph *guh3 = new TGraph(bins, xvaluesc, yvalues);
        guh3->SetLineColor(kCyan);
        guh3->Draw("same");

        c1->Modified();
        c1->Modified();
        c1->Update();

        std::cin.get();
        c1->GetListOfPrimitives()->Remove(guh0);
        c1->GetListOfPrimitives()->Remove(guh1);
        c1->GetListOfPrimitives()->Remove(guh2);
        c1->GetListOfPrimitives()->Remove(guh3);
    }
    app.Run();

    return 0;
}
