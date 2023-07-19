#include <vector>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TChain.h>
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

using std::string;
using std::to_string;
using std::vector;

typedef struct
{
    float eff;
    float etop;
    float ebot;
    float err;

} Efficiency;

vector<string> analysisfiles = {"reco_run988_analysis.root",
                                "reco_run525_analysis.root",
                                "reco_run526_analysis.root",
                                "reco_run527_analysis.root",
                                "reco_run528_analysis.root",
                                //"reco_run529_analysis.root",
                                "reco_run2850_analysis.root",
                                "reco_run2800_analysis.root",
                                "reco_run522_analysis.root"};
/*"reco_run522_analysis.root",
 "reco_run525_analysis.root",
 "reco_run526_analysis.root",
 "reco_run527_analysis.root",
 "reco_run528_analysis.root",
 "reco_run529_analysis.root",
 "reco_run2850_analysis.root"};*/

 //vector<string> analysisfiles = {"reco_run525_analysis.root"};

// Declaration of leaf types
Double_t t0s;
Double_t t0_errors;
TF1 *t0_fits;
vector<TF1 *> *r_functions;
Int_t event_num;

int efftopv;
int effbotv;
float layereff_top;
float layereff_bot;

// vector<float> *time;
vector<float> *charge;
vector<int> *chamber;
vector<int> *layer;
vector<int> *tube;
vector<int> *is_inliner;
Double_t a;
Double_t aerr;
Double_t b;
Double_t berr;
Double_t chisq;
Double_t a_sys_plus;
Double_t b_sys_plus;
Double_t a_sys_minus;
Double_t b_sys_minus;
Double_t eff;
Double_t eff_err_plus;
Double_t eff_err_minus;
Double_t eff_sys_plus;
Double_t eff_sys_minus;

// List of branches
TBranch *b_t0s;
TBranch *b_t0_errors;
TBranch *b_t0_fits;
TBranch *b_r_functions;
TBranch *b_event_num;
// TBranch *b_time;
TBranch *b_charge;
TBranch *b_chamber;
TBranch *b_layer;
TBranch *b_tube;
TBranch *b_is_inliner;
TBranch *b_a;
TBranch *b_aerr;
TBranch *b_b;
TBranch *b_berr;
TBranch *b_chisq;
TBranch *b_a_sys_plus;
TBranch *b_b_sys_plus;
TBranch *b_a_sys_minus;
TBranch *b_b_sys_minus;
TBranch *b_eff;
TBranch *b_eff_err_plus;
TBranch *b_eff_err_minus;
TBranch *b_eff_sys_plus;
TBranch *b_eff_sys_minus;

TBranch *b_efftop;
TBranch *b_effbot;

TBranch *b_layereff_top;
TBranch *b_layereff_bot;

int main()
{

    vector<int> voltages;

    int initialvoltage = 3100;

    for (int i = 0; i < analysisfiles.size(); i++)
        voltages.push_back(initialvoltage - i * 50);

    vector<vector<Efficiency>> efficiencies_all;

    for (string filename : analysisfiles)
    {

        TFile *f = new TFile((filename.c_str()));
        TTree *tree = (TTree *)f->Get("myTree");

        tree->SetBranchAddress("t0s", &t0s, &b_t0s);
        tree->SetBranchAddress("t0_errors", &t0_errors, &b_t0_errors);
        tree->SetBranchAddress("t0_fits", &t0_fits, &b_t0_fits);
        tree->SetBranchAddress("r_functions", &r_functions, &b_r_functions);
        tree->SetBranchAddress("event_num", &event_num, &b_event_num);
        //  tree->SetBranchAddress("time", &time, &b_time);
        tree->SetBranchAddress("charge", &charge, &b_charge);
        tree->SetBranchAddress("chamber", &chamber, &b_chamber);
        tree->SetBranchAddress("layer", &layer, &b_layer);
        tree->SetBranchAddress("tube", &tube, &b_tube);
        tree->SetBranchAddress("is_inliner", &is_inliner, &b_is_inliner);
        tree->SetBranchAddress("a", &a, &b_a);
        tree->SetBranchAddress("aerr", &aerr, &b_aerr);
        tree->SetBranchAddress("b", &b, &b_b);
        tree->SetBranchAddress("berr", &berr, &b_berr);
        tree->SetBranchAddress("chisq", &chisq, &b_chisq);
        tree->SetBranchAddress("a_sys_plus", &a_sys_plus, &b_a_sys_plus);
        tree->SetBranchAddress("b_sys_plus", &b_sys_plus, &b_b_sys_plus);
        tree->SetBranchAddress("a_sys_minus", &a_sys_minus, &b_a_sys_minus);
        tree->SetBranchAddress("b_sys_minus", &b_sys_minus, &b_b_sys_minus);
        tree->SetBranchAddress("eff", &eff, &b_eff);
        tree->SetBranchAddress("eff_err_plus", &eff_err_plus, &b_eff_err_plus);
        tree->SetBranchAddress("eff_err_minus", &eff_err_minus, &b_eff_err_minus);
        tree->SetBranchAddress("eff_sys_plus", &eff_sys_plus, &b_eff_sys_plus);
        tree->SetBranchAddress("eff_sys_minus", &eff_sys_minus, &b_eff_sys_minus);

        tree->SetBranchAddress("eff_top", &efftopv, &b_efftop);
        tree->SetBranchAddress("eff_bot", &effbotv, &b_effbot);

        tree->SetBranchAddress("layer_eff_top", &layereff_top, &b_layereff_top);
        tree->SetBranchAddress("layer_eff_bot", &layereff_bot, &b_layereff_bot);

        long int entries = 144;

        std::cout << "entries is " << entries << std::endl;

        vector<Efficiency> efficiencies;

        for (int i = 0; i < entries; i++)
        {
            tree->GetEntry(i);

            Efficiency efficiency;

            float eff_val = (float)layereff_top / (float)layereff_bot;

            efficiency.eff = isnan(eff_val) ? 0 : eff_val; // If the efficiency is nan, set it to 0

            efficiency.etop = isnan(layereff_top) ? 0 : layereff_top;
            efficiency.ebot = isnan(layereff_bot) ? 0 : layereff_bot;
            efficiency.err =  sqrt(efficiency.eff*(1-efficiency.eff)/efficiency.ebot);

            efficiencies.push_back(efficiency);


        }

        efficiencies_all.push_back(efficiencies);

        delete tree;

        f->Close();

        delete f;
    }

    // vector<vector<float>> tubeeff; // filenames x 9 size array
    // vector<vector<float>> tubeerrs;
    // for (int i = 0; i < analysisfiles.size(); ++i)
    // {
    //     vector<float> wow;
    //     tubeeff.push_back(wow);
    //     tubeerrs.push_back(wow);
    // }

    // for (int i = 0; i < analysisfiles.size(); ++i)
    // {
    //     vector<float> tubeeffs;
    //     vector<float> errs;

    //     for (int j = 0; j < efficiencies_all.at(i).size(); ++j)
    //     {
    //         float tubeeffn = (float)efficiencies_all.at(i).at(j).etop / (float)efficiencies_all.at(i).at(j).ebot;

    //         tubeeffs.push_back(tubeeffn);

    //         float asdf = ((tubeeffn) * (1. - tubeeffn)) / (float)efficiencies_all.at(i).at(j).ebot;

    //         errs.push_back(asdf);
    //     }

    //     tubeeff.at(i) = tubeeffs;
    //     tubeerrs.at(i) = errs;
    // }

    float layereffs[analysisfiles.size()][9];
    float layererrs[analysisfiles.size()][9];

    for (int i = 0; i < analysisfiles.size(); ++i)
        for (int j = 0; j < 9; ++j)
        {
            layereffs[i][j] = efficiencies_all.at(i).at(j).eff;
            layererrs[i][j] = efficiencies_all.at(i).at(j).err;
        }

    // for (int i = 0; i < analysisfiles.size(); ++i)
    //     for (int j = 0; j < 9; ++j)
    //     {
    //         layereffs[i][j] = 0;
    //         layererrs[i][j] = 0;
    //     }

    // int counter[analysisfiles.size()][9];
    // for (int i = 0; i < analysisfiles.size(); ++i)
    //     for (int j = 0; j < 9; ++j)
    //         counter[i][j] = 0;

    // for (int i = 0; i < analysisfiles.size(); ++i)
    // {
    //     for (int j = 0; j < 144; ++j)
    //     {
    //         if (tubeeff.at(i).at(j) > 0. && tubeeff.at(i).at(j) < 1.1)
    //         {
    //             counter[i][j % 9] += 1;
    //             layereffs[i][j % 9] += tubeeff.at(i).at(j);
    //             layererrs[i][j % 9] += tubeerrs.at(i).at(j) * tubeerrs.at(i).at(j);
    //         }
    //     }
    // }

    // for (int i = 0; i < analysisfiles.size(); ++i)
    //     for (int j = 0; j < 9; ++j)
    //     {
    //         layereffs[i][j] = layereffs[i][j] / (float)counter[i][j];
    //         layererrs[i][j] = layereffs[i][j] / (float)counter[i][j];
    //     }

    for (int i = 0; i < analysisfiles.size(); ++i)
    {
        std::cout << "File " << analysisfiles.at(i) << " (Voltage: " << voltages.at(i) << ")" << std::endl;
        for (int j = 0; j < 9; ++j)
        {
            std::cout << "\t" << layereffs[i][j] << " +- " << layererrs[i][j] << std::endl;
        }
    }

    TCanvas *c2 = new TCanvas("c2", "c2", 1000, 1000);
    c2->Divide(3, 3);

    for (int i = 0; i < 9; ++i)
    {
        c2->cd(i + 1);
        vector<double> x;
        vector<double> y;
        vector<double> x_err;
        vector<double> y_err_positive;
        vector<double> y_err_negative;

        for (int j = 0; j < analysisfiles.size(); ++j)
        {
            x.push_back(3100 - 50 * j);
            y.push_back(layereffs[j][i]);
            x_err.push_back(0);
            y_err_positive.push_back(layererrs[j][i]);
            y_err_negative.push_back(layererrs[j][i]);
        }

        TGraphAsymmErrors *gr = new TGraphAsymmErrors(x.size(), &x[0], &y[0], &x_err[0], &x_err[0], &y_err_negative[0], &y_err_positive[0]);
        string title = "Layer " + to_string(i + 1) + ";Chamber 2 Voltage (V);Efficiency";
        gr->SetTitle(title.c_str());
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(1);
        gr->SetMarkerColor(kBlack);
        gr->Draw("AP");
    }

    c2->SaveAs("efficiencies_layers.pdf");

    return 0;
}