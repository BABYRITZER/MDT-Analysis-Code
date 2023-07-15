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
    float eff_err_plus;
    float eff_err_minus;
    float eff_sys_plus;
    float eff_sys_minus;

    float etop;
    float ebot;

} Efficiency;

vector<string> analysisfiles ={"reco_run3100_analysis.root",
                                "reco_run525_analysis.root",
                                "reco_run526_analysis.root"};//,
                                //"reco_run527_analysis.root",
                                //"reco_run528_analysis.root",
                                //"reco_run529_analysis.root",
                               // "reco_run2800_analysis.root"};


/*{"reco_run522_analysis.root",
                                "reco_run525_analysis.root",
                                "reco_run526_analysis.root",
                                "reco_run527_analysis.root",
                                "reco_run528_analysis.root",
                                "reco_run529_analysis.root"};*/

// Declaration of leaf types
Double_t t0s;
Double_t t0_errors;
TF1 *t0_fits;
vector<TF1 *> *r_functions;
Int_t event_num;

int efftopv;
int effbotv;

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


        long int entries = 144;

        std::cout << "entries is " << entries << std::endl;

        vector<Efficiency> efficiencies;

        for (int i = 0; i < entries; i++)
        {
            tree->GetEntry(i);

            Efficiency efficiency;
            
            efficiency.eff = isnan(eff) ? 0 : eff; // If the efficiency is nan, set it to 0
            efficiency.eff_err_plus = isnan(eff_err_plus) ? 0 : eff_err_plus;
            efficiency.eff_err_minus = isnan(eff_err_minus) ? 0 : eff_err_minus;
            efficiency.eff_sys_plus = isnan(eff_sys_plus) ? 0 : eff_sys_plus;
            efficiency.eff_sys_minus = isnan(eff_sys_minus) ? 0 : eff_sys_minus;

            efficiency.etop = isnan(efftopv) ? 0 : efftopv;
            efficiency.etop = isnan(effbotv) ? 0 : effbotv;

            efficiencies.push_back(efficiency);
        }

        efficiencies_all.push_back(efficiencies);

        delete tree;

        f->Close();

        delete f;
    }

    // Average the efficiencies for each chamber
    vector<vector<Efficiency>> average_efficiencies_all_chambers; // Filesize x 3 (chambers)

    for (int i = 0; i < efficiencies_all.size(); ++i)
    {
        Efficiency average_efficiency_chamb1{0, 0, 0, 0, 0};
        Efficiency average_efficiency_chamb2{0, 0, 0, 0, 0};
        Efficiency average_efficiency_chamb3{0, 0, 0, 0, 0};

        std::cout << "File: " << analysisfiles.at(i) << std::endl;
        for (int j = 0; j < 144; ++j)
        {
            
            if (j < 48)
            {
                average_efficiency_chamb1.eff += efficiencies_all[i][j].eff;
                average_efficiency_chamb1.eff_err_plus += sqrt(efficiencies_all[i][j].eff * (1.-efficiencies_all[i][j].eff) / efficiencies_all[i][j].ebot);  ///efficiencies_all[i][j].eff_err_plus * efficiencies_all[i][j].eff_err_plus;
                average_efficiency_chamb1.eff_err_minus += sqrt(efficiencies_all[i][j].eff * (1.-efficiencies_all[i][j].eff) / efficiencies_all[i][j].ebot); //efficiencies_all[i][j].eff_err_minus * efficiencies_all[i][j].eff_err_minus;
                average_efficiency_chamb1.eff_sys_plus += efficiencies_all[i][j].eff_sys_plus * efficiencies_all[i][j].eff_sys_plus;
                average_efficiency_chamb1.eff_sys_minus += efficiencies_all[i][j].eff_sys_minus * efficiencies_all[i][j].eff_sys_minus;
            }
            else if (j < 96)
            {
                average_efficiency_chamb2.eff += efficiencies_all[i][j].eff;
                average_efficiency_chamb2.eff_err_plus += (efficiencies_all[i][j].eff * (1.-efficiencies_all[i][j].eff) )/ efficiencies_all[i][j].ebot;//efficiencies_all[i][j].eff_err_plus * efficiencies_all[i][j].eff_err_plus;
                average_efficiency_chamb2.eff_err_minus += (efficiencies_all[i][j].eff * (1.-efficiencies_all[i][j].eff) )/ efficiencies_all[i][j].ebot;//efficiencies_all[i][j].eff_err_minus * efficiencies_all[i][j].eff_err_minus;
                average_efficiency_chamb2.eff_sys_plus += efficiencies_all[i][j].eff_sys_plus * efficiencies_all[i][j].eff_sys_plus;
                average_efficiency_chamb2.eff_sys_minus += efficiencies_all[i][j].eff_sys_minus * efficiencies_all[i][j].eff_sys_minus;
            }
            else
            {
                average_efficiency_chamb3.eff += efficiencies_all[i][j].eff;
                average_efficiency_chamb3.eff_err_plus += (efficiencies_all[i][j].eff * (1.-efficiencies_all[i][j].eff) )/ efficiencies_all[i][j].ebot;//efficiencies_all[i][j].eff_err_plus * efficiencies_all[i][j].eff_err_plus;
                average_efficiency_chamb3.eff_err_minus += (efficiencies_all[i][j].eff * (1.-efficiencies_all[i][j].eff) )/ efficiencies_all[i][j].ebot;//efficiencies_all[i][j].eff_err_minus * efficiencies_all[i][j].eff_err_minus;
                average_efficiency_chamb3.eff_sys_plus += efficiencies_all[i][j].eff_sys_plus * efficiencies_all[i][j].eff_sys_plus;
                average_efficiency_chamb3.eff_sys_minus += efficiencies_all[i][j].eff_sys_minus * efficiencies_all[i][j].eff_sys_minus;
            }
        }

        average_efficiency_chamb1.eff /= 48;
        std::cout << "\tAverage Chamber 1 eff is " << average_efficiency_chamb1.eff << std::endl;
        average_efficiency_chamb1.eff_err_plus = sqrt(average_efficiency_chamb1.eff_err_plus / 48);
        average_efficiency_chamb1.eff_err_minus = sqrt(average_efficiency_chamb1.eff_err_minus / 48);
        average_efficiency_chamb1.eff_sys_plus = sqrt(average_efficiency_chamb1.eff_sys_plus / 48);
        average_efficiency_chamb1.eff_sys_minus = sqrt(average_efficiency_chamb1.eff_sys_minus / 48);

        average_efficiency_chamb2.eff /= 48;
        std::cout << "\tAverage Chamber 2 eff is " << average_efficiency_chamb2.eff << std::endl;
        average_efficiency_chamb2.eff_err_plus = sqrt(average_efficiency_chamb2.eff_err_plus / 48);
        average_efficiency_chamb2.eff_err_minus = sqrt(average_efficiency_chamb2.eff_err_minus / 48);
        average_efficiency_chamb2.eff_sys_plus = sqrt(average_efficiency_chamb2.eff_sys_plus / 48);
        average_efficiency_chamb2.eff_sys_minus = sqrt(average_efficiency_chamb2.eff_sys_minus / 48);

        average_efficiency_chamb3.eff /= 48;
        std::cout << "\tAverage Chamber 3 eff is " << average_efficiency_chamb3.eff << std::endl;
        average_efficiency_chamb3.eff_err_plus = sqrt(average_efficiency_chamb3.eff_err_plus / 48);
        average_efficiency_chamb3.eff_err_minus = sqrt(average_efficiency_chamb3.eff_err_minus / 48);
        average_efficiency_chamb3.eff_sys_plus = sqrt(average_efficiency_chamb3.eff_sys_plus / 48);
        average_efficiency_chamb3.eff_sys_minus = sqrt(average_efficiency_chamb3.eff_sys_minus / 48);

        vector<Efficiency> average_efficiencies;
        average_efficiencies.push_back(average_efficiency_chamb1);
        average_efficiencies.push_back(average_efficiency_chamb2);
        average_efficiencies.push_back(average_efficiency_chamb3);

        average_efficiencies_all_chambers.push_back(average_efficiencies);
    }

    std::cout << "Average efficiencies for all chambers: " << std::endl;
    for (int k = 0; k < average_efficiencies_all_chambers.size(); ++k)
    {
        std::cout << "\tRun " << k << ":" << std::endl;
        for (int l = 0; l < average_efficiencies_all_chambers[k].size(); ++l)
        {
            std::cout << "\t\tChamber " << l << ":" << std::endl;
            std::cout << "\t\t\tEfficiency: " << average_efficiencies_all_chambers[k][l].eff << std::endl;
            std::cout << "\t\t\tEfficiency error plus: " << average_efficiencies_all_chambers[k][l].eff_err_plus << std::endl;
            std::cout << "\t\t\tEfficiency error minus: " << average_efficiencies_all_chambers[k][l].eff_err_minus << std::endl;
            std::cout << "\t\t\tEfficiency systematic error plus: " << average_efficiencies_all_chambers[k][l].eff_sys_plus << std::endl;
            std::cout << "\t\t\tEfficiency systematic error minus: " << average_efficiencies_all_chambers[k][l].eff_sys_minus << std::endl;
        }
    }

    // Average the efficiencies for each layer (16 tubes per layer)
    vector<vector<Efficiency>> average_efficiencies_all_layers;

    int num_layers = 9;

    for (int i = 0; i < efficiencies_all.size(); ++i)
    {

        vector<Efficiency> average_efficiencies_per_layer;
        for (int i = 0; i < num_layers; ++i)
        {
            average_efficiencies_per_layer.push_back(Efficiency{0, 0, 0, 0, 0});
        }

        for (int j = 0; j < 144; ++j)
        {
            int layer = j / 16;

            average_efficiencies_per_layer.at(layer).eff += efficiencies_all[i][j].eff;
            average_efficiencies_per_layer.at(layer).eff_err_plus += efficiencies_all[i][j].eff_err_plus * efficiencies_all[i][j].eff_err_plus;
            average_efficiencies_per_layer.at(layer).eff_err_minus += efficiencies_all[i][j].eff_err_minus * efficiencies_all[i][j].eff_err_minus;
            average_efficiencies_per_layer.at(layer).eff_sys_plus += efficiencies_all[i][j].eff_sys_plus * efficiencies_all[i][j].eff_sys_plus;
            average_efficiencies_per_layer.at(layer).eff_sys_minus += efficiencies_all[i][j].eff_sys_minus * efficiencies_all[i][j].eff_sys_minus;
        }

        for (int i = 0; i < num_layers; ++i)
        {
            average_efficiencies_per_layer.at(i).eff /= 16;
            average_efficiencies_per_layer.at(i).eff_err_plus = sqrt(average_efficiencies_per_layer.at(i).eff_err_plus / 16);
            average_efficiencies_per_layer.at(i).eff_err_minus = sqrt(average_efficiencies_per_layer.at(i).eff_err_minus / 16);
            average_efficiencies_per_layer.at(i).eff_sys_plus = sqrt(average_efficiencies_per_layer.at(i).eff_sys_plus / 16);
            average_efficiencies_per_layer.at(i).eff_sys_minus = sqrt(average_efficiencies_per_layer.at(i).eff_sys_minus / 16);
        }

        average_efficiencies_all_layers.push_back(average_efficiencies_per_layer);
    }

    // Graph the average efficiencies for each chamber across all files. So each graph is for a different chamber, and each point is from a different file.
    vector<vector<double>> x_values;
    vector<vector<double>> y_values;
    vector<vector<double>> x_err_values;
    vector<vector<double>> y_err_positive_values;
    vector<vector<double>> y_err_negative_values;

    // sqrt( [eps(1-eps)]/N )

    for (int i = 0; i < average_efficiencies_all_chambers.size(); ++i)
    {
        vector<double> x;
        vector<double> y;
        vector<double> x_err;
        vector<double> y_err_positive;
        vector<double> y_err_negative;

        for (int j = 0; j < average_efficiencies_all_chambers.at(i).size(); ++j)
        {
            x.push_back(voltages.at(i));
            y.push_back(average_efficiencies_all_chambers.at(i).at(j).eff);
            x_err.push_back(0.01);

            double y_e_p = average_efficiencies_all_chambers.at(i).at(j).eff - sqrt(average_efficiencies_all_chambers.at(i).at(j).eff_err_plus * average_efficiencies_all_chambers.at(i).at(j).eff_err_plus + average_efficiencies_all_chambers.at(i).at(j).eff_sys_plus * average_efficiencies_all_chambers.at(i).at(j).eff_sys_plus) / 2;
            //y_err_positive.push_back( abs(y_e_p) );
            y_err_positive.push_back( average_efficiencies_all_chambers.at(i).at(j).eff_err_plus );
//TODO:
            double y_e_n = average_efficiencies_all_chambers.at(i).at(j).eff - sqrt(average_efficiencies_all_chambers.at(i).at(j).eff_err_minus * average_efficiencies_all_chambers.at(i).at(j).eff_err_minus + average_efficiencies_all_chambers.at(i).at(j).eff_sys_minus * average_efficiencies_all_chambers.at(i).at(j).eff_sys_minus) /  2;
            //y_err_negative.push_back( abs(y_e_n) );
            y_err_negative.push_back( average_efficiencies_all_chambers.at(i).at(j).eff_err_plus );

        }

        x_values.push_back(x);
        y_values.push_back(y);
        x_err_values.push_back(x_err);
        y_err_positive_values.push_back(y_err_positive);
        y_err_negative_values.push_back(y_err_negative);
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000);
    c1->Divide(3);

    int num_chambers = 3;


    for (int i = 0; i < num_chambers; ++i)
    {
        c1->cd(i + 1);
        vector<double> x;
        vector<double> y;
        vector<double> x_err;
        vector<double> y_err_positive;
        vector<double> y_err_negative;

        for (int j = 0; j < analysisfiles.size(); ++j)
        {
            x.push_back(x_values.at(j).at(i));
            y.push_back(y_values.at(j).at(i));
            x_err.push_back(x_err_values.at(j).at(i));
            y_err_positive.push_back(y_err_positive_values.at(j).at(i));
            y_err_negative.push_back(y_err_negative_values.at(j).at(i));
        }

        TGraphAsymmErrors *gr = new TGraphAsymmErrors(x.size(), &x[0], &y[0], &x_err[0], &x_err[0], &y_err_negative[0], &y_err_positive[0]);
        string title = "Chamber " + to_string(i + 1) + ";Chamber 2 Voltage (V);Efficiency";
        gr->SetTitle(title.c_str());
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(1);
        gr->SetMarkerColor(kBlack);
        gr->Draw("AP");
    }

    c1->SaveAs("efficiencies_chambers.pdf");

    vector<vector<double>> x_values_layer;
    vector<vector<double>> y_values_layer;
    vector<vector<double>> x_err_values_layer;
    vector<vector<double>> y_err_positive_values_layer;
    vector<vector<double>> y_err_negative_values_layer;

    for(int i = 0; i < average_efficiencies_all_layers.size(); ++i)
    {
        vector<double> x;
        vector<double> y;
        vector<double> x_err;
        vector<double> y_err_positive;
        vector<double> y_err_negative;

        for(int j = 0; j < average_efficiencies_all_layers.at(i).size(); ++j)
        {
            x.push_back(voltages.at(i));
            y.push_back(average_efficiencies_all_layers.at(i).at(j).eff);
            x_err.push_back(0.01);

            double y_e_p = average_efficiencies_all_layers.at(i).at(j).eff - sqrt(average_efficiencies_all_layers.at(i).at(j).eff_err_plus * average_efficiencies_all_layers.at(i).at(j).eff_err_plus + average_efficiencies_all_layers.at(i).at(j).eff_sys_plus * average_efficiencies_all_layers.at(i).at(j).eff_sys_plus) / 2;
            //y_err_positive.push_back( abs(y_e_p) );
            y_err_positive.push_back ( average_efficiencies_all_layers.at(i).at(j).eff_err_plus );

            double y_e_n = average_efficiencies_all_layers.at(i).at(j).eff - sqrt(average_efficiencies_all_layers.at(i).at(j).eff_err_minus * average_efficiencies_all_layers.at(i).at(j).eff_err_minus + average_efficiencies_all_layers.at(i).at(j).eff_sys_minus * average_efficiencies_all_layers.at(i).at(j).eff_sys_minus) /  2;
            //y_err_negative.push_back( abs(y_e_n) );
            y_err_negative.push_back( average_efficiencies_all_layers.at(i).at(j).eff_err_plus );
        }

        x_values_layer.push_back(x);
        y_values_layer.push_back(y);
        x_err_values_layer.push_back(x_err);
        y_err_positive_values_layer.push_back(y_err_positive);
        y_err_negative_values_layer.push_back(y_err_negative);
    }

    // Graph the average efficiencies for each layer (9 layers) across all files. So each graph is for a different layer, and each point is from a different file.
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
            x.push_back(x_values_layer.at(j).at(i));
            y.push_back(y_values_layer.at(j).at(i));
            x_err.push_back(x_err_values_layer.at(j).at(i));
            y_err_positive.push_back(y_err_positive_values_layer.at(j).at(i));
            y_err_negative.push_back(y_err_negative_values_layer.at(j).at(i));
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
