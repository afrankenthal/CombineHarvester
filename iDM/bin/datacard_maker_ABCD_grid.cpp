#include <TCanvas.h>
#include <TKey.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>

#include <string>
#include <sys/stat.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "utils/json.hpp"
using json = nlohmann::json;

using namespace std;

int main(int argc, char * argv[]) {
    if (argc == 1) {
        cout << "Error! No config filename specfied (in json format). Exiting..." << endl;
        return 1;
    }
    std::ifstream cfg_file(argv[1]);
    if (!cfg_file) {
        cout << "Error! File not found. Exiting..." << endl;
        return 2;
    }

    vector<float> grid_x{0.2, 0.6, 1.0};
    vector<float> grid_y{400, 350, 300, 250, 200, 150, 100};

    json cfg;
    cfg_file >> cfg;

    TString in_filename_signal = TString(cfg["infilename_signal"].get<std::string>());
    TString in_filename_VR = TString(cfg["infilename_VR"].get<std::string>());
    TString in_filename_template = TString(cfg["infilename_template"].get<std::string>());
    TString in_filename_SR = TString(cfg["infilename_SR"].get<std::string>());
    TFile * in_file_signal = new TFile(in_filename_signal, "READ");
    TFile * in_file_VR = new TFile(in_filename_VR, "READ");
    TFile * in_file_template = new TFile(in_filename_template, "READ");
    TFile * in_file_SR = new TFile(in_filename_SR, "READ");

    ofstream out_file;
    out_file.open("discovery_signif_scan.txt");

    for (auto && keyAsObj : *in_file_VR->GetListOfKeys()) {
        auto key = (TKey*)keyAsObj;
        if (TString(key->GetClassName()) != "TCanvas") continue;
        TString canvas_name = TString(key->GetName());
        if (!canvas_name.Contains("36")) continue;
        if (!canvas_name.Contains("canvas2D")) continue;
        if (!canvas_name.Contains("-DATA")) continue;
        cout << "Processing " << canvas_name << ", class " << key->GetClassName() << endl;

        TString h_name = canvas_name;
        h_name.Remove(0, 9); // remove "Canvas2D_" to get TH2D name
        cout << h_name << endl;
        TCanvas * c_VR = (TCanvas*)in_file_VR->Get(canvas_name);
        TCanvas * c_SR = (TCanvas*)in_file_SR->Get(canvas_name);
        TH2D * h_VR = (TH2D*)c_VR->FindObject(h_name);
        TH2D * h_SR = (TH2D*)c_SR->FindObject(h_name);

        TH1D * template_px = (TH1D*)in_file_template->Get(h_name + TString("_px"));
        TH1D * template_py = (TH1D*)in_file_template->Get(h_name + TString("_py"));

        template_px->Scale(1/template_px->Integral());
        template_py->Scale(1/template_py->Integral());

        canvas_name.ReplaceAll("-DATA", "_sig_52p5_100");
        h_name = canvas_name;
        h_name.Remove(0, 9);
        TCanvas * c_sig = (TCanvas*)in_file_signal->Get(canvas_name);
        TH2D * h_sig = (TH2D*)c_sig->FindObject(h_name);
        //TH2D * h_map_sig = (TH2D*)h_sig->Clone();
        //h_map_sig->Reset();
        //h_map_sig->SetName("SensitivityMap_" + h_name);

        for (auto y : grid_y) {
            int bin_y = h_VR->GetYaxis()->FindBin(y);
            cout << "biny " << bin_y << endl;

            for (auto x : grid_x) {
                int bin_x = h_VR->GetXaxis()->FindBin(x);
                cout << "binx " << bin_x << endl;

                float A_SR = h_SR->Integral(0, bin_x-1, 0, bin_y-1);
                // Template method
                float c1 = template_px->Integral(bin_x, template_px->GetNbinsX()+1);
                float c2 = template_py->Integral(bin_y, template_py->GetNbinsX()+1);
                float B_SR_pred = A_SR * c1;
                float C_SR_pred = A_SR * c2;
                float D_SR_pred = A_SR * c1 * c2;

                // Closure error (estimate)
                float C_SR_clos_err = C_SR_pred * 0.14;
                float C_SR_err = C_SR_clos_err;

                float A_sig = h_sig->Integral(0, bin_x-1, 0, bin_y-1);
                float B_sig = h_sig->Integral(bin_x, h_sig->GetNbinsX()+1, 0, bin_y-1);
                float C_sig = h_sig->Integral(0, bin_x-1, bin_y, h_sig->GetNbinsY()+1);
                float D_sig = h_sig->Integral(bin_x, h_sig->GetNbinsX()+1, bin_y, h_sig->GetNbinsY()+1);

                // rescale to 2018+2017+2016 (137 / 60 = 2.30)
                float sig_sf = 2.30;
                A_sig *= sig_sf;
                B_sig *= sig_sf;
                C_sig *= sig_sf;
                D_sig *= sig_sf;

                // Compute asymptotic discovery significance with formula
                auto calc_ZA = [&](float s, float b, float b_err) {
                    return sqrt(2.0*((s+b)*log((s+b)*(b+b_err*b_err)/(b*b+(s+b)*b_err*b_err)) - (b*b/(b_err*b_err))*log(1+(b_err*b_err*s)/(b*(b+b_err*b_err)))));
                };
                float ZA_sig = calc_ZA(C_sig, C_SR_pred, C_SR_err);

                cout << "x = " << x << ", y = " << y << endl;
                cout << "Asymptotic discovery sign. = " << ZA_sig << endl;
                out_file << ZA_sig << " ";
                cout << "Background in A, B, C, D (predicted) = " << A_SR << ", " << B_SR_pred << ", " << C_SR_pred << ", " << D_SR_pred << endl;
                cout << "Signal yields in A, B, C, D = " << A_sig << ", " << B_sig << ", " << C_sig << ", " << D_sig << endl;

                // Now create combine datacard to compute asymptotic exclusion limit

                //! [part1]
                // Define four categories labelled A, B, C and D, and
                // set the observed yields in a map.
                ch::Categories cats = {
                    {0, "A"},
                    {1, "B"},
                    {2, "C"},
                    {3, "D"}
                };
                std::map<std::string, float> obs_rates = {
                    {"A", A_SR},
                    {"B", B_SR_pred},
                    {"C", C_SR_pred},
                    {"D", D_SR_pred}
                };
                std::map<std::string, float> sig_rates = {
                    {"A", A_sig},
                    {"B", B_sig},
                    {"C", C_sig},
                    {"D", D_sig}
                };
                //! [part1]

                //! [part2]
                ch::CombineHarvester cb;
                cb.SetVerbosity(0);

                cb.AddObservations({"*"}, {""}, {"13TeV"}, {""},          cats);
                cb.AddProcesses(   {"*"}, {""}, {"13TeV"}, {""}, {"sig"}, cats, true);
                cb.AddProcesses(   {"*"}, {""}, {"13TeV"}, {""}, {"bkg"}, cats, false);

                cb.ForEachObs([&](ch::Observation *x) {
                        x->set_rate(obs_rates[x->bin()]);
                        });
                cb.cp().backgrounds().ForEachProc([](ch::Process *x) {
                        x->set_rate(1);
                        });
                cb.cp().signals().ForEachProc([&](ch::Process *x) {
                        x->set_rate(sig_rates[x->bin()]);
                        });
                //! [part2]

                //! [part3]
                using ch::syst::SystMap;
                using ch::syst::SystMapFunc;
                using ch::syst::bin;

                // Add a traditional lnN systematic
                //cb.cp().bin({"D"}).AddSyst(cb, "DummySys", "lnN", SystMap<>::init(1.0001));

                // Create a unqiue floating parameter in each bin
                cb.cp().backgrounds().bin({"A", "B", "C", "D"}).AddSyst(cb, "bkgA_norm", "rateParam", SystMap<>::init(A_SR));
                cb.cp().backgrounds().bin({"B", "D"}).AddSyst(cb, "c1", "rateParam", SystMap<>::init(c1));
                cb.cp().backgrounds().bin({"C", "D"}).AddSyst(cb, "c2", "rateParam", SystMap<>::init(c2));

                //! [part3]

                //! [part4]
                //cb.PrintAll();
                std::stringstream ss;
                ss << "datacard_" << std::fixed << std::setprecision(1) << x << std::setprecision(0) << "_" << y << ".txt";
                cb.WriteDatacard(ss.str());
                //! [part4]
            }
            out_file << "\n";
        }
    }

    in_file_signal->Close();
    in_file_VR->Close();
    in_file_template->Close();
    in_file_SR->Close();
    out_file.close();
}





