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
        cout << "Error! No data filename specfied (in json format). Exiting..." << endl;
        return 1;
    }
    std::ifstream data_file(argv[1]);
    if (!data_file) {
        cout << "Error! File not found. Exiting..." << endl;
        return 2;
    }

    json data;
    data_file >> data;

    for (auto & [name, yields] : data.items()) {

        float A_sig = yields["A_sig"].get<float>();
        float B_sig = yields["B_sig"].get<float>();
        float C_sig = yields["C_sig"].get<float>();
        float D_sig = yields["D_sig"].get<float>();

        float A_bkg = yields["A_bkg"].get<float>();
        float B_bkg = yields["B_bkg"].get<float>();
        float C_bkg = yields["C_bkg"].get<float>();
        float D_bkg = yields["D_bkg"].get<float>();

        float c1 = B_bkg / A_bkg;
        float c2 = C_bkg / A_bkg;

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
            {"A", A_bkg},
            {"B", B_bkg},
            {"C", C_bkg},
            {"D", D_bkg}
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

        cb.cp().ForEachObs([&](ch::Observation *x) {
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
        cb.cp().backgrounds().bin({"A", "B", "C", "D"}).AddSyst(cb, "bkgA_norm", "rateParam", SystMap<>::init(A_bkg));
        cb.cp().backgrounds().bin({"B", "D"}).AddSyst(cb, "c1", "rateParam", SystMap<>::init(c1));
        cb.cp().backgrounds().bin({"C", "D"}).AddSyst(cb, "c2", "rateParam", SystMap<>::init(c2));

        //! [part3]

        //! [part4]
        //cb.PrintAll();

        cout << ">> Writing datacard for hist: " << name << "\n";
        mkdir(name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        cb.WriteDatacard(name + "/datacard.txt");
        //! [part4]
    }
}





