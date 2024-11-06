#include "../../Grouper/include/Detectors.hh"
#include "CLHEP/Vector/ThreeVector.h"

using namespace CLHEP;

int Verbosee = 0;

TFile *ANALYSIS_File;

TTreeReader *Reader;
TTreeReaderValue<int> *Tree_PDG;
TTreeReaderValue<int> *Tree_Event;
TTreeReaderValue<double> *Tree_time;
TTreeReaderValue<double> *Tree_Ex;
TTreeReaderValue<double> *Tree_E0;
TTreeReaderValue<double> *Tree_px;
TTreeReaderValue<double> *Tree_py;
TTreeReaderValue<double> *Tree_pz;

map<int, TH1D *> H_E0;
map<int, TH1D *> H_Ex;
map<int, TH1D *> H_Time;
map<int, TH1D *> H_px;
map<int, TH1D *> H_py;
map<int, TH1D *> H_pz;
map<int , string> ParticleName;
map<int, TDirectory *> dir_Particle;

string SearchName(int PDG)
{
    if (PDG == 0)
    {
        return "All";
    }

    if (PDG > 1000000000)
    {
        PDG = (PDG/10)*10;
    }
    for (const auto &pair : NametoCode_map)
    {
        if (pair.second == PDG)
        {
            return pair.first;
        }
    }

    Error("No name found for PDG: " + to_string(PDG));
    return "";
}


void VerifyPDG(int PDG)
{
    if (H_E0.find(PDG) == H_E0.end())
    {
        ParticleName[PDG] = SearchName(PDG);

        dir_Particle[PDG] = ANALYSIS_File->mkdir(ParticleName[PDG].c_str());
        H_E0[PDG] = new TH1D(("H_E0_" + to_string(PDG)).c_str(), ("H_E0_" + to_string(PDG)).c_str(), 10000, 0, 10000);
        H_Ex[PDG] = new TH1D(("H_Ex_" + to_string(PDG)).c_str(), ("H_Ex_" + to_string(PDG)).c_str(), 10000, 0, 10000);
        H_Time[PDG] = new TH1D(("H_Time_" + to_string(PDG)).c_str(), ("H_Time_" + to_string(PDG)).c_str(), 1000, 0, 10000);
        H_px[PDG] = new TH1D(("H_px_" + to_string(PDG)).c_str(), ("H_px_" + to_string(PDG)).c_str(), 100, 0, 100);
        H_py[PDG] = new TH1D(("H_py_" + to_string(PDG)).c_str(), ("H_py_" + to_string(PDG)).c_str(), 100, 0, 100);
        H_pz[PDG] = new TH1D(("H_pz_" + to_string(PDG)).c_str(), ("H_pz_" + to_string(PDG)).c_str(), 100, 0, 100);

        
    }

    
}

void InitHistograms()
{
    Info("Initializing Histograms");

    // loop on map
    
}




void WriteHistograms()
{   

    Info("Writing Histograms");
    

    //adding for over all
    ANALYSIS_File->cd();

    for (auto &pair : H_E0)
    {
        int PDG = pair.first;

        dir_Particle[PDG]->cd();
        H_E0[PDG]->Write();
        H_Ex[PDG]->Write();
        H_Time[PDG]->Write();
        H_px[PDG]->Write();
        H_py[PDG]->Write();
        H_pz[PDG]->Write();
    }
}