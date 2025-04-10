
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TFile.h"
#include "TH1D.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TString.h"

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <fstream> 
#include <dirent.h>
#include <random>
#include <ctime>
#include <string>
#include <algorithm>
#include <gsl/gsl_statistics.h>
#include <filesystem>

using namespace std;

int Verbosee = 0;

TFile *ANALYSIS_File;

TTreeReader *Reader;
TTreeReaderArray<int> *Tree_PDG;
TTreeReaderArray<int> *Tree_Event;
TTreeReaderArray<double> *Tree_time;
TTreeReaderArray<double> *Tree_Ex;
TTreeReaderArray<double> *Tree_E0;
TTreeReaderArray<double> *Tree_px;
TTreeReaderArray<double> *Tree_py;
TTreeReaderArray<double> *Tree_pz;

map<int, TH1D *> H_E0;
map<int, TH1D *> H_Ex;
map<int, TH1D *> H_Time;
map<int, TH1D *> H_px;
map<int, TH1D *> H_py;
map<int, TH1D *> H_pz;
map<int , string> ParticleName;
map<int, TDirectory *> dir_Particle;

void ProgressBar(int cEntry, int TotalEntries, clock_t start, clock_t Current, string Prefix = "")
{
  if (cEntry % 10000 == 0 && cEntry > 2 * 1000)
  {
    Current = clock();
    const Char_t *Color;
    Double_t Frac = 1.0 * cEntry / TotalEntries;
    Double_t Timeclock = ((double)(Current - start) / CLOCKS_PER_SEC);
    Double_t TimeLeft = Timeclock * (1 / Frac - 1.);
    Color = "\e[1;31m";

    cout << Form(("\r"+Prefix+" Entries : ").c_str())
         << TotalEntries
         << " --- "
         << Form("%4.2f", 100. * cEntry / TotalEntries) << " %"
         << " --- "
         << " Time Left : " << Form("%2d min ", (int)TimeLeft / 60)
         << Form("%02d sec", (int)TimeLeft % 60)
         << flush;
  }

  if (cEntry == TotalEntries-1)
  {
    Current = clock();
    const Char_t *Color;
    Double_t Frac = 1.0 * cEntry / TotalEntries;
    Double_t Timeclock = ((double)(Current - start) / CLOCKS_PER_SEC);
    Double_t TimeLeft = Timeclock * (1 / Frac - 1.);
    Color = "\e[1;31m";
    cout << Form(("\r"+Prefix+" Entries : ").c_str())
         << TotalEntries
         << " --- "
         << Form("%4.2f", 100. * cEntry / TotalEntries) << " %"
         << " --- "
         << " Time Left : " << Form("%2d min ", (int)TimeLeft / 60)
         << Form("%02d sec", (int)TimeLeft % 60)
         << flush;
    cout << endl;
  }
}

int CRADLE_READER()
{
    string name = "/run/media/local1/Disque_Dur/2024_DATA/SIMULATED_DATA/04-01/test";
    // name = "241Am_700nm_width";
    TFile *SIMULATED_File = new TFile((name + ".root").c_str(), "READ");
    ANALYSIS_File = new TFile((name + "_analysed.root").c_str(), "RECREATE");
    TTree *tree = (TTree *)SIMULATED_File->Get("ParticleTree");

    Reader = new TTreeReader(tree);
    Tree_PDG = new TTreeReaderArray<int>(*Reader, "code");
    Tree_E0 = new TTreeReaderArray<double>(*Reader, "energy");
    Tree_time = new TTreeReaderArray<double>(*Reader, "time");
    Tree_Ex = new TTreeReaderArray<double>(*Reader, "excitation_energy");
    Tree_px = new TTreeReaderArray<double>(*Reader, "px");
    Tree_py = new TTreeReaderArray<double>(*Reader, "py");
    Tree_pz = new TTreeReaderArray<double>(*Reader, "pz");


    int Entries = tree->GetEntries();
    clock_t start = clock(), Current;

    int Event_MAX = 1e6;

    TH1D *H_E0[100];
    H_E0[11] = new TH1D("H_E0_e-", "H_E0_e-", 10000, 0, 10000);
    H_E0[11]->GetXaxis()->SetTitle("Energy [keV]");
    H_E0[11]->GetYaxis()->SetTitle("Counts");
    H_E0[11]->GetXaxis()->CenterTitle();
    H_E0[11]->GetYaxis()->CenterTitle();

    while (Reader->Next() && Reader->GetCurrentEntry() < Event_MAX)
    {
        ProgressBar(Reader->GetCurrentEntry(), Event_MAX, start, Current, "Reading Tree");

        // cout << "Event: " << **Tree_Event << endl;
        for (int part_i = 0; part_i < Tree_PDG->GetSize(); part_i++)
        {
            if (Tree_PDG->At(part_i) == 11)
            {
                int PDG = Tree_PDG->At(part_i);
                // cout << "Particle: " << **Tree_PDG << " Energy: " << **Tree_E0 << " Ex: " << **Tree_Ex << " Time: " << **Tree_time << " px: " << **Tree_px << " py: " << **Tree_py << " pz: " << **Tree_pz << endl;
                H_E0[PDG]->Fill(Tree_E0->At(part_i));
                // H_Ex[PDG]->Fill(Tree_Ex->At(part_i));
                // H_Time[PDG]->Fill(Tree_time->At(part_i));
                // H_px[PDG]->Fill(Tree_px->At(part_i));
                // H_py[PDG]->Fill(Tree_py->At(part_i));
                // H_pz[PDG]->Fill(Tree_pz->At(part_i));

                break; // taking only the first e-
            }
        }
    }

    H_E0[11]->Write();


    SIMULATED_File->Close();
    ANALYSIS_File->Close();
}