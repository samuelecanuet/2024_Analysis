#ifndef MERGER_HH
#define MERGER_HH

#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <boost/filesystem.hpp>
#include <array>

#include "TFile.h"
#include <TStyle.h>
#include "TCanvas.h"
#include <TKey.h>
#include <TMath.h>
#include "TF1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TRandom.h"
#include "TPaveText.h"
#include "TLegend.h"

#include "TMinuit.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"

#include "../../../lib/SignalDict/Signal.h"
#include "Detectors.hh"

#include <thread>
#include <chrono>

using namespace std;
namespace fs = boost::filesystem;
using namespace ROOT::Math;
random_device rd;
mt19937 gen(rd());


TFile *MERGED_File;
TTree *MERGED_Tree;
TFile *MATCHED_File;
TFile *MATCHED_SiPM_FILE;
string GROUPED_filename;
TFile *GROUPED_File;
TF1* Matching_function[SIGNAL_MAX];
TF1* SiPM_ONOFF[SIGNAL_MAX];
TF1* SiPM_HighLow[SIGNAL_MAX];
TF1* SiPM_Gain[SIGNAL_MAX];
TDirectory *dir_det[SIGNAL_MAX];
TH1D *H[SIGNAL_MAX];
TGraph *G_RearStrip[SIGNAL_MAX];
TH2D *H_RearStrip[SIGNAL_MAX];
TH2D *H_RearStrip_Matched[SIGNAL_MAX];
TH1D *H_Strip_Matched[SIGNAL_MAX];
int counter_graph[SIGNAL_MAX] = {0};
TF1 *MatchingRearStrip[SIGNAL_MAX];

TH1D *H_SiPM_Left[SIGNAL_MAX];
TH1D *H_SiPM_Right[SIGNAL_MAX];
TH1D *H_SiPM_Center[SIGNAL_MAX];
TH1D *H_SiPM_Full[SIGNAL_MAX];
TH1D *H_SiPM_LeftLeft[SIGNAL_MAX];

TH1D* H_S[SIGNAL_MAX];
TGraph *G_N[SIGNAL_MAX];
TGraph *G_Fake;

vector<Signal> MERGED_Tree_Silicon;
vector<vector<pair<Signal, Signal>>> MERGED_Tree_SiPMGroup;
Signal MERGED_Tree_HRS;

TTreeReaderArray<Signal> *Silicon;
TTreeReaderValue<vector<vector<pair<Signal, Signal>>>> *SiPM_Groups;
TTreeReaderValue<Signal> *HRS;
double Channel;

TH1D* H_Sum;

// for ISOLDE proton pulses
map<string, double> Calibration_Shift_Proton_Pulses;
map<string, TTree*> Tree_ISOLDE_Proton; 
map<string, TTreeReader*> Reader_ISOLDE_Proton;
map<string, TTreeReaderValue<double>*> Reader_ISOLDE_Proton_Time;



void InitGraph()
{
  H_Sum = new TH1D("H_SuM", "H_SuM", winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
  for (int i = 0; i <= SIGNAL_MAX; i++)
  {
    if (IsDetectorSiliStrip(i))
    {
      G_RearStrip[i] = new TGraph();
      G_RearStrip[i]->SetName(("G_RearStrip_" + detectorName[i]).c_str());
      G_RearStrip[i]->SetTitle(("G_RearStrip_" + detectorName[i]).c_str());
      G_RearStrip[i]->GetXaxis()->SetTitle("Rear Channel");
      G_RearStrip[i]->GetYaxis()->SetTitle("Strip Channel");
      G_RearStrip[i]->GetXaxis()->CenterTitle();
      G_RearStrip[i]->GetYaxis()->CenterTitle();

      H_Strip_Matched[i] = new TH1D(("H_Strip_Matched_" + detectorName[i]).c_str(), ("H_Strip_Matched_" + detectorName[i]).c_str(), eSiliN/10, 0, eSiliMax);
      H_Strip_Matched[i]->GetXaxis()->SetTitle("Strip Channel");
      H_Strip_Matched[i]->GetYaxis()->SetTitle("Counts");
      H_Strip_Matched[i]->GetXaxis()->CenterTitle();
      H_Strip_Matched[i]->GetYaxis()->CenterTitle();
    }

    if (IsDetectorSiliBack(i))
    {
      H_RearStrip[i] = new TH2D(("H_RearStrip_" + detectorName[i]).c_str(), ("H_RearStrip_" + detectorName[i]).c_str(), eSiliN/10, 0, eSiliMax, eSiliN/10, 0, eSiliMax);
      H_RearStrip[i]->GetXaxis()->SetTitle("Rear Channel");
      H_RearStrip[i]->GetYaxis()->SetTitle("Strip Channel");
      H_RearStrip[i]->GetXaxis()->CenterTitle();
      H_RearStrip[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Matched[i] = new TH2D(("H_RearStrip_Matched_" + detectorName[i]).c_str(), ("H_RearStrip_Matched_" + detectorName[i]).c_str(), eSiliN/10, 0, eSiliMax, eSiliN/10, 0, eSiliMax);
      H_RearStrip_Matched[i]->GetXaxis()->SetTitle("Rear Channel");
      H_RearStrip_Matched[i]->GetYaxis()->SetTitle("Strip Channel");
      H_RearStrip_Matched[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Matched[i]->GetYaxis()->CenterTitle();

      
    }
  }

  for (int mul = 1; mul <= BETA_SIZE; mul++)
  {
    H_SiPM_Left[mul] = new TH1D(("H_SiPM_Left_" + to_string(mul)).c_str(), ("H_SiPM_Left_" + to_string(mul)).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
    H_SiPM_Left[mul]->GetXaxis()->SetTitle("Time (ns)");
    H_SiPM_Left[mul]->GetYaxis()->SetTitle("Counts");
    H_SiPM_Left[mul]->GetXaxis()->CenterTitle();
    H_SiPM_Left[mul]->GetYaxis()->CenterTitle();

    H_SiPM_Right[mul] = new TH1D(("H_SiPM_Right_" + to_string(mul)).c_str(), ("H_SiPM_Right_" + to_string(mul)).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
    H_SiPM_Right[mul]->GetXaxis()->SetTitle("Time (ns)");
    H_SiPM_Right[mul]->GetYaxis()->SetTitle("Counts");
    H_SiPM_Right[mul]->GetXaxis()->CenterTitle();
    H_SiPM_Right[mul]->GetYaxis()->CenterTitle();

    H_SiPM_Center[mul] = new TH1D(("H_SiPM_Center_" + to_string(mul)).c_str(), ("H_SiPM_Center_" + to_string(mul)).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
    H_SiPM_Center[mul]->GetXaxis()->SetTitle("Time (ns)");
    H_SiPM_Center[mul]->GetYaxis()->SetTitle("Counts");
    H_SiPM_Center[mul]->GetXaxis()->CenterTitle();
    H_SiPM_Center[mul]->GetYaxis()->CenterTitle();

    H_SiPM_Full[mul] = new TH1D(("H_SiPM_Full_" + to_string(mul)).c_str(), ("H_SiPM_Full_" + to_string(mul)).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
    H_SiPM_Full[mul]->GetXaxis()->SetTitle("Time (ns)");
    H_SiPM_Full[mul]->GetYaxis()->SetTitle("Counts");
    H_SiPM_Full[mul]->GetXaxis()->CenterTitle();
    H_SiPM_Full[mul]->GetYaxis()->CenterTitle();

    H_SiPM_LeftLeft[mul] = new TH1D(("H_SiPM_LeftLeft_" + to_string(mul)).c_str(), ("H_SiPM_LeftLeft_" + to_string(mul)).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
    H_SiPM_LeftLeft[mul]->GetXaxis()->SetTitle("Time (ns)");
    H_SiPM_LeftLeft[mul]->GetYaxis()->SetTitle("Counts");
    H_SiPM_LeftLeft[mul]->GetXaxis()->CenterTitle();
    H_SiPM_LeftLeft[mul]->GetYaxis()->CenterTitle();

    H_S[mul] = new TH1D(("H_S_" + to_string(mul)).c_str(), ("H_S_" + to_string(mul)).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
    H_S[mul]->GetXaxis()->SetTitle("Time (ns)");
    H_S[mul]->GetYaxis()->SetTitle("Counts");
    H_S[mul]->GetXaxis()->CenterTitle();
    H_S[mul]->GetYaxis()->CenterTitle();

    G_N[mul] = new TGraph();
    
  }
  G_Fake = new TGraph();
}

void InitProtonPulse()
{
  if (YEAR != 2024)
    return;

  // load ISOLDE proton pulses
  TFile *fISOLDE = MyTFile(DIR_DATA_ISOLDE + "ISOLDE_pulses.root", "READ");
  for (const auto &pair : Map_RunFiles)
  {
    string NUCLEUS = pair.first;
    vector<string> Runs = pair.second;
    for (const auto &Run : Runs)
    {
      Tree_ISOLDE_Proton[Run] = (TTree*)fISOLDE->Get(("Tree_ISOLDE_Proton_" + Run).c_str());
      if (Tree_ISOLDE_Proton[Run] == NULL)
      {
        Error("Tree_ISOLDE_Proton_" + Run + " not found in ISOLDE_pulses.root");
      }
      else
      {
        Reader_ISOLDE_Proton[Run] = new TTreeReader(Tree_ISOLDE_Proton[Run]);
        Reader_ISOLDE_Proton_Time[Run] = new TTreeReaderValue<double>(*Reader_ISOLDE_Proton[Run], "Time");
        Reader_ISOLDE_Proton[Run]->Next();
        Info("Tree_ISOLDE_Proton_" + Run + " loaded successfully", 2);
      }
    }
  }

  // correction factor
  TFile *f = MyTFile((DIR_ROOT_DATA_MATCHED + "Time_Pulse_FitParameters_" + to_string(YEAR) + " (copy).root").c_str(), "READ");
  TGraphErrors *G_Calibration_Shift_Proton_Pulse = (TGraphErrors*)f->Get("G_Time_Pulse_Calibration");
  for (int i = 0; i < G_Calibration_Shift_Proton_Pulse->GetN(); i++)
  {
    double Run, Shift;
    G_Calibration_Shift_Proton_Pulse->GetPoint(i, Run, Shift);
    string Run_str = Run < 100 ? "0" + to_string((int)Run) : to_string((int)Run);
    Calibration_Shift_Proton_Pulses[Run_str] = Shift;
  }
  f->Close();
  
}

void LoadMatchingFunction(int Run)
{
  // Info("Loading matching functions for run " + to_string(Run));
  // MATCHING SILICON
  for (int i = 0; i <= SIGNAL_MAX; i++)
  {
    if (IsDetectorSiliStrip(i))
    {
      Matching_function[i] = (TF1*)MATCHED_File->Get(("poll1" + to_string(Run) + to_string(i)).c_str()); 

      if (Matching_function[i] == NULL)
      {
        int run_past = Run;
        // int run_future = Run;
        TF1 *fpast = NULL;
        // TF1 *future = NULL;
        // looking for past run 
        while (fpast == NULL)
        {
          run_past--;
          fpast = (TF1*)MATCHED_File->Get(("poll1" + to_string(run_past) + to_string(i)).c_str());
        } 
        // looking for future run
        // while (future == NULL)
        // {
        //   run_future++;
        //   future = (TF1*)MATCHED_File->Get(("poll1" + to_string(run_future) + to_string(i)).c_str());
        // }
        // average the two functions
        // double a1 = fpast->GetParameter(1);
        // double b1 = fpast->GetParameter(0);
        // double a2 = future->GetParameter(1);
        // double b2 = future->GetParameter(0);
        // double a = (a1 + a2) / 2;
        // double b = (b1 + b2) / 2;
        Matching_function[i] = (TF1*)fpast->Clone(("poll1" + to_string(Run) + to_string(i)).c_str());
        // Matching_function[i]->SetParameter(0, b);
        // Matching_function[i]->SetParameter(1, a);
        Warning("Matching function for detector " + detectorName[i] + " not found for run " + to_string(Run) + ", using averaged function from run " + to_string(run_past));

        // Matching_function[i] = new TF1(("poll1" + to_string(Run) + to_string(i)).c_str(), "x", 0, 10000); 
        // Warning("Matching function for detector " + detectorName[i] + " not found for run " + to_string(Run) + ", using default function");
      }
    }
  }

  // MATCHING SIPMs
  for (int det = 0; det < SIGNAL_MAX; det++)
  {
    if (IsDetectorBeta(det))
    {
      if (Run < 80)
        SiPM_ONOFF[det] = ((TF1 *)MATCHED_SiPM_FILE->Get(("SiPM_ONOFF_" + detectorName[det]).c_str()));
      else
        SiPM_ONOFF[det] = new TF1(("SiPM_ONOFF_" + detectorName[det]).c_str(), "x", 0, 100e6);
      if (SiPM_ONOFF[det] == nullptr)
      {
        Error("SiPM ONOFF Function not found for detector: " + detectorName[det]);
      }
    }
    if (IsDetectorBetaLow(det))
    {
      SiPM_HighLow[det] = ((TF1 *)MATCHED_SiPM_FILE->Get(("SiPM_HighLow_" + to_string(GetDetectorChannel(det))).c_str()));
      if (SiPM_HighLow[det] == nullptr)
      {
        Error("SiPM HighLow Function not found for detector: " + detectorName[det]);
      }
    }
    if (IsDetectorBeta(det))
    {
      SiPM_Gain[det] = InvertFunction((TF1 *)MATCHED_SiPM_FILE->Get(("SiPM_Gain_" + detectorName[det]).c_str()));
      if (SiPM_Gain[det] == nullptr)
      {
        Error("SiPM Gain Function not found for detector: " + detectorName[det]);
      }
    }
  }

  // Info("Matching functions loaded");
}

void InsertProtonPulse(string Run, double Time)
{
  if (YEAR != 2024)
    return;

  double ISOLDE_Time = (**Reader_ISOLDE_Proton_Time[Run] - Calibration_Shift_Proton_Pulses[Run] + 0.5) * 1e9;
  // cout << "ISOLDE Time : " << ISOLDE_Time*1e-9 << " s" << endl;
  // cout << "Current Time : " << Time*1e-9 << " s" << endl;
  if (ISOLDE_Time < Time)
  {
    Signal HRS = Signal(99, ISOLDE_Time, 0);
    MERGED_Tree_HRS = HRS;
    MERGED_Tree->Fill();
    
    // cout << "        Inserted proton pulse at " << ISOLDE_Time << " s" << endl;

    // next proton pulse
    Reader_ISOLDE_Proton[Run]->Next();

    // initialize the MERGED_Tree for current detector data
    MERGED_Tree_HRS = Signal();
    MERGED_Tree_SiPMGroup = vector<vector<pair<Signal, Signal>>>();
    MERGED_Tree_Silicon = vector<Signal>();

    InsertProtonPulse(Run, Time);
  }
}

#endif