#ifndef GROUPER_HH
#define GROUPER_HH

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <map>

#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TRandom.h"
#include "TLatex.h"

#include "../../../lib/SignalDict/Signal.h"
#include "Detectors.hh"

#include "/home/local1/Documents/lib/GTools1.0/include/GString.hh"

using namespace std;

int Verbose = 0;
int Run;

pair<double, double> detectorCleaning[SIGNAL_MAX];

double Event = 0;
int counterraph[BETA_SIZE + 1];

string ROOT_filename;
string ROOT_basefilename;

TFile *GROUPED_File;
TTree *GROUPED_Tree;

vector<Signal> Tree_Silicon;
vector<Signal> Tree_SiPMHigh;
vector<Signal> Tree_SiPMLow;

//////////////GROUPED////////////////
TH2D *hist_channel[SIGNAL_MAX];
TH2D *hist_label;
TH1D *hist_time[SIGNAL_MAX];
/// FOR ALL
TH1D *H_Channel_RAW[SIGNAL_MAX];
TH1D *H_Channel_ASSOCIATED[SIGNAL_MAX];
TH1D *H_Channel_A[SIGNAL_MAX];
TH1D *H_Channel_B[SIGNAL_MAX];
TH1D *H_Channel_C[SIGNAL_MAX];
TH1D *H_Channel_D[SIGNAL_MAX];
TH1D *H_Channel_E[SIGNAL_MAX];
TH1D *H_Channel_Cutted[SIGNAL_MAX];

/// Silicon
TH1D *H_Strip_Mulitplicity_RAW;
TH1D *H_Strip_Mulitplicity_ASSOCIATED;

/// Rear
TH1D *H_Rear_Mulitplicity_RAW;
TH1D *H_Rear_Mulitplicity_ASSOCIATED;

/// Rear-Strip
TCanvas *C_Strip_Channel_DiffRear[SIGNAL_MAX];
TGraphErrors *spread_up[SIGNAL_MAX];
TGraphErrors *spread_down[SIGNAL_MAX];
TH2D *H_RearStrip_ASSOCIATED;
TH2D *H_RearStrip_NON_ASSOCIATED;
TH2D *H_RearStrip_Channel_ASSOCIATED[SIGNAL_MAX];
TH2D *H_RearStrip_Channel_A[SIGNAL_MAX];
TH2D *H_RearStrip_Channel_C[SIGNAL_MAX];
TH2D *H_RearStrip_Channel_D[SIGNAL_MAX];
TH2D *H_RearStrip_Channel_E;
TH2D *H_RearStrip_Channel_Cutted[SIGNAL_MAX];
TH2D *H_Rear2Strip_Channel_C[SIGNAL_MAX];
TH2D *H_Rear2Strip_Channel_D[SIGNAL_MAX];
TH1D *H_Strip_Channel_DiffRear_A[SIGNAL_MAX];
TH1D *H_Strip_Channel_DiffRear_Cutted[SIGNAL_MAX];
TH2D *H_Strip_Channel_DiffRearvsStrip_A[SIGNAL_MAX];
TH2D *H_Strip_Channel_DiffRearvsStrip_Cutted[SIGNAL_MAX];
TH2D *H_2Strip_Channel_E;
TH1D *H_RearStrip_Time_ASSOCIATED[SIGNAL_MAX];
TH1D *H_RearStrip_Time_A[SIGNAL_MAX];
TH1D *H_RearStrip_Time_B[SIGNAL_MAX];
TH1D *H_RearStrip_Time_C[SIGNAL_MAX];
TH1D *H_RearStrip_Time_D[SIGNAL_MAX];
TH1D *H_RearStrip_Time_E;
TH1D *H_RearStrip_Time_Cutted[SIGNAL_MAX];
TH1D *H_2Strip_Time_E;
TH2D *H_2Strip_Label_E;

/// SiPM High
TH1D *H_SiPMHigh_Mulitplicity_RAW;
TH1D *H_SiPMHigh_Mulitplicity_Cutted;

/// SiPM Low
TH1D *H_SiPMLow_Mulitplicity_RAW;
TH1D *H_SiPMLow_Mulitplicity_Cutted;

TDirectory *dir_Strip;
TDirectory *dir_Strip_Channel;
TDirectory *dir_Strip_Channel_Write;
TDirectory *dir_Strip_Time;
TDirectory *dir_Strip_Time_Write;
TDirectory *dir_Strip_Channel_Time;
TDirectory *dir_Strip_Channel_DiffRear;
TDirectory *dir_Strip_Multiplicities;

TDirectory *dir_Rear;
TDirectory *dir_Rear_Channel;
TDirectory *dir_Rear_Channel_Write;
TDirectory *dir_Rear_Multiplicity;

TDirectory *dir_Strip_Rear;
TDirectory *dir_Rear_Strip_Channel;
TDirectory *dir_Strip_Rear_Time;

TDirectory *dir_SiPM_High;
TDirectory *dir_SiPMHigh_Channel;
TDirectory *dir_SiPMHigh_Channel_Write;
TDirectory *dir_SiPMHigh_Time;
TDirectory *dir_SiPMHighRear_TimeChannel;
TDirectory *dir_SiPMHigh_Multiplicity;

TDirectory *dir_SiPMLow;
TDirectory *dir_SiPMLow_Channel;
TDirectory *dir_SiPMLow_Channel_Write;
TDirectory *dir_SiPMLow_Time;
TDirectory *dir_SiPMLowRear_TimeChannel;
TDirectory *dir_SiPMLow_Multiplicity;

TDirectory *dir_SiPM_Multiplicities;

//////////////////////////////////////
/// CUTTING THE SIGNALS IN GROUP /////
//////////////////////////////////////
TF1 *SpreadAcceptance[SIGNAL_MAX];
double MeanAcceptance[SIGNAL_MAX];

//////////////////////////////////////

inline int InitHistograms_Grouped()
{
  dir_Strip = GROUPED_File->mkdir("Strip");
  dir_Strip_Channel = dir_Strip->mkdir("Strip_Channel");
  dir_Strip_Channel_Write = dir_Strip_Channel->mkdir("Strip_Channel_Write");
  dir_Strip_Time = dir_Strip->mkdir("Strip_Time");
  dir_Strip_Time_Write = dir_Strip_Time->mkdir("Strip_Time_Write");
  dir_Strip_Channel_Time = dir_Strip->mkdir("Strip_Channel_Time");
  dir_Strip_Channel_DiffRear = dir_Strip->mkdir("Strip_Channel_DiffRear");
  dir_Strip_Multiplicities = dir_Strip->mkdir("Strip_Multiplicities");

  ////REAR
  dir_Rear = GROUPED_File->mkdir("Rear");
  dir_Rear_Channel = dir_Rear->mkdir("Rear_Channel");
  dir_Rear_Channel_Write = dir_Rear_Channel->mkdir("Rear_Channel_Write");
  dir_Rear_Multiplicity = dir_Rear->mkdir("Rear_Multiplicity");

  /// STRIP REAR
  dir_Strip_Rear = GROUPED_File->mkdir("Strip_Rear");
  dir_Rear_Strip_Channel = dir_Strip_Rear->mkdir("Rear_Strip_Channel");
  dir_Strip_Rear_Time = dir_Strip_Rear->mkdir("Strip_Rear_Time");

  ////SiPM High
  dir_SiPM_High = GROUPED_File->mkdir("SiPM_High");
  dir_SiPMHigh_Channel = dir_SiPM_High->mkdir("SiPMHigh_Channel");
  dir_SiPMHigh_Channel_Write = dir_SiPM_High->mkdir("SiPMHigh_Channel_Write");
  dir_SiPMHigh_Time = dir_SiPM_High->mkdir("SiPMHigh_Time");
  dir_SiPMHighRear_TimeChannel = dir_SiPM_High->mkdir("SiPMHighRear_TimeChannel");
  dir_SiPMHigh_Multiplicity = dir_SiPM_High->mkdir("SiPMHigh_Multiplicity");

  ////SiPM Low
  dir_SiPMLow = GROUPED_File->mkdir("SiPMLow");
  dir_SiPMLow_Channel = dir_SiPMLow->mkdir("SiPMLow_Channel");
  dir_SiPMLow_Channel_Write = dir_SiPMLow_Channel->mkdir("SiPMLow_Channel_Write");
  dir_SiPMLow_Time = dir_SiPMLow->mkdir("SiPMLow_Time");
  dir_SiPMLowRear_TimeChannel = dir_SiPMLow->mkdir("SiPMLowRear_TimeChannel");
  dir_SiPMLow_Multiplicity = dir_SiPMLow->mkdir("SiPMLow_Multiplicity");

  ////SiPM
  dir_SiPM_Multiplicities = GROUPED_File->mkdir("SiPM_Multiplicities");

  for (size_t i = 0; i < detectorNum; ++i)
  {
    if (IsDetectorSiliStrip(i))
    {
      H_Channel_RAW[i] = new TH1D(("Channel_RAW_" + detectorName[i]).c_str(), ("Channel_RAW_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_RAW[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_RAW[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_RAW[i]->GetXaxis()->CenterTitle();
      H_Channel_RAW[i]->GetYaxis()->CenterTitle();

      H_Channel_ASSOCIATED[i] = new TH1D(("Channel_ASSOCIATED_" + detectorName[i]).c_str(), ("Channel_ASSOCIATED_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_ASSOCIATED[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_ASSOCIATED[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_ASSOCIATED[i]->GetXaxis()->CenterTitle();
      H_Channel_ASSOCIATED[i]->GetYaxis()->CenterTitle();

      H_Channel_A[i] = new TH1D(("Channel_A_" + detectorName[i]).c_str(), ("Channel_A_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_A[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_A[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_A[i]->GetXaxis()->CenterTitle();
      H_Channel_A[i]->GetYaxis()->CenterTitle();

      H_Channel_B[i] = new TH1D(("Channel_B_" + detectorName[i]).c_str(), ("Channel_B_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_B[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_B[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_B[i]->GetXaxis()->CenterTitle();
      H_Channel_B[i]->GetYaxis()->CenterTitle();

      H_Channel_C[i] = new TH1D(("Channel_C_" + detectorName[i]).c_str(), ("Channel_C_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_C[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_C[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_C[i]->GetXaxis()->CenterTitle();
      H_Channel_C[i]->GetYaxis()->CenterTitle();

      H_Channel_D[i] = new TH1D(("Channel_D_" + detectorName[i]).c_str(), ("Channel_D_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_D[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_D[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_D[i]->GetXaxis()->CenterTitle();
      H_Channel_D[i]->GetYaxis()->CenterTitle();

      H_Channel_E[i] = new TH1D(("Channel_E_" + detectorName[i]).c_str(), ("Channel_E_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_E[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_E[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_E[i]->GetXaxis()->CenterTitle();
      H_Channel_E[i]->GetYaxis()->CenterTitle();

      H_Channel_Cutted[i] = new TH1D(("Channel_Cutted_" + detectorName[i]).c_str(), ("Channel_Cutted_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_Cutted[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_Cutted[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_Cutted[i]->GetXaxis()->CenterTitle();
      H_Channel_Cutted[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Channel_A[i] = new TH2D(("RearStrip_Channel_A_" + detectorName[i]).c_str(), ("RearStrip_Channel_A_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearStrip_Channel_A[i]->GetXaxis()->SetTitle("Rear Channel");
      H_RearStrip_Channel_A[i]->GetYaxis()->SetTitle("Strip Channel");
      H_RearStrip_Channel_A[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Channel_A[i]->GetYaxis()->CenterTitle();

      H_Strip_Channel_DiffRear_A[i] = new TH1D(("Strip_Channel_DiffRear_A_" + detectorName[i]).c_str(), ("Strip_Channel_DiffRear_A_" + detectorName[i]).c_str(), 2000, -1, 1);
      H_Strip_Channel_DiffRear_A[i]->GetXaxis()->SetTitle("#DeltaE/E");
      H_Strip_Channel_DiffRear_A[i]->GetYaxis()->SetTitle("Counts");
      H_Strip_Channel_DiffRear_A[i]->GetXaxis()->CenterTitle();
      H_Strip_Channel_DiffRear_A[i]->GetYaxis()->CenterTitle();

      H_Strip_Channel_DiffRearvsStrip_A[i] = new TH2D(("Strip_Channel_DiffRearvsStrip_A_" + detectorName[i]).c_str(), ("Strip_Channel_DiffRearvsStrip_A_" + detectorName[i]).c_str(), 2000, -1, 1, eSiliN / 10, eSiliMin, eSiliMax);
      H_Strip_Channel_DiffRearvsStrip_A[i]->GetXaxis()->SetTitle("#DeltaE/E");
      H_Strip_Channel_DiffRearvsStrip_A[i]->GetYaxis()->SetTitle("Strip Channel");
      H_Strip_Channel_DiffRearvsStrip_A[i]->GetXaxis()->CenterTitle();
      H_Strip_Channel_DiffRearvsStrip_A[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Channel_Cutted[i] = new TH2D(("RearStrip_Channel_Cutted_" + detectorName[i]).c_str(), ("RearStrip_Channel_Cutted_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearStrip_Channel_Cutted[i]->GetXaxis()->SetTitle("Rear Channel");
      H_RearStrip_Channel_Cutted[i]->GetYaxis()->SetTitle("Strip Channel");
      H_RearStrip_Channel_Cutted[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Channel_Cutted[i]->GetYaxis()->CenterTitle();

      H_Strip_Channel_DiffRear_Cutted[i] = new TH1D(("Strip_Channel_DiffRear_Cutted_" + detectorName[i]).c_str(), ("Strip_Channel_DiffRear_Cutted_" + detectorName[i]).c_str(), 2000, -1, 1);
      H_Strip_Channel_DiffRear_Cutted[i]->GetXaxis()->SetTitle("#DeltaE/E");
      H_Strip_Channel_DiffRear_Cutted[i]->GetYaxis()->SetTitle("Counts");
      H_Strip_Channel_DiffRear_Cutted[i]->GetXaxis()->CenterTitle();
      H_Strip_Channel_DiffRear_Cutted[i]->GetYaxis()->CenterTitle();

      H_Strip_Channel_DiffRearvsStrip_Cutted[i] = new TH2D(("Strip_Channel_DiffRearvsStrip_Cutted_" + detectorName[i]).c_str(), ("Strip_Channel_DiffRearvsStrip_Cutted_" + detectorName[i]).c_str(), 2000, -1, 1, eSiliN / 10, eSiliMin, eSiliMax);
      H_Strip_Channel_DiffRearvsStrip_Cutted[i]->GetXaxis()->SetTitle("#DeltaE/E");
      H_Strip_Channel_DiffRearvsStrip_Cutted[i]->GetYaxis()->SetTitle("Strip Channel");
      H_Strip_Channel_DiffRearvsStrip_Cutted[i]->GetXaxis()->CenterTitle();
      H_Strip_Channel_DiffRearvsStrip_Cutted[i]->GetYaxis()->CenterTitle();

      hist_channel[i] = new TH2D(("hist_channel_" + detectorName[i]).c_str(), ("hist_channel_" + detectorName[i]).c_str(), eSiliN/10, eSiliMin, eSiliMax, eSiliN/10, eSiliMin, eSiliMax);
      hist_channel[i]->GetXaxis()->SetTitle("Channel (grouped)");
      hist_channel[i]->GetYaxis()->SetTitle("Channel");
      hist_channel[i]->GetXaxis()->CenterTitle();
      hist_channel[i]->GetYaxis()->CenterTitle();

      hist_time[i] = new TH1D(("hist_time_" + detectorName[i]).c_str(), ("hist_time_" + detectorName[i]).c_str(), winGroupN, winGroupMin, winGroupMax);
      hist_time[i]->GetXaxis()->SetTitle("Time (ns)");
      hist_time[i]->GetYaxis()->SetTitle("Counts");
      hist_time[i]->GetXaxis()->CenterTitle();
      hist_time[i]->GetYaxis()->CenterTitle();

    }

    if (IsDetectorSiliBack(i))
    {
      H_Channel_RAW[i] = new TH1D(("Channel_RAW_" + detectorName[i]).c_str(), ("Channel_RAW_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_RAW[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_RAW[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_RAW[i]->GetXaxis()->CenterTitle();
      H_Channel_RAW[i]->GetYaxis()->CenterTitle();

      H_Channel_ASSOCIATED[i] = new TH1D(("Channel_ASSOCIATED_" + detectorName[i]).c_str(), ("Channel_ASSOCIATED_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_ASSOCIATED[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_ASSOCIATED[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_ASSOCIATED[i]->GetXaxis()->CenterTitle();
      H_Channel_ASSOCIATED[i]->GetYaxis()->CenterTitle();

      H_Channel_A[i] = new TH1D(("Channel_A_" + detectorName[i]).c_str(), ("Channel_A_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_A[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_A[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_A[i]->GetXaxis()->CenterTitle();
      H_Channel_A[i]->GetYaxis()->CenterTitle();

      H_Channel_B[i] = new TH1D(("Channel_B_" + detectorName[i]).c_str(), ("Channel_B_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_B[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_B[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_B[i]->GetXaxis()->CenterTitle();
      H_Channel_B[i]->GetYaxis()->CenterTitle();

      H_Channel_C[i] = new TH1D(("Channel_C_" + detectorName[i]).c_str(), ("Channel_C_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_C[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_C[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_C[i]->GetXaxis()->CenterTitle();
      H_Channel_C[i]->GetYaxis()->CenterTitle();

      H_Channel_D[i] = new TH1D(("Channel_D_" + detectorName[i]).c_str(), ("Channel_D_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_D[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_D[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_D[i]->GetXaxis()->CenterTitle();
      H_Channel_D[i]->GetYaxis()->CenterTitle();

      H_Channel_E[i] = new TH1D(("Channel_E_" + detectorName[i]).c_str(), ("Channel_E_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_E[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_E[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_E[i]->GetXaxis()->CenterTitle();
      H_Channel_E[i]->GetYaxis()->CenterTitle();

      H_Channel_Cutted[i] = new TH1D(("Channel_Cutted_" + detectorName[i]).c_str(), ("Channel_Cutted_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_Cutted[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_Cutted[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_Cutted[i]->GetXaxis()->CenterTitle();
      H_Channel_Cutted[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Channel_ASSOCIATED[i] = new TH2D(("RearStrip_Channel_ASSOCIATED_" + detectorName[i]).c_str(), ("RearStrip_Channel_ASSOCIATED_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearStrip_Channel_ASSOCIATED[i]->GetXaxis()->SetTitle("Rear Channel");
      H_RearStrip_Channel_ASSOCIATED[i]->GetYaxis()->SetTitle("Strip Channel");
      H_RearStrip_Channel_ASSOCIATED[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Channel_ASSOCIATED[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Channel_A[i] = new TH2D(("RearStrip_Channel_A_" + detectorName[i]).c_str(), ("RearStrip_Channel_A_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearStrip_Channel_A[i]->GetXaxis()->SetTitle("Rear Channel");
      H_RearStrip_Channel_A[i]->GetYaxis()->SetTitle("Strip Channel");
      H_RearStrip_Channel_A[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Channel_A[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Channel_C[i] = new TH2D(("RearStrip_Channel_C_" + detectorName[i]).c_str(), ("RearStrip_Channel_C_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearStrip_Channel_C[i]->GetXaxis()->SetTitle("Rear Channel");
      H_RearStrip_Channel_C[i]->GetYaxis()->SetTitle("Strip Channel");
      H_RearStrip_Channel_C[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Channel_C[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Channel_D[i] = new TH2D(("RearStrip_Channel_D_" + detectorName[i]).c_str(), ("RearStrip_Channel_D_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearStrip_Channel_D[i]->GetXaxis()->SetTitle("Rear Channel");
      H_RearStrip_Channel_D[i]->GetYaxis()->SetTitle("Strip Channel");
      H_RearStrip_Channel_D[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Channel_D[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Channel_Cutted[i] = new TH2D(("RearStrip_Channel_Cutted_" + detectorName[i]).c_str(), ("RearStrip_Channel_Cutted_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearStrip_Channel_Cutted[i]->GetXaxis()->SetTitle("Rear Channel");
      H_RearStrip_Channel_Cutted[i]->GetYaxis()->SetTitle("Strip Channel");
      H_RearStrip_Channel_Cutted[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Channel_Cutted[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Time_ASSOCIATED[i] = new TH1D(("RearStrip_Time_ASSOCIATED_" + detectorName[i]).c_str(), ("RearStrip_Time_ASSOCIATED_" + detectorName[i]).c_str(), winGroupN, winGroupMin, winGroupMax);
      H_RearStrip_Time_ASSOCIATED[i]->GetXaxis()->SetTitle("Rear-Strip Time (ns)");
      H_RearStrip_Time_ASSOCIATED[i]->GetYaxis()->SetTitle("Counts");
      H_RearStrip_Time_ASSOCIATED[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Time_ASSOCIATED[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Time_A[i] = new TH1D(("RearStrip_Time_A_" + detectorName[i]).c_str(), ("RearStrip_Time_A_" + detectorName[i]).c_str(), winGroupN, winGroupMin, winGroupMax);
      H_RearStrip_Time_A[i]->GetXaxis()->SetTitle("Rear-Strip Time (ns)");
      H_RearStrip_Time_A[i]->GetYaxis()->SetTitle("Counts");
      H_RearStrip_Time_A[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Time_A[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Time_B[i] = new TH1D(("RearStrip_Time_B_" + detectorName[i]).c_str(), ("RearStrip_Time_B_" + detectorName[i]).c_str(), winGroupN, winGroupMin, winGroupMax);
      H_RearStrip_Time_B[i]->GetXaxis()->SetTitle("Rear-Strip Time (ns)");
      H_RearStrip_Time_B[i]->GetYaxis()->SetTitle("Counts");
      H_RearStrip_Time_B[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Time_B[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Time_C[i] = new TH1D(("RearStrip_Time_C_" + detectorName[i]).c_str(), ("RearStrip_Time_C_" + detectorName[i]).c_str(), winGroupN, winGroupMin, winGroupMax);
      H_RearStrip_Time_C[i]->GetXaxis()->SetTitle("Rear-Strip Time (ns)");
      H_RearStrip_Time_C[i]->GetYaxis()->SetTitle("Counts");
      H_RearStrip_Time_C[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Time_C[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Time_D[i] = new TH1D(("RearStrip_Time_D_" + detectorName[i]).c_str(), ("RearStrip_Time_D_" + detectorName[i]).c_str(), winGroupN, winGroupMin, winGroupMax);
      H_RearStrip_Time_D[i]->GetXaxis()->SetTitle("Rear-Strip Time (ns)");
      H_RearStrip_Time_D[i]->GetYaxis()->SetTitle("Counts");
      H_RearStrip_Time_D[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Time_D[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Time_Cutted[i] = new TH1D(("RearStrip_Time_Cutted_" + detectorName[i]).c_str(), ("RearStrip_Time_Cutted_" + detectorName[i]).c_str(), winGroupN, winGroupMin, winGroupMax);
      H_RearStrip_Time_Cutted[i]->GetXaxis()->SetTitle("Rear-Strip Time (ns)");
      H_RearStrip_Time_Cutted[i]->GetYaxis()->SetTitle("Counts");
      H_RearStrip_Time_Cutted[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Time_Cutted[i]->GetYaxis()->CenterTitle();

      H_Rear2Strip_Channel_C[i] = new TH2D(("Rear2Strip_Channel_C_" + detectorName[i]).c_str(), ("Rear2Strip_Channel_C_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_Rear2Strip_Channel_C[i]->GetXaxis()->SetTitle("Rear Channel");
      H_Rear2Strip_Channel_C[i]->GetYaxis()->SetTitle("2 Sum Strip Channel");
      H_Rear2Strip_Channel_C[i]->GetXaxis()->CenterTitle();
      H_Rear2Strip_Channel_C[i]->GetYaxis()->CenterTitle();

      H_Rear2Strip_Channel_D[i] = new TH2D(("Rear2Strip_Channel_D_" + detectorName[i]).c_str(), ("Rear2Strip_Channel_D_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_Rear2Strip_Channel_D[i]->GetXaxis()->SetTitle("Rear Channel");
      H_Rear2Strip_Channel_D[i]->GetYaxis()->SetTitle("2 Sum Strip Channel");
      H_Rear2Strip_Channel_D[i]->GetXaxis()->CenterTitle();
      H_Rear2Strip_Channel_D[i]->GetYaxis()->CenterTitle();
    }

    if (IsDetectorBetaHigh(i))
    {
      H_Channel_RAW[i] = new TH1D(("Channel_RAW_" + detectorName[i]).c_str(), ("Channel_RAW_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
      H_Channel_RAW[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_RAW[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_RAW[i]->GetXaxis()->CenterTitle();
      H_Channel_RAW[i]->GetYaxis()->CenterTitle();

      H_Channel_Cutted[i] = new TH1D(("Channel_Cutted_" + detectorName[i]).c_str(), ("Channel_Cutted_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
      H_Channel_Cutted[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_Cutted[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_Cutted[i]->GetXaxis()->CenterTitle();
      H_Channel_Cutted[i]->GetYaxis()->CenterTitle();
    }

    if (IsDetectorBetaLow(i))
    {
      H_Channel_RAW[i] = new TH1D(("Channel_RAW_" + detectorName[i]).c_str(), ("Channel_RAW_" + detectorName[i]).c_str(), eLowN, eLowMin, eLowMax);
      H_Channel_RAW[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_RAW[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_RAW[i]->GetXaxis()->CenterTitle();
      H_Channel_RAW[i]->GetYaxis()->CenterTitle();

      H_Channel_Cutted[i] = new TH1D(("Channel_Cutted_" + detectorName[i]).c_str(), ("Channel_Cutted_" + detectorName[i]).c_str(), eLowN, eLowMin, eLowMax);
      H_Channel_Cutted[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_Cutted[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_Cutted[i]->GetXaxis()->CenterTitle();
      H_Channel_Cutted[i]->GetYaxis()->CenterTitle();
    }
  }

  //////////////////////////////
  /////// FROM RAW GROUP ///////
  //////////////////////////////

  hist_label = new TH2D("hist_label", "hist_label", 80, 10, 90, 80, 10, 90);
  hist_label->GetXaxis()->SetTitle("Strip (grouped )");
  hist_label->GetYaxis()->SetTitle("Strip");
  hist_label->GetXaxis()->CenterTitle();
  hist_label->GetYaxis()->CenterTitle();

  /// HOW MANY STRIP IN THE RAW GROUP
  H_Strip_Mulitplicity_RAW = new TH1D(("Strip_Mulitplicity_RAW"), ("Strip_Mulitplicity_RAW"), 10, 0, 10);
  H_Strip_Mulitplicity_RAW->GetXaxis()->SetTitle("Multiplicity");
  H_Strip_Mulitplicity_RAW->GetYaxis()->SetTitle("Counts");
  H_Strip_Mulitplicity_RAW->GetXaxis()->CenterTitle();
  H_Strip_Mulitplicity_RAW->GetYaxis()->CenterTitle();

  /// HOW MANY REAR IN THE RAW GROUP
  H_Rear_Mulitplicity_RAW = new TH1D(("Rear_Mulitplicity_RAW"), ("Rear_Mulitplicity_RAW"), 10, 0, 10);
  H_Rear_Mulitplicity_RAW->GetXaxis()->SetTitle("Multiplicity");
  H_Rear_Mulitplicity_RAW->GetYaxis()->SetTitle("Counts");
  H_Rear_Mulitplicity_RAW->GetXaxis()->CenterTitle();
  H_Rear_Mulitplicity_RAW->GetYaxis()->CenterTitle();

  /// HOW MANY SiPM High IN THE RAW GROUP
  H_SiPMHigh_Mulitplicity_RAW = new TH1D(("SiPMHigh_Mulitplicity_RAW"), ("SiPMHigh_Mulitplicity_RAW"), 20, 0, 20);
  H_SiPMHigh_Mulitplicity_RAW->GetXaxis()->SetTitle("Multiplicity");
  H_SiPMHigh_Mulitplicity_RAW->GetYaxis()->SetTitle("Counts");
  H_SiPMHigh_Mulitplicity_RAW->GetXaxis()->CenterTitle();
  H_SiPMHigh_Mulitplicity_RAW->GetYaxis()->CenterTitle();

  /// HOW MANY SiPM Low IN THE RAW GROUP
  H_SiPMLow_Mulitplicity_RAW = new TH1D(("SiPMLow_Mulitplicity_RAW"), ("SiPMLow_Mulitplicity_RAW"), 20, 0, 20);
  H_SiPMLow_Mulitplicity_RAW->GetXaxis()->SetTitle("Multiplicity");
  H_SiPMLow_Mulitplicity_RAW->GetYaxis()->SetTitle("Counts");
  H_SiPMLow_Mulitplicity_RAW->GetXaxis()->CenterTitle();
  H_SiPMLow_Mulitplicity_RAW->GetYaxis()->CenterTitle();

  //////////////////////////////////////////
  /////// FROM ASSOCIATED REAR-STRIP ///////
  //////////////////////////////////////////

  /// HOW MANY STRIP IN THE ASSOCIATED GROUP
  H_Strip_Mulitplicity_ASSOCIATED = new TH1D(("Strip_Mulitplicity_ASSOCIATED"), ("Strip_Mulitplicity_ASSOCIATED"), 10, 0, 10);
  H_Strip_Mulitplicity_ASSOCIATED->GetXaxis()->SetTitle("Multiplicity");
  H_Strip_Mulitplicity_ASSOCIATED->GetYaxis()->SetTitle("Counts");
  H_Strip_Mulitplicity_ASSOCIATED->GetXaxis()->CenterTitle();
  H_Strip_Mulitplicity_ASSOCIATED->GetYaxis()->CenterTitle();

  /// HOW MANY REAR IN THE ASSOCIATED GROUP
  H_Rear_Mulitplicity_ASSOCIATED = new TH1D(("Rear_Mulitplicity_ASSOCIATED"), ("Rear_Mulitplicity_ASSOCIATED"), 10, 0, 10);
  H_Rear_Mulitplicity_ASSOCIATED->GetXaxis()->SetTitle("Multiplicity");
  H_Rear_Mulitplicity_ASSOCIATED->GetYaxis()->SetTitle("Counts");
  H_Rear_Mulitplicity_ASSOCIATED->GetXaxis()->CenterTitle();
  H_Rear_Mulitplicity_ASSOCIATED->GetYaxis()->CenterTitle();

  /// WHICH ASSOCIATED
  H_RearStrip_ASSOCIATED = new TH2D("RearStrip_ASSOCIATED", "RearStrip_ASSOCIATED", 80, 10, 90, 8, 1, 9);
  H_RearStrip_ASSOCIATED->GetXaxis()->SetTitle("Strip");
  H_RearStrip_ASSOCIATED->GetYaxis()->SetTitle("Rear");
  H_RearStrip_ASSOCIATED->GetXaxis()->CenterTitle();
  H_RearStrip_ASSOCIATED->GetYaxis()->CenterTitle();

  /// WHICH NON-ASSOCIATED
  H_RearStrip_NON_ASSOCIATED = new TH2D("RearStrip_NON_ASSOCIATED", "RearStrip_NON_ASSOCIATED", 80, 10, 90, 8, 1, 9);
  H_RearStrip_NON_ASSOCIATED->GetXaxis()->SetTitle("Strip");
  H_RearStrip_NON_ASSOCIATED->GetYaxis()->SetTitle("Rear");
  H_RearStrip_NON_ASSOCIATED->GetXaxis()->CenterTitle();
  H_RearStrip_NON_ASSOCIATED->GetYaxis()->CenterTitle();

  //////////////////////////////////////////////////////
  /////// FOR CASE E : MULTIPLE DETECTOR DETECTION /////
  //////////////////////////////////////////////////////

  // TIME BETWEEN REAR AND STRIP
  H_RearStrip_Time_E = new TH1D("RearStrip_Time_E", "RearStrip_Time_E", winGroupN, winGroupMin, winGroupMax);
  H_RearStrip_Time_E->GetXaxis()->SetTitle("Rear-Strip Time (ns)");
  H_RearStrip_Time_E->GetYaxis()->SetTitle("Counts");
  H_RearStrip_Time_E->GetXaxis()->CenterTitle();
  H_RearStrip_Time_E->GetYaxis()->CenterTitle();

  /// ENERGY SHARING BETWEEN THE 2 STRIPS
  H_2Strip_Channel_E = new TH2D("2Strip_Channel_E", "2Strip_Channel_E", eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
  H_2Strip_Channel_E->GetXaxis()->SetTitle("Channel");
  H_2Strip_Channel_E->GetYaxis()->SetTitle("Channel");
  H_2Strip_Channel_E->GetXaxis()->CenterTitle();
  H_2Strip_Channel_E->GetYaxis()->CenterTitle();

  /// ENERGY BETWEEN STRIP AND REAR
  H_RearStrip_Channel_E = new TH2D("RearStrip_Channel_E", "RearStrip_Channel_E", eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
  H_RearStrip_Channel_E->GetXaxis()->SetTitle("Rear Channel");
  H_RearStrip_Channel_E->GetYaxis()->SetTitle("Strip Channel");
  H_RearStrip_Channel_E->GetXaxis()->CenterTitle();
  H_RearStrip_Channel_E->GetYaxis()->CenterTitle();

  /// Label BETWEEN THE 2 STRIPS
  H_2Strip_Label_E = new TH2D("2Strip_Label_E", "2Strip_Label_E", 80, 10, 90, 80, 10, 90);
  H_2Strip_Label_E->GetXaxis()->SetTitle("Label");
  H_2Strip_Label_E->GetYaxis()->SetTitle("Label");
  H_2Strip_Label_E->GetXaxis()->CenterTitle();
  H_2Strip_Label_E->GetYaxis()->CenterTitle();

  /// TIME BETWEEN THE 2 STRIPS
  H_2Strip_Time_E = new TH1D("2Strip_Time_E", "2Strip_Time_E", winGroupN, winGroupMin, winGroupMax);
  H_2Strip_Time_E->GetXaxis()->SetTitle("Time between 2 signals (ns)");
  H_2Strip_Time_E->GetYaxis()->SetTitle("Counts");
  H_2Strip_Time_E->GetXaxis()->CenterTitle();
  H_2Strip_Time_E->GetYaxis()->CenterTitle();

  //////////////////////////////////////////
  /////// FROM CUTTED //////////////////////
  //////////////////////////////////////////

  H_SiPMLow_Mulitplicity_Cutted = new TH1D(("SiPMLow_Mulitplicity_Cutted"), ("SiPMLow_Mulitplicity_Cutted"), 20, 0, 20);
  H_SiPMLow_Mulitplicity_Cutted->GetXaxis()->SetTitle("Multiplicity");
  H_SiPMLow_Mulitplicity_Cutted->GetYaxis()->SetTitle("Counts");
  H_SiPMLow_Mulitplicity_Cutted->GetXaxis()->CenterTitle();
  H_SiPMLow_Mulitplicity_Cutted->GetYaxis()->CenterTitle();

  H_SiPMHigh_Mulitplicity_Cutted = new TH1D(("SiPMHigh_Mulitplicity_Cutted"), ("SiPMHigh_Mulitplicity_Cutted"), 20, 0, 20);
  H_SiPMHigh_Mulitplicity_Cutted->GetXaxis()->SetTitle("Multiplicity");
  H_SiPMHigh_Mulitplicity_Cutted->GetYaxis()->SetTitle("Counts");
  H_SiPMHigh_Mulitplicity_Cutted->GetXaxis()->CenterTitle();
  H_SiPMHigh_Mulitplicity_Cutted->GetYaxis()->CenterTitle();

  return 0;
}

inline int InitTree_Grouped()
{
  GROUPED_Tree = new TTree("Tree", "Tree");
  GROUPED_Tree->Branch("Tree_Silicon", &Tree_Silicon);
  GROUPED_Tree->Branch("Tree_SiPMHigh", &Tree_SiPMHigh);
  GROUPED_Tree->Branch("Tree_SiPMLow", &Tree_SiPMLow);
  return 0;
}

inline int WriteHistograms_Grouped()
{
  GROUPED_File->cd();
  hist_label->Write();

  for (size_t i = 0; i < detectorNum; ++i)
  {
    GROUPED_File->cd();
    if (IsDetectorSiliStrip(i))
    {
      hist_channel[i]->Write();
      hist_time[i]->Write();
      

      // CHANNEL
      dir_Strip_Channel_Write->cd();
      H_Channel_RAW[i]->Write();
      H_Channel_ASSOCIATED[i]->Write();
      H_Channel_A[i]->Write();
      H_Channel_B[i]->Write();
      H_Channel_C[i]->Write();
      H_Channel_D[i]->Write();
      H_Channel_E[i]->Write();

      dir_Strip_Channel->cd();
      TCanvas *C_RAW_ASSOCIATED = new TCanvas(("C_RAW_ASSOCIATED_" + detectorName[i]).c_str(), ("C_RAW_ASSOCIATED_" + detectorName[i]).c_str(), 800, 400);
      C_RAW_ASSOCIATED->Divide(1, 2);
      C_RAW_ASSOCIATED->cd(1);
      H_Channel_RAW[i]->Draw();
      H_Channel_ASSOCIATED[i]->SetLineColor(kRed);
      H_Channel_ASSOCIATED[i]->Draw("SAME");
      C_RAW_ASSOCIATED->cd(2);
      TH1D *H_Channel_RAW_ASSOCIATED_SUB = (TH1D *)H_Channel_RAW[i]->Clone(("H_Channel_RAW_ASSOCIATED_SUB" + detectorName[i]).c_str());
      H_Channel_RAW_ASSOCIATED_SUB->Add(H_Channel_ASSOCIATED[i], -1);
      H_Channel_RAW_ASSOCIATED_SUB->Draw();
      C_RAW_ASSOCIATED->Write();

      TCanvas *C_ASSOCIATED_CASE_A = new TCanvas(("C_ASSOCIATED_CASE_A_" + detectorName[i]).c_str(), ("C_ASSOCIATED_CASE_A_" + detectorName[i]).c_str(), 800, 400);
      C_ASSOCIATED_CASE_A->Divide(1, 2);
      C_ASSOCIATED_CASE_A->cd(1);
      H_Channel_ASSOCIATED[i]->SetLineColor(kBlue);
      H_Channel_ASSOCIATED[i]->Draw();
      H_Channel_A[i]->SetLineColor(kRed);
      H_Channel_A[i]->Draw("SAME");
      C_ASSOCIATED_CASE_A->cd(2);
      TH1D *H_Channel_ASSOCIATED_A_SUB = (TH1D *)H_Channel_ASSOCIATED[i]->Clone(("H_Channel_ASSOCIATED_A_SUB" + detectorName[i]).c_str());
      H_Channel_ASSOCIATED_A_SUB->Add(H_Channel_A[i], -1);
      H_Channel_ASSOCIATED_A_SUB->Draw();
      C_ASSOCIATED_CASE_A->Write();

      /// frac strip rear

      dir_Strip_Channel_DiffRear->cd();

      double x_center = H_Strip_Channel_DiffRear_A[i]->GetBinCenter(H_Strip_Channel_DiffRear_A[i]->GetMaximumBin());
      int counter = 0;

      TGraphErrors *graph_spread = new TGraphErrors();
      graph_spread->SetTitle("Spread");
      graph_spread->GetXaxis()->SetTitle("Strip Channel");
      graph_spread->GetYaxis()->SetTitle("Spread (#sigma of #DeltaE/E)");
      graph_spread->GetXaxis()->CenterTitle();
      graph_spread->GetYaxis()->CenterTitle();

      // graph mean
      TGraphErrors *graph_mean = new TGraphErrors();
      graph_mean->SetTitle("Mean");
      graph_mean->GetXaxis()->SetTitle("Channel");
      graph_mean->GetYaxis()->SetTitle("Mean of #DeltaE/E");
      graph_mean->GetXaxis()->CenterTitle();
      graph_mean->GetYaxis()->CenterTitle();

      TH1D *H_Strip_Channel_DiffRear_A_YProjection = H_Strip_Channel_DiffRearvsStrip_A[i]->ProjectionY(("H_Strip_Channel_DiffRearvsStrip_A_YProjection_" + detectorName[i]).c_str());
      for (int bin = 0; bin < H_Strip_Channel_DiffRear_A[i]->GetNbinsX(); bin += 10)
      {
        TH1D *H_Strip_Channel_DiffRear_A_XProjection = H_Strip_Channel_DiffRearvsStrip_A[i]->ProjectionX(("H_Strip_Channel_DiffRearvsStrip_A_XProjection_" + detectorName[i] + "_" + to_string(bin)).c_str(), bin, bin + 10);
        H_Strip_Channel_DiffRear_A_XProjection->GetXaxis()->SetRangeUser(x_center - 0.05, x_center + 0.05);
        TF1 *gauss = new TF1(("gauss_" + detectorName[i] + "_" + to_string(bin)).c_str(), "gaus", x_center - 0.05, x_center + 0.05);
        gauss->SetParameter(1, x_center);
        gauss->SetParameter(2, 0.005);
        if (H_Strip_Channel_DiffRear_A_XProjection->Integral() < 20)
          continue; /// avoid empty data for fit

        H_Strip_Channel_DiffRear_A_XProjection->Fit(gauss, "QRN");

        if (gauss->GetParError(2) / gauss->GetParameter(2) > 0.2)
          continue; /// avoid weird data points due to low statistics

        graph_spread->SetPoint(counter, H_Strip_Channel_DiffRear_A_YProjection->GetBinCenter(bin + 5), gauss->GetParameter(2));
        graph_spread->SetPointError(counter, H_Strip_Channel_DiffRear_A_YProjection->GetBinCenter(bin + 5) - H_Strip_Channel_DiffRear_A_YProjection->GetBinCenter(bin), gauss->GetParError(2));

        graph_mean->SetPoint(counter, H_Strip_Channel_DiffRear_A_YProjection->GetBinCenter(bin + 5), gauss->GetParameter(1));
        graph_mean->SetPointError(counter, H_Strip_Channel_DiffRear_A_YProjection->GetBinCenter(bin + 5) - H_Strip_Channel_DiffRear_A_YProjection->GetBinCenter(bin), gauss->GetParError(1));
        counter++;
      }

      TF1 *linear = new TF1("linear", "[0]");
      graph_mean->Fit(linear, "Q");
      double x_mean = linear->GetParameter(0);
      MeanAcceptance[i] = x_mean;

      C_Strip_Channel_DiffRear[i] = new TCanvas(("C_Strip_Channel_DiffRear_" + detectorName[i]).c_str(), ("C_Strip_Channel_DiffRear_" + detectorName[i]).c_str(), 800, 400);

      ////1st
      C_Strip_Channel_DiffRear[i]->cd();
      TPad *P_RearStrip = new TPad("RearStrip", "RearStrip", 0.05, 0.5, 0.5, 1.0);
      P_RearStrip->Draw();
      P_RearStrip->cd();
      H_RearStrip_Channel_A[i]->Draw("COLZ");

      /////3rd
      C_Strip_Channel_DiffRear[i]->cd();
      TPad *P_SpreadRearStripvsE = new TPad("SpreadRearStripvsE", "SpreadRearStripvsE", 0.05, 0.05, 0.5, 0.5);
      P_SpreadRearStripvsE->Draw();
      P_SpreadRearStripvsE->cd();
      graph_spread->SetMarkerStyle(20);
      graph_spread->SetMarkerSize(0.5);
      graph_spread->Draw("AP");
      SpreadAcceptance[i] = new TF1("spread_fit", "[0]/x");
      SpreadAcceptance[i]->SetParameter(0, 100);
      graph_spread->Fit(SpreadAcceptance[i], "Q");

      double yfullrange = graph_spread->GetYaxis()->GetXmax() - graph_spread->GetYaxis()->GetXmin();
      TLatex *formula = new TLatex(60000, yfullrange * 0.85, "#bf{#color[2]{f(x) = #frac{a}{x}}}");
      formula->Draw();
      TLatex *a = new TLatex(60000, yfullrange * 0.7, ("#bf{#color[2]{a = " + formatValueWithError(SpreadAcceptance[i]->GetParameter(0), SpreadAcceptance[i]->GetParError(0)) + "}}").c_str());
      a->Draw();

      ////2nd
      C_Strip_Channel_DiffRear[i]->cd();
      TPad *P_DiffRearStripvsStrip = new TPad("DiffRearStripvsStrip", "DiffRearStripvsStrip", 0.55, 0.50, 1.0, 1.0);
      P_DiffRearStripvsStrip->Draw();
      P_DiffRearStripvsStrip->cd();
      H_Strip_Channel_DiffRearvsStrip_A[i]->Draw("COLZ");
      spread_up[i] = new TGraphErrors();
      spread_down[i] = new TGraphErrors();
      for (int bin = 0; bin < H_Strip_Channel_DiffRearvsStrip_A[i]->GetNbinsY(); bin++)
      {
        double spread = SpreadAcceptance[i]->Eval(H_Strip_Channel_DiffRearvsStrip_A[i]->GetYaxis()->GetBinCenter(bin));
        spread_up[i]->SetPoint(bin, x_mean + 3 * spread, H_Strip_Channel_DiffRearvsStrip_A[i]->GetYaxis()->GetBinCenter(bin));
        spread_down[i]->SetPoint(bin, x_mean - 3 * spread, H_Strip_Channel_DiffRearvsStrip_A[i]->GetYaxis()->GetBinCenter(bin));
      }

      ////4th
      C_Strip_Channel_DiffRear[i]->cd();
      TPad *P_DiffRearStrip = new TPad("DiffRearStrip", "DiffRearStrip", 0.55, 0.05, 1.0, 0.5);
      P_DiffRearStrip->Draw();
      P_DiffRearStrip->cd();
      P_DiffRearStrip->SetLogy();
      H_Strip_Channel_DiffRear_A[i]->Draw();
      C_Strip_Channel_DiffRear[i]->cd();

      ///////deleting
      delete C_RAW_ASSOCIATED;
      delete H_Channel_RAW_ASSOCIATED_SUB;
      delete H_Channel_ASSOCIATED[i];
      delete H_Channel_RAW[i];
      delete H_Channel_A[i];
      delete H_Channel_B[i];
      delete H_Channel_C[i];
      delete H_Channel_D[i];
      delete H_Channel_E[i];
      delete C_ASSOCIATED_CASE_A;
      delete H_Channel_ASSOCIATED_A_SUB;
    }

    if (IsDetectorSiliBack(i))
    {

      // CHANNEL
      dir_Rear_Channel_Write->cd();
      H_Channel_RAW[i]->Write();
      H_Channel_ASSOCIATED[i]->Write();
      H_Channel_A[i]->Write();
      H_Channel_B[i]->Write();
      H_Channel_C[i]->Write();
      H_Channel_D[i]->Write();
      H_Channel_E[i]->Write();

      dir_Rear_Channel->cd();
      TCanvas *C_RAW_ASSOCIATED = new TCanvas(("C_RAW_ASSOCIATED_" + detectorName[i]).c_str(), ("C_RAW_ASSOCIATED_" + detectorName[i]).c_str(), 800, 400);
      C_RAW_ASSOCIATED->Divide(1, 2);
      C_RAW_ASSOCIATED->cd(1);
      H_Channel_RAW[i]->Draw();
      H_Channel_ASSOCIATED[i]->SetLineColor(kRed);
      H_Channel_ASSOCIATED[i]->Draw("SAME");
      C_RAW_ASSOCIATED->cd(2);
      TH1D *H_Channel_RAW_ASSOCIATED_SUB = (TH1D *)H_Channel_RAW[i]->Clone(("H_Channel_RAW_ASSOCIATED_SUB" + detectorName[i]).c_str());
      H_Channel_RAW_ASSOCIATED_SUB->SetTitle(("H_Channel_RAW_ASSOCIATED_SUB" + detectorName[i]).c_str());
      H_Channel_RAW_ASSOCIATED_SUB->Add(H_Channel_ASSOCIATED[i], -1);
      H_Channel_RAW_ASSOCIATED_SUB->Draw();
      C_RAW_ASSOCIATED->Write();

      TCanvas *C_ASSOCIATED_CASE_A = new TCanvas(("C_ASSOCIATED_CASE_A_" + detectorName[i]).c_str(), ("C_ASSOCIATED_CASE_A_" + detectorName[i]).c_str(), 800, 400);
      C_ASSOCIATED_CASE_A->Divide(1, 2);
      C_ASSOCIATED_CASE_A->cd(1);
      H_Channel_ASSOCIATED[i]->SetLineColor(kBlue);
      H_Channel_ASSOCIATED[i]->Draw();
      H_Channel_A[i]->SetLineColor(kRed);
      H_Channel_A[i]->Draw("SAME");
      C_ASSOCIATED_CASE_A->cd(2);
      TH1D *H_Channel_ASSOCIATED_A_SUB = (TH1D *)H_Channel_ASSOCIATED[i]->Clone(("H_Channel_ASSOCIATED_A_SUB" + detectorName[i]).c_str());
      H_Channel_ASSOCIATED_A_SUB->SetTitle(("H_Channel_ASSOCIATED_A_SUB" + detectorName[i]).c_str());
      H_Channel_ASSOCIATED_A_SUB->Add(H_Channel_A[i], -1);
      H_Channel_ASSOCIATED_A_SUB->Draw();
      C_ASSOCIATED_CASE_A->Write();

      // TIME

      /////// deleting
      delete C_RAW_ASSOCIATED;
      delete H_Channel_RAW_ASSOCIATED_SUB;
      delete H_Channel_RAW[i];
      delete H_Channel_ASSOCIATED[i];

      //// REAR - STRIP

      // CHANNEL
      dir_Strip_Rear->cd();
      dir_Rear_Strip_Channel->cd();
      H_RearStrip_Channel_ASSOCIATED[i]->Write();
      H_RearStrip_Channel_A[i]->Write();
      H_RearStrip_Channel_C[i]->Write();
      H_Rear2Strip_Channel_C[i]->Write();
      H_RearStrip_Channel_D[i]->Write();
      H_Rear2Strip_Channel_D[i]->Write();

      // TIME
      dir_Strip_Time_Write->cd();
      H_RearStrip_Time_ASSOCIATED[i]->Write();
      H_RearStrip_Time_A[i]->Write();
      H_RearStrip_Time_B[i]->Write();
      H_RearStrip_Time_C[i]->Write();
      H_RearStrip_Time_D[i]->Write();

      dir_Strip_Time->cd();
      TCanvas *C_TIME_ASSOCIATED_CASES = new TCanvas(("C_TIME_ASSOCIATED_CASES_" + detectorName[i]).c_str(), ("C_TIME_ASSOCIATED_CASES_" + detectorName[i]).c_str(), 800, 400);
      C_TIME_ASSOCIATED_CASES->cd();
      H_RearStrip_Time_ASSOCIATED[i]->SetTitle(("Time_ASSOCIATED_CASES " + detectorName[i]).c_str());
      H_RearStrip_Time_ASSOCIATED[i]->SetLineColor(kBlack);
      H_RearStrip_Time_ASSOCIATED[i]->Draw();
      H_RearStrip_Time_A[i]->SetLineColor(kRed);
      H_RearStrip_Time_A[i]->Draw("SAME");
      H_RearStrip_Time_B[i]->SetLineColor(kGreen);
      H_RearStrip_Time_B[i]->Draw("SAME");
      H_RearStrip_Time_C[i]->SetLineColor(kBlue);
      H_RearStrip_Time_C[i]->Draw("SAME");
      H_RearStrip_Time_D[i]->SetLineColor(kViolet);
      H_RearStrip_Time_D[i]->Draw("SAME");
      TLegend *legend = new TLegend(0.1, 0.7, 0.25, 0.9);
      legend->AddEntry(H_RearStrip_Time_ASSOCIATED[i], "ASSOCIATED", "l");
      legend->AddEntry(H_RearStrip_Time_A[i], "Case A", "l");
      legend->AddEntry(H_RearStrip_Time_B[i], "Case B", "l");
      legend->AddEntry(H_RearStrip_Time_C[i], "Case C", "l");
      legend->AddEntry(H_RearStrip_Time_D[i], "Case D", "l");
      legend->Draw();
      // C_TIME_ASSOCIATED_CASES_CUTTED->Write();
    }

    if (IsDetectorBetaHigh(i))
    {
      dir_SiPM_High->cd();

      dir_SiPMHigh_Channel_Write->cd();
      H_Channel_RAW[i]->Write();
    }

    if (IsDetectorBetaLow(i))
    {
      dir_SiPM_High->cd();

      dir_SiPMLow_Channel_Write->cd();
      H_Channel_RAW[i]->Write();
    }
  }

  /// STRIP/////////////////////////////////////////////
  dir_Strip_Multiplicities->cd();
  H_Strip_Mulitplicity_RAW->Write();
  H_Strip_Mulitplicity_ASSOCIATED->Write();

  TCanvas *C_Strip_Multiplicities = new TCanvas("Strip_Multiplicities_RAW_ASSOCIATED", "Strip_Multiplicities_RAW_ASSOCIATED", 800, 400);
  H_Strip_Mulitplicity_RAW->Draw();
  H_Strip_Mulitplicity_ASSOCIATED->SetLineColor(kRed);
  H_Strip_Mulitplicity_ASSOCIATED->Draw("SAME");
  C_Strip_Multiplicities->Write();

  dir_Strip_Time->cd();
  H_RearStrip_Time_E->Write();

  // deleting
  delete C_Strip_Multiplicities;
  delete H_Strip_Mulitplicity_RAW;
  delete H_Strip_Mulitplicity_ASSOCIATED;

  /// REAR/////////////////////////////////////////////
  dir_Rear_Multiplicity->cd();
  H_Rear_Mulitplicity_RAW->Write();
  H_Rear_Mulitplicity_ASSOCIATED->Write();

  TCanvas *C_Rear_Multiplicities = new TCanvas("Rear_Multiplicities_RAW_ASSOCIATED", "Rear_Multiplicities_RAW_ASSOCIATED", 800, 400);
  H_Rear_Mulitplicity_RAW->Draw();
  H_Rear_Mulitplicity_ASSOCIATED->SetLineColor(kRed);
  H_Rear_Mulitplicity_ASSOCIATED->Draw("SAME");
  C_Rear_Multiplicities->Write();

  // deleting
  delete C_Rear_Multiplicities;
  delete H_Rear_Mulitplicity_RAW;
  delete H_Rear_Mulitplicity_ASSOCIATED;

  /// STRIP-REAR/////////////////////////////////////////////
  dir_Strip_Rear->cd();
  H_RearStrip_ASSOCIATED->Write();
  H_RearStrip_NON_ASSOCIATED->Write();

  dir_Rear_Strip_Channel->cd();
  H_2Strip_Channel_E->Write();
  H_RearStrip_Channel_E->Write();
  dir_Strip_Rear->cd();
  H_2Strip_Label_E->Write();
  dir_Strip_Rear_Time->cd();
  H_2Strip_Time_E->Write();

  /// SiPM Low/////////////////////////////////////////////
  dir_SiPMLow_Multiplicity->cd();
  H_SiPMLow_Mulitplicity_RAW->Write();

  /// SiPM High/////////////////////////////////////////////
  dir_SiPM_Multiplicities->cd();
  H_SiPMHigh_Mulitplicity_RAW->Write();

  return 0;
}

void ProcessCuttingGroups()
{
}

inline void WriteHistograms_Cutted()
{
  GROUPED_File->cd();

  for (int i = 0; i < detectorNum; i++)
  {
    if (IsDetectorSiliStrip(i))
    {
      TPad *P_DiffRearStrip = (TPad *)C_Strip_Channel_DiffRear[i]->GetPrimitive("DiffRearStrip");
      P_DiffRearStrip->cd();
      H_Strip_Channel_DiffRear_Cutted[i]->SetLineColor(kRed);
      H_Strip_Channel_DiffRear_Cutted[i]->Draw("SAME");

      TPad *P_DiffRearStripvsStrip = (TPad *)C_Strip_Channel_DiffRear[i]->GetPrimitive("DiffRearStripvsStrip");
      P_DiffRearStripvsStrip->cd();
      spread_up[i]->SetLineColor(kRed);
      spread_down[i]->SetLineColor(kRed);
      spread_up[i]->Draw("SAME");
      spread_down[i]->Draw("SAME");

      dir_Strip_Channel_DiffRear->cd();
      C_Strip_Channel_DiffRear[i]->Write();
    }

    if (IsDetectorBetaHigh(i))
    {
      /// SIPM High
      dir_SiPMHigh_Channel_Write->cd();
      H_Channel_Cutted[i]->Write();
      dir_SiPMHigh_Channel->cd();
      TCanvas *C_SiPMHigh_Channel = new TCanvas(("C_SiPMHigh_Channel_RAW_CUTTED" + detectorName[i]).c_str(), ("C_SiPMHigh_Channel_" + detectorName[i]).c_str(), 800, 400);
      H_Channel_RAW[i]->Draw();
      H_Channel_Cutted[i]->SetLineColor(kRed);
      H_Channel_Cutted[i]->Draw("SAME");
      C_SiPMHigh_Channel->Write();
    }

    if (IsDetectorBetaLow(i))
    {
      /// SIPM Low
      dir_SiPMLow_Channel_Write->cd();
      H_Channel_Cutted[i]->Write();
      dir_SiPMLow_Channel->cd();
      TCanvas *C_SiPMLow_Channel = new TCanvas(("C_SiPMLow_Channel_RAW_CUTTED" + detectorName[i]).c_str(), ("C_SiPMLow_Channel_" + detectorName[i]).c_str(), 800, 400);
      H_Channel_RAW[i]->Draw();
      H_Channel_Cutted[i]->SetLineColor(kRed);
      H_Channel_Cutted[i]->Draw("SAME");
      C_SiPMLow_Channel->Write();
    }
  }

  dir_SiPMHigh_Multiplicity->cd();
    H_SiPMHigh_Mulitplicity_Cutted->Write();
    TCanvas *C_SiPMHigh_Multiplicities = new TCanvas("SiPMHigh_Multiplicities_RAW_CUTTED", "SiPMHigh_Multiplicities_RAW_CUTTED", 800, 400);
    H_SiPMHigh_Mulitplicity_RAW->Draw();
    H_SiPMHigh_Mulitplicity_Cutted->SetLineColor(kRed);
    H_SiPMHigh_Mulitplicity_Cutted->Draw("SAME");
    C_SiPMHigh_Multiplicities->Write();

    dir_SiPMLow_Multiplicity->cd();
    H_SiPMLow_Mulitplicity_Cutted->Write();
    TCanvas *C_SiPMLow_Multiplicities = new TCanvas("SiPMLow_Multiplicities_RAW_CUTTED", "SiPMLow_Multiplicities_RAW_CUTTED", 800, 400);
    H_SiPMLow_Mulitplicity_RAW->Draw();
    H_SiPMLow_Mulitplicity_Cutted->SetLineColor(kRed);
    H_SiPMLow_Mulitplicity_Cutted->Draw("SAME");
    C_SiPMLow_Multiplicities->Write();
}

inline int WriteTree_Grouped()
{
  GROUPED_File->cd();
  GROUPED_Tree->Write();
  delete GROUPED_Tree;
  return 0;
}

void SearchForCoincidence(TTreeReaderArray<Signal> &signals)
{
  if (Verbose > 0)
    cout << "### Starting a Coincidence group on " << detectorName[signals[0].Label] << "  TIME : " << signals[0].Time << endl;

  vector<int> Rear_Position;
  vector<int> Strip_Position;
  vector<int> SiPM_High_Position;
  vector<int> SiPM_Low_Position;

  int Rear_Multiplicity_RAW = 0;
  int Strip_Multiplicity_RAW = 0;
  int SiPM_High_Multiplicity_RAW = 0;
  int SiPM_Low_Multiplicity_RAW = 0;

  for (int index = 0; index < signals.GetSize(); index++)
  {
    int current_label = signals[index].Label;
    if (IsDetectorSiliBack(current_label))
    {
      Rear_Position.push_back(index);
      Rear_Multiplicity_RAW++;
      H_Channel_RAW[current_label]->Fill(signals[index].Channel);
    }
    else if (IsDetectorSiliStrip(current_label))
    {
      Strip_Position.push_back(index);
      Strip_Multiplicity_RAW++;
      H_Channel_RAW[current_label]->Fill(signals[index].Channel);
    }
    else if (IsDetectorBetaHigh(current_label))
    {
      SiPM_High_Position.push_back(index);
      SiPM_High_Multiplicity_RAW++;
      H_Channel_RAW[current_label]->Fill(signals[index].Channel);
    }
    else if (IsDetectorBetaLow(current_label))
    {
      SiPM_Low_Position.push_back(index);
      SiPM_Low_Multiplicity_RAW++;
      H_Channel_RAW[current_label]->Fill(signals[index].Channel);
    }
  }

  ////// RAW //////
  H_Strip_Mulitplicity_RAW->Fill(Strip_Multiplicity_RAW);
  H_Rear_Mulitplicity_RAW->Fill(Rear_Multiplicity_RAW);
  H_SiPMHigh_Mulitplicity_RAW->Fill(SiPM_High_Multiplicity_RAW);
  H_SiPMLow_Mulitplicity_RAW->Fill(SiPM_Low_Multiplicity_RAW);
  /////////////////

  ///////////////////////////////////
  ///// ASSOCIATED REAR - STRIP /////
  ///////////////////////////////////
  // Determining Rear strip couples
  vector<pair<Signal, Signal>> RearStrip_ASSOCIATED;
  vector<pair<Signal, Signal>> RearStrip_NONASSOCIATED;
  vector<bool> Rear_Associated = vector<bool>(Rear_Position.size(), false);
  vector<bool> Strip_Associated = vector<bool>(Strip_Position.size(), false);
  for (int index_rear = 0; index_rear < Rear_Position.size(); index_rear++)
  {
    for (int index_strip = 0; index_strip < Strip_Position.size(); index_strip++)
    {
      if (IsSameSiliDetector(signals[Strip_Position[index_strip]].Label, signals[Rear_Position[index_rear]].Label))
      {
        RearStrip_ASSOCIATED.push_back(make_pair(signals[Rear_Position[index_rear]], signals[Strip_Position[index_strip]]));
        Strip_Associated[index_strip] = true;
        Rear_Associated[index_rear] = true;
      }
    }
  }

  // Saving Rear and Strip associated or non-associated
  int Rear_Multiplicity_ASSOCIATED = 0;
  int Strip_Multiplicity_ASSOCIATED = 0;
  for (int index_rear = 0; index_rear < Rear_Position.size(); index_rear++)
  {
    if (Rear_Associated[index_rear])
    {
      Rear_Multiplicity_ASSOCIATED++;
      H_Channel_ASSOCIATED[signals[Rear_Position[index_rear]].Label]->Fill(signals[Rear_Position[index_rear]].Channel);
    }
    else
    {
      for (int index_strip = 0; index_strip < Strip_Position.size(); index_strip++)
      {
        if (!Strip_Associated[index_strip])
        {
          RearStrip_NONASSOCIATED.push_back(make_pair(signals[Rear_Position[index_rear]], signals[Strip_Position[index_strip]]));
        }
      }
    }
  }

  for (int index_strip = 0; index_strip < Strip_Position.size(); index_strip++)
  {
    if (Strip_Associated[index_strip])
    {
      Strip_Multiplicity_ASSOCIATED++;
      H_Channel_ASSOCIATED[signals[Strip_Position[index_strip]].Label]->Fill(signals[Strip_Position[index_strip]].Channel);
    }
  }

  H_Strip_Mulitplicity_ASSOCIATED->Fill(Strip_Multiplicity_ASSOCIATED);
  H_Rear_Mulitplicity_ASSOCIATED->Fill(Rear_Multiplicity_ASSOCIATED);
  //////////////////////////////

  // Saving Rear and Strip  associated and NON-associated
  for (auto &pair : RearStrip_ASSOCIATED)
  {
    H_RearStrip_Channel_ASSOCIATED[pair.first.Label]->Fill(pair.first.Channel, pair.second.Channel);
    H_RearStrip_Time_ASSOCIATED[pair.first.Label]->Fill(pair.first.Time - pair.second.Time);
    H_RearStrip_ASSOCIATED->Fill(pair.second.Label, GetDetector(pair.first.Label));
  }
  for (auto &pair : RearStrip_NONASSOCIATED) // trying to couple the non-associated rear-strip
  {
    H_RearStrip_NON_ASSOCIATED->Fill(pair.second.Label, GetDetector(pair.first.Label));
  }

  ///////////////////////////////////
  ///////// CASES A B C D E /////////
  ///////////////////////////////////

  // Case A : SINGLE Rear and Strip associated (most common case)
  if (RearStrip_ASSOCIATED.size() == 1 && Rear_Position.size() == 1 && Strip_Position.size() == 2)
  {
    H_Channel_A[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[0].first.Channel);   // REAR
    H_Channel_A[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].second.Channel); // STRIP

    H_RearStrip_Channel_A[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[0].first.Channel, RearStrip_ASSOCIATED[0].second.Channel);                                                                                                  // REAR-STRIP for rear plot
    H_RearStrip_Channel_A[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].first.Channel, RearStrip_ASSOCIATED[0].second.Channel);                                                                                                 // REAR-STRIP for strip plot
    H_Strip_Channel_DiffRear_A[RearStrip_ASSOCIATED[0].second.Label]->Fill((RearStrip_ASSOCIATED[0].first.Channel - RearStrip_ASSOCIATED[0].second.Channel) / RearStrip_ASSOCIATED[0].second.Channel);                                                // STRIP/REAR
    H_Strip_Channel_DiffRearvsStrip_A[RearStrip_ASSOCIATED[0].second.Label]->Fill((RearStrip_ASSOCIATED[0].first.Channel - RearStrip_ASSOCIATED[0].second.Channel) / RearStrip_ASSOCIATED[0].second.Channel, RearStrip_ASSOCIATED[0].second.Channel); // STRIP/REAR vs Strip channel
    H_RearStrip_Time_A[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[0].first.Time - RearStrip_ASSOCIATED[0].second.Time);

    for (int index_strip : Strip_Position)
    {
      if (signals[index_strip].Label != RearStrip_ASSOCIATED[0].second.Label && RearStrip_ASSOCIATED[0].second.Channel > 26000 && RearStrip_ASSOCIATED[0].second.Channel < 32000 )
      {
        hist_channel[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].second.Channel, signals[index_strip].Channel);
        hist_label->Fill(RearStrip_ASSOCIATED[0].second.Label, signals[index_strip].Label); 
        hist_time[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].second.Time - signals[index_strip].Time);   
        break;
      }
    }
    
  }

  // Case B : MULTIPLE SAME Rear and Strip associated
  //  do the code
  if (RearStrip_ASSOCIATED.size() > 1)
  {
    for (int index1 = 0; index1 < RearStrip_ASSOCIATED.size(); index1++)
    {
      for (int index2 = index1 + 1; index2 < RearStrip_ASSOCIATED.size(); index2++)
      {
        if (RearStrip_ASSOCIATED[index1].first.Label == RearStrip_ASSOCIATED[index2].first.Label && RearStrip_ASSOCIATED[index1].second.Label == RearStrip_ASSOCIATED[index2].second.Label)
        {
          H_Channel_B[RearStrip_ASSOCIATED[index1].first.Label]->Fill(RearStrip_ASSOCIATED[index1].first.Channel);   // REAR
          H_Channel_B[RearStrip_ASSOCIATED[index1].second.Label]->Fill(RearStrip_ASSOCIATED[index1].second.Channel); // STRIP

          H_Channel_B[RearStrip_ASSOCIATED[index2].first.Label]->Fill(RearStrip_ASSOCIATED[index2].first.Channel);   // REAR
          H_Channel_B[RearStrip_ASSOCIATED[index2].second.Label]->Fill(RearStrip_ASSOCIATED[index2].second.Channel); // STRIP

          // H_RearStrip_Channel_B[RearStrip_ASSOCIATED[index1].first.Label]->Fill(RearStrip_ASSOCIATED[index1].first.Channel, RearStrip_ASSOCIATED[index1].second.Channel);  //REAR-STRIP
          // H_RearStrip_Channel_B[RearStrip_ASSOCIATED[index2].first.Label]->Fill(RearStrip_ASSOCIATED[index2].first.Channel, RearStrip_ASSOCIATED[index2].second.Channel);  //REAR-STRIP
        }
      }
    }
  }

  // Case C : SAME Rear and DIFFERENT NEIGHBOURG Strip associated (intertrip)
  if (RearStrip_ASSOCIATED.size() > 1)
  {
    for (int index1 = 0; index1 < RearStrip_ASSOCIATED.size(); index1++)
    {
      for (int index2 = index1 + 1; index2 < RearStrip_ASSOCIATED.size(); index2++)
      {
        if (RearStrip_ASSOCIATED[index1].first.Label == RearStrip_ASSOCIATED[index2].first.Label)
        {
          if (IsDetectorSiliInterStrip(RearStrip_ASSOCIATED[index1].second.Label, RearStrip_ASSOCIATED[index2].second.Label))
          {
            H_Channel_C[RearStrip_ASSOCIATED[index1].first.Label]->Fill(RearStrip_ASSOCIATED[index1].first.Channel);   // REAR
            H_Channel_C[RearStrip_ASSOCIATED[index1].second.Label]->Fill(RearStrip_ASSOCIATED[index1].second.Channel); // STRIP
            H_Channel_C[RearStrip_ASSOCIATED[index2].second.Label]->Fill(RearStrip_ASSOCIATED[index2].second.Channel); // STRIP

            H_RearStrip_Channel_C[RearStrip_ASSOCIATED[index1].first.Label]->Fill(RearStrip_ASSOCIATED[index1].first.Channel, RearStrip_ASSOCIATED[index1].second.Channel);                                                // REAR-STRIP
            H_RearStrip_Channel_C[RearStrip_ASSOCIATED[index1].first.Label]->Fill(RearStrip_ASSOCIATED[index1].first.Channel, RearStrip_ASSOCIATED[index2].second.Channel);                                                // REAR-STRIP
            H_Rear2Strip_Channel_C[RearStrip_ASSOCIATED[index1].first.Label]->Fill(RearStrip_ASSOCIATED[index2].first.Channel, RearStrip_ASSOCIATED[index1].second.Channel + RearStrip_ASSOCIATED[index2].second.Channel); // REAR-STRIP+STRIP

            H_RearStrip_Time_C[RearStrip_ASSOCIATED[index1].first.Label]->Fill(RearStrip_ASSOCIATED[index1].first.Time - RearStrip_ASSOCIATED[index1].second.Time);
            H_RearStrip_Time_C[RearStrip_ASSOCIATED[index1].first.Label]->Fill(RearStrip_ASSOCIATED[index1].first.Time - RearStrip_ASSOCIATED[index2].second.Time);
          }
        }
      }
    }
  }

  // Case D : SAME Rear and DIFFERENT non-NEIGHBOURG Strip associated (non-intertrip)
  if (RearStrip_ASSOCIATED.size() > 1)
  {
    for (int index1 = 0; index1 < RearStrip_ASSOCIATED.size(); index1++)
    {
      for (int index2 = index1 + 1; index2 < RearStrip_ASSOCIATED.size(); index2++)
      {
        if (RearStrip_ASSOCIATED[index1].first.Label == RearStrip_ASSOCIATED[index2].first.Label)
        {
          if (!IsDetectorSiliInterStrip(RearStrip_ASSOCIATED[index1].second.Label, RearStrip_ASSOCIATED[index2].second.Label))
          {
            H_Channel_D[RearStrip_ASSOCIATED[index1].first.Label]->Fill(RearStrip_ASSOCIATED[index1].first.Channel);   // REAR
            H_Channel_D[RearStrip_ASSOCIATED[index1].second.Label]->Fill(RearStrip_ASSOCIATED[index1].second.Channel); // STRIP
            H_Channel_D[RearStrip_ASSOCIATED[index2].second.Label]->Fill(RearStrip_ASSOCIATED[index2].second.Channel); // STRIP

            H_RearStrip_Channel_D[RearStrip_ASSOCIATED[index1].first.Label]->Fill(RearStrip_ASSOCIATED[index1].first.Channel, RearStrip_ASSOCIATED[index1].second.Channel);                                                // REAR-STRIP
            H_RearStrip_Channel_D[RearStrip_ASSOCIATED[index1].first.Label]->Fill(RearStrip_ASSOCIATED[index1].first.Channel, RearStrip_ASSOCIATED[index2].second.Channel);                                                // REAR-STRIP
            H_Rear2Strip_Channel_D[RearStrip_ASSOCIATED[index1].first.Label]->Fill(RearStrip_ASSOCIATED[index2].first.Channel, RearStrip_ASSOCIATED[index1].second.Channel + RearStrip_ASSOCIATED[index2].second.Channel); // REAR-STRIP+STRIP

            H_RearStrip_Time_D[RearStrip_ASSOCIATED[index1].first.Label]->Fill(RearStrip_ASSOCIATED[index1].first.Time - RearStrip_ASSOCIATED[index1].second.Time);
            H_RearStrip_Time_D[RearStrip_ASSOCIATED[index1].first.Label]->Fill(RearStrip_ASSOCIATED[index1].first.Time - RearStrip_ASSOCIATED[index2].second.Time);
          }
        }
      }
    }
  }

  // Case E : DIFFERENT Rear and DIFFERENT Strip associated (multiple detector detection)
  if (RearStrip_ASSOCIATED.size() > 1)
  {
    for (int index1 = 0; index1 < RearStrip_ASSOCIATED.size(); index1++)
    {
      for (int index2 = index1 + 1; index2 < RearStrip_ASSOCIATED.size(); index2++)
      {
        if (GetDetector(RearStrip_ASSOCIATED[index1].first.Label) != GetDetector(RearStrip_ASSOCIATED[index2].first.Label))
        {
          H_Channel_E[RearStrip_ASSOCIATED[index1].first.Label]->Fill(RearStrip_ASSOCIATED[index1].first.Channel);   // REAR
          H_Channel_E[RearStrip_ASSOCIATED[index1].second.Label]->Fill(RearStrip_ASSOCIATED[index1].second.Channel); // STRIP
          H_Channel_E[RearStrip_ASSOCIATED[index2].first.Label]->Fill(RearStrip_ASSOCIATED[index2].first.Channel);   // REAR
          H_Channel_E[RearStrip_ASSOCIATED[index2].second.Label]->Fill(RearStrip_ASSOCIATED[index2].second.Channel); // STRIP

          H_RearStrip_Channel_E->Fill(RearStrip_ASSOCIATED[index1].first.Channel, RearStrip_ASSOCIATED[index1].second.Channel); // REAR-STRIP
          H_RearStrip_Channel_E->Fill(RearStrip_ASSOCIATED[index2].first.Channel, RearStrip_ASSOCIATED[index2].second.Channel); // REAR-STRIP
          H_2Strip_Channel_E->Fill(RearStrip_ASSOCIATED[index1].second.Channel, RearStrip_ASSOCIATED[index2].second.Channel);   // STRIP-STRIP
          H_2Strip_Label_E->Fill(RearStrip_ASSOCIATED[index1].second.Label, RearStrip_ASSOCIATED[index2].second.Label);         // STRIP-STRIP

          H_RearStrip_Time_E->Fill(RearStrip_ASSOCIATED[index1].first.Time - RearStrip_ASSOCIATED[index1].second.Time);
          H_RearStrip_Time_E->Fill(RearStrip_ASSOCIATED[index2].first.Time - RearStrip_ASSOCIATED[index2].second.Time);

          H_2Strip_Time_E->Fill(RearStrip_ASSOCIATED[index1].second.Time - RearStrip_ASSOCIATED[index2].second.Time);
        }
      }
    }
  }
}

void CuttingGroups(TTreeReaderArray<Signal> &signals)
{
  vector<int> Rear_Position;
  vector<int> Strip_Position;
  vector<int> SiPM_High_Position;
  vector<int> SiPM_Low_Position;

  vector<pair<Signal, Signal>> RearStrip_ASSOCIATED;
  vector<bool> Rear_Associated = vector<bool>(Rear_Position.size(), false);
  vector<bool> Strip_Associated = vector<bool>(Strip_Position.size(), false);

  for (int index = 0; index < signals.GetSize(); index++)
  {
    int current_label = signals[index].Label;
    if (IsDetectorSiliBack(current_label))
    {
      Rear_Position.push_back(index);
    }
    else if (IsDetectorSiliStrip(current_label))
    {
      Strip_Position.push_back(index);
    }
    else if (IsDetectorBetaHigh(current_label))
    {
      SiPM_High_Position.push_back(index);
    }
    else if (IsDetectorBetaLow(current_label))
    {
      SiPM_Low_Position.push_back(index);
    }
  }

  for (int index_rear = 0; index_rear < Rear_Position.size(); index_rear++)
  {
    for (int index_strip = 0; index_strip < Strip_Position.size(); index_strip++)
    {
      if (IsSameSiliDetector(signals[Strip_Position[index_strip]].Label, signals[Rear_Position[index_rear]].Label))
      {
        RearStrip_ASSOCIATED.push_back(make_pair(signals[Rear_Position[index_rear]], signals[Strip_Position[index_strip]]));
      }
    }
  }

  if (RearStrip_ASSOCIATED.size() == 1 && Rear_Position.size() == 1 && Strip_Position.size() == 1)
  {
    // Cutting
    double diff = (RearStrip_ASSOCIATED[0].first.Channel - RearStrip_ASSOCIATED[0].second.Channel) / RearStrip_ASSOCIATED[0].second.Channel;
    double channel = RearStrip_ASSOCIATED[0].second.Channel;
    double spread = SpreadAcceptance[RearStrip_ASSOCIATED[0].second.Label]->Eval(channel);
    double mean = MeanAcceptance[RearStrip_ASSOCIATED[0].second.Label];

    if (diff > mean - 3 * spread && diff < mean + 3 * spread)
    {
      H_Channel_Cutted[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[0].first.Channel);   // REAR
      H_Channel_Cutted[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].second.Channel); // STRIP

      H_RearStrip_Channel_Cutted[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[0].first.Channel, RearStrip_ASSOCIATED[0].second.Channel);  // REAR-STRIP for rear plot
      H_RearStrip_Channel_Cutted[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].first.Channel, RearStrip_ASSOCIATED[0].second.Channel); // REAR-STRIP for strip plot

      H_Strip_Channel_DiffRear_Cutted[RearStrip_ASSOCIATED[0].second.Label]->Fill((RearStrip_ASSOCIATED[0].first.Channel - RearStrip_ASSOCIATED[0].second.Channel) / RearStrip_ASSOCIATED[0].second.Channel);                                                // STRIP/REAR
      H_Strip_Channel_DiffRearvsStrip_Cutted[RearStrip_ASSOCIATED[0].second.Label]->Fill((RearStrip_ASSOCIATED[0].first.Channel - RearStrip_ASSOCIATED[0].second.Channel) / RearStrip_ASSOCIATED[0].second.Channel, RearStrip_ASSOCIATED[0].second.Channel); // STRIP/REAR vs Strip channel

      H_RearStrip_Time_Cutted[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[0].first.Time - RearStrip_ASSOCIATED[0].second.Time);

      for (int index_h = 0; index_h < SiPM_High_Position.size(); index_h++)
      {
        H_Channel_Cutted[signals[SiPM_High_Position[index_h]].Label]->Fill(signals[SiPM_High_Position[index_h]].Channel);
      }
      for (int index_l = 0; index_l < SiPM_Low_Position.size(); index_l++)
      {
        H_Channel_Cutted[signals[SiPM_Low_Position[index_l]].Label]->Fill(signals[SiPM_Low_Position[index_l]].Channel);
      }

      H_SiPMLow_Mulitplicity_Cutted->Fill(SiPM_Low_Position.size());
      H_SiPMHigh_Mulitplicity_Cutted->Fill(SiPM_High_Position.size());
    }
  }
}

#endif