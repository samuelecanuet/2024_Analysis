// #ifndef GROUPER_HH
// #define GROUPER_HH

#include "../../../lib/SignalDict/Signal.h"
#include "Detectors.hh"

//
using namespace std;

int VERBOSE = 0;
int Run;
bool FULL = false;

/// FILE ///
string ROOT_filename;
string ROOT_basefilename;
TFile *GROUPED_File;

/// TREE ///
TTree *GROUPED_Tree;
vector<Signal> GROUPED_Tree_Silicon;
vector<Signal> GROUPED_Tree_SiPMHigh;
vector<Signal> GROUPED_Tree_SiPMLow;

TTree *CUTTED_Tree;
vector<Signal> CUTTED_Tree_Silicon;
vector<Signal> CUTTED_Tree_SiPMHigh;
vector<Signal> CUTTED_Tree_SiPMLow;

TTree *CLEANED_Tree;
TTree *CLEANED_Tree_detector[SIGNAL_MAX];
vector<Signal> CLEANED_Tree_Silicon;
vector<vector<pair<Signal, Signal>>> CLEANED_Tree_SiPMGroup;
Signal CLEANED_Tree_HRS;
double Tree_Channel_detector;

TTreeReaderArray<Signal> *Silicon;
TTreeReaderArray<Signal> *SiPM_High;
TTreeReaderArray<Signal> *SiPM_Low;

////////////// HISTOGRAMS ////////////////

// Group
TH1D *H_Group_Time[SIGNAL_MAX];
TH2D *H_Group_ChannelTime[SIGNAL_MAX];
TH1D *H_Channel_Group[SIGNAL_MAX];
// RAW
TH1D *H_Strip_Mulitplicity_RAW;
TH1D *H_Rear_Mulitplicity_RAW;
TH1D *H_Channel_RAW[SIGNAL_MAX];
TH2D *H_RearStrip_Channel_RAW[SIGNAL_MAX];
TH1D *H_SiPMHigh_Mulitplicity_RAW;
TH1D *H_SiPMLow_Mulitplicity_RAW;

// PileUp
TH1D *H_Channel_PileUp[SIGNAL_MAX];

// CUTTING
TH1D *H_Channel_Cutted[SIGNAL_MAX];
TH2D *H_RearStrip_ChannelTime_Cutted[SIGNAL_MAX];
TH1D *H_RearSiPM_Time_Cutted[SIGNAL_MAX];
TH2D *H_RearSiPM_ChannelTime_Cutted[SIGNAL_MAX];
TH2D *H_SiPM_ChannelTime_Cutted[SIGNAL_MAX];
TH2D *H_RearSiPM_MulitplicityTime_Cutted[SIGNAL_MAX];
TH2D *H_SiPM_MultiplicityTime_Cutted[SIGNAL_MAX];
TH2D *H_RearStrip_Channel_Cutted[SIGNAL_MAX];
TGraph *P_RearStrip_Channel_Cutted[SIGNAL_MAX];
pair<double, double> RearStrip_Matching[SIGNAL_MAX];
TH1D *H_Strip_Channel_DiffRear_Cutted[SIGNAL_MAX];
TH2D *H_Strip_Channel_DiffRearvsStrip_Cutted[SIGNAL_MAX];
TH1D *H_RearStrip_Time_Cutted[SIGNAL_MAX];
TH1D *H_SiPMLow_Mulitplicity_Cutted;
TH1D *H_SiPMHigh_Mulitplicity_Cutted;
TCanvas *C_Strip_Channel_DiffRear[SIGNAL_MAX];
TGraphErrors *G_RearStrip_Spread[SIGNAL_MAX];

// Walk_Silicon
TCanvas *C_RearSiPM_ChannelTime_Walk_Silicon[SIGNAL_MAX];
TH2D *H_RearSiPM_ChannelTime_Walk_Silicon[SIGNAL_MAX];
TGraphErrors *G_Rear_Mean_Walk_Silicon[SIGNAL_MAX];

// Walk_SiPM
TCanvas *C_SiPM_ChannelTime_Walk_SiPM[SIGNAL_MAX];
TH2D *H_SiPM_ChannelTime_Walk_SiPM[SIGNAL_MAX];

// Walk Corrected
TH2D *H_RearSiPM_ChannelTime_Walk_Corrected[SIGNAL_MAX];
TH2D *H_SiPM_ChannelTime_Walk_Corrected[SIGNAL_MAX];

// Cleaned
TH1D *H_Channel_Cleaned[SIGNAL_MAX];
TH1D *H_RearMatched_Channel_Cleaned[SIGNAL_MAX];
TH1D *H_RearMatchedMean_Channel_Cleaned[SIGNAL_MAX];
TH2D *H_RearSiPM_ChannelTime_Cleaned[SIGNAL_MAX];
TH2D *H_SiPM_ChannelTime_Cleaned[SIGNAL_MAX];
TH1D *H_SiPMHigh_Mulitplicity_Cleaned;
TH1D *H_SiPMLow_Mulitplicity_Cleaned;
pair<int, int> Rear_IAS[9] = {make_pair(0, 0), make_pair(40000, 44500), make_pair(36000, 39500), make_pair(36000, 39500), make_pair(33000, 36500), make_pair(36000, 39000), make_pair(33000, 36500), make_pair(35500, 38500), make_pair(35000, 38500)};
pair<int, int> peaks_window_F[SIGNAL_MAX];
TH1D *H_Rear_Channel_IAS_Cleaned[SIGNAL_MAX];
TH1D *H_Rear_Channel_IAScoinc_Cleaned[SIGNAL_MAX];
TCanvas *C_IAS_Channel_Cleaned;
TH1D *H_SiPMLowHigh_Time_One[SIGNAL_MAX];
TH1D *H_SiPMLowHigh_Time_Multiple[SIGNAL_MAX];
TH2D *H_SiPMLowHigh_TimeChannel_One[SIGNAL_MAX];
TH1D *H_RearSiPM_Time_Nearest[SIGNAL_MAX];
TH1D *H_RearSiPM_Time_Nearest_Nearest[SIGNAL_MAX][SIGNAL_MAX];
TH2D *H_RearSiPM_Channel_Nearset_Nearest[SIGNAL_MAX][SIGNAL_MAX];
TH2D *H_RearSiPM_Channel_Nearest_Group[SIGNAL_MAX][SIGNAL_MAX];
TH1D *H_2SiPM_Time_Next[SIGNAL_MAX][SIGNAL_MAX];
TH1D *H_SiPM_Pair_Time[SIGNAL_MAX];
TH1D *H_RearSiPM_Time_New_Nearest[SIGNAL_MAX];
TH1D *H_RearSiPM_Time_New_Mean[SIGNAL_MAX];
TH1D *H_SiPM_Channel_Alone[SIGNAL_MAX];
TH1D *H_SiPM_Channel_Couple[SIGNAL_MAX];
TH2D *H_2SiPM_Channel_Couple[SIGNAL_MAX];
TH1D *H_Multiplicity_Groups;
TH1D *H_Multiplicity_SubGroups;
TH1D *H_Multiplicity_Groups_IAS;
TH1D* H_Multiplicity_SubGroups_IAS;
TH1D* H_Time_SubGroups[SIGNAL_MAX];

// DT
TH2D *H_RearSiPM_ChannelTime_dt[SIGNAL_MAX];
TH2D *H_RearSiPM_ChannelTime_dt_det[SIGNAL_MAX][SIGNAL_MAX];


// FITS
TCanvas *C_Strip_Cutting_Fits;
TCanvas *C_Rear_Walk_Silicon_Fits;
TCanvas *C_SiPM_Walk_SiPM_Fits;

// CASES
//  A
TH1D *H_Channel_A[SIGNAL_MAX];
TH2D *H_RearStrip_Channel_A[SIGNAL_MAX];
TH2D *H_RearStrip_ChannelTime_A[SIGNAL_MAX];
TH1D *H_Strip_Channel_DiffRear_A[SIGNAL_MAX];
TH2D *H_Strip_Channel_DiffRearvsStrip_A[SIGNAL_MAX];
TH1D *H_RearStrip_Time_A[SIGNAL_MAX];
// B
TH1D *H_Channel_B[SIGNAL_MAX];
TH2D *H_RearStrip_Channel_B[SIGNAL_MAX];
TH2D *H_Rear2Strip_Channel_B[SIGNAL_MAX];
TH2D *H_2Strip_Channel_B[SIGNAL_MAX];
TH1D *H_RearStrip_Time_B[SIGNAL_MAX];
// C
TH1D *H_Channel_C[SIGNAL_MAX];
TH2D *H_RearStrip_Channel_C[SIGNAL_MAX];
TH2D *H_2Strip_Channel_C[SIGNAL_MAX];
TH1D *H_RearStrip_Time_C[SIGNAL_MAX];
TH1D *H_2Strip_Time_C[SIGNAL_MAX];
TH2D *H_2Strip_Label_C;
// D
TH1D *H_Channel_D[SIGNAL_MAX];
TH2D *H_2Strip_Channel_D[SIGNAL_MAX];
TH1D *H_2Strip_Time_D[SIGNAL_MAX];
TH2D *H_2Strip_Label_D;
// E
TH1D *H_Channel_E[SIGNAL_MAX];
TH2D *H_2Strip_Channel_E[SIGNAL_MAX];
TH1D *H_2Strip_Time_E[SIGNAL_MAX];
TH2D *H_2Strip_Label_E;
// F
TH1D *H_Channel_F[SIGNAL_MAX];
// G
TH1D *H_Channel_G[SIGNAL_MAX];
TH2D *H_RearStrip_Channel_G;
TH1D *H_RearStrip_Time_G;
TH2D *H_RearStrip_Label_G;
// H
TH1D *H_Channel_H[SIGNAL_MAX];
// check
TH2D *H_2Strip_Label_Check;
//////////////////////////////////////

TGraphErrors *spread_up[SIGNAL_MAX];
TGraphErrors *spread_down[SIGNAL_MAX];
TGraphErrors *maximum[SIGNAL_MAX];

///// DIRECTORY //////////////////
// Group
TDirectory *dir_Group;
TDirectory *dir_Group_Strip;
TDirectory *dir_Group_Rear;
// Detector
TDirectory *dir_Strip;
TDirectory *dir_Strip_Channel;
TDirectory *dir_Strip_Time;
TDirectory *dir_Strip_Multiplicities;
TDirectory *dir_Strip_Channel_Detector[SIGNAL_MAX];
TDirectory *dir_Strip_Time_Detector[SIGNAL_MAX];

TDirectory *dir_Rear;
TDirectory *dir_Rear_Channel;
TDirectory *dir_Rear_Time;
TDirectory *dir_Rear_Multiplicity;
TDirectory *dir_Rear_Channel_Detector[SIGNAL_MAX];
TDirectory *dir_Rear_Time_Detector[SIGNAL_MAX];

TDirectory *dir_SiPM_High;
TDirectory *dir_SiPMHigh_Channel;
TDirectory *dir_SiPMHigh_Time;
TDirectory *dir_SiPMHigh_Multiplicity;
TDirectory *dir_SiPMHigh_Channel_Detector[SIGNAL_MAX];
TDirectory *dir_SiPMHigh_Time_Detector[SIGNAL_MAX];

TDirectory *dir_SiPMLow;
TDirectory *dir_SiPMLow_Channel;
TDirectory *dir_SiPMLow_Time;
TDirectory *dir_SiPMLow_Multiplicity;
TDirectory *dir_SiPMLow_Channel_Detector[SIGNAL_MAX];
TDirectory *dir_SiPMLow_Time_Detector[SIGNAL_MAX];

// Selection step
TDirectory *dir_CASES;
TDirectory *dir_CASES_Channel;

TDirectory *dir_Cutting;
TDirectory *dir_Cutting_Fitting;
TDirectory *dir_Cutting_SiPM;

TDirectory *dir_Walk_Silicon;

TDirectory *dir_Walk_SiPM;

TDirectory *dir_Cleaned;
TDirectory *dir_Cleaned_Multiplicity;
TDirectory *dir_Cleaned_HighSiPM_TimeCoincidence_Mulitplicity[SIGNAL_MAX];
TDirectory *dir_Cleaned_LowSiPM_TimeCoincidence_Mulitplicity[SIGNAL_MAX];
TDirectory *dir_Matching_Cleaned;

TDirectory *dir_Fits;
TDirectory *dirCoincidences;

// FASTER GROUP //
double lastGroupTime[SIGNAL_MAX];
double LastrealSignalTime[SIGNAL_MAX];
Signal LastrealSignal;

///////// FUCNTION ARRAY /////////////
TF1 *SpreadAcceptance[SIGNAL_MAX];
double MeanAcceptance[SIGNAL_MAX];
TF1 *MeanAcceptance_Walk_Silicon[SIGNAL_MAX];
TF1 *MeanAcceptance_Walk_SiPM[SIGNAL_MAX];

//////// calibration from Calibration.cc  ////////
double ManualCalibFitted[SIGNAL_MAX][3];

// COUNTER //
int counter_all_IAS = 0;
int counter_condition = 0;


///////// FITS FUNCTIONS /////////////
double skewedgauss(double *x, double *p)
{
  double xi = p[0];
  double omega = p[1];
  double alpha = p[2];
  double arg = (x[0] - xi) / omega;
  double smallphi = TMath::Gaus(arg, 0.0, 1.0, true);
  double bigphi = 0.5 * (1 + erf(alpha * arg / sqrt(2)));
  return p[4] + p[3] * 2. / omega * smallphi * bigphi;
}

double gausslorrentz(double *x, double *p)
{
  double amp_gauss = p[0];
  double mean_gauss = p[1];
  double sigma_gauss = p[2];
  double amp_lorentz = p[3];
  double mean_lorentz = p[4];
  double gamma_lorentz = p[5];
  double BKG = p[6];

  double gauss = amp_gauss * TMath::Gaus(x[0], mean_gauss, sigma_gauss);
  double lorentz = amp_lorentz * TMath::BreitWigner(x[0], mean_lorentz, gamma_lorentz);

  return gauss + lorentz + BKG;
}

void InitPeakWindow(string nuclei)
{
  string CalibFileName;

  CalibFileName = "./Config_Files/Win_32Ar_Catcher1_14.txt";

  ifstream fileF(CalibFileName);

  for (auto &value : peaks_window_F)
  {
    value = make_pair(0, 0);
  }

  if (fileF.is_open())
  {
    Info("Window file found");

    string line;
    while (getline(fileF, line))
    {
      istringstream iss(line);
      int min;
      int max;
      int min1 = 0;
      int max1 = 0;
      string DetName;
      iss >> DetName >> min >> max;

      for (size_t i = 0; i < detectorNum; ++i)
      {
        if (DetName == detectorName[i])
        {
          peaks_window_F[i] = make_pair(min, max);
          break;
        }
      }
    }
  }
  else
  {
    Error("No Window file found");
  }
}

void InitCalibration()
{
  string CalibFileName;

  CalibFileName = "Config_Files/Calibration_" + to_string(YEAR) + ".txt";

  ifstream file(CalibFileName);

  if (file.is_open())
  {
    Info("Calibration file found");

    string line;
    while (getline(file, line))
    {
      istringstream iss(line);
      int det;
      double a;
      double b;
      double c;
      iss >> det >> a >> b >> c;
      ManualCalibFitted[det][0] = a;
      ManualCalibFitted[det][1] = b;
      ManualCalibFitted[det][2] = c;
    }
  }
  else
  {
    Error("No Calibration file found");
  }
}

void SaveFitParameters()
{ 
  string add = "_2024";
  if (FLAG2021) add = "_2021";
  if (FLAG2025) add = "_2025";
  TFile *file = MyTFile((DIR_ROOT_DATA_GROUPED + "Grouping_FitParameters"+add+".root").c_str(), "RECREATE");

  TGraph *MeanAcceptance_Graph = new TGraph();
  for (int i = 0; i < detectorNum; i++)
  {
    if (IsDetectorBeta(i))
    {
      MeanAcceptance_Walk_SiPM[i]->SetName(("MeanAcceptance_Walk_SiPM_" + detectorName[i]).c_str());
      MeanAcceptance_Walk_SiPM[i]->Write();
    }
    else if (IsDetectorSiliBack(i))
    {
      MeanAcceptance_Walk_Silicon[i]->SetName(("MeanAcceptance_Walk_Silicon_" + detectorName[i]).c_str());
      MeanAcceptance_Walk_Silicon[i]->Write();
    }
    else if (IsDetectorSiliStrip(i))
    {
      SpreadAcceptance[i]->SetName(("SpreadAcceptance_" + detectorName[i]).c_str());
      SpreadAcceptance[i]->Write();

      MeanAcceptance_Graph->SetPoint(i, i, MeanAcceptance[i]);
    }
  }
  MeanAcceptance_Graph->SetName("MeanAcceptance");
  MeanAcceptance_Graph->Write();

  file->Close();
}

void LoadFitParameters()
{
  TFile *file = MyTFile((DIR_ROOT_DATA_GROUPED + "Grouping_FitParameters_"+ to_string(YEAR) +".root").c_str(), "READ");
  for (int i = 0; i < detectorNum; i++)
  {
    if (IsDetectorBeta(i))
    {
      MeanAcceptance_Walk_SiPM[i] = (TF1 *)file->Get(("MeanAcceptance_Walk_SiPM_" + detectorName[i]).c_str());
      if (MeanAcceptance_Walk_SiPM[i] == nullptr)
      {
        Error(("MeanAcceptance_Walk_SiPM_" + detectorName[i]));
      }
    }
    else if (IsDetectorSiliBack(i))
    {
      MeanAcceptance_Walk_Silicon[i] = (TF1 *)file->Get(("MeanAcceptance_Walk_Silicon_" + detectorName[i]).c_str());
      if (MeanAcceptance_Walk_Silicon[i] == nullptr)
      {
        Error(("MeanAcceptance_Walk_Silicon_" + detectorName[i]));
      }
    }
    else if (IsDetectorSiliStrip(i))
    {
      SpreadAcceptance[i] = (TF1 *)file->Get(("SpreadAcceptance_" + detectorName[i]).c_str());
      MeanAcceptance[i] = ((TGraph *)file->Get("MeanAcceptance"))->GetY()[i];
      if (SpreadAcceptance[i] == nullptr)
      {
        Error(("SpreadAcceptance_" + detectorName[i]).c_str());
      }
      if (MeanAcceptance[i] == 0.)
      {
        Warning(("MeanAcceptance_" + detectorName[i]).c_str());
      }
    }
  }
  file->Close();
  GROUPED_File->cd();
  Success("Fit parameters loaded");
}

//////////////////////////////////////

inline int InitHistograms_Grouped()
{
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);

  // Group
  dir_Group = GROUPED_File->mkdir("Group");
  dir_Group_Strip = dir_Group->mkdir("Strip");
  dir_Group_Rear = dir_Group->mkdir("Rear");
  /// STRIP
  dir_Strip = GROUPED_File->mkdir("Strip");
  dir_Strip_Channel = dir_Strip->mkdir("Strip_Channel");
  dir_Strip_Time = dir_Strip->mkdir("Strip_Time");
  dir_Strip_Multiplicities = dir_Strip->mkdir("Strip_Multiplicities");

  ////REAR
  dir_Rear = GROUPED_File->mkdir("Rear");
  dir_Rear_Channel = dir_Rear->mkdir("Rear_Channel");
  dir_Rear_Time = dir_Rear->mkdir("Rear_Time");
  dir_Rear_Multiplicity = dir_Rear->mkdir("Rear_Multiplicity");

  ////SiPM High
  dir_SiPM_High = GROUPED_File->mkdir("SiPM_High");
  dir_SiPMHigh_Channel = dir_SiPM_High->mkdir("SiPMHigh_Channel");
  dir_SiPMHigh_Time = dir_SiPM_High->mkdir("SiPMHigh_Time");
  dir_SiPMHigh_Multiplicity = dir_SiPM_High->mkdir("SiPMHigh_Multiplicity");

  ////SiPM Low
  dir_SiPMLow = GROUPED_File->mkdir("SiPMLow");
  dir_SiPMLow_Channel = dir_SiPMLow->mkdir("SiPMLow_Channel");
  dir_SiPMLow_Time = dir_SiPMLow->mkdir("SiPMLow_Time");
  dir_SiPMLow_Multiplicity = dir_SiPMLow->mkdir("SiPMLow_Multiplicity");

  /// CASES
  dir_CASES = GROUPED_File->mkdir("CASES");
  dir_CASES_Channel = dir_CASES->mkdir("CASES_Channel");

  /// CUTTING
  dir_Cutting = GROUPED_File->mkdir("Cutting");
  dir_Cutting_Fitting = dir_Cutting->mkdir("Fitting");
  dir_Cutting_SiPM = dir_Cutting->mkdir("SiPM");

  /// Walk_Silicon
  dir_Walk_Silicon = GROUPED_File->mkdir("Walk_Silicon");

  /// Walk_SiPM
  dir_Walk_SiPM = GROUPED_File->mkdir("Walk_SiPM");

  // Cleaned
  dir_Cleaned = GROUPED_File->mkdir("Cleaned");
  dir_Matching_Cleaned = dir_Cleaned->mkdir("Matching");
  C_IAS_Channel_Cleaned = new TCanvas("C_IAS_Channel_Cleaned", "C_IAS_Channel_Cleaned", 800, 800);
  C_IAS_Channel_Cleaned->Divide(4, 2);

  // FIT RESULTS
  dir_Fits = GROUPED_File->mkdir("Fits");
  C_Strip_Cutting_Fits = new TCanvas("C_Strip_Cutting_Fits", "C_Strip_Cutting_Fits", 800, 800);
  C_Strip_Cutting_Fits->Divide(SILI_SIZE - 1, SILI_NUM);
  C_Rear_Walk_Silicon_Fits = new TCanvas("C_Rear_Walk_Silicon_Fits", "C_Rear_Walk_Silicon_Fits", 800, 800);
  C_Rear_Walk_Silicon_Fits->Divide(4, 2);
  C_SiPM_Walk_SiPM_Fits = new TCanvas("C_SiPM_Walk_SiPM_Fits", "C_SiPM_Walk_SiPM_Fits", 800, 800);
  C_SiPM_Walk_SiPM_Fits->Divide(4, 2);

  // Coincidences
  dirCoincidences = GROUPED_File->mkdir("Coincidences");

  for (size_t i = 0; i < detectorNum; ++i)
  {
    if (IsDetectorSiliStrip(i))
    {
      dir_Strip_Channel_Detector[i] = dir_Strip_Channel->mkdir(detectorName[i].c_str());
      dir_Strip_Time_Detector[i] = dir_Strip_Time->mkdir(detectorName[i].c_str());

      H_Group_Time[i] = new TH1D(("Group_Time_" + detectorName[i]).c_str(), ("Group_Time_" + detectorName[i]).c_str(), 1e4, 0, 1e9); // 1 ms
      H_Group_Time[i]->GetXaxis()->SetTitle("Time (ns)");
      H_Group_Time[i]->GetYaxis()->SetTitle("Counts");
      H_Group_Time[i]->GetXaxis()->CenterTitle();
      H_Group_Time[i]->GetYaxis()->CenterTitle();

      H_Group_ChannelTime[i] = new TH2D(("Group_ChannelTime_" + detectorName[i]).c_str(), ("Group_ChannelTime_" + detectorName[i]).c_str(), 1e4, 0, 1e8, 1.5e3, 35e3, 50e3);
      H_Group_ChannelTime[i]->GetXaxis()->SetTitle("Time (ns)");
      H_Group_ChannelTime[i]->GetYaxis()->SetTitle("Channel");
      H_Group_ChannelTime[i]->GetXaxis()->CenterTitle();
      H_Group_ChannelTime[i]->GetYaxis()->CenterTitle();

      H_Channel_Group[i] = new TH1D(("Channel_Group_" + detectorName[i]).c_str(), ("Channel_Group_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_Group[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_Group[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_Group[i]->GetXaxis()->CenterTitle();
      H_Channel_Group[i]->GetYaxis()->CenterTitle();

      H_Channel_RAW[i] = new TH1D(("Channel_RAW_" + detectorName[i]).c_str(), ("Channel_RAW_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_RAW[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_RAW[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_RAW[i]->GetXaxis()->CenterTitle();
      H_Channel_RAW[i]->GetYaxis()->CenterTitle();

      H_Channel_PileUp[i] = new TH1D(("Channel_PileUp_" + detectorName[i]).c_str(), ("Channel_PileUp_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_PileUp[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_PileUp[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_PileUp[i]->GetXaxis()->CenterTitle();
      H_Channel_PileUp[i]->GetYaxis()->CenterTitle();

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

      H_Channel_G[i] = new TH1D(("Channel_G_" + detectorName[i]).c_str(), ("Channel_G_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_G[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_G[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_G[i]->GetXaxis()->CenterTitle();
      H_Channel_G[i]->GetYaxis()->CenterTitle();

      H_Channel_H[i] = new TH1D(("Channel_H_" + detectorName[i]).c_str(), ("Channel_H_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_H[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_H[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_H[i]->GetXaxis()->CenterTitle();
      H_Channel_H[i]->GetYaxis()->CenterTitle();

      H_Channel_Cleaned[i] = new TH1D(("Channel_Cleaned_" + detectorName[i]).c_str(), ("Channel_Cleaned_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_Cleaned[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_Cleaned[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_Cleaned[i]->GetXaxis()->CenterTitle();
      H_Channel_Cleaned[i]->GetYaxis()->CenterTitle();

      H_Channel_Cutted[i] = new TH1D(("Channel_Cutted_" + detectorName[i]).c_str(), ("Channel_Cutted_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_Cutted[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_Cutted[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_Cutted[i]->GetXaxis()->CenterTitle();
      H_Channel_Cutted[i]->GetYaxis()->CenterTitle();

      H_RearMatched_Channel_Cleaned[i] = new TH1D(("RearMatched_Channel_Cleaned_" + detectorName[i]).c_str(), ("RearMatched_Channel_Cleaned_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_RearMatched_Channel_Cleaned[i]->GetXaxis()->SetTitle("Channel");
      H_RearMatched_Channel_Cleaned[i]->GetYaxis()->SetTitle("Counts");
      H_RearMatched_Channel_Cleaned[i]->GetXaxis()->CenterTitle();
      H_RearMatched_Channel_Cleaned[i]->GetYaxis()->CenterTitle();

      H_RearMatchedMean_Channel_Cleaned[i] = new TH1D(("RearMatchedMean_Channel_Cleaned_" + detectorName[i]).c_str(), ("RearMatchedMean_Channel_Cleaned_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_RearMatchedMean_Channel_Cleaned[i]->GetXaxis()->SetTitle("Channel");
      H_RearMatchedMean_Channel_Cleaned[i]->GetYaxis()->SetTitle("Counts");
      H_RearMatchedMean_Channel_Cleaned[i]->GetXaxis()->CenterTitle();
      H_RearMatchedMean_Channel_Cleaned[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Channel_RAW[i] = new TH2D(("RearStrip_Channel_RAW_" + detectorName[i]).c_str(), ("RearStrip_Channel_RAW_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearStrip_Channel_RAW[i]->GetXaxis()->SetTitle("Rear Channel");
      H_RearStrip_Channel_RAW[i]->GetYaxis()->SetTitle("Strip Channel");
      H_RearStrip_Channel_RAW[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Channel_RAW[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Channel_A[i] = new TH2D(("RearStrip_Channel_A_" + detectorName[i]).c_str(), ("RearStrip_Channel_A_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearStrip_Channel_A[i]->GetXaxis()->SetTitle("Rear Channel");
      H_RearStrip_Channel_A[i]->GetYaxis()->SetTitle("Strip Channel");
      H_RearStrip_Channel_A[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Channel_A[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Channel_B[i] = new TH2D(("RearStrip_Channel_B_" + detectorName[i]).c_str(), ("RearStrip_Channel_B_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearStrip_Channel_B[i]->GetXaxis()->SetTitle("Rear Channel");
      H_RearStrip_Channel_B[i]->GetYaxis()->SetTitle("Strip Channel");
      H_RearStrip_Channel_B[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Channel_B[i]->GetYaxis()->CenterTitle();

      H_2Strip_Channel_B[i] = new TH2D(("2Strip_Channel_B_" + detectorName[i]).c_str(), ("2Strip_Channel_B_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_2Strip_Channel_B[i]->GetXaxis()->SetTitle("Rear Channel");
      H_2Strip_Channel_B[i]->GetYaxis()->SetTitle("Strip Channel");
      H_2Strip_Channel_B[i]->GetXaxis()->CenterTitle();
      H_2Strip_Channel_B[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Channel_C[i] = new TH2D(("RearStrip_Channel_C_" + detectorName[i]).c_str(), ("RearStrip_Channel_C_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearStrip_Channel_C[i]->GetXaxis()->SetTitle("Rear Channel");
      H_RearStrip_Channel_C[i]->GetYaxis()->SetTitle("Strip Channel");
      H_RearStrip_Channel_C[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Channel_C[i]->GetYaxis()->CenterTitle();

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

      H_RearStrip_Time_A[i] = new TH1D(("RearStrip_Time_A_" + detectorName[i]).c_str(), ("RearStrip_Time_A_" + detectorName[i]).c_str(), winGroupN_Sili, winGroupMin_Sili, winGroupMax_Sili);
      H_RearStrip_Time_A[i]->GetXaxis()->SetTitle("Rear-Strip Time (ns)");
      H_RearStrip_Time_A[i]->GetYaxis()->SetTitle("Counts");
      H_RearStrip_Time_A[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Time_A[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Time_B[i] = new TH1D(("RearStrip_Time_B_" + detectorName[i]).c_str(), ("RearStrip_Time_B_" + detectorName[i]).c_str(), winGroupN_Sili, winGroupMin_Sili, winGroupMax_Sili);
      H_RearStrip_Time_B[i]->GetXaxis()->SetTitle("Rear-Strip Time (ns)");
      H_RearStrip_Time_B[i]->GetYaxis()->SetTitle("Counts");
      H_RearStrip_Time_B[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Time_B[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Time_C[i] = new TH1D(("RearStrip_Time_C_" + detectorName[i]).c_str(), ("RearStrip_Time_C_" + detectorName[i]).c_str(), winGroupN_Sili, winGroupMin_Sili, winGroupMax_Sili);
      H_RearStrip_Time_C[i]->GetXaxis()->SetTitle("Rear-Strip Time (ns)");
      H_RearStrip_Time_C[i]->GetYaxis()->SetTitle("Counts");
      H_RearStrip_Time_C[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Time_C[i]->GetYaxis()->CenterTitle();

      H_2Strip_Time_C[i] = new TH1D(("2Strip_Time_C_" + detectorName[i]).c_str(), ("2Strip_Time_C_" + detectorName[i]).c_str(), winGroupN_Sili, winGroupMin_Sili, winGroupMax_Sili);
      H_2Strip_Time_C[i]->GetXaxis()->SetTitle("Time between 2 strips (ns)");
      H_2Strip_Time_C[i]->GetYaxis()->SetTitle("Counts");
      H_2Strip_Time_C[i]->GetXaxis()->CenterTitle();
      H_2Strip_Time_C[i]->GetYaxis()->CenterTitle();

      H_2Strip_Time_D[i] = new TH1D(("2Strip_Time_D_" + detectorName[i]).c_str(), ("2Strip_Time_D_" + detectorName[i]).c_str(), winGroupN_Sili, winGroupMin_Sili, winGroupMax_Sili);
      H_2Strip_Time_D[i]->GetXaxis()->SetTitle("Time between 2 strips (ns)");
      H_2Strip_Time_D[i]->GetYaxis()->SetTitle("Counts");
      H_2Strip_Time_D[i]->GetXaxis()->CenterTitle();
      H_2Strip_Time_D[i]->GetYaxis()->CenterTitle();

      H_2Strip_Time_E[i] = new TH1D(("2Strip_Time_E_" + detectorName[i]).c_str(), ("2Strip_Time_E_" + detectorName[i]).c_str(), winGroupN_Sili, winGroupMin_Sili, winGroupMax_Sili);
      H_2Strip_Time_E[i]->GetXaxis()->SetTitle("Time between 2 strips (ns)");
      H_2Strip_Time_E[i]->GetYaxis()->SetTitle("Counts");
      H_2Strip_Time_E[i]->GetXaxis()->CenterTitle();
      H_2Strip_Time_E[i]->GetYaxis()->CenterTitle();

      H_2Strip_Channel_C[i] = new TH2D(("2Strip_Channel_C_" + detectorName[i]).c_str(), ("2Strip_Channel_C_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_2Strip_Channel_C[i]->GetXaxis()->SetTitle("Strip Channel");
      H_2Strip_Channel_C[i]->GetYaxis()->SetTitle("Strip Channel");
      H_2Strip_Channel_C[i]->GetXaxis()->CenterTitle();
      H_2Strip_Channel_C[i]->GetYaxis()->CenterTitle();

      H_2Strip_Channel_D[i] = new TH2D(("2Strip_Channel_D_" + detectorName[i]).c_str(), ("2Strip_Channel_D_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_2Strip_Channel_D[i]->GetXaxis()->SetTitle("Strip grouped Channel");
      H_2Strip_Channel_D[i]->GetYaxis()->SetTitle("Strip Channel");
      H_2Strip_Channel_D[i]->GetXaxis()->CenterTitle();
      H_2Strip_Channel_D[i]->GetYaxis()->CenterTitle();

      H_2Strip_Channel_E[i] = new TH2D(("2Strip_Channel_E_" + detectorName[i]).c_str(), ("2Strip_Channel_E_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_2Strip_Channel_E[i]->GetXaxis()->SetTitle("Strip grouped Channel");
      H_2Strip_Channel_E[i]->GetYaxis()->SetTitle("Strip Channel");
      H_2Strip_Channel_E[i]->GetXaxis()->CenterTitle();
      H_2Strip_Channel_E[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Time_Cutted[i] = new TH1D(("RearStrip_Time_Cutted_" + detectorName[i]).c_str(), ("RearStrip_Time_Cutted_" + detectorName[i]).c_str(), winGroupN_Sili, winGroupMin_Sili, winGroupMax_Sili);
      H_RearStrip_Time_Cutted[i]->GetXaxis()->SetTitle("Rear-Strip Time (ns)");
      H_RearStrip_Time_Cutted[i]->GetYaxis()->SetTitle("Counts");
      H_RearStrip_Time_Cutted[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Time_Cutted[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Channel_Cutted[i] = new TH2D(("RearStrip_Channel_Cutted_" + detectorName[i]).c_str(), ("RearStrip_Channel_Cutted_" + detectorName[i]).c_str(), 1000, 30e3, 50e3, 1000, 40e3, 50e3);
      H_RearStrip_Channel_Cutted[i]->GetXaxis()->SetTitle("Rear Channel");
      H_RearStrip_Channel_Cutted[i]->GetYaxis()->SetTitle("Strip Channel");
      H_RearStrip_Channel_Cutted[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Channel_Cutted[i]->GetYaxis()->CenterTitle();

      P_RearStrip_Channel_Cutted[i] = new TGraph();
      P_RearStrip_Channel_Cutted[i]->GetXaxis()->SetTitle("Rear Channel");
      P_RearStrip_Channel_Cutted[i]->GetYaxis()->SetTitle("Strip Channel");
      P_RearStrip_Channel_Cutted[i]->GetXaxis()->CenterTitle();
      P_RearStrip_Channel_Cutted[i]->GetYaxis()->CenterTitle();

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

      H_RearStrip_ChannelTime_A[i] = new TH2D(("RearStrip_ChannelTime_A_" + detectorName[i]).c_str(), ("RearStrip_ChannelTime_A_" + detectorName[i]).c_str(), winGroupN_Sili, winGroupMin_Sili, winGroupMax_Sili, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearStrip_ChannelTime_A[i]->GetXaxis()->SetTitle("Rear-Strip Time (ns)");
      H_RearStrip_ChannelTime_A[i]->GetYaxis()->SetTitle("Strip Channel");
      H_RearStrip_ChannelTime_A[i]->GetXaxis()->CenterTitle();
      H_RearStrip_ChannelTime_A[i]->GetYaxis()->CenterTitle();

      H_RearStrip_ChannelTime_Cutted[i] = new TH2D(("RearStrip_ChannelTime_Cutted_" + detectorName[i]).c_str(), ("RearStrip_ChannelTime_Cutted_" + detectorName[i]).c_str(), winGroupN_Sili, winGroupMin_Sili, winGroupMax_Sili, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearStrip_ChannelTime_Cutted[i]->GetXaxis()->SetTitle("Rear-Strip Time (ns)");
      H_RearStrip_ChannelTime_Cutted[i]->GetYaxis()->SetTitle("Strip Channel");
      H_RearStrip_ChannelTime_Cutted[i]->GetXaxis()->CenterTitle();
      H_RearStrip_ChannelTime_Cutted[i]->GetYaxis()->CenterTitle();

      H_RearSiPM_ChannelTime_dt[i] = new TH2D(("StripSiPM_ChannelTime_dt_" + detectorName[i]).c_str(), ("StripSiPM_ChannelTime_dt_" + detectorName[i]).c_str(), 2e3, 6e3, 10e3, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearSiPM_ChannelTime_dt[i]->GetXaxis()->SetTitle("Strip SiPM Time (ns)");
      H_RearSiPM_ChannelTime_dt[i]->GetYaxis()->SetTitle("Strip Channel");
      H_RearSiPM_ChannelTime_dt[i]->GetXaxis()->CenterTitle();
      H_RearSiPM_ChannelTime_dt[i]->GetYaxis()->CenterTitle();
    }

    if (IsDetectorSiliBack(i))
    {

      dir_Rear_Channel_Detector[i] = dir_Rear_Channel->mkdir(detectorName[i].c_str());
      dir_Rear_Time_Detector[i] = dir_Rear_Time->mkdir(detectorName[i].c_str());

      H_Group_Time[i] = new TH1D(("Group_Time_" + detectorName[i]).c_str(), ("Group_Time_" + detectorName[i]).c_str(), 1e4, 0, 1e9); // 1 ms            H_Group_Time[i]->GetXaxis()->SetTitle("Time (ns)");
      H_Group_Time[i]->GetYaxis()->SetTitle("Counts");
      H_Group_Time[i]->GetXaxis()->CenterTitle();
      H_Group_Time[i]->GetYaxis()->CenterTitle();

      H_Group_ChannelTime[i] = new TH2D(("Group_ChannelTime_" + detectorName[i]).c_str(), ("Group_ChannelTime_" + detectorName[i]).c_str(), 1e4, 0, 1e9, eSiliN / 10, 20e3, 60e3);
      H_Group_ChannelTime[i]->GetYaxis()->SetTitle("Channel");
      H_Group_ChannelTime[i]->GetXaxis()->SetTitle("Time (ns)");
      H_Group_ChannelTime[i]->GetXaxis()->CenterTitle();
      H_Group_ChannelTime[i]->GetYaxis()->CenterTitle();

      H_Channel_Group[i] = new TH1D(("Channel_Group_" + detectorName[i]).c_str(), ("Channel_Group_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_Group[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_Group[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_Group[i]->GetXaxis()->CenterTitle();
      H_Channel_Group[i]->GetYaxis()->CenterTitle();

      H_Channel_RAW[i] = new TH1D(("Channel_RAW_" + detectorName[i]).c_str(), ("Channel_RAW_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_RAW[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_RAW[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_RAW[i]->GetXaxis()->CenterTitle();
      H_Channel_RAW[i]->GetYaxis()->CenterTitle();

      H_Channel_PileUp[i] = new TH1D(("Channel_PileUp_" + detectorName[i]).c_str(), ("Channel_PileUp_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_PileUp[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_PileUp[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_PileUp[i]->GetXaxis()->CenterTitle();
      H_Channel_PileUp[i]->GetYaxis()->CenterTitle();

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

      H_Channel_F[i] = new TH1D(("Channel_F_" + detectorName[i]).c_str(), ("Channel_F_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_F[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_F[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_F[i]->GetXaxis()->CenterTitle();
      H_Channel_F[i]->GetYaxis()->CenterTitle();

      H_Channel_G[i] = new TH1D(("Channel_G_" + detectorName[i]).c_str(), ("Channel_G_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_G[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_G[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_G[i]->GetXaxis()->CenterTitle();
      H_Channel_G[i]->GetYaxis()->CenterTitle();

      H_Channel_H[i] = new TH1D(("Channel_H_" + detectorName[i]).c_str(), ("Channel_H_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_H[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_H[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_H[i]->GetXaxis()->CenterTitle();
      H_Channel_H[i]->GetYaxis()->CenterTitle();

      H_Channel_Cutted[i] = new TH1D(("Channel_Cutted_" + detectorName[i]).c_str(), ("Channel_Cutted_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_Cutted[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_Cutted[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_Cutted[i]->GetXaxis()->CenterTitle();
      H_Channel_Cutted[i]->GetYaxis()->CenterTitle();

      H_Channel_Cleaned[i] = new TH1D(("Channel_Cleaned_" + detectorName[i]).c_str(), ("Channel_Cleaned_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      H_Channel_Cleaned[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_Cleaned[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_Cleaned[i]->GetXaxis()->CenterTitle();
      H_Channel_Cleaned[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Channel_RAW[i] = new TH2D(("RearStrip_Channel_RAW_" + detectorName[i]).c_str(), ("RearStrip_Channel_RAW_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearStrip_Channel_RAW[i]->GetXaxis()->SetTitle("Rear Channel");
      H_RearStrip_Channel_RAW[i]->GetYaxis()->SetTitle("Strip Channel");
      H_RearStrip_Channel_RAW[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Channel_RAW[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Channel_A[i] = new TH2D(("RearStrip_Channel_A_" + detectorName[i]).c_str(), ("RearStrip_Channel_A_" + detectorName[i]).c_str(), eSiliN / 10, 30e3, 40e3, eSiliN / 10, 30e3, 50e3);
      H_RearStrip_Channel_A[i]->GetXaxis()->SetTitle("Rear Channel");
      H_RearStrip_Channel_A[i]->GetYaxis()->SetTitle("Strip Channel");
      H_RearStrip_Channel_A[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Channel_A[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Channel_B[i] = new TH2D(("RearStrip_Channel_B_" + detectorName[i]).c_str(), ("RearStrip_Channel_B_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearStrip_Channel_B[i]->GetXaxis()->SetTitle("Rear Channel");
      H_RearStrip_Channel_B[i]->GetYaxis()->SetTitle("Strip Channel");
      H_RearStrip_Channel_B[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Channel_B[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Channel_C[i] = new TH2D(("RearStrip_Channel_C_" + detectorName[i]).c_str(), ("RearStrip_Channel_C_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearStrip_Channel_C[i]->GetXaxis()->SetTitle("Rear Channel");
      H_RearStrip_Channel_C[i]->GetYaxis()->SetTitle("Strip Channel");
      H_RearStrip_Channel_C[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Channel_C[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Channel_Cutted[i] = new TH2D(("RearStrip_Channel_Cutted_" + detectorName[i]).c_str(), ("RearStrip_Channel_Cutted_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearStrip_Channel_Cutted[i]->GetXaxis()->SetTitle("Rear Channel");
      H_RearStrip_Channel_Cutted[i]->GetYaxis()->SetTitle("Strip Channel");
      H_RearStrip_Channel_Cutted[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Channel_Cutted[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Time_A[i] = new TH1D(("RearStrip_Time_A_" + detectorName[i]).c_str(), ("RearStrip_Time_A_" + detectorName[i]).c_str(), winGroupN_Sili, winGroupMin_Sili, winGroupMax_Sili);
      H_RearStrip_Time_A[i]->GetXaxis()->SetTitle("Rear-Strip Time (ns)");
      H_RearStrip_Time_A[i]->GetYaxis()->SetTitle("Counts");
      H_RearStrip_Time_A[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Time_A[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Time_B[i] = new TH1D(("RearStrip_Time_B_" + detectorName[i]).c_str(), ("RearStrip_Time_B_" + detectorName[i]).c_str(), winGroupN_Sili, winGroupMin_Sili, winGroupMax_Sili);
      H_RearStrip_Time_B[i]->GetXaxis()->SetTitle("Rear-Strip Time (ns)");
      H_RearStrip_Time_B[i]->GetYaxis()->SetTitle("Counts");
      H_RearStrip_Time_B[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Time_B[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Time_C[i] = new TH1D(("RearStrip_Time_C_" + detectorName[i]).c_str(), ("RearStrip_Time_C_" + detectorName[i]).c_str(), winGroupN_Sili, winGroupMin_Sili, winGroupMax_Sili);
      H_RearStrip_Time_C[i]->GetXaxis()->SetTitle("Rear-Strip Time (ns)");
      H_RearStrip_Time_C[i]->GetYaxis()->SetTitle("Counts");
      H_RearStrip_Time_C[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Time_C[i]->GetYaxis()->CenterTitle();

      H_RearStrip_Time_Cutted[i] = new TH1D(("RearStrip_Time_Cutted_" + detectorName[i]).c_str(), ("RearStrip_Time_Cutted_" + detectorName[i]).c_str(), winGroupN_Sili, winGroupMin_Sili, winGroupMax_Sili);
      H_RearStrip_Time_Cutted[i]->GetXaxis()->SetTitle("Rear-Strip Time (ns)");
      H_RearStrip_Time_Cutted[i]->GetYaxis()->SetTitle("Counts");
      H_RearStrip_Time_Cutted[i]->GetXaxis()->CenterTitle();
      H_RearStrip_Time_Cutted[i]->GetYaxis()->CenterTitle();

      H_Rear2Strip_Channel_B[i] = new TH2D(("Rear2Strip_Channel_B_" + detectorName[i]).c_str(), ("Rear2Strip_Channel_B_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_Rear2Strip_Channel_B[i]->GetXaxis()->SetTitle("Rear Channel");
      H_Rear2Strip_Channel_B[i]->GetYaxis()->SetTitle("Sum 2 Strip Channel");
      H_Rear2Strip_Channel_B[i]->GetXaxis()->CenterTitle();
      H_Rear2Strip_Channel_B[i]->GetYaxis()->CenterTitle();

      H_2Strip_Channel_D[i] = new TH2D(("2Strip_Channel_D_" + detectorName[i]).c_str(), ("2Strip_Channel_D_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_2Strip_Channel_D[i]->GetXaxis()->SetTitle("Strip grouped Channel");
      H_2Strip_Channel_D[i]->GetYaxis()->SetTitle("Strip Channel");
      H_2Strip_Channel_D[i]->GetXaxis()->CenterTitle();
      H_2Strip_Channel_D[i]->GetYaxis()->CenterTitle();

      H_2Strip_Channel_E[i] = new TH2D(("2Strip_Channel_E_" + detectorName[i]).c_str(), ("2Strip_Channel_E_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      H_2Strip_Channel_E[i]->GetXaxis()->SetTitle("Strip grouped Channel");
      H_2Strip_Channel_E[i]->GetYaxis()->SetTitle("Strip Channel");
      H_2Strip_Channel_E[i]->GetXaxis()->CenterTitle();
      H_2Strip_Channel_E[i]->GetYaxis()->CenterTitle();

      H_2Strip_Time_D[i] = new TH1D(("2Strip_Time_D_" + detectorName[i]).c_str(), ("2Strip_Time_D_" + detectorName[i]).c_str(), winGroupN_Sili, winGroupMin_Sili, winGroupMax_Sili);
      H_2Strip_Time_D[i]->GetXaxis()->SetTitle("Time between 2 strips (ns)");
      H_2Strip_Time_D[i]->GetYaxis()->SetTitle("Counts");
      H_2Strip_Time_D[i]->GetXaxis()->CenterTitle();
      H_2Strip_Time_D[i]->GetYaxis()->CenterTitle();

      H_2Strip_Time_E[i] = new TH1D(("2Strip_Time_E_" + detectorName[i]).c_str(), ("2Strip_Time_E_" + detectorName[i]).c_str(), winGroupN_Sili, winGroupMin_Sili, winGroupMax_Sili);
      H_2Strip_Time_E[i]->GetXaxis()->SetTitle("Time between 2 strips (ns)");
      H_2Strip_Time_E[i]->GetYaxis()->SetTitle("Counts");
      H_2Strip_Time_E[i]->GetXaxis()->CenterTitle();
      H_2Strip_Time_E[i]->GetYaxis()->CenterTitle();

      H_RearSiPM_Time_Cutted[i] = new TH1D(("RearSiPM_Time_Cutted_" + detectorName[i]).c_str(), ("RearSiPM_Time_Cutted_" + detectorName[i]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
      H_RearSiPM_Time_Cutted[i]->GetXaxis()->SetTitle("Rear-SiPM Time (ns)");
      H_RearSiPM_Time_Cutted[i]->GetYaxis()->SetTitle("Counts");
      H_RearSiPM_Time_Cutted[i]->GetXaxis()->CenterTitle();
      H_RearSiPM_Time_Cutted[i]->GetYaxis()->CenterTitle();

      H_RearSiPM_ChannelTime_Cutted[i] = new TH2D(("RearSiPM_ChannelTime_Cutted_" + detectorName[i]).c_str(), ("RearSiPM_ChannelTime_Cutted_" + detectorName[i]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearSiPM_ChannelTime_Cutted[i]->GetXaxis()->SetTitle("Rear-SiPM Time (ns)");
      H_RearSiPM_ChannelTime_Cutted[i]->GetYaxis()->SetTitle("Rear Channel");
      H_RearSiPM_ChannelTime_Cutted[i]->GetXaxis()->CenterTitle();
      H_RearSiPM_ChannelTime_Cutted[i]->GetYaxis()->CenterTitle();

      H_RearSiPM_MulitplicityTime_Cutted[i] = new TH2D(("RearSiPM_MulitplicityTime_Cutted_" + detectorName[i]).c_str(), ("RearSiPM_MulitplicityTime_Cutted_" + detectorName[i]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta, 20, 0, 20);
      H_RearSiPM_MulitplicityTime_Cutted[i]->GetXaxis()->SetTitle("Rear-SiPM Time (ns)");
      H_RearSiPM_MulitplicityTime_Cutted[i]->GetYaxis()->SetTitle("Multiplicity");
      H_RearSiPM_MulitplicityTime_Cutted[i]->GetXaxis()->CenterTitle();
      H_RearSiPM_MulitplicityTime_Cutted[i]->GetYaxis()->CenterTitle();

      H_SiPM_MultiplicityTime_Cutted[i] = new TH2D(("SiPM_MultiplicityTime_Cutted_" + detectorName[i]).c_str(), ("SiPM_MultiplicityTime_Cutted_" + detectorName[i]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta, 20, 0, 20);
      H_SiPM_MultiplicityTime_Cutted[i]->GetXaxis()->SetTitle("SiPM #Deltat (ns)");
      H_SiPM_MultiplicityTime_Cutted[i]->GetYaxis()->SetTitle("Multiplicity");
      H_SiPM_MultiplicityTime_Cutted[i]->GetXaxis()->CenterTitle();
      H_SiPM_MultiplicityTime_Cutted[i]->GetYaxis()->CenterTitle();

      H_RearSiPM_ChannelTime_Walk_Silicon[i] = new TH2D(("RearSiPM_ChannelTime_Walk_Silicon_" + detectorName[i]).c_str(), ("RearSiPM_ChannelTime_Walk_Silicon_" + detectorName[i]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearSiPM_ChannelTime_Walk_Silicon[i]->GetXaxis()->SetTitle("Rear-SiPM Time (ns)");
      H_RearSiPM_ChannelTime_Walk_Silicon[i]->GetYaxis()->SetTitle("Rear Channel");
      H_RearSiPM_ChannelTime_Walk_Silicon[i]->GetXaxis()->CenterTitle();
      H_RearSiPM_ChannelTime_Walk_Silicon[i]->GetYaxis()->CenterTitle();

      H_RearSiPM_ChannelTime_Walk_Corrected[i] = new TH2D(("RearSiPM_ChannelTime_Walk_Corrected_" + detectorName[i]).c_str(), ("RearSiPM_ChannelTime_Walk_Corrected_" + detectorName[i]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearSiPM_ChannelTime_Walk_Corrected[i]->GetXaxis()->SetTitle("Rear-SiPM Time (ns)");
      H_RearSiPM_ChannelTime_Walk_Corrected[i]->GetYaxis()->SetTitle("Rear Channel");
      H_RearSiPM_ChannelTime_Walk_Corrected[i]->GetXaxis()->CenterTitle();
      H_RearSiPM_ChannelTime_Walk_Corrected[i]->GetYaxis()->CenterTitle();

      H_RearSiPM_ChannelTime_Cleaned[i] = new TH2D(("RearSiPM_ChannelTime_Cleaned_" + detectorName[i]).c_str(), ("RearSiPM_ChannelTime_Cleaned_" + detectorName[i]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearSiPM_ChannelTime_Cleaned[i]->GetXaxis()->SetTitle("Rear-SiPM Time (ns)");
      H_RearSiPM_ChannelTime_Cleaned[i]->GetYaxis()->SetTitle("Rear Channel");
      H_RearSiPM_ChannelTime_Cleaned[i]->GetXaxis()->CenterTitle();
      H_RearSiPM_ChannelTime_Cleaned[i]->GetYaxis()->CenterTitle();

      H_Rear_Channel_IAS_Cleaned[i] = new TH1D(("Rear_Channel_IAS_Cleaned_" + detectorName[i]).c_str(), ("Rear_Channel_IAS_Cleaned_" + detectorName[i]).c_str(), (Rear_IAS[GetDetector(i)].first + Rear_IAS[GetDetector(i)].second) / 10, Rear_IAS[GetDetector(i)].first, Rear_IAS[GetDetector(i)].second);
      H_Rear_Channel_IAS_Cleaned[i]->GetXaxis()->SetTitle("Rear Channel");
      H_Rear_Channel_IAS_Cleaned[i]->GetYaxis()->SetTitle("Counts");
      H_Rear_Channel_IAS_Cleaned[i]->GetXaxis()->CenterTitle();
      H_Rear_Channel_IAS_Cleaned[i]->GetYaxis()->CenterTitle();

      H_Rear_Channel_IAScoinc_Cleaned[i] = new TH1D(("Rear_Channel_IAScoinc_Cleaned_" + detectorName[i]).c_str(), ("Rear_Channel_IAScoinc_Cleaned_" + detectorName[i]).c_str(), (Rear_IAS[GetDetector(i)].first + Rear_IAS[GetDetector(i)].second) / 10, Rear_IAS[GetDetector(i)].first, Rear_IAS[GetDetector(i)].second);
      H_Rear_Channel_IAScoinc_Cleaned[i]->GetXaxis()->SetTitle("Rear Channel");
      H_Rear_Channel_IAScoinc_Cleaned[i]->GetYaxis()->SetTitle("Counts");
      H_Rear_Channel_IAScoinc_Cleaned[i]->GetXaxis()->CenterTitle();
      H_Rear_Channel_IAScoinc_Cleaned[i]->GetYaxis()->CenterTitle();

      H_RearSiPM_Time_Nearest[i] = new TH1D(("RearSiPM_Time_Nearest_" + detectorName[i]).c_str(), ("RearSiPM_Time_Nearest_" + detectorName[i]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
      H_RearSiPM_Time_Nearest[i]->GetXaxis()->SetTitle("Rear-SiPM Time (ns)");
      H_RearSiPM_Time_Nearest[i]->GetYaxis()->SetTitle("Counts");
      H_RearSiPM_Time_Nearest[i]->GetXaxis()->CenterTitle();
      H_RearSiPM_Time_Nearest[i]->GetYaxis()->CenterTitle();

      H_RearSiPM_ChannelTime_dt[i] = new TH2D(("RearSiPM_ChannelTime_dt_" + detectorName[i]).c_str(), ("RearSiPM_ChannelTime_dt_" + detectorName[i]).c_str(), 2e3, 6e3, 10e3, eSiliN / 10, eSiliMin, eSiliMax);
      H_RearSiPM_ChannelTime_dt[i]->GetXaxis()->SetTitle("Rear(dt)-SiPM Time (ns)");
      H_RearSiPM_ChannelTime_dt[i]->GetYaxis()->SetTitle("Rear Channel");
      H_RearSiPM_ChannelTime_dt[i]->GetXaxis()->CenterTitle();
      H_RearSiPM_ChannelTime_dt[i]->GetYaxis()->CenterTitle();

    }

    if (IsDetectorBetaHigh(i))
    {
      dir_SiPMHigh_Channel_Detector[i] = dir_SiPMHigh_Channel->mkdir(detectorName[i].c_str());
      dir_SiPMHigh_Time_Detector[i] = dir_SiPMHigh_Time->mkdir(detectorName[i].c_str());

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

      H_SiPM_ChannelTime_Cutted[i] = new TH2D(("SiPM_ChannelTime_Cutted_" + detectorName[i]).c_str(), ("SiPM_ChannelTime_Cutted_" + detectorName[i]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta, eHighN / 10, eHighMin, eHighMax);
      H_SiPM_ChannelTime_Cutted[i]->GetXaxis()->SetTitle("Rear-SiPM Time (ns)");
      H_SiPM_ChannelTime_Cutted[i]->GetYaxis()->SetTitle("SiPM Channel");
      H_SiPM_ChannelTime_Cutted[i]->GetXaxis()->CenterTitle();
      H_SiPM_ChannelTime_Cutted[i]->GetYaxis()->CenterTitle();

      H_SiPM_ChannelTime_Walk_SiPM[i] = new TH2D(("SiPM_ChannelTime_Walk_SiPM_" + detectorName[i]).c_str(), ("SiPM_ChannelTime_Walk_SiPM_" + detectorName[i]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta, eHighN / 10, eHighMin, eHighMax);
      H_SiPM_ChannelTime_Walk_SiPM[i]->GetXaxis()->SetTitle("SiPM Time (ns)");
      H_SiPM_ChannelTime_Walk_SiPM[i]->GetYaxis()->SetTitle("SiPM Channel");
      H_SiPM_ChannelTime_Walk_SiPM[i]->GetXaxis()->CenterTitle();
      H_SiPM_ChannelTime_Walk_SiPM[i]->GetYaxis()->CenterTitle();

      H_SiPM_ChannelTime_Walk_Corrected[i] = new TH2D(("SiPM_ChannelTime_Walk_Corrected_" + detectorName[i]).c_str(), ("SiPM_ChannelTime_Walk_Corrected_" + detectorName[i]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta, eHighN / 10, eHighMin, eHighMax);
      H_SiPM_ChannelTime_Walk_Corrected[i]->GetXaxis()->SetTitle("Rear-SiPM Time (ns)");
      H_SiPM_ChannelTime_Walk_Corrected[i]->GetYaxis()->SetTitle("SiPM Channel");
      H_SiPM_ChannelTime_Walk_Corrected[i]->GetXaxis()->CenterTitle();
      H_SiPM_ChannelTime_Walk_Corrected[i]->GetYaxis()->CenterTitle();

      H_SiPM_ChannelTime_Cleaned[i] = new TH2D(("SiPM_ChannelTime_Cleaned_" + detectorName[i]).c_str(), ("SiPM_ChannelTime_Cleaned_" + detectorName[i]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta, eHighN / 10, eHighMin, eHighMax);
      H_SiPM_ChannelTime_Cleaned[i]->GetXaxis()->SetTitle("Rear-SiPM Time (ns)");
      H_SiPM_ChannelTime_Cleaned[i]->GetYaxis()->SetTitle("SiPM Channel");
      H_SiPM_ChannelTime_Cleaned[i]->GetXaxis()->CenterTitle();
      H_SiPM_ChannelTime_Cleaned[i]->GetYaxis()->CenterTitle();

      H_SiPMLowHigh_Time_One[i] = new TH1D(("SiPMLowHigh_Time_One_" + detectorName[i]).c_str(), ("SiPMLowHigh_Time_One_" + detectorName[i]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
      H_SiPMLowHigh_Time_One[i]->GetXaxis()->SetTitle("Time between LowHigh (ns)");
      H_SiPMLowHigh_Time_One[i]->GetYaxis()->SetTitle("Counts");
      H_SiPMLowHigh_Time_One[i]->GetXaxis()->CenterTitle();
      H_SiPMLowHigh_Time_One[i]->GetYaxis()->CenterTitle();

      H_SiPMLowHigh_TimeChannel_One[i] = new TH2D(("SiPMLowHigh_TimeChannel_One_" + detectorName[i]).c_str(), ("SiPMLowHigh_TimeChannel_One_" + detectorName[i]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta, eHighN, eHighMin, eHighMax);
      H_SiPMLowHigh_TimeChannel_One[i]->GetXaxis()->SetTitle("Time between LowHigh (ns)");
      H_SiPMLowHigh_TimeChannel_One[i]->GetYaxis()->SetTitle("SiPM High Channel");
      H_SiPMLowHigh_TimeChannel_One[i]->GetXaxis()->CenterTitle();
      H_SiPMLowHigh_TimeChannel_One[i]->GetYaxis()->CenterTitle();

      H_SiPMLowHigh_Time_Multiple[i] = new TH1D(("SiPMLowHigh_Time_Multiple_" + detectorName[i]).c_str(), ("SiPMLowHigh_Time_Multiple_" + detectorName[i]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
      H_SiPMLowHigh_Time_Multiple[i]->GetXaxis()->SetTitle("Time between LowHigh (ns)");
      H_SiPMLowHigh_Time_Multiple[i]->GetYaxis()->SetTitle("Counts");
      H_SiPMLowHigh_Time_Multiple[i]->GetXaxis()->CenterTitle();
      H_SiPMLowHigh_Time_Multiple[i]->GetYaxis()->CenterTitle();

      H_SiPM_Channel_Alone[i] = new TH1D(("SiPM_Channel_Alone_" + detectorName[i]).c_str(), ("SiPM_Channel_Alone_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
      H_SiPM_Channel_Alone[i]->GetXaxis()->SetTitle("SiPM High Channel");
      H_SiPM_Channel_Alone[i]->GetYaxis()->SetTitle("Counts");
      H_SiPM_Channel_Alone[i]->GetXaxis()->CenterTitle();
      H_SiPM_Channel_Alone[i]->GetYaxis()->CenterTitle();

      H_SiPM_Channel_Couple[i] = new TH1D(("SiPM_Channel_Couple_" + detectorName[i]).c_str(), ("SiPM_Channel_Couple_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
      H_SiPM_Channel_Couple[i]->GetXaxis()->SetTitle("SiPM High Channel");
      H_SiPM_Channel_Couple[i]->GetYaxis()->SetTitle("Counts");
      H_SiPM_Channel_Couple[i]->GetXaxis()->CenterTitle();
      H_SiPM_Channel_Couple[i]->GetYaxis()->CenterTitle();

      H_2SiPM_Channel_Couple[i] = new TH2D(("2SiPM_Channel_Couple_" + detectorName[i]).c_str(), ("2SiPM_Channel_Couple_" + detectorName[i]).c_str(), eLowN / 10, eLowMin, eLowMax, eHighN / 10, eHighMin, eHighMax);
      H_2SiPM_Channel_Couple[i]->GetXaxis()->SetTitle("SiPM Low Channel");
      H_2SiPM_Channel_Couple[i]->GetYaxis()->SetTitle("SiPM HIgh Channel");
      H_2SiPM_Channel_Couple[i]->GetXaxis()->CenterTitle();
      H_2SiPM_Channel_Couple[i]->GetYaxis()->CenterTitle();

      H_Channel_Cleaned[i] = new TH1D(("Channel_Cleaned_" + detectorName[i]).c_str(), ("Channel_Cleaned_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
      H_Channel_Cleaned[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_Cleaned[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_Cleaned[i]->GetXaxis()->CenterTitle();
      H_Channel_Cleaned[i]->GetYaxis()->CenterTitle();

      H_SiPM_Pair_Time[i] = new TH1D(("SiPM_Pair_Time_SiPM" + to_string(i)).c_str(), ("SiPM_Pair_Time_SiPM" + to_string(i)).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
      H_SiPM_Pair_Time[i]->GetXaxis()->SetTitle("Time between Low and High (ns)");
      H_SiPM_Pair_Time[i]->GetYaxis()->SetTitle("Counts");
      H_SiPM_Pair_Time[i]->GetXaxis()->CenterTitle();
      H_SiPM_Pair_Time[i]->GetYaxis()->CenterTitle();

      for (int j = i; j < SIGNAL_MAX; j++)
      {
        if (IsDetectorBetaHigh(j))
        {
          H_RearSiPM_Time_Nearest_Nearest[i][j] = new TH1D(("RearSiPM_Time_Nearest_Nearest_" + detectorName[i] + "_" + detectorName[j]).c_str(), ("RearSiPM_Time_Nearest_Nearest_" + detectorName[i] + "_" + detectorName[j]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
          H_RearSiPM_Time_Nearest_Nearest[i][j]->GetXaxis()->SetTitle("Rear-SiPM Time (ns)");
          H_RearSiPM_Time_Nearest_Nearest[i][j]->GetYaxis()->SetTitle("Counts");
          H_RearSiPM_Time_Nearest_Nearest[i][j]->GetXaxis()->CenterTitle();
          H_RearSiPM_Time_Nearest_Nearest[i][j]->GetYaxis()->CenterTitle();

          H_RearSiPM_Channel_Nearset_Nearest[i][j] = new TH2D(("RearSiPM_Channel_Nearset_Nearest_" + detectorName[i] + "_" + detectorName[j]).c_str(), ("RearSiPM_Channel_Nearset_Nearest_" + detectorName[i] + "_" + detectorName[j]).c_str(), eHighN / 10, eHighMin, eHighMax, eHighN / 10, eHighMin, eHighMax);
          H_RearSiPM_Channel_Nearset_Nearest[i][j]->GetXaxis()->SetTitle("SiPM Channel");
          H_RearSiPM_Channel_Nearset_Nearest[i][j]->GetYaxis()->SetTitle("SiPM Channel");
          H_RearSiPM_Channel_Nearset_Nearest[i][j]->GetXaxis()->CenterTitle();
          H_RearSiPM_Channel_Nearset_Nearest[i][j]->GetYaxis()->CenterTitle();

          H_RearSiPM_Channel_Nearest_Group[i][j] = new TH2D(("RearSiPM_Channel_Nearest_Grouped_" + detectorName[i] + "_" + detectorName[j]).c_str(), ("RearSiPM_Channel_Nearest_Grouped_" + detectorName[i] + "_" + detectorName[j]).c_str(), eHighN / 10, eHighMin, eHighMax, eHighN / 10, eHighMin, eHighMax);
          H_RearSiPM_Channel_Nearest_Group[i][j]->GetXaxis()->SetTitle("SiPM Channel");
          H_RearSiPM_Channel_Nearest_Group[i][j]->GetYaxis()->SetTitle("SiPM Channel");
          H_RearSiPM_Channel_Nearest_Group[i][j]->GetXaxis()->CenterTitle();
          H_RearSiPM_Channel_Nearest_Group[i][j]->GetYaxis()->CenterTitle();

          H_2SiPM_Time_Next[i][j] = new TH1D(("2SiPM_Time_Next_" + detectorName[i] + "_" + detectorName[j]).c_str(), ("2SiPM_Time_Next_" + detectorName[i] + "_" + detectorName[j]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
          H_2SiPM_Time_Next[i][j]->GetXaxis()->SetTitle("Time between 2 SiPMs (ns)");
          H_2SiPM_Time_Next[i][j]->GetYaxis()->SetTitle("Counts");
          H_2SiPM_Time_Next[i][j]->GetXaxis()->CenterTitle();
          H_2SiPM_Time_Next[i][j]->GetYaxis()->CenterTitle();
        }
      }

      for (int det = 0; det < SIGNAL_MAX; det++)
      {
        if (IsDetectorSiliBack(det))
        {
          H_RearSiPM_ChannelTime_dt_det[det][i] = new TH2D(("RearSiPM_ChannelTime_dt_" + detectorName[det] + "_" + detectorName[i]).c_str(), ("RearSiPM_ChannelTime_dt_" + detectorName[det] + "_" + detectorName[i]).c_str(), 2e3, 6e3, 10e3, eHighN / 10, eHighMin, eHighMax);
          H_RearSiPM_ChannelTime_dt_det[det][i]->GetXaxis()->SetTitle("Rear(dt)-SiPM Time (ns)");
          H_RearSiPM_ChannelTime_dt_det[det][i]->GetYaxis()->SetTitle("Rear Channel");
          H_RearSiPM_ChannelTime_dt_det[det][i]->GetXaxis()->CenterTitle();
          H_RearSiPM_ChannelTime_dt_det[det][i]->GetYaxis()->CenterTitle();
        }
      }

      
    }

    if (IsDetectorBetaLow(i))
    {
      dir_SiPMLow_Channel_Detector[i] = dir_SiPMLow_Channel->mkdir(detectorName[i].c_str());
      dir_SiPMLow_Time_Detector[i] = dir_SiPMLow_Time->mkdir(detectorName[i].c_str());

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

      H_SiPM_ChannelTime_Cutted[i] = new TH2D(("SiPM_ChannelTime_Cutted_" + detectorName[i]).c_str(), ("SiPM_ChannelTime_Cutted_" + detectorName[i]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta, eLowN / 10, eLowMin, eLowMax);
      H_SiPM_ChannelTime_Cutted[i]->GetXaxis()->SetTitle("Rear-SiPM Time (ns)");
      H_SiPM_ChannelTime_Cutted[i]->GetYaxis()->SetTitle("SiPM Channel");
      H_SiPM_ChannelTime_Cutted[i]->GetXaxis()->CenterTitle();
      H_SiPM_ChannelTime_Cutted[i]->GetYaxis()->CenterTitle();

      H_SiPM_ChannelTime_Walk_SiPM[i] = new TH2D(("SiPM_ChannelTime_Walk_SiPM_" + detectorName[i]).c_str(), ("SiPM_ChannelTime_Walk_SiPM_" + detectorName[i]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta, eLowN / 10, eLowMin, eLowMax);
      H_SiPM_ChannelTime_Walk_SiPM[i]->GetXaxis()->SetTitle("SiPM Time (ns)");
      H_SiPM_ChannelTime_Walk_SiPM[i]->GetYaxis()->SetTitle("SiPM Channel");
      H_SiPM_ChannelTime_Walk_SiPM[i]->GetXaxis()->CenterTitle();
      H_SiPM_ChannelTime_Walk_SiPM[i]->GetYaxis()->CenterTitle();

      H_SiPM_ChannelTime_Walk_Corrected[i] = new TH2D(("SiPM_ChannelTime_Walk_Corrected_" + detectorName[i]).c_str(), ("SiPM_ChannelTime_Walk_Corrected_" + detectorName[i]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta, eLowN / 10, eLowMin, eLowMax);
      H_SiPM_ChannelTime_Walk_Corrected[i]->GetXaxis()->SetTitle("SiPM Time (ns)");
      H_SiPM_ChannelTime_Walk_Corrected[i]->GetYaxis()->SetTitle("SiPM Channel");
      H_SiPM_ChannelTime_Walk_Corrected[i]->GetXaxis()->CenterTitle();
      H_SiPM_ChannelTime_Walk_Corrected[i]->GetYaxis()->CenterTitle();

      H_SiPM_ChannelTime_Cleaned[i] = new TH2D(("SiPM_ChannelTime_Cleaned_" + detectorName[i]).c_str(), ("SiPM_ChannelTime_Cleaned_" + detectorName[i]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta, eLowN / 10, eLowMin, eLowMax);
      H_SiPM_ChannelTime_Cleaned[i]->GetXaxis()->SetTitle("SiPM Time (ns)");
      H_SiPM_ChannelTime_Cleaned[i]->GetYaxis()->SetTitle("SiPM Channel");
      H_SiPM_ChannelTime_Cleaned[i]->GetXaxis()->CenterTitle();
      H_SiPM_ChannelTime_Cleaned[i]->GetYaxis()->CenterTitle();

      H_SiPM_Channel_Alone[i] = new TH1D(("SiPM_Channel_Alone_" + detectorName[i]).c_str(), ("SiPM_Channel_Alone_" + detectorName[i]).c_str(), eLowN, eLowMin, eLowMax);
      H_SiPM_Channel_Alone[i]->GetXaxis()->SetTitle("SiPM Low Channel");
      H_SiPM_Channel_Alone[i]->GetYaxis()->SetTitle("Counts");
      H_SiPM_Channel_Alone[i]->GetXaxis()->CenterTitle();
      H_SiPM_Channel_Alone[i]->GetYaxis()->CenterTitle();

      H_SiPM_Channel_Couple[i] = new TH1D(("SiPM_Channel_Couple_" + detectorName[i]).c_str(), ("SiPM_Channel_Couple_" + detectorName[i]).c_str(), eLowN, eLowMin, eLowMax);
      H_SiPM_Channel_Couple[i]->GetXaxis()->SetTitle("SiPM Low Channel");
      H_SiPM_Channel_Couple[i]->GetYaxis()->SetTitle("Counts");
      H_SiPM_Channel_Couple[i]->GetXaxis()->CenterTitle();
      H_SiPM_Channel_Couple[i]->GetYaxis()->CenterTitle();

      H_Channel_Cleaned[i] = new TH1D(("Channel_Cleaned_" + detectorName[i]).c_str(), ("Channel_Cleaned_" + detectorName[i]).c_str(), eLowN, eLowMin, eLowMax);
      H_Channel_Cleaned[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_Cleaned[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_Cleaned[i]->GetXaxis()->CenterTitle();
      H_Channel_Cleaned[i]->GetYaxis()->CenterTitle();

      for (int j = i; j < SIGNAL_MAX; j++)
      {
        if (IsDetectorBetaLow(j))
        {
          H_2SiPM_Time_Next[i][j] = new TH1D(("2SiPM_Time_Next_" + detectorName[i] + "_" + detectorName[j]).c_str(), ("2SiPM_Time_Next_" + detectorName[i] + "_" + detectorName[j]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
          H_2SiPM_Time_Next[i][j]->GetXaxis()->SetTitle("Time between 2 SiPMs (ns)");
          H_2SiPM_Time_Next[i][j]->GetYaxis()->SetTitle("Counts");
          H_2SiPM_Time_Next[i][j]->GetXaxis()->CenterTitle();
          H_2SiPM_Time_Next[i][j]->GetYaxis()->CenterTitle();
        }
      }
    }
  }

  for (int mul = 0; mul < SIGNAL_MAX; mul++)
  {
    H_RearSiPM_Time_New_Nearest[mul] = new TH1D(("RearSiPM_Time_New_Nearest_M" + to_string(mul)).c_str(), ("RearSiPM_Time_New_Nearest_M" + to_string(mul)).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
    H_RearSiPM_Time_New_Nearest[mul]->GetXaxis()->SetTitle("Rear-SiPM Time (ns)");
    H_RearSiPM_Time_New_Nearest[mul]->GetYaxis()->SetTitle("Counts");
    H_RearSiPM_Time_New_Nearest[mul]->GetXaxis()->CenterTitle();
    H_RearSiPM_Time_New_Nearest[mul]->GetYaxis()->CenterTitle();

    H_RearSiPM_Time_New_Mean[mul] = new TH1D(("RearSiPM_Time_New_Mean_M" + to_string(mul)).c_str(), ("RearSiPM_Time_New_Mean_M" + to_string(mul)).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
    H_RearSiPM_Time_New_Mean[mul]->GetXaxis()->SetTitle("Rear-SiPM Time (ns)");
    H_RearSiPM_Time_New_Mean[mul]->GetYaxis()->SetTitle("Counts");
    H_RearSiPM_Time_New_Mean[mul]->GetXaxis()->CenterTitle();
    H_RearSiPM_Time_New_Mean[mul]->GetYaxis()->CenterTitle();

    H_Time_SubGroups[mul] = new TH1D(("Time_SubGroups_M" + to_string(mul)).c_str(), ("Time_SubGroups_M" + to_string(mul)).c_str(), 300, 0, 600);
    H_Time_SubGroups[mul]->GetXaxis()->SetTitle("SubGroups Time Size (ns)");
    H_Time_SubGroups[mul]->GetYaxis()->SetTitle("Counts");
    H_Time_SubGroups[mul]->GetXaxis()->CenterTitle();
    H_Time_SubGroups[mul]->GetYaxis()->CenterTitle();
  }

  //////////////////////////////
  /////// FROM RAW GROUP ///////
  //////////////////////////////

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

  //////////////////////////////////////////////////////
  /////// FOR CASE C /////
  //////////////////////////////////////////////////////

  H_2Strip_Label_C = new TH2D("2Strip_Label_C", "2Strip_Label_C", 80, 10, 90, 80, 10, 90);
  H_2Strip_Label_C->GetXaxis()->SetTitle("Strip Label");
  H_2Strip_Label_C->GetYaxis()->SetTitle("Strip Label");
  H_2Strip_Label_C->GetXaxis()->CenterTitle();
  H_2Strip_Label_C->GetYaxis()->CenterTitle();

  //////////////////////////////////////////////////////
  /////// FOR CASE D&E : MULTIPLE DETECTOR DETECTION /////
  //////////////////////////////////////////////////////

  /// Label BETWEEN THE 2 STRIPS
  H_2Strip_Label_D = new TH2D("2Strip_Label_D", "2Strip_Label_D", 80, 10, 90, 80, 10, 90);
  H_2Strip_Label_D->GetXaxis()->SetTitle("Strip associated Label");
  H_2Strip_Label_D->GetYaxis()->SetTitle("Label");
  H_2Strip_Label_D->GetXaxis()->CenterTitle();
  H_2Strip_Label_D->GetYaxis()->CenterTitle();

  H_2Strip_Label_E = new TH2D("2Strip_Label_E", "2Strip_Label_E", 80, 10, 90, 80, 10, 90);
  H_2Strip_Label_E->GetXaxis()->SetTitle("Strip associated Label");
  H_2Strip_Label_E->GetYaxis()->SetTitle("Label");
  H_2Strip_Label_E->GetXaxis()->CenterTitle();
  H_2Strip_Label_E->GetYaxis()->CenterTitle();

  /////////////////////////////////////////////////////
  /////// FOR CASE G : MULTIPLE DETECTOR DETECTION /////
  /////////////////////////////////////////////////////

  H_RearStrip_Channel_G = new TH2D("RearStrip_Channel_G", "RearStrip_Channel_G", eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
  H_RearStrip_Channel_G->GetXaxis()->SetTitle("Rear Channel");
  H_RearStrip_Channel_G->GetYaxis()->SetTitle("Strip Channel");
  H_RearStrip_Channel_G->GetXaxis()->CenterTitle();
  H_RearStrip_Channel_G->GetYaxis()->CenterTitle();

  H_RearStrip_Time_G = new TH1D("RearStrip_Time_G", "RearStrip_Time_G", winGroupN_Sili, winGroupMin_Sili, winGroupMax_Sili);
  H_RearStrip_Time_G->GetXaxis()->SetTitle("Rear-Strip Time (ns)");
  H_RearStrip_Time_G->GetYaxis()->SetTitle("Counts");
  H_RearStrip_Time_G->GetXaxis()->CenterTitle();
  H_RearStrip_Time_G->GetYaxis()->CenterTitle();

  H_RearStrip_Label_G = new TH2D("RearStrip_Label_G", "RearStrip_Label_G", 80, 10, 90, 80, 10, 90);
  H_RearStrip_Label_G->GetXaxis()->SetTitle("Strip Label");
  H_RearStrip_Label_G->GetYaxis()->SetTitle("Rear Label");
  H_RearStrip_Label_G->GetXaxis()->CenterTitle();
  H_RearStrip_Label_G->GetYaxis()->CenterTitle();

  ///////////////////////////////////////////
  /////// checking faster configuration /////
  ///////////////////////////////////////////

  H_2Strip_Label_Check = new TH2D("2Strip_Label_Check", "2Strip_Label_Check", 80, 10, 90, 80, 10, 90);
  H_2Strip_Label_Check->GetXaxis()->SetTitle("Strip Label");
  H_2Strip_Label_Check->GetYaxis()->SetTitle("Strip Label");
  H_2Strip_Label_Check->GetXaxis()->CenterTitle();
  H_2Strip_Label_Check->GetYaxis()->CenterTitle();

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

  //////////////////////////////////////////
  /////// FROM CLEANED /////////////////////
  //////////////////////////////////////////

  H_SiPMLow_Mulitplicity_Cleaned = new TH1D(("SiPMLow_Mulitplicity_Cleaned"), ("SiPMLow_Mulitplicity_Cleaned"), 20, 0, 20);
  H_SiPMLow_Mulitplicity_Cleaned->GetXaxis()->SetTitle("Multiplicity");
  H_SiPMLow_Mulitplicity_Cleaned->GetYaxis()->SetTitle("Counts");
  H_SiPMLow_Mulitplicity_Cleaned->GetXaxis()->CenterTitle();
  H_SiPMLow_Mulitplicity_Cleaned->GetYaxis()->CenterTitle();

  H_SiPMHigh_Mulitplicity_Cleaned = new TH1D(("SiPMHigh_Mulitplicity_Cleaned"), ("SiPMHigh_Mulitplicity_Cleaned"), 20, 0, 20);
  H_SiPMHigh_Mulitplicity_Cleaned->GetXaxis()->SetTitle("Multiplicity");
  H_SiPMHigh_Mulitplicity_Cleaned->GetYaxis()->SetTitle("Counts");
  H_SiPMHigh_Mulitplicity_Cleaned->GetXaxis()->CenterTitle();
  H_SiPMHigh_Mulitplicity_Cleaned->GetYaxis()->CenterTitle();

  //////////////////////////////////////////
  /////// FROM SIPM GROUPS /////////////////////
  //////////////////////////////////////////

  H_Multiplicity_Groups = new TH1D(("Multiplicity_Groups"), ("Multiplicity_Groups"), 20, 0, 20);
  H_Multiplicity_Groups->GetXaxis()->SetTitle("# of SubGroups");
  H_Multiplicity_Groups->GetYaxis()->SetTitle("Counts");
  H_Multiplicity_Groups->GetXaxis()->CenterTitle();
  H_Multiplicity_Groups->GetYaxis()->CenterTitle();

  H_Multiplicity_SubGroups = new TH1D(("Multiplicity_SubGroups"), ("Multiplicity_SubGroups"), 20, 0, 20);
  H_Multiplicity_SubGroups->GetXaxis()->SetTitle("Multiplicity");
  H_Multiplicity_SubGroups->GetYaxis()->SetTitle("Counts");
  H_Multiplicity_SubGroups->GetXaxis()->CenterTitle();
  H_Multiplicity_SubGroups->GetYaxis()->CenterTitle();

  H_Multiplicity_Groups_IAS = new TH1D(("Multiplicity_Groups_IAS"), ("Multiplicity_Groups_IAS"), 20, 0, 20);
  H_Multiplicity_Groups_IAS->GetXaxis()->SetTitle("# of SubGroups");
  H_Multiplicity_Groups_IAS->GetYaxis()->SetTitle("Counts");
  H_Multiplicity_Groups_IAS->GetXaxis()->CenterTitle();
  H_Multiplicity_Groups_IAS->GetYaxis()->CenterTitle();

  H_Multiplicity_SubGroups_IAS = new TH1D(("Multiplicity_SubGroups_IAS"), ("Multiplicity_SubGroups_IAS"), 20, 0, 20);
  H_Multiplicity_SubGroups_IAS->GetXaxis()->SetTitle("Multiplicity");
  H_Multiplicity_SubGroups_IAS->GetYaxis()->SetTitle("Counts");
  H_Multiplicity_SubGroups_IAS->GetXaxis()->CenterTitle();
  H_Multiplicity_SubGroups_IAS->GetYaxis()->CenterTitle();

  //////////////////////////////////////////////////////

  return 0;
}

inline int WriteHistograms_Grouped()
{
  GROUPED_File->cd();
  for (size_t i = 0; i < detectorNum; ++i)
  {
    ProgressCounter(i, detectorNum, " <INFO>  Writing Histograms");
    GROUPED_File->cd();
    Verbose("Writing " + detectorName[i], VERBOSE, 1);
    if (IsDetectorSiliStrip(i))
    {
      ///////////// CHANNEL //////////////////////////////////////////////////////
      Verbose("Channel", VERBOSE, 1);
      // WRTTING //
      dir_Group_Strip->cd();
      dir_Strip_Channel_Detector[i]->cd();
      H_Channel_RAW[i]->Write();
      H_Channel_PileUp[i]->Write();
      H_Channel_A[i]->Write();
      H_Channel_B[i]->Write();
      H_Channel_C[i]->Write();
      H_Channel_D[i]->Write();
      H_Channel_E[i]->Write();
      H_Channel_G[i]->Write();
      H_Channel_H[i]->Write();
      H_RearStrip_Channel_RAW[i]->Write();
      H_RearStrip_Channel_A[i]->Write();
      H_RearStrip_Channel_B[i]->Write();
      H_2Strip_Channel_B[i]->Write();
      H_RearStrip_Channel_C[i]->Write();
      H_2Strip_Channel_C[i]->Write();
      H_2Strip_Channel_D[i]->Write();
      H_2Strip_Channel_E[i]->Write();
      /////////////

      dir_CASES_Channel->cd();
      TCanvas *C_RAW_CASE_A = new TCanvas(("C_RAW_CASE_A_" + detectorName[i]).c_str(), ("C_RAW_CASE_A_" + detectorName[i]).c_str(), 800, 400);
      C_RAW_CASE_A->Divide(1, 2);
      C_RAW_CASE_A->cd(1);
      H_Channel_RAW[i]->SetLineColor(kBlue);
      H_Channel_RAW[i]->Draw();
      H_Channel_A[i]->SetLineColor(kRed);
      H_Channel_A[i]->Draw("SAME");
      C_RAW_CASE_A->cd(2);
      TH1D *H_Channel_RAW_A_SUB = (TH1D *)H_Channel_RAW[i]->Clone(("H_Channel_RAW_A_SUB" + detectorName[i]).c_str());
      H_Channel_RAW_A_SUB->Add(H_Channel_A[i], -1);
      H_Channel_RAW_A_SUB->Draw();
      C_RAW_CASE_A->Write();

      TCanvas *C_RAW_PileUp = new TCanvas(("C_RAW_PileUp_" + detectorName[i]).c_str(), ("C_RAW_PileUp_" + detectorName[i]).c_str(), 800, 400);
      C_RAW_PileUp->Divide(1, 2);
      C_RAW_PileUp->cd(1);
      H_Channel_RAW[i]->SetLineColor(kBlue);
      H_Channel_RAW[i]->Draw();
      H_Channel_PileUp[i]->SetLineColor(kRed);
      H_Channel_PileUp[i]->Draw("SAME");
      C_RAW_PileUp->cd(2);
      TH1D *H_Channel_RAW_PileUp_SUB = (TH1D *)H_Channel_RAW[i]->Clone(("H_Channel_RAW_PileUp_SUB" + detectorName[i]).c_str());
      H_Channel_RAW_PileUp_SUB->Add(H_Channel_PileUp[i], -1);
      H_Channel_RAW_PileUp_SUB->Draw();
      C_RAW_PileUp->Write();

      ///////////// TIME //////////////////////////////////////////////////////
      Verbose("Time", VERBOSE, 1);
      // WRTTING //
      dir_Strip_Time_Detector[i]->cd();
      H_RearStrip_Time_A[i]->Write();
      H_RearStrip_Time_B[i]->Write();
      H_RearStrip_Time_C[i]->Write();
      H_2Strip_Time_C[i]->Write();
      H_2Strip_Time_D[i]->Write();
      H_2Strip_Time_E[i]->Write();
      /////////////

      ////////////////////////////////////////////////////////////////////////////
      ///////////// CUTTING //////////////////////////////////////////////////////
      Verbose("Cutting", VERBOSE, 1);
      // Writting //
      dir_Strip_Channel_Detector[i]->cd();
      H_Strip_Channel_DiffRear_A[i]->Write();
      H_Strip_Channel_DiffRearvsStrip_A[i]->Write();

      dir_Cutting_Fitting->cd();

      double x_center = H_Strip_Channel_DiffRear_A[i]->GetBinCenter(H_Strip_Channel_DiffRear_A[i]->GetMaximumBin());
      int counter = 0;

      G_RearStrip_Spread[i] = new TGraphErrors();
      G_RearStrip_Spread[i]->SetTitle(("Spread_" + detectorName[i]).c_str());
      G_RearStrip_Spread[i]->GetXaxis()->SetTitle("Strip Channel");
      G_RearStrip_Spread[i]->GetYaxis()->SetTitle("Spread (#sigma of #DeltaE/E)");
      G_RearStrip_Spread[i]->GetXaxis()->CenterTitle();
      G_RearStrip_Spread[i]->GetYaxis()->CenterTitle();

      // graph mean
      TGraphErrors *graph_mean = new TGraphErrors();
      graph_mean->SetTitle("Mean");
      graph_mean->GetXaxis()->SetTitle("Channel");
      graph_mean->GetYaxis()->SetTitle("Mean of #DeltaE/E");
      graph_mean->GetXaxis()->CenterTitle();
      graph_mean->GetYaxis()->CenterTitle();

      int step = 10;

      TH1D *H_Strip_Channel_DiffRear_A_YProjection = H_Strip_Channel_DiffRearvsStrip_A[i]->ProjectionY(("H_Strip_Channel_DiffRearvsStrip_A_YProjection_" + detectorName[i]).c_str());
      for (int bin = 0; bin < H_Strip_Channel_DiffRear_A[i]->GetNbinsX(); bin += step)
      {
        TH1D *H_Strip_Channel_DiffRear_A_XProjection = H_Strip_Channel_DiffRearvsStrip_A[i]->ProjectionX(("H_Strip_Channel_DiffRearvsStrip_A_XProjection_" + detectorName[i] + "_" + to_string(bin)).c_str(), bin, bin + step);
        H_Strip_Channel_DiffRear_A_XProjection->GetXaxis()->SetRangeUser(x_center - 0.05, x_center + 0.05);
        TF1 *gauss = new TF1(("gauss_" + detectorName[i] + "_" + to_string(bin)).c_str(), "gaus", x_center - 0.05, x_center + 0.05);
        gauss->SetParameter(1, x_center);
        gauss->SetParameter(2, 0.005);
        if (H_Strip_Channel_DiffRear_A_XProjection->Integral() < 20)
          continue; /// avoid empty data for fit

        H_Strip_Channel_DiffRear_A_XProjection->Fit(gauss, "QRN");

        if (gauss->GetParError(2) / gauss->GetParameter(2) > 0.2)
          continue; /// avoid weird data points due to low statistics

        G_RearStrip_Spread[i]->SetPoint(counter, H_Strip_Channel_DiffRear_A_YProjection->GetBinCenter(bin + step / 2), gauss->GetParameter(2));
        G_RearStrip_Spread[i]->SetPointError(counter, H_Strip_Channel_DiffRear_A_YProjection->GetBinCenter(bin + step / 2) - H_Strip_Channel_DiffRear_A_YProjection->GetBinCenter(bin), gauss->GetParError(2));

        graph_mean->SetPoint(counter, H_Strip_Channel_DiffRear_A_YProjection->GetBinCenter(bin + step / 2), gauss->GetParameter(1));
        graph_mean->SetPointError(counter, H_Strip_Channel_DiffRear_A_YProjection->GetBinCenter(bin + step / 2) - H_Strip_Channel_DiffRear_A_YProjection->GetBinCenter(bin), gauss->GetParError(1));
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
      G_RearStrip_Spread[i]->SetMarkerStyle(20);
      G_RearStrip_Spread[i]->SetMarkerSize(0.5);
      G_RearStrip_Spread[i]->Draw("AP");
      SpreadAcceptance[i] = new TF1("spread_fit", "[0]/x");
      SpreadAcceptance[i]->SetParameter(0, 100);
      int status = 0;
      if (G_RearStrip_Spread[i]->GetN() > 5)
      {
        if (SpreadAcceptance[i] == nullptr)
          status = G_RearStrip_Spread[i]->Fit(SpreadAcceptance[i], "QR");

        if (status != 0)
          Warning("Fit failed for " + detectorName[i] + " during Rear/Strip spread correction.");
      }

      double yfullrange = G_RearStrip_Spread[i]->GetYaxis()->GetXmax() - G_RearStrip_Spread[i]->GetYaxis()->GetXmin();
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
      delete H_Channel_RAW[i];
      delete H_Channel_A[i];
      // delete H_Channel_B[i];
      // delete H_Channel_C[i];
      // delete H_Channel_D[i];
      // delete H_Channel_E[i];
      // delete H_Channel_G[i];
      // delete C_RAW_CASE_A;
    }

    if (IsDetectorSiliBack(i))
    {
      // GROUP
      dir_Group_Rear->cd();
      // CHANNEL
      dir_Rear_Channel_Detector[i]->cd();
      H_Channel_RAW[i]->Write();
      H_Channel_PileUp[i]->Write();
      H_Channel_A[i]->Write();
      H_Channel_B[i]->Write();
      H_Channel_C[i]->Write();
      H_Channel_D[i]->Write();
      H_Channel_E[i]->Write();
      H_Channel_F[i]->Write();
      H_Channel_G[i]->Write();
      H_Channel_H[i]->Write();
      H_RearStrip_Channel_RAW[i]->Write();
      H_RearStrip_Channel_A[i]->Write();
      H_RearStrip_Channel_B[i]->Write();
      H_RearStrip_Channel_C[i]->Write();
      H_Rear2Strip_Channel_B[i]->Write();

      dir_CASES_Channel->cd();
      TCanvas *C_RAW_CASE_A = new TCanvas(("C_RAW_CASE_A_" + detectorName[i]).c_str(), ("C_RAW_CASE_A_" + detectorName[i]).c_str(), 800, 400);
      C_RAW_CASE_A->Divide(1, 2);
      C_RAW_CASE_A->cd(1);
      H_Channel_RAW[i]->SetLineColor(kBlue);
      H_Channel_RAW[i]->Draw();
      H_Channel_A[i]->SetLineColor(kRed);
      H_Channel_A[i]->Draw("SAME");
      C_RAW_CASE_A->cd(2);
      TH1D *H_Channel_RAW_A_SUB = (TH1D *)H_Channel_RAW[i]->Clone(("H_Channel_RAW_A_SUB" + detectorName[i]).c_str());
      H_Channel_RAW_A_SUB->Add(H_Channel_A[i], -1);
      H_Channel_RAW_A_SUB->Draw();
      C_RAW_CASE_A->Write();

      TCanvas *C_RAW_PileUp = new TCanvas(("C_RAW_PileUp_" + detectorName[i]).c_str(), ("C_RAW_PileUp_" + detectorName[i]).c_str(), 800, 400);
      C_RAW_PileUp->Divide(1, 2);
      C_RAW_PileUp->cd(1);
      H_Channel_RAW[i]->SetLineColor(kBlue);
      H_Channel_RAW[i]->Draw();
      H_Channel_PileUp[i]->SetLineColor(kRed);
      H_Channel_PileUp[i]->Draw("SAME");
      C_RAW_PileUp->cd(2);
      TH1D *H_Channel_RAW_PileUp_SUB = (TH1D *)H_Channel_RAW[i]->Clone(("H_Channel_RAW_PileUp_SUB" + detectorName[i]).c_str());
      H_Channel_RAW_PileUp_SUB->Add(H_Channel_PileUp[i], -1);
      H_Channel_RAW_PileUp_SUB->Draw();
      C_RAW_PileUp->Write();

      // TIME
      // Writing
      dir_Rear_Time_Detector[i]->cd();
      H_RearStrip_Time_A[i]->Write();
      H_RearStrip_Time_B[i]->Write();
      H_RearStrip_Time_C[i]->Write();
      H_2Strip_Time_D[i]->Write();
      H_2Strip_Time_E[i]->Write();

      /////// deleting
      delete C_RAW_CASE_A;
      // delete H_Channel_RAW[i];
    }

    if (IsDetectorBetaHigh(i))
    {
      dir_SiPMHigh_Channel_Detector[i]->cd();
      H_Channel_RAW[i]->Write();
    }

    if (IsDetectorBetaLow(i))
    {
      dir_SiPMLow_Channel_Detector[i]->cd();
      H_Channel_RAW[i]->Write();
    }
  }

  /// STRIP/////////////////////////////////////////////
  dir_Strip_Multiplicities->cd();
  H_Strip_Mulitplicity_RAW->Write();
  H_2Strip_Label_C->Write();
  H_2Strip_Label_D->Write();
  H_2Strip_Label_E->Write();
  H_RearStrip_Label_G->Write();
  H_2Strip_Label_Check->Write();

  dir_Strip_Time->cd();
  H_RearStrip_Time_G->Write();

  // deleting
  delete H_Strip_Mulitplicity_RAW;

  /// REAR/////////////////////////////////////////////
  dir_Rear_Multiplicity->cd();
  H_Rear_Mulitplicity_RAW->Write();

  // deleting
  delete H_Rear_Mulitplicity_RAW;

  /// SiPM Low/////////////////////////////////////////////
  dir_SiPMLow_Multiplicity->cd();
  H_SiPMLow_Mulitplicity_RAW->Write();

  /// SiPM High/////////////////////////////////////////////
  dir_SiPMHigh_Multiplicity->cd();
  H_SiPMHigh_Mulitplicity_RAW->Write();

  GROUPED_File->cd();
  return 0;
}

inline void WriteHistograms_Cutted()
{
  GROUPED_File->cd();

  for (int i = 0; i < detectorNum; i++)
  {
    ProgressCounter(i, detectorNum, " <INFO>  Writing Histograms");
    if (IsDetectorSiliStrip(i))
    {
      dir_Strip_Channel_Detector[i]->cd();
      H_Channel_Cutted[i]->Write();
      H_RearStrip_Channel_Cutted[i]->Write();
      H_Strip_Channel_DiffRear_Cutted[i]->Write();
      H_Strip_Channel_DiffRearvsStrip_Cutted[i]->Write();
      dir_Strip_Time_Detector[i]->cd();
      H_RearStrip_Time_Cutted[i]->Write();
      H_RearStrip_ChannelTime_Cutted[i]->Write();

      // UPDATE CANVAS
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

      dir_Cutting_Fitting->cd();
      C_Strip_Channel_DiffRear[i]->Write();

      // general fir result
      C_Strip_Cutting_Fits->cd(ConvertBase5ToBase10(i) - 5);
      G_RearStrip_Spread[i]->Draw("AP");
      SpreadAcceptance[i]->Draw("SAME");

      // TIME CANVAS
      dir_Strip_Time->cd();
      TCanvas *C_Strip_Time = new TCanvas(("C_Strip_Time_A_CUTTED_" + detectorName[i]).c_str(), ("C_Strip_Time_A_CUTTED_" + detectorName[i]).c_str(), 800, 400);
      H_RearStrip_Time_Cutted[i]->SetLineColor(kRed);
      H_RearStrip_Time_Cutted[i]->Draw("HIST");
      H_RearStrip_Time_A[i]->SetLineColor(kBlue);
      H_RearStrip_Time_A[i]->Draw("SAME");
      TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);
      legend->AddEntry(H_RearStrip_Time_A[i], "RAW", "l");
      legend->AddEntry(H_RearStrip_Time_Cutted[i], "CUTTED", "l");
      legend->Draw();
      C_Strip_Time->Write();

      // Matching Channel Strip-Rear
      TF1 *pol1 = new TF1("pol1", "pol1", eSiliMin, eSiliMax);
      P_RearStrip_Channel_Cutted[i]->Fit(pol1, "QN");
      RearStrip_Matching[i] = make_pair(pol1->GetParameter(1), pol1->GetParameter(0));
      TCanvas *C_RearStrip_Channel_Cutted = new TCanvas(("C_RearStrip_Channel_Matching_" + detectorName[i]).c_str(), ("C_RearStrip_Channel_Matching_" + detectorName[i]).c_str(), 800, 400);
      P_RearStrip_Channel_Cutted[i]->Draw("AP");
      pol1->SetLineColor(kRed);
      pol1->Draw("SAME");
      TPaveText *pt = new TPaveText(0.1, 0.7, 0.48, 0.9, "NDC");
      pt->AddText(("f(x) = " + to_string(pol1->GetParameter(1)) + "x + " + to_string(pol1->GetParameter(0))).c_str());
      pt->AddText(("#chi^{2}/NDF = " + to_string(pol1->GetChisquare() / pol1->GetNDF())).c_str());
      pt->Draw("SAME");
      C_RearStrip_Channel_Cutted->Write();
    }

    if (IsDetectorSiliBack(i))
    {
      dir_Rear_Channel_Detector[i]->cd();
      H_Channel_Cutted[i]->Write();
      dir_Rear_Time_Detector[i]->cd();
      H_RearStrip_Time_Cutted[i]->Write();
      H_RearSiPM_Time_Cutted[i]->Write();
      H_RearSiPM_ChannelTime_Cutted[i]->Write();
      H_RearSiPM_MulitplicityTime_Cutted[i]->Write();
      H_SiPM_MultiplicityTime_Cutted[i]->Write();

      /////// CLEANING SIPMs ON REAR ///////
      dir_Walk_Silicon->cd();

      int counter = 0;
      // graph mean
      G_Rear_Mean_Walk_Silicon[i] = new TGraphErrors();
      G_Rear_Mean_Walk_Silicon[i]->SetTitle(("Mean_#Deltat_" + detectorName[i]).c_str());
      G_Rear_Mean_Walk_Silicon[i]->GetXaxis()->SetTitle("Rear Channel");
      G_Rear_Mean_Walk_Silicon[i]->GetYaxis()->SetTitle("Mean Rear-SiPM Time (ns)");
      G_Rear_Mean_Walk_Silicon[i]->GetXaxis()->CenterTitle();
      G_Rear_Mean_Walk_Silicon[i]->GetYaxis()->CenterTitle();

      int step = 10;
      TH1D *H_RearSiPM_ChannelTime_Cutted_YProjection = H_RearSiPM_ChannelTime_Cutted[i]->ProjectionY(("H_RearSiPM_ChannelTime_Cutted_YProjection_" + detectorName[i]).c_str());
      for (int bin = 0; bin < H_RearSiPM_ChannelTime_Cutted[i]->GetNbinsY(); bin += step)
      {
        TH1D *H_RearSiPM_ChannelTime_Cutted_XProjection = H_RearSiPM_ChannelTime_Cutted[i]->ProjectionX(("H_RearSiPM_ChannelTime_Cutted_XProjection_" + detectorName[i] + "_" + to_string(bin)).c_str(), bin, bin + step);
        double x_center = H_RearSiPM_ChannelTime_Cutted_XProjection->GetBinCenter(H_RearSiPM_ChannelTime_Cutted_XProjection->GetMaximumBin());
        H_RearSiPM_ChannelTime_Cutted_XProjection->GetXaxis()->SetRangeUser(x_center - 100, x_center + 100);
        TF1 *gauss = new TF1(("skewgauss_" + detectorName[i] + "_" + to_string(bin)).c_str(), skewedgauss, x_center - 100, x_center + 100, 5);
        gauss->SetParameter(3, H_RearSiPM_ChannelTime_Cutted_XProjection->GetBinContent(H_RearSiPM_ChannelTime_Cutted_XProjection->GetMaximumBin()));
        gauss->SetParameter(0, x_center);
        gauss->SetParameter(1, 10);
        gauss->SetParLimits(2, -100, 0);
        if (H_RearSiPM_ChannelTime_Cutted_XProjection->Integral() < 100)
          continue; /// avoid empty data for fit

        H_RearSiPM_ChannelTime_Cutted_XProjection->Fit(gauss, "QRN");

        if (gauss->GetParError(0) / gauss->GetParameter(0) > 0.2)
          continue; /// avoid weird data points due to low statistics

        if ((x_center - gauss->GetMaximumX()) / gauss->GetMaximumX() > 0.2)
          continue; /// avoid weird data points due fit failling

        G_Rear_Mean_Walk_Silicon[i]->SetPoint(counter, H_RearSiPM_ChannelTime_Cutted_YProjection->GetBinCenter(bin + step / 2), abs(gauss->GetMaximumX()));
        G_Rear_Mean_Walk_Silicon[i]->SetPointError(counter, H_RearSiPM_ChannelTime_Cutted_YProjection->GetBinCenter(bin + step / 2) - H_RearSiPM_ChannelTime_Cutted_YProjection->GetBinCenter(bin), gauss->GetParError(0));
        counter++;

        TCanvas *c = new TCanvas(("C_RearSiPM_ChannelTime_Cutted_" + detectorName[i] + "_" + to_string(bin)).c_str(), ("C_RearSiPM_ChannelTime_Cutted_" + detectorName[i] + "_" + to_string(bin)).c_str(), 800, 400);
        H_RearSiPM_ChannelTime_Cutted_XProjection->Draw("HIST");
        gauss->Draw("SAME");
        // c->Write();
      }

      C_RearSiPM_ChannelTime_Walk_Silicon[i] = new TCanvas(("C_RearSiPM_ChannelTime_Walk_Silicon_" + detectorName[i]).c_str(), ("C_RearSiPM_ChannelTime_Walk_Silicon_" + detectorName[i]).c_str(), 800, 400);

      ////3rd
      C_RearSiPM_ChannelTime_Walk_Silicon[i]->cd();
      TPad *P_DiffRearStripvsStrip = new TPad("MeanRearSiPMvsT", "MeanRearSiPMvsT", 0.05, 0.05, 0.5, 0.5);
      P_DiffRearStripvsStrip->Draw();
      P_DiffRearStripvsStrip->cd();
      G_Rear_Mean_Walk_Silicon[i]->SetMarkerStyle(20);
      G_Rear_Mean_Walk_Silicon[i]->SetMarkerSize(0.5);
      G_Rear_Mean_Walk_Silicon[i]->Draw("AP");
      MeanAcceptance_Walk_Silicon[i] = new TF1("inverse", "[0] + [1]/(x+[2])");
      MeanAcceptance_Walk_Silicon[i]->SetParameters(150, -1e6, 2000);
      MeanAcceptance_Walk_Silicon[i]->SetParLimits(0, 20, 150);
      MeanAcceptance_Walk_Silicon[i]->SetParLimits(1, -1e7, -1e5);
      MeanAcceptance_Walk_Silicon[i]->SetParLimits(2, -150000, 10000);

      int status = 1;
      if (G_Rear_Mean_Walk_Silicon[i]->GetN() != 0)
      {
        status = G_Rear_Mean_Walk_Silicon[i]->Fit(MeanAcceptance_Walk_Silicon[i], "QW");
        if (status != 0)
          Warning("Fit failed for " + detectorName[i] + " during walk Silicon correction.");
      }
      // cout << "     " << MeanAcceptance_Walk_Silicon[i]->GetChisquare () << endl;

      ////1st
      C_RearSiPM_ChannelTime_Walk_Silicon[i]->cd();
      TPad *P_RearStrip = new TPad("RearSiPM_ChannelTime", "RearSiPM_ChannelTime", 0.05, 0.5, 0.5, 1.0);
      P_RearStrip->Draw();
      P_RearStrip->cd();
      P_RearStrip->SetLogz();
      H_RearSiPM_ChannelTime_Cutted[i]->Draw("COLZ");

      maximum[i] = new TGraphErrors();
      for (int bin = 0; bin < H_RearSiPM_ChannelTime_Cutted[i]->GetNbinsY(); bin++)
      {

        double x_mean = -MeanAcceptance_Walk_Silicon[i]->Eval(H_RearSiPM_ChannelTime_Cutted[i]->GetYaxis()->GetBinCenter(bin));
        maximum[i]->SetPoint(bin, x_mean, H_RearSiPM_ChannelTime_Cutted[i]->GetYaxis()->GetBinCenter(bin));
      }
      maximum[i]->SetLineColor(kRed);
      maximum[i]->Draw("SAME");
    }

    if (IsDetectorBetaHigh(i))
    {
      /// SIPM High
      dir_SiPMHigh_Channel_Detector[i]->cd();
      H_Channel_Cutted[i]->Write();
      dir_SiPMHigh_Time_Detector[i]->cd();
      H_SiPM_ChannelTime_Cutted[i]->Write();

      dir_Cutting_SiPM->cd();
      TCanvas *C_SiPMHigh_Channel = new TCanvas(("C_SiPMHigh_Channel_RAW_CUTTED_" + detectorName[i]).c_str(), ("C_SiPMHigh_Channel_" + detectorName[i]).c_str(), 800, 400);
      H_Channel_RAW[i]->Draw();
      H_Channel_Cutted[i]->SetLineColor(kRed);
      H_Channel_Cutted[i]->Draw("SAME");
      C_SiPMHigh_Channel->Write();
    }

    if (IsDetectorBetaLow(i))
    {
      /// SIPM Low
      dir_SiPMLow_Channel_Detector[i]->cd();
      H_Channel_Cutted[i]->Write();
      dir_SiPMLow_Time_Detector[i]->cd();
      H_SiPM_ChannelTime_Cutted[i]->Write();

      dir_Cutting_SiPM->cd();
      TCanvas *C_SiPMLow_Channel = new TCanvas(("C_SiPMLow_Channel_RAW_CUTTED_" + detectorName[i]).c_str(), ("C_SiPMLow_Channel_" + detectorName[i]).c_str(), 800, 400);
      H_Channel_RAW[i]->Draw();
      H_Channel_Cutted[i]->SetLineColor(kRed);
      H_Channel_Cutted[i]->Draw("SAME");
      C_SiPMLow_Channel->Write();
    }
  }

  dir_SiPMHigh_Multiplicity->cd();
  H_SiPMHigh_Mulitplicity_Cutted->Write();
  dir_Cutting_SiPM->cd();
  TCanvas *C_SiPMHigh_Multiplicities = new TCanvas("SiPMHigh_Multiplicities_RAW_CUTTED", "SiPMHigh_Multiplicities_RAW_CUTTED", 800, 400);
  H_SiPMHigh_Mulitplicity_RAW->Draw();
  H_SiPMHigh_Mulitplicity_Cutted->SetLineColor(kRed);
  H_SiPMHigh_Mulitplicity_Cutted->Draw("SAME");
  C_SiPMHigh_Multiplicities->Write();

  dir_SiPMLow_Multiplicity->cd();
  H_SiPMLow_Mulitplicity_Cutted->Write();
  dir_Cutting_SiPM->cd();
  H_SiPMLow_Mulitplicity_RAW->Draw();
  H_SiPMLow_Mulitplicity_Cutted->SetLineColor(kRed);
  H_SiPMLow_Mulitplicity_Cutted->Draw("SAME");
  // C_SiPMLow_Multiplicities->Write();
}

inline void WriteHistograms_SiliconWalk()
{
  for (int i = 0; i < detectorNum; i++)
  {
    ProgressCounter(i, detectorNum, " <INFO>  Writing Histograms");
    if (IsDetectorSiliBack(i))
    {
      dir_Walk_Silicon->cd();

      C_RearSiPM_ChannelTime_Walk_Silicon[i]->cd();
      TPad *P_RearSiPM_ChannelTime_Walk_Silicon = new TPad("RearSiPM_ChannelTime_Walk_Silicon", "RearSiPM_ChannelTime_Walk_Silicon", 0.55, 0.50, 1.0, 1.0);
      P_RearSiPM_ChannelTime_Walk_Silicon->Draw();
      P_RearSiPM_ChannelTime_Walk_Silicon->cd();
      P_RearSiPM_ChannelTime_Walk_Silicon->SetLogz();
      H_RearSiPM_ChannelTime_Walk_Silicon[i]->Draw("COLZ");

      C_RearSiPM_ChannelTime_Walk_Silicon[i]->cd();
      TPad *P_DiffRearStrip = new TPad("RearSiPM_ChannelTime_Walk_Silicon_proj", "RearSiPM_ChannelTime_Walk_Silicon_proj", 0.55, 0.05, 1.0, 0.5);
      P_DiffRearStrip->Draw();
      P_DiffRearStrip->cd();
      P_DiffRearStrip->SetLogy();
      TH1D *proj = (TH1D *)H_RearSiPM_ChannelTime_Walk_Silicon[i]->ProjectionX(("H_RearSiPM_ChannelTime_Walk_Silicon_projx" + detectorName[i]).c_str());
      proj->Draw("HIST");

      C_RearSiPM_ChannelTime_Walk_Silicon[i]->Write();

      C_Rear_Walk_Silicon_Fits->cd(GetDetector(i));
      G_Rear_Mean_Walk_Silicon[i]->Draw("AP");
      MeanAcceptance_Walk_Silicon[i]->Draw("SAME");
    }

    if (IsDetectorBeta(i))
    {

      int step = 0;
      int mini = 0;
      int maxi = 0;
      int suppar2 = 0;
      if (IsDetectorBetaHigh(i))
      {
        step = 2;
        mini = 0;
        maxi = 2e6;
        suppar2 = 1e5;

        dir_SiPMHigh_Time_Detector[i]->cd();
        H_SiPM_ChannelTime_Walk_SiPM[i]->Write();
      }
      else
      {
        step = 1;
        mini = 4e4;
        maxi = 4e5;
        suppar2 = 0;

        dir_SiPMLow_Time_Detector[i]->cd();
        H_SiPM_ChannelTime_Walk_SiPM[i]->Write();
      }

      int counter = 0;
      // graph mean
      TGraphErrors *graph_mean = new TGraphErrors();
      graph_mean->SetTitle("Maximum");
      graph_mean->GetXaxis()->SetTitle("Channel");
      graph_mean->GetYaxis()->SetTitle("Maximum");
      graph_mean->GetXaxis()->CenterTitle();
      graph_mean->GetYaxis()->CenterTitle();

      vector<double> x;
      vector<double> y;

      TH1D *H_SiPM_ChannelTime_Walk_SiPM_YProjection = H_SiPM_ChannelTime_Walk_SiPM[i]->ProjectionY(("H_SiPM_ChannelTime_Walk_SiPM_YProjection_" + detectorName[i]).c_str());
      for (int bin = 0; bin < H_SiPM_ChannelTime_Walk_SiPM[i]->GetNbinsY(); bin += step)
      {
        TH1D *H_SiPM_ChannelTime_Walk_SiPM_XProjection = H_SiPM_ChannelTime_Walk_SiPM[i]->ProjectionX(("H_SiPM_ChannelTime_Walk_SiPM_XProjection_" + detectorName[i] + "_" + to_string(bin)).c_str(), bin, bin + step);
        double x_center = H_SiPM_ChannelTime_Walk_SiPM_XProjection->GetBinCenter(H_SiPM_ChannelTime_Walk_SiPM_XProjection->GetMaximumBin());
        H_SiPM_ChannelTime_Walk_SiPM_XProjection->GetXaxis()->SetRangeUser(x_center - 100, x_center + 100);
        TF1 *gauss = new TF1(("skewgauss_" + detectorName[i] + "_" + to_string(bin)).c_str(), skewedgauss, x_center - 100, x_center + 100, 5);
        gauss->SetParameter(3, H_SiPM_ChannelTime_Walk_SiPM_XProjection->GetBinContent(H_SiPM_ChannelTime_Walk_SiPM_XProjection->GetMaximumBin()));
        gauss->SetParameter(0, x_center);
        gauss->SetParameter(1, 10);
        gauss->SetParLimits(2, -100, 0);
        if (H_SiPM_ChannelTime_Walk_SiPM_XProjection->Integral() < 100)
          continue; /// avoid empty data for fit

        H_SiPM_ChannelTime_Walk_SiPM_XProjection->Fit(gauss, "QRN");

        if (gauss->GetParError(0) / abs(gauss->GetParameter(0)) > 0.2)
          continue; /// avoid weird data points due to low statistics

        x.push_back(H_SiPM_ChannelTime_Walk_SiPM_YProjection->GetBinCenter(bin + step / 2));
        y.push_back(gauss->GetMaximumX());
        graph_mean->SetPoint(counter, H_SiPM_ChannelTime_Walk_SiPM_YProjection->GetBinCenter(bin + step / 2), -gauss->GetMaximumX());
        graph_mean->SetPointError(counter, H_SiPM_ChannelTime_Walk_SiPM_YProjection->GetBinCenter(bin + step / 2) - H_SiPM_ChannelTime_Walk_SiPM_YProjection->GetBinCenter(bin), gauss->GetParError(0));
        counter++;

        TCanvas *c = new TCanvas(("C_SiPM_ChannelTime_Walk_SiPM_" + detectorName[i] + "_" + to_string(bin)).c_str(), ("C_SiPM_ChannelTime_Walk_SiPM_" + detectorName[i] + "_" + to_string(bin)).c_str(), 800, 400);
        H_SiPM_ChannelTime_Walk_SiPM_XProjection->Draw("HIST");
        gauss->Draw("SAME");
        // c->Write();
      }

      C_SiPM_ChannelTime_Walk_SiPM[i] = new TCanvas(("C_SiPM_ChannelTime_Walk_SiPM_" + detectorName[i]).c_str(), ("C_SiPM_ChannelTime_Walk_SiPM_" + detectorName[i]).c_str(), 800, 400);

      ////3rd
      C_SiPM_ChannelTime_Walk_SiPM[i]->cd();
      TPad *P_DiffRearStripvsStrip = new TPad("MeanRearSiPMvsT", "MeanRearSiPMvsT", 0.05, 0.05, 0.5, 0.5);
      P_DiffRearStripvsStrip->Draw();
      P_DiffRearStripvsStrip->cd();
      graph_mean->SetMarkerStyle(20);
      graph_mean->SetMarkerSize(0.5);
      graph_mean->Draw("AP");

      MeanAcceptance_Walk_SiPM[i] = new TF1("inverse", "[0] + [1]/sqrt(x+[2])", mini, maxi);
      MeanAcceptance_Walk_SiPM[i]->SetParameters(-15, 3000, 0);
      MeanAcceptance_Walk_SiPM[i]->SetParLimits(0, -20, 20);
      MeanAcceptance_Walk_SiPM[i]->SetParLimits(1, 0, 20000);
      MeanAcceptance_Walk_SiPM[i]->SetParLimits(2, -mini, suppar2);

      int status = 1;
      if (graph_mean->GetN() != 0)
      {
        status = graph_mean->Fit(MeanAcceptance_Walk_SiPM[i], "WQR");
        if (status != 0)
          Warning("Fit failed for " + detectorName[i] + " during walk SiPM correction.");
      }
      //
      // genereal fit canvas
      if (IsDetectorBetaLow(i))
      {
        C_SiPM_Walk_SiPM_Fits->cd(GetDetectorChannel(i));
        graph_mean->Draw("AP");
        MeanAcceptance_Walk_SiPM[i]->Draw("SAME");
      }

      // ////1st
      C_SiPM_ChannelTime_Walk_SiPM[i]->cd();
      TPad *P_RearStrip = new TPad("RearSiPM_ChannelTime", "RearSiPM_ChannelTime", 0.05, 0.5, 0.5, 1.0);
      P_RearStrip->Draw();
      P_RearStrip->cd();
      P_RearStrip->SetLogz();
      H_SiPM_ChannelTime_Walk_SiPM[i]->Draw("COLZ");

      maximum[i] = new TGraphErrors();
      for (int bin = 0; bin < H_SiPM_ChannelTime_Walk_SiPM[i]->GetNbinsY(); bin++)
      {
        double x_mean = -MeanAcceptance_Walk_SiPM[i]->Eval(H_SiPM_ChannelTime_Walk_SiPM[i]->GetYaxis()->GetBinCenter(bin));
        maximum[i]->SetPoint(bin, x_mean, H_SiPM_ChannelTime_Walk_SiPM[i]->GetYaxis()->GetBinCenter(bin));
      }
      maximum[i]->SetLineColor(kRed);
      maximum[i]->Draw("SAME");
    }
  }
}

inline void WriteHistograms_SiPMWalk()
{
  for (int i = 0; i < detectorNum; i++)
  {
    ProgressCounter(i, detectorNum, " <INFO>  Writing Histograms");
    if (IsDetectorBeta(i))
    {
      dir_Walk_SiPM->cd();

      // 2nd
      C_SiPM_ChannelTime_Walk_SiPM[i]->cd();
      TPad *P_SiPM_ChannelTime_Walk_SiPM = new TPad("SiPM_ChannelTime_Walk_SiPM", "SiPM_ChannelTime_Walk_SiPM", 0.55, 0.50, 1.0, 1.0);
      P_SiPM_ChannelTime_Walk_SiPM->Draw();
      P_SiPM_ChannelTime_Walk_SiPM->cd();
      P_SiPM_ChannelTime_Walk_SiPM->SetLogz();
      H_SiPM_ChannelTime_Walk_Corrected[i]->Draw("COLZ");
      C_SiPM_ChannelTime_Walk_SiPM[i]->Write();
    }
  }
}

inline void WriteHistograms_Cleaned()
{
  GROUPED_File->cd();
  // H_RearSiPM_MeanTime_Nearest_Group->Write();
  for (int i = 0; i < detectorNum; i++)
  {
    ProgressCounter(i, detectorNum, " <INFO>  Writing Histograms");
    if (IsDetectorSiliStrip(i))
    {
      // if (!FULL)
      // {
      dir_Group_Strip->cd();
      H_Group_Time[i]->Write();
      H_Group_ChannelTime[i]->Write();
      H_Channel_Group[i]->Write();
      TCanvas *c = new TCanvas(("C_Group_ChannelTime_" + detectorName[i]).c_str(), ("C_Group_ChannelTime_" + detectorName[i]).c_str(), 800, 400);
      H_Channel_RAW[i]->Draw("HIST");
      H_Channel_Group[i]->SetLineColor(kRed);
      H_Channel_Group[i]->Draw("SAME");
      c->Write();
      // }

      dir_Strip_Channel_Detector[i]->cd();
      H_Channel_Cleaned[i]->Write();
      H_RearMatched_Channel_Cleaned[i]->Write();
      H_RearMatchedMean_Channel_Cleaned[i]->Write();

      dir_Cleaned->cd();
      TCanvas *C_Strip_Cleaned = new TCanvas(("C_Strip_Cleaned_" + detectorName[i]).c_str(), ("C_Strip_Cleaned_" + detectorName[i]).c_str(), 800, 400);
      H_Channel_RAW[i]->Draw("HIST");
      H_Channel_Cleaned[i]->SetLineColor(kRed);
      H_Channel_Cleaned[i]->Draw("SAME");
      C_Strip_Cleaned->Write();

      dir_Strip_Time_Detector[i]->cd();
      H_RearSiPM_ChannelTime_dt[i]->Write();

      dir_Matching_Cleaned->cd();
      TCanvas *C_Strip_Matching_Cleaned = new TCanvas(("C_Strip_Matching_Cleaned_" + detectorName[i]).c_str(), ("C_Strip_Matching_Cleaned_" + detectorName[i]).c_str(), 800, 400);
      H_Channel_Cleaned[i]->Draw("HIST");
      H_RearMatched_Channel_Cleaned[i]->SetLineColor(kRed);
      H_RearMatched_Channel_Cleaned[i]->Draw("SAME");
      H_RearMatchedMean_Channel_Cleaned[i]->SetLineColor(kGreen);
      H_RearMatchedMean_Channel_Cleaned[i]->Draw("SAME");
      TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);
      legend->AddEntry(H_Channel_Cleaned[i], "Strip", "l");
      legend->AddEntry(H_RearMatched_Channel_Cleaned[i], "Rear", "l");
      legend->AddEntry(H_RearMatchedMean_Channel_Cleaned[i], "Mean", "l");
      legend->Draw("SAME");
      C_Strip_Matching_Cleaned->Write();
    }
    if (IsDetectorSiliBack(i))
    {
      // if (!FULL)
      // {
      dir_Group_Rear->cd();
      H_Group_Time[i]->Write();
      H_Group_ChannelTime[i]->Write();
      H_Channel_Group[i]->Write();
      TCanvas *c = new TCanvas(("C_Group_ChannelTime_" + detectorName[i]).c_str(), ("C_Group_ChannelTime_" + detectorName[i]).c_str(), 800, 400);
      H_Channel_RAW[i]->Draw("HIST");
      H_Channel_Group[i]->SetLineColor(kRed);
      H_Channel_Group[i]->Draw("SAME");
      c->Write();
      // }
      dir_Rear_Channel_Detector[i]->cd();
      H_Channel_Cleaned[i]->Write();
      H_Rear_Channel_IAS_Cleaned[i]->Write();
      H_Rear_Channel_IAScoinc_Cleaned[i]->Write();
      dir_Rear_Time_Detector[i]->cd();
      H_RearSiPM_ChannelTime_Walk_Corrected[i]->Write();
      H_RearSiPM_ChannelTime_Cleaned[i]->Write();
      H_RearSiPM_Time_Nearest[i]->Write();
      H_RearSiPM_ChannelTime_dt[i]->Write();

      C_IAS_Channel_Cleaned->cd(GetDetector(i));
      H_Rear_Channel_IAS_Cleaned[i]->SetFillColor(0);
      H_Rear_Channel_IAS_Cleaned[i]->SetLineColor(kBlack);
      H_Rear_Channel_IAS_Cleaned[i]->Rebin(10);
      H_Rear_Channel_IAS_Cleaned[i]->Draw("HIST");
      H_Rear_Channel_IAScoinc_Cleaned[i]->SetFillColor(0);
      H_Rear_Channel_IAScoinc_Cleaned[i]->SetLineColor(kRed);
      H_Rear_Channel_IAScoinc_Cleaned[i]->Rebin(10);
      H_Rear_Channel_IAScoinc_Cleaned[i]->Draw("SAME");
      TH1D *H_Rear_Channel_IASnocoinc_Cleaned_Clone = (TH1D *)H_Rear_Channel_IAS_Cleaned[i]->Clone(("H_Rear_Channel_IASnocoinc_Cleaned_Clone_" + detectorName[i]).c_str());
      H_Rear_Channel_IASnocoinc_Cleaned_Clone->Add(H_Rear_Channel_IAScoinc_Cleaned[i], -1);
      H_Rear_Channel_IASnocoinc_Cleaned_Clone->SetFillColor(0);
      H_Rear_Channel_IASnocoinc_Cleaned_Clone->SetLineColor(kBlue);
      H_Rear_Channel_IASnocoinc_Cleaned_Clone->Draw("SAME");
    }

    if (IsDetectorBetaHigh(i))
    {
      dir_SiPMHigh_Channel_Detector[i]->cd();
      H_SiPM_Channel_Alone[i]->Write();
      H_SiPM_Channel_Couple[i]->Write();
      H_2SiPM_Channel_Couple[i]->Write();
      H_Channel_Cleaned[i]->Write();

      dir_SiPMHigh_Time_Detector[i]->cd();
      H_SiPM_ChannelTime_Walk_Corrected[i]->Write();
      H_SiPM_ChannelTime_Cleaned[i]->Write();
      H_SiPM_Pair_Time[i]->Write();
      for (int det = 0; det < SIGNAL_MAX; det++)
      {
        if (IsDetectorSiliBack(det))
        {
          H_RearSiPM_ChannelTime_dt_det[det][i]->Write();
        }
      }

      H_SiPMLowHigh_Time_One[i]->Write();
      H_SiPMLowHigh_TimeChannel_One[i]->Write();
      H_SiPMLowHigh_Time_Multiple[i]->Write();

      for (int j = i + 1; j < SIGNAL_MAX; j++)
      {
        if (IsDetectorBetaHigh(j))
        {
          H_RearSiPM_Time_Nearest_Nearest[i][j]->Write();
          H_RearSiPM_Channel_Nearset_Nearest[i][j]->Write();
          H_RearSiPM_Channel_Nearest_Group[i][j]->Write();
          H_2SiPM_Time_Next[i][j]->Write();
        }
      }
    }

    if (IsDetectorBetaLow(i))
    {
      dir_SiPMLow_Channel_Detector[i]->cd();
      H_SiPM_Channel_Alone[i]->Write();
      H_SiPM_Channel_Couple[i]->Write();
      H_Channel_Cleaned[i]->Write();

      dir_SiPMLow_Time_Detector[i]->cd();
      H_SiPM_ChannelTime_Walk_Corrected[i]->Write();
      H_SiPM_ChannelTime_Cleaned[i]->Write();

      for (int j = i + 1; j < SIGNAL_MAX; j++)
      {
        if (IsDetectorBetaLow(j))
        {
          H_2SiPM_Time_Next[i][j]->Write();
        }
      }
    }
  }

  dir_SiPMHigh_Time->cd();
  TCanvas *C_SiPMHigh_Time = new TCanvas("C_2SiPMHigh_Time_Next", "C_2SiPMHigh_Time_Next", 800, 400);
  TCanvas *C_SiPMLow_Time = new TCanvas("C_2SiPMLow_Time_Next", "C_2SiPMLow_Time_Next", 800, 400);
  for (int i = 0; i < detectorNum; i++)
  {
    for (int j = i; j < detectorNum; j++)
    {
      if (IsDetectorBetaHigh(i) && IsDetectorBetaHigh(j))
      {
        C_SiPMHigh_Time->cd();
        H_2SiPM_Time_Next[i][j]->SetLineColor(i + j);
        H_2SiPM_Time_Next[i][j]->Draw("HIST SAME");
      }
      if (IsDetectorBetaLow(i) && IsDetectorBetaLow(j))
      {
        C_SiPMLow_Time->cd();
        H_2SiPM_Time_Next[i][j]->SetLineColor(i + j - 20);
        H_2SiPM_Time_Next[i][j]->Draw("HIST SAME");
      }
    }
  }
  C_SiPMHigh_Time->Write();
  C_SiPMLow_Time->Write();

  TCanvas *C_SiPM_Pair_Time = new TCanvas("C_SiPM_Pair_Time", "C_SiPM_Pair_Time", 800, 400);
  for (int i = 0; i < detectorNum; i++)
  {
    if (IsDetectorBetaHigh(i))
    {
      H_SiPM_Pair_Time[i]->SetLineColor(2 * i);
      H_SiPM_Pair_Time[i]->Draw("HIST SAME");
    }
  }
  C_SiPM_Pair_Time->Write();

  dir_SiPMHigh_Multiplicity->cd();
  H_SiPMHigh_Mulitplicity_Cleaned->Write();
  dir_SiPMLow_Multiplicity->cd();
  H_SiPMLow_Mulitplicity_Cleaned->Write();

  dir_Cleaned->cd();
  TCanvas *C_Multiplicity_Groups = new TCanvas("C_Multiplicity_Groups", "C_Multiplicity_Groups", 800, 400);
  H_Multiplicity_Groups->Draw("HIST");
  H_Multiplicity_Groups_IAS->SetLineColor(kRed);
  H_Multiplicity_Groups_IAS->Draw("SAME");
  C_Multiplicity_Groups->Write();

  TCanvas *C_Multiplicity_SubGroups = new TCanvas("C_Multiplicity_SubGroups", "C_Multiplicity_SubGroups", 800, 400);
  H_Multiplicity_SubGroups->Draw("HIST");
  H_Multiplicity_SubGroups_IAS->SetLineColor(kRed);
  H_Multiplicity_SubGroups_IAS->Draw("SAME");
  C_Multiplicity_SubGroups->Write();

  TCanvas *C_Time_SubGroups = new TCanvas("C_Time_SubGroups", "C_Time_SubGroups", 800, 400);  
  TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);
  for (int mul = 2; mul < BETA_SIZE+1; mul++)
  {
    H_Time_SubGroups[mul]->SetLineColor(mul);
    H_Time_SubGroups[mul]->Draw("HIST SAME");
    legend->AddEntry(H_Time_SubGroups[mul], to_string(mul).c_str(), "l");
  }
  legend->Draw("SAME");
  C_Time_SubGroups->Write();

  /// GENERAL FITS CANVAS ///
  dir_Fits->cd();
  C_Strip_Cutting_Fits->Write();
  C_Rear_Walk_Silicon_Fits->Write();
  C_SiPM_Walk_SiPM_Fits->Write();
}

inline int WriteTree_Grouped()
{
  GROUPED_File->cd();
  CLEANED_Tree->Write();
  for (int i = 0; i < detectorNum; i++)
  {
    if (IsDetectorSiliStrip(i))
      CLEANED_Tree_detector[i]->Write();
  }

  delete CLEANED_Tree;
  return 0;
}

void SearchForCoincidence(TTreeReaderArray<Signal> &signals)
{
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
    if (FLAG2021) signals[index].Label = YearConverter(signals[index].Label);
    int current_label = signals[index].Label;
    if (IsDetectorSiliBack(current_label))
    {
      Rear_Position.push_back(index);
      Rear_Multiplicity_RAW++;
      H_Channel_RAW[current_label]->Fill(signals[index].Channel);
      if (signals[index].Pileup == 1)
      {
        H_Channel_PileUp[current_label]->Fill(signals[index].Channel);
      }
    }
    else if (IsDetectorSiliStrip(current_label))
    {
      Strip_Position.push_back(index);
      Strip_Multiplicity_RAW++;
      H_Channel_RAW[current_label]->Fill(signals[index].Channel);
      if (signals[index].Pileup == 1)
      {
        H_Channel_PileUp[current_label]->Fill(signals[index].Channel);
      }
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
        H_RearStrip_Channel_RAW[signals[Strip_Position[index_strip]].Label]->Fill(signals[Rear_Position[index_rear]].Channel, signals[Strip_Position[index_strip]].Channel);
        H_RearStrip_Channel_RAW[signals[Rear_Position[index_rear]].Label]->Fill(signals[Rear_Position[index_rear]].Channel, signals[Strip_Position[index_strip]].Channel);
        Strip_Associated[index_strip] = true;
        Rear_Associated[index_rear] = true;
      }
    }
  }

  Signal THIRD_ALPHA_PEAK_FLAG = Signal();

  bool IAS = false;
  if (Rear_Position.size() != 0)
  {
    int Rear_Channel = signals[Rear_Position[0]].Channel;
    int Rear_Label = signals[Rear_Position[0]].Label;

    if (Rear_Channel > Rear_IAS[GetDetector(Rear_Label)].first && Rear_Channel < Rear_IAS[GetDetector(Rear_Label)].second)
    {
      IAS = true;
    }
  }

  // SAVING ALL THE PAIRS
  for (auto it : RearStrip_ASSOCIATED)
  {
    //////////////////////// TREE ////////////////////////
    GROUPED_Tree_Silicon.push_back(it.first);
    GROUPED_Tree_Silicon.push_back(it.second);
    for (int index_h = 0; index_h < SiPM_High_Position.size(); index_h++)
    {
      GROUPED_Tree_SiPMHigh.push_back(signals[SiPM_High_Position[index_h]]);
    }
    for (int index_l = 0; index_l < SiPM_Low_Position.size(); index_l++)
    {
      GROUPED_Tree_SiPMLow.push_back(signals[SiPM_Low_Position[index_l]]);
    }
    GROUPED_Tree->Fill();
    GROUPED_Tree_Silicon.clear();
    GROUPED_Tree_SiPMHigh.clear();
    GROUPED_Tree_SiPMLow.clear();
    ///////////////////////////////////////////////////////
  }

  ///////////////////////////////////////
  ///////// CASES A B C D E F G /////////
  ///////////////////////////////////////

  for (auto it : RearStrip_ASSOCIATED)
  {
    // Case A : SINGLE Rear and Strip associated (most common case)
    // if (RearStrip_ASSOCIATED.size() == 1 && Rear_Position.size() == 1 && Strip_Position.size() == 1)
    // {
      // cout << "CASE A" << endl;
      H_Channel_A[it.first.Label]->Fill(it.first.Channel);   // REAR
      H_Channel_A[it.second.Label]->Fill(it.second.Channel); // STRIP

      H_RearStrip_Channel_A[it.first.Label]->Fill(it.first.Channel, it.second.Channel);  // REAR-STRIP for rear plot
      H_RearStrip_Channel_A[it.second.Label]->Fill(it.first.Channel, it.second.Channel); // REAR-STRIP for strip plot

      H_Strip_Channel_DiffRear_A[it.second.Label]->Fill((it.first.Channel - it.second.Channel) / it.second.Channel);                                                // STRIP/REAR
      H_Strip_Channel_DiffRearvsStrip_A[it.second.Label]->Fill((it.first.Channel - it.second.Channel) / it.second.Channel, it.second.Channel); // STRIP/REAR vs Strip channel

      H_RearStrip_Time_A[it.first.Label]->Fill(it.first.Time - it.second.Time);
      H_RearStrip_Time_A[it.second.Label]->Fill(it.first.Time - it.second.Time);
    // }
  }

    // Case B : SAME Rear and DIFFERENT NEIGHBOURG Strip associated (intertrip)
   if (RearStrip_ASSOCIATED.size() == 2 && IsDetectorSiliInterStrip(RearStrip_ASSOCIATED[0].second.Label, RearStrip_ASSOCIATED[1].second.Label))
    {
      // cout << "CASE B" << endl;
      H_Channel_B[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[0].first.Channel);   // REAR
      H_Channel_B[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].second.Channel); // STRIP
      H_Channel_B[RearStrip_ASSOCIATED[1].second.Label]->Fill(RearStrip_ASSOCIATED[1].second.Channel); // STRIP

      H_RearStrip_Channel_B[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[0].first.Channel, RearStrip_ASSOCIATED[0].second.Channel);                                           // REAR-STRIP for strip
      H_RearStrip_Channel_B[RearStrip_ASSOCIATED[1].first.Label]->Fill(RearStrip_ASSOCIATED[1].first.Channel, RearStrip_ASSOCIATED[1].second.Channel);                                           // REAR-STRIP for strip
      H_RearStrip_Channel_B[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].first.Channel, RearStrip_ASSOCIATED[0].second.Channel);                                          // REAR-STRIP for rear
      H_RearStrip_Channel_B[RearStrip_ASSOCIATED[1].second.Label]->Fill(RearStrip_ASSOCIATED[1].first.Channel, RearStrip_ASSOCIATED[1].second.Channel);                                          // REAR-STRIP for rear
      H_Rear2Strip_Channel_B[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[1].first.Channel, RearStrip_ASSOCIATED[0].second.Channel + RearStrip_ASSOCIATED[1].second.Channel); // REAR-STRIP+STRIP
      H_2Strip_Channel_B[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].second.Channel, RearStrip_ASSOCIATED[1].second.Channel);                                            // STRIP-STRIP
      H_2Strip_Channel_B[RearStrip_ASSOCIATED[1].second.Label]->Fill(RearStrip_ASSOCIATED[1].second.Channel, RearStrip_ASSOCIATED[0].second.Channel);                                            // STRIP-STRIP

      H_RearStrip_Time_B[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[0].first.Time - RearStrip_ASSOCIATED[0].second.Time);
      H_RearStrip_Time_B[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[0].first.Time - RearStrip_ASSOCIATED[1].second.Time);
      H_RearStrip_Time_B[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].first.Time - RearStrip_ASSOCIATED[0].second.Time);
      H_RearStrip_Time_B[RearStrip_ASSOCIATED[1].second.Label]->Fill(RearStrip_ASSOCIATED[0].first.Time - RearStrip_ASSOCIATED[1].second.Time);
    }

    // Case C : DIFFERENT Rear and DIFFERENT Strip associated (NOT intertrip)

    else if (RearStrip_ASSOCIATED.size() == 2 && !IsDetectorSiliInterStrip(RearStrip_ASSOCIATED[0].second.Label, RearStrip_ASSOCIATED[1].second.Label))
    {
      // cout << "CASE C" << endl;
      H_Channel_C[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[0].first.Channel);   // REAR
      H_Channel_C[RearStrip_ASSOCIATED[1].first.Label]->Fill(RearStrip_ASSOCIATED[1].first.Channel);   // REAR
      H_Channel_C[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].second.Channel); // STRIP
      H_Channel_C[RearStrip_ASSOCIATED[1].second.Label]->Fill(RearStrip_ASSOCIATED[1].second.Channel); // STRIP

      H_RearStrip_Channel_C[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[0].first.Channel, RearStrip_ASSOCIATED[0].second.Channel);  // REAR-STRIP for strip
      H_RearStrip_Channel_C[RearStrip_ASSOCIATED[1].first.Label]->Fill(RearStrip_ASSOCIATED[1].first.Channel, RearStrip_ASSOCIATED[1].second.Channel);  // REAR-STRIP for strip
      H_RearStrip_Channel_C[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].first.Channel, RearStrip_ASSOCIATED[0].second.Channel); // REAR-STRIP for rear
      H_RearStrip_Channel_C[RearStrip_ASSOCIATED[1].second.Label]->Fill(RearStrip_ASSOCIATED[1].first.Channel, RearStrip_ASSOCIATED[1].second.Channel); // REAR-STRIP for rear
      H_2Strip_Channel_C[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].second.Channel, RearStrip_ASSOCIATED[1].second.Channel);   // STRIP-STRIP

      H_RearStrip_Time_C[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[0].first.Time - RearStrip_ASSOCIATED[0].second.Time);
      H_RearStrip_Time_C[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].first.Time - RearStrip_ASSOCIATED[1].second.Time);
      H_2Strip_Time_C[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].first.Time - RearStrip_ASSOCIATED[1].first.Time);

      H_2Strip_Label_C->Fill(RearStrip_ASSOCIATED[0].second.Label, RearStrip_ASSOCIATED[1].second.Label);
    }

    // Case D : SINGLE Rear and Strip + other strip

    else if (RearStrip_ASSOCIATED.size() == 1 && Rear_Position.size() == 1 && Strip_Position.size() == 2 && GetDetector(signals[Strip_Position[0]].Label) != GetDetector(signals[Strip_Position[1]].Label))
    {
      // cout << "CASE D" << endl;
      int index_other_strip = -1;
      if (signals[Strip_Position[0]].Label != RearStrip_ASSOCIATED[0].second.Label)
      {
        index_other_strip = Strip_Position[0];
      }
      else
      {
        index_other_strip = Strip_Position[1];
      }
      if (index_other_strip != -1)
      {
        if (RearStrip_ASSOCIATED[0].second.Channel > 20000 && RearStrip_ASSOCIATED[0].second.Channel < 36000)
        {
          THIRD_ALPHA_PEAK_FLAG = signals[index_other_strip];
          H_Channel_D[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].second.Channel);                                           // STRIP
          H_Channel_D[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[0].first.Channel);                                             // Rear
          H_2Strip_Channel_D[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[0].second.Channel, signals[index_other_strip].Channel); // STRIP-STRIP
          H_2Strip_Time_D[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[0].second.Time - signals[index_other_strip].Time);         // STRIP-STRIP
          H_2Strip_Time_D[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].second.Time - signals[index_other_strip].Time);        // STRIP-STRIP
          H_2Strip_Label_D->Fill(RearStrip_ASSOCIATED[0].second.Label, signals[index_other_strip].Label);                                            // STRIP-STRIP
        }
        else
        {
          H_Channel_E[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].second.Channel);                                           // STRIP
          H_Channel_E[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[0].first.Channel);                                             // Rear
          H_2Strip_Channel_E[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[0].second.Channel, signals[index_other_strip].Channel); // STRIP-STRIP
          H_2Strip_Time_E[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[0].second.Time - signals[index_other_strip].Time);         // STRIP-STRIP
          H_2Strip_Time_E[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].second.Time - signals[index_other_strip].Time);        // STRIP-STRIP
          H_2Strip_Label_E->Fill(RearStrip_ASSOCIATED[0].second.Label, signals[index_other_strip].Label);                                            // STRIP-STRIP
        }
      }
    }

    // Case F : SINGLE Rear
    else if (RearStrip_ASSOCIATED.size() == 0 && Rear_Position.size() == 1 && Strip_Position.size() == 0)
    {
      H_Channel_F[signals[Rear_Position[0]].Label]->Fill(signals[Rear_Position[0]].Channel); // REAR
    }

    // Case G : Not associated REAR and STRIP

    else if (RearStrip_ASSOCIATED.size() == 0 && Rear_Position.size() >= 1 && Strip_Position.size() >= 1)
    {
      // cout << "CASE G" << endl;
      for (int index_rear : Rear_Position)
      {
        for (int index_strip : Strip_Position)
        {
          H_Channel_G[signals[index_rear].Label]->Fill(signals[index_rear].Channel);              // REAR
          H_Channel_G[signals[index_strip].Label]->Fill(signals[index_strip].Channel);            // STRIP
          H_RearStrip_Channel_G->Fill(signals[index_rear].Channel, signals[index_strip].Channel); // REAR-STRIP
          H_RearStrip_Label_G->Fill(signals[index_strip].Label, signals[index_rear].Label);       // REAR-STRIP
          H_RearStrip_Time_G->Fill(signals[index_rear].Time - signals[index_strip].Time);
        }
      }
    }

    else
    {
      for (int index_strip : Strip_Position)
      {
        H_Channel_H[signals[index_strip].Label]->Fill(signals[index_strip].Channel); // STRIP
      }
      for (int index_rear : Rear_Position)
      {
        H_Channel_H[signals[index_rear].Label]->Fill(signals[index_rear].Channel); // REAR
      }
    }

    // checking if strip are well defined in faster
    if (RearStrip_ASSOCIATED.size() == 2 && GetDetector(RearStrip_ASSOCIATED[0].second.Label) == GetDetector(RearStrip_ASSOCIATED[1].second.Label))
    {
      if (RearStrip_ASSOCIATED[0].second.Label < RearStrip_ASSOCIATED[1].second.Label)
        H_2Strip_Label_Check->Fill(RearStrip_ASSOCIATED[0].second.Label, RearStrip_ASSOCIATED[1].second.Label);
      else
        H_2Strip_Label_Check->Fill(RearStrip_ASSOCIATED[1].second.Label, RearStrip_ASSOCIATED[0].second.Label);
    }
  
}

void CuttingGroups()
{
  // Cutting
  //// Strip
  double Strip_Channel = (*Silicon)[1].Channel;
  int Strip_Label = (*Silicon)[1].Label;
  double Strip_Time = (*Silicon)[1].Time;

  ////Rear
  double Rear_Channel = (*Silicon)[0].Channel;
  int Rear_Label = (*Silicon)[0].Label;
  double Rear_Time = (*Silicon)[0].Time;

  /// SiPM
  double SiPM_Time = 0;
  int SiPM_High_Size = SiPM_High->GetSize();
  int SiPM_Low_Size = SiPM_Low->GetSize();

  /// Selection
  double diff = (Rear_Channel - Strip_Channel) / Strip_Channel;
  double spread = SpreadAcceptance[Strip_Label]->Eval(Strip_Channel);
  double mean = MeanAcceptance[Strip_Label];

  if (diff > mean - 3 * spread && diff < mean + 3 * spread)
  {
    H_Channel_Cutted[Rear_Label]->Fill(Rear_Channel);   // REAR
    H_Channel_Cutted[Strip_Label]->Fill(Strip_Channel); // STRIP

    H_RearStrip_Channel_Cutted[Rear_Label]->Fill(Rear_Channel, Strip_Channel);      // REAR-STRIP for rear plot
    H_RearStrip_Channel_Cutted[Strip_Label]->Fill(Rear_Channel, Strip_Channel);     // REAR-STRIP for strip plot
    P_RearStrip_Channel_Cutted[Strip_Label]->AddPoint(Rear_Channel, Strip_Channel); // REAR-STRIP for strip plot TPROFILE for fit matching

    H_RearStrip_Time_Cutted[Rear_Label]->Fill(Rear_Time - Strip_Time);
    H_RearStrip_Time_Cutted[Strip_Label]->Fill(Rear_Time - Strip_Time);

    H_RearStrip_ChannelTime_Cutted[Strip_Label]->Fill(Rear_Time - Strip_Time, Strip_Channel);

    H_Strip_Channel_DiffRear_Cutted[Strip_Label]->Fill(diff);                       // STRIP/REAR
    H_Strip_Channel_DiffRearvsStrip_Cutted[Strip_Label]->Fill(diff, Strip_Channel); // STRIP/REAR vs Strip channel

    H_RearStrip_Time_Cutted[Rear_Label]->Fill(Rear_Time - Strip_Time);

    for (int i = 0; i < SiPM_High_Size; i++)
    {
      SiPM_Time = (*SiPM_High)[i].Time;
      H_Channel_Cutted[(*SiPM_High)[i].Label]->Fill((*SiPM_High)[i].Channel);
      H_RearSiPM_Time_Cutted[Rear_Label]->Fill(Rear_Time - SiPM_Time);
      H_RearSiPM_ChannelTime_Cutted[Rear_Label]->Fill(Rear_Time - SiPM_Time, Rear_Channel);
      H_SiPM_ChannelTime_Cutted[(*SiPM_High)[i].Label]->Fill(Rear_Time - SiPM_Time, (*SiPM_High)[i].Channel);
      H_RearSiPM_MulitplicityTime_Cutted[Rear_Label]->Fill(Rear_Time - SiPM_Time, SiPM_High_Size);
    }

    for (int i = 0; i < SiPM_Low_Size; i++)
    {
      SiPM_Time = (*SiPM_Low)[i].Time;
      H_Channel_Cutted[(*SiPM_Low)[i].Label]->Fill((*SiPM_Low)[i].Channel);
      H_RearSiPM_Time_Cutted[Rear_Label]->Fill(Rear_Time - SiPM_Time);
      H_RearSiPM_ChannelTime_Cutted[Rear_Label]->Fill(Rear_Time - SiPM_Time, Rear_Channel);
      H_SiPM_ChannelTime_Cutted[(*SiPM_Low)[i].Label]->Fill(Rear_Time - SiPM_Time, (*SiPM_Low)[i].Channel);
      H_RearSiPM_MulitplicityTime_Cutted[Rear_Label]->Fill(Rear_Time - SiPM_Time, SiPM_Low_Size);
    }

    if (Rear_Channel > Rear_IAS[GetDetector(Rear_Label)].first && Rear_Channel < Rear_IAS[GetDetector(Rear_Label)].second)
    {
      H_SiPMLow_Mulitplicity_Cutted->Fill(SiPM_Low_Size);
      H_SiPMHigh_Mulitplicity_Cutted->Fill(SiPM_High_Size);
    }

    for (auto &signal : *Silicon)
    {
      CUTTED_Tree_Silicon.push_back(signal);
    }
    for (auto &signal : *SiPM_High)
    {
      CUTTED_Tree_SiPMHigh.push_back(signal);
    }
    for (auto &signal : *SiPM_Low)
    {
      CUTTED_Tree_SiPMLow.push_back(signal);
    }
    CUTTED_Tree->Fill();
    CUTTED_Tree_Silicon.clear();
    CUTTED_Tree_SiPMHigh.clear();
    CUTTED_Tree_SiPMLow.clear();
  }
}

void SiliconWalkCorrection()
{

  int Rear_Label = (*Silicon)[0].Label;
  double Rear_Time = (*Silicon)[0].Time;
  double Rear_Channel = (*Silicon)[0].Channel;

  double mean_time = MeanAcceptance_Walk_Silicon[Rear_Label]->Eval(Rear_Channel);

  for (int i = 0; i < SiPM_High->GetSize(); i++)
  {
    double diff_time = Rear_Time - (*SiPM_High)[i].Time;
    H_RearSiPM_ChannelTime_Walk_Silicon[Rear_Label]->Fill(diff_time + mean_time, Rear_Channel);
    H_SiPM_ChannelTime_Walk_SiPM[(*SiPM_High)[i].Label]->Fill(diff_time + mean_time, (*SiPM_High)[i].Channel);
  }

  for (int i = 0; i < SiPM_Low->GetSize(); i++)
  {
    double diff_time = Rear_Time - (*SiPM_Low)[i].Time;
    H_RearSiPM_ChannelTime_Walk_Silicon[Rear_Label]->Fill(diff_time + mean_time, Rear_Channel);
    H_SiPM_ChannelTime_Walk_SiPM[(*SiPM_Low)[i].Label]->Fill(diff_time + mean_time, (*SiPM_Low)[i].Channel);
  }
}

inline void SiPMWalkCorrection()
{

  int Rear_Label = (*Silicon)[0].Label;
  double Rear_Time = (*Silicon)[0].Time;
  double Rear_Channel = (*Silicon)[0].Channel;

  double mean_time_silicon = MeanAcceptance_Walk_Silicon[Rear_Label]->Eval(Rear_Channel);

  for (int i = 0; i < SiPM_High->GetSize(); i++)
  {
    double diff_time = Rear_Time - (*SiPM_High)[i].Time;
    double mean_time_SiPM = MeanAcceptance_Walk_SiPM[(*SiPM_High)[i].Label]->Eval((*SiPM_High)[i].Channel);
    H_RearSiPM_ChannelTime_Walk_Corrected[Rear_Label]->Fill(diff_time + mean_time_silicon + mean_time_SiPM, Rear_Channel);
    H_SiPM_ChannelTime_Walk_Corrected[(*SiPM_High)[i].Label]->Fill(diff_time + mean_time_silicon + mean_time_SiPM, (*SiPM_High)[i].Channel);
  }

  for (int i = 0; i < SiPM_Low->GetSize(); i++)
  {
    double diff_time = Rear_Time - (*SiPM_Low)[i].Time;
    double mean_time_SiPM = MeanAcceptance_Walk_SiPM[(*SiPM_Low)[i].Label]->Eval((*SiPM_Low)[i].Channel);
    H_RearSiPM_ChannelTime_Walk_Corrected[Rear_Label]->Fill(diff_time + mean_time_silicon + mean_time_SiPM, Rear_Channel);
    H_SiPM_ChannelTime_Walk_Corrected[(*SiPM_Low)[i].Label]->Fill(diff_time + mean_time_silicon + mean_time_SiPM, (*SiPM_Low)[i].Channel);
  }
}

inline void CleaningGroups(TTreeReaderArray<Signal> &signals, int verbose=0)
{

  int TIME_BETWEEN_SAME_GAIN = 20;
  int TIME_BETWEEN_DIFFERENT_GAIN = 30;

  ////////////////////####################////////////////////
  /// Selecting CASE A ///
  vector<int> Rear_Position;
  vector<int> Strip_Position;

  vector<Signal> Silicon;
  vector<Signal> SiPM_High;
  vector<Signal> SiPM_Low;

  for (int index = 0; index < signals.GetSize(); index++)
  {
    if (FLAG2021) signals[index].Label = YearConverter(signals[index].Label);
    int current_label = signals[index].Label;
    if (IsDetectorSiliBack(current_label))
    {
      signals[index].Channel = signals[index].Channel - 0.5 + gRandom->Rndm();
      signals[index].Time = signals[index].Time + (-0.5 + gRandom->Rndm()) * 8;
      Rear_Position.push_back(index);
    }
    else if (IsDetectorSiliStrip(current_label))
    {
      signals[index].Channel = signals[index].Channel - 0.5 + gRandom->Rndm();
      signals[index].Time = signals[index].Time + (-0.5 + gRandom->Rndm()) * 8;
      Strip_Position.push_back(index);
    }
    else if (IsDetectorBetaHigh(current_label))
    {
      signals[index].Channel = signals[index].Channel - 0.5 + gRandom->Rndm();
      signals[index].Time = signals[index].Time + (-0.5 + gRandom->Rndm()) * 2;
      SiPM_High.push_back(signals[index]);
    }
    else if (IsDetectorBetaLow(current_label))
    {
      signals[index].Channel = signals[index].Channel - 0.5 + gRandom->Rndm();
      signals[index].Time = signals[index].Time + (-0.5 + gRandom->Rndm()) * 2;
      SiPM_Low.push_back(signals[index]);
    }
    if (IsHRSProton(current_label))
    {
      CLEANED_Tree_HRS = signals[index];
      CLEANED_Tree->Fill();
      CLEANED_Tree_HRS = Signal();
      return;
    }    
  }

  if (verbose == 1) Info("Silicon Stuff", 1);
  // Determining Rear strip couples
  vector<pair<Signal, Signal>> RearStrip_ASSOCIATED;
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

  if (RearStrip_ASSOCIATED.size() < 1)
    return;

  vector<bool> RearStrip_ASSOCIATED_ACCEPTED = vector<bool>(RearStrip_ASSOCIATED.size(), true);
  // ###### Cuts ####### //
  for (int index_pair = 0; index_pair < RearStrip_ASSOCIATED.size(); index_pair++)
  {
    /// RearStrip energy cut (3 sigmas)
    double diff = (RearStrip_ASSOCIATED[index_pair].first.Channel - RearStrip_ASSOCIATED[index_pair].second.Channel) / RearStrip_ASSOCIATED[index_pair].second.Channel;
    double spread = SpreadAcceptance[RearStrip_ASSOCIATED[index_pair].second.Label]->Eval(RearStrip_ASSOCIATED[index_pair].second.Channel);
    double mean = MeanAcceptance[RearStrip_ASSOCIATED[index_pair].second.Label];
    if (diff > mean + 3 * spread || diff < mean - 3 * spread) 
    {
      RearStrip_ASSOCIATED_ACCEPTED[index_pair] = false;
    }

    // PileUp rejection from faster flag
    else if (RearStrip_ASSOCIATED[index_pair].first.Pileup == 1 || RearStrip_ASSOCIATED[index_pair].second.Pileup == 1)
    {
      RearStrip_ASSOCIATED_ACCEPTED[index_pair] = false;
    }

    // PileUp rejection between groups (500µs cut on strip)
    else if (RearStrip_ASSOCIATED[0].second.Time - lastGroupTime[RearStrip_ASSOCIATED[0].second.Label] < 500e3)
    {
      lastGroupTime[RearStrip_ASSOCIATED[0].second.Label] = RearStrip_ASSOCIATED[0].second.Time;
      RearStrip_ASSOCIATED_ACCEPTED[index_pair] = false;
    }

    // RearStrip time cut (100 ns)
    else if (abs(RearStrip_ASSOCIATED[index_pair].first.Time - RearStrip_ASSOCIATED[index_pair].second.Time) > 100)
    {
      RearStrip_ASSOCIATED_ACCEPTED[index_pair] = false;
    }

  }
  
  // Only one pair in the group
  int counter_Accepted = 0;
  for (int index_pair = 0; index_pair < RearStrip_ASSOCIATED.size(); index_pair++)
  {
    if (RearStrip_ASSOCIATED_ACCEPTED[index_pair])
    {
      counter_Accepted++;
    }
  }

  if (counter_Accepted > 1 || counter_Accepted == 0)
  {
    // cout << "Rejected " << endl;
    return;
  }

  Silicon.push_back(RearStrip_ASSOCIATED[0].first);
  Silicon.push_back(RearStrip_ASSOCIATED[0].second);
  ////////////////////####################////////////////////


  /// Cutting data Rear/Strip ///
  //// Strip
  double Strip_Channel = (Silicon)[1].Channel;
  int Strip_Label = (Silicon)[1].Label;
  double Strip_Time = (Silicon)[1].Time;

  ////Rear
  double Rear_Channel = (Silicon)[0].Channel;
  int Rear_Label = (Silicon)[0].Label;
  double Rear_Time = (Silicon)[0].Time;

  /// SiPM
  double SiPM_Time = 0;
  int SiPM_High_Size = SiPM_High.size();
  int SiPM_Low_Size = SiPM_Low.size();

  /// Selection
  double diff = (Rear_Channel - Strip_Channel) / Strip_Channel;
  double spread = SpreadAcceptance[Strip_Label]->Eval(Strip_Channel);
  double mean = MeanAcceptance[Strip_Label];

  if (!FULL)
  {
    H_Channel_RAW[Strip_Label]->Fill(Strip_Channel);
    H_Channel_RAW[Rear_Label]->Fill(Rear_Channel);
  }

  int number = Strip_Label;

  bool IAS = false;

  if (Rear_Channel > Rear_IAS[GetDetector(Rear_Label)].first && Rear_Channel < Rear_IAS[GetDetector(Rear_Label)].second)
  {
    IAS = true;
  }

  if (IAS)
  {
    H_Group_Time[Strip_Label]->Fill(Strip_Time - lastGroupTime[number]);
    H_Channel_Group[Strip_Label]->Fill(Strip_Channel);
    H_Group_ChannelTime[Strip_Label]->Fill(Strip_Time - lastGroupTime[number], Strip_Channel);

    H_Group_Time[Rear_Label]->Fill(Rear_Time - lastGroupTime[Rear_Label]);
    H_Channel_Group[Rear_Label]->Fill(Rear_Channel);
    H_Group_ChannelTime[Rear_Label]->Fill(Rear_Time - lastGroupTime[Rear_Label], Rear_Channel);
  }

  lastGroupTime[number] = Rear_Time;

  if (!FULL)
  {
    H_Channel_Cutted[Strip_Label]->Fill(Strip_Channel);
    H_Channel_Cutted[Rear_Label]->Fill(Rear_Channel);
  }

  ////////////////////####################////////////////////

  ////////////////////####################////////////////////
  /// CONTINUOUS VALUE OF CHANNELS ///
  for (auto &signal : Silicon)
  {
    H_Channel_Cleaned[signal.Label]->Fill(signal.Channel);
    CLEANED_Tree_Silicon.push_back(signal);
  }

  for (auto &signal : SiPM_High)
  {
    // signal.Time = CLEANED_Tree_Silicon[0].Time - signal.Time + MeanAcceptance_Walk_SiPM[signal.Label]->Eval(signal.Channel) + MeanAcceptance_Walk_Silicon[CLEANED_Tree_Silicon[0].Label]->Eval(CLEANED_Tree_Silicon[0].Channel);
  }

  for (auto &signal : SiPM_Low)
  {
    // signal.Time = CLEANED_Tree_Silicon[0].Time - signal.Time + MeanAcceptance_Walk_SiPM[signal.Label]->Eval(signal.Channel) + MeanAcceptance_Walk_Silicon[CLEANED_Tree_Silicon[0].Label]->Eval(CLEANED_Tree_Silicon[0].Channel);
  }

  /////////////////////////////////////

  double mean_time_silicon = MeanAcceptance_Walk_Silicon[Rear_Label]->Eval(Rear_Channel); // silicon walk correction

  vector<Signal> SiPM_H[10]; // vector of High per Label
  vector<Signal> SiPM_L[10]; // vector of Low per Label
  //////////// NEW ANALYSIS ///////////
  

  // #######################
  // ######### (A) #########
  if (verbose == 1) Info("A", 1);
  // # Searching Time difference between 2 SiPMs High (IAS) #
  // # - plotting
  // # - grouping SiPM High per Label
  for (int i = 0; i < SiPM_High.size(); i++)
  {
    SiPM_H[GetDetectorChannel((SiPM_High)[i].Label)].push_back((SiPM_High)[i]);
    if (i > 1)
    {
      if ((SiPM_High)[i].Label < (SiPM_High)[i - 1].Label)
      {
        if (IAS)
          H_2SiPM_Time_Next[(SiPM_High)[i].Label][(SiPM_High)[i - 1].Label]->Fill((SiPM_High)[i].Time - (SiPM_High)[i - 1].Time);
      }
      else
      {
        if (IAS)
          H_2SiPM_Time_Next[(SiPM_High)[i - 1].Label][(SiPM_High)[i].Label]->Fill((SiPM_High)[i].Time - (SiPM_High)[i - 1].Time);
      }
    }
  }

  // # Searching Time difference between 2 SiPMs Low (IAS) #
  for (int i = 0; i < SiPM_Low.size(); i++)
  {
    SiPM_L[GetDetectorChannel((SiPM_Low)[i].Label)].push_back((SiPM_Low)[i]);
    if (i > 1)
    {
      if ((SiPM_Low)[i].Label < (SiPM_Low)[i - 1].Label)
      {
        if (IAS)
          H_2SiPM_Time_Next[(SiPM_Low)[i].Label][(SiPM_Low)[i - 1].Label]->Fill((SiPM_Low)[i].Time - (SiPM_Low)[i - 1].Time);
      }
      else
      {
        if (IAS)
          H_2SiPM_Time_Next[(SiPM_Low)[i - 1].Label][(SiPM_Low)[i].Label]->Fill((SiPM_Low)[i].Time - (SiPM_Low)[i - 1].Time);
      }
    }
  }
  // #######################
  // #######################

  // #######################
  // ######### (B) #########
  if (verbose == 1) Info("B", 1);
  // # Grouping SiPM High in subgroup with the time difference deduced in (A) #
  // # - Grouped in *SiPM_HGroups* #
  vector<vector<Signal>> SiPM_HGroups;
  for (int i = 0; i < SiPM_High.size(); i++)
  {
    // Initialisation
    if (i == 0)
    {
      SiPM_HGroups.push_back({(SiPM_High)[i]});
    }
    else
    {
      if ((abs(SiPM_HGroups[SiPM_HGroups.size() - 1][(SiPM_HGroups[SiPM_HGroups.size() - 1]).size() - 1].Time - (SiPM_High)[i].Time)) <= TIME_BETWEEN_SAME_GAIN)
      {
        // Same Group : adding
        SiPM_HGroups[SiPM_HGroups.size() - 1].push_back((SiPM_High)[i]);
      }
      else
      {
        // New Group : creating
        SiPM_HGroups.push_back({(SiPM_High)[i]});
      }
    }
  }
  // #######################
  // #######################

  // #######################
  // ######### (C) #########
  if (verbose == 1) Info("C", 1);
  // # Associating Low annd High in pair keeping the made groups #
  // # - Grouped in *SiPM_Groups* #

  vector<vector<pair<Signal, Signal>>> SiPM_Groups;

  // Loop on group
  for (int i_group = 0; i_group < SiPM_HGroups.size(); i_group++)
  {
    SiPM_Groups.push_back({});
    // Loop on High in the subgroup
    for (int i_high = 0; i_high < SiPM_HGroups[i_group].size(); i_high++)
    {
      SiPM_Groups[i_group].push_back(make_pair(SiPM_HGroups[i_group][i_high], Signal()));
      // Loop on Low to find the nearest in time with the same label
      double min_diff = 1000000;
      for (int i_low = 0; i_low < SiPM_Low.size(); i_low++)
      {
        if (GetDetectorChannel(SiPM_Low[i_low].Label) == GetDetectorChannel(SiPM_HGroups[i_group][i_high].Label))
        {
          // Correcting walk effect
          double high = -SiPM_HGroups[i_group][i_high].Time + MeanAcceptance_Walk_SiPM[SiPM_HGroups[i_group][i_high].Label]->Eval(SiPM_HGroups[i_group][i_high].Channel);
          double low = -SiPM_Low[i_low].Time + MeanAcceptance_Walk_SiPM[SiPM_Low[i_low].Label]->Eval(SiPM_Low[i_low].Channel);
          double diff = abs(low - high);
          if (min_diff > diff)
          {
            min_diff = diff;
            if (diff <= TIME_BETWEEN_DIFFERENT_GAIN)
            {
              SiPM_Groups[i_group][i_high].second = SiPM_Low[i_low];
            }
          }
        }
      }
    }
  }

  // # Plotting Alone High Signal #
  for (int i_group = 0; i_group < SiPM_Groups.size(); i_group++)
  {
    for (int i_pair = 0; i_pair < SiPM_Groups[i_group].size(); i_pair++)
    {
      if (!SiPM_Groups[i_group][i_pair].second.isValid)
      {
        H_SiPM_Channel_Alone[SiPM_Groups[i_group][i_pair].first.Label]->Fill(SiPM_Groups[i_group][i_pair].first.Channel);
      }
    }
  }

  // # Get index of associated Low #
  vector<bool> Associated_Low(SiPM_Low.size(), false);
  for (int i_group = 0; i_group < SiPM_Groups.size(); i_group++)
  {
    for (int i_pair = 0; i_pair < SiPM_Groups[i_group].size(); i_pair++)
    {
      if (SiPM_Groups[i_group][i_pair].second.isValid)
      {
        Associated_Low[SiPM_Groups[i_group][i_pair].second.Label] = true;
      }
    }
  }

  
  // # Get Alone Low to create a new pair in the group at the end #
  for (int i_group = 0; i_group < SiPM_Groups.size(); i_group++)
  {
    double first = 10e6;
    double last = -10e6;
    for (int i_pair = 0; i_pair < SiPM_Groups[i_group].size(); i_pair++)
    {
      double diff_time = Rear_Time - SiPM_Groups[i_group][i_pair].first.Time;
      double mean_time_SiPM = MeanAcceptance_Walk_SiPM[SiPM_Groups[i_group][i_pair].first.Label]->Eval(SiPM_Groups[i_group][i_pair].first.Channel);
      double time = -(diff_time + mean_time_silicon + mean_time_SiPM);

      if (first > time || first == 10e6)
      {
        first = time;
      }
      if (last < time || last == -10e6)
      {
        last = time;
      }
    }

    H_Time_SubGroups[SiPM_Groups[i_group].size()]->Fill(SiPM_Groups[i_group][SiPM_Groups[i_group].size()-1].first.Time-SiPM_Groups[i_group][0].first.Time);

    // Loop on All low
    for (int i_low = 0; i_low < SiPM_Low.size(); i_low++)
    {
      // Only if Alone
      if (!Associated_Low[SiPM_Low[i_low].Label])
      {
        double diff_time = Rear_Time - SiPM_Low[i_low].Time;
        double mean_time_SiPM = MeanAcceptance_Walk_SiPM[SiPM_Low[i_low].Label]->Eval(SiPM_Low[i_low].Channel);
        double low_alone = -(diff_time + mean_time_silicon + mean_time_SiPM);
        if (low_alone > first - TIME_BETWEEN_DIFFERENT_GAIN && low_alone < last + TIME_BETWEEN_DIFFERENT_GAIN)
        {
          for (int i_pair = 0; i_pair < SiPM_Groups[i_group].size(); i_pair++)
          {
            if (GetDetectorChannel(SiPM_Groups[i_group][i_pair].first.Label) == GetDetectorChannel(SiPM_Low[i_low].Label))
            {
              SiPM_Groups[i_group][i_pair].second = SiPM_Low[i_low];
              Associated_Low[SiPM_Low[i_low].Label] = true;
            }
          }
          if (!Associated_Low[SiPM_Low[i_low].Label])
          {
            SiPM_Groups[i_group].push_back(make_pair(Signal(), SiPM_Low[i_low]));
            Associated_Low[SiPM_Low[i_low].Label] = true;
          }
        }
      }
    }
  }


  // # Plotting High Low in groups #
  // - time difference
  // - channel (1D and 2D) both and alone
  for (int i_group = 0; i_group < SiPM_Groups.size(); i_group++)
  {
    for (int i_pair = 0; i_pair < SiPM_Groups[i_group].size(); i_pair++)
    {
      // BOTH
      if (SiPM_Groups[i_group][i_pair].second.isValid && SiPM_Groups[i_group][i_pair].first.isValid)
      {
        double high = -SiPM_Groups[i_group][i_pair].first.Time + MeanAcceptance_Walk_SiPM[SiPM_Groups[i_group][i_pair].first.Label]->Eval(SiPM_Groups[i_group][i_pair].first.Channel);
        double low = -SiPM_Groups[i_group][i_pair].second.Time + MeanAcceptance_Walk_SiPM[SiPM_Groups[i_group][i_pair].second.Label]->Eval(SiPM_Groups[i_group][i_pair].second.Channel);
        double diff = low - high;

        H_SiPM_Pair_Time[SiPM_Groups[i_group][i_pair].first.Label]->Fill(diff);
        H_2SiPM_Channel_Couple[SiPM_Groups[i_group][i_pair].first.Label]->Fill(SiPM_Groups[i_group][i_pair].second.Channel, SiPM_Groups[i_group][i_pair].first.Channel);
      }

      // Alone Low
      if (!SiPM_Groups[i_group][i_pair].first.isValid)
      {
        H_SiPM_Channel_Alone[SiPM_Groups[i_group][i_pair].second.Label]->Fill(SiPM_Groups[i_group][i_pair].second.Channel);
      }
      // High 
      else
      {
        H_Channel_Cleaned[SiPM_Groups[i_group][i_pair].first.Label]->Fill(SiPM_Groups[i_group][i_pair].first.Channel);
      }
      
      // Alone High
      if (!SiPM_Groups[i_group][i_pair].second.isValid)
      {
        H_SiPM_Channel_Alone[SiPM_Groups[i_group][i_pair].first.Label]->Fill(SiPM_Groups[i_group][i_pair].first.Channel);
      }
      // Low
      else
      {
        H_Channel_Cleaned[SiPM_Groups[i_group][i_pair].second.Label]->Fill(SiPM_Groups[i_group][i_pair].second.Channel);
      }
    }
  }

  // # Plotting MULTIPLICITY of groups, and subgroups #
  H_Multiplicity_Groups->Fill(SiPM_Groups.size());
  if (IAS) H_Multiplicity_Groups_IAS->Fill(SiPM_Groups.size());
  
  for (int i_group = 0; i_group < SiPM_Groups.size(); i_group++)
  {
    H_Multiplicity_SubGroups->Fill(SiPM_Groups[i_group].size());
    if (IAS) H_Multiplicity_SubGroups_IAS->Fill(SiPM_Groups[i_group].size());
  }

  // #####################################
  // ######### DISPLAYING GROUPS #########

  int Valid_Pair_Counter = 0;
  for (int i_group = 0; i_group < SiPM_Groups.size(); i_group++)
  {
    for (int i_pair = 0; i_pair < SiPM_Groups[i_group].size(); i_pair++)
    {
      if (SiPM_Groups[i_group][i_pair].second.isValid && SiPM_Groups[i_group][i_pair].first.isValid)
      {
        Valid_Pair_Counter++;
      }
    }
  }

  bool Verbose = false;
  bool Found = false;

  if (IAS)
  {
    counter_all_IAS += 1;
  }

  // # TIME FASTER SORTED #
  if (IAS && SiPM_Groups.size() > 0)
  {
     if (Verbose)
    {
      cout << BLUE << endl;
      cout << "----------------------------------------------------------------------------" << endl;
      cout << "High Signals : " << SiPM_High_Size << "     Low Signals : " << SiPM_Low_Size << endl;
      cout << "SiPM Valid Pair : " << Valid_Pair_Counter << "     (High Alone : " << SiPM_High_Size - Valid_Pair_Counter << "     Low Alone : " << SiPM_Low_Size - Valid_Pair_Counter << ")" << endl;
      cout << "SiPM Group : " << SiPM_Groups.size() << endl;
      cout << "----------------------------------------------------------------------------" << WHITE << endl;
      for (int i_signal = 0; i_signal < signals.GetSize(); i_signal++)
      {
        Found = false;
        if (!IsDetectorBeta(signals[i_signal].Label))
          continue;

        // correcting time
        double diff_time = Rear_Time - signals[i_signal].Time;
        double mean_time_SiPM = MeanAcceptance_Walk_SiPM[signals[i_signal].Label]->Eval(signals[i_signal].Channel);

        for (int i_group = 0; i_group < SiPM_Groups.size(); i_group++)
        {
          for (int i_pair = 0; i_pair < SiPM_Groups[i_group].size(); i_pair++)
          {
            if (signals[i_signal] == SiPM_Groups[i_group][i_pair].first || signals[i_signal] == SiPM_Groups[i_group][i_pair].second)
            {
              if (!SiPM_Groups[i_group][i_pair].first.isValid)
                break;
              cout << left << GREEN
                   << setw(2) << "(G" << i_group << ") "
                   << setw(2) << "(P" << i_pair << ") "
                   << setw(8) << signals[i_signal].Label
                   << setw(10) << "Channel : "
                   << setw(15) << signals[i_signal].Channel
                   << setw(6) << "Time : "
                   << setw(15) << signals[i_signal].Time - signals[Rear_Position[0]].Time
                   << setw(6) << "Corrected Time : "
                   << setw(15)
                   << -(diff_time + mean_time_silicon + mean_time_SiPM)
                   << WHITE << endl;
              Found = true;
              break;
            }
          }
        }
        if (!Found)
        {
          cout << left << RED << "          "
               << setw(8) << signals[i_signal].Label
               << setw(10) << "Channel : "
               << setw(15) << signals[i_signal].Channel
               << setw(6) << "Time : "
               << setw(15) << signals[i_signal].Time - signals[Rear_Position[0]].Time
               << setw(6) << "Corrected Time : "
               << setw(15)
               << -(diff_time + mean_time_silicon + mean_time_SiPM)
               << WHITE << endl;
        }
      }

      cout << BLUE << "----------------------------------------------------------------------------" << WHITE << endl;

      // # GROUP PAIR SORTED #
      for (int i_group = 0; i_group < SiPM_Groups.size(); i_group++)
      {
        cout << left << BLUE
             << setw(2) << "Group n°" << i_group << "   " << "Group Size : " << SiPM_Groups[i_group].size();

        // Mean time group (on High)
        double mean_time = 0;
        int counter = 0;
        for (int i_pair = 0; i_pair < SiPM_Groups[i_group].size(); i_pair++)
        {
          if (!SiPM_Groups[i_group][i_pair].first.isValid)
            continue;
          double diff_time = Rear_Time - SiPM_Groups[i_group][i_pair].first.Time;
          double mean_time_SiPM = MeanAcceptance_Walk_SiPM[SiPM_Groups[i_group][i_pair].first.Label]->Eval(SiPM_Groups[i_group][i_pair].first.Channel);
          double t = -(diff_time + mean_time_silicon + mean_time_SiPM);
          // t different than nan
          if (!isnan(t))
          {
            counter++;
            mean_time += t;
          }
        }

        cout << "     " << setw(4) << "t = "
             << setw(15) << mean_time / counter << WHITE << endl;

        for (int i_pair = 0; i_pair < SiPM_Groups[i_group].size(); i_pair++)
        {

          if (SiPM_Groups[i_group][i_pair].first.isValid)
          {
            cout << GREEN << left
                 << "-- " << setw(4) << GetDetectorChannel(SiPM_Groups[i_group][i_pair].first.Label)
                 << setw(10) << "Channel High : "
                 << setw(15) << SiPM_Groups[i_group][i_pair].first.Channel;

            if (SiPM_Groups[i_group][i_pair].second.isValid)
            {

              cout << setw(10) << "Channel Low : "
                   << setw(15) << SiPM_Groups[i_group][i_pair].second.Channel;

              double high = -SiPM_Groups[i_group][i_pair].first.Time + MeanAcceptance_Walk_SiPM[SiPM_Groups[i_group][i_pair].first.Label]->Eval(SiPM_Groups[i_group][i_pair].first.Channel);
              double low = -SiPM_Groups[i_group][i_pair].second.Time + MeanAcceptance_Walk_SiPM[SiPM_Groups[i_group][i_pair].second.Label]->Eval(SiPM_Groups[i_group][i_pair].second.Channel);
              double diff = abs(low - high);

              cout << setw(6) << "ΔT : "
                   << setw(15) << diff
                   << WHITE << endl;
            }
            else
            {
              cout << YELLOW << "......................................... "
                   << WHITE << endl;
            }
          }
          else
          {
            cout << YELLOW << left
                 << "-- " << setw(34) << GetDetectorChannel(SiPM_Groups[i_group][i_pair].second.Label)
                 << setw(10) << "Channel Low : "
                 << setw(15) << SiPM_Groups[i_group][i_pair].second.Channel
                 << WHITE << endl;
          }
        }
      }
    }
  }

  /////////////////////////////////////

  // #######################
  // ######### (D) #########
  // # Time correction applied on SiPM #
  // # Saving in Tree #
  if (IAS)
  {
    H_RearSiPM_ChannelTime_dt[Rear_Label]->Fill((Silicon)[0].dt, Rear_Channel);
    H_RearSiPM_ChannelTime_dt[Strip_Label]->Fill((Silicon)[1].dt-(Silicon)[0].dt+8100, Strip_Channel);
  }


  for (int i_group = 0; i_group < SiPM_Groups.size(); i_group++)
  {
    for (int i_pair = 0; i_pair < SiPM_Groups[i_group].size(); i_pair++)
    {
      if (SiPM_Groups[i_group][i_pair].first.isValid)
      {
        double diff_time = Rear_Time - SiPM_Groups[i_group][i_pair].first.Time;
        double mean_time_SiPM = MeanAcceptance_Walk_SiPM[SiPM_Groups[i_group][i_pair].first.Label]->Eval(SiPM_Groups[i_group][i_pair].first.Channel);
        double t = -(diff_time + mean_time_silicon + mean_time_SiPM);
        SiPM_Groups[i_group][i_pair].first.Time = t;

        if (IAS)
        {
          H_RearSiPM_ChannelTime_dt_det[Rear_Label][SiPM_Groups[i_group][i_pair].first.Label]->Fill(t + (Silicon)[0].dt, SiPM_Groups[i_group][i_pair].first.Channel);
        }

        // if ((t + (Silicon)[0].dt) < 8000 || (t + (Silicon)[0].dt) > 8200);
        // {
        //   SiPM_Groups[i_group][i_pair].first.isValid = false;
        //   SiPM_Groups[i_group][i_pair].second.isValid = false;
        // }
      }
      if (SiPM_Groups[i_group][i_pair].second.isValid)
      {
        double diff_time = Rear_Time - SiPM_Groups[i_group][i_pair].second.Time;
        double mean_time_SiPM = MeanAcceptance_Walk_SiPM[SiPM_Groups[i_group][i_pair].second.Label]->Eval(SiPM_Groups[i_group][i_pair].second.Channel);
        double t = -(diff_time + mean_time_silicon + mean_time_SiPM);
        SiPM_Groups[i_group][i_pair].second.Time = t;
      }
    }
  }

  CLEANED_Tree_SiPMGroup = SiPM_Groups;
  
 

  CLEANED_Tree->Fill();
  Tree_Channel_detector = (Silicon)[1].Channel;
  CLEANED_Tree_detector[(Silicon)[1].Label]->Fill();
  CLEANED_Tree_Silicon.clear();
  CLEANED_Tree_SiPMGroup.clear();

  //////////////////////////////////////
}