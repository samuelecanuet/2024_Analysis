#ifndef GROUPER_HH
#define GROUPER_HH

#include "../../../lib/SignalDict/Signal.h"
#include "Detectors.hh"


using namespace std;

int VERBOSE = 0;
int Run;

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

TTree *SILICON_WALK_Tree;
vector<Signal> SILICON_WALK_Tree_Silicon;
vector<Signal> SILICON_WALK_Tree_SiPMHigh;
vector<Signal> SILICON_WALK_Tree_SiPMLow;

TTree *CLEANED_Tree;
TTree *CLEANED_Tree_detector[SIGNAL_MAX];
vector<Signal> CLEANED_Tree_Silicon;
vector<Signal> CLEANED_Tree_SiPMHigh;
vector<Signal> CLEANED_Tree_SiPMLow;
double Tree_Channel_detector;

TTreeReaderArray<Signal> *Silicon;
TTreeReaderArray<Signal> *SiPM_High;
TTreeReaderArray<Signal> *SiPM_Low;

////////////// HISTOGRAMS ////////////////
// RAW
TH1D *H_Strip_Mulitplicity_RAW;
TH1D *H_Rear_Mulitplicity_RAW;
TH1D *H_Channel_RAW[SIGNAL_MAX];
TH2D *H_RearStrip_Channel_RAW[SIGNAL_MAX];
TH1D *H_SiPMHigh_Mulitplicity_RAW;
TH1D *H_SiPMLow_Mulitplicity_RAW;

//PileUp
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
TH1D *H_SiPM_Mulitplicity_Cleaned[SIGNAL_MAX];
pair<int, int> Rear_IAS[9] = {make_pair(0, 0), make_pair(40000, 44500), make_pair(36000, 39500), make_pair(36000, 39500), make_pair(33000, 36500), make_pair(36000, 39000), make_pair(33000, 36500), make_pair(35500, 38500), make_pair(35000, 38500)};
pair<int, int> peaks_window_F[SIGNAL_MAX];
TH1D *H_Rear_Channel_IAS_Cleaned[SIGNAL_MAX];
TH1D *H_Rear_Channel_IAScoinc_Cleaned[SIGNAL_MAX];
vector<int> CoincidenceTime;
int CoincidenceTime_MAX = 26;
int CoincidenceTime_BIN = 26;
int CoincidenceTime_Left_Window = -50;
int CoincidenceTime_Right_Window = 50;
TH1D *H_HighSiPM_Time_IAS_Cleaned[99][99][SIGNAL_MAX];
TH1D *H_HighSiPM_Time_IAS_UnCleaned[99][99][SIGNAL_MAX];
TH1D *H_Energy_IAS;
TH1D *H_Energy_IAS_Cleaned[99][99][SIGNAL_MAX];
TH1D *H_Energy_IAS_UnCleaned[99][99][SIGNAL_MAX];
TH1D *H_LowSiPM_Time_IAS_Cleaned[99][99][SIGNAL_MAX];
TH1D *H_LowSiPM_Time_IAS_UnCleaned[99][99][SIGNAL_MAX];
TCanvas *C_IAS_Channel_Cleaned;
TH2D *H_2SiPM_Channel_Cleaned[SIGNAL_MAX];
// TGraphErrors *G_HighSiPM_FakeEvent_Multiplicity[SIGNAL_MAX];
// TGraphErrors *G_LowSiPM_FakeEvent_Multiplicity[SIGNAL_MAX];
// TGraphErrors *G_HighSiPM_Losses_Multiplicity[SIGNAL_MAX];
// TGraphErrors *G_LowSiPM_Losses_Multiplicity[SIGNAL_MAX];
TH2D *H_HighSiPM_FakeEvent_Multiplicity[SIGNAL_MAX];
TH2D *H_LowSiPM_FakeEvent_Multiplicity[SIGNAL_MAX];
TH2D *H_HighSiPM_Losses_Multiplicity[SIGNAL_MAX];
TH2D *H_LowSiPM_Losses_Multiplicity[SIGNAL_MAX];
int counter_graph_fake_events = 0;
TH1D* H_2SiPM_Time_IAS_Cleaned[SIGNAL_MAX];
TH1D* H_SiPMLowHigh_Time_One[SIGNAL_MAX];
TH1D* H_SiPMLowHigh_Time_Multiple[SIGNAL_MAX];
TH1D* H_2SiPM_Time_One[SIGNAL_MAX][SIGNAL_MAX];
TH2D* H_SiPMLowHigh_TimeChannel_One[SIGNAL_MAX];
TH1D* H_2SiPM_Time_Multiple[SIGNAL_MAX][SIGNAL_MAX];
TH1D* H_RearSiPM_Time_Nearest[SIGNAL_MAX];
TH1D* H_RearSiPM_Time_Nearest_Nearest[SIGNAL_MAX][SIGNAL_MAX];
TH2D* H_RearSiPM_Channel_Nearset_Nearest[SIGNAL_MAX][SIGNAL_MAX];
TH1D* H_RearSiPM_MeanTime_Nearest_Group;
TH2D* H_RearSiPM_Channel_Nearest_Group[SIGNAL_MAX][SIGNAL_MAX];
TH1D* H_2SiPM_Time_Next[SIGNAL_MAX][SIGNAL_MAX];
TH1D* H_2SiPM_Time_MeanFurthest[SIGNAL_MAX];
TH1D* H_RearSiPM_Time_New_Nearest[SIGNAL_MAX];
TH1D* H_RearSiPM_Time_New_Mean[SIGNAL_MAX];
int counter_graph_fake_events_multiplicity = 0;


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
//////////////////////////////////////

TGraphErrors *spread_up[SIGNAL_MAX];
TGraphErrors *spread_down[SIGNAL_MAX];
TGraphErrors *maximum[SIGNAL_MAX];

///// DIRECTORY //////////////////
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
TDirectory *dir_SiPMHighRear_TimeChannel;
TDirectory *dir_SiPMHigh_Multiplicity;
TDirectory *dir_SiPMHigh_Channel_Detector[SIGNAL_MAX];
TDirectory *dir_SiPMHigh_Time_Detector[SIGNAL_MAX];

TDirectory *dir_SiPMLow;
TDirectory *dir_SiPMLow_Channel;
TDirectory *dir_SiPMLow_Time;
TDirectory *dir_SiPMLowRear_TimeChannel;
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

///////// FUCNTION ARRAY /////////////
TF1 *SpreadAcceptance[SIGNAL_MAX];
double MeanAcceptance[SIGNAL_MAX];
TF1 *MeanAcceptance_Walk_Silicon[SIGNAL_MAX];
TF1 *MeanAcceptance_Walk_SiPM[SIGNAL_MAX];
pair<double, double> SpreadAcceptance_Walk_SiPM[SIGNAL_MAX];

//////// calibration from Calibration.cc  ////////
double ManualCalibFitted[SIGNAL_MAX][3];

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
      int min1=0;
      int max1=0;
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

  CalibFileName = "./Config_Files/Calibration.txt";

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

//////////////////////////////////////

inline int InitHistograms_Grouped()
{

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
  dir_SiPMHighRear_TimeChannel = dir_SiPM_High->mkdir("SiPMHighRear_TimeChannel");
  dir_SiPMHigh_Multiplicity = dir_SiPM_High->mkdir("SiPMHigh_Multiplicity");

  ////SiPM Low
  dir_SiPMLow = GROUPED_File->mkdir("SiPMLow");
  dir_SiPMLow_Channel = dir_SiPMLow->mkdir("SiPMLow_Channel");
  dir_SiPMLow_Time = dir_SiPMLow->mkdir("SiPMLow_Time");
  dir_SiPMLowRear_TimeChannel = dir_SiPMLow->mkdir("SiPMLowRear_TimeChannel");
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
  
  for (int i = 0; i <= MAX_MULTIPLICTY; i++)
  {
    dir_Cleaned_HighSiPM_TimeCoincidence_Mulitplicity[i] = dir_Cleaned->mkdir(("HighSiPM_TimeCoincidence_M" + to_string(i)).c_str());
    dir_Cleaned_LowSiPM_TimeCoincidence_Mulitplicity[i] = dir_Cleaned->mkdir(("LowSiPM_TimeCoincidence_M" + to_string(i)).c_str());
  }

  for (int i = CoincidenceTime_BIN; i <= CoincidenceTime_MAX; i += CoincidenceTime_BIN) {
    CoincidenceTime.push_back(i);
  }

  // FIT RESULTS
  dir_Fits = GROUPED_File->mkdir("Fits");
  C_Strip_Cutting_Fits = new TCanvas("C_Strip_Cutting_Fits", "C_Strip_Cutting_Fits", 800, 800);
  C_Strip_Cutting_Fits->Divide(SILI_SIZE - 1, SILI_NUM);
  C_Rear_Walk_Silicon_Fits = new TCanvas("C_Rear_Walk_Silicon_Fits", "C_Rear_Walk_Silicon_Fits", 800, 800);
  C_Rear_Walk_Silicon_Fits->Divide(4, 2);
  C_SiPM_Walk_SiPM_Fits = new TCanvas("C_SiPM_Walk_SiPM_Fits", "C_SiPM_Walk_SiPM_Fits", 800, 800);
  C_SiPM_Walk_SiPM_Fits->Divide(4, 2);


  for (size_t i = 0; i < detectorNum; ++i)
  {
    if (IsDetectorSiliStrip(i))
    {
      dir_Strip_Channel_Detector[i] = dir_Strip_Channel->mkdir(detectorName[i].c_str());
      dir_Strip_Time_Detector[i] = dir_Strip_Time->mkdir(detectorName[i].c_str());

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

      H_RearStrip_Channel_Cutted[i] = new TH2D(("RearStrip_Channel_Cutted_" + detectorName[i]).c_str(), ("RearStrip_Channel_Cutted_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
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
    }

    if (IsDetectorSiliBack(i))
    {

      dir_Rear_Channel_Detector[i] = dir_Rear_Channel->mkdir(detectorName[i].c_str());
      dir_Rear_Time_Detector[i] = dir_Rear_Time->mkdir(detectorName[i].c_str());

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

      H_SiPM_Mulitplicity_Cleaned[i] = new TH1D(("SiPM_Mulitplicity_Cleaned_" + detectorName[i]).c_str(), ("SiPM_Mulitplicity_Cleaned_" + detectorName[i]).c_str(), 20, 0, 20);
      H_SiPM_Mulitplicity_Cleaned[i]->GetXaxis()->SetTitle("Multiplicity");
      H_SiPM_Mulitplicity_Cleaned[i]->GetYaxis()->SetTitle("Counts");
      H_SiPM_Mulitplicity_Cleaned[i]->GetXaxis()->CenterTitle();
      H_SiPM_Mulitplicity_Cleaned[i]->GetYaxis()->CenterTitle();

      H_2SiPM_Channel_Cleaned[i] = new TH2D(("2SiPM_Channel_Cleaned_" + detectorName[i]).c_str(), ("2SiPM_Channel_Cleaned_" + detectorName[i]).c_str(), eHighN / 10, eHighMin, eHighMax, eHighN / 10, eHighMin, eHighMax);
      H_2SiPM_Channel_Cleaned[i]->GetXaxis()->SetTitle("SiPM 1 Channel");
      H_2SiPM_Channel_Cleaned[i]->GetYaxis()->SetTitle(("SiPM " + to_string(GetDetectorChannel(i)) + " Channel").c_str());
      H_2SiPM_Channel_Cleaned[i]->GetXaxis()->CenterTitle();
      H_2SiPM_Channel_Cleaned[i]->GetYaxis()->CenterTitle();

      H_2SiPM_Time_IAS_Cleaned[i] = new TH1D(("2SiPM_Time_IAS_Cleaned_" + detectorName[i]).c_str(), ("2SiPM_Time_IAS_Cleaned_" + detectorName[i]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
      H_2SiPM_Time_IAS_Cleaned[i]->GetXaxis()->SetTitle("Time between 2 SiPMs (ns)");
      H_2SiPM_Time_IAS_Cleaned[i]->GetYaxis()->SetTitle("Counts");
      H_2SiPM_Time_IAS_Cleaned[i]->GetXaxis()->CenterTitle();
      H_2SiPM_Time_IAS_Cleaned[i]->GetYaxis()->CenterTitle();

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

      for (int j = i; j < SIGNAL_MAX; j++)
      {
        if (IsDetectorBetaHigh(j))
        {
          H_2SiPM_Time_One[i][j] = new TH1D(("2SiPM_Time_One_" + detectorName[i] + "_" + detectorName[j]).c_str(), ("2SiPM_Time_One_" + detectorName[i] + "_" + detectorName[j]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
          H_2SiPM_Time_One[i][j]->GetXaxis()->SetTitle("Time between 2 SiPMs (ns)");
          H_2SiPM_Time_One[i][j]->GetYaxis()->SetTitle("Counts");
          H_2SiPM_Time_One[i][j]->GetXaxis()->CenterTitle();
          H_2SiPM_Time_One[i][j]->GetYaxis()->CenterTitle();

          H_2SiPM_Time_Multiple[i][j] = new TH1D(("2SiPM_Time_Multiple_" + detectorName[i] + "_" + detectorName[j]).c_str(), ("2SiPM_Time_Multiple_" + detectorName[i] + "_" + detectorName[j]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
          H_2SiPM_Time_Multiple[i][j]->GetXaxis()->SetTitle("Time between 2 SiPMs (ns)");
          H_2SiPM_Time_Multiple[i][j]->GetYaxis()->SetTitle("Counts");
          H_2SiPM_Time_Multiple[i][j]->GetXaxis()->CenterTitle();
          H_2SiPM_Time_Multiple[i][j]->GetYaxis()->CenterTitle();

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

      H_SiPM_Mulitplicity_Cleaned[i] = new TH1D(("SiPM_Mulitplicity_Cleaned_" + detectorName[i]).c_str(), ("SiPM_Mulitplicity_Cleaned_" + detectorName[i]).c_str(), 20, 0, 20);
      H_SiPM_Mulitplicity_Cleaned[i]->GetXaxis()->SetTitle("Multiplicity");
      H_SiPM_Mulitplicity_Cleaned[i]->GetYaxis()->SetTitle("Counts");
      H_SiPM_Mulitplicity_Cleaned[i]->GetXaxis()->CenterTitle();
      H_SiPM_Mulitplicity_Cleaned[i]->GetYaxis()->CenterTitle();
    }

    
  }

  // for each coinc time and multiplicity HIGH
  for (int index_time = 0; index_time < CoincidenceTime.size(); index_time++)
  {
    for (int index_time1 = 0; index_time1 < CoincidenceTime.size(); index_time1++)
    {
      for (int mul = 0; mul < SIGNAL_MAX; mul++)
      {
        H_HighSiPM_Time_IAS_Cleaned[index_time][index_time1][mul] = new TH1D(("HighSiPM_Time_IAS_Cleaned_" + to_string(CoincidenceTime[index_time]) + "ns_" + to_string(CoincidenceTime[index_time1]) + "ns_M" + to_string(mul)).c_str(), ("HighSiPM_Time_IAS_Cleaned_" + to_string(CoincidenceTime[index_time]) + "ns_" + to_string(CoincidenceTime[index_time1]) + "ns_M" + to_string(mul)).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
        H_HighSiPM_Time_IAS_Cleaned[index_time][index_time1][mul]->GetXaxis()->SetTitle("Rear-SiPM Time (ns)");
        H_HighSiPM_Time_IAS_Cleaned[index_time][index_time1][mul]->GetYaxis()->SetTitle("Counts");
        H_HighSiPM_Time_IAS_Cleaned[index_time][index_time1][mul]->GetXaxis()->CenterTitle();
        H_HighSiPM_Time_IAS_Cleaned[index_time][index_time1][mul]->GetYaxis()->CenterTitle();

        H_HighSiPM_Time_IAS_UnCleaned[index_time][index_time1][mul] = new TH1D(("HighSiPM_Time_IAS_UnCleaned_" + to_string(CoincidenceTime[index_time]) + "ns_" + to_string(CoincidenceTime[index_time1]) + "ns_M" + to_string(mul)).c_str(), ("HighSiPM_Time_IAS_UnCleaned_" + to_string(CoincidenceTime[index_time]) + "ns_" + to_string(CoincidenceTime[index_time1]) + "ns_M" + to_string(mul)).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
        H_HighSiPM_Time_IAS_UnCleaned[index_time][index_time1][mul]->GetXaxis()->SetTitle("Rear-SiPM Time (ns)");
        H_HighSiPM_Time_IAS_UnCleaned[index_time][index_time1][mul]->GetYaxis()->SetTitle("Counts");
        H_HighSiPM_Time_IAS_UnCleaned[index_time][index_time1][mul]->GetXaxis()->CenterTitle();
        H_HighSiPM_Time_IAS_UnCleaned[index_time][index_time1][mul]->GetYaxis()->CenterTitle();

        H_Energy_IAS_Cleaned[index_time][index_time1][mul] = new TH1D(("Energy_IAS_Cleaned_" + to_string(CoincidenceTime[index_time]) + "ns_" + to_string(CoincidenceTime[index_time1]) + "ns_M" + to_string(mul)).c_str(), ("Energy_IAS_Cleaned_" + to_string(CoincidenceTime[index_time]) + "ns_" + to_string(CoincidenceTime[index_time1]) + "ns_M" + to_string(mul)).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
        H_Energy_IAS_Cleaned[index_time][index_time1][mul]->GetXaxis()->SetTitle("Energy (keV)");
        H_Energy_IAS_Cleaned[index_time][index_time1][mul]->GetYaxis()->SetTitle("Counts");
        H_Energy_IAS_Cleaned[index_time][index_time1][mul]->GetXaxis()->CenterTitle();
        H_Energy_IAS_Cleaned[index_time][index_time1][mul]->GetYaxis()->CenterTitle();

        H_Energy_IAS_UnCleaned[index_time][index_time1][mul] = new TH1D(("Energy_IAS_UnCleaned_" + to_string(CoincidenceTime[index_time]) + "ns_" + to_string(CoincidenceTime[index_time1]) + "ns_M" + to_string(mul)).c_str(), ("Energy_IAS_UnCleaned_" + to_string(CoincidenceTime[index_time]) + "ns_" + to_string(CoincidenceTime[index_time1]) + "ns_M" + to_string(mul)).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
        H_Energy_IAS_UnCleaned[index_time][index_time1][mul]->GetXaxis()->SetTitle("Energy (keV)");
        H_Energy_IAS_UnCleaned[index_time][index_time1][mul]->GetYaxis()->SetTitle("Counts");
        H_Energy_IAS_UnCleaned[index_time][index_time1][mul]->GetXaxis()->CenterTitle();
        H_Energy_IAS_UnCleaned[index_time][index_time1][mul]->GetYaxis()->CenterTitle();
      }
    }
  }

  // for each coinc time and multiplicity LOW
  for (int index_time = 0; index_time < CoincidenceTime.size(); index_time++)
  {
    for (int index_time1 = 0; index_time1 < CoincidenceTime.size(); index_time1++)
    {
      for (int mul = 0; mul < SIGNAL_MAX; mul++)
      {
        H_LowSiPM_Time_IAS_Cleaned[index_time][index_time1][mul] = new TH1D(("LowSiPM_Time_IAS_Cleaned_" + to_string(CoincidenceTime[index_time]) + "ns_" + to_string(CoincidenceTime[index_time1]) + "ns_M" + to_string(mul)).c_str(), ("LowSiPM_Time_IAS_Cleaned_" + to_string(CoincidenceTime[index_time]) + "ns_" + to_string(CoincidenceTime[index_time1]) + "ns_M" + to_string(mul)).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
        H_LowSiPM_Time_IAS_Cleaned[index_time][index_time1][mul]->GetXaxis()->SetTitle("Rear-SiPM Time (ns)");
        H_LowSiPM_Time_IAS_Cleaned[index_time][index_time1][mul]->GetYaxis()->SetTitle("Counts");
        H_LowSiPM_Time_IAS_Cleaned[index_time][index_time1][mul]->GetXaxis()->CenterTitle();
        H_LowSiPM_Time_IAS_Cleaned[index_time][index_time1][mul]->GetYaxis()->CenterTitle();

        H_LowSiPM_Time_IAS_UnCleaned[index_time][index_time1][mul] = new TH1D(("LowSiPM_Time_IAS_UnCleaned_" + to_string(CoincidenceTime[index_time]) + "ns_" + to_string(CoincidenceTime[index_time1]) + "ns_M" + to_string(mul)).c_str(), ("LowSiPM_Time_IAS_UnCleaned_" + to_string(CoincidenceTime[index_time]) + "ns_" + to_string(CoincidenceTime[index_time1]) + "ns_M" + to_string(mul)).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
        H_LowSiPM_Time_IAS_UnCleaned[index_time][index_time1][mul]->GetXaxis()->SetTitle("Rear-SiPM Time (ns)");
        H_LowSiPM_Time_IAS_UnCleaned[index_time][index_time1][mul]->GetYaxis()->SetTitle("Counts");
        H_LowSiPM_Time_IAS_UnCleaned[index_time][index_time1][mul]->GetXaxis()->CenterTitle();
        H_LowSiPM_Time_IAS_UnCleaned[index_time][index_time1][mul]->GetYaxis()->CenterTitle();
      }
    }
  }
    

  for (int mul = 0; mul < SIGNAL_MAX; mul++)
    {
      
      H_HighSiPM_FakeEvent_Multiplicity[mul] = new TH2D(("HighSiPM_FakeEvent_Multiplicity_M" + to_string(mul)).c_str(), ("HighSiPM_FakeEvent_Multiplicity_M" + to_string(mul)).c_str(), (CoincidenceTime_MAX)/CoincidenceTime_BIN+1, -CoincidenceTime_MAX-CoincidenceTime_BIN/2, 0+CoincidenceTime_BIN/2, (CoincidenceTime_MAX)/CoincidenceTime_BIN+1, 0-CoincidenceTime_BIN/2, CoincidenceTime_MAX+CoincidenceTime_BIN/2);
      H_HighSiPM_FakeEvent_Multiplicity[mul]->GetXaxis()->SetTitle("Rear-SiPM Time Left (ns)");
      H_HighSiPM_FakeEvent_Multiplicity[mul]->GetYaxis()->SetTitle("Rear-SiPM Time Right (ns)");
      H_HighSiPM_FakeEvent_Multiplicity[mul]->GetXaxis()->CenterTitle();
      H_HighSiPM_FakeEvent_Multiplicity[mul]->GetYaxis()->CenterTitle();

      H_LowSiPM_FakeEvent_Multiplicity[mul] = new TH2D(("LowSiPM_FakeEvent_Multiplicity_M" + to_string(mul)).c_str(), ("LowSiPM_FakeEvent_Multiplicity_M" + to_string(mul)).c_str(), (CoincidenceTime_MAX)/CoincidenceTime_BIN+1, -CoincidenceTime_MAX-CoincidenceTime_BIN/2, 0+CoincidenceTime_BIN/2, (CoincidenceTime_MAX)/CoincidenceTime_BIN+1, 0-CoincidenceTime_BIN/2, CoincidenceTime_MAX+CoincidenceTime_BIN/2);
      H_LowSiPM_FakeEvent_Multiplicity[mul]->GetXaxis()->SetTitle("Rear-SiPM Time Left (ns)");
      H_LowSiPM_FakeEvent_Multiplicity[mul]->GetYaxis()->SetTitle("Rear-SiPM Time Right (ns)");
      H_LowSiPM_FakeEvent_Multiplicity[mul]->GetXaxis()->CenterTitle();
      H_LowSiPM_FakeEvent_Multiplicity[mul]->GetYaxis()->CenterTitle();

      H_HighSiPM_Losses_Multiplicity[mul] = new TH2D(("HighSiPM_Losses_Multiplicity_M" + to_string(mul)).c_str(), ("HighSiPM_Losses_Multiplicity_M" + to_string(mul)).c_str(), (CoincidenceTime_MAX)/CoincidenceTime_BIN+1, -CoincidenceTime_MAX-CoincidenceTime_BIN/2, 0+CoincidenceTime_BIN/2, (CoincidenceTime_MAX)/CoincidenceTime_BIN+1, 0-CoincidenceTime_BIN/2, CoincidenceTime_MAX+CoincidenceTime_BIN/2);
      H_HighSiPM_Losses_Multiplicity[mul]->GetXaxis()->SetTitle("Rear-SiPM Time Left (ns)");
      H_HighSiPM_Losses_Multiplicity[mul]->GetYaxis()->SetTitle("Rear-SiPM Time Right (ns)");
      H_HighSiPM_Losses_Multiplicity[mul]->GetXaxis()->CenterTitle();
      H_HighSiPM_Losses_Multiplicity[mul]->GetYaxis()->CenterTitle();

      H_LowSiPM_Losses_Multiplicity[mul] = new TH2D(("LowSiPM_Losses_Multiplicity_M" + to_string(mul)).c_str(), ("LowSiPM_Losses_Multiplicity_M" + to_string(mul)).c_str(), (CoincidenceTime_MAX)/CoincidenceTime_BIN+1, -CoincidenceTime_MAX-CoincidenceTime_BIN/2, 0+CoincidenceTime_BIN/2, (CoincidenceTime_MAX)/CoincidenceTime_BIN+1, 0-CoincidenceTime_BIN/2, CoincidenceTime_MAX+CoincidenceTime_BIN/2);
      H_LowSiPM_Losses_Multiplicity[mul]->GetXaxis()->SetTitle("Rear-SiPM Time Left (ns)");
      H_LowSiPM_Losses_Multiplicity[mul]->GetYaxis()->SetTitle("Rear-SiPM Time Right (ns)");
      H_LowSiPM_Losses_Multiplicity[mul]->GetXaxis()->CenterTitle();
      H_LowSiPM_Losses_Multiplicity[mul]->GetYaxis()->CenterTitle();

      H_2SiPM_Time_MeanFurthest[mul] = new TH1D(("2SiPM_Time_MeanFurthest_M" + to_string(mul)).c_str(), ("2SiPM_Time_MeanFurthest_M" + to_string(mul)).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
      H_2SiPM_Time_MeanFurthest[mul]->GetXaxis()->SetTitle("Time between 2 SiPMs (ns)");
      H_2SiPM_Time_MeanFurthest[mul]->GetYaxis()->SetTitle("Counts");
      H_2SiPM_Time_MeanFurthest[mul]->GetXaxis()->CenterTitle();
      H_2SiPM_Time_MeanFurthest[mul]->GetYaxis()->CenterTitle();

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

  //////////////////////////////////////////////////////
  H_Energy_IAS = new TH1D("Energy_IAS", "Energy_IAS", eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
  H_Energy_IAS->GetXaxis()->SetTitle("Energy (keV)");
  H_Energy_IAS->GetYaxis()->SetTitle("Counts");
  H_Energy_IAS->GetXaxis()->CenterTitle();
  H_Energy_IAS->GetYaxis()->CenterTitle();


  // new 
  H_RearSiPM_MeanTime_Nearest_Group = new TH1D("RearSiPM_MeanTime_Nearest_Group", "RearSiPM_MeanTime_Nearest_Group", winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
  H_RearSiPM_MeanTime_Nearest_Group->GetXaxis()->SetTitle("Rear-SiPM Time (ns)");
  H_RearSiPM_MeanTime_Nearest_Group->GetYaxis()->SetTitle("Counts");
  H_RearSiPM_MeanTime_Nearest_Group->GetXaxis()->CenterTitle();
  H_RearSiPM_MeanTime_Nearest_Group->GetYaxis()->CenterTitle();

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
      int status = G_RearStrip_Spread[i]->Fit(SpreadAcceptance[i], "Q");

      if (status != 0)
        Warning("Fit failed for " + detectorName[i] + " during Rear/Strip spread correction.");

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
      MeanAcceptance_Walk_Silicon[i]->SetParameters(150, -1e6, -1000);
      MeanAcceptance_Walk_Silicon[i]->SetParLimits(0, 20, 150);
      MeanAcceptance_Walk_Silicon[i]->SetParLimits(1, -1e7, -1e5);
      MeanAcceptance_Walk_Silicon[i]->SetParLimits(2, -150000, 0);

      int status = G_Rear_Mean_Walk_Silicon[i]->Fit(MeanAcceptance_Walk_Silicon[i], "Q");
      if (status != 0)
        Warning("Fit failed for " + detectorName[i] + " during walk Silicon correction.");

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

      int status = graph_mean->Fit(MeanAcceptance_Walk_SiPM[i], "WQR");
      if (status != 0)
        Warning("Fit failed for " + detectorName[i] + " during walk SiPM correction.");

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

      // 4th
      //  C_SiPM_ChannelTime_Walk_SiPM[i]->cd();
      //  TPad *P_SiPM_ChannelTime_Walk_SiPM_proj = new TPad("P_SiPM_ChannelTime_Walk_SiPM_proj", "P_SiPM_ChannelTime_Walk_SiPM_proj", 0.55, 0.05, 1., 0.5);
      //  P_SiPM_ChannelTime_Walk_SiPM_proj->Draw();
      //  P_SiPM_ChannelTime_Walk_SiPM_proj->cd();
      //  P_SiPM_ChannelTime_Walk_SiPM_proj->SetLogy();
      //  TH1D *proj = (TH1D *)H_SiPM_ChannelTime_Walk_Corrected[i]->ProjectionX(("H_SiPM_ChannelTime_Walk_SiPM_projx" + detectorName[i]).c_str());
      //  proj->Draw("HIST");

      // proj->GetXaxis()->SetRangeUser(-200, 200);
      // TH1D *BKG = (TH1D*)proj->ShowBackground(20, "SAME");

      // TH1D *Distribution = (TH1D*)proj->Clone(("Distribution_" + detectorName[i]).c_str());
      // Distribution->Add(BKG, -1);

      // LIMIT CALCULATION 3-SIGMA
      double range[2];
      // double total = 0;
      // for (int bin = 0; bin < Distribution->GetNbinsX(); bin++)
      // {
      //   if (Distribution->GetBinCenter(bin) > -200 && Distribution->GetBinCenter(bin) < 200)
      //   {
      //     total += Distribution->GetBinContent(bin);
      //   }
      // }
      // double first = total * 0.00135;
      // double last = total * 0.99865;
      // double integral = 0;
      // for (int bin = 0; bin < Distribution->GetNbinsX(); bin++)
      // {
      //   if (Distribution->GetBinCenter(bin) > -200 && Distribution->GetBinCenter(bin) < 200)
      //   {
      //     if ((integral < first) && ((integral + Distribution->GetBinContent(bin)) > first))
      //     {
      //       range[0] = Distribution->GetBinCenter(bin);
      //     }
      //     else if ((integral < last) && ((integral + Distribution->GetBinContent(bin)) > last))
      //     {
      //       range[1] = Distribution->GetBinCenter(bin);
      //     }
      //     integral += Distribution->GetBinContent(bin);
      //   }
      // }
      /////////////////////////////////

      /////// MANUAL CUTS ////////
      range[0] = -30;
      range[1] = 20;
      ////////////////////////////

      // DISPLAY SELECTION
      // double integral_tot = 0;
      // double integral_BKG = 0;
      // for (int bin = 0; bin < proj->GetNbinsX(); bin++)
      // {
      //   if (proj->GetBinCenter(bin) < range[0] || proj->GetBinCenter(bin) > range[1])
      //   {
      //     proj->SetBinContent(bin, 0);
      //   }
      //   if (proj->GetBinCenter(bin) > range[0] && proj->GetBinCenter(bin) < range[1])
      //   {
      //     integral_tot += proj->GetBinContent(bin);
      //     integral_BKG += BKG->GetBinContent(bin);
      //   }
      // }
      // TH1D *H_Selection = (TH1D *)proj->DrawCopy("SAME");
      // H_Selection->GetXaxis()->SetRangeUser(range[0], range[1]);
      // H_Selection->SetFillColor(kRed);
      // H_Selection->SetFillStyle(3244);

      // Saving for cuts
      // if (IsDetectorBetaLow(i))
      SpreadAcceptance_Walk_SiPM[GetDetectorChannel(i)] = make_pair(range[0], range[1]);

      // cout << " Signal/Noise: " << setprecision(10) << (integral_tot-integral_BKG)/integral_tot << "        Ranges:" << range[0] << "  " << range[1] << endl;

      C_SiPM_ChannelTime_Walk_SiPM[i]->Write();
    }
  }
}

inline void WriteHistograms_Cleaned()
{
  GROUPED_File->cd();
  H_RearSiPM_MeanTime_Nearest_Group->Write();
  for (int i = 0; i < detectorNum; i++)
  {
    ProgressCounter(i, detectorNum, " <INFO>  Writing Histograms");
    if (IsDetectorSiliStrip(i))
    {
      dir_Strip_Channel_Detector[i]->cd();
      H_Channel_Cleaned[i]->Write();
      H_RearMatched_Channel_Cleaned[i]->Write();
      H_RearMatchedMean_Channel_Cleaned[i]->Write();

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
      dir_Rear_Channel_Detector[i]->cd();
      H_Channel_Cleaned[i]->Write();
      H_Rear_Channel_IAS_Cleaned[i]->Write();
      H_Rear_Channel_IAScoinc_Cleaned[i]->Write();
      dir_Rear_Time_Detector[i]->cd();
      H_RearSiPM_ChannelTime_Walk_Corrected[i]->Write();
      H_RearSiPM_ChannelTime_Cleaned[i]->Write();
      H_RearSiPM_Time_Nearest[i]->Write();

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
      dir_SiPMHigh_Time_Detector[i]->cd();
      H_SiPM_ChannelTime_Walk_Corrected[i]->Write();
      H_SiPM_ChannelTime_Cleaned[i]->Write();
      H_SiPM_Mulitplicity_Cleaned[i]->Write();
      H_2SiPM_Channel_Cleaned[i]->Write();
      H_2SiPM_Time_IAS_Cleaned[i]->Write();
      
      H_SiPMLowHigh_Time_One[i]->Write();
      H_SiPMLowHigh_TimeChannel_One[i]->Write();
      H_SiPMLowHigh_Time_Multiple[i]->Write();

      for (int j = i; j < SIGNAL_MAX; j++)
      {
        if (IsDetectorBetaHigh(j))
        {
          H_2SiPM_Time_One[i][j]->Write();
          H_2SiPM_Time_Multiple[i][j]->Write();
          H_RearSiPM_Time_Nearest_Nearest[i][j]->Write();
          H_RearSiPM_Channel_Nearset_Nearest[i][j]->Write();
          H_RearSiPM_Channel_Nearest_Group[i][j]->Write();
          H_2SiPM_Time_Next[i][j]->Write();
        }
      }

    }

    if (IsDetectorBetaLow(i))
    {
      dir_SiPMLow_Time_Detector[i]->cd();
      H_SiPM_ChannelTime_Walk_Corrected[i]->Write();
      H_SiPM_ChannelTime_Cleaned[i]->Write();
      H_SiPM_Mulitplicity_Cleaned[i]->Write();
    }
  }

  //// TOTAL AS REF FOR LOSSES /////
  int vector_size = CoincidenceTime.size() - 1;
  double total_integral_maxtimeHigh[BETA_SIZE + 1];
  double total_integral_maxtimeLow[BETA_SIZE + 1];
  for (int mul = 1; mul <= BETA_SIZE; mul++)
  {
    // H_HighSiPM_Time_IAS_UnCleaned[vector_size][vector_size][mul]->GetXaxis()->SetRangeUser(-200, 300);
    // TH1D *BKG = (TH1D *)H_HighSiPM_Time_IAS_UnCleaned[vector_size][vector_size][mul]->ShowBackground(20);
    //// COMPUTING SIGNAL/NOISE
    double integral_tot = H_HighSiPM_Time_IAS_Cleaned[vector_size][vector_size][mul]->Integral();
    double integral_rl_window = H_HighSiPM_Time_IAS_UnCleaned[vector_size][vector_size][mul]->Integral()/2;
    // for (int bin = 0; bin < H_HighSiPM_Time_IAS_Cleaned[vector_size][vector_size][mul]->GetNbinsX(); bin++)
    // {
    //   if (H_HighSiPM_Time_IAS_Cleaned[vector_size][vector_size][mul]->GetBinCenter(bin) > -CoincidenceTime[vector_size] && H_HighSiPM_Time_IAS_Cleaned[vector_size][vector_size][mul]->GetBinCenter(bin) < CoincidenceTime[vector_size])
    //   {
    //     integral_tot += H_HighSiPM_Time_IAS_Cleaned[vector_size][vector_size][mul]->GetBinContent(bin);
    //     integral_BKG += BKG->GetBinContent(bin);
    //   }
    // }
    total_integral_maxtimeHigh[mul] = integral_tot-integral_rl_window;

    // H_LowSiPM_Time_IAS_UnCleaned[vector_size][vector_size][mul]->GetXaxis()->SetRangeUser(-200, 300);
    // TH1D *BKG_Low = (TH1D *)H_LowSiPM_Time_IAS_UnCleaned[vector_size][vector_size][mul]->ShowBackground(20);
    //// COMPUTING SIGNAL/NOISE
    double integral_tot_Low = H_LowSiPM_Time_IAS_Cleaned[vector_size][vector_size][mul]->Integral();
    double integral_rl_window_Low = H_LowSiPM_Time_IAS_UnCleaned[vector_size][vector_size][mul]->Integral()/2;
    // for (int bin = 0; bin < H_LowSiPM_Time_IAS_UnCleaned[vector_size][vector_size][mul]->GetNbinsX(); bin++)
    // {
    //   if (H_LowSiPM_Time_IAS_Cleaned[vector_size][vector_size][mul]->GetBinCenter(bin) > -CoincidenceTime[vector_size] && H_LowSiPM_Time_IAS_Cleaned[vector_size][vector_size][mul]->GetBinCenter(bin) < CoincidenceTime[vector_size])
    //   {
    //     integral_tot_Low += H_LowSiPM_Time_IAS_Cleaned[vector_size][vector_size][mul]->GetBinContent(bin);
    //     integral_BKG_Low += BKG_Low->GetBinContent(bin);
    //   }
    // }
    total_integral_maxtimeLow[mul] = integral_tot_Low - integral_rl_window_Low;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////
  for (int Time_index_l = 0; Time_index_l < CoincidenceTime.size(); Time_index_l++)
    {
      for (int Time_index_r = 0; Time_index_r < CoincidenceTime.size(); Time_index_r++)
    {
      for (int mul = 1; mul <= BETA_SIZE; mul++)
      {
        // HIGH
        dir_Cleaned_HighSiPM_TimeCoincidence_Mulitplicity[mul]->cd();
        // H_HighSiPM_Time_IAS_Cleaned[Time_index][mul]->Write();
        // H_HighSiPM_Time_IAS_UnCleaned[Time_index][mul]->Write();

        dir_Cleaned_HighSiPM_TimeCoincidence_Mulitplicity[mul]->cd();
        TCanvas *C_IAS_Time_Cleaned_High = new TCanvas(("C_IAS_Time_Cleaned_High_" + to_string(-CoincidenceTime[Time_index_l]) + "ns_" + to_string(CoincidenceTime[Time_index_r]) +"ns_M" + to_string(mul)).c_str(), ("C_IAS_Time_Cleaned_" + to_string(CoincidenceTime[Time_index_l]) + "ns_" + to_string(CoincidenceTime[Time_index_r]) +"ns_M" + to_string(mul)).c_str(), 800, 400);
        
        C_IAS_Time_Cleaned_High->Divide(3, 1);
        C_IAS_Time_Cleaned_High->cd(1);
        C_IAS_Time_Cleaned_High->SetLogy();
        H_HighSiPM_Time_IAS_Cleaned[Time_index_l][Time_index_r][mul]->GetYaxis()->SetRangeUser(1, -1111);
        H_HighSiPM_Time_IAS_Cleaned[Time_index_l][Time_index_r][mul]->SetFillStyle(3244);
        H_HighSiPM_Time_IAS_Cleaned[Time_index_l][Time_index_r][mul]->SetFillColor(kRed);
        H_HighSiPM_Time_IAS_Cleaned[Time_index_l][Time_index_r][mul]->Draw("HIST");
        H_HighSiPM_Time_IAS_UnCleaned[Time_index_l][Time_index_r][mul]->SetFillColor(3244);
        H_HighSiPM_Time_IAS_UnCleaned[Time_index_l][Time_index_r][mul]->SetLineColor(kBlack);
        H_HighSiPM_Time_IAS_UnCleaned[Time_index_l][Time_index_r][mul]->Draw("SAME");
        

        // H_HighSiPM_Time_IAS_UnCleaned[Time_index_l][Time_index_r][mul]->GetXaxis()->SetRangeUser(-200, 300);
        // TH1D *BKG = (TH1D *)H_HighSiPM_Time_IAS_UnCleaned[Time_index_l][Time_index_r][mul]->ShowBackground(20);

        //// COMPUTING SIGNAL/NOISE
        double integral_tot = H_HighSiPM_Time_IAS_Cleaned[Time_index_l][Time_index_r][mul]->Integral();
        double integral_rl_window = H_HighSiPM_Time_IAS_UnCleaned[Time_index_l][Time_index_r][mul]->Integral()/2;
        
        TLatex *text = new TLatex();
        text->SetNDC();
        text->SetTextSize(0.03);
        text->DrawLatex(0.2, 0.8, ("Fake events: " + to_string(integral_rl_window * 100 / integral_tot) + " %").c_str());
        text->Draw("SAME");

        TLatex *text_loss = new TLatex();
        text_loss->SetNDC();
        text_loss->SetTextSize(0.03);
        text_loss->DrawLatex(0.2, 0.7, ("Losses: " + to_string((total_integral_maxtimeHigh[mul]-(integral_tot-integral_rl_window)) * 100 / total_integral_maxtimeHigh[mul]) + " %").c_str());
        text_loss->Draw("SAME");

        H_HighSiPM_FakeEvent_Multiplicity[mul]->Fill(-CoincidenceTime[Time_index_l], CoincidenceTime[Time_index_r], integral_rl_window * 100 / integral_tot);
        H_HighSiPM_Losses_Multiplicity[mul]->Fill(-CoincidenceTime[Time_index_l], CoincidenceTime[Time_index_r], (total_integral_maxtimeHigh[mul]-(integral_tot-integral_rl_window)) * 100 / total_integral_maxtimeHigh[mul]);
        
        C_IAS_Time_Cleaned_High->cd(2);
        C_IAS_Time_Cleaned_High->SetLogy();
        H_Energy_IAS->GetXaxis()->SetRangeUser(3300, 3400);
        H_Energy_IAS->GetYaxis()->SetRangeUser(1, -1111);
        H_Energy_IAS->SetLineColor(kBlue);
        H_Energy_IAS->Draw("HIST");
        H_Energy_IAS_Cleaned[Time_index_l][Time_index_r][mul]->SetLineColor(kRed);
        H_Energy_IAS_Cleaned[Time_index_l][Time_index_r][mul]->Draw("HIST SAME");
        H_Energy_IAS_UnCleaned[Time_index_l][Time_index_r][mul]->SetLineColor(kBlack);
        H_Energy_IAS_UnCleaned[Time_index_l][Time_index_r][mul]->Draw("HIST SAME");
        TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);
        legend->AddEntry(H_Energy_IAS, "Single", "l");
        legend->AddEntry(H_Energy_IAS_Cleaned[Time_index_l][Time_index_r][mul], "Coincidence", "l");
        legend->AddEntry(H_Energy_IAS_UnCleaned[Time_index_l][Time_index_r][mul], "Excluded", "l");
        legend->Draw("SAME");

        C_IAS_Time_Cleaned_High->cd(3);
        C_IAS_Time_Cleaned_High->SetLogy();
        TH1D* sub = (TH1D*)H_Energy_IAS->Clone(("sub_" + to_string(Time_index_l) + "_" + to_string(Time_index_r) + "_" + to_string(mul)).c_str());
        sub->GetXaxis()->SetRangeUser(3300, 3400);
        sub->GetYaxis()->SetRangeUser(1, -1111);
        sub->Add(H_Energy_IAS_Cleaned[Time_index_l][Time_index_r][mul], -1);
        sub->SetLineColor(kGreen);
        sub->Scale(H_Energy_IAS_UnCleaned[Time_index_l][Time_index_r][mul]->Integral()/sub->Integral());
        sub->Draw("HIST");
        H_Energy_IAS_UnCleaned[Time_index_l][Time_index_r][mul]->Draw("HIST SAME");
               

        C_IAS_Time_Cleaned_High->Write();

        // LOW
        dir_Cleaned_LowSiPM_TimeCoincidence_Mulitplicity[mul]->cd();
        // H_LowSiPM_Time_IAS_Cleaned[Time_index][mul]->Write();
        // H_LowSiPM_Time_IAS_UnCleaned[Time_index][mul]->Write();

        dir_Cleaned_LowSiPM_TimeCoincidence_Mulitplicity[mul]->cd();
        TCanvas *C_IAS_Time_Cleaned_low = new TCanvas(("C_IAS_Time_Cleaned_Low_" + to_string(-CoincidenceTime[Time_index_l]) + "ns_" + to_string(CoincidenceTime[Time_index_r]) +"ns_M" + to_string(mul)).c_str(), ("C_IAS_Time_Cleaned_" + to_string(CoincidenceTime[Time_index_l]) + "ns_" + to_string(CoincidenceTime[Time_index_r]) +"ns_M" + to_string(mul)).c_str(), 800, 400);
        C_IAS_Time_Cleaned_low->cd();
        C_IAS_Time_Cleaned_low->SetLogy();
        H_LowSiPM_Time_IAS_Cleaned[Time_index_l][Time_index_r][mul]->SetFillStyle(3244);
        H_LowSiPM_Time_IAS_Cleaned[Time_index_l][Time_index_r][mul]->SetFillColor(kRed);
        H_LowSiPM_Time_IAS_Cleaned[Time_index_l][Time_index_r][mul]->Draw("HIST");
        H_LowSiPM_Time_IAS_UnCleaned[Time_index_l][Time_index_r][mul]->SetFillColor(3244);
       
        H_LowSiPM_Time_IAS_UnCleaned[Time_index_l][Time_index_r][mul]->SetLineColor(kBlack);
        H_LowSiPM_Time_IAS_UnCleaned[Time_index_l][Time_index_r][mul]->Draw("SAME");
        

        // H_LowSiPM_Time_IAS_UnCleaned[Time_index_l][Time_index_r][mul]->GetXaxis()->SetRangeUser(-200, 300);
        // TH1D *BKG_low = (TH1D *)H_LowSiPM_Time_IAS_UnCleaned[Time_index_l][Time_index_r][mul]->ShowBackground(20);

        //// COMPUTING SIGNAL/NOISE
        double integral_tot_low = H_LowSiPM_Time_IAS_Cleaned[Time_index_l][Time_index_r][mul]->Integral();
        double integral_rl_window_Low = H_LowSiPM_Time_IAS_UnCleaned[Time_index_l][Time_index_r][mul]->Integral()/2;
        // for (int bin = 0; bin < H_LowSiPM_Time_IAS_UnCleaned[Time_index_l][Time_index_r][mul]->GetNbinsX(); bin++)
        // {
        //   if (H_LowSiPM_Time_IAS_Cleaned[Time_index_l][Time_index_r][mul]->GetBinCenter(bin) > -CoincidenceTime[Time_index_l] && H_LowSiPM_Time_IAS_Cleaned[Time_index_l][Time_index_r][mul]->GetBinCenter(bin) < CoincidenceTime[Time_index_r])
        //   {
        //     integral_tot_low += H_LowSiPM_Time_IAS_Cleaned[Time_index_l][Time_index_r][mul]->GetBinContent(bin);
        //     integral_BKG_low += BKG_low->GetBinContent(bin);
        //   }
        // }

        TLatex *text_low = new TLatex();
        text_low->SetNDC();
        text_low->SetTextSize(0.03);
        text_low->DrawLatex(0.2, 0.8, ("Fake events: " + to_string(integral_rl_window_Low * 100 / integral_tot_low) + " %").c_str());
        text_low->Draw("SAME");

        TLatex *text_low_loss = new TLatex();
        text_low_loss->SetNDC();
        text_low_loss->SetTextSize(0.03);
        text_low_loss->DrawLatex(0.2, 0.7, ("Losses: " + to_string((total_integral_maxtimeLow[mul]-(integral_tot_low-integral_rl_window_Low)) * 100 / total_integral_maxtimeLow[mul]) + " %").c_str());
        text_low_loss->Draw("SAME");

        H_LowSiPM_FakeEvent_Multiplicity[mul]->Fill(-CoincidenceTime[Time_index_l], CoincidenceTime[Time_index_r], integral_rl_window_Low * 100 / integral_tot_low);
        H_LowSiPM_Losses_Multiplicity[mul]->Fill(-CoincidenceTime[Time_index_l], CoincidenceTime[Time_index_r], (total_integral_maxtimeLow[mul]-(integral_tot_low-integral_rl_window_Low)) * 100 / total_integral_maxtimeLow[mul]);
        C_IAS_Time_Cleaned_low->Write();
      }
    }
    }

  dir_Cleaned->cd();

  //HIGH
  for (int mul = 1; mul <= MAX_MULTIPLICTY; mul++)
  {
    TCanvas *c = new TCanvas(("C_HighSiPM_FakeMultiplicity_Cleaned_M"+to_string(mul)).c_str(), ("C_HighSiPM_FakeMultiplicity_Cleaned_M"+to_string(mul)).c_str(), 800, 400);
    H_HighSiPM_FakeEvent_Multiplicity[mul]->Draw("COLZ");
    c->Write();

    TCanvas *c1 = new TCanvas(("C_HighSiPM_Losses_Multiplicity_Cleaned_M"+to_string(mul)).c_str(), ("C_HighSiPM_Losses_Multiplicity_Cleaned_M"+to_string(mul)).c_str(), 800, 400);
    H_HighSiPM_Losses_Multiplicity[mul]->Draw("COLZ");
    c1->Write();
  }

  //LOW
  for (int mul = 1; mul <= MAX_MULTIPLICTY; mul++)
  {
    TCanvas *c = new TCanvas(("C_LowSiPM_FakeMultiplicity_Cleaned_M"+to_string(mul)).c_str(), ("C_LowSiPM_FakeMultiplicity_Cleaned_M"+to_string(mul)).c_str(), 800, 400);
    H_LowSiPM_FakeEvent_Multiplicity[mul]->Draw("COLZ");
    c->Write();

    TCanvas *c1 = new TCanvas(("C_LowSiPM_Losses_Multiplicity_Cleaned_M"+to_string(mul)).c_str(), ("C_LowSiPM_Losses_Multiplicity_Cleaned_M"+to_string(mul)).c_str(), 800, 400);
    H_LowSiPM_Losses_Multiplicity[mul]->Draw("COLZ");
    c1->Write();
  }


  dir_SiPMHigh_Multiplicity->cd();
  H_SiPMHigh_Mulitplicity_Cleaned->Write();
  dir_SiPMLow_Multiplicity->cd();
  H_SiPMLow_Mulitplicity_Cleaned->Write();

  dir_Cleaned->cd();
  C_IAS_Channel_Cleaned->Write();

  /// GENERAL FITS CANVAS ///
  dir_Fits->cd();
  C_Strip_Cutting_Fits->Write();
  C_Rear_Walk_Silicon_Fits->Write();
  C_SiPM_Walk_SiPM_Fits->Write();

  GROUPED_File->cd();
  for (int mul = 1; mul <= BETA_SIZE; mul++)
  {
    H_2SiPM_Time_MeanFurthest[mul]->Write();
    H_RearSiPM_Time_New_Nearest[mul]->Write();
    H_RearSiPM_Time_New_Mean[mul]->Write();
  }
}

inline void WriteIASLosses()
{

  GROUPED_File->cd();
  /// FOR SILICON STRIPS ///
  TCanvas *C_Losses = new TCanvas("C_Losses", "C_Losses", 800, 400);
  TH1D* H_Losses_all = new TH1D("H_Losses_all", "H_Losses_all", 6, 0, 3);
  TH1D* H_Losses_IAS = new TH1D("H_Losses_IAS", "H_Losses_IAS", 6, 0, 3);

  for (int i = 0; i < detectorNum; i++)
  {
    if (IsDetectorSiliStrip(i))
    {
      H_Losses_all->SetBinContent(0, H_Losses_all->GetBinContent(0) + H_Channel_RAW[i]->Integral()); 
      H_Losses_all->SetBinContent(2, H_Losses_all->GetBinContent(2) + H_Channel_A[i]->Integral());   
      H_Losses_all->SetBinContent(5, H_Losses_all->GetBinContent(5) + H_Channel_Cleaned[i]->Integral());  
      
      H_Channel_RAW[i]->GetXaxis()->SetRangeUser(peaks_window_F[i].first, peaks_window_F[i].second);
      H_Channel_A[i]->GetXaxis()->SetRangeUser(peaks_window_F[i].first, peaks_window_F[i].second);
      H_Channel_Cleaned[i]->GetXaxis()->SetRangeUser(peaks_window_F[i].first, peaks_window_F[i].second);

      H_Losses_IAS->SetBinContent(0, H_Losses_IAS->GetBinContent(0) + H_Channel_RAW[i]->Integral());
      H_Losses_IAS->SetBinContent(2, H_Losses_IAS->GetBinContent(2) + H_Channel_A[i]->Integral());
      H_Losses_IAS->SetBinContent(5, H_Losses_IAS->GetBinContent(5) + H_Channel_Cleaned[i]->Integral());

      H_Channel_RAW[i]->GetXaxis()->SetRangeUser(-1111, -1111);
      H_Channel_A[i]->GetXaxis()->SetRangeUser(-1111, -1111);
      H_Channel_Cleaned[i]->GetXaxis()->SetRangeUser(-1111, -1111);
    }
  }
  H_Losses_all->SetLineColor(kBlack);
  // H_Losses_all->SetBinLabel(1, "RAW");
  // H_Losses_all->SetBinLabel(3, "SELECTED");
  // H_Losses_all->SetBinLabel(6, "CLEANED");
  H_Losses_all->Draw("HIST");
  H_Losses_IAS->SetLineColor(kRed);
  H_Losses_IAS->Draw("SAME");
  C_Losses->Write();


  // /// FOR SIPM
  // TCanvas *C_Losses_SiPM = new TCanvas("C_Losses_SiPM", "C_Losses_SiPM", 800, 400);
  // TH1D* H_Losses_all_SiPM = new TH1D("H_Losses_all_SiPM", "H_Losses_all_SiPM", 6, 0, 3);
  // TH1D* H_Losses_IAS_SiPM = new TH1D("H_Losses_IAS_SiPM", "H_Losses_IAS_SiPM", 6, 0, 3);

  // for (int i = 0; i < detectorNum; i++)
  // {
  //   if (IsDetectorBetaHigh(i))
  //   {
  //     H_Losses_all_SiPM->SetBinContent(0, H_Losses_all_SiPM->GetBinContent(0) + H_SiPMHigh_Mulitplicity_RAW[i]->Integral()); 
  //     // H_Losses_all_SiPM->SetBinContent(2, H_Losses_all_SiPM->GetBinContent(2) + H_SiPMHigh_Mulitplicity_A[i]->Integral());   
  //     H_Losses_all_SiPM->SetBinContent(5, H_Losses_all_SiPM->GetBinContent(5) + H_SiPM_Mulitplicity_Cleaned[i]->Integral());  
      
  //     H_SiPMHigh_Mulitplicity_RAW[i]->GetXaxis()->SetRangeUser(peaks_window_F[i].first, peaks_window_F[i].second);
  //     // H_SiPMHigh_Mulitplicity_A[i]->GetXaxis()->SetRangeUser(peaks_window_F[i].first, peaks_window_F[i].second);
  //     H_SiPM_Mulitplicity_Cleaned[i]->GetXaxis()->SetRangeUser(peaks_window_F[i].first, peaks_window_F[i].second);

  //     H_Losses_IAS_SiPM->SetBinContent(0, H_Losses_IAS_SiPM->GetBinContent(0) + H_SiPMHigh_Mulitplicity_RAW[i]->Integral());
  //     // H_Losses_IAS_SiPM->SetBinContent(2, H_Losses_IAS_SiPM->GetBinContent(2) + H_SiPMHigh_Mulitplicity_A[i]->Integral());
  //     H_Losses_IAS_SiPM->SetBinContent(5, H_Losses_IAS_SiPM->GetBinContent(5) + H_SiPM_Mulitplicity_Cleaned[i]->Integral());

  //     H_SiPMHigh_Mulitplicity_RAW[i]->GetXaxis()->SetRangeUser(-1111, 1111);
  //     // H_SiPMHigh_Mulitplicity_A[i]->GetXaxis()->SetRangeUser(-1111, 1111);
  //     H_SiPM_Mulitplicity_Cleaned[i]->GetXaxis()->SetRangeUser(-1111, 1111);
  //   }
  // }

  // H_Losses_all_SiPM->SetLineColor(kBlack);
  // // H_Losses_all_SiPM->SetBinLabel(1, "RAW");
  // // H_Losses_all_SiPM->SetBinLabel(3, "SELECTED");
  // // H_Losses_all_SiPM->SetBinLabel(6, "CLEANED");
  // H_Losses_all_SiPM->Draw();
  // H_Losses_IAS_SiPM->SetLineColor(kRed);
  // H_Losses_IAS_SiPM->Draw("SAME");
  // C_Losses_SiPM->Write();
}

inline int WriteTree_Grouped()
{
  GROUPED_File->cd();
  CLEANED_Tree->Write();
  for (int i = 0 ; i< detectorNum; i++)
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
  ///////////////////////////////////////
  ///////// CASES A B C D E F G /////////
  ///////////////////////////////////////

  // Case A : SINGLE Rear and Strip associated (most common case)

  if (RearStrip_ASSOCIATED.size() == 1 && Rear_Position.size() == 1 && Strip_Position.size() == 1)
  {
    // RearStrip_ASSOCIATED[0].first.Channel = 1.0911*RearStrip_ASSOCIATED[0].first.Channel;
    // cout << "CASE A" << endl;
    H_Channel_A[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[0].first.Channel);   // REAR
    H_Channel_A[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].second.Channel); // STRIP

    H_RearStrip_Channel_A[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[0].first.Channel, RearStrip_ASSOCIATED[0].second.Channel);  // REAR-STRIP for rear plot
    H_RearStrip_Channel_A[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].first.Channel, RearStrip_ASSOCIATED[0].second.Channel); // REAR-STRIP for strip plot

    H_Strip_Channel_DiffRear_A[RearStrip_ASSOCIATED[0].second.Label]->Fill((RearStrip_ASSOCIATED[0].first.Channel - RearStrip_ASSOCIATED[0].second.Channel) / RearStrip_ASSOCIATED[0].second.Channel);                                                // STRIP/REAR
    H_Strip_Channel_DiffRearvsStrip_A[RearStrip_ASSOCIATED[0].second.Label]->Fill((RearStrip_ASSOCIATED[0].first.Channel - RearStrip_ASSOCIATED[0].second.Channel) / RearStrip_ASSOCIATED[0].second.Channel, RearStrip_ASSOCIATED[0].second.Channel); // STRIP/REAR vs Strip channel

    H_RearStrip_Time_A[RearStrip_ASSOCIATED[0].first.Label]->Fill(RearStrip_ASSOCIATED[0].first.Time - RearStrip_ASSOCIATED[0].second.Time);
    H_RearStrip_Time_A[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].first.Time - RearStrip_ASSOCIATED[0].second.Time);

    //////////////////////// TREE ////////////////////////
    GROUPED_Tree_Silicon.push_back(RearStrip_ASSOCIATED[0].first);
    GROUPED_Tree_Silicon.push_back(RearStrip_ASSOCIATED[0].second);
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

  // Case B : SAME Rear and DIFFERENT NEIGHBOURG Strip associated (intertrip)
  else if (RearStrip_ASSOCIATED.size() == 2 && IsDetectorSiliInterStrip(RearStrip_ASSOCIATED[0].second.Label, RearStrip_ASSOCIATED[1].second.Label))
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
    H_2Strip_Channel_B[RearStrip_ASSOCIATED[0].second.Label]->Fill(RearStrip_ASSOCIATED[0].second.Channel, RearStrip_ASSOCIATED[1].second.Channel);                                             // STRIP-STRIP
    H_2Strip_Channel_B[RearStrip_ASSOCIATED[1].second.Label]->Fill(RearStrip_ASSOCIATED[1].second.Channel, RearStrip_ASSOCIATED[0].second.Channel);                                             // STRIP-STRIP


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

    //////////////////////// TREE ////////////////////////
    GROUPED_Tree_Silicon.push_back(RearStrip_ASSOCIATED[0].first);
    GROUPED_Tree_Silicon.push_back(RearStrip_ASSOCIATED[0].second);
    if (THIRD_ALPHA_PEAK_FLAG.isValid)
      GROUPED_Tree_Silicon.push_back(THIRD_ALPHA_PEAK_FLAG);
    
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

  if (diff > mean - 3 * spread && diff < mean + 3 * spread && (*Silicon)[1].Pileup == 0 && (*Silicon)[0].Pileup == 0)
  {
    H_Channel_Cutted[Rear_Label]->Fill(Rear_Channel);   // REAR
    H_Channel_Cutted[Strip_Label]->Fill(Strip_Channel); // STRIP

    H_RearStrip_Channel_Cutted[Rear_Label]->Fill(Rear_Channel, Strip_Channel);  // REAR-STRIP for rear plot
    H_RearStrip_Channel_Cutted[Strip_Label]->Fill(Rear_Channel, Strip_Channel); // REAR-STRIP for strip plot
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

    H_SiPMLow_Mulitplicity_Cutted->Fill(SiPM_Low_Size);
    H_SiPMHigh_Mulitplicity_Cutted->Fill(SiPM_High_Size);

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

inline void CleaningGroups()
{

  /// CONTINUOUS VALUE OF CHANNELS ///
  for (auto &signal : *Silicon)
  {
    signal.Channel = signal.Channel - 0.5 + gRandom->Rndm();    
    // signal.Time = signal.Time + (- 0.5 + gRandom->Rndm()) * 8; 
    CLEANED_Tree_Silicon.push_back(signal);
  }

  if ((*Silicon).GetSize() == 2)
    {
      double matched = RearStrip_Matching[(*Silicon)[1].Label].first * (*Silicon)[0].Channel + RearStrip_Matching[(*Silicon)[1].Label].second;
      H_RearMatched_Channel_Cleaned[(*Silicon)[1].Label]->Fill(matched);
      H_RearMatchedMean_Channel_Cleaned[(*Silicon)[1].Label]->Fill(( matched + (*Silicon)[1].Channel) /2);
    }

  for (auto &signal : *SiPM_High)
  {
    signal.Channel = signal.Channel - 0.5 + gRandom->Rndm();
    signal.Time = signal.Time + (- 0.5 + gRandom->Rndm()) * 2;
  }

  for (auto &signal : *SiPM_Low)
  {
    signal.Channel = signal.Channel - 0.5 + gRandom->Rndm();
    signal.Time = signal.Time + (- 0.5 + gRandom->Rndm()) * 2;
  }
  /////////////////////////////////////

  
  int Rear_Label = (*Silicon)[0].Label;
  double Rear_Time = (*Silicon)[0].Time;
  double Rear_Channel = (*Silicon)[0].Channel;

  double mean_time_silicon = MeanAcceptance_Walk_Silicon[Rear_Label]->Eval(Rear_Channel); // silicon walk correction

  vector<Signal> SiPM_H[10];
  vector<Signal> SiPM_L[10];
  //////////// NEW ANALYSIS ///////////
  if (Rear_Channel > Rear_IAS[GetDetector(Rear_Label)].first && Rear_Channel < Rear_IAS[GetDetector(Rear_Label)].second)
  {
    Signal nearest = Signal();
    nearest.Time = 1000000;
    for (int i = 0; i < SiPM_High->GetSize(); i++)
    {
      double diff_time = Rear_Time - (*SiPM_High)[i].Time;
      double mean_time_SiPM = MeanAcceptance_Walk_SiPM[(*SiPM_High)[i].Label]->Eval((*SiPM_High)[i].Channel);
      Signal signal = (*SiPM_High)[i];
      signal.Time = diff_time + mean_time_silicon + mean_time_SiPM;
      SiPM_H[GetDetectorChannel((*SiPM_High)[i].Label)].push_back(signal);

      // if ((*SiPM_High)[0].Time > (*SiPM_High)[i].Time)
      // {
      //   cout << "ERROR" << endl;
      // }

      // saving the time difference between Rear and SiPMs nearest to 0
      if (abs(nearest.Time) > abs(signal.Time))
      {
        nearest = signal;
      }

      // time between Two SiPMs next to each other
      if (i > 0)
      {
        if ((*SiPM_High)[i].Label < (*SiPM_High)[i-1].Label)
        {
          H_2SiPM_Time_Next[(*SiPM_High)[i].Label][(*SiPM_High)[i-1].Label]->Fill((*SiPM_High)[i].Time - (*SiPM_High)[i-1].Time);
        }
        else
        {
          H_2SiPM_Time_Next[(*SiPM_High)[i-1].Label][(*SiPM_High)[i].Label]->Fill((*SiPM_High)[i].Time - (*SiPM_High)[i-1].Time);
        }
      }
    }

    for (int i = 0; i < SiPM_Low->GetSize(); i++)
    {
      double diff_time = Rear_Time - (*SiPM_Low)[i].Time;
      double mean_time_SiPM = MeanAcceptance_Walk_SiPM[(*SiPM_Low)[i].Label]->Eval((*SiPM_Low)[i].Channel);
      Signal signal = (*SiPM_Low)[i];
      signal.Time = diff_time + mean_time_silicon + mean_time_SiPM;
      SiPM_L[GetDetectorChannel((*SiPM_Low)[i].Label)].push_back(signal);
    }

    // Plotting only nearsest for high for each rear
    if (nearest.Time < 1000000)
    {
      H_RearSiPM_Time_Nearest[Rear_Label]->Fill(nearest.Time);
    }

    if (nearest.Time != 1000000)
    {
      // Plotting Time between nearsest SiPMs from the rear and other SIPM
      for (int SiPMi = 0; SiPMi < 10; SiPMi++)
      {
        for (int i = 0; i < SiPM_H[SiPMi].size(); i++)
        {
          if (nearest != SiPM_H[SiPMi][i])
          {
            if (nearest.Label < SiPM_H[SiPMi][i].Label)
            {
              H_RearSiPM_Time_Nearest_Nearest[nearest.Label][SiPM_H[SiPMi][i].Label]->Fill(nearest.Time - SiPM_H[SiPMi][i].Time);
              H_RearSiPM_Channel_Nearset_Nearest[nearest.Label][SiPM_H[SiPMi][i].Label]->Fill(nearest.Channel, SiPM_H[SiPMi][i].Channel);
            }
            else
            {
              H_RearSiPM_Time_Nearest_Nearest[SiPM_H[SiPMi][i].Label][nearest.Label]->Fill(nearest.Time - SiPM_H[SiPMi][i].Time);
              H_RearSiPM_Channel_Nearset_Nearest[SiPM_H[SiPMi][i].Label][nearest.Label]->Fill(SiPM_H[SiPMi][i].Channel, nearest.Channel);
            }
          }
        }
      }
    
      // grouping beta event within the faster group in subgroup
      vector<vector<Signal>> SiPM_HGroups;
      for (int i = 0; i < SiPM_High->GetSize(); i++)
      {
        double diff_time = Rear_Time - (*SiPM_High)[i].Time;
        double mean_time_SiPM = MeanAcceptance_Walk_SiPM[(*SiPM_High)[i].Label]->Eval((*SiPM_High)[i].Channel);
        (*SiPM_High)[i].Time = diff_time + mean_time_silicon + mean_time_SiPM;

        if (i == 0)
        {
          SiPM_HGroups.push_back({(*SiPM_High)[i]});
        }
        else
        {
          if ((abs(SiPM_HGroups[SiPM_HGroups.size()-1][(SiPM_HGroups[SiPM_HGroups.size()-1]).size()-1].Time - (*SiPM_High)[i].Time)) <= 20)
          {
            // same group
            SiPM_HGroups[SiPM_HGroups.size()-1].push_back((*SiPM_High)[i]);
          }
          else
          { 
            // new group
            SiPM_HGroups.push_back({(*SiPM_High)[i]});
          }
        }
      }

      vector<double> SiPM_HGroupsMeanTime;
      double lowest_mean_time = 1000000;
      double nearest_group_time = 1000000;
      int nearest_index = -1; 
      int lowest_index = -1;
      // loop on subgroup
      for (int i = 0; i < SiPM_HGroups.size(); i++)
      {
        double mean_time = 0;
        //loop on event of the subgroup
        for (int j = 0; j < SiPM_HGroups[i].size(); j++)
        {
          // mean of subgroup
          mean_time += SiPM_HGroups[i][j].Time;


          // nearest event whithin subgroups
          if (abs(SiPM_HGroups[i][j].Time) < abs(nearest_group_time))
          {
            nearest_group_time = SiPM_HGroups[i][j].Time;
            nearest_index = i;
          }
        }
        mean_time = mean_time / SiPM_HGroups[i].size();
        SiPM_HGroupsMeanTime.push_back(mean_time);

        if (abs(mean_time) < abs(lowest_mean_time))
        {
          lowest_mean_time = mean_time;
          lowest_index = i;
        }
      }

      // Plotting SiPM channel between SiPM in the nearest subgroup
      if (SiPM_HGroups[nearest_index].size() > 0)
      {
        H_RearSiPM_MeanTime_Nearest_Group->Fill(nearest_group_time);
        for (int i = 0; i < SiPM_HGroups[nearest_index].size(); i++)
        {
          for (int j = i; j < SiPM_HGroups[nearest_index].size(); j++)
          {
            if (SiPM_HGroups[nearest_index][i].Label < SiPM_HGroups[nearest_index][j].Label)
            {
              H_RearSiPM_Channel_Nearest_Group[SiPM_HGroups[nearest_index][i].Label][SiPM_HGroups[nearest_index][j].Label]->Fill(SiPM_HGroups[nearest_index][i].Channel, SiPM_HGroups[nearest_index][j].Channel);
            }
            else
            {
              H_RearSiPM_Channel_Nearest_Group[SiPM_HGroups[nearest_index][j].Label][SiPM_HGroups[nearest_index][i].Label]->Fill(SiPM_HGroups[nearest_index][j].Channel, SiPM_HGroups[nearest_index][i].Channel);
            }
          }
        }
      }

      // plotting mean subgroup time vs fartest SiPM in group time for each multiplicity

      // find furthest event in the group from the subgroup mean time
      int furthest_index = -1;
      double furthest_time = -1;
      for (int j = 0; j < SiPM_HGroups[nearest_index].size(); j++)
      {
        double mean = SiPM_HGroupsMeanTime[nearest_index];
        double time_SiPM = SiPM_HGroups[nearest_index][j].Time;

        if ( abs(mean - time_SiPM) > abs(furthest_time))
        {
          furthest_time = mean - time_SiPM;
          furthest_index = j;
        }
      }
      H_2SiPM_Time_MeanFurthest[SiPM_HGroups[nearest_index].size()]->Fill(furthest_time);

      
      H_RearSiPM_Time_New_Nearest[SiPM_HGroups[nearest_index].size()]->Fill(nearest_group_time);
      H_RearSiPM_Time_New_Mean[SiPM_HGroups[nearest_index].size()]->Fill(SiPM_HGroupsMeanTime[nearest_index]);
    }
    /////

    // plotting time difference between SiPMs
    for (int SiPMi = 0; SiPMi < 10; SiPMi++)
    {
      for (int SiPMj = SiPMi; SiPMj < 10; SiPMj++)
      {
        if (SiPM_H[SiPMi].size() == 1 && SiPM_H[SiPMj].size() == 1)
          H_2SiPM_Time_One[SiPM_H[SiPMi][0].Label][SiPM_H[SiPMj][0].Label]->Fill(SiPM_H[SiPMi][0].Time - SiPM_H[SiPMj][0].Time);
        for (int i = 0; i < SiPM_H[SiPMi].size(); i++)
        {
          for (int j = 0; j < SiPM_H[SiPMj].size(); j++)
          {
            H_2SiPM_Time_Multiple[SiPM_H[SiPMi][i].Label][SiPM_H[SiPMj][j].Label]->Fill(SiPM_H[SiPMi][i].Time - SiPM_H[SiPMj][j].Time);
          }
        }
      }
    }

    /// Plotting Time between Low and High for each SiPM
    for (int SiPMi = 1; SiPMi < 10; SiPMi++)
    {
      if (SiPM_H[SiPMi].size() == 1 && SiPM_L[SiPMi].size() == 1)
      {
        H_SiPMLowHigh_Time_One[SiPM_H[SiPMi][0].Label]->Fill(SiPM_H[SiPMi][0].Time - SiPM_L[SiPMi][0].Time);
        H_SiPMLowHigh_TimeChannel_One[SiPM_H[SiPMi][0].Label]->Fill(SiPM_H[SiPMi][0].Time - SiPM_L[SiPMi][0].Time, SiPM_H[SiPMi][0].Channel);
      }
      for (int i = 0; i < SiPM_H[SiPMi].size(); i++)
      {
        for (int j = 0; j < SiPM_L[SiPMi].size(); j++)
        {
          H_SiPMLowHigh_Time_Multiple[SiPM_H[SiPMi][i].Label]->Fill(SiPM_H[SiPMi][i].Time - SiPM_L[SiPMi][j].Time);
        }
      }
    }

  }

    

  // // Plotting Mean time Difference between Rear and SiPMs for each multiplciity
  // for (int SiPMi = 0; SiPMi < 10; SiPMi++)
  // {
  //   if (SiPM_H[SiPMi].size() == 1)
  //     H_RearSiPM_Time_One[SiPMi]->Fill(Rear_Time - SiPM_H[SiPMi][0].Time);
  //   for (int i = 0; i < SiPM_H[SiPMi].size(); i++)
  //   {
  //     H_RearSiPM_Time_Multiple[SiPMi]->Fill(Rear_Time - SiPM_H[SiPMi][i].Time);
  //   }
  // }


  /////////////////////////////////////

  int Multiplicity_High = 0;
  int Multiplicity_SiPMHigh[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int Multiplicity_Low = 0;
  int Multiplicity_SiPMLow[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  for (int i = 0; i < SiPM_High->GetSize(); i++)
  {
    double diff_time = Rear_Time - (*SiPM_High)[i].Time;
    double mean_time_SiPM = MeanAcceptance_Walk_SiPM[(*SiPM_High)[i].Label]->Eval((*SiPM_High)[i].Channel);
    pair<double, double> spread_time_SiPM = SpreadAcceptance_Walk_SiPM[GetDetectorChannel((*SiPM_High)[i].Label)];

    if (diff_time + mean_time_silicon + mean_time_SiPM > spread_time_SiPM.first && diff_time + mean_time_silicon + mean_time_SiPM < spread_time_SiPM.second)
    {
      H_RearSiPM_ChannelTime_Cleaned[Rear_Label]->Fill(diff_time, Rear_Channel);
      H_SiPM_ChannelTime_Cleaned[(*SiPM_High)[i].Label]->Fill(diff_time, (*SiPM_High)[i].Channel);
      Multiplicity_High += 1;
      Multiplicity_SiPMHigh[GetDetectorChannel((*SiPM_High)[i].Label)] += 1;

      CLEANED_Tree_SiPMHigh.push_back((*SiPM_High)[i]);
    }
  }

  for (int i = 0; i < SiPM_Low->GetSize(); i++)
  {
    double diff_time = Rear_Time - (*SiPM_Low)[i].Time;
    double mean_time_SiPM = MeanAcceptance_Walk_SiPM[(*SiPM_Low)[i].Label]->Eval((*SiPM_Low)[i].Channel);
    pair<double, double> spread_time_SiPM = SpreadAcceptance_Walk_SiPM[GetDetectorChannel((*SiPM_Low)[i].Label)];

    if (diff_time + mean_time_silicon + mean_time_SiPM > spread_time_SiPM.first && diff_time + mean_time_silicon + mean_time_SiPM < spread_time_SiPM.second)
    {
      H_RearSiPM_ChannelTime_Cleaned[Rear_Label]->Fill(diff_time, Rear_Channel);
      H_SiPM_ChannelTime_Cleaned[(*SiPM_Low)[i].Label]->Fill(diff_time, (*SiPM_Low)[i].Channel);
      Multiplicity_Low += 1;
      Multiplicity_SiPMLow[GetDetectorChannel((*SiPM_Low)[i].Label)] += 1;

      CLEANED_Tree_SiPMLow.push_back((*SiPM_Low)[i]);
    }
  }
  H_SiPMHigh_Mulitplicity_Cleaned->Fill(Multiplicity_High);
  H_SiPMLow_Mulitplicity_Cleaned->Fill(Multiplicity_Low);

  for (int i = 0; i < detectorNum; i++)
  {
    if (IsDetectorBetaHigh(i))
    {
      H_SiPM_Mulitplicity_Cleaned[i]->Fill(Multiplicity_SiPMHigh[GetDetectorChannel(i)]);
    }
    else if (IsDetectorBetaLow(i))
    {
      H_SiPM_Mulitplicity_Cleaned[i]->Fill(Multiplicity_SiPMLow[GetDetectorChannel(i)]);
    }
  }

  if ((*Silicon)[1].Pileup == 0 && (*Silicon)[0].Pileup == 0)
  {
    CLEANED_Tree->Fill();
    Tree_Channel_detector = (*Silicon)[1].Channel;
    CLEANED_Tree_detector[(*Silicon)[1].Label]->Fill();
  }
  CLEANED_Tree_Silicon.clear();
  CLEANED_Tree_SiPMHigh.clear();
  CLEANED_Tree_SiPMLow.clear();

  /////// IAS CHECK //////

  // counting mutiplicity for each time coinc
  int Multiplicity_High_IAS_in[99][99] = {0};
  int Multiplicity_High_IAS_out_l[99][99] = {0};
  int Multiplicity_High_IAS_out_r[99][99] = {0};
  int Multiplicity_Low_IAS_in[99][99] = {0};
  int Multiplicity_Low_IAS_out_l[99][99] = {0};
  int Multiplicity_Low_IAS_out_r[99][99] = {0};

  if (Rear_Channel > Rear_IAS[GetDetector(Rear_Label)].first && Rear_Channel < Rear_IAS[GetDetector(Rear_Label)].second)
  {
    for (int i = 0; i < SiPM_High->GetSize(); i++)
    {
      int SiPM_Label = (*SiPM_High)[i].Label;
      double diff_time = Rear_Time - (*SiPM_High)[i].Time;
      double mean_time_SiPM = MeanAcceptance_Walk_SiPM[SiPM_Label]->Eval((*SiPM_High)[i].Channel);

      for (int Time_Coinc_index_l = 0; Time_Coinc_index_l < CoincidenceTime.size(); Time_Coinc_index_l++)
      {
        for (int Time_Coinc_index_r = 0; Time_Coinc_index_r < CoincidenceTime.size(); Time_Coinc_index_r++)
        {
          if (diff_time + mean_time_silicon + mean_time_SiPM > -CoincidenceTime[Time_Coinc_index_l] && diff_time + mean_time_silicon + mean_time_SiPM < CoincidenceTime[Time_Coinc_index_r])
          {
            Multiplicity_High_IAS_in[Time_Coinc_index_l][Time_Coinc_index_r] += 1;
          }
          if (diff_time + mean_time_silicon + mean_time_SiPM < CoincidenceTime_Left_Window && diff_time + mean_time_silicon + mean_time_SiPM > CoincidenceTime_Left_Window-(CoincidenceTime[Time_Coinc_index_r]+CoincidenceTime[Time_Coinc_index_l]))
          {
            Multiplicity_High_IAS_out_l[Time_Coinc_index_l][Time_Coinc_index_r] += 1;
          }
          if (diff_time + mean_time_silicon + mean_time_SiPM > CoincidenceTime_Right_Window && diff_time + mean_time_silicon + mean_time_SiPM < CoincidenceTime_Right_Window+(CoincidenceTime[Time_Coinc_index_r]+CoincidenceTime[Time_Coinc_index_l]))
          {
            Multiplicity_High_IAS_out_r[Time_Coinc_index_l][Time_Coinc_index_r] += 1;
          }
        }
      }
    }
  }

  if (Rear_Channel > Rear_IAS[GetDetector(Rear_Label)].first && Rear_Channel < Rear_IAS[GetDetector(Rear_Label)].second)
  {
    for (int i = 0; i < SiPM_Low->GetSize(); i++)
    {
      int SiPM_Label = (*SiPM_Low)[i].Label;
      double diff_time = Rear_Time - (*SiPM_Low)[i].Time;
      double mean_time_SiPM = MeanAcceptance_Walk_SiPM[SiPM_Label]->Eval((*SiPM_Low)[i].Channel);
      for (int Time_Coinc_index_l = 0; Time_Coinc_index_l < CoincidenceTime.size(); Time_Coinc_index_l++)
      {
        for (int Time_Coinc_index_r = 0; Time_Coinc_index_r < CoincidenceTime.size(); Time_Coinc_index_r++)
        {
          if (diff_time + mean_time_silicon + mean_time_SiPM > -CoincidenceTime[Time_Coinc_index_l] && diff_time + mean_time_silicon + mean_time_SiPM < CoincidenceTime[Time_Coinc_index_r])
          {
            Multiplicity_Low_IAS_in[Time_Coinc_index_l][Time_Coinc_index_r] += 1;
          }
          if (diff_time + mean_time_silicon + mean_time_SiPM < CoincidenceTime_Left_Window && diff_time + mean_time_silicon + mean_time_SiPM > CoincidenceTime_Left_Window-(CoincidenceTime[Time_Coinc_index_r]+CoincidenceTime[Time_Coinc_index_l]))
          {
            Multiplicity_Low_IAS_out_l[Time_Coinc_index_l][Time_Coinc_index_r] += 1;
          }
          if (diff_time + mean_time_silicon + mean_time_SiPM > CoincidenceTime_Right_Window && diff_time + mean_time_silicon + mean_time_SiPM < CoincidenceTime_Right_Window+(CoincidenceTime[Time_Coinc_index_r]+CoincidenceTime[Time_Coinc_index_l]))
          {
            Multiplicity_Low_IAS_out_r[Time_Coinc_index_l][Time_Coinc_index_r] += 1;
          }
        }
      }
    }
  }

  //

  double ref = 0;
  if (Rear_Channel > Rear_IAS[GetDetector(Rear_Label)].first && Rear_Channel < Rear_IAS[GetDetector(Rear_Label)].second)
  {
    H_Rear_Channel_IAS_Cleaned[Rear_Label]->Fill(Rear_Channel);
    double energy_calib = 0;
    double channel = (*Silicon)[1].Channel/1000;
    if ((*Silicon)[1].Label > 50)
    {
        energy_calib = ManualCalibFitted[(*Silicon)[1].Label][0] + ManualCalibFitted[(*Silicon)[1].Label][1] * channel + ManualCalibFitted[(*Silicon)[1].Label][2] * pow(channel, 2);
        H_Energy_IAS->Fill(energy_calib);
    }

    for (int i = 0; i < SiPM_High->GetSize(); i++)
    {
      int SiPM_Label = (*SiPM_High)[i].Label;
      double diff_time = Rear_Time - (*SiPM_High)[i].Time;
      double mean_time_SiPM = MeanAcceptance_Walk_SiPM[SiPM_Label]->Eval((*SiPM_High)[i].Channel);

      for (int Time_Coinc_index_l = 0; Time_Coinc_index_l < CoincidenceTime.size(); Time_Coinc_index_l++)
      {
        for (int Time_Coinc_index_r = 0; Time_Coinc_index_r < CoincidenceTime.size(); Time_Coinc_index_r++)
        {
          if (diff_time + mean_time_silicon + mean_time_SiPM > -CoincidenceTime[Time_Coinc_index_l] && diff_time + mean_time_silicon + mean_time_SiPM < CoincidenceTime[Time_Coinc_index_r])
          {
            for (int multiplicity = 1; multiplicity <= Multiplicity_High_IAS_in[Time_Coinc_index_l][Time_Coinc_index_r]; multiplicity++)
            {
              H_HighSiPM_Time_IAS_Cleaned[Time_Coinc_index_l][Time_Coinc_index_r][multiplicity]->Fill(diff_time + mean_time_silicon + mean_time_SiPM);
              if ((*Silicon)[1].Label > 50)
              {
                H_Energy_IAS_Cleaned[Time_Coinc_index_l][Time_Coinc_index_r][multiplicity]->Fill(energy_calib, (double)1./Multiplicity_High_IAS_in[Time_Coinc_index_l][Time_Coinc_index_r]);
              }
            }           
            if (ref != 0 && GetDetectorChannel(SiPM_Label) != 7)
            {
              H_2SiPM_Time_IAS_Cleaned[SiPM_Label]->Fill(diff_time + mean_time_silicon + mean_time_SiPM - ref);
            }
            if (ref == 0 && GetDetectorChannel(SiPM_Label) == 7)
            {
              ref = diff_time + mean_time_silicon + mean_time_SiPM;
            }
          }
          if (diff_time + mean_time_silicon + mean_time_SiPM < CoincidenceTime_Left_Window && diff_time + mean_time_silicon + mean_time_SiPM > CoincidenceTime_Left_Window-(CoincidenceTime[Time_Coinc_index_r]+CoincidenceTime[Time_Coinc_index_l]))
          {
            for (int multiplicity = 1; multiplicity <= Multiplicity_High_IAS_out_l[Time_Coinc_index_l][Time_Coinc_index_r]; multiplicity++)
            {
              H_HighSiPM_Time_IAS_UnCleaned[Time_Coinc_index_l][Time_Coinc_index_r][multiplicity]->Fill(diff_time + mean_time_silicon + mean_time_SiPM);
              if ((*Silicon)[1].Label > 50)
              {
                H_Energy_IAS_UnCleaned[Time_Coinc_index_l][Time_Coinc_index_r][multiplicity]->Fill(energy_calib, (double)1./Multiplicity_High_IAS_out_l[Time_Coinc_index_l][Time_Coinc_index_r]);
              }
            }
          }
          if (diff_time + mean_time_silicon + mean_time_SiPM > CoincidenceTime_Right_Window && diff_time + mean_time_silicon + mean_time_SiPM < CoincidenceTime_Right_Window+(CoincidenceTime[Time_Coinc_index_r]+CoincidenceTime[Time_Coinc_index_l]))
          {
            for (int multiplicity = 1; multiplicity <= Multiplicity_High_IAS_out_r[Time_Coinc_index_l][Time_Coinc_index_r]; multiplicity++)
            {
              H_HighSiPM_Time_IAS_UnCleaned[Time_Coinc_index_l][Time_Coinc_index_r][multiplicity]->Fill(diff_time + mean_time_silicon + mean_time_SiPM);
              if ((*Silicon)[1].Label > 50)
              {
                H_Energy_IAS_UnCleaned[Time_Coinc_index_l][Time_Coinc_index_r][multiplicity]->Fill(energy_calib, (double)1./Multiplicity_High_IAS_out_r[Time_Coinc_index_l][Time_Coinc_index_r]);
              }
            }
          }
        }
      }
    }

    for (int i = 0; i < SiPM_Low->GetSize(); i++)
    {
      int SiPM_Label = (*SiPM_Low)[i].Label;
      double diff_time = Rear_Time - (*SiPM_Low)[i].Time;
      double mean_time_SiPM = MeanAcceptance_Walk_SiPM[SiPM_Label]->Eval((*SiPM_Low)[i].Channel);


           
      for (int Time_Coinc_index_l = 0; Time_Coinc_index_l < CoincidenceTime.size(); Time_Coinc_index_l++)
      {
        for (int Time_Coinc_index_r = 0; Time_Coinc_index_r < CoincidenceTime.size(); Time_Coinc_index_r++)
        {
          if (diff_time + mean_time_silicon + mean_time_SiPM > -CoincidenceTime[Time_Coinc_index_l] && diff_time + mean_time_silicon + mean_time_SiPM < CoincidenceTime[Time_Coinc_index_r])
          {
            for (int multiplicity = 1; multiplicity <= Multiplicity_Low_IAS_in[Time_Coinc_index_l][Time_Coinc_index_r]; multiplicity++)
            {
              H_LowSiPM_Time_IAS_Cleaned[Time_Coinc_index_l][Time_Coinc_index_r][multiplicity]->Fill(diff_time + mean_time_silicon + mean_time_SiPM);
            }
          }
          if (diff_time + mean_time_silicon + mean_time_SiPM < CoincidenceTime_Left_Window && diff_time + mean_time_silicon + mean_time_SiPM > CoincidenceTime_Left_Window-(CoincidenceTime[Time_Coinc_index_r]+CoincidenceTime[Time_Coinc_index_l]))
          {
            for (int multiplicity = 1; multiplicity <= Multiplicity_Low_IAS_out_l[Time_Coinc_index_l][Time_Coinc_index_r]; multiplicity++)
            {
              H_LowSiPM_Time_IAS_UnCleaned[Time_Coinc_index_l][Time_Coinc_index_r][multiplicity]->Fill(diff_time + mean_time_silicon + mean_time_SiPM);
            }
          }
          if (diff_time + mean_time_silicon + mean_time_SiPM > CoincidenceTime_Right_Window && diff_time + mean_time_silicon + mean_time_SiPM < CoincidenceTime_Right_Window+(CoincidenceTime[Time_Coinc_index_r]+CoincidenceTime[Time_Coinc_index_l]))
          {
            for (int multiplicity = 1; multiplicity <= Multiplicity_Low_IAS_out_r[Time_Coinc_index_l][Time_Coinc_index_r]; multiplicity++)
            {
              H_LowSiPM_Time_IAS_UnCleaned[Time_Coinc_index_l][Time_Coinc_index_r][multiplicity]->Fill(diff_time + mean_time_silicon + mean_time_SiPM);
            }
          }
        }
      }
    }

    if (Multiplicity_High >= 9)
    {
      H_Rear_Channel_IAScoinc_Cleaned[Rear_Label]->Fill(Rear_Channel);
    }
    ///////////////////////
  }
}

void SaveFitParameters()
{
  TFile *file = new TFile((DIR_ROOT_DATA_GROUPED + "Grouping_FitParameters.root").c_str(), "RECREATE");

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
  }

  file->Close();
}

void LoadFitParameters()
{
  TFile *file = new TFile((DIR_ROOT_DATA_GROUPED + "Grouping_FitParameters.root").c_str(), "READ");
  for (int i = 0; i < detectorNum; i++)
  {
    if (IsDetectorBeta(i))
    {
      MeanAcceptance_Walk_SiPM[i] = (TF1 *)file->Get(("MeanAcceptance_Walk_SiPM_" + detectorName[i]).c_str());
    }
    else if (IsDetectorSiliBack(i))
    {
      MeanAcceptance_Walk_Silicon[i] = (TF1 *)file->Get(("MeanAcceptance_Walk_Silicon_" + detectorName[i]).c_str());
    }
  }
  file->Close();
  GROUPED_File->cd();
}

#endif