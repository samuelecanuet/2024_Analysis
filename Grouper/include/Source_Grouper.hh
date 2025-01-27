// #ifndef GROUPER_HH
// #define GROUPER_HH

#include "../../../lib/SignalDict/Signal.h"
#include "Detectors.hh"


using namespace std;

int VERBOSE = 0;
int Run;
bool FULL = false;

/// FILE ///
string ROOT_filename;
string ROOT_basefilename;
TFile *GROUPED_File;

/// TREE ///
TTree *CLEANED_Tree;
TTree *CLEANED_Tree_detector[SIGNAL_MAX];
vector<Signal> CLEANED_Tree_Silicon;
vector<vector<pair<Signal, Signal>>> CLEANED_Tree_SiPMGroup;
double Tree_Channel_detector;

////////////// HISTOGRAMS ////////////////

// RAW
TH1D *H_Channel_RAW[SIGNAL_MAX];
TH1D *H_SiPMHigh_Mulitplicity_RAW;
TH1D *H_SiPMLow_Mulitplicity_RAW;

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
TCanvas *C_IAS_Channel_Cleaned;
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

// FITS
TCanvas *C_Strip_Cutting_Fits;
TCanvas *C_Rear_Walk_Silicon_Fits;
TCanvas *C_SiPM_Walk_SiPM_Fits;


///// DIRECTORY //////////////////
// Group
TDirectory *dir_Group;
TDirectory *dir_Group_Strip;
TDirectory *dir_Group_Rear;
// Detector
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

TDirectory *dir_Cleaned;
TDirectory *dir_Cleaned_Multiplicity;
TDirectory *dir_Cleaned_HighSiPM_TimeCoincidence_Mulitplicity[SIGNAL_MAX];
TDirectory *dir_Cleaned_LowSiPM_TimeCoincidence_Mulitplicity[SIGNAL_MAX];
TDirectory *dir_Matching_Cleaned;

///////// FUCNTION ARRAY /////////////
TF1 *SpreadAcceptance[SIGNAL_MAX];
double MeanAcceptance[SIGNAL_MAX];
TF1 *MeanAcceptance_Walk_Silicon[SIGNAL_MAX];
TF1 *MeanAcceptance_Walk_SiPM[SIGNAL_MAX];

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

void SaveFitParameters()
{
  TFile *file = new TFile((DIR_ROOT_DATA_GROUPED + "Grouping_FitParameters.root").c_str(), "RECREATE");

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
  TFile *file = new TFile((DIR_ROOT_DATA_GROUPED + "Grouping_FitParameters.root").c_str(), "READ");
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
        Error(("MeanAcceptance_" + detectorName[i]).c_str());
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

  for (size_t i = 0; i < detectorNum; ++i)
  {
    if (IsDetectorBetaHigh(i))
    {
      dir_SiPMHigh_Channel_Detector[i] = dir_SiPMHigh_Channel->mkdir(detectorName[i].c_str());
      dir_SiPMHigh_Time_Detector[i] = dir_SiPMHigh_Time->mkdir(detectorName[i].c_str());

      H_Channel_RAW[i] = new TH1D(("Channel_RAW_" + detectorName[i]).c_str(), ("Channel_RAW_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
      H_Channel_RAW[i]->GetXaxis()->SetTitle("Channel");
      H_Channel_RAW[i]->GetYaxis()->SetTitle("Counts");
      H_Channel_RAW[i]->GetXaxis()->CenterTitle();
      H_Channel_RAW[i]->GetYaxis()->CenterTitle();

      H_SiPM_ChannelTime_Cleaned[i] = new TH2D(("SiPM_ChannelTime_Cleaned_" + detectorName[i]).c_str(), ("SiPM_ChannelTime_Cleaned_" + detectorName[i]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta, eHighN / 10, eHighMin, eHighMax);
      H_SiPM_ChannelTime_Cleaned[i]->GetXaxis()->SetTitle("Rear-SiPM Time (ns)");
      H_SiPM_ChannelTime_Cleaned[i]->GetYaxis()->SetTitle("SiPM Channel");
      H_SiPM_ChannelTime_Cleaned[i]->GetXaxis()->CenterTitle();
      H_SiPM_ChannelTime_Cleaned[i]->GetYaxis()->CenterTitle();

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

      H_2SiPM_Channel_Couple[i] = new TH2D(("2SiPM_Channel_Couple_" + detectorName[i]).c_str(), ("2SiPM_Channel_Couple_" + detectorName[i]).c_str(), eLowN, eLowMin, eLowMax, eHighN, eHighMin, eHighMax);
      H_2SiPM_Channel_Couple[i]->GetXaxis()->SetTitle("SiPM Low Channel");
      H_2SiPM_Channel_Couple[i]->GetYaxis()->SetTitle("SiPM High Channel");
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

  return 0;
}

inline void WriteHistograms_Cleaned()
{
  GROUPED_File->cd();
  // H_RearSiPM_MeanTime_Nearest_Group->Write();
  for (int i = 0; i < detectorNum; i++)
  {
    ProgressCounter(i, detectorNum, " <INFO>  Writing Histograms");

    if (IsDetectorBetaHigh(i))
    {
      dir_SiPMHigh_Channel_Detector[i]->cd();
      H_SiPM_Channel_Alone[i]->Write();
      H_SiPM_Channel_Couple[i]->Write();
      H_2SiPM_Channel_Couple[i]->Write();
      H_Channel_Cleaned[i]->Write();

      dir_SiPMHigh_Time_Detector[i]->cd();
      H_SiPM_ChannelTime_Cleaned[i]->Write();
      H_SiPM_Pair_Time[i]->Write();

      for (int j = i + 1; j < SIGNAL_MAX; j++)
      {
        if (IsDetectorBetaHigh(j))
        {
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

  /// GENERAL FITS CANVAS ///
  C_Strip_Cutting_Fits->Write();
  C_Rear_Walk_Silicon_Fits->Write();
  C_SiPM_Walk_SiPM_Fits->Write();
}

inline int WriteTree_Grouped()
{
  GROUPED_File->cd();
  CLEANED_Tree->Write();

  delete CLEANED_Tree;
  return 0;
}

inline void CleaningGroups(TTreeReaderArray<Signal> &signals)
{

  int TIME_BETWEEN_SAME_GAIN = 20;
  int TIME_BETWEEN_DIFFERENT_GAIN = 30;

  ////////////////////####################////////////////////

  vector<Signal> Silicon;
  vector<Signal> SiPM_High;
  vector<Signal> SiPM_Low;

  double Rear_Time = signals[0].Time;

  for (int index = 0; index < signals.GetSize(); index++)
  {
    int current_label = signals[index].Label;
    // cout << current_label << endl;
    if (IsDetectorBetaHigh(current_label))
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
  }

  ////////////////////####################////////////////////

  /// SiPM
  double SiPM_Time = 0;
  int SiPM_High_Size = SiPM_High.size();
  int SiPM_Low_Size = SiPM_Low.size();

  ////////////////////####################////////////////////

  ////////////////////####################////////////////////
  /// CONTINUOUS VALUE OF CHANNELS ///
  for (auto &signal : SiPM_High)
  {
    // signal.Time = CLEANED_Tree_Silicon[0].Time - signal.Time + MeanAcceptance_Walk_SiPM[signal.Label]->Eval(signal.Channel) + MeanAcceptance_Walk_Silicon[CLEANED_Tree_Silicon[0].Label]->Eval(CLEANED_Tree_Silicon[0].Channel);
  }

  for (auto &signal : SiPM_Low)
  {
    // signal.Time = CLEANED_Tree_Silicon[0].Time - signal.Time + MeanAcceptance_Walk_SiPM[signal.Label]->Eval(signal.Channel) + MeanAcceptance_Walk_Silicon[CLEANED_Tree_Silicon[0].Label]->Eval(CLEANED_Tree_Silicon[0].Channel);
  }

  /////////////////////////////////////
  vector<Signal> SiPM_H[10]; // vector of High per Label
  vector<Signal> SiPM_L[10]; // vector of Low per Label
  //////////// NEW ANALYSIS ///////////

  // #######################
  // ######### (A) #########
  // cout << "A" << endl;
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
          H_2SiPM_Time_Next[(SiPM_High)[i].Label][(SiPM_High)[i - 1].Label]->Fill((SiPM_High)[i].Time - (SiPM_High)[i - 1].Time);
      }
      else
      {
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
          H_2SiPM_Time_Next[(SiPM_Low)[i].Label][(SiPM_Low)[i - 1].Label]->Fill((SiPM_Low)[i].Time - (SiPM_Low)[i - 1].Time);
      }
      else
      {
          H_2SiPM_Time_Next[(SiPM_Low)[i - 1].Label][(SiPM_Low)[i].Label]->Fill((SiPM_Low)[i].Time - (SiPM_Low)[i - 1].Time);
      }
    }
  }
  // #######################
  // #######################

  // #######################
  // ######### (B) #########
  // cout << "B" << endl;
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
  // cout << "C" << endl;
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
      double time = -(diff_time + mean_time_SiPM);

      if (first > time || first == 10e6)
      {
        first = time;
      }
      if (last < time || last == -10e6)
      {
        last = time;
      }
    }

    // Loop on All low
    for (int i_low = 0; i_low < SiPM_Low.size(); i_low++)
    {
      // Only if Alone
      if (!Associated_Low[SiPM_Low[i_low].Label])
      {
        double diff_time = Rear_Time - SiPM_Low[i_low].Time;
        double mean_time_SiPM = MeanAcceptance_Walk_SiPM[SiPM_Low[i_low].Label]->Eval(SiPM_Low[i_low].Channel);
        double low_alone = -(diff_time + mean_time_SiPM);
        if (low_alone > first - TIME_BETWEEN_DIFFERENT_GAIN && low_alone < last + TIME_BETWEEN_DIFFERENT_GAIN)
        {
          SiPM_Groups[i_group].push_back(make_pair(Signal(), SiPM_Low[i_low]));
          Associated_Low[SiPM_Low[i_low].Label] = true;
        }
      }
    }
  }


  // # Plotting High Low in groups #
  // - time difference
  // - channel (1D and 2D) both and alone
  if (SiPM_Groups.size() == 1)
    {
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
    }

  // #######################
  // ######### (C) #########
  // cout << "D" << endl;
  // # Saving in Tree #
  CLEANED_Tree_SiPMGroup = SiPM_Groups;

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

  // # TIME FASTER SORTED #
  if (SiPM_Groups.size() > 0 && Verbose && SiPM_High.size() > 9)
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
                 << setw(15) << signals[i_signal].Time - Rear_Time
                 << setw(6) << "Corrected Time : "
                 << setw(15)
                 << -(diff_time + mean_time_SiPM)
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
             << setw(15) << signals[i_signal].Time - Rear_Time
             << setw(6) << "Corrected Time : "
             << setw(15)
             << -(diff_time + mean_time_SiPM)
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
        double t = -(diff_time + mean_time_SiPM);
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
  

  /////////////////////////////////////

  CLEANED_Tree->Fill();
  CLEANED_Tree_SiPMGroup.clear();

  //////////////////////////////////////
}

// #endif