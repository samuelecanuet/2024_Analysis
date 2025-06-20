#ifndef DETECTORS_HH
#define DETECTORS_HH

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

//  ROOT includes
//// files ////
#include "TSystem.h"
#include "TFile.h"
//// trees ////
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TNtuple.h"
//// data containers ////
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"	
#include "TGraph.h"	
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
//// plotting ////
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"
#include "TStyle.h"
#include "TArrow.h" 
#include "TMarker3DBox.h"
#include "TPolyLine3D.h"
#include "TView.h"
//// functions ////
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TFitResult.h"
//// others ////
#include "TRandom.h"
#include "TMinuit.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include <TVectorD.h>
#include <TMatrixDSym.h>

#include "Messenger.hh"
#include "Utilities.hh"


#define FDATA_MAX 80  ///< Maximum coder label
#define SIGNAL_MAX 190 ///< Maximum number of signals
#define BETA_NUM 2    ///< Number of beta (SiPM) detector groups
#define SILI_NUM 8    ///< Number of silicon detector groups
#define BETA_SIZE 9   ///< Number of channels in beta (SiPM) detector groups
#define SILI_SIZE 6   ///< Number of channels in silicon detector groups
#define MAX_MULTIPLICTY 9

#define BETA_HI 1
#define BETA_LO 2

using namespace std;
using namespace ROOT::Math;

/// FOLDERS ///
string DIR_FAST_DATA;
string DIR_ROOT_DATA;
string DIR_ROOT_DATA_GROUPED;
string DIR_ROOT_DATA_CLEANED;
string DIR_ROOT_DATA_MATCHED;
string DIR_ROOT_DATA_MERGED;
string DIR_ROOT_DATA_ANALYSED;
string DIR_ROOT_DATA_CALIBRATED;
string DIR_ROOT_DATA_MCP;
string DIR_ROOT_DATA_RATE;
string DIR_ROOT_DATA_SIMULATED;
string DIR_DATA_ISOLDE;
string DIR_DATA_HDD;
///////////////


// YEARS
bool FLAG2021 = false;
bool FLAG2024 = false;
bool FLAG2025 = false;
int YEAR = -1;
int REFERENCE_RUN = -1;
int MATCHING_RUN = -1;
map<string, vector<string>> Map_RunFiles;
// 

string detectorFileName; ///< Detectors definition file name
size_t detectorNum = SIGNAL_MAX;      ///< Number of defined detectors channels
string detectorName[SIGNAL_MAX];

map<string, int> NucleusColor;

// DATA
map<string, double> IAS;
map<string, map<double, vector<double>>> PeakData; ///< Peak data for each nucleus and peak type
map<string, map<double, double>> WindowsBetaMap;
map<string, double> Qbeta; ///< Beta energy for each nucleus
map<string, vector<double>> PeakList; ///< List of peaks for each nucleus
map<string, map<double, pair<double, double>[SIGNAL_MAX]>> WindowsMap; ///< Energy windows for each nucleus and peak type
map<string, pair<int, int>> CanvasMap; ///< Canvas for each nucleus

/// Group
int mini = -500;
int maxi = 500;

int n_Sili = 8;
int n_Beta = 2;

double winGroupMin_Sili = mini - n_Sili + mini % n_Sili;
double winGroupMax_Sili = maxi + n_Sili - maxi % n_Sili;
int winGroupN_Sili = (abs(winGroupMin_Sili) + winGroupMax_Sili) / n_Sili;

double winGroupMin_Beta = mini - n_Beta + mini % n_Beta;
double winGroupMax_Beta = maxi + n_Beta - maxi % n_Beta;
int winGroupN_Beta = (abs(winGroupMin_Beta) + winGroupMax_Beta) / n_Beta;

/// Silicon ///
double winSiliMin = -50;
double winSiliMax = 50;
int winSiliN = (abs(winSiliMin) + winSiliMax) / 8;
double eSiliMin = 0;
double eSiliMax = 100000;
int eSiliN = 10000;

double eSiliMin_cal = 0;
double eSiliMax_cal = 10000;
int eSiliN_cal = 100000;

/// SiPM ///
/// High
double winHighMin = 0;
double winHighMax = 250;
int winHighN = (abs(winHighMin) + winHighMax) / 2;
double eHighMin = 0;
double eHighMax = 6000000;
int eHighN = 3000;

/// Low
double winLowMin = 0;
double winLowMax = 250;
int winLowN = (abs(winLowMin) + winLowMax) / 2;
double eLowMin = 0;
double eLowMax = 1200000;
int eLowN = 3000;

// limit
double eHighLowLimit = 2.e6;

/// TOTAL WINDOW ///
double tabMIN[3] = {winSiliMin, winHighMin, winLowMin};
double winTotalMin = gsl_stats_min(tabMIN, 1, 3);
double tabMAX[3] = {winSiliMax, winHighMax, winLowMax};
double winTotalMax = gsl_stats_max(tabMAX, 1, 3);

int YearConverter(int det)
{

  // SiPM LOW
  if (det == 1) return 111;
  if (det == 2) return 112;
  if (det == 3) return 113;
  if (det == 4) return 114;
  if (det == 5) return 115;
  if (det == 6) return 116;
  if (det == 7) return 117;
  if (det == 8) return 118;
  if (det == 10) return 119;

  // SiPM HIGH
  if (det == 11) return 101;
  if (det == 12) return 102;
  if (det == 13) return 103;
  if (det == 14) return 104;
  if (det == 15) return 105;
  if (det == 16) return 106;
  if (det == 17) return 107;
  if (det == 18) return 108;
  if (det == 9) return 109;

  // Silicon
  if (det == 19) return 15;
  if (det == 20) return 14;
  if (det == 21) return 13;
  if (det == 22) return 11;
  if (det == 23) return 12;
  if (det == 24) return 16;

  if (det == 25) return 21;
  if (det == 26) return 22;
  if (det == 27) return 23;
  if (det == 28) return 24;
  if (det == 29) return 25;
  if (det == 30) return 26;

  if (det == 31) return 35;
  if (det == 32) return 34;
  if (det == 33) return 33;
  if (det == 34) return 31;
  if (det == 35) return 32;
  if (det == 36) return 36;

  if (det == 37) return 45;
  if (det == 38) return 44;
  if (det == 39) return 43;
  if (det == 40) return 41;
  if (det == 41) return 42;
  if (det == 42) return 46;

  if (det == 43) return 55;
  if (det == 44) return 51;
  if (det == 45) return 53;
  if (det == 46) return 54;
  if (det == 47) return 52;
  if (det == 48) return 56;

  if (det == 49) return 65;
  if (det == 50) return 64;
  if (det == 51) return 63;
  if (det == 52) return 61;
  if (det == 53) return 62;
  if (det == 54) return 66;

  if (det == 55) return 75;
  if (det == 56) return 74;
  if (det == 57) return 73;
  if (det == 58) return 71;
  if (det == 59) return 72;
  if (det == 60) return 76;

  if (det == 61) return 85;
  if (det == 62) return 84;
  if (det == 63) return 83;
  if (det == 64) return 81;
  if (det == 65) return 82;
  if (det == 66) return 86;

  return det;

}

int YearConverterInverter(int det)
{
    if (YEAR != 2021) return det;

    // SiPM LOW
    if (det == 111) return 1;
    if (det == 112) return 2;
    if (det == 113) return 3;
    if (det == 114) return 4;
    if (det == 115) return 5;
    if (det == 116) return 6;
    if (det == 117) return 7;
    if (det == 118) return 8;
    if (det == 119) return 10;
  
    // SiPM HIGH
    if (det == 101) return 11;
    if (det == 102) return 12;
    if (det == 103) return 13;
    if (det == 104) return 14;
    if (det == 105) return 15;
    if (det == 106) return 16;
    if (det == 107) return 17;
    if (det == 108) return 18;
    if (det == 109) return 9;
  
    // Silicon
    if (det == 15) return 19;
    if (det == 14) return 20;
    if (det == 13) return 21;
    if (det == 11) return 22;
    if (det == 12) return 23;
    if (det == 16) return 24;
  
    if (det == 21) return 25;
    if (det == 22) return 26;
    if (det == 23) return 27;
    if (det == 24) return 28;
    if (det == 25) return 29;
    if (det == 26) return 30;

    if (det == 35) return 31;
    if (det == 34) return 32;
    if (det == 33) return 33;
    if (det == 31) return 34;
    if (det == 32) return 35;
    if (det == 36) return 36;

    if (det == 45) return 37;
    if (det == 44) return 38;
    if (det == 43) return 39;
    if (det == 41) return 40;
    if (det == 42) return 41;
    if (det == 46) return 42;

    if (det == 55) return 43;
    if (det == 51) return 44;
    if (det == 53) return 45;
    if (det == 54) return 46;
    if (det == 52) return 47;
    if (det == 56) return 48;

    if (det == 65) return 49;
    if (det == 64) return 50;
    if (det == 63) return 51;
    if (det == 61) return 52;
    if (det == 62) return 53;
    if (det == 66) return 54;

    if (det == 75) return 55;
    if (det == 74) return 56;
    if (det == 73) return 57;
    if (det == 71) return 58;
    if (det == 72) return 59;
    if (det == 76) return 60;

    if (det == 85) return 61;
    if (det == 84) return 62;
    if (det == 83) return 63;
    if (det == 81) return 64;
    if (det == 82) return 65;
    if (det == 86) return 66;

    return 0;
}



inline bool IsDetectorBetaLow(int det)
{
  return (det >= 111 && det <= 119);
}

inline bool IsDetectorBetaHigh(int det)
{
  return (det >= 101 && det <= 109);
}

inline bool IsDetectorBeta(int det)
{
  return (IsDetectorBetaHigh(det) || IsDetectorBetaLow(det));
}

inline int GetDetectorChannel(int det)
{
  return (det % 10);
}

inline int GetDetector(int det)
{
  return (det / 10);
}

inline bool IsDetectorSili(int det)
{
  int num = GetDetector(det);
  int channel = GetDetectorChannel(det);
  return (num >= 1 && num <= 8 && channel >= 1 && channel <= 6);
}

//----------------------------------------------------------------------
inline bool IsDetectorSiliBack(int det)
{
  return (IsDetectorSili(det) && (GetDetectorChannel(det) == 6));
}

inline bool IsDetectorSiliStrip(int det)
{
  return (IsDetectorSili(det) && (GetDetectorChannel(det) != 6));
}

inline bool IsSameSiliDetector(int det1, int det2)
{
  return (GetDetector(det1) == GetDetector(det2));
}

inline int HighLowLabelConversion(int det)
{

  if (!IsDetectorBeta(det))
  {
    Error("Using HighLowLabelConversion for non-beta detector");
  }
  if (IsDetectorBetaLow(det))
  {
    return (det - 10);
  }
  else if (IsDetectorBetaHigh(det))
  {
    return (det + 10);
  }
}

inline bool IsDetectorSiliInterStrip(int det1, int det2)
{
  if (IsDetectorSiliStrip(det1) && IsDetectorSiliStrip(det2))
  {
    if (abs(GetDetectorChannel(det1) - GetDetectorChannel(det2)) == 1)
    {
      if (IsSameSiliDetector(det1, det2))
      {
        return true;
      }
    }
  }
  return false;
}

inline bool IsHRSProton(int label)
{
  if (label == 99)
  {
    return true;
  }
  return false;
}

vector<int> Dir2Det(string dir, int strip)
{
    vector<int> detectors;
    for (int i = 0; i < detectorNum; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            if (GetDetectorChannel(i) == strip && dir == "Up" && GetDetector(i) <= 4)
            {
                detectors.push_back(i);
            }
            else if (GetDetectorChannel(i) == strip && dir == "Down" && GetDetector(i) >= 5)
            {
                detectors.push_back(i);
            }
        }
    }
    return detectors;
}

void InitPeakData(string NUCLEUS, string Source = "ENSDF")
{
    Qbeta["32Ar"] = 11134-1200;
    Qbeta["33Ar"] = 11620-2200;
    Qbeta["32Cl"] = 12680-8800;

    string path = "/home/local1/Documents/2024_Analysis/Grouper/Config_Files/";
    std::ifstream file;
    if (NUCLEUS.find("32Ar") != string::npos)
    {
        file.open(path + "32Ar_protons_" + Source + ".txt");
        if (!file.is_open())
        {
            Warning("Impossible to open 32Ar_protons_" + Source + ".txt");
        }
    }
    else if (NUCLEUS.find("33Ar") != string::npos)
    {
        file.open(path + "33Ar_protons_" + Source + ".txt");
        if (!file.is_open())
        {
            Warning("Impossible to open 33Ar_protons_" + Source + ".txt");
        }
    }
    else if (NUCLEUS.find("18N") != string::npos)
    {
        file.open(path + "18N_protons_" + Source + ".txt");
        if (!file.is_open())
        {
            Warning("Impossible to open 18N_protons_" + Source + ".txt");
        }
    }
    else if (NUCLEUS.find("32Cl") != string::npos)
    {
        file.open(path + "32Cl_delayed_" + Source + ".txt");
        if (!file.is_open())
        {
            Warning("Impossible to open 32Cl_delayed_" + Source + ".txt");
        }
    }
    else
    {
        Warning("No energy error file for " + NUCLEUS);
        return;
    }

    string line;
    double error;
    double energy, BR, BR_error;
    double peak;
    int counter = 0;
    while (getline(file, line))
    {
        counter++;
        stringstream ss(line);
        ss >> peak >> energy >> error >> BR >> BR_error;
        PeakData[NUCLEUS][peak] = vector<double>({energy, error, BR, BR_error});
        WindowsBetaMap[NUCLEUS][peak] = Qbeta[NUCLEUS] - energy;    
    }

    file.close();
}

void InitRuns()
{
  string path = "/home/local1/Documents/2024_Analysis/Grouper/Config_Files/";
  ifstream file((path + to_string(YEAR) + "/Runs_" + to_string(YEAR) + ".txt").c_str());
  if (!file.is_open())
  {
    Error("Could not open the file Runs_" + to_string(YEAR) + ".txt");
  }

  string line;
  string nucleus;
  while (getline(file, line))
  {
    if (line.empty())
      continue;

    if (line[0] == '#')
      nucleus = line.substr(1);

    else
    {
      stringstream ss(line);
      string number;
      while (ss >> number)
      {
        Map_RunFiles[nucleus].push_back(number);
      }
    }
  }

  file.close();

  Info("Runs loaded");
  for (const auto &pair : Map_RunFiles)
  {
    string nucleus = pair.first;
    vector<string> runs = pair.second;
    Info("Nucleus : " + nucleus, 1);
    string runstring = "";
    for (const auto &run : runs)
    {
      runstring += run + " ";
    }
    Info(runstring, 2);
  }
}

void InitWindows()
{

    CanvasMap["32Ar"] = make_pair(8, 5);
    CanvasMap["32Ar_thick"] = make_pair(8, 5);
    CanvasMap["33Ar"] = make_pair(10, 5);
    CanvasMap["33Ar_thick"] = make_pair(10, 5);
    CanvasMap["18N"] = make_pair(3, 3);
    
    string direction[2] = {"Up", "Down"};
    string path = "/home/local1/Documents/2024_Analysis/Grouper/Config_Files/";
    for (auto dir : direction)
    {
        for (int strip = 1; strip <= 5; strip++)
        {
            ifstream file(path + "Detector_Window/" + dir + "_" + to_string(strip) + ".txt");
            if (!file.is_open())
            {
                Error("Impossible to open " + dir + "_" + to_string(strip) + ".txt");
            }

            string line;
            double energy_low;
            double energy_high;
            double number;
            string nuclei;
            while (getline(file, line))
            {
                energy_high = -1;
                energy_low = -1;

                if (line.empty())
                {
                    continue;
                }

                if (line.find("#") != string::npos)
                {
                    nuclei = line.substr(1);
                    continue;
                }
                stringstream ss(line);
                ss >> number >> energy_low >> energy_high;

                for (int i : Dir2Det(dir, strip))
                {
                    WindowsMap[nuclei][number][i] = make_pair(energy_low, energy_high);
                    // cout << "Nuclei : " << nuclei << " Number : " << number << " Detector : " << detectorName[i] << " Energy Low : " << energy_low << " Energy High : " << energy_high << endl;
                }

                if (dir == "Up" && strip == 1)
                {
                    PeakList[nuclei].push_back(number);
                }
            }
        }
    }
}


inline void InitDetectors(const string &fname)
{
  if (FLAG2021) YEAR = 2021;
  else if (FLAG2024) YEAR = 2024;
  else if (FLAG2025) YEAR = 2025;
  else
  {
    Error("Impossible to determine the year");
  }
  Info("---------  YEAR : " + to_string(YEAR) + " ---------");
  ifstream file(fname.c_str());
  if (!file.is_open())
  {
    Error("Impossible to open " + fname + "file");
  }

  string line;
  while (getline(file, line))
  {

    if (line.find("#") != string::npos)
      continue;

    istringstream iss(line);

    size_t lastDelimiterPos = line.rfind(':');
    size_t firstDelimiterPos = line.find(':');

    int label = stoi(line.substr(0, firstDelimiterPos));
    string name = line.substr(lastDelimiterPos + 1);
    
    detectorName[label] = name;

  }

  detectorName[100] = "SimSiPM";
  detectorName[99] = "HRS";


  // INIT PATHS
  if (FLAG2025)
  {
    DIR_FAST_DATA = "/run/media/local1/Disque_Dur/2025_DATA/DETECTOR_DATA/DATA/";
    DIR_ROOT_DATA = "/mnt/hgfs/shared-2/2025_DATA/DETECTOR_DATA/ROOT/";
    DIR_ROOT_DATA_GROUPED = "/mnt/hgfs/shared-2/2025_DATA/DETECTOR_DATA/GROUPED/";
    DIR_ROOT_DATA_CLEANED = "/mnt/hgfs/shared-2/2025_DATA/DETECTOR_DATA/CLEANED/";
    DIR_ROOT_DATA_MATCHED = "/mnt/hgfs/shared-2/2025_DATA/DETECTOR_DATA/MATCHED/";
    DIR_ROOT_DATA_MERGED = "/mnt/hgfs/shared-2/2025_DATA/DETECTOR_DATA/MERGED/";
    DIR_ROOT_DATA_ANALYSED = "/mnt/hgfs/shared-2/2025_DATA/DETECTOR_DATA/ANALYSED/";
    DIR_ROOT_DATA_CALIBRATED = "/mnt/hgfs/shared-2/2025_DATA/DETECTOR_DATA/CALIBRATED/";
    DIR_ROOT_DATA_MCP = "/mnt/hgfs/shared-2/2025_DATA/MCP_DATA/";
    DIR_ROOT_DATA_RATE = "/mnt/hgfs/shared-2/2025_DATA/DETECTOR_DATA/RATE/";
    DIR_ROOT_DATA_SIMULATED = "/mnt/hgfs/shared-2/2024_DATA/SIMULATED_DATA/";
    DIR_DATA_ISOLDE = "/mnt/hgfs/shared-2/2025_DATA/ISOLDE_DATA/";
    DIR_DATA_HDD = "/run/media/local1/Disque_Dur/2025_DATA/DETECTOR_DATA/";
    REFERENCE_RUN = 98;
    MATCHING_RUN = 69;
  }
  else if (FLAG2024)
  {
    DIR_FAST_DATA = "/run/media/local1/Disque_Dur/2024_DATA/DETECTOR_DATA/DATA/";
    DIR_ROOT_DATA = "/mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/ROOT/";
    DIR_ROOT_DATA_GROUPED = "/mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/GROUPED/";
    DIR_ROOT_DATA_CLEANED = "/mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/CLEANED/";
    DIR_ROOT_DATA_MATCHED = "/mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/MATCHED/";
    DIR_ROOT_DATA_MERGED = "/mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/MERGED/";
    DIR_ROOT_DATA_ANALYSED = "/mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/ANALYSED/";
    DIR_ROOT_DATA_CALIBRATED = "/mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/CALIBRATED/";
    DIR_ROOT_DATA_MCP = "/mnt/hgfs/shared-2/2024_DATA/MCP_DATA/";
    DIR_ROOT_DATA_RATE = "/mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/RATE/";
    DIR_ROOT_DATA_SIMULATED = "/mnt/hgfs/shared-2/2024_DATA/SIMULATED_DATA/";
    DIR_DATA_ISOLDE = "/mnt/hgfs/shared-2/2024_DATA/ISOLDE_DATA/";
    DIR_DATA_HDD = "/run/media/local1/Disque_Dur/2024_DATA/DETECTOR_DATA/";
    REFERENCE_RUN = 114;
    MATCHING_RUN = 77;
  }
  else if (FLAG2021)
  {
    DIR_FAST_DATA = "/run/media/local1/Disque_Dur/2021_DATA/DETECTOR_DATA/DATA/";
    DIR_ROOT_DATA = "/mnt/hgfs/shared-2/2021_DATA/DETECTOR_DATA/ROOT/";
    DIR_ROOT_DATA_GROUPED = "/mnt/hgfs/shared-2/2021_DATA/DETECTOR_DATA/GROUPED/";
    DIR_ROOT_DATA_CLEANED = "/mnt/hgfs/shared-2/2021_DATA/DETECTOR_DATA/CLEANED/";
    DIR_ROOT_DATA_MATCHED = "/mnt/hgfs/shared-2/2021_DATA/DETECTOR_DATA/MATCHED/";
    DIR_ROOT_DATA_MERGED = "/mnt/hgfs/shared-2/2021_DATA/DETECTOR_DATA/MERGED/";
    DIR_ROOT_DATA_ANALYSED = "/mnt/hgfs/shared-2/2021_DATA/DETECTOR_DATA/ANALYSED/";
    DIR_ROOT_DATA_CALIBRATED = "/mnt/hgfs/shared-2/2021_DATA/DETECTOR_DATA/CALIBRATED/";
    // DIR_ROOT_DATA_MCP = "/mnt/hgfs/shared-2/2021_DATA/MCP_DATA/";
    DIR_ROOT_DATA_RATE = "/mnt/hgfs/shared-2/2021_DATA/DETECTOR_DATA/RATE/";
    DIR_ROOT_DATA_SIMULATED = "/mnt/hgfs/shared-2/2021_DATA/SIMULATED_DATA/";
    DIR_DATA_ISOLDE = "/mnt/hgfs/shared-2/2021_DATA/ISOLDE_DATA/";
    DIR_DATA_HDD = "/run/media/local1/Disque_Dur/2021_DATA/DETECTOR_DATA/";
    REFERENCE_RUN = 16;
    MATCHING_RUN = 36;
  }

  IAS["32Ar"] = 14.;
  IAS["32Ar_thick"] = 14.;
  IAS["33Ar"] = 21.;
  IAS["33Ar_thick"] = 21.;

  // Nuclei
  NucleusColor["32Ar"] = kCyan;
  NucleusColor["32Ar_thick"] = kCyan + 2;
  NucleusColor["33Ar"] = kGreen + 1;
  NucleusColor["33Ar_thick"] = kGreen + 3;
  NucleusColor["32Cl"] = kMagenta + 1;
  NucleusColor["32Cl_thick"] = kMagenta + 3;

}



#endif