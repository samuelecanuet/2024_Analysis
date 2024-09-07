#ifndef DETECTORS_HH
#define DETECTORS_HH

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream> 
#include <dirent.h>
#include <gsl/gsl_statistics.h>
#include "TFile.h"

#include "Messenger.hh"

#include "/home/local1/Documents/lib/GTools1.0/include/GString.hh"


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

/// FOLDERS ///
string DIR_FAST_DATA = "../../../../../run/media/local1/Disque_Dur/2024_DATA/DETECTOR_DATA/DATA/";
string DIR_ROOT_DATA = "../../../../../mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/ROOT/";
string DIR_ROOT_DATA_GROUPED = "../../../../../mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/GROUPED/";
string DIR_ROOT_DATA_CLEANED = "../../../../../mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/CLEANED/";
string DIR_ROOT_DATA_MATCHED = "../../../../../mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/MATCHED/";
string DIR_ROOT_DATA_MERGED  = "../../../../../mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/MERGED/";
string DIR_ROOT_DATA_SIMULATED = "../../../../../mnt/hgfs/shared-2/2024_DATA/SIMULATED_DATA/";
string DIR_DATA_ISOLDE = "../../../../../mnt/hgfs/shared-2/2024_DATA/ISOLDE_DATA/";
///////////////

string detectorFileName; ///< Detectors definition file name
size_t detectorNum = SIGNAL_MAX;      ///< Number of defined detectors channels
string detectorName[SIGNAL_MAX];
int detectorCoder[SIGNAL_MAX]; ///< Coder label of all signals
int detectorInfo[SIGNAL_MAX];  ///< Detector information of all signals (type identifier)

int detBeta[BETA_NUM][BETA_SIZE]; ///< Detector number of beta signals
int detSili[SILI_NUM][SILI_SIZE]; ///< Detector number of silicon signals

int coderDetector[FDATA_MAX]; ///< Detector number of all coders

/// Group
int mini = -300;
int maxi = 350;

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
int eSiliN_cal = 10000;

/// SiPM ///
/// High
double winHighMin = 0;
double winHighMax = 250;
int winHighN = (abs(winHighMin) + winHighMax) / 2;
double eHighMin = 0;
double eHighMax = 6000000;
int eHighN = 6000;

/// Low
double winLowMin = 0;
double winLowMax = 250;
int winLowN = (abs(winLowMin) + winLowMax) / 2;
double eLowMin = 0;
double eLowMax = 1200000;
int eLowN = 2400;

/// TOTAL WINDOW ///
double tabMIN[3] = {winSiliMin, winHighMin, winLowMin};
double winTotalMin = gsl_stats_min(tabMIN, 1, 3);
double tabMAX[3] = {winSiliMax, winHighMax, winLowMax};
double winTotalMax = gsl_stats_max(tabMAX, 1, 3);

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

inline void InitDetectors(const string &fname)
{
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
}

void WriteTime(TFile* from_File, TFile* to_file)
{
  from_File->cd();
  TObject* start = from_File->Get("Start_Time");
  TObject* stop = from_File->Get("End_Time");

  to_file->cd();
  TNamed("Start_Time", start->GetTitle()).Write();
  TNamed("End_Time", stop->GetTitle()).Write();
}

pair<string, string> GetTime(TFile* File)
{
  TObject* start = File->Get("Start_Time");
  TObject* stop = File->Get("End_Time");

  string str_start = start->GetTitle();
  string str_stop = stop->GetTitle(); 

  return make_pair(str_start, str_stop);
}

string SearchFiles(const string &directory, const string &run_number)
{
  DIR *dir;
  struct dirent *ent;
  if ((dir = opendir(directory.c_str())) != NULL)
  {
    while ((ent = readdir(dir)) != NULL)
    {
      string filename = ent->d_name;
      if (filename.find(("run_" + run_number).c_str()) == 0)
      {
        Info(("Found file: " + filename).c_str());
        closedir(dir);
        return filename;
      }
    }
    Error(("No file found for run: " + run_number).c_str());
    closedir(dir);
    return "";
  }
  else
  {
    Error(("Could not open directory: " + directory).c_str());
    return "";
  }
}

int ConvertBase5ToBase10(int base5_int) {
    string base5 = to_string(base5_int);
    int num = 0;
    for (char digit : base5) {
        num = num * 5 + (digit - '0');
    }
    return num;
}
#endif