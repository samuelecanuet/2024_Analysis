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
#include "/home/local1/Documents/lib/GTools1.0/include/GString.hh"

#include <thread>
#include <chrono>

using namespace std;
namespace fs = boost::filesystem;
using namespace ROOT::Math;
random_device rd;
mt19937 gen(rd());

vector<string> NUCLEI = {"32Ar", "32Ar_thick", "33Ar"};
map<string, vector<string>> Map_RunFiles;
TFile *MERGED_File;
TFile *MATCHED_File;
string GROUPED_filename;
TFile *GROUPED_File;
TF1* Matching_function[SIGNAL_MAX];
TDirectory *dir_det[SIGNAL_MAX];
TH1D *H[SIGNAL_MAX];

vector<Signal> MERGED_Tree_Silicon;
vector<Signal> MERGED_Tree_SiPMHigh;
vector<Signal> MERGED_Tree_SiPMLow;


TTreeReaderArray<Signal> *Silicon;
TTreeReaderArray<Signal> *SiPM_High;
TTreeReaderArray<Signal> *SiPM_Low;
double Channel;

void Init()
{
  Map_RunFiles["32Ar"] = {"057", "058", "059",
                        "061", "062", "064", "065", "066", "067", "068", "069",
                        "070", "071", "072", "074", "077",
                        
                        
                        
                        "112", "113", "114", "115", "116", "118"};

  Map_RunFiles["32Ar_thick"] = {"075", "076"};

  Map_RunFiles["33Ar"] = {"078"};
}

void LoadMatchingFunction(int Run)
{
  for (int i = 0; i <= SIGNAL_MAX; i++)
  {
    if (IsDetectorSiliStrip(i))
    {
      Matching_function[i] = (TF1*)MATCHED_File->Get(("poll1" + to_string(Run) + to_string(i)).c_str()); 

      if (Matching_function[i] == NULL)
      {
        Matching_function[i] = new TF1(("poll1" + to_string(Run) + to_string(i)).c_str(), "x", 0, 10000); 
      }
    }
  }
}



#endif