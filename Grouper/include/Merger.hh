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

vector<string> NUCLEI = {"32Ar", "32Ar_thick", "33Ar"};
map<string, vector<string>> Map_RunFiles;
TFile *MERGED_File;
TFile *MATCHED_File;
string GROUPED_filename;
TFile *GROUPED_File;
TF1* Matching_function[SIGNAL_MAX];
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
// vector<Signal> MERGED_Tree_SiPMHigh;
// vector<Signal> MERGED_Tree_SiPMLow;
vector<vector<pair<Signal, Signal>>> MERGED_Tree_SiPMGroup;


TTreeReaderArray<Signal> *Silicon;
// TTreeReaderArray<Signal> *SiPM_High;
// TTreeReaderArray<Signal> *SiPM_Low;
TTreeReaderValue<vector<vector<pair<Signal, Signal>>>> *SiPM_Groups;
double Channel;

TH1D* H_Sum;

void Init()
{
  Map_RunFiles["32Ar"] = {"057", "058", "059",
                         "061", "062", "064", "065", "066", "067", "068", "069",
                         "070", "071", "072", "074", "077",
                        
                        
                        
                         "112", "113", "114", "115", "116", "118"};

  Map_RunFiles["32Ar_thick"] = {"075", "076"};

  Map_RunFiles["33Ar"] = {"078"};
}

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

void LoadMatchingFunction(int Run)
{
  for (int i = 0; i <= SIGNAL_MAX; i++)
  {
    if (IsDetectorSili(i))
    {
      Matching_function[i] = (TF1*)MATCHED_File->Get(("poll1" + to_string(Run) + to_string(i)).c_str()); 

      if (Matching_function[i] == NULL)
      {
        Matching_function[i] = new TF1(("poll1" + to_string(Run) + to_string(i)).c_str(), "x", 0, 10000); 
      }
    }
  }
}

void ComputeFakeCoincidences()
{
  for (int mul = 1; mul <= BETA_SIZE; mul++)
  {

    H_SiPM_LeftLeft[mul] = (TH1D*)H_SiPM_Full[mul]->Clone(("H_SiPM_LeftLeft_" + to_string(mul)).c_str());
    for (int bin = 0; bin < H_SiPM_LeftLeft[mul]->GetNbinsX(); bin++)
    {
      if (H_SiPM_LeftLeft[mul]->GetBinCenter(bin) > -240 || H_SiPM_LeftLeft[mul]->GetBinCenter(bin) < -300)
      {
        H_SiPM_LeftLeft[mul]->SetBinContent(bin, 0);
      }
    }

    double integral_left_left = H_SiPM_LeftLeft[mul]->Integral();
    double mean_height_left_left = integral_left_left/30;

    H_SiPM_Center[mul] = (TH1D*)H_SiPM_Full[mul]->Clone(("H_SiPM_Center_" + to_string(mul)).c_str());
    for (int bin = 0; bin < H_SiPM_Center[mul]->GetNbinsX(); bin++)
    {
      if (H_SiPM_Center[mul]->GetBinCenter(bin) > 230 || H_SiPM_Center[mul]->GetBinCenter(bin) < -240)
      {
        H_SiPM_Center[mul]->SetBinContent(bin, 0);
      }
    }

    TH1D* H_Fortuitous = (TH1D*)H_SiPM_Full[mul]->Clone("_Fortuitous_");
    H_Fortuitous->Reset();
    for (int bin = 0; bin < H_SiPM_Center[mul]->GetNbinsX(); bin++)
    {
      if (H_SiPM_Center[mul]->GetBinCenter(bin) > -240 && H_SiPM_Center[mul]->GetBinCenter(bin) < 130)
      {
        H_Fortuitous->SetBinContent(bin, mean_height_left_left);
      }
      if (H_SiPM_Center[mul]->GetBinCenter(bin) > 130 || H_SiPM_Center[mul]->GetBinCenter(bin) < -240)
      {
        H_Fortuitous->SetBinContent(bin, H_SiPM_Full[mul]->GetBinContent(bin));
      }
    }

    TH1D* H_Real = (TH1D*)H_SiPM_Center[mul]->Clone("_Real_");
    for (int bin = 0; bin < H_SiPM_Center[mul]->GetNbinsX(); bin++)
    {
      H_Real->SetBinContent(bin, H_SiPM_Center[mul]->GetBinContent(bin) - mean_height_left_left );
    }


  
    double integral_center = H_SiPM_Center[mul]->Integral();


    TCanvas *C_FakeCoincidence = new TCanvas(("C_Fake_" + to_string(mul)).c_str(), ("C_Fake_" + to_string(mul)).c_str(), 800, 400);
    H_SiPM_Full[mul]->SetFillColor(kBlack);
    H_SiPM_Full[mul]->SetFillStyle(3244);
    H_SiPM_Full[mul]->Draw("HIST");
    H_Fortuitous->SetFillStyle(3244);
    H_Fortuitous->SetFillColor(kRed);
    H_Fortuitous->Draw("SAME");

    TLatex *text_low = new TLatex();
    text_low->SetNDC();
    text_low->SetTextSize(0.03);
    H_SiPM_Full[mul]->GetXaxis()->SetRangeUser(-240, 130);
    double pourcentage = mean_height_left_left/2*370 * 100 / H_SiPM_Full[mul]->Integral();
    H_SiPM_Full[mul]->GetXaxis()->SetRangeUser(-1111, -1111);
    cout << "Real : " << H_Real->Integral() << "    Fake : " << mean_height_left_left/2*370 << endl;
    text_low->DrawLatex(0.15, 0.8, ("Fake coincidences: " + to_string(pourcentage) + " %").c_str());
    text_low->Draw("SAME");
    C_FakeCoincidence->Write(); 

    TCanvas *C_Windows = new TCanvas(("C_Windows_" + to_string(mul)).c_str(), ("C_Windows_" + to_string(mul)).c_str(), 800, 400);
    H_SiPM_Center[mul]->Draw("HIST");
    H_SiPM_Left[mul]->SetLineColor(kRed);
    H_SiPM_Left[mul]->Draw("SAME");
    H_SiPM_Right[mul]->SetLineColor(kRed);
    H_SiPM_Right[mul]->Draw("SAME");
    C_Windows->Write();

    G_Fake->AddPoint(mul, pourcentage);
  }

  TCanvas *C_Fake = new TCanvas("C_Fake", "C_Fake", 800, 400);
  G_Fake->GetXaxis()->SetTitle("Multiplicity");
  G_Fake->GetYaxis()->SetTitle("Fake coincidences (%)");
  G_Fake->SetMarkerStyle(20);
  G_Fake->SetMarkerSize(1);
  G_Fake->Draw("AP");
  C_Fake->Write();
}





#endif