#ifndef LIVETIME_HH
#define LIVETIME_HH

#include "../../../lib/SignalDict/Signal.h"
#include "Detectors.hh"

default_random_engine generator;

TFile *GROUPED_File;
TFile *fTime;

vector<pair<double, double>> peaks_window_F[SIGNAL_MAX];

int run_number;

pair<time_t, time_t> RunDate[1000];
TH1D* H_Data[1000];

vector<string> NUCLEI_LIST = {"18N", "32Ar", "33Ar"};
string run_number_str;

map<string, vector<double>> MAP_Delta_Time;
map<string, vector<int>> MAP_Run_Number;
map<string, double> MAP_Start;
map<string, pair<double, double>> MAP_CT;
map<string, double> MAP_LifeTime;
map<string, vector<string>> MAP_Filename;
map<string, TH1D*> H_Time_All;
map<string, TH1D*[5]> H_Time_All_Group;
map<string, TH1D*> H_Data_A;
map<string, TH1D*> H_Data_B;
map<string, TH1D*> H_Data_C;
map<string, TH1D*> H_Data_D;
map<string, TH1D*[5]> H_Decay;
map<string, TH1D*[SIGNAL_MAX]> H_Exp;
map<string, TH1D*[SIGNAL_MAX]> H_Time;
map<string, TDirectory*> dir;
map<string, TDirectory*> dir_det;
map<string, TF1*> fReleaseCurve;


TH1D* H_Channel18N[SIGNAL_MAX];
TH1D* H_Channel18N_coinc[SIGNAL_MAX];
TH1D* H_Channel18N_3a[SIGNAL_MAX];
TH1D* H_Time_All_3a;
TH1D* H_Decay_3a;

double Intensity_Threshold = 100;
double Pulse_Period = 1.2;
double Neighboor_A_Threshold = 1.5;
double Neighboor_B_Threshold = 2.7;
double Neighboor_C_Threshold = 3.9;


void Init()
{
  MAP_Filename["33Ar"] = {"run_078_multifast_33Ar_grouped.root"};
  MAP_Filename["32Ar"] = {"run_114_multifast_32Ar_grouped.root"};
  MAP_Filename["18N"] = {"run_114_multifast_32Ar_grouped.root"};

  MAP_Delta_Time["33Ar"] = {0.27};
  MAP_Delta_Time["32Ar"] = {0.62};
  MAP_Delta_Time["18N"] = {0.82};

  MAP_Run_Number["33Ar"] = {78};
  MAP_Run_Number["32Ar"] = {114};
  MAP_Run_Number["18N"] = {114};

  MAP_Start["33Ar"] = 0.;
  MAP_Start["32Ar"] = 0.;
  MAP_Start["18N"] = 0.;

  MAP_CT["33Ar"] = make_pair(44, 50);
  MAP_CT["32Ar"] = make_pair(44, 50);
  MAP_CT["18N"] = make_pair(44, 50);

  MAP_LifeTime["33Ar"] = 0.173;
  MAP_LifeTime["32Ar"] = 0.098;
  MAP_LifeTime["18N"] = 0.619;

  H_Time_All_3a = (TH1D *)H_Data[MAP_Run_Number["18N"][0]]->Clone("H_Time_All_3a");
  H_Time_All_3a->Reset();
  H_Decay_3a = new TH1D("H_Decay_3a", "H_Decay_3a", Pulse_Period*4*1000, 0, Pulse_Period*4);

  for (auto NUCLEUS : NUCLEI_LIST)
  {

    dir[NUCLEUS] = fTime->mkdir(NUCLEUS.c_str());
    dir_det[NUCLEUS] = dir[NUCLEUS]->mkdir("Detectors");

    H_Time_All[NUCLEUS] = (TH1D *)H_Data[MAP_Run_Number[NUCLEUS][0]]->Clone(("H_Time_All_" + NUCLEUS).c_str());
    H_Time_All[NUCLEUS]->Reset();

    H_Time_All_Group[NUCLEUS][0] = (TH1D *)H_Time_All[NUCLEUS]->Clone(("H_Time_All_Group_" + NUCLEUS + "_0").c_str());
    H_Time_All_Group[NUCLEUS][1] = (TH1D *)H_Time_All[NUCLEUS]->Clone(("H_Time_All_Group_" + NUCLEUS + "_A").c_str());
    H_Time_All_Group[NUCLEUS][2] = (TH1D *)H_Time_All[NUCLEUS]->Clone(("H_Time_All_Group_" + NUCLEUS + "_B").c_str());
    H_Time_All_Group[NUCLEUS][3] = (TH1D *)H_Time_All[NUCLEUS]->Clone(("H_Time_All_Group_" + NUCLEUS + "_C").c_str());
    H_Time_All_Group[NUCLEUS][4] = (TH1D *)H_Time_All[NUCLEUS]->Clone(("H_Time_All_Group_" + NUCLEUS + "_D").c_str());

    H_Data_A[NUCLEUS] = (TH1D *)H_Time_All[NUCLEUS]->Clone(("H_Data_A_" + NUCLEUS).c_str());
    H_Data_B[NUCLEUS] = (TH1D *)H_Time_All[NUCLEUS]->Clone(("H_Data_B_" + NUCLEUS).c_str());
    H_Data_C[NUCLEUS] = (TH1D *)H_Time_All[NUCLEUS]->Clone(("H_Data_C_" + NUCLEUS).c_str());
    H_Data_D[NUCLEUS] = (TH1D *)H_Time_All[NUCLEUS]->Clone(("H_Data_D_" + NUCLEUS).c_str());

    H_Decay[NUCLEUS][0] = new TH1D(("H_Decay_" + NUCLEUS + "_0").c_str(), ("H_Decay_" + NUCLEUS + "_0").c_str(), Pulse_Period*4*1000, 0, Pulse_Period*4);
    H_Decay[NUCLEUS][1] = new TH1D(("H_Decay_" + NUCLEUS + "_A").c_str(), ("H_Decay_" + NUCLEUS + "_A").c_str(), Pulse_Period*4*1000, 0, Pulse_Period*4);
    H_Decay[NUCLEUS][2] = new TH1D(("H_Decay_" + NUCLEUS + "_B").c_str(), ("H_Decay_" + NUCLEUS + "_B").c_str(), Pulse_Period*4*1000, 0, Pulse_Period*4);
    H_Decay[NUCLEUS][3] = new TH1D(("H_Decay_" + NUCLEUS + "_C").c_str(), ("H_Decay_" + NUCLEUS + "_C").c_str(), Pulse_Period*4*1000, 0, Pulse_Period*4);
    H_Decay[NUCLEUS][4] = new TH1D(("H_Decay_" + NUCLEUS + "_D").c_str(), ("H_Decay_" + NUCLEUS + "_D").c_str(), Pulse_Period*4*1000, 0, Pulse_Period*4);

  
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
      if (!IsDetectorSiliStrip(i))
        continue;
      
      H_Exp[NUCLEUS][i] = new TH1D(("H_Exp_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_Exp_" + NUCLEUS + "_" + detectorName[i]).c_str(), 10000, 0, 100000);
      H_Exp[NUCLEUS][i]->GetXaxis()->SetTitle("Channel");
      H_Exp[NUCLEUS][i]->GetYaxis()->SetTitle("Counts");
      H_Exp[NUCLEUS][i]->GetXaxis()->CenterTitle();
      H_Exp[NUCLEUS][i]->GetYaxis()->CenterTitle();

      H_Time[NUCLEUS][i] = (TH1D *)H_Data[MAP_Run_Number[NUCLEUS][0]]->Clone(("H_Time_" + NUCLEUS + "_" + detectorName[i]).c_str());
      H_Time[NUCLEUS][i]->Reset();
      H_Time[NUCLEUS][i]->GetXaxis()->SetTitle("Time [s]");
      H_Time[NUCLEUS][i]->GetYaxis()->SetTitle("Counts / ms");
      H_Time[NUCLEUS][i]->GetXaxis()->CenterTitle();
      H_Time[NUCLEUS][i]->GetYaxis()->CenterTitle();

      H_Channel18N[i] = new TH1D(("H_Channel18N_" + detectorName[i]).c_str(), ("H_Channel18N_" + detectorName[i]).c_str(), 10000, 0, 100000);
      H_Channel18N[i]->GetXaxis()->SetTitle("Channel");
      H_Channel18N[i]->GetYaxis()->SetTitle("Counts");
      H_Channel18N[i]->GetXaxis()->CenterTitle();
      H_Channel18N[i]->GetYaxis()->CenterTitle();

      H_Channel18N_coinc[i] = new TH1D(("H_Channel18N_coinc_" + detectorName[i]).c_str(), ("H_Channel18N_coinc_" + detectorName[i]).c_str(), 10000, 0, 100000);
      H_Channel18N_coinc[i]->GetXaxis()->SetTitle("Channel");
      H_Channel18N_coinc[i]->GetYaxis()->SetTitle("Counts");
      H_Channel18N_coinc[i]->GetXaxis()->CenterTitle();
      H_Channel18N_coinc[i]->GetYaxis()->CenterTitle();

      H_Channel18N_3a[i] = new TH1D(("H_Channel18N_3a_" + detectorName[i]).c_str(), ("H_Channel18N_3a_" + detectorName[i]).c_str(), 10000, 0, 100000);
      H_Channel18N_3a[i]->GetXaxis()->SetTitle("Channel");
      H_Channel18N_3a[i]->GetYaxis()->SetTitle("Counts");
      H_Channel18N_3a[i]->GetXaxis()->CenterTitle();
      H_Channel18N_3a[i]->GetYaxis()->CenterTitle();

    }
  }
}

void WriteHistograms()
{
  fTime->cd();
  for (auto NUCLEUS : NUCLEI_LIST)
  {
    dir[NUCLEUS]->cd();
    H_Time_All[NUCLEUS]->Write();
    H_Data[MAP_Run_Number[NUCLEUS][0]]->Write();
    H_Data_A[NUCLEUS]->Write();
    H_Data_B[NUCLEUS]->Write();
    H_Data_C[NUCLEUS]->Write();
    H_Data_D[NUCLEUS]->Write();

    TCanvas *c = new TCanvas(("ISOLDE_DETECTOR_" + NUCLEUS).c_str(), ("ISOLDE_DETECTOR_" + NUCLEUS).c_str(), 800, 600);
    H_Time_All[NUCLEUS]->Draw("HIST");
    H_Data[MAP_Run_Number[NUCLEUS][0]]->SetLineColor(kRed);
    H_Data[MAP_Run_Number[NUCLEUS][0]]->Draw("HIST SAME");
    c->Write();

    TCanvas *c_CrossTalk = new TCanvas(("CrossTalk_" + NUCLEUS).c_str(), ("CrossTalk_" + NUCLEUS).c_str(), 800, 600);
    TH1D* H_CT = (TH1D *)H_Time_All[NUCLEUS]->Clone(("H_CT_" + NUCLEUS).c_str());
    H_CT->Rebin(10);
    H_CT->GetXaxis()->SetRangeUser(MAP_CT[NUCLEUS].first, MAP_CT[NUCLEUS].second);
    H_CT->Draw("HIST");
    c_CrossTalk->Write();


    //// A
    TCanvas *cA = new TCanvas(("cA_" + NUCLEUS).c_str(), ("cA_" + NUCLEUS).c_str(), 800, 600);
    H_Data[MAP_Run_Number[NUCLEUS][0]]->Draw("HIST");
    H_Data_A[NUCLEUS]->SetLineColor(kBlue);
    H_Data_A[NUCLEUS]->Draw("HIST SAME");
    cA->Write();

    TCanvas *cA_Time = new TCanvas(("cA_Time_" + NUCLEUS).c_str(), ("cA_Time_" + NUCLEUS).c_str(), 800, 600);
    H_Time_All_Group[NUCLEUS][1]->Draw("HIST");
    cA_Time->Write();

    TCanvas *cA_Decay = new TCanvas(("cA_Decay_" + NUCLEUS).c_str(), ("cA_Decay_" + NUCLEUS).c_str(), 800, 600);
    H_Decay[NUCLEUS][1]->Draw("HIST");
    cA_Decay->Write();

    //// B
    TCanvas *cB = new TCanvas(("cB_" + NUCLEUS).c_str(), ("cB_" + NUCLEUS).c_str(), 800, 600);
    H_Data[MAP_Run_Number[NUCLEUS][0]]->Draw("HIST");
    H_Data_B[NUCLEUS]->SetLineColor(kRed);
    H_Data_B[NUCLEUS]->Draw("HIST SAME");
    cB->Write();

    TCanvas *cB_Time = new TCanvas(("cB_Time_" + NUCLEUS).c_str(), ("cB_Time_" + NUCLEUS).c_str(), 800, 600);
    H_Time_All_Group[NUCLEUS][2]->Draw("HIST");
    cB_Time->Write();

    TCanvas *cB_Decay = new TCanvas(("cB_Decay_" + NUCLEUS).c_str(), ("cB_Decay_" + NUCLEUS).c_str(), 800, 600);
    H_Decay[NUCLEUS][2]->Draw("HIST");
    cB_Decay->Write();
    
    //// C
    TCanvas *cC = new TCanvas(("cC_" + NUCLEUS).c_str(), ("cC_" + NUCLEUS).c_str(), 800, 600);
    H_Data[MAP_Run_Number[NUCLEUS][0]]->Draw("HIST");
    H_Data_C[NUCLEUS]->SetLineColor(kGreen);
    H_Data_C[NUCLEUS]->Draw("HIST SAME");
    cC->Write();

    TCanvas *cC_Time = new TCanvas(("cC_Time_" + NUCLEUS).c_str(), ("cC_Time_" + NUCLEUS).c_str(), 800, 600);
    H_Time_All_Group[NUCLEUS][3]->Draw("HIST");
    cC_Time->Write();

    TCanvas *cC_Decay = new TCanvas(("cC_Decay_" + NUCLEUS).c_str(), ("cC_Decay_" + NUCLEUS).c_str(), 800, 600);
    H_Decay[NUCLEUS][3]->Draw("HIST");
    cC_Decay->Write();

    //// D
    TCanvas *cD = new TCanvas(("cD_" + NUCLEUS).c_str(), ("cD_" + NUCLEUS).c_str(), 800, 600);
    H_Data[MAP_Run_Number[NUCLEUS][0]]->Draw("HIST");
    H_Data_D[NUCLEUS]->SetLineColor(kBlack);
    H_Data_D[NUCLEUS]->Draw("HIST SAME");
    cD->Write();

    TCanvas *cD_Time = new TCanvas(("cD_Time_" + NUCLEUS).c_str(), ("cD_Time_" + NUCLEUS).c_str(), 800, 600);
    H_Time_All_Group[NUCLEUS][4]->Draw("HIST");
    cD_Time->Write();

    TCanvas *cD_Decay = new TCanvas(("cD_Decay_" + NUCLEUS).c_str(), ("cD_Decay_" + NUCLEUS).c_str(), 800, 600);
    H_Decay[NUCLEUS][4]->Draw("HIST");
    cD_Decay->Write();

    TCanvas *cABCD = new TCanvas(("cABCD_" + NUCLEUS).c_str(), ("cABCD_" + NUCLEUS).c_str(), 800, 600);
    H_Data[MAP_Run_Number[NUCLEUS][0]]->Draw("HIST");
    H_Data_A[NUCLEUS]->SetLineColor(kBlue);
    H_Data_A[NUCLEUS]->Draw("HIST SAME");
    H_Data_B[NUCLEUS]->SetLineColor(kRed);
    H_Data_B[NUCLEUS]->Draw("HIST SAME");
    H_Data_C[NUCLEUS]->SetLineColor(kGreen);
    H_Data_C[NUCLEUS]->Draw("HIST SAME");
    H_Data_D[NUCLEUS]->SetLineColor(kBlack);
    H_Data_D[NUCLEUS]->Draw("HIST SAME");
    cABCD->Write();

    TCanvas *cABCD_Time = new TCanvas(("cABCD_Time_" + NUCLEUS).c_str(), ("cABCD_Time_" + NUCLEUS).c_str(), 800, 600);  
    H_Time_All_Group[NUCLEUS][1]->SetLineColor(kBlue);
    H_Time_All_Group[NUCLEUS][1]->Draw("HIST SAME");
    H_Time_All_Group[NUCLEUS][2]->SetLineColor(kRed);
    H_Time_All_Group[NUCLEUS][2]->Draw("HIST SAME");
    H_Time_All_Group[NUCLEUS][3]->SetLineColor(kGreen);
    H_Time_All_Group[NUCLEUS][3]->Draw("HIST SAME");

    cABCD_Time->Write();

    TCanvas *cABCD_Decay = new TCanvas(("cABCD_Decay_" + NUCLEUS).c_str(), ("cABCD_Decay_" + NUCLEUS).c_str(), 800, 600);
    TLegend *lABCD_Decay = new TLegend(0.1, 0.7, 0.3, 0.9);
    H_Decay[NUCLEUS][0]->GetXaxis()->SetRangeUser(0, Pulse_Period);
    H_Decay[NUCLEUS][1]->GetXaxis()->SetRangeUser(0, Pulse_Period);
    H_Decay[NUCLEUS][2]->GetXaxis()->SetRangeUser(0, Pulse_Period);
    H_Decay[NUCLEUS][3]->GetXaxis()->SetRangeUser(0, Pulse_Period);
    H_Decay[NUCLEUS][4]->GetXaxis()->SetRangeUser(0, Pulse_Period);


    vector<double> integrals = {H_Decay[NUCLEUS][0]->Integral(), H_Decay[NUCLEUS][1]->Integral(), H_Decay[NUCLEUS][2]->Integral(), H_Decay[NUCLEUS][3]->Integral(), H_Decay[NUCLEUS][4]->Integral()};
    auto integral = *min_element(integrals.begin(), integrals.end());
    H_Decay[NUCLEUS][0]->Scale(integral/H_Decay[NUCLEUS][0]->Integral());
    H_Decay[NUCLEUS][1]->Scale(integral/H_Decay[NUCLEUS][1]->Integral());
    H_Decay[NUCLEUS][2]->Scale(integral/H_Decay[NUCLEUS][2]->Integral());
    H_Decay[NUCLEUS][3]->Scale(integral/H_Decay[NUCLEUS][3]->Integral());
    H_Decay[NUCLEUS][4]->Scale(integral/H_Decay[NUCLEUS][4]->Integral());

    H_Decay[NUCLEUS][0]->GetXaxis()->SetRangeUser(-1111, -1111);
    H_Decay[NUCLEUS][1]->GetXaxis()->SetRangeUser(-1111, -1111);
    H_Decay[NUCLEUS][2]->GetXaxis()->SetRangeUser(-1111, -1111);
    H_Decay[NUCLEUS][3]->GetXaxis()->SetRangeUser(-1111, -1111);
    H_Decay[NUCLEUS][4]->GetXaxis()->SetRangeUser(-1111, -1111);

    H_Decay[NUCLEUS][1]->SetLineColor(kBlue);
    H_Decay[NUCLEUS][1]->Draw("HIST");
    H_Decay[NUCLEUS][1]->GetXaxis()->SetTitle("Time [s]");
    lABCD_Decay->AddEntry(H_Decay[NUCLEUS][1], "1.2s", "l");
    H_Decay[NUCLEUS][2]->SetLineColor(kRed);
    H_Decay[NUCLEUS][2]->Draw("HIST SAME");
    lABCD_Decay->AddEntry(H_Decay[NUCLEUS][2], "2.4s", "l");
    H_Decay[NUCLEUS][3]->SetLineColor(kGreen);
    H_Decay[NUCLEUS][3]->Draw("HIST SAME");
    lABCD_Decay->AddEntry(H_Decay[NUCLEUS][3], "3.6s", "l");
    H_Decay[NUCLEUS][4]->SetLineColor(kBlack);
    H_Decay[NUCLEUS][4]->Draw("HIST SAME");
    lABCD_Decay->AddEntry(H_Decay[NUCLEUS][4], "4.8s", "l");
    lABCD_Decay->Draw("SAME");
    cABCD_Decay->Write();

    

    dir_det[NUCLEUS]->cd();
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
      if (!IsDetectorSiliStrip(i))
        continue;
      
      H_Exp[NUCLEUS][i]->Write();
      H_Time[NUCLEUS][i]->Write();

      if (NUCLEUS == "18N")
      {
        TCanvas *c = new TCanvas(("Channel_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("ISOLDE_DETECTOR_" + NUCLEUS + "_" + detectorName[i]).c_str(), 800, 600);
        H_Channel18N[i]->Draw("HIST");
        H_Channel18N_coinc[i]->SetLineColor(kRed);
        H_Channel18N_coinc[i]->Draw("HIST SAME");
        c->Write();

        H_Channel18N_3a[i]->Write();

      }
    }
  }


  //////////// PLOTTING DECAY NUCLEI COMPARAISON WITH PULSE PERIOD ////////////
  fTime->cd();
  TCanvas *cA_all_Decay = new TCanvas("cA_all_Decay", "cA_all_Decay", 800, 600);
  TLegend *lA_all_Decay = new TLegend(0.1, 0.7, 0.3, 0.9);
  TCanvas *cB_all_Decay = new TCanvas("cB_all_Decay", "cB_all_Decay", 800, 600);
  TLegend *lB_all_Decay = new TLegend(0.1, 0.7, 0.3, 0.9);
  TCanvas *cC_all_Decay = new TCanvas("cC_all_Decay", "cC_all_Decay", 800, 600);
  TLegend *lC_all_Decay = new TLegend(0.1, 0.7, 0.3, 0.9);


  //rescale
  vector<double> integrals1;
  vector<double> integrals2;
  vector<double> integrals3;

  int bin = 0;
  for (auto NUCLEUS : NUCLEI_LIST)
  {
    H_Decay[NUCLEUS][1]->GetXaxis()->SetRangeUser(0, 0.1);
    integrals1.push_back(H_Decay[NUCLEUS][1]->Integral());
    H_Decay[NUCLEUS][2]->GetXaxis()->SetRangeUser(0, 0.1);
    integrals2.push_back(H_Decay[NUCLEUS][2]->Integral());
    H_Decay[NUCLEUS][3]->GetXaxis()->SetRangeUser(0, 0.1);
    integrals3.push_back(H_Decay[NUCLEUS][3]->Integral());
  }

  double integral1 = *min_element(integrals1.begin(), integrals1.end());
  double integral2 = *min_element(integrals2.begin(), integrals2.end());
  double integral3 = *min_element(integrals3.begin(), integrals3.end());

  int counter = 0;
  for (auto NUCLEUS : NUCLEI_LIST)
  {
    
    H_Decay[NUCLEUS][1]->Scale(integral1/integrals1[counter]);
    H_Decay[NUCLEUS][1]->GetXaxis()->SetRangeUser(-1111, -1111);
    H_Decay[NUCLEUS][2]->Scale(integral2/integrals2[counter]);
    H_Decay[NUCLEUS][2]->GetXaxis()->SetRangeUser(-1111, -1111);
    H_Decay[NUCLEUS][3]->Scale(integral3/integrals3[counter]);
    H_Decay[NUCLEUS][3]->GetXaxis()->SetRangeUser(-1111, -1111);
    counter++;
  }

  //plotting
  counter = 0;
  for (auto NUCLEUS : NUCLEI_LIST)
  {
    counter++ ;
    cA_all_Decay->cd();
    lA_all_Decay->AddEntry(H_Decay[NUCLEUS][1], NUCLEUS.c_str(), "l");
    H_Decay[NUCLEUS][1]->SetLineColor(kBlue+counter);
    H_Decay[NUCLEUS][1]->Draw("HIST SAME");

    cB_all_Decay->cd();
    lB_all_Decay->AddEntry(H_Decay[NUCLEUS][2], NUCLEUS.c_str(), "l");
    H_Decay[NUCLEUS][2]->SetLineColor(kBlue+counter);
    H_Decay[NUCLEUS][2]->Draw("HIST SAME");

    cC_all_Decay->cd();
    lC_all_Decay->AddEntry(H_Decay[NUCLEUS][3], NUCLEUS.c_str(), "l");
    H_Decay[NUCLEUS][3]->SetLineColor(kBlue+counter);
    H_Decay[NUCLEUS][3]->Draw("HIST SAME");
  }

  cA_all_Decay->cd();
  lA_all_Decay->Draw("SAME");
  cA_all_Decay->Write();
  cB_all_Decay->cd();
  lB_all_Decay->Draw("SAME");
  cB_all_Decay->Write();
  cC_all_Decay->cd();
  lC_all_Decay->Draw("SAME");
  cC_all_Decay->Write();

  
}

void InitPeakWindow(string nuclei)
{
  string CalibFileName;

  CalibFileName = "./Config_Files/Win_"+nuclei+"_Catcher1_14.txt";
  
  
  ifstream fileF(CalibFileName);

  // for (auto &value : peaks_window_F)
  // {
  //   value = make_pair(0, 0);
  // }

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
      iss >> DetName >> min >> max >> min1 >> max1;

      for (size_t i = 0; i < detectorNum; ++i)
      {
        if (DetName == detectorName[i])
        {
          peaks_window_F[i].clear();
          peaks_window_F[i].push_back(make_pair(min, max));

          if (min1 != 0 && max1 != 0)
          {
            peaks_window_F[i].push_back(make_pair(min1, max1));
          }
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

void ReadAllRunsDate()
{

  for (int i = 0; i < 1000; i++)
  {
    RunDate[i] = make_pair(0, 0);
  }

  DIR *dir;
  struct dirent *ent;
  if ((dir = opendir(DIR_ROOT_DATA.c_str())) != NULL)
  {
    while ((ent = readdir(dir)) != NULL)
    {
      string filename = ent->d_name;
      if (filename.find("run_") == 0)
      {
        string run_number = filename.substr(4, 3);
        int i = stoi(run_number);
        TFile *f = new TFile((DIR_ROOT_DATA + filename).c_str(), "READ");
        if (!f->IsOpen())
        {
          Error("Impossible to open ROOT file");
        }
        TObject *start = f->Get("Start_Time");
        TObject *stop = f->Get("End_Time");
        if (!start || !stop)
        {
          Warning("Impossible to find Start_Time or End_Time for run " + run_number);
        }

        else
        {
          Info("Found for run " + run_number);
          string start_time_str = start->GetTitle();
          string stop_time_str = stop->GetTitle();

          struct tm tm_start = {};
          struct tm tm_stop = {};

          istringstream ss1(start_time_str);
          ss1 >> std::get_time(&tm_start, "%y-%m-%d at %Hh-%Mm-%Ss");
          tm_start.tm_year += 100;

          istringstream ss2(stop_time_str);
          ss2 >> std::get_time(&tm_stop, "%y-%m-%d at %Hh-%Mm-%Ss");
          tm_stop.tm_year += 100;


          time_t start_time_t = mktime(&tm_start);
          time_t stop_time_t = mktime(&tm_stop);
          int start_time = static_cast<int>(start_time_t);
          int stop_time = static_cast<int>(stop_time_t);  
          double diff= stop_time-start_time;

          H_Data[i] = new TH1D(Form("H_Data_run_%d", i), Form("H_Data_run_%d", i), diff*1000, 0, diff);
          H_Data[i]->GetXaxis()->SetTitle("Time [s]");
          H_Data[i]->GetYaxis()->SetTitle("Intensity");
          H_Data[i]->GetXaxis()->CenterTitle();
          H_Data[i]->GetYaxis()->CenterTitle();
          RunDate[i] = make_pair(start_time_t, stop_time_t);
        }
      }
    }
    closedir(dir);
  }
  else
  {
    Error("Impossible to open ROOT data directory");
  }
}

void ReadISOLDE()
{
  string ISOLDE_filename = DIR_DATA_ISOLDE + "2024_05_ISOLDE_data.csv";
  ifstream file(ISOLDE_filename.c_str());

  if (!file.is_open())
  {
    Error("Impossible to open ISOLDE file");
  }

  string line;
  int counter = 0;
  int first = 56;
  double milliseconds;
  while (getline(file, line))
  {
    counter++;
    if (counter < 4)
      continue;
    
    ProgressBar(counter, 2000000, 0, 0, "ISOLDE data : ");

    istringstream iss(line);

    string data_time, int_str;
    getline(iss, data_time, ',');
    getline(iss, int_str);

    double intensity = 0.0;
    try
    {
      intensity = std::stod(int_str);
    }
    catch (const std::invalid_argument &ia)
    {
      // Warning(("Invalid argument"));
      continue;
    }
    catch (const std::out_of_range &oor)
    {
      // Warning(("Out of range"));
      continue;
    }

    std::tm tm_data = {};

    std::istringstream ss1(data_time);
    ss1 >> std::get_time(&tm_data, "%Y-%m-%d %H:%M:%S"); 
    ss1.ignore();
    ss1 >> milliseconds;

    // cout << "Data : " << mktime(&tm_data) << "    " << data_time << "    " << milliseconds << endl;

    for (int i = first; i < 130; i++)
    {
      if (RunDate[i].first == 0)
        continue;


      // // cout all the date 
      // cout << "Run : " << i << " Start : " << RunDate[i].first << " Stop : " << RunDate[i].second << endl;
      // // cout tm_data
      // cout << "Data : " << mktime(&tm_data) << endl;

      
      if (i != 78)
        tm_data.tm_hour += 2;

      if (mktime(&tm_data) > RunDate[i].first && mktime(&tm_data) < RunDate[i].second)
      {        
        double time = mktime(&tm_data)-RunDate[i].first+milliseconds/1000;
        H_Data[i]->SetBinContent(H_Data[i]->FindBin(time), intensity);
        first = i;
        milliseconds = 0;
        break;
      }
    }
  }
}

double expo(double *x, double *par)
{
  double lambda = par[0];
  double sigma = par[1];
  double T0 = par[2];
  double t0 = par[3];

  return sqrt(2*M_PI)*sigma*exp(lambda/2*(lambda*sigma*sigma+2*T0-2*x[0]+2*t0));
}

double MyReleaseCurve(double *x, double *par)
{
  double lambda1 = par[0];
  double lambda2 = par[1];
  double lambda3 = par[2];
  double alpha = par[3];
  double t0 = par[4];
  double A = par[5];
  double sigma_T = par[6];
  double T0 = par[7];

  double parC1[4] = {lambda2, sigma_T, T0, t0};
  double C1 = alpha * expo(x, par);

  double parC2[4] = {lambda3, sigma_T, T0, t0};
  double C2 = (1 - alpha) * expo(x, par);

  double parC3[4] = {lambda1+lambda2, sigma_T, T0, t0};
  double C3 = - alpha * expo(x, par);

  double parC4[4] = {lambda1+lambda3, sigma_T, T0, t0};
  double C4 = - (1 - alpha) * expo(x, par);

  return 1/A * (C1 + C2 + C3 + C4);
}

double ReleaseCurve(double *x, double *par)
{
  double A = par[0];
  double alpha = par[1];
  double lambda1 = log(2)/par[2];
  double lambda2 = log(2)/par[3];
  double lambda3 = log(2)/par[4];
  double t0 = par[5];

  double main = 1/A * (1-exp(-lambda1*(x[0]-t0))) *(alpha * exp(-lambda2*(x[0]-t0)) + (1-alpha) * exp(-lambda3*(x[0]-t0)));
  t0 -= 2.4;
  double pre = 1/A * (1-exp(-lambda1*(x[0]-t0))) *(alpha * exp(-lambda2*(x[0]-t0)) + (1-alpha) * exp(-lambda3*(x[0]-t0)));

  return main+pre;
}

int Group_Determining(int i, string nucleus)
{
  if (H_Data_A[nucleus]->GetBinContent(i) > 0)
  {
    return 1;
  }
  else if (H_Data_B[nucleus]->GetBinContent(i) > 0)
  {
    return 2;
  }
  else if (H_Data_C[nucleus]->GetBinContent(i) > 0)
  {
    return 3;
  }
  else if (H_Data_D[nucleus]->GetBinContent(i) > 0)
  {
    return 4;
  }
  else
  {
    return 0;
  }
}

#endif