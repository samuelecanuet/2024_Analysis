#ifndef PULSE_HH
#define PULSE_HH

#include "Detectors.hh"

default_random_engine generator;

TFile *GROUPED_File;
TFile *fTime;

TFile *fISOLDE;
map<string, TTree*> Tree;

int run_number;
int VERBOSE = 1;

map<string, pair<time_t, time_t>> RunDate;
map<string, TH1D*> H_Data;

TF1* Calibration[SIGNAL_MAX];
int start_gate = -10;
int end_gate = 40;

map<string, vector<double>> MAP_Delta_Time;
map<string, vector<int>> MAP_Run_Number;
map<string, double> MAP_Start;
map<string, pair<double, double>> MAP_CT;
map<string, double> MAP_LifeTime;
map<string, vector<string>> MAP_Filename;
map<string, TH1D*> H_Time_All;
map<string, TH1D*> H_Time_All_Coinc;
map<string, TH1D*[5]> H_Time_All_Group;
map<string, TH1D*> H_Data_A;
map<string, TH1D*> H_Data_B;
map<string, TH1D*> H_Data_C;
map<string, TH1D*> H_Data_D;
map<string, TH1D*[5]> H_Decay;
map<string, TH1D*[5]> H_Decay_All;
map<string, TH1D*[5]> H_Decay_Coinc;
map<string, TH1D*[SIGNAL_MAX]> H_Exp;
map<string, TH1D*[SIGNAL_MAX]> H_Exp_Coinc;
map<string, TH1D*[SIGNAL_MAX]> H_Time;
map<string, TDirectory*> dir;
map<string, TDirectory*[SIGNAL_MAX]> dir_det;
map<string, TF1*> fReleaseCurve;

map<string, TGraphErrors*> G_TCalibration;
double StepCalibration = 0.002; // s
TGraphErrors *G_TimeShift = new TGraphErrors();


double Intensity_Threshold = 100;
double Pulse_Period = 1.2;
double Neighboor_A_Threshold = 1.5;
double Neighboor_B_Threshold = 2.7;
double Neighboor_C_Threshold = 3.9;

void InitHistograms(string Run)
{
    Info("Initialising Histograms for Run " + Run, 2);

    dir[Run] = fTime->mkdir(Run.c_str());

    H_Time_All_Coinc[Run] = (TH1D *)H_Data[Run]->Clone(("H_Time_All_Coinc_" + Run).c_str());
    H_Time_All_Coinc[Run]->Reset();

    H_Time_All[Run] = (TH1D *)H_Data[Run]->Clone(("H_Time_All_" + Run).c_str());
    H_Time_All[Run]->Reset();

    H_Decay[Run][0] = new TH1D(("H_Decay_" + Run + "_0").c_str(), ("H_Decay_" + Run + "_0").c_str(), Pulse_Period*4*1000, 0, Pulse_Period*4);
    H_Decay[Run][1] = new TH1D(("H_Decay_" + Run + "_A").c_str(), ("H_Decay_" + Run + "_A").c_str(), Pulse_Period*4*1000, 0, Pulse_Period*4);
    H_Decay[Run][2] = new TH1D(("H_Decay_" + Run + "_B").c_str(), ("H_Decay_" + Run + "_B").c_str(), Pulse_Period*4*1000, 0, Pulse_Period*4);
    H_Decay[Run][3] = new TH1D(("H_Decay_" + Run + "_C").c_str(), ("H_Decay_" + Run + "_C").c_str(), Pulse_Period*4*1000, 0, Pulse_Period*4);
    H_Decay[Run][4] = new TH1D(("H_Decay_" + Run + "_D").c_str(), ("H_Decay_" + Run + "_D").c_str(), Pulse_Period*4*1000, 0, Pulse_Period*4);

    H_Decay_Coinc[Run][0] = (TH1D *)H_Decay[Run][0]->Clone(("H_Decay_Coinc_" + Run + "_0").c_str());
    H_Decay_Coinc[Run][1] = (TH1D *)H_Decay[Run][1]->Clone(("H_Decay_Coinc_" + Run + "_A").c_str());
    H_Decay_Coinc[Run][2] = (TH1D *)H_Decay[Run][2]->Clone(("H_Decay_Coinc_" + Run + "_B").c_str());
    H_Decay_Coinc[Run][3] = (TH1D *)H_Decay[Run][3]->Clone(("H_Decay_Coinc_" + Run + "_C").c_str());
    H_Decay_Coinc[Run][4] = (TH1D *)H_Decay[Run][4]->Clone(("H_Decay_Coinc_" + Run + "_D").c_str());

    H_Decay_All[Run][0] = (TH1D *)H_Decay[Run][0]->Clone(("H_Decay_All_" + Run + "_0").c_str());
    H_Decay_All[Run][1] = (TH1D *)H_Decay[Run][1]->Clone(("H_Decay_All_" + Run + "_A").c_str());
    H_Decay_All[Run][2] = (TH1D *)H_Decay[Run][2]->Clone(("H_Decay_All_" + Run + "_B").c_str());
    H_Decay_All[Run][3] = (TH1D *)H_Decay[Run][3]->Clone(("H_Decay_All_" + Run + "_C").c_str());
    H_Decay_All[Run][4] = (TH1D *)H_Decay[Run][4]->Clone(("H_Decay_All_" + Run + "_D").c_str());
   
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
      if (!IsDetectorSiliStrip(i))
        continue;
      
      // if (VERBOSE == 1) Info("Detector : " + detectorName[i], 3);

      dir_det[Run][i] = dir[Run]->mkdir(detectorName[i].c_str());
      
      H_Exp[Run][i] = new TH1D(("H_Exp_" + Run + "_" + detectorName[i]).c_str(), ("H_Exp_" + Run + "_" + detectorName[i]).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
      H_Exp[Run][i]->GetXaxis()->SetTitle("Energy [keV]");
      H_Exp[Run][i]->GetYaxis()->SetTitle("Counts");
      H_Exp[Run][i]->GetXaxis()->CenterTitle();
      H_Exp[Run][i]->GetYaxis()->CenterTitle();

      H_Exp_Coinc[Run][i] = new TH1D(("H_Exp_Coinc_" + Run + "_" + detectorName[i]).c_str(), ("H_Exp_Coinc_" + Run + "_" + detectorName[i]).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
      H_Exp_Coinc[Run][i]->GetXaxis()->SetTitle("Energy [keV]");
      H_Exp_Coinc[Run][i]->GetYaxis()->SetTitle("Counts");
      H_Exp_Coinc[Run][i]->GetXaxis()->CenterTitle();
      H_Exp_Coinc[Run][i]->GetYaxis()->CenterTitle();
    }
}
void LoadHistograms(string Run)
{
  TFile *f = MyTFile((DIR_ROOT_DATA_MATCHED + "Time_Pulse_Calibration_" + to_string(YEAR) + ".root").c_str(), "READ");

  H_Decay_Coinc[Run][1] = (TH1D *)f->Get((Run + "/H_Decay_Coinc_" + Run + "_A").c_str());
  if (H_Decay_Coinc[Run][1] == NULL)
  {
    f->Close();
    Error("H_Decay_Coinc_" + Run + "_A not found in Time_Pulse_Calibration_" + to_string(YEAR) + ".root");
  }
  
  string REFERENCE_RUN_str = REFERENCE_RUN > 100 ? to_string(REFERENCE_RUN) : "0" + to_string(REFERENCE_RUN);
  H_Decay_Coinc[REFERENCE_RUN_str][1] = (TH1D *)f->Get((REFERENCE_RUN_str + "/H_Decay_Coinc_" + REFERENCE_RUN_str + "_A").c_str());
  if (H_Decay_Coinc[REFERENCE_RUN_str][1] == NULL)
  {
    f->Close();
    Error("H_Decay_Coinc_" + REFERENCE_RUN_str + "_A not found in Time_Pulse_Calibration_" + to_string(YEAR) + ".root");
  }
}

void InitCalib()
{
  TFile *CALIBRATED_File = MyTFile((DIR_ROOT_DATA_CALIBRATED + "Calibrated_" + to_string(YEAR)+ ".root").c_str(), "READ");
  for (int i = 0; i < SIGNAL_MAX; i++)
  {
    if (IsDetectorSiliStrip(i))
    {
      Calibration[i] = (TF1 *)CALIBRATED_File->Get(("Calibration_" + detectorName[i]).c_str());

      if (Calibration[i] == NULL)
      {
        Calibration[i] = new TF1(("Calibration_" + detectorName[i]).c_str(), "x", 0, 10000);
        Warning("No calibration found for " + detectorName[i]);
      }
    }
  }
}

string GetNextRun(int next_run)
{
  string stop_time_str = " ";
  string next_run_str = next_run > 100 ? to_string(next_run) : "0" + to_string(next_run);
  string next_run_filename = SearchFiles(DIR_ROOT_DATA_GROUPED, next_run_str);
  TFile *next_run_file = MyTFile((DIR_ROOT_DATA_GROUPED + next_run_filename).c_str(), "READ");
  if (next_run_file == nullptr)
  {
    return stop_time_str;
  }
  TObject *next_start = next_run_file->Get("Start_Time");
  if (!next_start)
  {
    Warning("Impossible to find Start_Time for run " + next_run_str);
  }
  else if (next_start->GetTitle() == " ")
  {
    Warning("No start time found for run " + next_run_str);
  }
  else
  {
    stop_time_str = next_start->GetTitle();
    Success("Using start time of run " + next_run_str , 1);
  }
  next_run_file->Close();

  return stop_time_str;
}

void ReadRunDate(string Run, TFile *f)
{
  if (!f->IsOpen())
  {
    Error("Impossible to open ROOT file");
  }
  TObject *start = f->Get("Start_Time");
  TObject *stop = f->Get("End_Time");
  if (!start)
  {
    Warning("Impossible to find Start_Time for run " + Run);
  }
  else if (!stop)
  {
    Warning("Impossible to find End_Time for run " + Run);
  }
  else
  {
    Info("Found for run " + Run);
    string start_time_str = start->GetTitle();
    string stop_time_str = stop->GetTitle();

    // stop time not written in the .setup file ?
    int next_run = stoi(Run);
    while(stop_time_str == " ")
    {
      Warning("No end time found for run " + to_string(next_run));
      next_run++;
      stop_time_str = GetNextRun(next_run);
      
    }

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
    double diff = stop_time - start_time;

    H_Data[Run] = new TH1D(("H_Data_run_" + Run).c_str(), ("H_Data_run_" + Run).c_str(), diff * 1000, 0, diff);
    H_Data[Run]->GetXaxis()->SetTitle("Time [s]");
    H_Data[Run]->GetYaxis()->SetTitle("Intensity");
    H_Data[Run]->GetXaxis()->CenterTitle();
    H_Data[Run]->GetYaxis()->CenterTitle();
    RunDate[Run] = make_pair(start_time_t, stop_time_t);
  }
  // f->Close();
}

void ReadISOLDE()
{
  fISOLDE = MyTFile((DIR_DATA_ISOLDE+"ISOLDE_pulses.root").c_str(), "RECREATE");
  
  double time;

  for (auto pairs : Map_RunFiles)
  {
    string NUCLEUS = pairs.first;
    vector<string> Runs = pairs.second;
    for (auto Run : Runs)
    {
      string GROUPED_filename = SearchFiles(DIR_ROOT_DATA_GROUPED, Run);
      GROUPED_File = MyTFile((DIR_ROOT_DATA_GROUPED + GROUPED_filename).c_str(), "READ");
      ReadRunDate(Run, GROUPED_File);

      fISOLDE->cd();
      Tree[Run] = new TTree(("Tree_ISOLDE_Proton_" + Run).c_str(), ("Tree_ISOLDE_Proton_" + Run).c_str());
      Tree[Run]->Branch("Time", &time);
    }
  }

  string ISOLDE_filename;
  if (YEAR != 2021)
    ISOLDE_filename = DIR_DATA_ISOLDE + "ISOLDE_data.csv";
  else 
    Error("ISOLDE data not available for 2021");

  ifstream file(ISOLDE_filename.c_str());

  if (!file.is_open())
  {
    Error("Impossible to open ISOLDE file");
  }

  string line;
  int counter = 0;
  int first = 0;
  double milliseconds;
  Start("Reading ISOLDE data");
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

    // tm_data.tm_year += 100; // Adjust for years after 2000

    time_t timefile = mktime(&tm_data);

    // continue if < to 2025/04/17
    if (YEAR == 2025)
    {
      time_t start = RunDate["043"].first;
      if (start > timefile)
        continue;

      time_t stop = RunDate["108"].second;
      if (stop < timefile)
        break;
    }

    for (auto pairs : Map_RunFiles)
    {
      string NUCLEUS = pairs.first;
      vector<string> Runs = pairs.second;
      for (auto Run : Runs)
      {
        if (timefile > RunDate[Run].first && timefile < RunDate[Run].second)
        {
          time = timefile - RunDate[Run].first + milliseconds / 1000;
          H_Data[Run]->SetBinContent(H_Data[Run]->FindBin(time), intensity);
          Tree[Run]->Fill();
          milliseconds = 0;
          break;
        }
      }
    }
  }
}

void WriteISOLDE()
{
  Info("Writing ISOLDE data to file", 2);
  fISOLDE->cd();

  for (auto pairs : Map_RunFiles)
  {
    string NUCLEUS = pairs.first;
    vector<string> Runs = pairs.second;
    for (auto Run : Runs)
    {
      H_Data[Run]->Write();
      Tree[Run]->Write();
    }
  }

  fISOLDE->Close();
}

void LoadISOLDE(string run)
{
  Info("Loading ISOLDE data for run " + run, 2);
  TFile *f = MyTFile(DIR_DATA_ISOLDE+"ISOLDE_pulses.root", "READ");
  H_Data[run] = (TH1D *)f->Get(("H_Data_run_" + run).c_str());
}

int Group_Determining(int i, string run)
{
  if (H_Data[run]->GetBinError(i) < 5)
    return H_Data[run]->GetBinError(i);
  else
    return 0;
}

bool IsCoincidence(double P_Time, vector<vector<pair<Signal, Signal>>> SiPM_Groups)
{

  double NEAREST = 1e9;
  int NEAREST_GROUP_INDEX = -1;
  vector<double> nearest_group_time;
  vector<double> mean_group_time = vector<double>(10, 0);
  if ((SiPM_Groups).size() > 0)
  {
    // cout << "Size: " << (SiPM_Groups).size() << endl;
    // Lopping on subgroups
    for (int i_group = 0; i_group < (SiPM_Groups).size(); i_group++)
    {
      nearest_group_time.push_back(1e9);
      int counter_mean_group_time = 0;
      // if ((SiPM_Groups)[i_group].size() > 9)
      // {
      //     for (int i_pair = 0; i_pair < (SiPM_Groups)[i_group].size(); i_pair++)
      //     {
      //         cout << (SiPM_Groups)[i_group][i_pair].first << endl;
      //     }
      // }
      // Looping on pair of the subgroup
      for (int i_pair = 0; i_pair < (SiPM_Groups)[i_group].size(); i_pair++)
      {

        // Taking valid element of pair with High priority
        if ((SiPM_Groups)[i_group][i_pair].first.isValid && !isnan((SiPM_Groups)[i_group][i_pair].first.Time))
        {
          if (abs(nearest_group_time[i_group]) > abs((SiPM_Groups)[i_group][i_pair].first.Time))
          {
            nearest_group_time[i_group] = (SiPM_Groups)[i_group][i_pair].first.Time;
          }

          mean_group_time[i_group] += (SiPM_Groups)[i_group][i_pair].first.Time;
          counter_mean_group_time++;
        }
        else if ((SiPM_Groups)[i_group][i_pair].second.isValid && !isnan((SiPM_Groups)[i_group][i_pair].second.Time))
        {
          if (abs(nearest_group_time[i_group]) > abs((SiPM_Groups)[i_group][i_pair].second.Time))
          {
            nearest_group_time[i_group] = (SiPM_Groups)[i_group][i_pair].second.Time;
          }
          mean_group_time[i_group] += (SiPM_Groups)[i_group][i_pair].first.Time;
          counter_mean_group_time++;
        }
      }
      mean_group_time[i_group] /= counter_mean_group_time;
      // cout << "Mean Time: " << mean_group_time[i_group] << endl;
    }

    // Saving the nearest group in time and index group
    for (int i_group = 0; i_group < nearest_group_time.size(); i_group++)
    {
      if (abs(NEAREST) > abs(nearest_group_time[i_group]))
      {
        NEAREST_GROUP_INDEX = i_group;
        NEAREST = nearest_group_time[i_group];
      }
    }
  }

  // cout << "Nearest: " << NEAREST << "    " << NEAREST_GROUP_INDEX << endl;

  int Multiplicity;
  if (NEAREST_GROUP_INDEX == -1)
    return false;
  else
    Multiplicity = (SiPM_Groups)[NEAREST_GROUP_INDEX].size();
  

  if (Multiplicity >= 3 && (NEAREST > start_gate && NEAREST < end_gate))
  {
    return true;
  }

  return false;
}

pair<double, double> ComputingChi2(TH1D *H_Ref, TH1D *H, string Run)
{
  double Best_Chi2 = 1e10;
  double Best_Tau = -1;

  double start = -0.6;
  double end = 0.6;

  G_TCalibration[Run] = new TGraphErrors();

  //scanning tau
  for (double tau = start; tau <= end; tau+= StepCalibration)
  {
    ProgressCounter( (tau-start)/StepCalibration, (int)((end - start)/StepCalibration), "Computing Chi2 for " + Run, 1, 2);

    TH1D *H_tau = (TH1D *)H->Clone(("H_tau_" + to_string(tau)).c_str());
    H_tau->Reset(); 

    int bin_tau = tau/abs(tau) * H->FindBin(abs(tau));

    // shifting the histogram by tau
    for (int i = 1; i <= H->GetNbinsX(); i++)
    {
      double bin_center = H->GetBinCenter(i);
      double bin_content = H->GetBinContent(i);
      double bin_error = H->GetBinError(i);

      if (i + bin_tau < 1)
        continue;
      H_tau->SetBinContent(i + bin_tau, bin_content);
      H_tau->SetBinError(i + bin_tau, bin_error);
    }

    gErrorIgnoreLevel = kError;
    double min = tau < 0 ? 0 : tau;
    H_tau->GetXaxis()->SetRangeUser(min, 1.2+min);
    H_Ref->GetXaxis()->SetRangeUser(min, 1.2+min);
    H_tau->Scale(H_Ref->Integral() / H_tau->Integral());
    double Chi2 = H_Ref->Chi2Test(H_tau, "UW CHI2/NDF");
    gErrorIgnoreLevel = kInfo;
    H_Ref->GetXaxis()->SetRangeUser(-1111, -1111);

    G_TCalibration[Run]->AddPoint(-tau, Chi2);

    // save best value
    if (Chi2 < Best_Chi2)
    {
      Best_Chi2 = Chi2;
      Best_Tau = -tau;
    }
  } 
  return make_pair(Best_Tau, Best_Chi2);
}

void WriteHistogram(string Run, string nucleus)
{
  Info("Writing histograms for " + Run, 2);
  fTime->cd();
  dir[Run]->cd();
  H_Time_All_Coinc[Run]->Write();
  H_Time_All[Run]->Write();
  H_Data[Run]->Write();
  H_Decay_Coinc[Run][1]->Write();
  H_Decay_All[Run][1]->Write();

  for (int i = 0; i < SIGNAL_MAX; i++)
  {
    if (!IsDetectorSiliStrip(i))
      continue;

    // up histograms comparison and residuals between then
    dir_det[Run][i]->cd();
    TCanvas *c = new TCanvas(("Channel_" + Run + "_" + detectorName[i]).c_str(), ("ISOLDE_DETECTOR_" + Run + "_" + detectorName[i]).c_str(), 800, 600);
    H_Exp[Run][i]->GetXaxis()->SetRangeUser(WindowsMap[nucleus][IAS[nucleus]][i].first, WindowsMap[nucleus][IAS[nucleus]][i].second);
    H_Exp[Run][i]->Draw("HIST");
    H_Exp_Coinc[Run][i]->SetLineColor(kRed);
    H_Exp_Coinc[Run][i]->Draw("HIST SAME");
    c->Write();
  }

  dir[Run]->cd();

  //// A

  // up histograms comparison and residuals between then
  TCanvas *cA_Decay = new TCanvas(("cA_Decay_" + Run).c_str(), ("cA_Decay_" + Run).c_str(), 800, 600);
  TPad *pad1A = new TPad("pad1A", "pad1A", 0, 0.3, 1, 1);
  pad1A->Draw();
  pad1A->cd();
  H_Decay[Run][1]->Draw("HIST");
  H_Decay_Coinc[Run][1]->SetLineColor(kRed);
  H_Decay_Coinc[Run][1]->Draw("HIST SAME");
  cA_Decay->cd();
  TPad *pad2A = new TPad("pad2A", "pad2A", 0, 0, 1, 0.3);
  pad2A->Draw();
  pad2A->cd();
  TH1D *hResidualsA = (TH1D *)H_Decay[Run][1]->Clone();
  hResidualsA->Add(H_Decay_Coinc[Run][1], -1);
  hResidualsA->Divide(H_Decay[Run][1]);
  hResidualsA->Draw();
  cA_Decay->Write();

  ////// B
  TCanvas *cB_Decay = new TCanvas(("cB_Decay_" + Run).c_str(), ("cB_Decay_" + Run).c_str(), 800, 600);
  TPad *pad1B = new TPad("pad1B", "pad1B", 0, 0.3, 1, 1);
  pad1B->Draw();
  pad1B->cd();
  H_Decay[Run][2]->Draw("HIST");
  H_Decay_Coinc[Run][2]->SetLineColor(kRed);
  H_Decay_Coinc[Run][2]->Draw("HIST SAME");
  cB_Decay->cd();
  TPad *pad2B = new TPad("pad2B", "pad2B", 0, 0, 1, 0.3);
  pad2B->Draw();
  pad2B->cd();
  TH1D *hResidualsB = (TH1D *)H_Decay[Run][2]->Clone();
  hResidualsB->Add(H_Decay_Coinc[Run][2], -1);
  hResidualsB->Divide(H_Decay[Run][2]);
  hResidualsB->Draw();
  cB_Decay->Write();

  ////// C
  TCanvas *cC_Decay = new TCanvas(("cC_Decay_" + Run).c_str(), ("cC_Decay_" + Run).c_str(), 800, 600);
  TPad *pad1C = new TPad("pad1C", "pad1C", 0, 0.3, 1, 1);
  pad1C->Draw();
  pad1C->cd();
  H_Decay[Run][3]->Draw("HIST");
  H_Decay_Coinc[Run][3]->SetLineColor(kRed);
  H_Decay_Coinc[Run][3]->Draw("HIST SAME");
  cC_Decay->cd();
  TPad *pad2C = new TPad("pad2C", "pad2C", 0, 0, 1, 0.3);
  pad2C->Draw();
  pad2C->cd();
  TH1D *hResidualsC = (TH1D *)H_Decay[Run][3]->Clone();
  hResidualsC->Add(H_Decay_Coinc[Run][3], -1);
  hResidualsC->Divide(H_Decay[Run][3]);
  hResidualsC->Draw();  
  cC_Decay->Write();

  ////// D
  TCanvas *cD_Decay = new TCanvas(("cD_Decay_" + Run).c_str(), ("cD_Decay_" + Run).c_str(), 800, 600);
  TPad *pad1D = new TPad("pad1D", "pad1D", 0, 0.3, 1, 1);
  pad1D->Draw();
  pad1D->cd();
  H_Decay[Run][4]->Draw("HIST");
  H_Decay_Coinc[Run][4]->SetLineColor(kRed);
  H_Decay_Coinc[Run][4]->Draw("HIST SAME");
  cD_Decay->cd();
  TPad *pad2D = new TPad("pad2D", "pad2D", 0, 0, 1, 0.3);
  pad2D->Draw();
  pad2D->cd();
  TH1D *hResidualsD = (TH1D *)H_Decay[Run][4]->Clone();
  hResidualsD->Add(H_Decay_Coinc[Run][4], -1);
  hResidualsD->Divide(H_Decay[Run][4]);
  hResidualsD->Draw();
  cD_Decay->Write();

  delete H_Time_All_Coinc[Run];
  delete H_Time_All[Run];
  delete H_Data[Run];
  for (int i = 0; i < SIGNAL_MAX; i++)
  {
    if (!IsDetectorSiliStrip(i))
      continue;

    delete H_Exp[Run][i];
    delete H_Exp_Coinc[Run][i];
  }
  delete H_Decay[Run][0];
  delete H_Decay[Run][1];
  delete H_Decay[Run][2];
  delete H_Decay[Run][3];
  delete H_Decay[Run][4];
  
}

#endif