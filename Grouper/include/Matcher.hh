#ifndef MATCHER_HH
#define MATCHER_HH

#include "../../../lib/SignalDict/Signal.h"
#include "Detectors.hh"

default_random_engine generator;

/// FILE ///
string GROUPED_filename;
string GROUPED_basefilename;
TFile *GROUPED_File;
string REFERENCE_filename;
TFile *REFERENCE_File;
TFile *MATCHED_File;

map<string, vector<string>> Map_RunFiles;   

int Run;
vector<string> Runs;
map<int, int> REFERENCE_Run;

TTree *GROUPED_Tree_Detectors[SIGNAL_MAX];
TTreeReader *Reader;
TTreeReaderArray<Signal> *Silicon;
double Channel;

//HISTOGRAMS
TCanvas *C_Det[SIGNAL_MAX];
TCanvas *C_Det_corr[SIGNAL_MAX];
TCanvas *C_Det_vs[SIGNAL_MAX];
TH1D *H_Run[999][SIGNAL_MAX];
TH1D *H_Run_corr[999][SIGNAL_MAX];
TH1D *H_Run_Ref[SIGNAL_MAX];
TF1 *fpol1[999][SIGNAL_MAX];
pair<double, double> ManualCalib[100][SIGNAL_MAX];
TF1 *linearfit[SIGNAL_MAX];
TGraph *G_ManualCalib[SIGNAL_MAX];
TGraphErrors *G_Mean_Run[99][SIGNAL_MAX];
TGraphErrors *G_Mean_Run_corr[99][SIGNAL_MAX];
TGraphErrors *G_resPar0[SIGNAL_MAX];
TGraphErrors *G_resPar1[SIGNAL_MAX];
TLine *line[99][SIGNAL_MAX];
TH1D* H_Run_BEFORE[SIGNAL_MAX];
TH1D* H_Run_AFTER[SIGNAL_MAX];
int counter_graph[SIGNAL_MAX] = {0};
pair<double, double> WindowsMap[100][SIGNAL_MAX];

int current_detector = 0;
double CHI2;
double parameter[2];

vector<int> Peaks = {5, 8, 14, 23};

void InitHistograms(string Run_string)
{
    int Run = atoi(Run_string.c_str());
    for (int i = 0; i < detectorNum; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            H_Run[Run][i] = new TH1D(("H_Run_" + Run_string + "_" + detectorName[i]).c_str(), ("H_Run_" + Run_string + "_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
            H_Run[Run][i]->GetXaxis()->SetTitle("Channel");
            H_Run[Run][i]->GetYaxis()->SetTitle("Counts");
            H_Run[Run][i]->GetXaxis()->CenterTitle();
            H_Run[Run][i]->GetYaxis()->CenterTitle();

            H_Run_corr[Run][i] = new TH1D(("H_Run_corr_" + Run_string + "_" + detectorName[i]).c_str(), ("H_Run_corr_" + Run_string + "_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
            H_Run_corr[Run][i]->GetXaxis()->SetTitle("Channel");
            H_Run_corr[Run][i]->GetYaxis()->SetTitle("Counts");
            H_Run_corr[Run][i]->GetXaxis()->CenterTitle();
            H_Run_corr[Run][i]->GetYaxis()->CenterTitle();
        }
    }
}

void InitRuns()
{
  ifstream file(("./Config_Files/" + to_string(YEAR) + "/Runs_" + to_string(YEAR) + ".txt").c_str());
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

    REFERENCE_Run[2021] = 36;
    REFERENCE_Run[2024] = 77;
    REFERENCE_Run[2025] = 69;
}

double FunctionToMinimize(const double *par)
{
    // cout << setprecision(5) << "PAR : " << par[0] << " " << par[1] << endl;

    Reader->Restart();
    if (par[0] == 0)
    {
        H_Run[Run][current_detector]->Reset();
        TTreeReaderValue<double> ChannelDet(*Reader, "Channel");
        while (Reader->Next())
        {
            H_Run[Run][current_detector]->Fill(par[0] + *ChannelDet * par[1]);
        }
    }
    else
    {
        H_Run_corr[Run][current_detector]->Reset();
        TTreeReaderValue<double> ChannelDet(*Reader, "Channel");
        while (Reader->Next())
        {
            H_Run_corr[Run][current_detector]->Fill(par[0] + *ChannelDet * par[1]);
        }
    }

    double chi2 = 0;
    //chi2 = H_Run[current_detector]->Chi2Test(H_Run_Ref[current_detector], "UW CHI2/NDF");

    if (CHI2 > chi2)
    {
        CHI2 = chi2;
        parameter[0] = par[0];
        parameter[1] = par[1];
    }

    // cout << "Chi2 : " << chi2 << endl;
    // cout << endl;
    
    return chi2;
}

void CHI2Minimization(int i)
{
    current_detector = i;

    Minimizer *minimizer = Factory::CreateMinimizer("Minuit2", "Migrad");
    ROOT::Math::Functor functor(&FunctionToMinimize, 2);
    minimizer->SetFunction(functor);
    minimizer->SetFixedVariable(0, "Offset", 0);
    minimizer->SetFixedVariable(1, "Proportionnality", 1);
    minimizer->SetPrecision(0.00000001);
    minimizer->SetTolerance(0.00000001);
    minimizer->SetMaxFunctionCalls(1000000);
    minimizer->SetMaxIterations(1000000);

    // minimizer->Minimize();
    // const double *par = minimizer->X();
    
    const double par[2] = {0, 1.0};
    // cout << par[0] << " " << par[1] << endl;
    FunctionToMinimize(par);

    C_Det[i]->cd();

    for (int peak : Peaks)
    {   
        H_Run[Run][i]->GetXaxis()->SetRangeUser(WindowsMap[peak][i].first, WindowsMap[peak][i].second);
        H_Run_Ref[i]->GetXaxis()->SetRangeUser(WindowsMap[peak][i].first, WindowsMap[peak][i].second);

        if (H_Run[Run][i]->GetMean() == 0)
            continue;
        G_Mean_Run[peak][i]->AddPoint(Run, H_Run[Run][i]->GetMean());
        G_Mean_Run[peak][i]->SetPointError(G_Mean_Run[peak][i]->GetN()-1, 0, H_Run[Run][i]->GetMeanError());
    }

    H_Run[Run][i]->GetXaxis()->SetRangeUser(WindowsMap[14][i].first, WindowsMap[14][i].second);
    H_Run_Ref[i]->GetXaxis()->SetRangeUser(WindowsMap[14][i].first, WindowsMap[14][i].second);
    // H_Run[Run][i]->Scale(H_Run_Ref[i]->Integral() / H_Run[Run][i]->Integral());
    H_Run[Run][i]->SetLineColor(Run);
    H_Run[Run][i]->Draw("HIST SAME");
    
    
    
    // //// correction
    
    TCanvas *cfit_corr = new TCanvas(("cfit_corr_Run_"+to_string(Run)+"_"+to_string(i)).c_str(), ("cfit_corr_Run_"+to_string(Run)+"_"+to_string(i)).c_str(), 1920, 1080);

    TGraphErrors *G = new TGraphErrors();
    fpol1[Run][i] = new TF1(("poll1"+to_string(Run)+to_string(i)).c_str(), "[0] + [1]*x", 0, 100000);

    // correction
    for (int peak : Peaks)
    {
        
        if (WindowsMap[peak][i].first == -1 || !WindowsMap[peak][i].first)
            continue;

        H_Run[Run][i]->GetXaxis()->SetRangeUser(WindowsMap[peak][i].first, WindowsMap[peak][i].second);
        H_Run_Ref[i]->GetXaxis()->SetRangeUser(WindowsMap[peak][i].first, WindowsMap[peak][i].second);

        if (H_Run[Run][i]->GetMean() < 1000 || H_Run_Ref[i]->GetMean() < 1000)
            continue;

        // Info("Peak : " + to_string(peak), 2);
        G->AddPoint(H_Run[Run][i]->GetMean(), H_Run_Ref[i]->GetMean());
        G->SetPointError(G->GetN()-1, H_Run[Run][i]->GetMeanError(), H_Run_Ref[i]->GetMeanError());
    }
    fpol1[Run][i]->SetParLimits(0, -1000, 1000);
    fpol1[Run][i]->SetParameter(0, 0);
    fpol1[Run][i]->SetParLimits(1, 0.98, 1.02);
    fpol1[Run][i]->SetParameter(1, 1);
    if (G->GetN() <= 1)
    {
        Warning("No points to fit for Run " + to_string(Run) + " Detector " + detectorName[i]);
        return;
    }

    TFitResultPtr s = G->Fit(fpol1[Run][i], "QSE");
    if (s->Status() != 0)
    {
        Warning("Fit failed for Run " + to_string(Run) + " Detector " + detectorName[i]);
    }
    cfit_corr->Divide(1, 2);
    //plotting
    cfit_corr->cd(1);
    G->Draw("AP");
    fpol1[Run][i]->Draw("SAME");
    // residus
    TGraphErrors *G_residus = (TGraphErrors*)G->Clone();
    for (int j = 0; j < G_residus->GetN(); j++)
    {
        double x, y;
        G_residus->GetPoint(j, x, y);
        G_residus->SetPoint(j, x, y - fpol1[Run][i]->Eval(x));
    }
    cfit_corr->cd(2);
    G_residus->Draw("AP");
    TLine *line = new TLine(G_residus->GetXaxis()->GetXmin(), 0, G_residus->GetXaxis()->GetXmax(), 0);
    line->SetLineColor(kRed);
    line->Draw("SAME");
    // cfit_corr->Write();

    C_Det_corr[i]->cd();



    const double par_corr[2] = {fpol1[Run][i]->GetParameter(0), fpol1[Run][i]->GetParameter(1)};
    FunctionToMinimize(par_corr);


    for (int peak : Peaks)
    {   
        H_Run_corr[Run][i]->GetXaxis()->SetRangeUser(WindowsMap[peak][i].first, WindowsMap[peak][i].second);
        H_Run_Ref[i]->GetXaxis()->SetRangeUser(WindowsMap[peak][i].first, WindowsMap[peak][i].second);

        G_Mean_Run_corr[peak][i]->AddPoint(Run, H_Run_corr[Run][i]->GetMean());
        G_Mean_Run_corr[peak][i]->SetPointError(G_Mean_Run_corr[peak][i]->GetN()-1, 0, H_Run_corr[Run][i]->GetMeanError());
    }

    H_Run_corr[Run][i]->GetXaxis()->SetRangeUser(WindowsMap[14][i].first, WindowsMap[14][i].second);
    H_Run_Ref[i]->GetXaxis()->SetRangeUser(WindowsMap[14][i].first, WindowsMap[14][i].second);
    // H_Run_corr[Run][i]->Scale(H_Run_Ref[i]->Integral() / H_Run_corr[Run][i]->Integral());
    H_Run_corr[Run][i]->SetLineColor(Run);
    H_Run_corr[Run][i]->Draw("HIST SAME");
}

void InitManualCalibration()
{
    ifstream file(("Config_Files/" + to_string(YEAR) + "/Manual_Calibration_" + to_string(YEAR) + ".txt").c_str());
    if (!file.is_open())
    {
        Error("Impossible to open Manual_Calibration_" + to_string(YEAR) + ".txt");
    }

    string line;
    double channel6 = -1, channel14 = -1, channel29 = -1;
    double energy6 = -1, energy14 = -1, energy29 = -1;
    int counter = 0;
    string detname;
    while (getline(file, line))
    {
        counter++;
        channel6 = -1, channel14 = -1, channel29 = -1;
        energy6 = -1, energy14 = -1, energy29 = -1;

        stringstream ss(line);

        if (YEAR == 2024)
        {
            ss >> detname >> channel6 >> channel14 >> channel29 >> energy6 >> energy14 >> energy29;
            for (int i = 0; i < detectorNum; i++)
            {
                if (IsDetectorSiliStrip(i))
                {
                    if (detname == detectorName[i])
                    {
                        ManualCalib[6][i] = make_pair(channel6 * 1000, energy6);
                        ManualCalib[14][i] = make_pair(channel14 * 1000, energy14);
                        ManualCalib[29][i] = make_pair(channel29 * 1000, energy29);

                        linearfit[i] = new TF1("linearfit", "[0] + [1]*x", 0, 10000);
                        G_ManualCalib[i] = new TGraph();
                        G_ManualCalib[i]->SetPoint(0, ManualCalib[6][i].second, ManualCalib[6][i].first);
                        G_ManualCalib[i]->SetPoint(1, ManualCalib[14][i].second, ManualCalib[14][i].first);
                        G_ManualCalib[i]->SetPoint(2, ManualCalib[29][i].second, ManualCalib[29][i].first);
                        G_ManualCalib[i]->Fit(linearfit[i], "Q");

                        break;
                    }
                }
            }
        }
        else if (YEAR == 2025 || YEAR == 2021)
        {
            ss >> detname >> channel6 >> channel14 >> energy6 >> energy14;
            for (int i = 0; i < detectorNum; i++)
            {
                if (IsDetectorSiliStrip(i))
                {
                    if (detname == detectorName[i])
                    {
                        ManualCalib[6][i] = make_pair(channel6 * 1000, energy6);
                        ManualCalib[14][i] = make_pair(channel14 * 1000, energy14);

                        linearfit[i] = new TF1("linearfit", "[0] + [1]*x", 0, 100000);
                        G_ManualCalib[i] = new TGraph();
                        G_ManualCalib[i]->SetPoint(0, ManualCalib[6][i].second, ManualCalib[6][i].first);
                        G_ManualCalib[i]->SetPoint(1, ManualCalib[14][i].second, ManualCalib[14][i].first);
                        G_ManualCalib[i]->Fit(linearfit[i], "Q");

                        

                        break;
                    }
                }
            }
        }
    }

    file.close();

    
}


void InitWindows()
{
    string direction[2] = {"Up", "Down"};
    for (auto dir : direction)
    {
        for (int strip = 1; strip <= 5; strip++)
        {
            ifstream file("Config_Files/Detector_Window/" + dir + "_" + to_string(strip) + ".txt");
            if (!file.is_open())
            {
                Error("Impossible to open " + dir + "_" + to_string(strip) + ".txt");
            }

            string line;
            double energy_low;
            double energy_high;
            int number;
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

                if (nuclei == "32Ar")
                {
                    for (int i : Dir2Det(dir, strip))
                    {
                        WindowsMap[number][i] = make_pair(linearfit[i]->Eval(energy_low), linearfit[i]->Eval(energy_high));
                        // cout << "Nuclei : " << nuclei << " Number : " << number << " Detector : " << detectorName[i] << " Energy Low : " << energy_low << " Energy High : " << energy_high << endl;
                    }
                }
            }
        }
    }
}



#endif