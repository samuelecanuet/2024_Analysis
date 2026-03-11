#ifndef SiPM_MATCHER_HH
#define SiPM_MATCHER_HH

#include "../../../lib/SignalDict/Signal.h"
#include "Detectors.hh"

/// FILE ///
map<string, TFile *> GROUPED_File;
map<string, TFile *> SIMULATED_File;
TFile *MATCHED_File;
string MATCHED_Filename;
TFile *f_tree;

TTree *Tree;
TTreeReader *Reader;
TTreeReaderArray<Signal> *Silicon;
TTreeReaderValue<vector<vector<pair<Signal, Signal>>>> *SiPM_Groups;
TTreeReaderValue<Signal> *HRS;  

TTree *Tree_SIMULATED;
TTreeReader *Reader_SIMULATED;
TTreeReaderValue<int> *Silicon_code;
TTreeReaderValue<double> *Silicon_energy;
TTreeReaderValue<double> *SiPM_energy;

TTree *Tree_MATCHED;
vector<Signal> SiPM;

int Entry_MAX = 10e8;
bool FULL = false;

// TREE PEAKS //
map<string, TTree *[50][SIGNAL_MAX]> Tree_Peaks;
Signal SiPMi;

/// HISTOGRAMS ///
// ONOFF
map<int, map<int, TH1D *>> H_SiPM_ONOFF_RUN;    
map<int, map<int, TTree *>>Tree_SiPM_ONOFF_RUN;
map<int, pair<TH1D *, TH1D *>> H_SiPM_ONOFF;
map<int, map<int, TF1*>> SiPM_ONOFF_Coef_RUN;
TH1D *H_SiPM_ONOFF_Corrected[SIGNAL_MAX];
TTree *Tree_SiPM_ONOFF[SIGNAL_MAX];
double Channel;
// single
map<string, TH1D *[SIGNAL_MAX]> H_SiPM_High;
map<string, TH1D *[SIGNAL_MAX]> H_SiPM_Low;
// matching low high
map<string, TH2D *[SIGNAL_MAX]> H_SiPM_HighLow;
map<string, TGraph *[SIGNAL_MAX]> G_SiPM_HighLow;
map<string, TH2D *[SIGNAL_MAX]> H_SiPM_Gain;
map<string, TGraph *[SIGNAL_MAX]> G_SiPM_Gain;
// Results
TGraphErrors *gHighLow_RUNS_a[SIGNAL_MAX];
TGraphErrors *gSiPM_High_RUNS_a[SIGNAL_MAX];
TGraphErrors *gSiPM_Low_RUNS_a[SIGNAL_MAX];
TGraphErrors *gHighLow_RUNS_b[SIGNAL_MAX];   
TGraphErrors *gSiPM_High_RUNS_b[SIGNAL_MAX];
TGraphErrors *gSiPM_Low_RUNS_b[SIGNAL_MAX];

/// DIRECTORY ///
map<string, TDirectory *> dir;
map<string, TDirectory *[SIGNAL_MAX]> dir_detector;
map<string, TDirectory *[SIGNAL_MAX]> dir_SiPM;
map<string, TDirectory *> dir_Nucleus;

// FUCNTIONS //
TF1 *F_SiliconCalibration[SIGNAL_MAX];
TF1 *f_linear = new TF1("f_linear", "[0]*x + [1]", 0, 10e6);
TF1 *f_linear_fixed = new TF1("f_linear_fixed", "[0]*x+[1]", 0, 10e6);
TF1 *f_erf = new TF1("f_erf", "x < [0] ? [1]*x + [2] : [3]*(erf((x-[4])/[5]))", 0, 10000e3);
TF1 *TF1_DoubleGaussian;
TF1 *TF1_SimpleGaussian;
pair<double, double> Range_SiPM_LowHigh = make_pair(50e3, 150e3);
int current_detector;
double SiPM_Window[10] = {0, 1500e3, 1500e3, 1500e3};
// FINAL FUNCTIONS
TF1 *SiPM_ONOFF[SIGNAL_MAX];   // fitted on low used for both
TF1 *SiPM_HighLow[SIGNAL_MAX]; // fitted on both used for low
TF1 *SiPM_Gain[SIGNAL_MAX];    // fitted on high used for both
// DATA //
string Nuclei[4] = {"90Sr", "207Bi", "32Ar", "33Ar"};
string TYPE;
string RUN;
string NUCLEI;
double SiliconCalibrationParameter[SIGNAL_MAX][3];
map<string, TF1 *[SIGNAL_MAX]> MatchingLowHigh;
map<string, TF1 *[SIGNAL_MAX]> MatchingSiPM;


void init()
{
    // Map_RunFiles["32Ar"] =
    //     {"057", "058", "059",
    //      "061", "062", "064", "065", "066", "067", "068", "069",
    //      "070", "071", "072", "074", "075", "076", 
    //     "077",

    //      "112", "113", "114", "115", "116", "118"};

    // Map_RunFiles["33Ar"] = {"078"};

    if (YEAR == 2024)
    {
        // Map_RunFiles["207Bi"] = {"137"};

        // Map_RunFiles["90Sr"] = {"133"};
        Map_RunFiles["32Ar"] = {"064", "065", "066", "067", 
                                "068", "069", "070", "071", "072",
                                 "074", "075", "076", "077",
                                 "112", "113", "114", 
                                 "115", 
                                 "116", "118"
                                };
        // Map_RunFiles["33Ar"] = {"078"};
    }
    else if (YEAR == 2025)
    {

        Map_RunFiles["32Ar"] = {"043", "044", 
        "045", "046", "047", "048", "049", "050", "051", "052", "053",
                                 "055", "056", "057", "058", "059", "060", "061", "062", "063", "064",
                                 "065", "069", "070", "071", "072", "073", "074", "075", "076",
                                 "083", "084", "085", "086", "087", "088", "089", "090", "091", "092",
                                 "093", "094", "095", "096", "097", "098", "099", "100", "101", "102",
                                 "103", "104", "105", "106", "107", 
                                 "108"};
        // Map_RunFiles["207Bi"] = {"116"};
        // Map_RunFiles["90Sr"] = {"122"};
    }
    else 
    {
        Error("Year not recognized or selected 2021");
    }
}

TF1 *InvertingLinear(TF1 *f)
{
    if (f->GetNpar() == 2)
    {
        TF1 *f_inv = new TF1("f_inv", "[0]*x + [1]", 0, 10000e3);
        double a = f->GetParameter(0);
        double b = f->GetParameter(1);
        f_inv->SetParameter(0, 1. / a);
        f_inv->SetParameter(1, -b / a);
        return f_inv;
    }

    else if (f->GetNpar() == 1)
    {
        TF1 *f_inv = new TF1("f_inv", "[0]*x", 0, 10000e3);
        double a = f->GetParameter(0);
        f_inv->SetParameter(0, 1. / a);
        return f_inv;
    }

    else
    {
        Error("Function not linear, cannot invert");
        return nullptr;
    }
}

double DELTA_PileUp_Range = 200e3;
pair<TF1 *, TF1*> PileUp_Function[SIGNAL_MAX];
void InitPileUp()
{
    Info("Init PileUp");

    if (YEAR == 2024)
    {
        // BEFORE GAIN CHANGE //
        TFile *f_pileup1 = MyTFile((DIR_ROOT_DATA_MATCHED + "PileUp_run_077.root").c_str(), "READ");

        if (f_pileup1 != nullptr)
        {
            for (int det = 0; det < SIGNAL_MAX  ; det++)
            {
                if (IsDetectorBeta(det))
                {
                    PileUp_Function[det].first = (TF1 *)f_pileup1->Get((detectorName[det]+"/PileUp_Correction_"+detectorName[det]).c_str());
                    
                    if (PileUp_Function[det].first == NULL)
                    {
                        Error("PileUp Function before gain changenot found for detector: " + detectorName[det]);
                    }
                }
            }
        }

        // AFTER GAIN CHANGE //
        TFile *f_pileup2 = MyTFile((DIR_ROOT_DATA_MATCHED + "PileUp_run_114.root").c_str(), "READ");

        if (f_pileup2 != nullptr)
        {
            for (int det = 0; det < SIGNAL_MAX  ; det++)
            {
                if (IsDetectorBeta(det))
                {
                    PileUp_Function[det].second = (TF1 *)f_pileup2->Get((detectorName[det]+"/PileUp_Correction_"+detectorName[det]).c_str());
                    
                    if (PileUp_Function[det].second == NULL)
                    {
                        Error("PileUp Function after gain changenot found for detector: " + detectorName[det]); 
                    }
                }
            }
        }
    }
}
bool PileUp(Signal SiPM)
{
    double coef = 1;
    if (SiPM.Label > 110)
        coef = 10;
    if (stoi(RUN) < 112)
    {
        if (abs(PileUp_Function[SiPM.Label].second->Eval(SiPM.Channel) - SiPM.Pileup) < DELTA_PileUp_Range/coef)
            return false;
        return true;
    }
    else
    {
        if (abs(PileUp_Function[SiPM.Label].first->Eval(SiPM.Channel) - SiPM.Pileup) < DELTA_PileUp_Range/coef)
            return false;
        return true;
    }
}

int WriteONOFFValues(bool first, bool dummy = false)
{   
    Info("Write matching ONOFF values to file");
    string option;
    if (first) option = "RECREATE";
    else option = "UPDATE";
    TFile *f = MyTFile((DIR_ROOT_DATA_MATCHED + "SiPM_Matching_values.root").c_str(), option.c_str());
    f->cd();
    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if(IsDetectorBeta(det))
        {
            for (auto &pairr : Map_RunFiles)
            {
                NUCLEI = pairr.first;
                for (string run : pairr.second)
                {
                    if (dummy || SiPM_ONOFF_Coef_RUN[det][stoi(run)] == nullptr)
                    {
                        SiPM_ONOFF_Coef_RUN[det][stoi(run)] = new TF1(("SiPM_ONOFF_Coef_RUN_" + detectorName[det] + "_run_" + run).c_str(), "[0]*x + [1]", 0, 100e6);
                        SiPM_ONOFF_Coef_RUN[det][stoi(run)]->SetParameter(0, 1.);
                        SiPM_ONOFF_Coef_RUN[det][stoi(run)]->SetParameter(1, 0.);
                    }
                    SiPM_ONOFF_Coef_RUN[det][stoi(run)]->SetName(("SiPM_ONOFF_Coef_RUN_" + detectorName[det] + "_run_" + run).c_str());
                    SiPM_ONOFF_Coef_RUN[det][stoi(run)]->Write();
                }
            }
        }
    }

    f->Close();
    MATCHED_File->cd();
    return 0;
}

int LoadONOFFValues()
{
    Warning("Read matching values from file");
    TFile *f = MyTFile((DIR_ROOT_DATA_MATCHED + "SiPM_Matching_values.root").c_str(), "READ");
    if (f == NULL)
    {
        return 1;
    }

    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorBeta(det))
        {

            for (auto &pairr : Map_RunFiles)
            {
                NUCLEI = pairr.first;
                for (string run : pairr.second)
                {
                    SiPM_ONOFF_Coef_RUN[det][stoi(run)] = ((TF1 *)f->Get(("SiPM_ONOFF_Coef_RUN_" + detectorName[det] + "_run_" + run).c_str()));
                }
            }
        }

        // f->Close();
        MATCHED_File->cd();
        return 0;
    }

    return 0;
}

void WriteValues(bool first)
{   
    string option;
    if (first) option = "RECREATE";
    else option = "UPDATE";

    Info("Write matching values to file");
    TFile *f = MyTFile((DIR_ROOT_DATA_MATCHED + "SiPM_Matching_values.root").c_str(), option.c_str());
    f->cd();
    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if(IsDetectorBeta(det))
        {
            MatchingSiPM[RUN][det]->Write(("MatchingSiPM_" + RUN + "_" + detectorName[det]).c_str());
            if (IsDetectorBetaHigh(det))
            {
                MatchingLowHigh[RUN][det]->Write(("MatchingLowHigh_" + RUN + "_" + detectorName[det]).c_str());
            }
        }
    }
    f->Close();
    MATCHED_File->cd();
}

int LoadValues(bool inverted_for_applying = true)
{
    Warning("Read matching values from file");
    // TFile *f = MyTFile((DIR_ROOT_DATA_MATCHED + "SiPM_Matching_values.root").c_str(), "READ");
    // if (f == NULL)
    // {
    //     return 1;
    // }
    // for (int det = 0; det < SIGNAL_MAX; det++)
    // {
    //     if (IsDetectorBeta(det))
    //     {
    //         if (inverted_for_applying)
    //             MatchingSiPM[RUN][det] = InvertingLinear((TF1 *)f->Get(("MatchingSiPM_" + RUN + "_" + detectorName[det]).c_str()));
    //         else
    //             MatchingSiPM[RUN][det] = (TF1 *)f->Get(("MatchingSiPM_" + RUN + "_" + detectorName[det]).c_str());
    //         if (IsDetectorBetaHigh(det))
    //         {
    //             MatchingLowHigh[RUN][det] = (TF1 *)f->Get(("MatchingLowHigh_" + RUN + "_" + detectorName[det]).c_str());
    //         }
    //     }
    // }

    // // f->Close();
    // MATCHED_File->cd();

    TFile *f = MyTFile((DIR_ROOT_DATA_MATCHED + "SiPM_Matching_Functions.root").c_str(), "READ");
    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorBeta(det))
        {
            SiPM_ONOFF[det] = ((TF1 *)f->Get(("SiPM_ONOFF_" + detectorName[det]).c_str()));
            if (SiPM_ONOFF[det] == nullptr)
            {
                Error("SiPM ONOFF Function not found for detector: " + detectorName[det]);
            }
        }
        if (IsDetectorBetaLow(det))
        {
            SiPM_HighLow[det] = ((TF1 *)f->Get(("SiPM_HighLow_" + to_string(GetDetectorChannel(det))).c_str()));
            if (SiPM_HighLow[det] == nullptr)
            {
                Error("SiPM HighLow Function not found for detector: " + detectorName[det]);
            }
        }
        if (IsDetectorBeta(det))
        {
            if (inverted_for_applying)
                SiPM_Gain[det] = InvertingLinear((TF1 *)f->Get(("SiPM_Gain_" + detectorName[det]).c_str()));
            else
                SiPM_Gain[det] = ((TF1 *)f->Get(("SiPM_Gain_" + detectorName[det]).c_str()));

            if (SiPM_Gain[det] == nullptr)
            {
                Error("SiPM Gain Function not found for detector: " + detectorName[det]);
            }
        }
    }

    MATCHED_File->cd();

    return 0;
}

void WriteTree()
{
    Info("Write Tree start");
    f_tree->cd();
    RUN = "32Ar";
    for (int peak = 1; peak <= 50; peak++)
    {
        if (WindowsMap[RUN][peak][11].first == -1 || !WindowsMap[RUN][peak][11].first)
            continue;
        for (int det = 1; det <= 9; det++)
        {
            Tree_Peaks[RUN][peak][det]->Write();
        }
    }
    for (int det = 1; det <= 9; det++)
    {
        Tree_Peaks["207Bi"][0][det]->Write();
        Tree_Peaks["90Sr"][0][det]->Write();
    }
    f_tree->Close();
    MATCHED_File->cd();
    Info("Write Tree end");
}

void InitHistogramsONOFF()
{
    MATCHED_File->cd();
    Info("Init Histograms ON OFF start");
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBetaHigh(i))
        {
            H_SiPM_ONOFF[i].first = new TH1D(("H_SiPM_ONOFF_" + detectorName[i] + "_before").c_str(), ("H_SiPM_ONOFF_" + detectorName[i] + "_before").c_str(), eHighN, eHighMin, eHighMax);
            H_SiPM_ONOFF[i].first->GetXaxis()->SetTitle("Channel");
            H_SiPM_ONOFF[i].first->GetYaxis()->SetTitle("Counts");
            H_SiPM_ONOFF[i].first->GetXaxis()->CenterTitle();
            H_SiPM_ONOFF[i].first->GetYaxis()->CenterTitle();

            H_SiPM_ONOFF[i].second = new TH1D(("H_SiPM_ONOFF_" + detectorName[i] + "_after").c_str(), ("H_SiPM_ONOFF_" + detectorName[i] + "_after").c_str(), eHighN, eHighMin, eHighMax);
            H_SiPM_ONOFF[i].second->GetXaxis()->SetTitle("Channel");
            H_SiPM_ONOFF[i].second->GetYaxis()->SetTitle("Counts");
            H_SiPM_ONOFF[i].second->GetXaxis()->CenterTitle();
            H_SiPM_ONOFF[i].second->GetYaxis()->CenterTitle();

            H_SiPM_ONOFF_Corrected[i] = new TH1D(("H_SiPM_ONOFF_before_Corrected_" + detectorName[i]).c_str(), ("H_SiPM_ONOFF_before_Corrected_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
            H_SiPM_ONOFF_Corrected[i]->GetXaxis()->SetTitle("Channel");
            H_SiPM_ONOFF_Corrected[i]->GetYaxis()->SetTitle("Counts");
            H_SiPM_ONOFF_Corrected[i]->GetXaxis()->CenterTitle();
            H_SiPM_ONOFF_Corrected[i]->GetYaxis()->CenterTitle();

            Tree_SiPM_ONOFF[i] = new TTree(("Tree_SiPM_ONOFF_before_" + detectorName[i]).c_str(), ("Tree_SiPM_ONOFF_before_" + detectorName[i]).c_str());
            Tree_SiPM_ONOFF[i]->Branch("Channel", &Channel);
        }

        if (IsDetectorBetaLow(i))
        {
            H_SiPM_ONOFF[i].first = new TH1D(("H_SiPM_ONOFF_" + detectorName[i] + "_before").c_str(), ("H_SiPM_ONOFF_" + detectorName[i] + "_before").c_str(), eLowN, eLowMin, eLowMax);
            H_SiPM_ONOFF[i].first->GetXaxis()->SetTitle("Channel");
            H_SiPM_ONOFF[i].first->GetYaxis()->SetTitle("Counts");
            H_SiPM_ONOFF[i].first->GetXaxis()->CenterTitle();
            H_SiPM_ONOFF[i].first->GetYaxis()->CenterTitle();

            H_SiPM_ONOFF[i].second = new TH1D(("H_SiPM_ONOFF_" + detectorName[i] + "_after").c_str(), ("H_SiPM_ONOFF_" + detectorName[i] + "_after").c_str(), eLowN, eLowMin, eLowMax);
            H_SiPM_ONOFF[i].second->GetXaxis()->SetTitle("Channel");
            H_SiPM_ONOFF[i].second->GetYaxis()->SetTitle("Counts");
            H_SiPM_ONOFF[i].second->GetXaxis()->CenterTitle();
            H_SiPM_ONOFF[i].second->GetYaxis()->CenterTitle();

            H_SiPM_ONOFF_Corrected[i] = new TH1D(("H_SiPM_ONOFF_before_Corrected_" + detectorName[i]).c_str(), ("H_SiPM_ONOFF_before_Corrected_" + detectorName[i]).c_str(), eLowN, eLowMin, eLowMax);
            H_SiPM_ONOFF_Corrected[i]->GetXaxis()->SetTitle("Channel");
            H_SiPM_ONOFF_Corrected[i]->GetYaxis()->SetTitle("Counts");
            H_SiPM_ONOFF_Corrected[i]->GetXaxis()->CenterTitle();
            H_SiPM_ONOFF_Corrected[i]->GetYaxis()->CenterTitle();

            Tree_SiPM_ONOFF[i] = new TTree(("Tree_SiPM_ONOFF_before_" + detectorName[i]).c_str(), ("Tree_SiPM_ONOFF_before_" + detectorName[i]).c_str());
            Tree_SiPM_ONOFF[i]->Branch("Channel", &Channel);
        }
    }    
}

void InitHistograms(string run, int Verbose = 0, bool matched = false)
{
    string add;
    if (matched) 
    {
        add = "_Matched";
        eLowMax = eHighMax;
        eLowMin = eHighMin;
        eLowN = eHighN;
    }
    else add = "";

    Info("Init Histograms start");

    if (Verbose == 1)
        Info(run, 1);
    if (!matched)
        dir[run] = MATCHED_File->mkdir(run.c_str());
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBetaHigh(i))
        {
            if (Verbose == 1)
                Info(detectorName[i], 2);

            if (!matched)
                dir_detector[run][i] = dir[run]->mkdir(("SiPM_" + to_string(GetDetectorChannel(i))).c_str());

            H_SiPM_High[run][i] = new TH1D(("H_SiPM_High_" + run + "_" + detectorName[i] + add).c_str(), ("H_SiPM_High_" + run + "_" + detectorName[i] + add).c_str(), eHighN, eHighMin, eHighMax);
            H_SiPM_High[run][i]->GetXaxis()->SetTitle("Channel");
            H_SiPM_High[run][i]->GetYaxis()->SetTitle("Counts");
            H_SiPM_High[run][i]->GetXaxis()->CenterTitle();
            H_SiPM_High[run][i]->GetYaxis()->CenterTitle();
            H_SiPM_High[run][i]->SetStats(0);

            H_SiPM_HighLow[run][i] = new TH2D(("H_SiPM_HighLow_" + run + "_" + detectorName[i] + add).c_str(), ("H_SiPM_HighLow_" + run + "_" + detectorName[i] + add).c_str(), eLowN, eLowMin, eLowMax, eHighN, eHighMin, eHighMax);
            H_SiPM_HighLow[run][i]->GetXaxis()->SetTitle("Channel");
            H_SiPM_HighLow[run][i]->GetYaxis()->SetTitle("Channel");
            H_SiPM_HighLow[run][i]->GetXaxis()->CenterTitle();
            H_SiPM_HighLow[run][i]->GetYaxis()->CenterTitle();
            H_SiPM_HighLow[run][i]->SetStats(0);

            G_SiPM_HighLow[run][i] = new TGraph();

            H_SiPM_Gain[run][i] = new TH2D(("H_SiPM_Gain_" + run + "_" + detectorName[i] + add).c_str(), ("H_SiPM_Gain_" + run + "_" + detectorName[i] + add).c_str(), eHighN, eHighMin, eHighMax, eHighN, eHighMin, eHighMax);
            H_SiPM_Gain[run][i]->GetXaxis()->SetTitle("Channel");
            H_SiPM_Gain[run][i]->GetYaxis()->SetTitle("SiPM 4 Channel");
            H_SiPM_Gain[run][i]->GetXaxis()->CenterTitle();
            H_SiPM_Gain[run][i]->GetYaxis()->CenterTitle();
            H_SiPM_Gain[run][i]->SetStats(0);

            G_SiPM_Gain[run][i] = new TGraph();
        }

        if (IsDetectorBetaLow(i))
        {
            if (Verbose == 1)
                Info(detectorName[i], 2);
            H_SiPM_Low[run][i] = new TH1D(("H_SiPM_Low_" + run + "_" + detectorName[i] + add).c_str(), ("H_SiPM_Low_" + run + "_" + detectorName[i] + add).c_str(), eLowN, eLowMin, eLowMax);
            H_SiPM_Low[run][i]->GetXaxis()->SetTitle("Channel");
            H_SiPM_Low[run][i]->GetYaxis()->SetTitle("Counts");
            H_SiPM_Low[run][i]->GetXaxis()->CenterTitle();
            H_SiPM_Low[run][i]->GetYaxis()->CenterTitle();

            H_SiPM_Gain[run][i] = new TH2D(("H_SiPM_Gain_" + run + "_" + detectorName[i] + add).c_str(), ("H_SiPM_Gain_" + run + "_" + detectorName[i] + add).c_str(), eLowN, eLowMin, eLowMax, eLowN, eLowMin, eLowMax);
            H_SiPM_Gain[run][i]->GetXaxis()->SetTitle("Channel");
            H_SiPM_Gain[run][i]->GetYaxis()->SetTitle("SiPM 4 Channel");
            H_SiPM_Gain[run][i]->GetXaxis()->CenterTitle();
            H_SiPM_Gain[run][i]->GetYaxis()->CenterTitle();

            G_SiPM_Gain[run][i] = new TGraph();
        }
    }

    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBetaHigh(i))
        {
            gHighLow_RUNS_a[i] = new TGraphErrors();
            gHighLow_RUNS_b[i] = new TGraphErrors();
            gSiPM_High_RUNS_a[i] = new TGraphErrors();
            gSiPM_High_RUNS_b[i] = new TGraphErrors();
        }
        if (IsDetectorBetaLow(i))
        {
            gSiPM_Low_RUNS_a[i] = new TGraphErrors();
            gSiPM_Low_RUNS_b[i] = new TGraphErrors();
        }
    }

    Info("Init Histograms done");
}

void ReadDataONOFF(int current_run)
{
    // fillingHistograms for M>=3
    Info("Reading Data for Fitting ON-OFF", 1);
    clock_t start = clock(), Current;
    Reader->Restart();
    int MAXENTRIES = Reader->GetEntries();

    while (Reader->Next() && Reader->GetCurrentEntry() < Entry_MAX)
    {
        ProgressBar(Reader->GetCurrentEntry(), MAXENTRIES, start, Current, "Reading Tree");

        // In case of Silicon gate
        if (Silicon != nullptr)
        {
            if ((**HRS).isValid)
                continue;

            int sili_label = (*Silicon)[1].Label;
            double sili_energy = F_SiliconCalibration[sili_label]->Eval((*Silicon)[1].Channel / 1000.);
            if (sili_energy < WindowsMap[NUCLEI][IAS[NUCLEI]][sili_label].first || sili_energy > WindowsMap[NUCLEI][IAS[NUCLEI]][sili_label].second)
                continue;
        }

        for (int i_groups = 0; i_groups < (**SiPM_Groups).size(); i_groups++)
        {
            if ((**SiPM_Groups)[i_groups].size() <= 2)
                continue;
            for (int i_pair = 0; i_pair < (**SiPM_Groups)[i_groups].size(); i_pair++)
            {
                if ((**SiPM_Groups)[i_groups][i_pair].first.isValid)
                {
                    if (current_run < 80)
                    {
                        H_SiPM_ONOFF[(**SiPM_Groups)[i_groups][i_pair].first.Label].first->Fill((**SiPM_Groups)[i_groups][i_pair].first.Channel);
                        Channel = (**SiPM_Groups)[i_groups][i_pair].first.Channel;
                        Tree_SiPM_ONOFF[(**SiPM_Groups)[i_groups][i_pair].first.Label]->Fill();
                    }
                    else
                        H_SiPM_ONOFF[(**SiPM_Groups)[i_groups][i_pair].first.Label].second->Fill((**SiPM_Groups)[i_groups][i_pair].first.Channel);
                    
                    // per run
                    H_SiPM_ONOFF_RUN[(**SiPM_Groups)[i_groups][i_pair].first.Label][current_run]->Fill((**SiPM_Groups)[i_groups][i_pair].first.Channel);
                    Channel = (**SiPM_Groups)[i_groups][i_pair].first.Channel;
                    Tree_SiPM_ONOFF_RUN[(**SiPM_Groups)[i_groups][i_pair].first.Label][current_run]->Fill();
                }

                if ((**SiPM_Groups)[i_groups][i_pair].second.isValid)
                {
                    if (current_run < 80)
                    {
                        H_SiPM_ONOFF[(**SiPM_Groups)[i_groups][i_pair].second.Label].first->Fill((**SiPM_Groups)[i_groups][i_pair].second.Channel);
                        Channel = (**SiPM_Groups)[i_groups][i_pair].second.Channel;
                        Tree_SiPM_ONOFF[(**SiPM_Groups)[i_groups][i_pair].second.Label]->Fill();
                    }
                    else
                        H_SiPM_ONOFF[(**SiPM_Groups)[i_groups][i_pair].second.Label].second->Fill((**SiPM_Groups)[i_groups][i_pair].second.Channel);
                    
                    // per run
                    H_SiPM_ONOFF_RUN[(**SiPM_Groups)[i_groups][i_pair].second.Label][current_run]->Fill((**SiPM_Groups)[i_groups][i_pair].second.Channel);
                    Channel = (**SiPM_Groups)[i_groups][i_pair].second.Channel;
                    Tree_SiPM_ONOFF_RUN[(**SiPM_Groups)[i_groups][i_pair].second.Label][current_run]->Fill();
                }
            }
        }
    }
}

void ReadData()
{
    Info("Reading Data for Fitting Low-High", 1);
    clock_t start = clock(), Current;

    // Graph counter
    int counter[SIGNAL_MAX] = {0};
    int counterSiPM[SIGNAL_MAX] = {0};
    Reader->Restart();
    int MAXENTRIES = Reader->GetEntries();

    while (Reader->Next() && Reader->GetCurrentEntry() < Entry_MAX)
    {
        ProgressBar(Reader->GetCurrentEntry(), MAXENTRIES, start, Current, "Reading Tree");

        // In case of Silicon gate
        if (Silicon != nullptr)
        {
            if ((**HRS).isValid)
                continue;

            int sili_label = (*Silicon)[1].Label;
            double sili_energy = F_SiliconCalibration[sili_label]->Eval((*Silicon)[1].Channel / 1000.);
            // gate on proton peak
            // int peak=0;
            // if (NUCLEI == "32Ar") peak = 14;
            // if (NUCLEI == "33Ar") peak = 21;
            if (sili_energy < WindowsMap[NUCLEI][IAS[NUCLEI]][sili_label].first || sili_energy > WindowsMap[NUCLEI][IAS[NUCLEI]][sili_label].second)
                continue;
        }

        // Applying ONOFF
        for (int i_groups = 0; i_groups < (**SiPM_Groups).size(); i_groups++)
        {
            for (int i_pair = 0; i_pair < (**SiPM_Groups)[i_groups].size(); i_pair++)
            {
                
                if ((**SiPM_Groups)[i_groups][i_pair].first.isValid)
                    (**SiPM_Groups)[i_groups][i_pair].first.Channel = SiPM_ONOFF_Coef_RUN[(**SiPM_Groups)[i_groups][i_pair].first.Label+10][stoi(RUN)]->Eval((**SiPM_Groups)[i_groups][i_pair].first.Channel);
                if ((**SiPM_Groups)[i_groups][i_pair].second.isValid)
                    (**SiPM_Groups)[i_groups][i_pair].second.Channel = SiPM_ONOFF_Coef_RUN[(**SiPM_Groups)[i_groups][i_pair].second.Label][stoi(RUN)]->Eval((**SiPM_Groups)[i_groups][i_pair].second.Channel);
            }
        }

        // Looping on SiPM groups
        for (int i_groups = 0; i_groups < (**SiPM_Groups).size(); i_groups++)
        {
            if ((**SiPM_Groups)[i_groups].size() <= 2)
                continue;
            for (int i_pair = 0; i_pair < (**SiPM_Groups)[i_groups].size(); i_pair++)
            {
                // Compare High Low for each SiPM
                if ((**SiPM_Groups)[i_groups][i_pair].first.isValid && (**SiPM_Groups)[i_groups][i_pair].second.isValid)
                {
                    if ((**SiPM_Groups)[i_groups][i_pair].first.Channel > 0 && (**SiPM_Groups)[i_groups][i_pair].second.Channel > 0)
                    {
                    // if (!PileUp((**SiPM_Groups)[i_groups][i_pair].first) && !PileUp((**SiPM_Groups)[i_groups][i_pair].second)) // pileup high or low
                    // {
                        H_SiPM_HighLow[RUN][(**SiPM_Groups)[i_groups][i_pair].first.Label]->Fill((**SiPM_Groups)[i_groups][i_pair].second.Channel, (**SiPM_Groups)[i_groups][i_pair].first.Channel);
                        G_SiPM_HighLow[RUN][(**SiPM_Groups)[i_groups][i_pair].first.Label]->SetPoint(counter[(**SiPM_Groups)[i_groups][i_pair].first.Label], (**SiPM_Groups)[i_groups][i_pair].second.Channel, (**SiPM_Groups)[i_groups][i_pair].first.Channel);
                        counter[(**SiPM_Groups)[i_groups][i_pair].first.Label]++;
                    // }
                    }   
                }

                if ((**SiPM_Groups)[i_groups][i_pair].first.isValid)
                {
                    H_SiPM_High[RUN][(**SiPM_Groups)[i_groups][i_pair].first.Label]->Fill((**SiPM_Groups)[i_groups][i_pair].first.Channel);
                }

                if ((**SiPM_Groups)[i_groups][i_pair].second.isValid)
                {
                    H_SiPM_Low[RUN][(**SiPM_Groups)[i_groups][i_pair].second.Label]->Fill((**SiPM_Groups)[i_groups][i_pair].second.Channel);
                }

                // Compare SiPMS to 101
                if ((**SiPM_Groups)[i_groups][i_pair].first.Label == 101)
                {
                    for (int i_pair_compare = 0; i_pair_compare < (**SiPM_Groups)[i_groups].size(); i_pair_compare++)
                    {
                        if ((**SiPM_Groups)[i_groups][i_pair].first.isValid && (**SiPM_Groups)[i_groups][i_pair_compare].first.isValid)
                        {
                            if ((**SiPM_Groups)[i_groups][i_pair].first.Channel > 0 && (**SiPM_Groups)[i_groups][i_pair_compare].first.Channel > 0)
                            {
                                // if (!PileUp((**SiPM_Groups)[i_groups][i_pair].first) && !PileUp((**SiPM_Groups)[i_groups][i_pair_compare].first)) // pileup high or low
                                // {
                                H_SiPM_Gain[RUN][(**SiPM_Groups)[i_groups][i_pair_compare].first.Label]->Fill((**SiPM_Groups)[i_groups][i_pair].first.Channel, (**SiPM_Groups)[i_groups][i_pair_compare].first.Channel);
                                G_SiPM_Gain[RUN][(**SiPM_Groups)[i_groups][i_pair_compare].first.Label]->SetPoint(counterSiPM[(**SiPM_Groups)[i_groups][i_pair_compare].first.Label], (**SiPM_Groups)[i_groups][i_pair].first.Channel, (**SiPM_Groups)[i_groups][i_pair_compare].first.Channel);
                                counterSiPM[(**SiPM_Groups)[i_groups][i_pair_compare].first.Label]++;
                                // }
                            }
                        }

                        if ((**SiPM_Groups)[i_groups][i_pair].second.isValid && (**SiPM_Groups)[i_groups][i_pair_compare].second.isValid)
                        {
                            if ((**SiPM_Groups)[i_groups][i_pair].second.Channel > 0 && (**SiPM_Groups)[i_groups][i_pair_compare].second.Channel > 0)
                            {
                                // if (!PileUp((**SiPM_Groups)[i_groups][i_pair].second) && !PileUp((**SiPM_Groups)[i_groups][i_pair_compare].second)) // pileup high or low
                                // {
                                H_SiPM_Gain[RUN][(**SiPM_Groups)[i_groups][i_pair_compare].second.Label]->Fill((**SiPM_Groups)[i_groups][i_pair].second.Channel, (**SiPM_Groups)[i_groups][i_pair_compare].second.Channel);
                                G_SiPM_Gain[RUN][(**SiPM_Groups)[i_groups][i_pair_compare].second.Label]->SetPoint(counterSiPM[(**SiPM_Groups)[i_groups][i_pair_compare].second.Label], (**SiPM_Groups)[i_groups][i_pair].second.Channel, (**SiPM_Groups)[i_groups][i_pair_compare].second.Channel);
                                counterSiPM[(**SiPM_Groups)[i_groups][i_pair_compare].second.Label]++;
                                // }
                            }
                        }
                    }
                }
            }
        }
    }
}

void ReadDataWithCorrections()
{
    Info("Reading Data for applying corrections Low-High", 1);
    clock_t start = clock(), Current;

    // Graph counter
    int counter[SIGNAL_MAX] = {0};
    int counterSiPM[SIGNAL_MAX] = {0};
    Reader->Restart();

    

    while (Reader->Next() && Reader->GetCurrentEntry() < Entry_MAX)
    {
        ProgressBar(Reader->GetCurrentEntry(), Reader->GetEntries(), start, Current, "Reading Tree");

        // In case of Silicon gate
        if (Silicon)
        {
            if ((**HRS).isValid)
                continue;
            
            int sili_label = (*Silicon)[1].Label;
            double sili_energy = F_SiliconCalibration[sili_label]->Eval((*Silicon)[1].Channel / 1000.);
            // gate on proton peak
            if (sili_energy < WindowsMap[NUCLEI][IAS[NUCLEI]][sili_label].first || sili_energy > WindowsMap[NUCLEI][IAS[NUCLEI]][sili_label].second)
                continue;
        }

        // cout << "Applying correction" << endl;
        /// old
        // for (int i_groups = 0; i_groups < (**SiPM_Groups).size(); i_groups++)
        // {
        //     // cout << "Group " << i_groups << endl;
        //     for (int i_pair = 0; i_pair < (**SiPM_Groups)[i_groups].size(); i_pair++)
        //     {
        //         // APPLYING MATCHING ONOFF
        //         if ((**SiPM_Groups)[i_groups][i_pair].first.isValid)
        //             (**SiPM_Groups)[i_groups][i_pair].first.Channel = SiPM_ONOFF_Coef_RUN[(**SiPM_Groups)[i_groups][i_pair].first.Label+10][stoi(RUN)]->Eval((**SiPM_Groups)[i_groups][i_pair].first.Channel);
        //         if ((**SiPM_Groups)[i_groups][i_pair].second.isValid)
        //             (**SiPM_Groups)[i_groups][i_pair].second.Channel = SiPM_ONOFF_Coef_RUN[(**SiPM_Groups)[i_groups][i_pair].second.Label][stoi(RUN)]->Eval((**SiPM_Groups)[i_groups][i_pair].second.Channel);
               
        //         // APPLYING MATCHING High Low
        //         if ((**SiPM_Groups)[i_groups][i_pair].second.isValid)
        //         {
        //             // cout << "Applying on " << (**SiPM_Groups)[i_groups][i_pair].second.Label << endl;
        //             (**SiPM_Groups)[i_groups][i_pair].second.Channel = MatchingLowHigh[RUN][(**SiPM_Groups)[i_groups][i_pair].second.Label-10]->Eval((**SiPM_Groups)[i_groups][i_pair].second.Channel);
        //         }

        //         // APPLYING MATCHING SiPM
        //         if ((**SiPM_Groups)[i_groups][i_pair].first.isValid)
        //         {
        //             // cout << "Applying on " << (**SiPM_Groups)[i_groups][i_pair].first.Label << endl;
        //             (**SiPM_Groups)[i_groups][i_pair].first.Channel = MatchingSiPM[RUN][(**SiPM_Groups)[i_groups][i_pair].first.Label]->Eval((**SiPM_Groups)[i_groups][i_pair].first.Channel);
        //         }
        //         if ((**SiPM_Groups)[i_groups][i_pair].second.isValid)
        //         {
        //             // cout << "Applying on " << (**SiPM_Groups)[i_groups][i_pair].second.Label << endl;
        //             (**SiPM_Groups)[i_groups][i_pair].second.Channel = MatchingSiPM[RUN][(**SiPM_Groups)[i_groups][i_pair].second.Label-10]->Eval((**SiPM_Groups)[i_groups][i_pair].second.Channel);
        //         }
        //     }
        // }

        // new
        for (int i_groups = 0; i_groups < (**SiPM_Groups).size(); i_groups++)
        {
            // cout << "Group " << i_groups << endl;
            for (int i_pair = 0; i_pair < (**SiPM_Groups)[i_groups].size(); i_pair++)
            {
                // APPLYING MATCHING ONOFF
                if (stoi(RUN)<80)
                {
                    if ((**SiPM_Groups)[i_groups][i_pair].first.isValid)
                        (**SiPM_Groups)[i_groups][i_pair].first.Channel = SiPM_ONOFF[(**SiPM_Groups)[i_groups][i_pair].first.Label]->Eval((**SiPM_Groups)[i_groups][i_pair].first.Channel);
                    if ((**SiPM_Groups)[i_groups][i_pair].second.isValid)
                        (**SiPM_Groups)[i_groups][i_pair].second.Channel = SiPM_ONOFF[(**SiPM_Groups)[i_groups][i_pair].second.Label]->Eval((**SiPM_Groups)[i_groups][i_pair].second.Channel);
                }
                // APPLYING MATCHING High Low
                if ((**SiPM_Groups)[i_groups][i_pair].second.isValid)
                {
                    // cout << "Applying on " << (**SiPM_Groups)[i_groups][i_pair].second.Label << endl;
                    (**SiPM_Groups)[i_groups][i_pair].second.Channel = SiPM_HighLow[(**SiPM_Groups)[i_groups][i_pair].second.Label]->Eval((**SiPM_Groups)[i_groups][i_pair].second.Channel);
                }

                // APPLYING MATCHING SiPM
                if ((**SiPM_Groups)[i_groups][i_pair].first.isValid)
                {
                    // cout << "Applying on " << (**SiPM_Groups)[i_groups][i_pair].first.Label << endl;
                    (**SiPM_Groups)[i_groups][i_pair].first.Channel = SiPM_Gain[(**SiPM_Groups)[i_groups][i_pair].first.Label]->Eval((**SiPM_Groups)[i_groups][i_pair].first.Channel);
                }
                if ((**SiPM_Groups)[i_groups][i_pair].second.isValid)
                {
                    // cout << "Applying on " << (**SiPM_Groups)[i_groups][i_pair].second.Label << endl;
                    (**SiPM_Groups)[i_groups][i_pair].second.Channel = SiPM_Gain[(**SiPM_Groups)[i_groups][i_pair].second.Label]->Eval((**SiPM_Groups)[i_groups][i_pair].second.Channel);
                }
            }
        }
        // cout << "Correction applied" << endl;
        

        // cout << "Filling histograms" << endl;
        // Looping on SiPM groups
        for (int i_groups = 0; i_groups < (**SiPM_Groups).size(); i_groups++)
        {
            // cout << "Group " << i_groups << endl;
            if ((**SiPM_Groups)[i_groups].size() <= 2)
                continue;
            for (int i_pair = 0; i_pair < (**SiPM_Groups)[i_groups].size(); i_pair++)
            {
                // cout << "Pair " << i_pair << endl;
                // Compare High Low for each SiPM
                if ((**SiPM_Groups)[i_groups][i_pair].first.isValid && (**SiPM_Groups)[i_groups][i_pair].second.isValid)
                {
                    H_SiPM_HighLow[RUN][(**SiPM_Groups)[i_groups][i_pair].first.Label]->Fill((**SiPM_Groups)[i_groups][i_pair].second.Channel, (**SiPM_Groups)[i_groups][i_pair].first.Channel);
                    G_SiPM_HighLow[RUN][(**SiPM_Groups)[i_groups][i_pair].first.Label]->SetPoint(counter[(**SiPM_Groups)[i_groups][i_pair].first.Label], (**SiPM_Groups)[i_groups][i_pair].second.Channel, (**SiPM_Groups)[i_groups][i_pair].first.Channel);
                    counter[(**SiPM_Groups)[i_groups][i_pair].first.Label]++;
                }

                if ((**SiPM_Groups)[i_groups][i_pair].first.isValid)
                {
                    H_SiPM_High[RUN][(**SiPM_Groups)[i_groups][i_pair].first.Label]->Fill((**SiPM_Groups)[i_groups][i_pair].first.Channel);
                }

                if ((**SiPM_Groups)[i_groups][i_pair].second.isValid)
                {

                    H_SiPM_Low[RUN][(**SiPM_Groups)[i_groups][i_pair].second.Label]->Fill((**SiPM_Groups)[i_groups][i_pair].second.Channel);
                }

                // cout << "SiPM 104 comparison" << endl;
                // Compare SiPMS to 104
                if ((**SiPM_Groups)[i_groups][i_pair].first.Label == 101)
                {
                    for (int i_pair_compare = 0; i_pair_compare < (**SiPM_Groups)[i_groups].size(); i_pair_compare++)
                    {
                        // cout << "Pair compare " << i_pair_compare << endl;
                        if ((**SiPM_Groups)[i_groups][i_pair].first.isValid && (**SiPM_Groups)[i_groups][i_pair_compare].first.isValid)
                        {
                            H_SiPM_Gain[RUN][(**SiPM_Groups)[i_groups][i_pair_compare].first.Label]->Fill((**SiPM_Groups)[i_groups][i_pair].first.Channel , (**SiPM_Groups)[i_groups][i_pair_compare].first.Channel);
                            G_SiPM_Gain[RUN][(**SiPM_Groups)[i_groups][i_pair_compare].first.Label]->SetPoint(counterSiPM[(**SiPM_Groups)[i_groups][i_pair_compare].first.Label], (**SiPM_Groups)[i_groups][i_pair].first.Channel, (**SiPM_Groups)[i_groups][i_pair_compare].first.Channel);
                            counterSiPM[(**SiPM_Groups)[i_groups][i_pair_compare].first.Label]++;
                        }

                        if ((**SiPM_Groups)[i_groups][i_pair].second.isValid && (**SiPM_Groups)[i_groups][i_pair_compare].second.isValid)
                        {
                            H_SiPM_Gain[RUN][(**SiPM_Groups)[i_groups][i_pair_compare].second.Label]->Fill((**SiPM_Groups)[i_groups][i_pair].second.Channel , (**SiPM_Groups)[i_groups][i_pair_compare].second.Channel);
                            G_SiPM_Gain[RUN][(**SiPM_Groups)[i_groups][i_pair_compare].second.Label]->SetPoint(counterSiPM[(**SiPM_Groups)[i_groups][i_pair_compare].second.Label], (**SiPM_Groups)[i_groups][i_pair].second.Channel, (**SiPM_Groups)[i_groups][i_pair_compare].second.Channel);
                            counterSiPM[(**SiPM_Groups)[i_groups][i_pair_compare].second.Label]++;
                        }
                    }
                }
            }
        }
        // cout << "Histograms filled" << endl;
    }
}

void FittingONOFF(int Verbose = 0)
{
    Info("Fitting ON-OFF gain factor");
    int rebin = 4;

    double start = 200e3;
    double end = 3000e3;
    
    TFile *f = MyTFile(Form("test_ONOFF_%d.root", YEAR), "RECREATE");
    f->cd();
    
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBeta(i))
        {

            if (IsDetectorBetaLow(i))
            {
                start = 20e3;
                end = 600e3;
            }

            Info(detectorName[i], 1);
            TGraphErrors *G_chi2 = new TGraphErrors();
            H_SiPM_ONOFF[i].second->Rebin(rebin);
            for (double coef = 0.85; coef <= 1.2; coef += 0.002)
            {
                TH1D *H_temp = (TH1D *)H_SiPM_ONOFF[i].first->Clone("H_temp");
                H_temp->Reset();
                Reader = new TTreeReader(Tree_SiPM_ONOFF[i]);
                TTreeReaderValue<double> Channeli(*Reader, "Channel");
                while (Reader->Next())
                {
                    double corrected_channel = *Channeli * coef;
                    H_temp->Fill(corrected_channel);
                }

                // scaling 
                
                H_temp->Rebin(rebin);
                H_SiPM_ONOFF[i].second->GetXaxis()->SetRangeUser(start, end);
                H_temp->GetXaxis()->SetRangeUser(start, end);
                H_temp->Scale(H_SiPM_ONOFF[i].second->Integral() / H_temp->Integral());

                double chi2 = H_temp->Chi2Test(H_SiPM_ONOFF[i].second, "CHI2/NDF");

                // cout << "        Coef: " << coef << " Chi2: " << chi2 << endl;

                G_chi2->SetPoint(G_chi2->GetN(), coef, chi2);
            }
            
            TCanvas *cChi2 = new TCanvas(("cChi2_ONOFF_" + detectorName[i]).c_str(), ("cChi2_ONOFF_" + detectorName[i]).c_str(), 800, 600);
            cChi2->Divide(2, 1);
            cChi2->cd(1);
            G_chi2->SetTitle(("Chi2 Test ON-OFF SiPM " + detectorName[i]).c_str());
            G_chi2->GetXaxis()->SetTitle("Gain Factor");
            G_chi2->GetYaxis()->SetTitle("Chi2/NDF");
            

            double x0 = 0;
            double y0 = 1e6;
            for (int i = 0; i < G_chi2->GetN(); i++)
            {
                double x, y;
                G_chi2->GetPoint(i, x, y);
                if ((y < y0 || y0 == 1e6) && x < 1.1 && x > 0.9)
                {
                    x0 = x;
                    y0 = y;
                }
            }
           

            TF1 *f_pol2_2 = new TF1("AsymetricPol2_2", AsymetricPol2, x0 - 0.05, x0 + 0.05, 4);
            f_pol2_2->SetParameter(0, x0);
            f_pol2_2->SetParLimits(0, x0 - 0.05, x0 + 0.05);
            f_pol2_2->SetParLimits(1, 0, 10000);
            f_pol2_2->SetParLimits(2, 0, 10000);
            f_pol2_2->SetParameter(3, y0);
            f_pol2_2->SetParLimits(3, 0, 10*y0);
            
            TFitResultPtr r = G_chi2->Fit(f_pol2_2, "RS", "", x0 - 0.05, x0 + 0.05);
            G_chi2->SetMarkerStyle(20);
            G_chi2->Draw("AP");

            TLine *line1 = new TLine(x0, 0, x0, y0);
            line1->SetLineColor(kGray);
            line1->SetLineStyle(2);
            line1->Draw("same");

            if (IsDetectorBetaLow(i))
                x0 = f_pol2_2->GetParameter(0);

            for (string run : Map_RunFiles["32Ar"])
            {
                if (stoi(run) < 80)
                {
                    SiPM_ONOFF_Coef_RUN[i][stoi(run)] = new TF1("f_linear", "[0]*x", 0, 1e7);
                    SiPM_ONOFF_Coef_RUN[i][stoi(run)]->SetParameter(0, x0);
                }
            }


            TLine *line = new TLine(x0, 0, x0, f_pol2_2->Eval(x0));
            line->SetLineColor(kRed);
            line->SetLineStyle(2);
            line->Draw("same");

            cChi2->cd(2);
            TH1D *H_temp = (TH1D *)H_SiPM_ONOFF[i].first->Clone("H_temp");
            H_temp->Reset();
            Reader = new TTreeReader(Tree_SiPM_ONOFF[i]);
            TTreeReaderValue<double> Channeli(*Reader, "Channel");
            while (Reader->Next())
            {
                double corrected_channel = *Channeli * x0;
                H_temp->Fill(corrected_channel);
            }

            // scaling
            H_SiPM_ONOFF[i].first->Rebin(rebin);
            H_temp->Rebin(rebin);
            H_SiPM_ONOFF[i].first->GetXaxis()->SetRangeUser(start, end);
            H_SiPM_ONOFF[i].second->GetXaxis()->SetRangeUser(start, end);
            H_temp->GetXaxis()->SetRangeUser(start, end);
            H_temp->Scale(H_SiPM_ONOFF[i].second->Integral() / H_temp->Integral());
            H_SiPM_ONOFF[i].first->Scale(H_SiPM_ONOFF[i].second->Integral() / H_SiPM_ONOFF[i].first->Integral());

            H_SiPM_ONOFF[i].first->Draw("HIST");
            H_SiPM_ONOFF[i].second->SetLineColor(kBlack);
            H_SiPM_ONOFF[i].second->Draw("HIST SAME");
            H_temp->SetLineColor(kRed);
            H_temp->Draw("HIST SAME");

            cChi2->Write();
        }
    }
    f->Close();
}

void FittingONOFF_Runs(int Verbose = 0)
{
    Info("Fitting ON-OFF gain factor per Run");
    int rebin = 10;

    double start = 200e3;
    double end = 6000e3;
    
    TFile *f = MyTFile(Form("test_ONOFF_run_%d.root", YEAR), "RECREATE");
    f->cd();

    // higher stat run for comparaison - selecting with Hi1 detector
    double entries_max = 0;
    int run_ref = -1;
    TH1D * H_ref[SIGNAL_MAX];
    for (string run : Map_RunFiles["32Ar"])
    {
        double entries = H_SiPM_ONOFF_RUN[101][stoi(run)]->GetEntries();
        cout << "Run " << run << " Entries Hi1: " << entries << endl;
        if (entries > entries_max)
        {
            entries_max = entries;
            run_ref = stoi(run);

            for (int i = 0; i < SIGNAL_MAX; i++)
            {
                if (IsDetectorBeta(i))
                {
                    H_ref[i] = (TH1D *)H_SiPM_ONOFF_RUN[i][stoi(run)]->Clone("H_ref_ONOFF");
                    H_ref[i]->Rebin(rebin);
                }
            }
        }
    }

    Info("Reference run for ON-OFF fitting: " + to_string(run_ref));

    
    // looping over runs   
    TGraphErrors *G_factor[SIGNAL_MAX] = {nullptr};
    TDirectory *dir_run[SIGNAL_MAX] = {nullptr};
    double xmin = 0.8;
    double xmax = 1.6;
    for (string run : Map_RunFiles["32Ar"])
    { 
        int RUN = stoi(run);
        Info("Fitting Run " + run);
        dir_run[RUN] = f->mkdir(("Run_" + run).c_str());
        dir_run[RUN]->cd();
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorBeta(i))
            {
                Info(detectorName[i], 1);
                TGraphErrors *G_chi2 = new TGraphErrors();
                
                for (double coef = xmin; coef <= xmax; coef += 0.01)
                {
                    TH1D *H_temp = (TH1D *)H_SiPM_ONOFF_RUN[i][RUN]->Clone("H_temp");
                    H_temp->Reset();
                    Reader = new TTreeReader(Tree_SiPM_ONOFF_RUN[i][RUN]);
                    Reader->Restart();
                    TTreeReaderValue<double> Channeli(*Reader, "Channel");
                    while (Reader->Next())
                    {
                        double corrected_channel = *Channeli * coef;
                        H_temp->Fill(corrected_channel);
                    }

                    // scaling 
                    H_temp->Rebin(rebin);
                    H_ref[i]->GetXaxis()->SetRangeUser(start, end);
                    H_temp->GetXaxis()->SetRangeUser(start, end);
                    H_temp->Scale(H_ref[i]->Integral() / H_temp->Integral());

                    double chi2 = H_temp->Chi2Test(H_ref[i], "CHI2/NDF");

                    // cout << "        Coef: " << coef << " Chi2: " << chi2 << endl;

                    G_chi2->AddPoint(coef, chi2);
                }
                
                TCanvas *cChi2 = new TCanvas(("cChi2_ONOFF_" + detectorName[i] + "_Run" + run).c_str(), ("cChi2_ONOFF_" + detectorName[i] + "_Run" + run).c_str(), 800, 600);
                cChi2->Divide(2, 1);
                cChi2->cd(1);
                G_chi2->SetName(("G_chi2_ONOFF_" + detectorName[i] + "_Run" + run).c_str());
                G_chi2->SetTitle(("Chi2 Test ON-OFF SiPM " + detectorName[i] + " Run " + run).c_str());
                G_chi2->GetXaxis()->SetTitle("Gain Factor");
                G_chi2->GetYaxis()->SetTitle("Chi2/NDF");
                

                double x0 = 0;
                double y0 = 1e6;
                for (int i = 0; i < G_chi2->GetN(); i++)
                {
                    double x, y;
                    G_chi2->GetPoint(i, x, y);
                    if ((y < y0 || y0 == 1e6) && x < 1.1 && x > 0.9)
                    {
                        x0 = x;
                        y0 = y;
                    }
                }

                
                gErrorIgnoreLevel = kWarning;
                TF1 *f_pol2_2 = new TF1("AsymetricPol2_2", AsymetricPol2, xmin, xmax, 4);
                f_pol2_2->SetParameter(0, x0);
                f_pol2_2->SetParLimits(0, x0 - 0.1, x0 + 0.1);
                f_pol2_2->SetParLimits(1, 0, 10000);
                f_pol2_2->SetParLimits(2, 0, 10000);
                f_pol2_2->SetParameter(3, y0);
                f_pol2_2->SetParLimits(3, 0, 5*y0);
                
                G_chi2->Fit(f_pol2_2, "RQ", "", x0 - 0.1, x0 + 0.1);
                G_chi2->SetMarkerStyle(20);
                G_chi2->Draw("AP");

                TLine *line1 = new TLine(x0, 0, x0, y0);
                line1->SetLineColor(kGray);
                line1->SetLineStyle(2);
                line1->Draw("same");

                x0 = y0 == 0 ? 1.0 : f_pol2_2->GetParameter(0);
                SiPM_ONOFF_Coef_RUN[i][RUN] = new TF1("f_linear", "[0]*x", 0, 1e7);
                SiPM_ONOFF_Coef_RUN[i][RUN]->SetParameter(0, x0);
                
                TLine *line = new TLine(x0, 0, x0, f_pol2_2->Eval(x0));
                line->SetLineColor(kRed);
                line->SetLineStyle(2);
                line->Draw("same");

                // save value
                if (G_factor[i] == nullptr)
                    G_factor[i] = new TGraphErrors();
                G_factor[i]->AddPoint(RUN, x0);
                Info("Run " + to_string(RUN) + " Gain Factor: " + to_string(x0) + " +/- " + to_string(f_pol2_2->GetParError(0)), 1);
                G_factor[i]->SetPointError(G_factor[i]->GetN()-1, 0, f_pol2_2->GetParError(0));

                cChi2->cd(2);
                TH1D *H_temp = (TH1D *)H_SiPM_ONOFF_RUN[i][RUN]->Clone(Form("H_temp_Run%d_det%d", RUN, i));
                H_temp->Reset();
                Reader = new TTreeReader(Tree_SiPM_ONOFF_RUN[i][RUN]);
                Reader->Restart();
                TTreeReaderValue<double> Channeli(*Reader, "Channel");
                while (Reader->Next())
                {
                    double corrected_channel = *Channeli * x0;
                    H_temp->Fill(corrected_channel);
                }

                // scaling
                H_SiPM_ONOFF_RUN[i][RUN]->Rebin(rebin);
                H_temp->Rebin(rebin);
                H_SiPM_ONOFF_RUN[i][RUN]->GetXaxis()->SetRangeUser(start, end);
                H_ref[i]->GetXaxis()->SetRangeUser(start, end);
                H_temp->GetXaxis()->SetRangeUser(start, end);
                H_temp->Scale(H_ref[i]->Integral() / H_temp->Integral());
                H_SiPM_ONOFF_RUN[i][RUN]->Scale(H_ref[i]->Integral() / H_SiPM_ONOFF_RUN[i][RUN]->Integral());

                H_SiPM_ONOFF_RUN[i][RUN]->Draw("HIST");
                H_ref[i]->SetLineColor(kBlack);
                H_ref[i]->Draw("HIST SAME");
                H_temp->SetLineColor(kRed);
                H_temp->Draw("HIST SAME");

                // cChi2->Write();
            }
        }
    }

    f->cd();
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBeta(i))
        {
            TCanvas *cFactor = new TCanvas(("cFactor_ONOFF_" + detectorName[i]).c_str(), ("cFactor_ONOFF_" + detectorName[i]).c_str(), 800, 600);
            G_factor[i]->SetTitle(("Gain Factor per Run ON-OFF SiPM " + detectorName[i]).c_str());
            G_factor[i]->GetXaxis()->SetTitle("Run Number");
            G_factor[i]->GetYaxis()->SetTitle("Gain Factor");
            G_factor[i]->SetMarkerStyle(20);
            G_factor[i]->Draw("AP");
            // cFactor->Write();
        }
    }


    f->Close();
}

void FittingLowHigh(int Verbose = 0)
{
    Info("Fitting High-Low gain factor for run " + RUN);
    f_linear->SetParameters(10., 0.);
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBetaHigh(i))
        {
            if (Verbose == 1)
                Info(detectorName[i], 1);
            current_detector = i;
            if (G_SiPM_HighLow[RUN][i]->GetN() > 1e3)
            {
                TFitResultPtr r = G_SiPM_HighLow[RUN][i]->Fit(f_linear, "QRS", "", Range_SiPM_LowHigh.first, Range_SiPM_LowHigh.second);
                MatchingLowHigh[RUN][i] = G_SiPM_HighLow[RUN][i]->GetFunction("f_linear");
                MatchingLowHigh[RUN][i]->SetParameter(1, 0.);
                if (MatchingLowHigh[RUN][i] == nullptr || r != 0)
                {
                    MatchingLowHigh[RUN][i] = new TF1("dummy", "[0]*x", 0, 10e6);
                    MatchingLowHigh[RUN][i]->SetParameter(0, 10);
                    MatchingLowHigh[RUN][i]->SetParError(0, 0);
                    Warning("[Fit Failed] Default function assigned for High-Low SiPM " + detectorName[i] + " in run " + RUN);
                }
            }
            else
            {
                MatchingLowHigh[RUN][i] = new TF1("dummy", "[0]*x", 0, 10e6);
                MatchingLowHigh[RUN][i]->SetParameter(0, 10.);
                MatchingLowHigh[RUN][i]->SetParError(0, 0);
                Warning("[No data] Default function assigned for High-Low SiPM " + detectorName[i] + " in run " + RUN);
            }
        }
    }
}

void FittingSiPM(int Verbose = 0)
{
    Info("Fitting SiPM gain factor for run " + RUN);
    f_linear_fixed->FixParameter(1, 0.);
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBeta(i))
        {
            if (Verbose == 1)
                Info(detectorName[i], 1);
            current_detector = i;
            
            // no data, dummy function
            if (G_SiPM_Gain[RUN][i]->GetN() < 1e3)
            {
                MatchingSiPM[RUN][i] = new TF1("dummy", "[0]*x+[1]", 0, 100e6);
                MatchingSiPM[RUN][i]->SetParameter(0, 1.);
                MatchingSiPM[RUN][i]->SetParError(0, 0);
                MatchingSiPM[RUN][i]->SetParameter(0, 0.);
                Warning("[No data] Default function assigned for SiPM matching " + detectorName[i] + " in run " + RUN);
                continue;
            }

            TFitResultPtr r;
            if (IsDetectorBetaLow(i))
            {
                r = G_SiPM_Gain[RUN][i]->Fit(f_linear_fixed, "QS");
                MatchingSiPM[RUN][i] = G_SiPM_Gain[RUN][i]->GetFunction("f_linear_fixed");
            }
            else
            {
                r = G_SiPM_Gain[RUN][i]->Fit(f_linear_fixed, "QRS", "", 200e3, Range_SiPM_LowHigh.second*10);
                MatchingSiPM[RUN][i] = G_SiPM_Gain[RUN][i]->GetFunction("f_linear_fixed");
            }

            // wrong fit result, dummy function
            if (MatchingSiPM[RUN][i] == nullptr || r != 0)
            {
                MatchingSiPM[RUN][i] = new TF1("dummy", "[0]*x+[1]", 0, 100e6);
                MatchingSiPM[RUN][i]->SetParameter(0, 1.);
                MatchingSiPM[RUN][i]->SetParError(0, 0);
                MatchingSiPM[RUN][i]->SetParameter(0, 0.);
                Warning("[Fit Failed] Default function assigned for SiPMs matching " + detectorName[i] + " in run " + RUN);
                continue;
            }
        }
    }
}

void WriteHistogram(int Verbose = 0, bool correction = false)
{
    string add = "";
    if (correction) add = "_Matched";

    gStyle->SetOptStat(0);
    MATCHED_File->cd();
    Info("Write Histograms for run " + RUN + add);

    dir[RUN]->cd();
    TCanvas *cHighLowSiPM = new TCanvas(("2D_HighLowSiPM_Fitting" + RUN + add).c_str(), ("2D_HighLowSiPM_Fitting" + RUN + add).c_str(), 800, 800);
    cHighLowSiPM->Divide(3, 3);

    TCanvas *cHighLowSiPM_Fitting = new TCanvas(("Graph_HighLowSiPM_Fitting_" + RUN + add).c_str(), ("Graph_HighLowSiPM_Fitting_" + RUN + add).c_str(), 800, 800);
    cHighLowSiPM_Fitting->Divide(3, 3);

    TCanvas *cSiPM_GainH = new TCanvas(("2D_SiPM_Gain_High_Fitting" + RUN + add).c_str(), ("2D_SiPM_Gain_High_Fitting" + RUN + add).c_str(), 800, 800);
    cSiPM_GainH->Divide(3, 3);

    TCanvas *cSiPM_Gain_FittingH = new TCanvas(("Graph_SiPM_Gain_High_Fitting_" + RUN + add).c_str(), ("Graph_SiPM_Gain_High_Fitting_" + RUN + add).c_str(), 800, 800);
    cSiPM_Gain_FittingH->Divide(3, 3);

    TCanvas *cSiPM_GainL = new TCanvas(("2D_SiPM_Gain_Low_Fitting" + RUN + add).c_str(), ("2D_SiPM_Gain_Low_Fitting" + RUN + add).c_str(), 800, 800);
    cSiPM_GainL->Divide(3, 3);

    TCanvas *cSiPM_Gain_FittingL = new TCanvas(("Graph_SiPM_Gain_Low_Fitting_" + RUN + add).c_str(), ("Graph_SiPM_Gain_Low_Fitting_" + RUN + add).c_str(), 800, 800);
    cSiPM_Gain_FittingL->Divide(3, 3);

    TCanvas *c_High = new TCanvas(("All_High" + RUN + add).c_str(), ("All_High" + RUN + add).c_str(), 800, 800);
    TCanvas *c_Low = new TCanvas(("All_Low" + RUN + add).c_str(), ("All_Low" + RUN + add).c_str(), 800, 800);

    for (int i = 0; i < SIGNAL_MAX; i++)
    {

        if (IsDetectorBetaHigh(i))
        {
            if (Verbose == 1)
                Info(detectorName[i], 1);

            dir_detector[RUN][i]->cd();
            H_SiPM_High[RUN][i]->Write();
            H_SiPM_HighLow[RUN][i]->Write();
            // G_SiPM_HighLow[RUN][i]->Write();

            c_High->cd();
            H_SiPM_High[RUN][i]->Draw("HIST SAME");

            if (!correction)
            {
                MatchingLowHigh[RUN][i]->SetName(("MatchingLowHigh_" + RUN + "_" + detectorName[i]).c_str());
                MatchingLowHigh[RUN][i]->Write();
                MatchingSiPM[RUN][i]->SetName(("MatchingSiPM_" + RUN + "_" + detectorName[i]).c_str());
                MatchingSiPM[RUN][i]->Write();
            }

            /////////// ##### FIITING HIGH LOW ###### ///////////
            // TH2D and FIT
            cHighLowSiPM->cd(GetDetectorChannel(i));
            H_SiPM_HighLow[RUN][i]->Draw("COLZ");
            TLatex *text; 
            if (!correction)
            {
                MatchingLowHigh[RUN][i]->SetLineColor(kRed);
                MatchingLowHigh[RUN][i]->Draw("SAME");
                text = new TLatex();
                text->SetNDC();
                text->SetTextSize(0.1);
                std::ostringstream streamObj;
                streamObj << std::fixed << std::setprecision(2) << MatchingLowHigh[RUN][i]->GetParameter(0);
                std::string parameterStr = streamObj.str();
                text->DrawLatex(0.7, 0.8, ("#color[2]{" + parameterStr + "}").c_str());
                text->Draw("SAME");
            }

            // TGraph and FIT
            cHighLowSiPM_Fitting->cd(GetDetectorChannel(i));
            G_SiPM_HighLow[RUN][i]->SetTitle(("SiPM " + to_string(GetDetectorChannel(i))).c_str());
            G_SiPM_HighLow[RUN][i]->GetXaxis()->SetTitle("Low SiPM");
            G_SiPM_HighLow[RUN][i]->GetYaxis()->SetTitle("High SiPM");
            G_SiPM_HighLow[RUN][i]->GetXaxis()->CenterTitle();
            G_SiPM_HighLow[RUN][i]->GetYaxis()->CenterTitle();
            // G_SiPM_HighLow[RUN][i]->Draw("AP");
            if (!correction)
            {
                MatchingLowHigh[RUN][i]->SetLineColor(kRed);
                MatchingLowHigh[RUN][i]->Draw("SAME");

                text->Draw("SAME");
            }

            /////////// ##### FIITING SiPM Gain ###### ///////////
            // TH2D and FIT
            cSiPM_GainH->cd(GetDetectorChannel(i));
            H_SiPM_Gain[RUN][i]->Draw("COLZ");
            TLatex *textGain;
            if (!correction)
            {
                MatchingSiPM[RUN][i]->SetLineColor(kRed);
                MatchingSiPM[RUN][i]->Draw("SAME");
                textGain = new TLatex();
                textGain->SetNDC();
                textGain->SetTextSize(0.1);
                std::ostringstream streamObjGain;
                streamObjGain << std::fixed << std::setprecision(2) << MatchingSiPM[RUN][i]->GetParameter(0);
                std::string parameterStrGain = streamObjGain.str();
                textGain->DrawLatex(0.7, 0.8, ("#color[2]{" + parameterStrGain + "}").c_str());
                textGain->Draw("SAME");
            }
            // TGraph and FIT
            cSiPM_Gain_FittingH->cd(GetDetectorChannel(i));
            G_SiPM_Gain[RUN][i]->SetTitle(("SiPM " + to_string(GetDetectorChannel(i))).c_str());
            G_SiPM_Gain[RUN][i]->GetXaxis()->SetTitle("Low SiPM");
            G_SiPM_Gain[RUN][i]->GetYaxis()->SetTitle("High SiPM");
            G_SiPM_Gain[RUN][i]->GetXaxis()->CenterTitle();
            G_SiPM_Gain[RUN][i]->GetYaxis()->CenterTitle();
            G_SiPM_Gain[RUN][i]->Draw("AP");
            if (!correction)
            {
                MatchingSiPM[RUN][i]->SetLineColor(kRed);
                MatchingSiPM[RUN][i]->Draw("SAME");

                textGain->Draw("SAME");
            }
        }
        if (IsDetectorBetaLow(i))
        {
            
            if (Verbose == 1)
                Info(detectorName[i], 1);
            dir_detector[RUN][i - 10]->cd();
            H_SiPM_Low[RUN][i]->Write();

            c_Low->cd();
            H_SiPM_Low[RUN][i]->Draw("HIST SAME");

            if (!correction)
            {
                MatchingSiPM[RUN][i]->SetName(("MatchingSiPM_" + RUN + "_" + detectorName[i]).c_str());
                MatchingSiPM[RUN][i]->Write();
            }

            /////////// ##### FIITING SiPM Gain ###### ///////////
            // TH2D and FIT
            cSiPM_GainL->cd(GetDetectorChannel(i));
            H_SiPM_Gain[RUN][i]->Draw("COLZ");
            TLatex *textGain;
            if (!correction)
            {
                MatchingSiPM[RUN][i]->SetLineColor(kRed);
                MatchingSiPM[RUN][i]->Draw("SAME");
                textGain = new TLatex();
                textGain->SetNDC();
                textGain->SetTextSize(0.1);
                std::ostringstream streamObjGain;
                streamObjGain << std::fixed << std::setprecision(2) << MatchingSiPM[RUN][i]->GetParameter(0);
                std::string parameterStrGain = streamObjGain.str();
                textGain->DrawLatex(0.7, 0.8, ("#color[2]{" + parameterStrGain + "}").c_str());
                textGain->Draw("SAME");
            }

            // // TGraph and FIT
            cSiPM_Gain_FittingH->cd(GetDetectorChannel(i));
            G_SiPM_Gain[RUN][i]->SetTitle(("SiPM " + to_string(GetDetectorChannel(i))).c_str());
            G_SiPM_Gain[RUN][i]->GetXaxis()->SetTitle("Low SiPM");
            G_SiPM_Gain[RUN][i]->GetYaxis()->SetTitle("High SiPM");
            G_SiPM_Gain[RUN][i]->GetXaxis()->CenterTitle();
            G_SiPM_Gain[RUN][i]->GetYaxis()->CenterTitle();
            // G_SiPM_Gain[RUN][i]->Draw("AP");
            if (!correction)
            {
                MatchingSiPM[RUN][i]->SetLineColor(kRed);
                MatchingSiPM[RUN][i]->Draw("SAME");
                textGain->Draw("SAME");
            }
        }
    }

    dir[RUN]->cd();
    cHighLowSiPM->Write();
    cHighLowSiPM_Fitting->Write();
    cSiPM_GainH->Write();
    cSiPM_Gain_FittingH->Write();
    cSiPM_GainL->Write();
    cSiPM_Gain_FittingL->Write();
    c_High->Write();
    c_Low->Write();
    MATCHED_File->Close();
    delete cHighLowSiPM;
    delete cHighLowSiPM_Fitting;
    delete cSiPM_GainH;
    delete cSiPM_Gain_FittingH;
    delete cSiPM_GainL;
    delete cSiPM_Gain_FittingL;
    delete c_High;
    delete c_Low;
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBetaHigh(i))
        {
            delete G_SiPM_HighLow[RUN][i];
            delete G_SiPM_Gain[RUN][i];
            delete H_SiPM_High[RUN][i];
            delete H_SiPM_HighLow[RUN][i];
            delete H_SiPM_Gain[RUN][i];
        }
        if (IsDetectorBetaLow(i))
        {
            delete H_SiPM_Low[RUN][i];
            delete G_SiPM_Gain[RUN][i];
            delete H_SiPM_Gain[RUN][i];
        }

    }
    MATCHED_File = MyTFile(DIR_ROOT_DATA_MATCHED + MATCHED_Filename, "UPDATE", "Q");
    Info("Write Histograms done");
}

void WriteHistogramAfterCorrection(int Verbose = 0)
{
    bool correction = true;
    string add = "";
    if (correction) add = "_Matched";

    gStyle->SetOptStat(0);
    MATCHED_File->cd();
    Info("Write Histograms for run " + RUN + add);

    MATCHED_File->cd(RUN.c_str());
    TCanvas *cHighLowSiPM = new TCanvas(("2D_HighLowSiPM_Fitting" + RUN + add).c_str(), ("2D_HighLowSiPM_Fitting" + RUN + add).c_str(), 800, 800);
    cHighLowSiPM->Divide(3, 3);

    TCanvas *cHighLowSiPM_Fitting = new TCanvas(("Graph_HighLowSiPM_Fitting_" + RUN + add).c_str(), ("Graph_HighLowSiPM_Fitting_" + RUN + add).c_str(), 800, 800);
    cHighLowSiPM_Fitting->Divide(3, 3);

    TCanvas *cSiPM_GainH = new TCanvas(("2D_SiPM_Gain_High_Fitting" + RUN + add).c_str(), ("2D_SiPM_Gain_High_Fitting" + RUN + add).c_str(), 800, 800);
    cSiPM_GainH->Divide(3, 3);

    TCanvas *cSiPM_Gain_FittingH = new TCanvas(("Graph_SiPM_Gain_High_Fitting_" + RUN + add).c_str(), ("Graph_SiPM_Gain_High_Fitting_" + RUN + add).c_str(), 800, 800);
    cSiPM_Gain_FittingH->Divide(3, 3);

    TCanvas *cSiPM_GainL = new TCanvas(("2D_SiPM_Gain_Low_Fitting" + RUN + add).c_str(), ("2D_SiPM_Gain_Low_Fitting" + RUN + add).c_str(), 800, 800);
    cSiPM_GainL->Divide(3, 3);

    TCanvas *cSiPM_Gain_FittingL = new TCanvas(("Graph_SiPM_Gain_Low_Fitting_" + RUN + add).c_str(), ("Graph_SiPM_Gain_Low_Fitting_" + RUN + add).c_str(), 800, 800);
    cSiPM_Gain_FittingL->Divide(3, 3);

    TCanvas *c_High = new TCanvas(("All_High" + RUN + add).c_str(), ("All_High" + RUN + add).c_str(), 800, 800);
    TCanvas *c_Low = new TCanvas(("All_Low" + RUN + add).c_str(), ("All_Low" + RUN + add).c_str(), 800, 800);
    
    TLine *linediag = new TLine(0, 0, 100e6, 100e6);
    linediag->SetLineColor(kRed);
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBetaHigh(i))
        {
            if (Verbose == 1)
                Info(detectorName[i], 1);

            MATCHED_File->cd((RUN + "/SiPM_" + to_string(GetDetectorChannel(i))).c_str());
            H_SiPM_High[RUN][i]->Write();
            H_SiPM_HighLow[RUN][i]->Write();
            // G_SiPM_HighLow[RUN][i]->Write();

            c_High->cd();
            H_SiPM_High[RUN][i]->SetLineColor(GetDetectorChannel(i));
            H_SiPM_High[RUN][i]->Draw("HIST SAME");

            /////////// ##### FIITING HIGH LOW ###### ///////////
            // TH2D and FIT
            cHighLowSiPM->cd(GetDetectorChannel(i));
            H_SiPM_HighLow[RUN][i]->Draw("COLZ");
            linediag->Draw("SAME");
            // TGraph and FIT
            cHighLowSiPM_Fitting->cd(GetDetectorChannel(i));
            G_SiPM_HighLow[RUN][i]->SetTitle(("SiPM " + to_string(GetDetectorChannel(i))).c_str());
            G_SiPM_HighLow[RUN][i]->GetXaxis()->SetTitle("Low SiPM");
            G_SiPM_HighLow[RUN][i]->GetYaxis()->SetTitle("High SiPM");
            G_SiPM_HighLow[RUN][i]->GetXaxis()->CenterTitle();
            G_SiPM_HighLow[RUN][i]->GetYaxis()->CenterTitle();
            // G_SiPM_HighLow[RUN][i]->Draw("AP");

            /////////// ##### FIITING SiPM Gain ###### ///////////
            // TH2D and FIT
            cSiPM_GainH->cd(GetDetectorChannel(i));
            H_SiPM_Gain[RUN][i]->Draw("COLZ");
            linediag->Draw("SAME");
            // TGraph and FIT
            cSiPM_Gain_FittingH->cd(GetDetectorChannel(i));
            G_SiPM_Gain[RUN][i]->SetTitle(("SiPM " + to_string(GetDetectorChannel(i))).c_str());
            G_SiPM_Gain[RUN][i]->GetXaxis()->SetTitle("Low SiPM");
            G_SiPM_Gain[RUN][i]->GetYaxis()->SetTitle("High SiPM");
            G_SiPM_Gain[RUN][i]->GetXaxis()->CenterTitle();
            G_SiPM_Gain[RUN][i]->GetYaxis()->CenterTitle();
            // G_SiPM_Gain[RUN][i]->Draw("AP");
        }
        if (IsDetectorBetaLow(i))
        {
            
            if (Verbose == 1)
                Info(detectorName[i], 1);
            MATCHED_File->cd((RUN + "/SiPM_" + to_string(GetDetectorChannel(i))).c_str());
            H_SiPM_Low[RUN][i]->Write();

            c_Low->cd();
            H_SiPM_Low[RUN][i]->SetLineColor(GetDetectorChannel(i));
            H_SiPM_Low[RUN][i]->Draw("HIST SAME");

            /////////// ##### FIITING SiPM Gain ###### ///////////
            // TH2D and FIT
            cSiPM_GainL->cd(GetDetectorChannel(i));
            H_SiPM_Gain[RUN][i]->Draw("COLZ");
            linediag->Draw("SAME");

            // // TGraph and FIT
            cSiPM_Gain_FittingH->cd(GetDetectorChannel(i));
            G_SiPM_Gain[RUN][i]->SetTitle(("SiPM " + to_string(GetDetectorChannel(i))).c_str());
            G_SiPM_Gain[RUN][i]->GetXaxis()->SetTitle("Low SiPM");
            G_SiPM_Gain[RUN][i]->GetYaxis()->SetTitle("High SiPM");
            G_SiPM_Gain[RUN][i]->GetXaxis()->CenterTitle();
            G_SiPM_Gain[RUN][i]->GetYaxis()->CenterTitle();
            // G_SiPM_Gain[RUN][i]->Draw("AP");
        }
    }

    MATCHED_File->cd(RUN.c_str());
    cHighLowSiPM->Write();
    cHighLowSiPM_Fitting->Write();
    cSiPM_GainH->Write();
    cSiPM_Gain_FittingH->Write();
    cSiPM_GainL->Write();
    cSiPM_Gain_FittingL->Write();
    c_High->Write();
    c_Low->Write();
    MATCHED_File->Close();
    delete cHighLowSiPM;
    delete cHighLowSiPM_Fitting;
    delete cSiPM_GainH;
    delete cSiPM_Gain_FittingH;
    delete cSiPM_GainL;
    delete cSiPM_Gain_FittingL;
    delete c_High;
    delete c_Low;
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBetaHigh(i))
        {
            delete G_SiPM_HighLow[RUN][i];
            delete G_SiPM_Gain[RUN][i];
            delete H_SiPM_High[RUN][i];
            delete H_SiPM_HighLow[RUN][i];
            delete H_SiPM_Gain[RUN][i];
        }
        if (IsDetectorBetaLow(i))
        {
            delete H_SiPM_Low[RUN][i];
            delete G_SiPM_Gain[RUN][i];
            delete H_SiPM_Gain[RUN][i];
        }

    }
    MATCHED_File = MyTFile(DIR_ROOT_DATA_MATCHED + "SiPM_Matching_test.root", "UPDATE", "Q");
    Info("Write Histograms done");
}

void WriteHistograms(int Verbose = 0)
{
    gStyle->SetOptStat(0);
    MATCHED_File->cd();
    Info("Write Histograms over all runs");

    int counterHighLow[SIGNAL_MAX] = {0};
    int counterSiPM_H[SIGNAL_MAX] = {0};
    int counterSiPM_L[SIGNAL_MAX] = {0};

    for (auto &pairr : Map_RunFiles)
    {
        TCanvas *cNucleus_SiPM_High[SIGNAL_MAX];
        TCanvas *cNucleus_SiPM_Low[SIGNAL_MAX];

        for (int i = 0; i < SIGNAL_MAX; i++)
            {
                if (IsDetectorBetaHigh(i))
                {
                    cNucleus_SiPM_High[i] = new TCanvas((detectorName[i] + "_" + pairr.first).c_str(), (detectorName[i] + "_" + pairr.first).c_str(), 1920, 1080);

                }
                if (IsDetectorBetaLow(i))
                {
                    cNucleus_SiPM_Low[i] = new TCanvas((detectorName[i] + "_" + pairr.first).c_str(), (detectorName[i] + "_" + pairr.first).c_str(), 1920, 1080);
                }
            }       

        for (string run : pairr.second)
        {
            // string filename = SearchFiles(DIR_ROOT_DATA_GROUPED, run);
            // GROUPED_File[run] = MyTFile(DIR_ROOT_DATA_GROUPED + filename, "READ", "Q");
            // if (GROUPED_File[run] == NULL)
            //     continue;
            // else
            //     GROUPED_File[run]->Close();

            MATCHED_File->cd();
            if (Verbose == 1)
                Info(run, 1);

            for (int i = 0; i < SIGNAL_MAX; i++)
            {
                if (IsDetectorBetaHigh(i))
                {
                    if (Verbose == 1)
                        Info(detectorName[i], 2);
                    
                    // All HIgh/Low factor on the same canvas for all runs
                    MatchingLowHigh[run][i] = (TF1 *)MATCHED_File->Get((run + "/SiPM_" + to_string(GetDetectorChannel(i)) + "/MatchingLowHigh_" + run + "_" + detectorName[i]).c_str());
                    gHighLow_RUNS_a[i]->AddPoint(stoi(run), MatchingLowHigh[run][i]->GetParameter(0));
                    gHighLow_RUNS_a[i]->SetPointError(gHighLow_RUNS_a[i]->GetN()-1, 0, MatchingLowHigh[run][i]->GetParError(0));
                    gHighLow_RUNS_b[i]->AddPoint(stoi(run), MatchingLowHigh[run][i]->GetParameter(1));
                    gHighLow_RUNS_b[i]->SetPointError(gHighLow_RUNS_b[i]->GetN()-1, 0, MatchingLowHigh[run][i]->GetParError(1));

                    // All SiPM factor on the same canvas for all runs for High
                    MatchingSiPM[run][i] = (TF1 *)MATCHED_File->Get((run + "/SiPM_" + to_string(GetDetectorChannel(i)) + "/MatchingSiPM_" + run + "_" + detectorName[i]).c_str());
                    gSiPM_High_RUNS_a[i]->AddPoint(stoi(run), MatchingSiPM[run][i]->GetParameter(0));
                    gSiPM_High_RUNS_a[i]->SetPointError(gSiPM_High_RUNS_a[i]->GetN()-1, 0, MatchingSiPM[run][i]->GetParError(0));
                    gSiPM_High_RUNS_b[i]->AddPoint(stoi(run), MatchingSiPM[run][i]->GetParameter(1));
                    gSiPM_High_RUNS_b[i]->SetPointError(gSiPM_High_RUNS_b[i]->GetN()-1, 0, MatchingSiPM[run][i]->GetParError(1));

                    
                    H_SiPM_High[run][i] = (TH1D*)MATCHED_File->Get((run + "/SiPM_" + to_string(GetDetectorChannel(i)) + "/H_SiPM_High_" + run + "_" + detectorName[i]).c_str());
                    H_SiPM_High[run][i]->Scale(1. / H_SiPM_High[run][i]->Integral());
                    cNucleus_SiPM_High[i]->cd();
                    if (stoi(run) == 77)
                    {
                        H_SiPM_High[run][i]->SetLineColor(kRed);
                        H_SiPM_High[run][i]->Draw("HIST SAME");
                    }
                    else
                    {
                        H_SiPM_High[run][i]->SetLineColor(stoi(run));
                        H_SiPM_High[run][i]->Draw("HIST SAME");
                    }
                    
                }

                else if (IsDetectorBetaLow(i))
                {
                    if (Verbose == 1)
                        Info(detectorName[i], 2);

                    // All SiPM factor on the same canvas for all runs for High
                    MatchingSiPM[run][i] = (TF1 *)MATCHED_File->Get((run + "/SiPM_" + to_string(GetDetectorChannel(i)) + "/MatchingSiPM_" + run + "_" + detectorName[i]).c_str());
                    gSiPM_Low_RUNS_a[i]->AddPoint(stoi(run), MatchingSiPM[run][i]->GetParameter(0));
                    gSiPM_Low_RUNS_a[i]->SetPointError(gSiPM_Low_RUNS_a[i]->GetN()-1, 0, MatchingSiPM[run][i]->GetParError(0));
                    gSiPM_Low_RUNS_b[i]->AddPoint(stoi(run), MatchingSiPM[run][i]->GetParameter(1));
                    gSiPM_Low_RUNS_b[i]->SetPointError(gSiPM_Low_RUNS_b[i]->GetN()-1, 0, MatchingSiPM[run][i]->GetParError(1));

                    cNucleus_SiPM_Low[i]->cd();
                    H_SiPM_Low[run][i] = (TH1D*)MATCHED_File->Get((run + "/SiPM_" + to_string(GetDetectorChannel(i)) + "/H_SiPM_Low_" + run + "_" + detectorName[i]).c_str());
                    H_SiPM_Low[run][i]->Scale(1. / H_SiPM_Low[run][i]->Integral());
                    if (stoi(run) == 77)
                    {
                        H_SiPM_Low[run][i]->SetLineColor(kRed);
                        H_SiPM_Low[run][i]->Draw("SAME");
                    }
                    else
                    {
                        H_SiPM_Low[run][i]->SetLineColor(stoi(run));
                        H_SiPM_Low[run][i]->Draw("SAME");
                    }
                }
            }
        }

        MATCHED_File->cd();
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorBetaHigh(i))
            {
                cNucleus_SiPM_High[i]->Write();

            }
            if (IsDetectorBetaLow(i))
            {
                cNucleus_SiPM_Low[i]->Write();
            }
        }
    }

    MATCHED_File->cd();
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBetaHigh(i))
        {
            if (Verbose == 1)
                Info(detectorName[i], 1);

            if (Verbose == 1)
                Info("Fitting High / Low", 1);
            TCanvas *c = new TCanvas(("Matching_HighLow" + to_string(GetDetectorChannel(i))).c_str(), ("Matching_HighLow" + to_string(GetDetectorChannel(i))).c_str(), 1920, 1080);
            c->Divide(1, 2);
            c->cd(1);
            // seperating wrong value from the fit and the TGraphErrors
            TGraphErrors *gfit_a = new TGraphErrors();
            gfit_a->SetMarkerStyle(20);
            gfit_a->SetMarkerColor(kBlack);
            gfit_a->SetTitle(("Matching " + detectorName[i]).c_str());
            gfit_a->GetXaxis()->SetTitle("Run");
            gfit_a->GetYaxis()->SetTitle("Gain factor Low/High");
            TGraphErrors *gnofit_a = new TGraphErrors();
            gnofit_a->SetMarkerStyle(20);
            gnofit_a->SetMarkerColor(kGray);
            gnofit_a->SetTitle(("Matching NotWorking" + detectorName[i]).c_str());
            gnofit_a->GetXaxis()->SetTitle("Run");
            gnofit_a->GetYaxis()->SetTitle("Gain factor Low/High");

            for (int j = 0; j < gHighLow_RUNS_a[i]->GetN(); j++)
            {
                double x, y;
                gHighLow_RUNS_a[i]->GetPoint(j, x, y);
                double ey = gHighLow_RUNS_a[i]->GetErrorY(j);
                if (ey == 0)
                {
                    gnofit_a->AddPoint(x, y);
                    gnofit_a->SetPointError(gnofit_a->GetN() - 1, 0, ey);
                }
                else
                {
                    gfit_a->AddPoint(x, y);
                    gfit_a->SetPointError(gfit_a->GetN() - 1, 0, ey);
                }
            }
            
            //draw graphs
            TMultiGraph *mg_a = new TMultiGraph();    
            mg_a->Add(gfit_a, "AP");
            mg_a->Add(gnofit_a, "AP");
            mg_a->Draw("AP");
            //draw fit
            gfit_a->Fit("pol0", "Q");
            gfit_a->GetFunction("pol0")->SetLineColor(kRed);
            gfit_a->GetFunction("pol0")->Draw("SAME");
            gHighLow_RUNS_a[i] = gfit_a;

            c->cd(2);
            // seperating wrong value from the fit and the TGraphErrors
            TGraphErrors *gfit_b = new TGraphErrors();
            gfit_b->SetMarkerStyle(20);
            gfit_b->SetMarkerColor(kBlack);
            gfit_b->SetTitle(("Matching " + detectorName[i]).c_str());
            gfit_b->GetXaxis()->SetTitle("Run");
            gfit_b->GetYaxis()->SetTitle("Offset factor Low/High");
            TGraphErrors *gnofit_b = new TGraphErrors();
            gnofit_b->SetMarkerStyle(20);
            gnofit_b->SetMarkerColor(kGray);
            gnofit_b->SetTitle(("Matching NotWorking" + detectorName[i]).c_str());
            gnofit_b->GetXaxis()->SetTitle("Run");
            gnofit_b->GetYaxis()->SetTitle("Offset factor Low/High");
            for (int j = 0; j < gHighLow_RUNS_b[i]->GetN(); j++)
            {
                double x, y;
                gHighLow_RUNS_b[i]->GetPoint(j, x, y);
                double ey = gHighLow_RUNS_b[i]->GetErrorY(j);
                if (ey == 0)
                {
                    gnofit_b->AddPoint(x, y);
                    gnofit_b->SetPointError(gnofit_b->GetN() - 1, 0, ey);
                }
                else
                {
                    gfit_b->AddPoint(x, y);
                    gfit_b->SetPointError(gfit_b->GetN() - 1, 0, ey);
                }
            }
            //draw graphs
            TMultiGraph *mg_b = new TMultiGraph();
            mg_b->Add(gfit_b, "AP");
            mg_b->Add(gnofit_b, "AP");
            mg_b->Draw("AP");
            //draw fit
            if (gfit_b->GetN() > 1)
            {
                gfit_b->Fit("pol0", "Q");
                gfit_b->GetFunction("pol0");
                gfit_b->GetFunction("pol0")->SetLineColor(kRed);
                gfit_b->GetFunction("pol0")->Draw("SAME");
            }
            
            c->Write();
            
            if (Verbose == 1)
                Info("Fitting SiPM", 1);
            TCanvas *cSiPM = new TCanvas(("Matching_SiPM_High" + to_string(GetDetectorChannel(i))).c_str(), ("Matching_SiPM_High" + to_string(GetDetectorChannel(i))).c_str(), 1920, 1080);
            cSiPM->Divide(1, 2);
            cSiPM->cd(1);
            // seperating wrong value from the fit and the TGraphErrors
            TGraphErrors *gSiPM_fit_a = new TGraphErrors();
            gSiPM_fit_a->SetMarkerStyle(20);
            gSiPM_fit_a->SetMarkerColor(kBlack);
            gSiPM_fit_a->SetTitle(("Matching SiPM with" + detectorName[i]).c_str());
            gSiPM_fit_a->GetXaxis()->SetTitle("Run");
            gSiPM_fit_a->GetYaxis()->SetTitle("Gain factor SiPM");
            TGraphErrors *gSiPM_nofit_a = new TGraphErrors();
            gSiPM_nofit_a->SetMarkerStyle(20);
            gSiPM_nofit_a->SetMarkerColor(kGray);
            gSiPM_nofit_a->SetTitle(("Matching SiPM NotWorking" + detectorName[i]).c_str());
            gSiPM_nofit_a->GetXaxis()->SetTitle("Run");
            gSiPM_nofit_a->GetYaxis()->SetTitle("Gain factor SiPM");
            for (int j = 0; j < gSiPM_High_RUNS_a[i]->GetN(); j++)
            {
                double x, y;
                gSiPM_High_RUNS_a[i]->GetPoint(j, x, y);
                double ey = gSiPM_High_RUNS_a[i]->GetErrorY(j);
                if (ey == 0)
                {
                    gSiPM_nofit_a->AddPoint(x, y);
                    gSiPM_nofit_a->SetPointError(gSiPM_nofit_a->GetN() - 1, 0, ey);
                }
                else
                {
                    gSiPM_fit_a->AddPoint(x, y);
                    gSiPM_fit_a->SetPointError(gSiPM_fit_a->GetN() - 1, 0, ey);
                }
            }
            gSiPM_High_RUNS_a[i] = gSiPM_fit_a;
            //draw graphs
            TMultiGraph *mgSiPM = new TMultiGraph();
            mgSiPM->Add(gSiPM_fit_a, "AP");
            mgSiPM->Add(gSiPM_nofit_a, "AP");
            mgSiPM->Draw("AP");
            //draw fit
            gSiPM_fit_a->Fit("pol0", "Q");
            gSiPM_fit_a->GetFunction("pol0");
            gSiPM_fit_a->GetFunction("pol0")->SetLineColor(kRed);
            gSiPM_fit_a->GetFunction("pol0")->Draw("SAME");

            cSiPM->cd(2);
            // seperating wrong value from the fit and the TGraphErrors
            TGraphErrors *gSiPM_fit_b = new TGraphErrors();
            gSiPM_fit_b->SetMarkerStyle(20);
            gSiPM_fit_b->SetMarkerColor(kBlack);
            gSiPM_fit_b->SetTitle(("Matching SiPM with" + detectorName[i]).c_str());
            gSiPM_fit_b->GetXaxis()->SetTitle("Run");
            gSiPM_fit_b->GetYaxis()->SetTitle("Offset factor SiPM");
            TGraphErrors *gSiPM_nofit_b = new TGraphErrors();
            gSiPM_nofit_b->SetMarkerStyle(20);
            gSiPM_nofit_b->SetMarkerColor(kGray);
            gSiPM_nofit_b->SetTitle(("Matching SiPM NotWorking" + detectorName[i]).c_str());
            gSiPM_nofit_b->GetXaxis()->SetTitle("Run");
            gSiPM_nofit_b->GetYaxis()->SetTitle("Offset factor SiPM");
            for (int j = 0; j < gSiPM_High_RUNS_b[i]->GetN(); j++)
            {
                double x, y;
                gSiPM_High_RUNS_b[i]->GetPoint(j, x, y);
                double ey = gSiPM_High_RUNS_b[i]->GetErrorY(j);
                if (ey == 0)
                {
                    gSiPM_nofit_b->AddPoint(x, y);
                    gSiPM_nofit_b->SetPointError(gSiPM_nofit_b->GetN() - 1, 0, ey);
                }
                else
                {
                    gSiPM_fit_b->AddPoint(x, y);
                    gSiPM_fit_b->SetPointError(gSiPM_fit_b->GetN() - 1, 0, ey);
                }
            }
            //draw graphs
            TMultiGraph *mgSiPM_b = new TMultiGraph();
            mgSiPM_b->Add(gSiPM_fit_b, "AP");
            mgSiPM_b->Add(gSiPM_nofit_b, "AP");
            mgSiPM_b->Draw("AP");
            //draw fit
            if (gSiPM_fit_b->GetN() > 1)
            {
                gSiPM_fit_b->Fit("pol0", "Q");
                gSiPM_fit_b->GetFunction("pol0");
                gSiPM_fit_b->GetFunction("pol0")->SetLineColor(kRed);
                gSiPM_fit_b->GetFunction("pol0")->Draw("SAME");     
            }       

    
            cSiPM->Write();
        }

        if (IsDetectorBetaLow(i))
        {
            TCanvas *cSiPM = new TCanvas(("Matching_SiPM_Low" + to_string(GetDetectorChannel(i))).c_str(), ("Matching_SiPM_Low" + to_string(GetDetectorChannel(i))).c_str(), 1920, 1080);
            cSiPM->Divide(1, 2);
            cSiPM->cd(1);
            // seperating wrong value from the fit and the TGraphErrors
            TGraphErrors *gSiPM_fit_a = new TGraphErrors();
            gSiPM_fit_a->SetMarkerStyle(20);
            gSiPM_fit_a->SetMarkerColor(kBlack);
            gSiPM_fit_a->SetTitle(("Matching SiPM with" + detectorName[i]).c_str());
            gSiPM_fit_a->GetXaxis()->SetTitle("Run");
            gSiPM_fit_a->GetYaxis()->SetTitle("Gain factor SiPM");
            TGraphErrors *gSiPM_nofit_a = new TGraphErrors();
            gSiPM_nofit_a->SetMarkerStyle(20);
            gSiPM_nofit_a->SetMarkerColor(kGray);
            gSiPM_nofit_a->SetTitle(("Matching SiPM NotWorking" + detectorName[i]).c_str());
            gSiPM_nofit_a->GetXaxis()->SetTitle("Run");
            gSiPM_nofit_a->GetYaxis()->SetTitle("Gain factor SiPM");
            for (int j = 0; j < gSiPM_Low_RUNS_a[i]->GetN(); j++)
            {
                double x, y;
                gSiPM_Low_RUNS_a[i]->GetPoint(j, x, y);
                double ey = gSiPM_Low_RUNS_a[i]->GetErrorY(j);
                if (ey == 0)
                {
                    gSiPM_nofit_a->AddPoint(x, y);
                    gSiPM_nofit_a->SetPointError(gSiPM_nofit_a->GetN() - 1, 0, ey);
                }
                else
                {
                    gSiPM_fit_a->AddPoint(x, y);
                    gSiPM_fit_a->SetPointError(gSiPM_fit_a->GetN() - 1, 0, ey);
                }
            }
            //draw graphs
            TMultiGraph *mgSiPM = new TMultiGraph();
            mgSiPM->Add(gSiPM_fit_a, "AP");
            mgSiPM->Add(gSiPM_nofit_a, "AP");
            mgSiPM->Draw("AP");
            //draw fit
            gSiPM_fit_a->Fit("pol0", "Q");
            gSiPM_fit_a->GetFunction("pol0");
            gSiPM_fit_a->GetFunction("pol0")->SetLineColor(kRed);
            gSiPM_fit_a->GetFunction("pol0")->Draw("SAME");

            cSiPM->cd(2);
            // seperating wrong value from the fit and the TGraphErrors
            TGraphErrors *gSiPM_fit_b = new TGraphErrors();
            gSiPM_fit_b->SetMarkerStyle(20);
            gSiPM_fit_b->SetMarkerColor(kBlack);
            gSiPM_fit_b->SetTitle(("Matching SiPM with" + detectorName[i]).c_str());
            gSiPM_fit_b->GetXaxis()->SetTitle("Run");
            gSiPM_fit_b->GetYaxis()->SetTitle("Offset factor SiPM");
            TGraphErrors *gSiPM_nofit_b = new TGraphErrors();
            gSiPM_nofit_b->SetMarkerStyle(20);
            gSiPM_nofit_b->SetMarkerColor(kGray);
            gSiPM_nofit_b->SetTitle(("Matching SiPM NotWorking" + detectorName[i]).c_str());
            gSiPM_nofit_b->GetXaxis()->SetTitle("Run");
            gSiPM_nofit_b->GetYaxis()->SetTitle("Offset factor SiPM");
            for (int j = 0; j < gSiPM_Low_RUNS_b[i]->GetN(); j++)
            {
                double x, y;
                gSiPM_Low_RUNS_b[i]->GetPoint(j, x, y);
                double ey = gSiPM_Low_RUNS_b[i]->GetErrorY(j);
                if (ey == 0)
                {
                    gSiPM_nofit_b->AddPoint(x, y);
                    gSiPM_nofit_b->SetPointError(gSiPM_nofit_b->GetN() - 1, 0, ey);
                }
                else
                {
                    gSiPM_fit_b->AddPoint(x, y);
                    gSiPM_fit_b->SetPointError(gSiPM_fit_b->GetN() - 1, 0, ey);
                }
            }
            //draw graphs
            TMultiGraph *mgSiPM_b = new TMultiGraph();
            mgSiPM_b->Add(gSiPM_fit_b, "AP");
            mgSiPM_b->Add(gSiPM_nofit_b, "AP");
            mgSiPM_b->Draw("AP");
            //draw fit
            if (gSiPM_fit_b->GetN() > 1)
            {
                gSiPM_fit_b->Fit("pol0", "Q");
                gSiPM_fit_b->GetFunction("pol0");
                gSiPM_fit_b->GetFunction("pol0")->SetLineColor(kRed);
                gSiPM_fit_b->GetFunction("pol0")->Draw("SAME");
            }
                          
            
            cSiPM->Write();
        }
    }

    // FINAL FILE
    //  functions used in this program to show the effect of correction
    //  will be used in Merger.c++
    TFile *file_matching = MyTFile(DIR_ROOT_DATA_MATCHED + "SiPM_Matching_Functions.root", "RECREATE");
    SiPM_ONOFF[SIGNAL_MAX];        // fitted on low used for both
    SiPM_HighLow[SIGNAL_MAX];      // fitted on both used for low
    SiPM_Gain[SIGNAL_MAX];         // fitted on high used for both
    file_matching->cd();
    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        // ONOFF --> a*x (correction for run before 80) (correction for both SiPM)
        if (IsDetectorBeta(det))
        {
            Info("Writing SiPM ONOFF for " + detectorName[det]);
            SiPM_ONOFF[det] = SiPM_ONOFF_Coef_RUN[110 + GetDetectorChannel(det)][stoi(Map_RunFiles["32Ar"][0])];
            SiPM_ONOFF[det]->SetName(("SiPM_ONOFF_" + detectorName[det]).c_str());
            SiPM_ONOFF[det]->Write();
        }

        // HighLow --> a*x (correction for low SiPM)
        if (IsDetectorBetaLow(det))
        {
            Info("Writing SiPM HighLow for " + detectorName[det]);
            SiPM_HighLow[det] = new TF1(("SiPM_HighLow_" + to_string(GetDetectorChannel(det))).c_str(), "[0]*x", 0, 100e6);
            SiPM_HighLow[det]->SetParameter(0, gHighLow_RUNS_a[100 + GetDetectorChannel(det)]->GetFunction("pol0")->GetParameter(0));
            SiPM_HighLow[det]->Write();
        }

        // // SiPM Gain --> a*x (correction for both SiPM)
        if (IsDetectorBeta(det))
        {
            Info("Writing SiPM Gain for " + detectorName[det]);
            SiPM_Gain[det] = new TF1(("SiPM_Gain_" + detectorName[det]).c_str(), "[0]*x", 0, 100e6);
            SiPM_Gain[det]->SetParameter(0, gSiPM_High_RUNS_a[100 + GetDetectorChannel(det)]->GetFunction("pol0")->GetParameter(0));
            SiPM_Gain[det]->Write();
        }
    }
    file_matching->Close();

    MATCHED_File->Close();
    MATCHED_File = MyTFile(DIR_ROOT_DATA_MATCHED + MATCHED_Filename, "UPDATE", "Q");
    Info("Write Histograms stop");
}


void WriteHistogramsAfterCorrection(int Verbose = 0)
{
    Info("Write Histograms after correction");
    TH1D* H_SiPM_NoCorrected[SIGNAL_MAX] = {NULL};
    TH1D* H_SiPM_Corrected[SIGNAL_MAX] = {NULL};


    // getting hist to merge 
    TFile *fread = MyTFile(DIR_ROOT_DATA_MATCHED + MATCHED_Filename, "READ");
    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorBetaHigh(det))
        {
            for (string run : Map_RunFiles["32Ar"])
            {
                if (H_SiPM_NoCorrected[det] == NULL)
                    H_SiPM_NoCorrected[det] = (TH1D*)fread->Get((run + "/SiPM_" + to_string(GetDetectorChannel(det)) + "/H_SiPM_High_" + run + "_" + detectorName[det]).c_str());
                else
                    H_SiPM_NoCorrected[det]->Add((TH1D*)fread->Get((run + "/SiPM_" + to_string(GetDetectorChannel(det)) + "/H_SiPM_High_" + run + "_" + detectorName[det]).c_str()));

                if (H_SiPM_Corrected[det] == NULL)
                    H_SiPM_Corrected[det] = (TH1D*)fread->Get((run + "/SiPM_" + to_string(GetDetectorChannel(det)) + "/H_SiPM_High_" + run + "_" + detectorName[det] + "_Matched").c_str());
                else
                    H_SiPM_Corrected[det]->Add((TH1D*)fread->Get((run + "/SiPM_" + to_string(GetDetectorChannel(det)) + "/H_SiPM_High_" + run + "_" + detectorName[det] + "_Matched").c_str()));
            }
        }

        if (IsDetectorBetaLow(det))
        {
            for (string run : Map_RunFiles["32Ar"])
            {
                if (H_SiPM_NoCorrected[det] == NULL)
                    H_SiPM_NoCorrected[det] = (TH1D*)fread->Get((run + "/SiPM_" + to_string(GetDetectorChannel(det)) + "/H_SiPM_Low_" + run + "_" + detectorName[det]).c_str());
                else
                    H_SiPM_NoCorrected[det]->Add((TH1D*)fread->Get((run + "/SiPM_" + to_string(GetDetectorChannel(det)) + "/H_SiPM_Low_" + run + "_" + detectorName[det]).c_str()));

                if (H_SiPM_Corrected[det] == NULL)
                    H_SiPM_Corrected[det] = (TH1D*)fread->Get((run + "/SiPM_" + to_string(GetDetectorChannel(det)) + "/H_SiPM_Low_" + run + "_" + detectorName[det] + "_Matched").c_str());
                else
                    H_SiPM_Corrected[det]->Add((TH1D*)fread->Get((run + "/SiPM_" + to_string(GetDetectorChannel(det)) + "/H_SiPM_Low_" + run + "_" + detectorName[det] + "_Matched").c_str()));
            }
        }
    }

    fread = MyTFile(DIR_ROOT_DATA_MATCHED + MATCHED_Filename, "UPDATE", "Q");

    TCanvas *cCompareCorrectionHigh = new TCanvas("cCompareCorrectionHigh", "cCompareCorrectionHigh", 800, 800);
    cCompareCorrectionHigh->Divide(3, 3);

    TCanvas *cCompareCorrectionLow = new TCanvas("cCompareCorrectionLow", "cCompareCorrectionLow", 800, 800);
    cCompareCorrectionLow->Divide(3, 3);

    TCanvas *cCompareCorrectionHigh_SiPM = new TCanvas("cCompareCorrectionHigh_SiPM", "cCompareCorrectionHigh_SiPM", 800, 800);
    cCompareCorrectionHigh_SiPM->Divide(1, 2);

    TCanvas *cCompareCorrectionLow_SiPM = new TCanvas("cCompareCorrectionLow_SiPM", "cCompareCorrectionLow_SiPM", 800, 800);
    cCompareCorrectionLow_SiPM->Divide(1, 2);


    double Highmin = 100e3;
    double Highmax = 500e3;
    double Lowmin = 500e3;
    double Lowmax = 2000e3;

    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorBetaHigh(det))
        {
            Info("Writing Corrected High SiPM for " + detectorName[det]);

            cCompareCorrectionHigh->cd(GetDetectorChannel(det));
            TH1D * H_SiPM_NoCorrected_temp = (TH1D*)H_SiPM_NoCorrected[det]->Clone();
            H_SiPM_NoCorrected_temp->GetXaxis()->SetRangeUser(Highmin, Highmax);
            H_SiPM_NoCorrected_temp->SetLineColor(kBlack);
            H_SiPM_NoCorrected_temp->Draw("HIST");
            TH1D * H_SiPM_Corrected_temp = (TH1D*)H_SiPM_Corrected[det]->Clone();
            H_SiPM_Corrected_temp->GetXaxis()->SetRangeUser(Highmin, Highmax);
            H_SiPM_Corrected_temp->Scale(H_SiPM_NoCorrected_temp->Integral() / H_SiPM_Corrected_temp->Integral());
            H_SiPM_Corrected_temp->SetLineColor(kRed);
            H_SiPM_Corrected_temp->Draw("HIST SAME");


            cCompareCorrectionHigh_SiPM->cd(1);
            H_SiPM_NoCorrected[det]->GetXaxis()->SetRangeUser(Highmin, Highmax);
            H_SiPM_NoCorrected[det]->Scale(H_SiPM_NoCorrected[101]->Integral() / H_SiPM_NoCorrected[det]->Integral());
            H_SiPM_NoCorrected[det]->SetLineColor(GetDetectorChannel(det));
            H_SiPM_NoCorrected[det]->Draw("HIST SAME");
            cCompareCorrectionHigh_SiPM->cd(2);
            H_SiPM_Corrected[det]->GetXaxis()->SetRangeUser(Highmin, Highmax);
            H_SiPM_Corrected[det]->Scale(H_SiPM_NoCorrected[101]->Integral() / H_SiPM_Corrected[det]->Integral());
            H_SiPM_Corrected[det]->SetLineColor(GetDetectorChannel(det));
            H_SiPM_Corrected[det]->Draw("HIST SAME");
        }
        if (IsDetectorBetaLow(det))
        {
            Info("Writing Corrected Low SiPM for " + detectorName[det]);
            cCompareCorrectionLow->cd(GetDetectorChannel(det));
            TH1D * H_SiPM_NoCorrected_temp = (TH1D*)H_SiPM_NoCorrected[det]->Clone();
            H_SiPM_NoCorrected_temp->GetXaxis()->SetRangeUser(Lowmin/10, Lowmax/10);
            H_SiPM_NoCorrected_temp->SetLineColor(kBlack);
            H_SiPM_NoCorrected_temp->Draw("HIST");
            TH1D * H_SiPM_Corrected_temp = (TH1D*)H_SiPM_Corrected[det]->Clone();
            H_SiPM_Corrected_temp->GetXaxis()->SetRangeUser(Lowmin, Lowmax);
            H_SiPM_Corrected_temp->Scale(H_SiPM_NoCorrected_temp->Integral() / H_SiPM_Corrected_temp->Integral());
            H_SiPM_Corrected_temp->SetLineColor(kRed);
            H_SiPM_Corrected_temp->Draw("HIST SAME");

            cCompareCorrectionLow_SiPM->cd(1);
            H_SiPM_NoCorrected[det]->GetXaxis()->SetRangeUser(Lowmin/10, Lowmax/10);
            H_SiPM_NoCorrected[det]->Scale(H_SiPM_NoCorrected[111]->Integral() / H_SiPM_NoCorrected[det]->Integral());
            H_SiPM_NoCorrected[det]->SetLineColor(GetDetectorChannel(det));
            H_SiPM_NoCorrected[det]->Draw("HIST SAME");
            cCompareCorrectionLow_SiPM->cd(2);
            H_SiPM_Corrected[det]->GetXaxis()->SetRangeUser(Lowmin, Lowmax);
            H_SiPM_Corrected[det]->Scale(H_SiPM_NoCorrected[111]->Integral() / H_SiPM_Corrected[det]->Integral());
            H_SiPM_Corrected[det]->SetLineColor(GetDetectorChannel(det));
            H_SiPM_Corrected[det]->Draw("HIST SAME");
        }
    }

    fread->cd();
    cCompareCorrectionHigh->Write();
    cCompareCorrectionLow->Write();
    cCompareCorrectionHigh_SiPM->Write();
    cCompareCorrectionLow_SiPM->Write();
    fread->Close();
}

void InitSiliconCalibration()
{

    string CalibFileName;

    CalibFileName = "Config_Files/" + to_string(YEAR) + "/Calibration_" + to_string(YEAR) + ".txt";

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
            SiliconCalibrationParameter[det][0] = a;
            SiliconCalibrationParameter[det][1] = b;
            SiliconCalibrationParameter[det][2] = c;
        }
    }
    else
    {
        Error("No Calibration file found");
    }

    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            F_SiliconCalibration[i] = new TF1(("F_SiliconCalibration_" + detectorName[i]).c_str(), "[0] + [1]*x + [2]*x*x", eSiliMin_cal, eSiliMax_cal);
            F_SiliconCalibration[i]->SetParameters(SiliconCalibrationParameter[i][0], SiliconCalibrationParameter[i][1], SiliconCalibrationParameter[i][2]);
        }
    }
}

#endif