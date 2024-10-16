#ifndef SIPM_CALIBRATION_HH
#define SIPM_CALIBRATION_HH

#include "../../../lib/SignalDict/Signal.h"
#include "Detectors.hh"

/// FILE ///
map<string, TFile *> GROUPED_File;
map<string, TFile *> SIMULATED_File;
TFile *CALIBRATED_File;

TTree *Tree;
TTreeReader *Reader;
TTreeReaderArray<Signal> *SiPM_High;
TTreeReaderArray<Signal> *SiPM_Low;
TTreeReaderArray<Signal> *Silicon;
TTreeReaderArray<Signal> *signals;

TTree* Tree_SIMULATED;
TTreeReader *Reader_SIMULATED;  
TTreeReaderValue<int> *Silicon_code;
TTreeReaderValue<double> *Silicon_energy;
TTreeReaderValue<double> *SiPM_energy;


TTree *Tree_MATCHED;
vector<Signal> SiPM;

/// WINDOWS ///
map<string, pair<double, double>[100][SIGNAL_MAX]> WindowsMap;

// TREE PEAKS //
map<string, TTree *[50][SIGNAL_MAX]> Tree_Peaks;
Signal SiPMi;

/// HISTOGRAMS ///
// single
map<string, TH1D *[SIGNAL_MAX]> H_SiPM_High;
map<string, TH1D *[SIGNAL_MAX]> H_SiPM_Low;
// matching low high
map<string, TH2D *[SIGNAL_MAX]> H_SiPM_HighLow;
map<string, TGraph *[SIGNAL_MAX]> G_SiPM_HighLow;
map<string, TH1D *[SIGNAL_MAX]> H_SiPM_HighLow_Projected;
map<string, TGraph *[SIGNAL_MAX]> G_SiPM_HighLow_Projected;
map<string, TH1D *[SIGNAL_MAX]> H_SiPM_HighLow_Matched;
// matching SiPM
map<string, TGraph *[SIGNAL_MAX]> G_Matching_SiPM;
map<string, TH2D*[SIGNAL_MAX]> H_Matching_SiPM;
map<string, TGraph *[SIGNAL_MAX]> G_SiPM_HighLow_Matched;
map<string, TH2D*[SIGNAL_MAX]> H_SiPM_HighLow_Matched_diff;
map<string, TH1D *[SIGNAL_MAX]> H_SiPM_Merged;
map<string, TH1D *[SIGNAL_MAX]> H_SiPM_Matched;
map<string, TH1D *[SIGNAL_MAX]> H_SiPM_UnMatched;
map<string, TH1D *[6]> H_SiPM_LowHighTriggered_x;
map<string, TH1D *> H_SiPM_HighTriggered_x;
map<string, TH1D *> H_SiPM_LowTriggered_x;
map<string, TH1D *[SIGNAL_MAX]> H_Channel_SiPM_High_Alone;
map<string, TH1D *[SIGNAL_MAX]> H_Channel_SiPM_Low_Alone;
// final sipm for each proton peaks
map<string, TH1D *[50][SIGNAL_MAX]> H_SiPM;
map<string, TH1D *[50][SIGNAL_MAX]> H_SiPM_Calibrated;

map<string, TH1D *[50]> H_Sim;
map<string, TH1D *[50][SIGNAL_MAX]> H_Sim_Conv;
/// DIRECTORY ///
map<string, TDirectory *> dir; 
map<string, TDirectory *[SIGNAL_MAX]> dir_detector;
map<string, TDirectory *[SIGNAL_MAX]> dir_SiPM;

// FUCNTIONS //
TF1 *F_SiliconCalibration[SIGNAL_MAX];
TF1 *f_linear = new TF1("f_linear", "[0]*x + [1]", 0, 10000e3);
TF1 *f_erf = new TF1("f_erf", "x < [0] ? [1]*x + [2] : [3]*(erf((x-[4])/[5]))", 0, 10000e3);
pair<double, double> Range_SiPM_LowHigh = make_pair(50e3, 200e3);
TF1* gauss;
TF1* Threshold_f;
int current_detector;
double SiPM_Window[10] = {0, 1500e3, 1500e3, 1500e3};
vector<double> res;
// DATA //
string Nuclei[1] = {"32Ar"}; //, "32Ar_thick", "33Ar"};
string Nucleis[2] = {"32Ar", "207Bi"};
string NUCLEUS;
double SiliconCalibrationParameter[SIGNAL_MAX][3];
map<string, TF1 *[SIGNAL_MAX]> MatchingLowHigh;
map<string, TF1 *[SIGNAL_MAX]> MatchingLowHigh_erf;
map<string, TF1 *[SIGNAL_MAX]> MatchingSiPM;

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

                for (int i : Dir2Det(dir, strip))
                {
                    WindowsMap[nuclei][number][i] = make_pair(energy_low, energy_high);
                    for (int det = 1; det <= 9; det++)
                    {
                        Tree_Peaks[NUCLEUS][number][det] = new TTree(("Tree_Peaks_" + NUCLEUS + "_" + to_string(number) + "_" + to_string(det)).c_str(), ("Tree_Peaks_" + NUCLEUS + "_" + to_string(number) + "_" + to_string(det)).c_str());
                        Tree_Peaks[NUCLEUS][number][det]->Branch("SiPM", &SiPMi);
                    }
                }
            }
        }
    }

    for (int det = 1; det <= 9; det++)
    {
        Tree_Peaks["207Bi"][0][det] = new TTree(("Tree_Peaks_" + NUCLEUS + "_" + to_string(0) + "_" + to_string(det)).c_str(), ("Tree_Peaks_" + NUCLEUS + "_" + to_string(0) + "_" + to_string(det)).c_str());
        Tree_Peaks["207Bi"][0][det]->Branch("SiPM", &SiPMi);
    }
}

pair<double, double> PointProjection(double xA, double yA, TF1 *f)
{
    double a = f->GetParameter(0);
    double b = -1;
    double c = f->GetParameter(1);

    double xB = (b * (b * xA - a * yA) - a * c) / (a * a + b * b);
    double yB = (a * (-b * xA + a * yA) - b * c) / (a * a + b * b);

    return make_pair(xB, yB);
}

void WriteValues()
{
    Info("Write matching values to file");
    TFile *f = new TFile((DIR_ROOT_DATA_CALIBRATED + "Matching_SiPM_values.root").c_str(), "RECREATE");
    f->cd();
    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorBetaHigh(det))
        {
            MatchingLowHigh[NUCLEUS][det]->Write(("MatchingLowHigh_" + NUCLEUS + "_" + detectorName[det]).c_str());
            MatchingSiPM[NUCLEUS][det]->Write(("MatchingSiPM_" + NUCLEUS + "_" + detectorName[det]).c_str());
        }
    }
    f->Close();
    CALIBRATED_File->cd();
}

void ReadValues()
{
    Warning("Read matching values from file");
    TFile *f = new TFile((DIR_ROOT_DATA_CALIBRATED + "Matching_SiPM_values.root").c_str(), "READ");
    NUCLEUS = "32Ar";
    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorBetaHigh(det))
        {
            for (string Nucleus : Nucleis)
            {
                MatchingLowHigh[Nucleus][det] = (TF1 *)f->Get(("MatchingLowHigh_" + NUCLEUS + "_" + detectorName[det]).c_str());
                MatchingSiPM[Nucleus][det] = (TF1 *)f->Get(("MatchingSiPM_" + NUCLEUS + "_" + detectorName[det]).c_str());
            }
        }
    }

    f->Close();
    CALIBRATED_File->cd();
}

void WriteTree()
{
    TFile *f = new TFile((DIR_ROOT_DATA_CALIBRATED + "Matching_SiPM_trees.root").c_str(), "RECREATE");
    f->cd();
    Info("Write Tree");
    NUCLEUS="32Ar";
    for (int peak = 1; peak <= 50; peak++)
    {
        for (int det = 1; det <= 9; det++)
        {
            if (WindowsMap[NUCLEUS][peak][11].first == -1 || !WindowsMap[NUCLEUS][peak][11].first)
                        continue;
            Tree_Peaks[NUCLEUS][peak][det]->Write();
        }
    }
    cout << "207Bi" << endl;
    for (int det = 1; det <= 9; det++)
    {
        Tree_Peaks["207Bi"][0][det]->Write();
    }
    f->Close();
    delete f;
    CALIBRATED_File->cd();
}

void ReadTree()
{
    TFile *f = new TFile((DIR_ROOT_DATA_CALIBRATED + "Matching_SiPM_trees.root").c_str(), "READ");
    Info("Read Tree");
    NUCLEUS="32Ar";
    for (int peak = 0; peak <= 50; peak++)
    {
        for (int det = 1; det <= 9; det++)
        {
            if (WindowsMap[NUCLEUS][peak][11].first == -1 || !WindowsMap[NUCLEUS][peak][11].first)
                        continue;
            Tree_Peaks[NUCLEUS][peak][det] = (TTree *)f->Get(("Tree_Peaks_" + NUCLEUS + "_" + to_string(peak) + "_" + to_string(det)).c_str());
        }
    }
    for (int det = 1; det <= 9; det++)
    {
        Tree_Peaks["207Bi"][0][det] = (TTree *)f->Get(("Tree_Peaks_207Bi_" + to_string(0) + "_" + to_string(det)).c_str());
    }
    // f->Close();
    CALIBRATED_File->cd();
}

void InitHistograms()
{

    Info("Init Histograms");
    for (string Nucleus : Nucleis)
    {
        NUCLEUS = Nucleus;
        dir[NUCLEUS] = CALIBRATED_File->mkdir(NUCLEUS.c_str());
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorBetaHigh(i))
            {
                H_SiPM_High[NUCLEUS][i] = new TH1D(("H_SiPM_High_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_SiPM_High_" + NUCLEUS + "_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
                H_SiPM_High[NUCLEUS][i]->GetXaxis()->SetTitle("Channel");
                H_SiPM_High[NUCLEUS][i]->GetYaxis()->SetTitle("Counts");
                H_SiPM_High[NUCLEUS][i]->GetXaxis()->CenterTitle();
                H_SiPM_High[NUCLEUS][i]->GetYaxis()->CenterTitle();
                H_SiPM_High[NUCLEUS][i]->SetStats(0);

                dir_detector[NUCLEUS][i] = dir[NUCLEUS]->mkdir(("SiPM_" + to_string(GetDetectorChannel(i))).c_str());
                H_SiPM_HighLow[NUCLEUS][i] = new TH2D(("H_SiPM_HighLow_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_SiPM_HighLow_" + NUCLEUS + "_" + detectorName[i]).c_str(), eLowN, eLowMin, eLowMax, eHighN, eHighMin, eHighMax);
                H_SiPM_HighLow[NUCLEUS][i]->GetXaxis()->SetTitle("Channel Low");
                H_SiPM_HighLow[NUCLEUS][i]->GetYaxis()->SetTitle("Channel High");
                H_SiPM_HighLow[NUCLEUS][i]->GetXaxis()->CenterTitle();
                H_SiPM_HighLow[NUCLEUS][i]->GetYaxis()->CenterTitle();
                H_SiPM_HighLow[NUCLEUS][i]->SetStats(0);

                G_SiPM_HighLow[NUCLEUS][i] = new TGraph();
                G_SiPM_HighLow_Projected[NUCLEUS][i] = new TGraph();
                G_SiPM_HighLow_Matched[NUCLEUS][i] = new TGraph();

                H_SiPM_HighLow_Projected[NUCLEUS][i] = new TH1D(("H_SiPM_HighLow_Projected_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_SiPM_HighLow_Projected_" + NUCLEUS + "_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
                H_SiPM_HighLow_Projected[NUCLEUS][i]->GetXaxis()->SetTitle("Channel");
                H_SiPM_HighLow_Projected[NUCLEUS][i]->GetYaxis()->SetTitle("Counts");
                H_SiPM_HighLow_Projected[NUCLEUS][i]->GetXaxis()->CenterTitle();
                H_SiPM_HighLow_Projected[NUCLEUS][i]->GetYaxis()->CenterTitle();
                H_SiPM_HighLow_Projected[NUCLEUS][i]->SetStats(0);

                // H_SiPM_HighLow_Matched_diff[NUCLEUS][i] = new TH2D(("H_SiPM_HighLow_Matched_diff_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_SiPM_HighLow_Matched_diff_" + NUCLEUS + "_" + detectorName[i]).c_str(), 2000, -100000, 100000, eHighN/5, eHighMin/5, eHighMax/5);
                // H_SiPM_HighLow_Matched_diff[NUCLEUS][i]->GetXaxis()->SetTitle("Channel Low - Channel High");
                // H_SiPM_HighLow_Matched_diff[NUCLEUS][i]->GetYaxis()->SetTitle("Channel High");
                // H_SiPM_HighLow_Matched_diff[NUCLEUS][i]->GetXaxis()->CenterTitle();
                // H_SiPM_HighLow_Matched_diff[NUCLEUS][i]->GetYaxis()->CenterTitle();
                // H_SiPM_HighLow_Matched_diff[NUCLEUS][i]->SetStats(0);

                H_Matching_SiPM[NUCLEUS][GetDetectorChannel(i)] = new TH2D(("H_Matching_SiPM_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_Matching_SiPM_" + NUCLEUS + "_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax, eHighN, eHighMin, eHighMax);
                H_Matching_SiPM[NUCLEUS][GetDetectorChannel(i)]->GetXaxis()->SetTitle("Channel Low - Channel High");
                H_Matching_SiPM[NUCLEUS][GetDetectorChannel(i)]->GetYaxis()->SetTitle("Channel High");
                H_Matching_SiPM[NUCLEUS][GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
                H_Matching_SiPM[NUCLEUS][GetDetectorChannel(i)]->GetYaxis()->CenterTitle();
                H_Matching_SiPM[NUCLEUS][GetDetectorChannel(i)]->SetStats(0);

                G_Matching_SiPM[NUCLEUS][GetDetectorChannel(i)] = new TGraph();
                

                H_SiPM_HighLow_Matched[NUCLEUS][i] = new TH1D(("H_SiPM_HighLow_Matched_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_SiPM_HighLow_Matched_" + NUCLEUS + "_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
                H_SiPM_HighLow_Matched[NUCLEUS][i]->GetXaxis()->SetTitle("Channel");
                H_SiPM_HighLow_Matched[NUCLEUS][i]->GetYaxis()->SetTitle("Counts");
                H_SiPM_HighLow_Matched[NUCLEUS][i]->GetXaxis()->CenterTitle();
                H_SiPM_HighLow_Matched[NUCLEUS][i]->GetYaxis()->CenterTitle();
                H_SiPM_HighLow_Matched[NUCLEUS][i]->SetStats(0);

                H_SiPM_Merged[NUCLEUS][i] = new TH1D(("H_SiPM_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_SiPM_" + NUCLEUS + "_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
                H_SiPM_Merged[NUCLEUS][i]->GetXaxis()->SetTitle("Channel");
                H_SiPM_Merged[NUCLEUS][i]->GetYaxis()->SetTitle("Counts");
                H_SiPM_Merged[NUCLEUS][i]->GetXaxis()->CenterTitle();
                H_SiPM_Merged[NUCLEUS][i]->GetYaxis()->CenterTitle();
                H_SiPM_Merged[NUCLEUS][i]->SetStats(0);

                H_SiPM_Matched[NUCLEUS][i] = new TH1D(("H_SiPM_Matched_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_SiPM_Matched_" + NUCLEUS + "_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
                H_SiPM_Matched[NUCLEUS][i]->GetXaxis()->SetTitle("Channel");
                H_SiPM_Matched[NUCLEUS][i]->GetYaxis()->SetTitle("Counts");
                H_SiPM_Matched[NUCLEUS][i]->GetXaxis()->CenterTitle();
                H_SiPM_Matched[NUCLEUS][i]->GetYaxis()->CenterTitle();
                H_SiPM_Matched[NUCLEUS][i]->SetStats(0);

                H_SiPM_UnMatched[NUCLEUS][i] = new TH1D(("H_SiPM_UnMatched_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_SiPM_UnMatched_" + NUCLEUS + "_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
                H_SiPM_UnMatched[NUCLEUS][i]->GetXaxis()->SetTitle("Channel");
                H_SiPM_UnMatched[NUCLEUS][i]->GetYaxis()->SetTitle("Counts");
                H_SiPM_UnMatched[NUCLEUS][i]->GetXaxis()->CenterTitle();
                H_SiPM_UnMatched[NUCLEUS][i]->GetYaxis()->CenterTitle();
                H_SiPM_UnMatched[NUCLEUS][i]->SetStats(0);

                H_Channel_SiPM_High_Alone[NUCLEUS][i] = new TH1D(("H_Channel_SiPM_High_Alone_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_Channel_SiPM_High_Alone_" + NUCLEUS + "_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
                H_Channel_SiPM_High_Alone[NUCLEUS][i]->GetXaxis()->SetTitle("Channel");
                H_Channel_SiPM_High_Alone[NUCLEUS][i]->GetYaxis()->SetTitle("Counts");
                H_Channel_SiPM_High_Alone[NUCLEUS][i]->GetXaxis()->CenterTitle();
                H_Channel_SiPM_High_Alone[NUCLEUS][i]->GetYaxis()->CenterTitle();
                H_Channel_SiPM_High_Alone[NUCLEUS][i]->SetStats(0);

                
                for (int peak = 0; peak <= 50; peak++)
                {

                    if ((WindowsMap[NUCLEUS][peak][11].first == -1 || !WindowsMap[NUCLEUS][peak][11].first))
                        continue;

                    H_SiPM[NUCLEUS][peak][GetDetectorChannel(i)] = new TH1D(("H_SiPM_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), ("H_SiPM_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), eHighN, eHighMin, eHighMax);
                    H_SiPM[NUCLEUS][peak][GetDetectorChannel(i)]->GetXaxis()->SetTitle("Channel");
                    H_SiPM[NUCLEUS][peak][GetDetectorChannel(i)]->GetYaxis()->SetTitle("Counts");
                    H_SiPM[NUCLEUS][peak][GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
                    H_SiPM[NUCLEUS][peak][GetDetectorChannel(i)]->GetYaxis()->CenterTitle();


                    H_SiPM_Calibrated[NUCLEUS][peak][GetDetectorChannel(i)] = new TH1D(("H_SiPM_Calibrated_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), ("H_SiPM_Calibrated_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
                    H_SiPM_Calibrated[NUCLEUS][peak][GetDetectorChannel(i)]->GetXaxis()->SetTitle("Energy [keV]");
                    H_SiPM_Calibrated[NUCLEUS][peak][GetDetectorChannel(i)]->GetYaxis()->SetTitle("Counts");
                    H_SiPM_Calibrated[NUCLEUS][peak][GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
                    H_SiPM_Calibrated[NUCLEUS][peak][GetDetectorChannel(i)]->GetYaxis()->CenterTitle();
                }

                if (NUCLEUS == "207Bi")
                {
                    H_SiPM[NUCLEUS][0][GetDetectorChannel(i)] = new TH1D(("H_SiPM_" + NUCLEUS + "_Peak_" + to_string(0) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), ("H_SiPM_" + NUCLEUS + "_Peak_" + to_string(0) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), eHighN, eHighMin, eHighMax);
                    H_SiPM[NUCLEUS][0][GetDetectorChannel(i)]->GetXaxis()->SetTitle("Channel");
                    H_SiPM[NUCLEUS][0][GetDetectorChannel(i)]->GetYaxis()->SetTitle("Counts");
                    H_SiPM[NUCLEUS][0][GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
                    H_SiPM[NUCLEUS][0][GetDetectorChannel(i)]->GetYaxis()->CenterTitle();

                    H_SiPM_Calibrated[NUCLEUS][0][GetDetectorChannel(i)] = new TH1D(("H_SiPM_Calibrated_" + NUCLEUS + "_Peak_" + to_string(0) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), ("H_SiPM_Calibrated_" + NUCLEUS + "_Peak_" + to_string(0) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
                    H_SiPM_Calibrated[NUCLEUS][0][GetDetectorChannel(i)]->GetXaxis()->SetTitle("Energy [keV]");
                    H_SiPM_Calibrated[NUCLEUS][0][GetDetectorChannel(i)]->GetYaxis()->SetTitle("Counts");
                    H_SiPM_Calibrated[NUCLEUS][0][GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
                    H_SiPM_Calibrated[NUCLEUS][0][GetDetectorChannel(i)]->GetYaxis()->CenterTitle();
                }
            }

            if (IsDetectorBetaLow(i))
            {
                H_SiPM_Low[NUCLEUS][i] = new TH1D(("H_SiPM_Low_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_SiPM_Low_" + NUCLEUS + "_" + detectorName[i]).c_str(), eLowN, eLowMin, eLowMax);
                H_SiPM_Low[NUCLEUS][i]->GetXaxis()->SetTitle("Channel");
                H_SiPM_Low[NUCLEUS][i]->GetYaxis()->SetTitle("Counts");
                H_SiPM_Low[NUCLEUS][i]->GetXaxis()->CenterTitle();
                H_SiPM_Low[NUCLEUS][i]->GetYaxis()->CenterTitle();

                H_SiPM_HighLow_Matched[NUCLEUS][i] = new TH1D(("H_SiPM_HighLow_Matched_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_SiPM_HighLow_Matched_" + NUCLEUS + "_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
                H_SiPM_HighLow_Matched[NUCLEUS][i]->GetXaxis()->SetTitle("Channel");
                H_SiPM_HighLow_Matched[NUCLEUS][i]->GetYaxis()->SetTitle("Counts");
                H_SiPM_HighLow_Matched[NUCLEUS][i]->GetXaxis()->CenterTitle();
                H_SiPM_HighLow_Matched[NUCLEUS][i]->GetYaxis()->CenterTitle();

                H_Channel_SiPM_Low_Alone[NUCLEUS][i] = new TH1D(("H_Channel_SiPM_Low_Alone_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_Channel_SiPM_Low_Alone_" + NUCLEUS + "_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
                H_Channel_SiPM_Low_Alone[NUCLEUS][i]->GetXaxis()->SetTitle("Channel");
                H_Channel_SiPM_Low_Alone[NUCLEUS][i]->GetYaxis()->SetTitle("Counts");
                H_Channel_SiPM_Low_Alone[NUCLEUS][i]->GetXaxis()->CenterTitle();
                H_Channel_SiPM_Low_Alone[NUCLEUS][i]->GetYaxis()->CenterTitle();

                H_SiPM_Matched[NUCLEUS][i] = new TH1D(("H_SiPM_Matched_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_SiPM_Matched_" + NUCLEUS + "_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
                H_SiPM_Matched[NUCLEUS][i]->GetXaxis()->SetTitle("Channel");
                H_SiPM_Matched[NUCLEUS][i]->GetYaxis()->SetTitle("Counts");
                H_SiPM_Matched[NUCLEUS][i]->GetXaxis()->CenterTitle();
                H_SiPM_Matched[NUCLEUS][i]->GetYaxis()->CenterTitle();
            }
        }

        for (int peak = 0; peak <= 50; peak++)
        {
            if (WindowsMap[NUCLEUS][peak][11].first == -1 || !WindowsMap[NUCLEUS][peak][11].first)
                continue;

            H_Sim[NUCLEUS][peak] = new TH1D(("H_Sim_" + NUCLEUS + "_Peak_" + to_string(peak)).c_str(), ("H_Sim_" + NUCLEUS + "_Peak_" + to_string(peak)).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
            H_Sim[NUCLEUS][peak]->GetXaxis()->SetTitle("Energy [keV]");
            H_Sim[NUCLEUS][peak]->GetYaxis()->SetTitle("Counts");
            H_Sim[NUCLEUS][peak]->GetXaxis()->CenterTitle();
            H_Sim[NUCLEUS][peak]->GetYaxis()->CenterTitle();
            for (int det = 1; det <= 9; det++)
            {

                H_Sim_Conv[NUCLEUS][peak][det] = new TH1D(("H_Sim_Conv_" + NUCLEUS + "_Peak_" + to_string(peak) + detectorName[100+det]).c_str(), ("H_Sim_Conv_" + NUCLEUS + "_Peak_" + to_string(peak) + detectorName[100+det]).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
                H_Sim_Conv[NUCLEUS][peak][det]->GetXaxis()->SetTitle("Energy [keV]");
                H_Sim_Conv[NUCLEUS][peak][det]->GetYaxis()->SetTitle("Counts");
                H_Sim_Conv[NUCLEUS][peak][det]->GetXaxis()->CenterTitle();
                H_Sim_Conv[NUCLEUS][peak][det]->GetYaxis()->CenterTitle();
            }
        }

        if (NUCLEUS == "207Bi")
        {
            H_Sim[NUCLEUS][0] = new TH1D(("H_Sim_" + NUCLEUS + "_Peak_" + to_string(0)).c_str(), ("H_Sim_" + NUCLEUS + "_Peak_" + to_string(0)).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
            H_Sim[NUCLEUS][0]->GetXaxis()->SetTitle("Energy [keV]");
            H_Sim[NUCLEUS][0]->GetYaxis()->SetTitle("Counts");
            H_Sim[NUCLEUS][0]->GetXaxis()->CenterTitle();
            H_Sim[NUCLEUS][0]->GetYaxis()->CenterTitle();
            for (int det = 1; det <= 9; det++)
            {

                H_Sim_Conv[NUCLEUS][0][det] = new TH1D(("H_Sim_Conv_" + NUCLEUS + "_Peak_" + to_string(0) + detectorName[100+det]).c_str(), ("H_Sim_Conv_" + NUCLEUS + "_Peak_" + to_string(0) + detectorName[100+det]).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
                H_Sim_Conv[NUCLEUS][0][det]->GetXaxis()->SetTitle("Energy [keV]");
                H_Sim_Conv[NUCLEUS][0][det]->GetYaxis()->SetTitle("Counts");
                H_Sim_Conv[NUCLEUS][0][det]->GetXaxis()->CenterTitle();
                H_Sim_Conv[NUCLEUS][0][det]->GetYaxis()->CenterTitle();
            }
        }

        for (int det = 1; det <= 9; det++)
        {
            dir_SiPM[NUCLEUS][det] = CALIBRATED_File->mkdir(("SiPM_" + NUCLEUS + to_string(det)).c_str());
        }
    }
}

void FittingLowHigh()
{
    ////////////////////////// FITTING Low/High //////////////////////////
    Info("Fitting Low/High SiPMs");
    clock_t start = clock(), Current;
    int counter[SIGNAL_MAX] = {0};

    while (Reader->Next() && Reader->GetCurrentEntry() < 1000000)
    {

        ProgressBar(Reader->GetCurrentEntry(), 1000000, start, Current, "Reading Tree");

        int sili_label = (*Silicon)[1].Label;
        double sili_energy = F_SiliconCalibration[sili_label]->Eval((*Silicon)[1].Channel / 1000);

        // gate on proton peak
        if (sili_energy < WindowsMap[NUCLEUS][14][sili_label].first || sili_energy > WindowsMap[NUCLEUS][14][sili_label].second)
            continue;

        for (int i = 0; i < SiPM_High->GetSize(); i++)
        {
            H_SiPM_High[NUCLEUS][(*SiPM_High)[i].Label]->Fill((*SiPM_High)[i].Channel);
            for (int j = 0; j < SiPM_Low->GetSize(); j++)
            {
                if (GetDetectorChannel((*SiPM_High)[i].Label) == GetDetectorChannel((*SiPM_Low)[j].Label))
                {
                    H_SiPM_HighLow[NUCLEUS][(*SiPM_High)[i].Label]->Fill((*SiPM_Low)[j].Channel, (*SiPM_High)[i].Channel);
                    G_SiPM_HighLow[NUCLEUS][(*SiPM_High)[i].Label]->SetPoint(counter[(*SiPM_High)[i].Label], (*SiPM_Low)[j].Channel, (*SiPM_High)[i].Channel);
                    counter[(*SiPM_High)[i].Label]++;
                }
            }
        }

        for (int i = 0; i < SiPM_Low->GetSize(); i++)
        {
            H_SiPM_Low[NUCLEUS][(*SiPM_Low)[i].Label]->Fill((*SiPM_Low)[i].Channel);
        }
    }

    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBetaHigh(i))
        {
            G_SiPM_HighLow[NUCLEUS][i]->Fit(f_linear, "QR", "", Range_SiPM_LowHigh.first, Range_SiPM_LowHigh.second);
            MatchingLowHigh[NUCLEUS][i] = G_SiPM_HighLow[NUCLEUS][i]->GetFunction("f_linear");

            // f_erf->SetParameter(0, 225e3);
            // f_erf->SetParLimits(0, 150e3, 300e3); 
            // f_erf->SetParameter(1, 10);
            // f_erf->SetParLimits(1, 9.5, 12);
            // f_erf->SetParameter(2, -50000);
            // f_erf->SetParLimits(2, -100000, 0);
            // f_erf->SetParameter(3, 5000e3);
            // f_erf->SetParLimits(3, 4000e3, 6000e3);
            // f_erf->SetParameter(4, 0);
            // f_erf->SetParLimits(4, -100000, 100000);
            // f_erf->SetParameter(5, 280e3);
            // f_erf->SetParLimits(5, 0, 600e3);
            // G_SiPM_HighLow[NUCLEUS][i]->Fit(f_erf);
            // MatchingLowHigh_erf[NUCLEUS][i] = G_SiPM_HighLow[NUCLEUS][i]->GetFunction("f_erf");

        }
    }
}

void FittingLowHighBi207()
{
    ////////////////////////// FITTING Low/High //////////////////////////
    Info("Fitting Low/High SiPMs");
    clock_t start = clock(), Current;
    int counter[SIGNAL_MAX] = {0};
    Tree = (TTree *)GROUPED_File[NUCLEUS]->Get("Tree_Group");
    Reader = new TTreeReader(Tree);
    signals = new TTreeReaderArray<Signal>(*Reader, "Signal");
    Reader->Restart();
    while (Reader->Next() && Reader->GetCurrentEntry() < 10000000)
    {

        ProgressBar(Reader->GetCurrentEntry(), 1000000, start, Current, "Reading Tree");
        
        vector<Signal> SiPM_High;
        vector<Signal> SiPM_Low;
        for (int i = 0; i < signals->GetSize(); i++)
        {
            if (IsDetectorBetaHigh((*signals)[i].Label))
            {
                SiPM_High.push_back((*signals)[i]);
            }
            if (IsDetectorBetaLow((*signals)[i].Label))
            {
                SiPM_Low.push_back((*signals)[i]);
            }
        }

        for (int i = 0; i < SiPM_High.size(); i++)
        {
            H_SiPM_High[NUCLEUS][(SiPM_High)[i].Label]->Fill((SiPM_High)[i].Channel);
            for (int j = 0; j < SiPM_Low.size(); j++)
            {
                if (GetDetectorChannel((SiPM_High)[i].Label) == GetDetectorChannel((SiPM_Low)[j].Label))
                {
                    H_SiPM_HighLow[NUCLEUS][(SiPM_High)[i].Label]->Fill((SiPM_Low)[j].Channel, (SiPM_High)[i].Channel);
                    G_SiPM_HighLow[NUCLEUS][(SiPM_High)[i].Label]->SetPoint(counter[(SiPM_High)[i].Label], (SiPM_Low)[j].Channel, (SiPM_High)[i].Channel);
                    counter[(SiPM_High)[i].Label]++;
                }
            }
        }

        for (int i = 0; i < SiPM_Low.size(); i++)
        {
            H_SiPM_Low[NUCLEUS][(SiPM_Low)[i].Label]->Fill((SiPM_Low)[i].Channel);
        }
    }

    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBetaHigh(i))
        {
            G_SiPM_HighLow[NUCLEUS][i]->Fit(f_linear, "QR", "", Range_SiPM_LowHigh.first, Range_SiPM_LowHigh.second);
            MatchingLowHigh[NUCLEUS][i] = G_SiPM_HighLow[NUCLEUS][i]->GetFunction("f_linear");

            

            // f_erf->SetParameter(0, 225e3);
            // f_erf->SetParLimits(0, 150e3, 300e3); 
            // f_erf->SetParameter(1, 10);
            // f_erf->SetParLimits(1, 9.5, 12);
            // f_erf->SetParameter(2, -50000);
            // f_erf->SetParLimits(2, -100000, 0);
            // f_erf->SetParameter(3, 5000e3);
            // f_erf->SetParLimits(3, 4000e3, 6000e3);
            // f_erf->SetParameter(4, 0);
            // f_erf->SetParLimits(4, -100000, 100000);
            // f_erf->SetParameter(5, 280e3);
            // f_erf->SetParLimits(5, 0, 600e3);
            // G_SiPM_HighLow[NUCLEUS][i]->Fit(f_erf, "Q");
            // MatchingLowHigh_erf[NUCLEUS][i] = G_SiPM_HighLow[NUCLEUS][i]->GetFunction("f_erf");

        }
    }

}

void FittingSiPMs()
{
    ////////////////////////// FITTING SiPMs //////////////////////////
    Info("Matching Low/High SiPM & Fitting SiPMs");
    clock_t start = clock(), Current;
    Reader->Restart();
    int counter_graph[SIGNAL_MAX] = {0};
    int counter[SIGNAL_MAX] = {0};
    while (Reader->Next() && Reader->GetCurrentEntry() < 1000000)
    {
        ProgressBar(Reader->GetCurrentEntry(), 1000000, start, Current, "Reading Tree");

        int sili_label = (*Silicon)[1].Label;
        double sili_energy = F_SiliconCalibration[sili_label]->Eval((*Silicon)[1].Channel / 1000);

        // gate on proton peak
        if (sili_energy < WindowsMap[NUCLEUS][14][sili_label].first || sili_energy > WindowsMap[NUCLEUS][14][sili_label].second)
            continue;

        /// WRITTING MATCHED LOW HIGH
        for (int i = 0; i < SiPM_High->GetSize(); i++)
        {
            H_SiPM_HighLow_Matched[NUCLEUS][(*SiPM_High)[i].Label]->Fill((*SiPM_High)[i].Channel);
        }

        for (int j = 0; j < SiPM_Low->GetSize(); j++)
        {
            H_SiPM_HighLow_Matched[NUCLEUS][(*SiPM_Low)[j].Label]->Fill(MatchingLowHigh["32Ar"][(*SiPM_Low)[j].Label - 10]->Eval((*SiPM_Low)[j].Channel));
        }

        /// double line computing
        for (int i = 0; i < SiPM_High->GetSize(); i++)
        {
            for (int j = 0; j < SiPM_Low->GetSize(); j++)
            {
                if (GetDetectorChannel((*SiPM_High)[i].Label) == GetDetectorChannel((*SiPM_Low)[j].Label))
                {
                    double high = (*SiPM_High)[i].Channel;
                    double low = MatchingLowHigh[NUCLEUS][(*SiPM_High)[i].Label]->Eval((*SiPM_Low)[j].Channel);

                    G_SiPM_HighLow_Matched[NUCLEUS][(*SiPM_High)[i].Label]->SetPoint(counter[(*SiPM_High)[i].Label], low, high);
                    counter[(*SiPM_High)[i].Label]++;

                    // H_SiPM_HighLow_Matched_diff[NUCLEUS][(*SiPM_High)[i].Label]->Fill(MatchingLowHigh[NUCLEUS][(*SiPM_High)[i].Label]->Eval((*SiPM_Low)[j].Channel)-high, (*SiPM_Low)[j].Channel);
                }
            }
        }
        ///////////////////////////////

        /// MATCHING SiPMs
        double SiPM[10] = {0};
        for (int i = 0; i < SiPM_High->GetSize(); i++)
        {
            for (int j = 0; j < SiPM_Low->GetSize(); j++)
            {
                if (GetDetectorChannel((*SiPM_High)[i].Label) == GetDetectorChannel((*SiPM_Low)[j].Label))
                {
                    if ((*SiPM_High)[i].Channel < Range_SiPM_LowHigh.second)
                    {
                        SiPM[GetDetectorChannel((*SiPM_High)[i].Label)] = (*SiPM_High)[i].Channel;
                    }
                    else
                    {
                        SiPM[GetDetectorChannel((*SiPM_High)[i].Label)] = MatchingLowHigh[NUCLEUS][(*SiPM_High)[i].Label]->Eval((*SiPM_Low)[j].Channel);
                    }
                }
            }
        }

        if (SiPM[1] != 0)
        {
            for (int index = 1; index <= 9; index++)
            {
                if (SiPM[index] == 0.)
                    continue;
                G_Matching_SiPM[NUCLEUS][index]->SetPoint(counter_graph[index], SiPM[index], SiPM[1]);
                H_Matching_SiPM[NUCLEUS][index]->Fill(SiPM[index], SiPM[1]);
                counter_graph[index]++;
            }
        }
        //////////////////////
    }

    /// fitting sipms matching
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBetaHigh(i))
        {
            cout << "Fitting SiPM " << i << endl;
            f_linear->FixParameter(1, 0);
            G_Matching_SiPM[NUCLEUS][GetDetectorChannel(i)]->Fit(f_linear, "QR", "", Range_SiPM_LowHigh.first*10, Range_SiPM_LowHigh.second*10);
            MatchingSiPM[NUCLEUS][i] = G_Matching_SiPM[NUCLEUS][GetDetectorChannel(i)]->GetFunction("f_linear");

            // TGraphErrors* g = new TGraphErrors();
            // int counter = 0;
            // for (int bin = 0; bin < H_Matching_SiPM[NUCLEUS][GetDetectorChannel(i)]->GetNbinsY(); bin+=10)
            // {
            //     TH1D* h = (TH1D*)H_Matching_SiPM[NUCLEUS][GetDetectorChannel(i)]->ProjectionY(("ap" + to_string(i)).c_str(), bin, bin+10);
            //     if (h->GetEntries() < 100)
            //         continue;
            //     // h->Write(); 
            //     h->Fit("gaus", "Q");
            //     g->SetPoint(counter, H_Matching_SiPM[NUCLEUS][GetDetectorChannel(i)]->GetXaxis()->GetBinCenter(bin+50), h->GetFunction("gaus")->GetParameter(1));
            //     g->SetPointError(counter, H_Matching_SiPM[NUCLEUS][GetDetectorChannel(i)]->GetXaxis()->GetBinCenter(25), h->GetFunction("gaus")->GetParError(1));
            //     counter++;
            //     delete h;
            // }
            // f_linear->FixParameter(1, 0);
            // g->Fit(f_linear, "QR", "", Range_SiPM_LowHigh.first * 10 - 200e3, Range_SiPM_LowHigh.second * 10 - 400e3);
            // MatchingSiPM[NUCLEUS][i] = g->GetFunction("f_linear");
            // g->Write();

        }
    }

    MatchingSiPM[NUCLEUS][101] = new TF1("MatchingSiPM_101", "[0]*x + [1]", 0, 10000e3);
    MatchingSiPM[NUCLEUS][101]->SetParameters(1, 0);
}

void FittingSiPMsBi207()
{
    ////////////////////////// FITTING SiPMs //////////////////////////
    Info("Matching Low/High SiPM & Fitting SiPMs");
    clock_t start = clock(), Current;
    Reader->Restart();
    int counter_graph[SIGNAL_MAX] = {0};
    int counter[SIGNAL_MAX] = {0};
    Tree = (TTree *)GROUPED_File[NUCLEUS]->Get("Tree_Group");
    Reader = new TTreeReader(Tree);
    signals = new TTreeReaderArray<Signal>(*Reader, "Signal");
    Reader->Restart();
    while (Reader->Next() && Reader->GetCurrentEntry() < 10000000)
    {
        ProgressBar(Reader->GetCurrentEntry(), 1000000, start, Current, "Reading Tree");

        vector<Signal> SiPM_High;
        vector<Signal> SiPM_Low;
        double ref_time;
        for (int i = 0; i < signals->GetSize(); i++)
        {
            if (i == 0)
                ref_time = (*signals)[i].Time;
            if (IsDetectorBetaHigh((*signals)[i].Label))
            {
                SiPM_High.push_back((*signals)[i]);
            }
            if (IsDetectorBetaLow((*signals)[i].Label))
            {
                SiPM_Low.push_back((*signals)[i]);
            }
        }

        /// WRITTING MATCHED LOW HIGH
        for (int i = 0; i < SiPM_High.size(); i++)
        {
            H_SiPM_HighLow_Matched[NUCLEUS][(SiPM_High)[i].Label]->Fill((SiPM_High)[i].Channel);
        }

        for (int j = 0; j < SiPM_Low.size(); j++)
        {
            H_SiPM_HighLow_Matched[NUCLEUS][(SiPM_Low)[j].Label]->Fill(MatchingLowHigh["32Ar"][(SiPM_Low)[j].Label - 10]->Eval((SiPM_Low)[j].Channel));
        }

        /// double line computing
        for (int i = 0; i < SiPM_High.size(); i++)
        {
            for (int j = 0; j < SiPM_Low.size(); j++)
            {
                if (GetDetectorChannel((SiPM_High)[i].Label) == GetDetectorChannel((SiPM_Low)[j].Label))
                {
                    double high = (SiPM_High)[i].Channel;
                    double low = MatchingLowHigh["32Ar"][(SiPM_High)[i].Label]->Eval((SiPM_Low)[j].Channel);

                    G_SiPM_HighLow_Matched[NUCLEUS][(SiPM_High)[i].Label]->SetPoint(counter[(SiPM_High)[i].Label], low, high);
                    counter[(SiPM_High)[i].Label]++;

                    // H_SiPM_HighLow_Matched_diff[NUCLEUS][(SiPM_High)[i].Label]->Fill(MatchingLowHigh["32Ar"][(SiPM_High)[i].Label]->Eval((SiPM_Low)[j].Channel)-high, (SiPM_Low)[j].Channel);
                    continue;
                }
            }
        }
        ///////////////////////////////

        /// MATCHING SiPMs
        Signal SiPM[10] = {Signal()};
        for (int i = 0; i < SiPM_High.size(); i++)
        {
            // for (int j = 0; j < SiPM_Low.size(); j++)
            // {
                // if ((GetDetectorChannel((SiPM_High)[i].Label) == GetDetectorChannel((SiPM_Low)[j].Label)))
                // {
                    // if ((SiPM_High)[i].Channel < Range_SiPM_LowHigh.second)
                    // {
                        SiPM[GetDetectorChannel((SiPM_High)[i].Label)] = (SiPM_High)[i];
                    // }
                    // else
                    // {
                    //     SiPM[GetDetectorChannel((SiPM_High)[i].Label)] = MatchingLowHigh["32Ar"][(SiPM_High)[i].Label]->Eval((SiPM_Low)[j].Channel);
                    // }
                    // break;
                // }
            // }
        }

        if (SiPM[1].isValid)
        {
            for (int index = 1; index <= 9; index++)
            {
                if (!SiPM[index].isValid)
                    continue;
                if (abs((SiPM[1].Time - (SiPM)[index].Time)) < 10)
                {
                    G_Matching_SiPM[NUCLEUS][index]->SetPoint(counter_graph[index], SiPM[index].Channel, SiPM[1].Channel);
                    counter_graph[index]++;
                }
            }
        }
        //////////////////////
    }

    /// fitting sipms matching
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBetaHigh(i))
        {
            f_linear->FixParameter(1, 0);
            G_Matching_SiPM[NUCLEUS][GetDetectorChannel(i)]->Fit(f_linear, "Q");
            MatchingSiPM[NUCLEUS][i] = G_Matching_SiPM[NUCLEUS][GetDetectorChannel(i)]->GetFunction("f_linear");
        }
    }

    MatchingSiPM[NUCLEUS][101] = new TF1("MatchingSiPM_101", "[0]*x + [1]", 0, 10000e3);
    MatchingSiPM[NUCLEUS][101]->SetParameters(1, 0);
}

void MergingSiPMs()
{
    Info("Matching SiPMs");

    clock_t start = clock(), Current;
    
    Reader->Restart();
    while (Reader->Next())
    {
        if (Reader->GetCurrentEntry() < Reader->GetEntries() - 1000000)
            continue;
        ProgressBar(Reader->GetCurrentEntry(), Reader->GetEntries(), start, Current, "Reading Tree");

        int sili_label = (*Silicon)[1].Label;
        double sili_energy = F_SiliconCalibration[sili_label]->Eval((*Silicon)[1].Channel / 1000);

        int peak_number = 0;    
        for (int i = 1; i < 100; i++)
        {
            if (WindowsMap[NUCLEUS][i][sili_label].first == -1 || !WindowsMap[NUCLEUS][i][sili_label].first)
                continue;

            if (sili_energy > WindowsMap[NUCLEUS][i][sili_label].first && sili_energy < WindowsMap[NUCLEUS][i][sili_label].second)
            {
                peak_number = i;
                break;
            }
        }

        if (peak_number == 0)
            continue;

        Signal High[10] = {Signal()};
        Signal Low[10] = {Signal()};

        for (int i = 0; i < SiPM_High->GetSize(); i++)
        {
            (*SiPM_High)[i].Channel = MatchingSiPM[NUCLEUS][100 + GetDetectorChannel((*SiPM_High)[i].Label)]->Eval((*SiPM_High)[i].Channel);
            High[GetDetectorChannel((*SiPM_High)[i].Label)] = (*SiPM_High)[i];
            if (peak_number == 14)
                H_SiPM_Matched[NUCLEUS][(*SiPM_High)[i].Label]->Fill((*SiPM_High)[i].Channel);
        }
        for (int i = 0; i < SiPM_Low->GetSize(); i++)
        {
            (*SiPM_Low)[i].Channel = MatchingSiPM[NUCLEUS][100 + GetDetectorChannel((*SiPM_Low)[i].Label)]->Eval(MatchingLowHigh[NUCLEUS][100 + GetDetectorChannel((*SiPM_Low)[i].Label)]->Eval((*SiPM_Low)[i].Channel));
            Low[GetDetectorChannel((*SiPM_Low)[i].Label)] = (*SiPM_Low)[i];
            if (peak_number == 14)
                H_SiPM_Matched[NUCLEUS][(*SiPM_Low)[i].Label]->Fill((*SiPM_Low)[i].Channel);
        }

        for (int det = 1; det <= 9; det++)
        {
            
            if (High[det].isValid && Low[det].isValid) /// BOTH
            {
                if (High[det].Channel < 800e3)
                {
                    H_SiPM[NUCLEUS][peak_number][det]->Fill(High[det].Channel);
                    // SiPM.push_back(High[det]);
                    SiPMi = High[det];
                    Tree_Peaks[NUCLEUS][peak_number][GetDetectorChannel(SiPMi.Label)]->Fill();
                }
                else //(Low[det].Channel > 850e3)
                {
                    H_SiPM[NUCLEUS][peak_number][det]->Fill(Low[det].Channel);
                    // SiPM.push_back(Low[det]);
                    SiPMi = Low[det];
                    Tree_Peaks[NUCLEUS][peak_number][GetDetectorChannel(SiPMi.Label)]->Fill();
                }
                // / make the ratio betwwen them in the range 750 and 850 to get smooth transition
                // else
                // {
                //     double ratio = (Low[det].Channel - 750e3) / (850e3 - 750e3);
                //     H_SiPM[NUCLEUS][peak_number][det]->Fill(High[det].Channel * (1 - ratio) + Low[det].Channel * ratio);
                // }
            }
            else if (High[det].isValid) /// High whitout Low
            {
                H_Channel_SiPM_High_Alone[NUCLEUS][100+det]->Fill(High[det].Channel);
                if (High[det].Channel < 800e3)
                {
                    H_SiPM[NUCLEUS][peak_number][det]->Fill(High[det].Channel);
                    // SiPM.push_back(High[det]);
                    SiPMi = High[det];
                    Tree_Peaks[NUCLEUS][peak_number][GetDetectorChannel(SiPMi.Label)]->Fill();
                }
                else
                {
                    // Mean value of other low SiPMs
                    double sum = 0;
                    int counter = 0;
                    for (int i = 1; i <= 9; i++)
                    {
                        if (Low[i].isValid)
                        {
                            counter++;
                            sum += Low[i].Channel;
                        }
                    }
                    if (counter != 0)
                    {
                        High[det].Channel = sum / counter;
                        H_SiPM[NUCLEUS][peak_number][det]->Fill( sum/counter );
                        // SiPM.push_back(High[det]);
                        SiPMi = High[det];
                        Tree_Peaks[NUCLEUS][peak_number][GetDetectorChannel(SiPMi.Label)]->Fill();
                    }
                }
            }
            else if (Low[det].isValid) /// Low whitout High
            {
                H_Channel_SiPM_Low_Alone[NUCLEUS][110+det]->Fill(Low[det].Channel);
                if (Low[det].Channel > 800e3)
                {
                    H_SiPM[NUCLEUS][peak_number][det]->Fill(Low[det].Channel);
                    // SiPM.push_back(Low[det]);
                    SiPMi = Low[det];
                    Tree_Peaks[NUCLEUS][peak_number][GetDetectorChannel(SiPMi.Label)]->Fill();
                }
            }
        }
        
        // SiPM.clear();
    }

    // for general test display result on IAS
    for (int det = 1; det <= 9; det++)
    {
        H_SiPM_Merged[NUCLEUS][100+det] = (TH1D*)H_SiPM[NUCLEUS][14][det]->Clone();
        Tree_Peaks[NUCLEUS][14][det]->Write();
    }
    
}

void MergingSiPMsBi207()
{
    Info("Matching SiPMs for 207Bi");

    clock_t start = clock(), Current;
    Tree = (TTree *)GROUPED_File[NUCLEUS]->Get("Tree_Group");
    Reader = new TTreeReader(Tree);
    signals = new TTreeReaderArray<Signal>(*Reader, "Signal");
    Reader->Restart();
    while (Reader->Next())
    {
        if (Reader->GetCurrentEntry() < Reader->GetEntries() - 10000000)
            continue;
        ProgressBar(Reader->GetCurrentEntry(), Reader->GetEntries(), start, Current, "Reading Tree");
        int peak_number = 0;

        vector<Signal> SiPM_High;
        vector<Signal> SiPM_Low;
        for (int i = 0; i < signals->GetSize(); i++)
        {
            if (IsDetectorBetaHigh((*signals)[i].Label))
            {
                SiPM_High.push_back((*signals)[i]);
            }
            if (IsDetectorBetaLow((*signals)[i].Label))
            {
                SiPM_Low.push_back((*signals)[i]);
            }
        }

        Signal High[10] = {Signal()};
        Signal Low[10] = {Signal()};

        for (int i = 0; i < SiPM_High.size(); i++)
        {
            SiPM_High[i].Channel = MatchingSiPM["207Bi"][100 + GetDetectorChannel(SiPM_High[i].Label)]->Eval(SiPM_High[i].Channel);
            High[GetDetectorChannel(SiPM_High[i].Label)] = SiPM_High[i];
            H_SiPM_Matched[NUCLEUS][SiPM_High[i].Label]->Fill(SiPM_High[i].Channel);
        }
        for (int i = 0; i < SiPM_Low.size(); i++)
        {
            SiPM_Low[i].Channel = MatchingSiPM["207Bi"][100 + GetDetectorChannel(SiPM_Low[i].Label)]->Eval(MatchingLowHigh["32Ar"][100 + GetDetectorChannel(SiPM_Low[i].Label)]->Eval(SiPM_Low[i].Channel));
            Low[GetDetectorChannel(SiPM_Low[i].Label)] = SiPM_Low[i];
            H_SiPM_Matched[NUCLEUS][SiPM_Low[i].Label]->Fill(SiPM_Low[i].Channel);
        }

        for (int det = 1; det <= 9; det++)
        {
            
            if (High[det].isValid) /// High whitout Low
            {
                H_SiPM[NUCLEUS][peak_number][det]->Fill(High[det].Channel);
                // SiPM.push_back(High[det]);
                SiPMi = High[det];
                Tree_Peaks[NUCLEUS][0][GetDetectorChannel(SiPMi.Label)]->Fill();
            }
        }
        
        SiPM.clear();
    }

    // for general test display result on IAS
    for (int det = 1; det <= 9; det++)
    {
        H_SiPM_Merged[NUCLEUS][100 + det] = (TH1D *)H_SiPM[NUCLEUS][0][det]->Clone();
        Tree_Peaks[NUCLEUS][0][det]->Write();
    }
    
}

double Chi2TreeHist_conv(const double *par)
{
    for (string NUCLEUS : Nucleis)
    {
        if (NUCLEUS == "32Ar")
            continue;
        ////////////////////////////////////////////////////////////////
        // Init Histograms
        for (int peak = 0; peak <= 50; peak++)
        {
            if ((WindowsMap[NUCLEUS][peak][11].first == -1 || !WindowsMap[NUCLEUS][peak][11].first))
                continue;
            
            H_Sim_Conv[NUCLEUS][peak][current_detector]->Reset();


            H_SiPM_Calibrated[NUCLEUS][peak][current_detector]->Reset();
            
        }

        if (NUCLEUS == "207Bi")
        {
            H_Sim_Conv[NUCLEUS][0][current_detector]->Reset();


            H_SiPM_Calibrated[NUCLEUS][0][current_detector]->Reset();
            
        }
        ////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////
        // parameters
        double Calibration_OffSet = par[0];
        double Calibration = par[1];
        double Resolution_OffSet = par[2];
        double Resolution_SQRT = par[3];
        double Resolution_2 = par[4];
        double Threshold = par[5];
        double Threshold_STD = par[6];
        ////////////////////////////////////////////////////////////////

        cout << "calibration exp" << endl;
        // calibration exp

        int peak_number = 14;
        if (NUCLEUS == "207Bi")
        {
            peak_number = 0;
        }

        TTreeReader *Reader_IAS = new TTreeReader(Tree_Peaks[NUCLEUS][peak_number][current_detector]);
        TTreeReaderValue<Signal> *SiPM = new TTreeReaderValue<Signal>(*Reader_IAS, "SiPM");

        Reader_IAS->Restart();
        double energy;
        while (Reader_IAS->Next())
        {
            energy = Calibration * (**SiPM).Channel / 1000 + Calibration_OffSet;
            H_SiPM_Calibrated[NUCLEUS][peak_number][GetDetectorChannel((**SiPM).Label)]->Fill(energy);
        }
        
        cout << "convoluting resolution" << endl;
        // convoluting resolution on histogram
        // for (int i = 1; i <= H_Sim[NUCLEUS][peak_number]->GetNbinsX(); i++)
        // {
        //     if (H_Sim[NUCLEUS][peak_number]->GetBinContent(i) == 0)
        //         continue;
        //     energy = H_Sim[NUCLEUS][peak_number]->GetBinCenter(i);
        //     double sigma_resolution = sqrt(pow(Resolution_OffSet, 2) + pow(Resolution_SQRT * sqrt(energy), 2) + pow(Resolution_2 * pow(energy, 2), 2));
        //     // double sigma_resolution = sqrt(pow(coefficents[2].second, 2) + pow(coefficents[3].second * sqrt(energy), 2));

        //     gauss = new TF1("gauss", "gaus", 0, eSiliMax_cal);
        //     gauss->SetNpx(eSiliN_cal);
        //     gauss->SetParameters(H_Sim[NUCLEUS][peak_number]->GetBinContent(i) / (sigma_resolution * sqrt(2 * M_PI)), energy, sigma_resolution); /// weighted gaussian

        //     H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Add(gauss->GetHistogram());

        //     delete gauss;

        //     // for (int j = 1; j <= H_Sim_Conv[NUCLEUS][peak_number]->GetNbinsX(); j++)
        //     // {
        //     //     H_Sim_Conv[NUCLEUS][peak_number]->SetBinContent(i, H_Sim_Conv[NUCLEUS][peak_number]->GetBinContent(i) + gauss->Eval(H_Sim[NUCLEUS][peak_number]->GetBinCenter(j)));
        //     // }
        //     // delete gauss;
        // }
        ///////////////:
        // random convoluting resolution on histogram
        for (int i = 0; i < H_Sim[NUCLEUS][peak_number]->GetEntries() && i < 10000000; i++)
        {
            energy = H_Sim[NUCLEUS][peak_number]->GetRandom();
            double sigma_resolution = sqrt(pow(Resolution_OffSet, 2) + pow(Resolution_SQRT * sqrt(energy), 2) + pow(Resolution_2 * pow(energy, 2), 2));
            // shoot in a gaussian with ramdom engine
            H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Fill(gRandom->Gaus(energy, sigma_resolution));
        }

        cout << "convoluting threshold" << endl;

        if (NUCLEUS == "207Bi")
        {
            TF1 *fexp = new TF1("fexp", "expo", 0, 10000);
            fexp->SetNpx(10000);

            fexp->SetParLimits(0, 0, 200);
            fexp->SetParameter(0, 13);
            fexp->SetParLimits(1, -0.1, 0);
            fexp->SetParameter(1, -0.05);

            H_SiPM_Calibrated[NUCLEUS][0][current_detector]->Fit(fexp, "QR", "", 80, 150);

            TH1D *h = (TH1D *)fexp->GetHistogram()->Clone();
            H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Add(h);
            delete h;
            delete fexp;
        }

        ////////// CONVOLUTING THRESHOLD//////////
        Threshold_f = new TF1("Threshold", "0.5*(1+erf((x-[0])/[1]))", 0, 8000);
        Threshold_f->SetParameters(Threshold, Threshold_STD);

        for (int i = 1; i <= H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetNbinsX(); i++)
        {
            H_Sim_Conv[NUCLEUS][peak_number][current_detector]->SetBinContent(i, H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetBinContent(i) * Threshold_f->Eval(H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetBinCenter(i)));
        }

        delete Threshold_f;

        //////////////// CHI2 ////////////////
        double chi2 = H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector]->Chi2Test(H_Sim_Conv[NUCLEUS][peak_number][current_detector], "UW CHI2/NDF");
        H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector]->GetXaxis()->SetRangeUser(400, 10000);
        H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetXaxis()->SetRangeUser(400, 10000);
        H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector]->Scale(H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Integral() / H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector]->Integral());

        cout << chi2 << "  " << Calibration_OffSet << " " << Calibration << " " << Resolution_OffSet << " " << Resolution_SQRT << " " << Resolution_2 << " " << Threshold << " " << Threshold_STD << endl;
    }

    return 0;
}

void CalibrationSiPM()
{
    res = {0, 8.17, 6.5, 5.0, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5}; //#2 8.17
    cout << "Calibration SiPM : " << NUCLEUS << endl;
    Minimizer *minimizer = Factory::CreateMinimizer("Minuit2", "Minuit2");
    ROOT::Math::Functor functor(&Chi2TreeHist_conv, 7);
    vector<double> Par = {30.3, 1.8, 0, res[current_detector], 0.e-05, 70, 15.0801};
    minimizer->SetFunction(functor);
    minimizer->SetLimitedVariable(0, "Calibration_OffSet", Par[0], 10, 0, 80);
    minimizer->SetLimitedVariable(1, "Calibration", Par[1], 0.1, 0.5, 6);
    minimizer->SetFixedVariable(2, "Resolution_OffSet", Par[2]);
    minimizer->SetLimitedVariable(3, "Resolution_SQRT", Par[3], 0.1, 2., 10);
    minimizer->SetLimitedVariable(4, "Resolution_2", Par[4], 0.1e-5, 1e-6, 8e-5);
    minimizer->SetLimitedVariable(5, "Threshold", Par[5], 5, 20, 100);
    minimizer->SetLimitedVariable(6, "Threshold_STD", Par[6], 10, 2, 20);
    // minimizer->SetMaxFunctionCalls(10000000);
    // minimizer->SetMaxIterations(10000000);
    const double *bestPar = Par.data();
    // minimizer->Minimize();
    double chi2 = Chi2TreeHist_conv(bestPar);
}

void WriteHistograms()
{
    gStyle->SetOptStat(0);
    Info("Write Histograms");
    for (string NUCLEUS : Nucleis)
    {
        cout << "Writing Nucleus " << NUCLEUS << endl;
        dir[NUCLEUS]->cd();

        TCanvas *cHighLowSiPM = new TCanvas(("cHighLowSiPM_" + NUCLEUS).c_str(), ("cHighLowSiPM_" + NUCLEUS).c_str(), 800, 800);
        cHighLowSiPM->Divide(3, 3);

        TCanvas *cHighLowSiPM_Fitting = new TCanvas(("cHighLowSiPM_Fitting_" + NUCLEUS).c_str(), ("cHighLowSiPM_Fitting_" + NUCLEUS).c_str(), 800, 800);
        cHighLowSiPM_Fitting->Divide(3, 3);

        TCanvas *cHighLowSiPM_Diff = new TCanvas(("cHighLowSiPM_Diff_" + NUCLEUS).c_str(), ("cHighLowSiPM_Diff_" + NUCLEUS).c_str(), 800, 800);
        cHighLowSiPM_Diff->Divide(3, 3);

        TCanvas *cMatchingSiPM_Fitting = new TCanvas(("cMatchingSiPM_Fiiting_" + NUCLEUS).c_str(), ("cMatchingSiPM_Fitting_" + NUCLEUS).c_str(), 800, 800);
        cMatchingSiPM_Fitting->Divide(3, 3);

        TCanvas *cMatchingSiPM = new TCanvas(("cMatchingSiPM_" + NUCLEUS).c_str(), ("cMatchingSiPM_" + NUCLEUS).c_str(), 800, 800);
        cMatchingSiPM->Divide(3, 3);

        TCanvas *cSiPM = new TCanvas(("cSiPM_" + NUCLEUS).c_str(), ("cSiPM_" + NUCLEUS).c_str(), 800, 800);
        TLegend *lSiPM = new TLegend(0.1, 0.7, 0.48, 0.9);
        TCanvas *cSiPM_UnMatched = new TCanvas(("cSiPM_UnMatched_" + NUCLEUS).c_str(), ("cSiPM_UnMatched_" + NUCLEUS).c_str(), 800, 800);
        TLegend *lSiPM_UnMatched = new TLegend(0.1, 0.7, 0.48, 0.9);
        TCanvas *cSiPM_LowHighTriggered = new TCanvas(("cSiPM_LowHighTriggered_" + NUCLEUS).c_str(), ("cSiPM_LowHighTriggered_" + NUCLEUS).c_str(), 800, 800);
        TLegend *lSiPM_LowHighTriggered = new TLegend(0.1, 0.7, 0.48, 0.9);
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorBetaHigh(i))
            {
                cout << "Writing SiPM " << i << endl;
                dir_detector[NUCLEUS][i]->cd();
                
                    
                    // single histograms
                    H_SiPM_High[NUCLEUS][i]->Write();
                    H_Channel_SiPM_High_Alone[NUCLEUS][i]->Write();
                    H_SiPM_HighLow[NUCLEUS][i]->Write();


                        // 3x3 HighLow graph fitted for matching
                        cHighLowSiPM_Fitting->cd(GetDetectorChannel(i));
                        G_SiPM_HighLow[NUCLEUS][i]->SetTitle(("SiPM " + to_string(GetDetectorChannel(i))).c_str());
                        G_SiPM_HighLow[NUCLEUS][i]->GetXaxis()->SetTitle("Low SiPM");
                        G_SiPM_HighLow[NUCLEUS][i]->GetYaxis()->SetTitle("High SiPM");
                        G_SiPM_HighLow[NUCLEUS][i]->GetXaxis()->CenterTitle();
                        G_SiPM_HighLow[NUCLEUS][i]->GetYaxis()->CenterTitle();
                        G_SiPM_HighLow[NUCLEUS][i]->Draw("AP");
                        if (NUCLEUS == "32Ar")
                        {
                            TLatex *text = new TLatex();
                            text->SetNDC();
                            text->SetTextSize(0.1);
                            std::ostringstream streamObj;
                            streamObj << std::fixed << std::setprecision(2) << MatchingLowHigh[NUCLEUS][i]->GetParameter(0);
                            std::string parameterStr = streamObj.str();
                            text->DrawLatex(0.7, 0.8, ("#color[2]{" + parameterStr + "}").c_str());
                            text->Draw("SAME");
                        }

                        // 3x3 HighLow histograms projected on diff
                        cHighLowSiPM_Diff->cd(GetDetectorChannel(i));
                        // H_SiPM_HighLow_Matched_diff[NUCLEUS][i]->Draw("COLZ");

                        // 3x3 HighLow histograms matched
                        cHighLowSiPM->cd(GetDetectorChannel(i));
                        H_SiPM_HighLow_Matched[NUCLEUS][i]->Rebin(10);
                        H_SiPM_HighLow_Matched[NUCLEUS][i]->Draw("HIST");
                        H_SiPM_HighLow_Matched[NUCLEUS][i + 10]->Rebin(10);
                        H_SiPM_HighLow_Matched[NUCLEUS][i + 10]->SetLineColor(kRed);
                        H_SiPM_HighLow_Matched[NUCLEUS][i + 10]->Draw("HIST SAME");

                        // 3x3 SIPM graph fitted for matching
                        cMatchingSiPM_Fitting->cd(GetDetectorChannel(i));
                        G_SiPM_HighLow[NUCLEUS][i]->SetTitle(("SiPM " + to_string(GetDetectorChannel(i))).c_str());
                        G_SiPM_HighLow[NUCLEUS][i]->GetXaxis()->SetTitle("Low SiPM");
                        G_SiPM_HighLow[NUCLEUS][i]->GetYaxis()->SetTitle("High SiPM");
                        G_SiPM_HighLow[NUCLEUS][i]->GetXaxis()->CenterTitle();
                        G_SiPM_HighLow[NUCLEUS][i]->GetYaxis()->CenterTitle();
                        G_Matching_SiPM[NUCLEUS][GetDetectorChannel(i)]->Draw("AP");
                        if (NUCLEUS == "32Ar")
                        {
                            TLatex *textSIPM = new TLatex();
                            textSIPM->SetNDC();
                            textSIPM->SetTextSize(0.1);
                            std::ostringstream streamObjSipm;
                            streamObjSipm << std::fixed << std::setprecision(2) << MatchingSiPM[NUCLEUS][i]->GetParameter(0);
                            std::string parameterStrsipm = streamObjSipm.str();
                            textSIPM->DrawLatex(0.7, 0.8, ("#color[2]{" + parameterStrsipm + "}").c_str());
                            textSIPM->Draw("SAME");
                        }

                        // 3x3 SiPM histograms matched
                        cMatchingSiPM->cd(GetDetectorChannel(i));
                        // H_SiPM_Matched[NUCLEUS][i]->Rebin(10);
                        H_SiPM_Matched[NUCLEUS][i]->Draw("HIST");
                        // H_SiPM_Matched[NUCLEUS][i + 10]->Rebin(10);
                        H_SiPM_Matched[NUCLEUS][i + 10]->SetLineColor(kRed);
                        H_SiPM_Matched[NUCLEUS][i + 10]->Draw("HIST SAME");
                        // H_SiPM_Merged[NUCLEUS][i]->Rebin(10);
                        H_SiPM_Merged[NUCLEUS][i]->SetLineColor(kBlack);
                        H_SiPM_Merged[NUCLEUS][i]->Draw("HIST SAME");
                    
                    // superimposed SiPM matched merged sipms
                    cSiPM->cd();
                    H_SiPM_Merged[NUCLEUS][i]->Draw("HIST SAME");
                
                // Losses
                TCanvas *cLosses_High = new TCanvas(("cLosses_"+NUCLEUS+"_High_"+to_string(GetDetectorChannel(i))).c_str(), ("cLosses_"+NUCLEUS+"_High_"+to_string(GetDetectorChannel(i))).c_str(), 800, 800);
                TPad *pad_hist = new TPad("pad_hist", "pad_hist", 0, 0.3, 1, 1);
                TPad *pad_ratio = new TPad("pad_ratio", "pad_ratio", 0, 0, 1, 0.3);
                pad_hist->Draw();
                pad_ratio->Draw();
                pad_hist->cd();
                H_SiPM_Matched[NUCLEUS][i]->Draw("HIST");
                H_Channel_SiPM_High_Alone[NUCLEUS][i]->SetLineColor(kBlack);
                H_Channel_SiPM_High_Alone[NUCLEUS][i]->Draw("HIST SAME");
                TPaveText *pt = new TPaveText(0.7, 0.7, 0.9, 0.9);
                double integral1 = H_SiPM_Matched[NUCLEUS][i]->Integral();
                double integral2 = H_Channel_SiPM_High_Alone[NUCLEUS][i]->Integral();
                pt->AddText(("Losses : " + to_string(integral2/integral1*100) + " %").c_str());
                pt->Draw("SAME");
                pad_ratio->cd();
                H_Channel_SiPM_High_Alone[NUCLEUS][i]->Rebin(10);   
                TH1D *H = (TH1D*)H_SiPM_Matched[NUCLEUS][i]->Clone();
                H->Reset(); 
                for (int bin = 0; bin < H->GetNbinsX(); bin++)
                {
                    double value1 = H_SiPM_Matched[NUCLEUS][i]->GetBinContent(bin);
                    double value2 = H_Channel_SiPM_High_Alone[NUCLEUS][i]->GetBinContent(bin);
                    if (value1 == 0 || value2 == 0)
                        continue;
                    H->SetBinContent(bin, value2/value1);
                }
                H->GetYaxis()->SetRangeUser(0, 1);
                H->GetXaxis()->SetTitle("Channel");
                H->GetYaxis()->SetTitle("Loss probability");
                H->Draw("E0");
                cLosses_High->Write();
                
        
            }

            if (IsDetectorBetaLow(i))
            {
                dir_detector[NUCLEUS][i - 10]->cd();
                H_SiPM_Low[NUCLEUS][i]->Write();
                H_Channel_SiPM_Low_Alone[NUCLEUS][i]->Write();
            
                TCanvas *cLosses_Low = new TCanvas(("cLosses_"+NUCLEUS+"_Low_"+to_string(GetDetectorChannel(i))).c_str(), ("cLosses_Low_"+NUCLEUS+"_"+to_string(GetDetectorChannel(i))).c_str(), 800, 800);
                TPad *pad_hist = new TPad("pad_hist", "pad_hist", 0, 0.3, 1, 1);
                TPad *pad_ratio = new TPad("pad_ratio", "pad_ratio", 0, 0, 1, 0.3);
                pad_hist->Draw();
                pad_ratio->Draw();
                pad_hist->cd();
                H_SiPM_Matched[NUCLEUS][i]->Draw("HIST");
                H_Channel_SiPM_Low_Alone[NUCLEUS][i]->SetLineColor(kBlack);
                H_Channel_SiPM_Low_Alone[NUCLEUS][i]->Draw("HIST SAME");
                TPaveText *pt = new TPaveText(0.7, 0.7, 0.9, 0.9);
                double integral1 = H_SiPM_Matched[NUCLEUS][i]->Integral();  
                double integral2 = H_Channel_SiPM_Low_Alone[NUCLEUS][i]->Integral();
                pt->AddText(("Losses : " + to_string(integral2/integral1*100) + " %").c_str());
                pt->Draw("SAME");
                pad_ratio->cd();
                H_Channel_SiPM_Low_Alone[NUCLEUS][i]->Rebin(10);
                TH1D *H = (TH1D*)H_SiPM_Matched[NUCLEUS][i]->Clone();
                H->Reset();
                for (int bin = 0; bin < H->GetNbinsX(); bin++)
                {
                    double value1 = H_SiPM_Matched[NUCLEUS][i]->GetBinContent(bin);
                    double value2 = H_Channel_SiPM_Low_Alone[NUCLEUS][i]->GetBinContent(bin);
                    if (value1 == 0 || value2 == 0)
                        continue;
                    H->SetBinContent(bin, value2/value1);
                }
                H->GetYaxis()->SetRangeUser(0, 1);
                H->GetXaxis()->SetTitle("Channel");
                H->GetYaxis()->SetTitle("Loss probability");
                H->Draw("E1");
                cLosses_Low->Write();

            }
        }

        // for (int i = 0; i <= 5; i++)
        // {
        //     cSiPM_LowHighTriggered->cd();
        //     H_SiPM_LowHighTriggered_x[NUCLEUS][i]->SetLineColor(i);
        //     lSiPM_LowHighTriggered->AddEntry(H_SiPM_LowHighTriggered_x[NUCLEUS][i], ("Triggered x" + to_string(i)).c_str(), "l");
        //     H_SiPM_LowHighTriggered_x[NUCLEUS][i]->Draw("SAME");
        // }

        dir[NUCLEUS]->cd();
        // H_SiPM_HighTriggered_x[NUCLEUS]->Write();
        // H_SiPM_LowTriggered_x[NUCLEUS]->Write();

        cHighLowSiPM->Write();
        cHighLowSiPM_Fitting->Write();
        cHighLowSiPM_Diff->Write();

        cMatchingSiPM->Write();
        cMatchingSiPM_Fitting->Write();

        cSiPM->cd();
        // lSiPM->Draw("SAME");
        cSiPM->Write();

        // if (NUCLEUS == "32Ar")
        // {
        //     for (int det = 1; det <= 9; det++)
        //     {
        //         for (int peak = 1; peak <= 50; peak++)
        //         {
        //             if (WindowsMap[NUCLEUS][peak][11].first == -1 || !WindowsMap[NUCLEUS][peak][11].first)
        //                 continue;
        //             dir_SiPM[NUCLEUS][det]->cd();
        //             TCanvas *cExp_Sim1 = new TCanvas(("cExp_Sim1_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(det)).c_str(), ("cExp_Sim1_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(det)).c_str(), 800, 800);
        //             H_SiPM_Calibrated[NUCLEUS][peak][det]->Draw("HIST");
        //             // H_Sim[NUCLEUS][peak]->SetLineColor(kRed);
        //             // H_Sim[NUCLEUS][peak]->Draw("HIST SAME");
        //             H_Sim_Conv[NUCLEUS][peak]->SetLineColor(kRed);
        //             H_Sim_Conv[NUCLEUS][peak]->Draw("HIST SAME");
        //             cExp_Sim1->Write();
        //         }
        //         MatchingSiPM[NUCLEUS][100 + det]->SetName(("MatchingSiPM_" + NUCLEUS + "_SiPM" + to_string(det)).c_str());
        //         MatchingSiPM[NUCLEUS][100 + det]->Write();
        //     }
        // }
    }

    NUCLEUS="207Bi";
    dir["207Bi"]->cd();
    TCanvas *cMatchingSiPM = new TCanvas(("cMatchingSiPM_" + NUCLEUS).c_str(), ("cMatchingSiPM_" + NUCLEUS).c_str(), 800, 800);
    cMatchingSiPM->Divide(3, 3);
    TLegend *lSiPM = new TLegend(0.1, 0.7, 0.48, 0.9);
    for (int det = 1; det <= 9; det++)
    {
        cMatchingSiPM->cd(det);
        cMatchingSiPM->cd(det)->SetLogy();
        H_SiPM_Merged["207Bi"][100 + det]->GetXaxis()->SetRangeUser(0, 1000000);
        H_SiPM_Merged["207Bi"][100 + det]->SetLineColor(det);
        H_SiPM_Merged["207Bi"][100 + det]->Draw("HIST");
        lSiPM->AddEntry(H_SiPM_Merged["207Bi"][100 + det], ("SiPM " + to_string(det)).c_str(), "l");

        // dir_SiPM["207Bi"][det]->cd();
    }
    lSiPM->Draw("SAME");
    dir["207Bi"]->cd();
    cMatchingSiPM->Write();

    for (int det = 1; det <= 9; det++)
    {
        int peak = 0;
        dir_SiPM[NUCLEUS][det]->cd();
        TCanvas *cExp_Sim1 = new TCanvas(("cExp_Sim1_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(det)).c_str(), ("cExp_Sim1_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(det)).c_str(), 800, 800);
        H_SiPM_Calibrated[NUCLEUS][peak][det]->Draw("HIST");
        // H_Sim[NUCLEUS][peak]->SetLineColor(kRed);
        // H_Sim[NUCLEUS][peak]->Draw("HIST SAME");
        H_Sim_Conv[NUCLEUS][peak][det]->SetLineColor(kRed);
        H_Sim_Conv[NUCLEUS][peak][det]->Draw("HIST SAME");
        cExp_Sim1->Write();
    }
}

void InitSiliconCalibration()
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


void SiPM_Calibration()
{

}
#endif