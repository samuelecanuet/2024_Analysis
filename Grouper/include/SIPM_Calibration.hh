#ifndef SIPM_CALIBRATION_HH
#define SIPM_CALIBRATION_HH

#include "../../../lib/SignalDict/Signal.h"
#include "Detectors.hh"

/// FILE ///
map<string, TFile *> GROUPED_File;
TFile *CALIBRATED_File;

TTree *Tree;
TTreeReader *Reader;
TTreeReaderArray<Signal> *SiPM_High;
TTreeReaderArray<Signal> *SiPM_Low;
TTreeReaderArray<Signal> *Silicon;

TTree *Tree_MATCHED;
vector<Signal> SiPM;

/// WINDOWS ///
map<string, pair<double, double>[100][SIGNAL_MAX]> WindowsMap;

// TREE PEAKS //
map<int, TTree *> Tree_Peaks;
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

map<string, TH2D *> H_SiPM_coinc;
/// DIRECTORY ///
map<string, TDirectory *> dir;
map<string, TDirectory *[SIGNAL_MAX]> dir_detector;
map<string, TDirectory *[SIGNAL_MAX]> dir_SiPM;

// FUCNTIONS //
TF1 *F_SiliconCalibration[SIGNAL_MAX];
TF1 *f_linear = new TF1("f_linear", "[0]*x + [1]", 0, 10000e3);
TF1 *f_erf = new TF1("f_erf", "x < [0] ? [1]*x + [2] : [3]*(erf((x-[4])/[5]))", 0, 10000e3);
pair<double, double> Range_SiPM_LowHigh = make_pair(50e3, 200e3);

double SiPM_Window[10] = {0, 1500e3, 1500e3, 1500e3};

// DATA //
string Nuclei[1] = {"32Ar"}; //, "32Ar_thick", "33Ar"};
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
                    Tree_Peaks[number] = new TTree(("Tree_Peaks_" + to_string(number)).c_str(), ("Tree_Peaks_" + to_string(number)).c_str());
                    Tree_Peaks[number]->Branch("SiPM", &SiPM);
                }
            }
        }
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

void InitHistograms()
{
    Info("Init Histograms");
    for (string Nucleus : Nuclei)
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

                dir_detector[NUCLEUS][i] = dir[NUCLEUS]->mkdir(("SiPM_" + to_string(GetDetectorChannel(i))).c_str());
                H_SiPM_HighLow[NUCLEUS][i] = new TH2D(("H_SiPM_HighLow_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_SiPM_HighLow_" + NUCLEUS + "_" + detectorName[i]).c_str(), eLowN, eLowMin, eLowMax, eHighN, eHighMin, eHighMax);
                H_SiPM_HighLow[NUCLEUS][i]->GetXaxis()->SetTitle("Channel Low");
                H_SiPM_HighLow[NUCLEUS][i]->GetYaxis()->SetTitle("Channel High");
                H_SiPM_HighLow[NUCLEUS][i]->GetXaxis()->CenterTitle();
                H_SiPM_HighLow[NUCLEUS][i]->GetYaxis()->CenterTitle();

                G_SiPM_HighLow[NUCLEUS][i] = new TGraph();
                G_SiPM_HighLow_Projected[NUCLEUS][i] = new TGraph();
                G_SiPM_HighLow_Matched[NUCLEUS][i] = new TGraph();

                H_SiPM_HighLow_Projected[NUCLEUS][i] = new TH1D(("H_SiPM_HighLow_Projected_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_SiPM_HighLow_Projected_" + NUCLEUS + "_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
                H_SiPM_HighLow_Projected[NUCLEUS][i]->GetXaxis()->SetTitle("Channel");
                H_SiPM_HighLow_Projected[NUCLEUS][i]->GetYaxis()->SetTitle("Counts");
                H_SiPM_HighLow_Projected[NUCLEUS][i]->GetXaxis()->CenterTitle();
                H_SiPM_HighLow_Projected[NUCLEUS][i]->GetYaxis()->CenterTitle();

                H_SiPM_HighLow_Matched_diff[NUCLEUS][i] = new TH2D(("H_SiPM_HighLow_Matched_diff_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_SiPM_HighLow_Matched_diff_" + NUCLEUS + "_" + detectorName[i]).c_str(), 2000, -100000, 100000, eHighN, eHighMin, eHighMax);
                H_SiPM_HighLow_Matched_diff[NUCLEUS][i]->GetXaxis()->SetTitle("Channel Low - Channel High");
                H_SiPM_HighLow_Matched_diff[NUCLEUS][i]->GetYaxis()->SetTitle("Channel High");
                H_SiPM_HighLow_Matched_diff[NUCLEUS][i]->GetXaxis()->CenterTitle();
                H_SiPM_HighLow_Matched_diff[NUCLEUS][i]->GetYaxis()->CenterTitle();

                G_Matching_SiPM[NUCLEUS][GetDetectorChannel(i)] = new TGraph();
                

                H_SiPM_HighLow_Matched[NUCLEUS][i] = new TH1D(("H_SiPM_HighLow_Matched_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_SiPM_HighLow_Matched_" + NUCLEUS + "_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
                H_SiPM_HighLow_Matched[NUCLEUS][i]->GetXaxis()->SetTitle("Channel");
                H_SiPM_HighLow_Matched[NUCLEUS][i]->GetYaxis()->SetTitle("Counts");
                H_SiPM_HighLow_Matched[NUCLEUS][i]->GetXaxis()->CenterTitle();
                H_SiPM_HighLow_Matched[NUCLEUS][i]->GetYaxis()->CenterTitle();

                H_SiPM_Merged[NUCLEUS][i] = new TH1D(("H_SiPM_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_SiPM_" + NUCLEUS + "_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
                H_SiPM_Merged[NUCLEUS][i]->GetXaxis()->SetTitle("Channel");
                H_SiPM_Merged[NUCLEUS][i]->GetYaxis()->SetTitle("Counts");
                H_SiPM_Merged[NUCLEUS][i]->GetXaxis()->CenterTitle();
                H_SiPM_Merged[NUCLEUS][i]->GetYaxis()->CenterTitle();

                H_SiPM_Matched[NUCLEUS][i] = new TH1D(("H_SiPM_Matched_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_SiPM_Matched_" + NUCLEUS + "_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
                H_SiPM_Matched[NUCLEUS][i]->GetXaxis()->SetTitle("Channel");
                H_SiPM_Matched[NUCLEUS][i]->GetYaxis()->SetTitle("Counts");
                H_SiPM_Matched[NUCLEUS][i]->GetXaxis()->CenterTitle();
                H_SiPM_Matched[NUCLEUS][i]->GetYaxis()->CenterTitle();

                H_SiPM_UnMatched[NUCLEUS][i] = new TH1D(("H_SiPM_UnMatched_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_SiPM_UnMatched_" + NUCLEUS + "_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
                H_SiPM_UnMatched[NUCLEUS][i]->GetXaxis()->SetTitle("Channel");
                H_SiPM_UnMatched[NUCLEUS][i]->GetYaxis()->SetTitle("Counts");
                H_SiPM_UnMatched[NUCLEUS][i]->GetXaxis()->CenterTitle();
                H_SiPM_UnMatched[NUCLEUS][i]->GetYaxis()->CenterTitle();

                H_Channel_SiPM_High_Alone[NUCLEUS][i] = new TH1D(("H_Channel_SiPM_High_Alone_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_Channel_SiPM_High_Alone_" + NUCLEUS + "_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
                H_Channel_SiPM_High_Alone[NUCLEUS][i]->GetXaxis()->SetTitle("Channel");
                H_Channel_SiPM_High_Alone[NUCLEUS][i]->GetYaxis()->SetTitle("Counts");
                H_Channel_SiPM_High_Alone[NUCLEUS][i]->GetXaxis()->CenterTitle();
                H_Channel_SiPM_High_Alone[NUCLEUS][i]->GetYaxis()->CenterTitle();

                dir_SiPM[NUCLEUS][GetDetectorChannel(i)] = CALIBRATED_File->mkdir(("SiPM_" + to_string(GetDetectorChannel(i))).c_str());
                for (int peak = 1; peak <= 50; peak++)
                {
                    if (WindowsMap[NUCLEUS][peak][11].first == -1 || !WindowsMap[NUCLEUS][peak][11].first)
                        continue;
                    
                    H_SiPM[NUCLEUS][peak][GetDetectorChannel(i)] = new TH1D(("H_SiPM_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), ("H_SiPM_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), eHighN, eHighMin, eHighMax);
                    H_SiPM[NUCLEUS][peak][GetDetectorChannel(i)]->GetXaxis()->SetTitle("Channel");
                    H_SiPM[NUCLEUS][peak][GetDetectorChannel(i)]->GetYaxis()->SetTitle("Counts");
                    H_SiPM[NUCLEUS][peak][GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
                    H_SiPM[NUCLEUS][peak][GetDetectorChannel(i)]->GetYaxis()->CenterTitle();
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

        H_SiPM_coinc[NUCLEUS] = new TH2D(("H_SiPM_coinc_" + NUCLEUS).c_str(), ("H_SiPM_coinc_" + NUCLEUS).c_str(), 20, 0, 20, 20, 0, 20);
        H_SiPM_coinc[NUCLEUS]->GetXaxis()->SetTitle("Low");
        H_SiPM_coinc[NUCLEUS]->GetYaxis()->SetTitle("High");
        H_SiPM_coinc[NUCLEUS]->GetXaxis()->CenterTitle();
        H_SiPM_coinc[NUCLEUS]->GetYaxis()->CenterTitle();

        for (int i = 0; i <= 5; i++)
        {
            H_SiPM_LowHighTriggered_x[NUCLEUS][i] = new TH1D(("H_SiPM_LowHighTriggered_" + to_string(i)).c_str(), ("H_SiPM_LowHighTriggered_x_" + to_string(i)).c_str(), 9, 0.5, 9.5);
            H_SiPM_LowHighTriggered_x[NUCLEUS][i]->GetXaxis()->SetTitle("SiPM");
            H_SiPM_LowHighTriggered_x[NUCLEUS][i]->GetYaxis()->SetTitle("Counts");
            H_SiPM_LowHighTriggered_x[NUCLEUS][i]->GetXaxis()->CenterTitle();
            H_SiPM_LowHighTriggered_x[NUCLEUS][i]->GetYaxis()->CenterTitle();
        }

        H_SiPM_HighTriggered_x[NUCLEUS] = new TH1D("H_SiPM_HighTriggered_x", "H_SiPM_HighTriggered_x", 10, 0.5, 10.5);
        H_SiPM_HighTriggered_x[NUCLEUS]->GetXaxis()->SetTitle("SiPM");
        H_SiPM_HighTriggered_x[NUCLEUS]->GetYaxis()->SetTitle("Counts");
        H_SiPM_HighTriggered_x[NUCLEUS]->GetXaxis()->CenterTitle();
        H_SiPM_HighTriggered_x[NUCLEUS]->GetYaxis()->CenterTitle();

        H_SiPM_LowTriggered_x[NUCLEUS] = new TH1D("H_SiPM_LowTriggered_x", "H_SiPM_LowTriggered_x", 10, 0.5, 10.5);
        H_SiPM_LowTriggered_x[NUCLEUS]->GetXaxis()->SetTitle("SiPM");
        H_SiPM_LowTriggered_x[NUCLEUS]->GetYaxis()->SetTitle("Counts");
        H_SiPM_LowTriggered_x[NUCLEUS]->GetXaxis()->CenterTitle();
        H_SiPM_LowTriggered_x[NUCLEUS]->GetYaxis()->CenterTitle();
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

            f_erf->SetParameter(0, 225e3);
            f_erf->SetParLimits(0, 150e3, 300e3); 
            f_erf->SetParameter(1, 10);
            f_erf->SetParLimits(1, 9.5, 12);
            f_erf->SetParameter(2, -50000);
            f_erf->SetParLimits(2, -100000, 0);
            f_erf->SetParameter(3, 5000e3);
            f_erf->SetParLimits(3, 4000e3, 6000e3);
            f_erf->SetParameter(4, 0);
            f_erf->SetParLimits(4, -100000, 100000);
            f_erf->SetParameter(5, 280e3);
            f_erf->SetParLimits(5, 0, 600e3);
            G_SiPM_HighLow[NUCLEUS][i]->Fit(f_erf);
            MatchingLowHigh_erf[NUCLEUS][i] = G_SiPM_HighLow[NUCLEUS][i]->GetFunction("f_erf");

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
            H_SiPM_HighLow_Matched[NUCLEUS][(*SiPM_Low)[j].Label]->Fill(MatchingLowHigh[NUCLEUS][(*SiPM_Low)[j].Label - 10]->Eval((*SiPM_Low)[j].Channel));
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

                    H_SiPM_HighLow_Matched_diff[NUCLEUS][(*SiPM_High)[i].Label]->Fill(MatchingLowHigh_erf[NUCLEUS][(*SiPM_High)[i].Label]->Eval((*SiPM_Low)[j].Channel)-high, (*SiPM_Low)[j].Channel);
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
            f_linear->FixParameter(1, 0);
            G_Matching_SiPM[NUCLEUS][GetDetectorChannel(i)]->Fit(f_linear, "QR", "", 1e6, 2 * eHighMax);
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
            H_SiPM_Matched[NUCLEUS][(*SiPM_High)[i].Label]->Fill((*SiPM_High)[i].Channel, 0.5);
        }
        for (int i = 0; i < SiPM_Low->GetSize(); i++)
        {
            (*SiPM_Low)[i].Channel = MatchingSiPM[NUCLEUS][100 + GetDetectorChannel((*SiPM_Low)[i].Label)]->Eval(MatchingLowHigh[NUCLEUS][100 + GetDetectorChannel((*SiPM_Low)[i].Label)]->Eval((*SiPM_Low)[i].Channel));
            Low[GetDetectorChannel((*SiPM_Low)[i].Label)] = (*SiPM_Low)[i];
            H_SiPM_Matched[NUCLEUS][(*SiPM_Low)[i].Label]->Fill((*SiPM_Low)[i].Channel, 0.5);
        }

        for (int det = 1; det <= 9; det++)
        {
            
            if (High[det].isValid && Low[det].isValid) /// BOTH
            {
                if (High[det].Channel < 750e3)
                {
                    H_SiPM[NUCLEUS][peak_number][det]->Fill(High[det].Channel);
                    SiPM.push_back(High[det]);
                }
                else if (Low[det].Channel > 850e3)
                {
                    H_SiPM[NUCLEUS][peak_number][det]->Fill(Low[det].Channel);
                    SiPM.push_back(Low[det]);
                }
                // / make the ratio betwwen them in the range 750 and 850 to get smooth transition
                else
                {
                    double ratio = (Low[det].Channel - 750e3) / (850e3 - 750e3);
                    H_SiPM[NUCLEUS][peak_number][det]->Fill(High[det].Channel * (1 - ratio) + Low[det].Channel * ratio);
                }
            }
            else if (High[det].isValid) /// High whitout Low
            {
                H_Channel_SiPM_High_Alone[NUCLEUS][100+det]->Fill(High[det].Channel);
                if (High[det].Channel < 800e3)
                {
                    H_SiPM[NUCLEUS][peak_number][det]->Fill(High[det].Channel);
                    SiPM.push_back(High[det]);
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
                        SiPM.push_back(High[det]);
                    }
                }
            }
            else if (Low[det].isValid) /// Low whitout High
            {
                H_Channel_SiPM_Low_Alone[NUCLEUS][110+det]->Fill(Low[det].Channel);
                if (Low[det].Channel > 800e3)
                {
                    H_SiPM[NUCLEUS][peak_number][det]->Fill(Low[det].Channel);
                    SiPM.push_back(Low[det]);
                }
            }
        }
        Tree_Peaks[peak_number]->Fill();
        SiPM.clear();
    }

    // for general test display result on IAS
    for (int det = 1; det <= 9; det++)
    {
        H_SiPM_Merged[NUCLEUS][100+det] = (TH1D*)H_SiPM[NUCLEUS][14][det]->Clone();
    }
    
}


void WriteHistograms()
{
    Info("Write Histograms");
    for (string NUCLEUS : Nuclei)
    {
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
                dir_detector[NUCLEUS][i]->cd();
                // single histograms
                H_SiPM_High[NUCLEUS][i]->Write();
                H_Channel_SiPM_High_Alone[NUCLEUS][i]->Write();
                H_SiPM_HighLow[NUCLEUS][i]->Write();

                // 3x3 HighLow graph fitted for matching
                cHighLowSiPM_Fitting->cd(GetDetectorChannel(i));
                G_SiPM_HighLow[NUCLEUS][i]->Draw("AP");

                // 3x3 HighLow histograms projected on diff
                cHighLowSiPM_Diff->cd(GetDetectorChannel(i));
                H_SiPM_HighLow_Matched_diff[NUCLEUS][i]->Draw("COLZ");

                // 3x3 HighLow histograms matched
                cHighLowSiPM->cd(GetDetectorChannel(i));
                H_SiPM_HighLow_Matched[NUCLEUS][i]->Rebin(10);
                H_SiPM_HighLow_Matched[NUCLEUS][i]->Draw("HIST");
                H_SiPM_HighLow_Matched[NUCLEUS][i + 10]->Rebin(10);
                H_SiPM_HighLow_Matched[NUCLEUS][i + 10]->SetLineColor(kRed);
                H_SiPM_HighLow_Matched[NUCLEUS][i + 10]->Draw("HIST SAME");

                // 3x3 SIPM graph fitted for matching
                cMatchingSiPM_Fitting->cd(GetDetectorChannel(i));
                G_Matching_SiPM[NUCLEUS][GetDetectorChannel(i)]->Draw("AP");

                // 3x3 SiPM histograms matched
                cMatchingSiPM->cd(GetDetectorChannel(i));
                H_SiPM_Matched[NUCLEUS][i]->Rebin(10);
                H_SiPM_Matched[NUCLEUS][i]->Draw("HIST");
                H_SiPM_Matched[NUCLEUS][i + 10]->Rebin(10);
                H_SiPM_Matched[NUCLEUS][i + 10]->SetLineColor(kRed);
                H_SiPM_Matched[NUCLEUS][i + 10]->Draw("HIST SAME");
                H_SiPM_Merged[NUCLEUS][i]->Rebin(10);
                H_SiPM_Merged[NUCLEUS][i]->SetLineColor(kGreen);
                H_SiPM_Merged[NUCLEUS][i]->Draw("HIST SAME");
            }

            if (IsDetectorBetaLow(i))
            {
                dir_detector[NUCLEUS][i - 10]->cd();
                H_SiPM_Low[NUCLEUS][i]->Write();
                H_Channel_SiPM_Low_Alone[NUCLEUS][i]->Write();
            }
        }

        for (int i = 0; i <= 5; i++)
        {
            cSiPM_LowHighTriggered->cd();
            H_SiPM_LowHighTriggered_x[NUCLEUS][i]->SetLineColor(i);
            lSiPM_LowHighTriggered->AddEntry(H_SiPM_LowHighTriggered_x[NUCLEUS][i], ("Triggered x" + to_string(i)).c_str(), "l");
            H_SiPM_LowHighTriggered_x[NUCLEUS][i]->Draw("SAME");
        }

        dir[NUCLEUS]->cd();
        H_SiPM_HighTriggered_x[NUCLEUS]->Write();
        H_SiPM_LowTriggered_x[NUCLEUS]->Write();

        cHighLowSiPM->Write();
        cHighLowSiPM_Fitting->Write();
        cHighLowSiPM_Diff->Write();

        cMatchingSiPM->Write();
        cMatchingSiPM_Fitting->Write();

        cSiPM->cd();
        lSiPM->Draw("SAME");
        cSiPM->Write();
        cSiPM_UnMatched->cd();
        lSiPM_UnMatched->Draw("SAME");
        cSiPM_UnMatched->Write();
        cSiPM_LowHighTriggered->cd();
        lSiPM_LowHighTriggered->Draw("SAME");
        cSiPM_LowHighTriggered->Write();

        
        for (int det = 1; det <= 9; det++)
        {
            
            for (int peak = 1; peak <= 50; peak++)
            {
                if (WindowsMap[NUCLEUS][peak][11].first == -1 || !WindowsMap[NUCLEUS][peak][11].first)
                    continue;
                
                dir_SiPM[NUCLEUS][det]->cd();
                H_SiPM[NUCLEUS][peak][det]->Write();
            }
            MatchingSiPM[NUCLEUS][100 + det]->SetName(("MatchingSiPM_" + NUCLEUS + "_SiPM" + to_string(det)).c_str());
            MatchingSiPM[NUCLEUS][100 + det]->Write();
        }
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