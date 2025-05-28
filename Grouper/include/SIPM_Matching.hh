#ifndef SiPM_MATCHER_HH
#define SiPM_MATCHER_HH

#include "../../../lib/SignalDict/Signal.h"
#include "Detectors.hh"

/// FILE ///
map<string, TFile *> GROUPED_File;
map<string, TFile *> SIMULATED_File;
TFile *MATCHED_File;
TFile *f_tree;

TTree *Tree;
TTreeReader *Reader;
TTreeReaderArray<Signal> *Silicon;
TTreeReaderValue<vector<vector<pair<Signal, Signal>>>> *SiPM_Groups;

TTree *Tree_SIMULATED;
TTreeReader *Reader_SIMULATED;
TTreeReaderValue<int> *Silicon_code;
TTreeReaderValue<double> *Silicon_energy;
TTreeReaderValue<double> *SiPM_energy;

TTree *Tree_MATCHED;
vector<Signal> SiPM;

int Entry_MAX;
bool FULL = false;
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
map<string, TH2D *[SIGNAL_MAX]> H_SiPM_Gain;
map<string, TGraph *[SIGNAL_MAX]> G_SiPM_Gain;
// Results
TGraphErrors *gHighLow_RUNS[SIGNAL_MAX];
TGraphErrors *gSiPM_High_RUNS[SIGNAL_MAX];
TGraphErrors *gSiPM_Low_RUNS[SIGNAL_MAX];   

/// DIRECTORY ///
map<string, TDirectory *> dir;
map<string, TDirectory *[SIGNAL_MAX]> dir_detector;
map<string, TDirectory *[SIGNAL_MAX]> dir_SiPM;
map<string, TDirectory *> dir_Nucleus;

// FUCNTIONS //
TF1 *F_SiliconCalibration[SIGNAL_MAX];
TF1 *f_linear = new TF1("f_linear", "[0]*x + [1]", 0, 10e6);
TF1 *f_erf = new TF1("f_erf", "x < [0] ? [1]*x + [2] : [3]*(erf((x-[4])/[5]))", 0, 10000e3);
TF1 *TF1_DoubleGaussian;
TF1 *TF1_SimpleGaussian;
pair<double, double> Range_SiPM_LowHigh = make_pair(50e3, 150e3);
int current_detector;
double SiPM_Window[10] = {0, 1500e3, 1500e3, 1500e3};
// DATA //
string Nuclei[4] = {"90Sr", "207Bi", "32Ar", "33Ar"};
string TYPE;
string RUN;
string NUCLEI;
double SiliconCalibrationParameter[SIGNAL_MAX][3];
map<string, TF1 *[SIGNAL_MAX]> MatchingLowHigh;
map<string, TF1 *[SIGNAL_MAX]> MatchingSiPM;

map<string, vector<string>> Map_RunFiles;





void init()
{
    // Map_RunFiles["32Ar"] =
    //     {"057", "058", "059",
    //      "061", "062", "064", "065", "066", "067", "068", "069",
    //      "070", "071", "072", "074", "075", "076", 
    //     "077",

    //      "112", "113", "114", "115", "116", "118"};

    // Map_RunFiles["33Ar"] = {"078"};

    Map_RunFiles["207Bi"] = {"137"};

    Map_RunFiles["90Sr"] = {"133"};
    Map_RunFiles["32Ar"] = {"077"};

}

TF1 *InvertingLinear(TF1 *f)
{
    TF1 *f_inv = new TF1("f_inv", "[0]*x + [1]", 0, 10000e3);
    double a = f->GetParameter(0);
    double b = f->GetParameter(1);
    f_inv->SetParameter(0,1. / a);
    f_inv->SetParameter(1, -b / a);
    return f_inv;
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

                // init silicon detector window
                for (int i : Dir2Det(dir, strip))
                {
                    WindowsMap[nuclei][number][i] = make_pair(energy_low, energy_high);
                }
            }
        }
    }
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

int LoadValues()
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
            MatchingSiPM[RUN][det] = InvertingLinear((TF1 *)f->Get(("MatchingSiPM_" + RUN + "_" + detectorName[det]).c_str()));
            if (IsDetectorBetaHigh(det))
            {
                MatchingLowHigh[RUN][det] = (TF1 *)f->Get(("MatchingLowHigh_" + RUN + "_" + detectorName[det]).c_str());
            }
        }
    }

    // f->Close();
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
            gHighLow_RUNS[i] = new TGraphErrors();
            gSiPM_High_RUNS[i] = new TGraphErrors();
        }
        if (IsDetectorBetaLow(i))
        {
            gSiPM_Low_RUNS[i] = new TGraphErrors();
        }
    }

    Info("Init Histograms done");
}

void ReadData()
{
    Info("Reading Data for Fitting Low-High", 1);
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
            int sili_label = (*Silicon)[1].Label;
            double sili_energy = F_SiliconCalibration[sili_label]->Eval((*Silicon)[1].Channel / 1000);
            // gate on proton peak
            int peak=0;
            if (NUCLEI == "32Ar") peak = 14;
            if (NUCLEI == "33Ar") peak = 21;
            if (sili_energy < WindowsMap[NUCLEI][peak][sili_label].first || sili_energy > WindowsMap[NUCLEI][peak][sili_label].second)
                continue;

        }

        // Looping on SiPM groups
        for (int i_groups = 0; i_groups < (**SiPM_Groups).size(); i_groups++)
        {
            // if ((**SiPM_Groups)[i_groups].size() != 9)
            //     continue;
            for (int i_pair = 0; i_pair < (**SiPM_Groups)[i_groups].size(); i_pair++)
            {
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


                // Compare SiPMS to 104
                if ((**SiPM_Groups)[i_groups][i_pair].first.Label == 104)
                {
                    for (int i_pair_compare = 0; i_pair_compare < (**SiPM_Groups)[i_groups].size(); i_pair_compare++)
                    {
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
            int sili_label = (*Silicon)[1].Label;
            double sili_energy = F_SiliconCalibration[sili_label]->Eval((*Silicon)[1].Channel / 1000);
            // gate on proton peak
            int peak=0;
            if (NUCLEI == "32Ar") peak = 14;
            if (NUCLEI == "33Ar") peak = 21;
            if (sili_energy < WindowsMap[NUCLEI][peak][sili_label].first || sili_energy > WindowsMap[NUCLEI][peak][sili_label].second)
                continue;

        }

        // cout << "Applying correction" << endl;
        for (int i_groups = 0; i_groups < (**SiPM_Groups).size(); i_groups++)
        {
            // cout << "Group " << i_groups << endl;
            for (int i_pair = 0; i_pair < (**SiPM_Groups)[i_groups].size(); i_pair++)
            {
                // APPLYING MATCHING High Low
                // cout << "Applying correction HL" << endl;
                if ((**SiPM_Groups)[i_groups][i_pair].second.isValid)
                {
                    // cout << "Applying on " << (**SiPM_Groups)[i_groups][i_pair].second.Label << endl;
                    (**SiPM_Groups)[i_groups][i_pair].second.Channel = MatchingLowHigh[RUN][(**SiPM_Groups)[i_groups][i_pair].second.Label-10]->Eval((**SiPM_Groups)[i_groups][i_pair].second.Channel);
                }

                // cout << "Applying correction SiPM" << endl;
                // APPLYING MATCHING SiPM
                if ((**SiPM_Groups)[i_groups][i_pair].first.isValid)
                {
                    // cout << "Applying on " << (**SiPM_Groups)[i_groups][i_pair].first.Label << endl;
                    (**SiPM_Groups)[i_groups][i_pair].first.Channel = MatchingSiPM[RUN][(**SiPM_Groups)[i_groups][i_pair].first.Label]->Eval((**SiPM_Groups)[i_groups][i_pair].first.Channel);
                }
                if ((**SiPM_Groups)[i_groups][i_pair].second.isValid)
                {
                    // cout << "Applying on " << (**SiPM_Groups)[i_groups][i_pair].second.Label << endl;
                    (**SiPM_Groups)[i_groups][i_pair].second.Channel = MatchingSiPM[RUN][(**SiPM_Groups)[i_groups][i_pair].second.Label-10]->Eval((**SiPM_Groups)[i_groups][i_pair].second.Channel);
                }
            }
        }
        // cout << "Correction applied" << endl;
        

        // cout << "Filling histograms" << endl;
        // Looping on SiPM groups
        for (int i_groups = 0; i_groups < (**SiPM_Groups).size(); i_groups++)
        {
            // cout << "Group " << i_groups << endl;
            if ((**SiPM_Groups)[i_groups].size() != 9)
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
                if ((**SiPM_Groups)[i_groups][i_pair].first.Label == 104)
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

void FittingLowHigh(int Verbose = 0)
{
    Info("Fitting High-Low gain factor for run " + RUN);

    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBetaHigh(i))
        {
            if (Verbose == 1)
                Info(detectorName[i], 1);
            current_detector = i;
            G_SiPM_HighLow[RUN][i]->Fit(f_linear, "QR", "", Range_SiPM_LowHigh.first, Range_SiPM_LowHigh.second);
            MatchingLowHigh[RUN][i] = G_SiPM_HighLow[RUN][i]->GetFunction("f_linear");
        }
    }
}

void FittingSiPM(int Verbose = 0)
{
    Info("Fitting SiPM gain factor for run " + RUN);
    f_linear->FixParameter(1, 0);

    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBeta(i))
        {
            if (Verbose == 1)
                Info(detectorName[i], 1);
            current_detector = i;
            if (IsDetectorBetaLow(i))
                G_SiPM_Gain[RUN][i]->Fit(f_linear, "Q");
            else
                G_SiPM_Gain[RUN][i]->Fit(f_linear, "QR", "", 200e3, Range_SiPM_LowHigh.second*10);
            MatchingSiPM[RUN][i] = G_SiPM_Gain[RUN][i]->GetFunction("f_linear");
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
            G_SiPM_HighLow[RUN][i]->Write();

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
            G_SiPM_HighLow[RUN][i]->Draw("AP");
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
            G_SiPM_Gain[RUN][i]->Draw("AP");
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
    MATCHED_File = MyTFile(DIR_ROOT_DATA_MATCHED + "SiPM_Matching.root", "UPDATE", "Q");
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
    

    for (int i = 0; i < SIGNAL_MAX; i++)
    {

        if (IsDetectorBetaHigh(i))
        {
            if (Verbose == 1)
                Info(detectorName[i], 1);

            MATCHED_File->cd((RUN + "/SiPM_" + to_string(GetDetectorChannel(i))).c_str());
            H_SiPM_High[RUN][i]->Write();
            H_SiPM_HighLow[RUN][i]->Write();
            G_SiPM_HighLow[RUN][i]->Write();

            c_High->cd();
            H_SiPM_High[RUN][i]->Draw("HIST SAME");

            /////////// ##### FIITING HIGH LOW ###### ///////////
            // TH2D and FIT
            cHighLowSiPM->cd(GetDetectorChannel(i));
            H_SiPM_HighLow[RUN][i]->Draw("COLZ");
            // TGraph and FIT
            cHighLowSiPM_Fitting->cd(GetDetectorChannel(i));
            G_SiPM_HighLow[RUN][i]->SetTitle(("SiPM " + to_string(GetDetectorChannel(i))).c_str());
            G_SiPM_HighLow[RUN][i]->GetXaxis()->SetTitle("Low SiPM");
            G_SiPM_HighLow[RUN][i]->GetYaxis()->SetTitle("High SiPM");
            G_SiPM_HighLow[RUN][i]->GetXaxis()->CenterTitle();
            G_SiPM_HighLow[RUN][i]->GetYaxis()->CenterTitle();
            G_SiPM_HighLow[RUN][i]->Draw("AP");

            /////////// ##### FIITING SiPM Gain ###### ///////////
            // TH2D and FIT
            cSiPM_GainH->cd(GetDetectorChannel(i));
            H_SiPM_Gain[RUN][i]->Draw("COLZ");
            // TGraph and FIT
            cSiPM_Gain_FittingH->cd(GetDetectorChannel(i));
            G_SiPM_Gain[RUN][i]->SetTitle(("SiPM " + to_string(GetDetectorChannel(i))).c_str());
            G_SiPM_Gain[RUN][i]->GetXaxis()->SetTitle("Low SiPM");
            G_SiPM_Gain[RUN][i]->GetYaxis()->SetTitle("High SiPM");
            G_SiPM_Gain[RUN][i]->GetXaxis()->CenterTitle();
            G_SiPM_Gain[RUN][i]->GetYaxis()->CenterTitle();
            G_SiPM_Gain[RUN][i]->Draw("AP");
        }
        if (IsDetectorBetaLow(i))
        {
            
            if (Verbose == 1)
                Info(detectorName[i], 1);
            MATCHED_File->cd((RUN + "/SiPM_" + to_string(GetDetectorChannel(i))).c_str());
            H_SiPM_Low[RUN][i]->Write();

            c_Low->cd();
            H_SiPM_Low[RUN][i]->Draw("HIST SAME");

            /////////// ##### FIITING SiPM Gain ###### ///////////
            // TH2D and FIT
            cSiPM_GainL->cd(GetDetectorChannel(i));
            H_SiPM_Gain[RUN][i]->Draw("COLZ");

            // // TGraph and FIT
            cSiPM_Gain_FittingH->cd(GetDetectorChannel(i));
            G_SiPM_Gain[RUN][i]->SetTitle(("SiPM " + to_string(GetDetectorChannel(i))).c_str());
            G_SiPM_Gain[RUN][i]->GetXaxis()->SetTitle("Low SiPM");
            G_SiPM_Gain[RUN][i]->GetYaxis()->SetTitle("High SiPM");
            G_SiPM_Gain[RUN][i]->GetXaxis()->CenterTitle();
            G_SiPM_Gain[RUN][i]->GetYaxis()->CenterTitle();
            G_SiPM_Gain[RUN][i]->Draw("AP");
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
    MATCHED_File = MyTFile(DIR_ROOT_DATA_MATCHED + "SiPM_Matching.root", "UPDATE", "Q");
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
            string type = "multifast";
            if (pairr.first == "90Sr")
                type = "data";
            GROUPED_File[run] = MyTFile(DIR_ROOT_DATA_GROUPED + "run_" + run + "_" + type + "_" + pairr.first + "_grouped.root", "READ", "Q");
            if (GROUPED_File[run] == NULL)
                continue;
            else
                GROUPED_File[run]->Close();

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
                    gHighLow_RUNS[i]->SetPoint(counterHighLow[i], stoi(run), MatchingLowHigh[run][i]->GetParameter(0));
                    gHighLow_RUNS[i]->SetPointError(counterHighLow[i], 0, MatchingLowHigh[run][i]->GetParError(0));
                    counterHighLow[i]++;

                    // All SiPM factor on the same canvas for all runs for High
                    MatchingSiPM[run][i] = (TF1 *)MATCHED_File->Get((run + "/SiPM_" + to_string(GetDetectorChannel(i)) + "/MatchingSiPM_" + run + "_" + detectorName[i]).c_str());
                    gSiPM_High_RUNS[i]->SetPoint(counterSiPM_H[i], stoi(run), MatchingSiPM[run][i]->GetParameter(0));
                    gSiPM_High_RUNS[i]->SetPointError(counterSiPM_H[i], 0, MatchingSiPM[run][i]->GetParError(0));
                    counterSiPM_H[i]++;

                    
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
                    gSiPM_Low_RUNS[i]->SetPoint(counterSiPM_L[i], stoi(run), MatchingSiPM[run][i]->GetParameter(0));
                    gSiPM_Low_RUNS[i]->SetPointError(counterSiPM_L[i], 0, MatchingSiPM[run][i]->GetParError(0));
                    counterSiPM_L[i]++;

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
            TCanvas *c = new TCanvas(("Matching_HighLow" + to_string(GetDetectorChannel(i))).c_str(), ("Matching_HighLow" + to_string(GetDetectorChannel(i))).c_str(), 1920, 1080);
            gHighLow_RUNS[i]->SetMarkerStyle(20);
            gHighLow_RUNS[i]->SetMarkerColor(kBlack);
            gHighLow_RUNS[i]->SetTitle(("Matching " + detectorName[i]).c_str());
            gHighLow_RUNS[i]->GetXaxis()->SetTitle("Run");
            gHighLow_RUNS[i]->GetYaxis()->SetTitle("Gain factor Low/High");
            gHighLow_RUNS[i]->Draw("AP");

            gHighLow_RUNS[i]->Fit("pol0", "Q");
            gHighLow_RUNS[i]->GetFunction("pol0")->SetLineColor(kRed);
            gHighLow_RUNS[i]->GetFunction("pol0")->Draw("SAME");

            c->Write();

            TCanvas *cSiPM = new TCanvas(("Matching_SiPM_High" + to_string(GetDetectorChannel(i))).c_str(), ("Matching_SiPM_High" + to_string(GetDetectorChannel(i))).c_str(), 1920, 1080);
            gSiPM_High_RUNS[i]->SetMarkerStyle(20);
            gSiPM_High_RUNS[i]->SetMarkerColor(kBlack);
            gSiPM_High_RUNS[i]->SetTitle(("Matching SiPM with" + detectorName[i]).c_str());
            gSiPM_High_RUNS[i]->GetXaxis()->SetTitle("Run");
            gSiPM_High_RUNS[i]->GetYaxis()->SetTitle("Gain factor SiPM");
            gSiPM_High_RUNS[i]->Draw("AP");
            gSiPM_High_RUNS[i]->Fit("pol0", "Q");
            gSiPM_High_RUNS[i]->GetFunction("pol0")->SetLineColor(kBlack);
            gSiPM_High_RUNS[i]->GetFunction("pol0")->Draw("SAME");


            cSiPM->Write();
        }

        if (IsDetectorBetaLow(i))
        {
            TCanvas *cSiPM = new TCanvas(("Matching_SiPM_Low" + to_string(GetDetectorChannel(i))).c_str(), ("Matching_SiPM_Low" + to_string(GetDetectorChannel(i))).c_str(), 1920, 1080);
            gSiPM_Low_RUNS[i]->SetMarkerStyle(20);
            gSiPM_Low_RUNS[i]->SetMarkerColor(kBlack);
            gSiPM_Low_RUNS[i]->SetTitle(("Matching SiPM with" + detectorName[i]).c_str());
            gSiPM_Low_RUNS[i]->GetXaxis()->SetTitle("Run");
            gSiPM_Low_RUNS[i]->GetYaxis()->SetTitle("Gain factor SiPM");
            gSiPM_Low_RUNS[i]->Draw("AP");
            gSiPM_Low_RUNS[i]->Fit("pol0", "Q");
            gSiPM_Low_RUNS[i]->GetFunction("pol0")->SetLineColor(kBlack);
            gSiPM_Low_RUNS[i]->GetFunction("pol0")->Draw("SAME");

            cSiPM->Write();
        }
    }

    MATCHED_File->Close();
    MATCHED_File = MyTFile(DIR_ROOT_DATA_MATCHED + "SiPM_Matching.root", "UPDATE", "Q");

    Info("Write Histograms stop");
}

void InitSiliconCalibration()
{

    string CalibFileName;

    CalibFileName = "/Config_Files/" + to_string(YEAR) + "/Calibration_" + to_string(YEAR) + ".txt";

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