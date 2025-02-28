#ifndef SIPM_CALIBRATION_HH
#define SIPM_CALIBRATION_HH

#include "../../../lib/SignalDict/Signal.h"
#include "Detectors.hh"

int VERBOSE = 0;

/// FILE ///
map<string, TFile *> GROUPED_File;
map<string, TFile *> SIMULATED_File;
TFile *CALIBRATED_File;
TFile *f_tree;
TFile *f_simulated_tree;

TTree *Tree;
TTreeReader *Reader;
TTreeReaderValue<vector<vector<pair<Signal, Signal>>>> *SiPM_Groups;
TTreeReaderArray<Signal> *Silicon;
TTreeReaderArray<Signal> *signals;

TTree *Tree_SIMULATED;
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
map<string, TTree *[50]> Tree_Peaks_Simulated;
Signal SiPMi;
double SiPMEnergy;

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
map<string, TH2D *[SIGNAL_MAX]> H_Matching_SiPM;
map<string, TGraph *[SIGNAL_MAX]> G_SiPM_HighLow_Matched;
map<string, TH2D *[SIGNAL_MAX]> H_SiPM_HighLow_Matched_diff;
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
map<string, TDirectory *[50]> dir_peaks;

// FUCNTIONS //
TF1 *F_SiliconCalibration[SIGNAL_MAX];
TF1 *F_MatchingSiPM[SIGNAL_MAX];
TF1 *F_MatchingLowHigh[SIGNAL_MAX];
TF1 *f_linear = new TF1("f_linear", "[0]*x + [1]", 0, 10000e3);
TF1 *f_erf = new TF1("f_erf", "x < [0] ? [1]*x + [2] : [3]*(erf((x-[4])/[5]))", 0, 10000e3);
pair<double, double> Range_SiPM_LowHigh = make_pair(50e3, 200e3);
TF1 *gauss;
TF1 *Threshold_f;
int current_detector;
double SiPM_Window[10] = {0, 1500e3, 1500e3, 1500e3};
vector<double> res;
// DATA //
string Nuclei[1] = {"32Ar"}; //, "32Ar_thick", "33Ar"};
string Nucleis[3] = {"32Ar", "207Bi", "90Sr"};
string NUCLEUS;
double SiliconCalibrationParameter[SIGNAL_MAX][3];
map<string, TF1 *[SIGNAL_MAX]> MatchingLowHigh;
map<string, TF1 *[SIGNAL_MAX]> MatchingLowHigh_erf;
map<string, TF1 *[SIGNAL_MAX]> MatchingSiPM;
map<string, pair<double, double>> SiPM_Range;

void InitWindows(int verbose = 0, string addpath="")
{
    Info("Init Windows");
    string direction[2] = {"Up", "Down"};
    for (auto dir : direction)
    {
        if (verbose > 0)
            Info("Direction: " + dir, 1);
        for (int strip = 1; strip <= 5; strip++)
        {
            if (verbose > 0)
                Info("Strip: " + to_string(strip), 2);
            ifstream file(addpath + "Config_Files/Detector_Window/" + dir + "_" + to_string(strip) + ".txt");
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
    Info("Init Windows Done");
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

TF1* InvertingLinear(TF1* f)
{
    TF1* f_inv = new TF1("f_inv", "[0]*x + [1]", 0, 10000e3);
    double a = f->GetParameter(0);
    double b = f->GetParameter(1);
    f_inv->SetParameter(0, 1 / a);
    f_inv->SetParameter(1, -b / a);
    return f_inv;
}

void InitMatchingSiPM()
{
    Info("Init Matching SiPM");
    TFile *f = MyTFile((DIR_ROOT_DATA_MATCHED + "SiPM_Matching_values.root").c_str(), "READ");

    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBetaHigh(i))
        {
            F_MatchingSiPM[GetDetectorChannel(i)] = InvertingLinear((TF1 *)f->Get(("MatchingSiPM_077_" + detectorName[i]).c_str()));
            if (F_MatchingSiPM[GetDetectorChannel(i)] == nullptr)
            {
                Error("No MatchingSiPM_077_" + detectorName[i] + " found");
            }
        }

        if (IsDetectorBetaHigh(i))
        {
            F_MatchingLowHigh[GetDetectorChannel(i)] = (TF1 *)f->Get(("MatchingLowHigh_077_" + detectorName[i]).c_str());
            if (F_MatchingLowHigh[GetDetectorChannel(i)] == nullptr)
            {
                Error("No MatchingLowHigh_077_" + detectorName[i] + " found");
            }
        }
    }
}



void ReadTree()
{
    Info("Read Tree");
    
    for (int peak = 0; peak <= 50; peak++)
    {
        for (int det = 1; det <= 9; det++)
        {
            if (WindowsMap["32Ar"][peak][11].first == -1 || !WindowsMap["32Ar"][peak][11].first)
                continue;
            Tree_Peaks["32Ar"][peak][det] = (TTree *)f_tree->Get(("Tree_Peaks_32Ar_" + to_string(peak) + "_" + to_string(det)).c_str());
        }
    }
    for (int det = 1; det <= 9; det++)
    {
        Tree_Peaks["207Bi"][0][det] = (TTree *)f_tree->Get(("Tree_Peaks_207Bi_" + to_string(0) + "_" + to_string(det)).c_str());
        Tree_Peaks["90Sr"][0][det] = (TTree *)f_tree->Get(("Tree_Peaks_90Sr_" + to_string(0) + "_" + to_string(det)).c_str());
    }
}

void WriteTree()
{
    Info("Write Tree");
    f_tree->cd();
    for (int peak = 0; peak <= 50; peak++)
    {
        for (int det = 1; det <= 9; det++)
        {
            if (WindowsMap["32Ar"][peak][11].first == -1 || !WindowsMap["32Ar"][peak][11].first)
                continue;
            Tree_Peaks["32Ar"][peak][det]->Write();
        }
    }
    for (int det = 1; det <= 9; det++)
    {
        Tree_Peaks["207Bi"][0][det]->Write();
        Tree_Peaks["90Sr"][0][det]->Write();
    }
    f_tree->Close();
    f_tree = MyTFile((DIR_ROOT_DATA_CALIBRATED + "Matching_SiPM_trees.root").c_str(), "READ");
}

void ExtractingTree()
{
    InitMatchingSiPM();
    Info("Extracting Tree");
    for (string NUCLEUS : Nucleis)
    {
        Info("Nucleus: " + NUCLEUS);

        if (NUCLEUS == "32Ar")
        {
            Tree = (TTree *)GROUPED_File[NUCLEUS]->Get("MERGED_Tree");
            Reader = new TTreeReader(Tree);
            Silicon = new TTreeReaderArray<Signal>(*Reader, "MERGED_Tree_Silicon");
            SiPM_Groups = new TTreeReaderValue<vector<vector<pair<Signal, Signal>>>>(*Reader, "MERGED_Tree_SiPMGroup");
        }
        else
        {
            Tree = (TTree *)GROUPED_File[NUCLEUS]->Get("CLEANED_Tree");
            Reader = new TTreeReader(Tree);
            SiPM_Groups = new TTreeReaderValue<vector<vector<pair<Signal, Signal>>>>(*Reader, "CLEANED_Tree_SiPMGroup");
        }

        clock_t start = clock(), Current;
        int Entries = Reader->GetEntries();
        int Entry_MAX = 1e7;
        while (Reader->Next() && Reader->GetCurrentEntry() < Entry_MAX)
        {
            ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "Reading Tree");

            double energy;
            int sili_code;
            int peak = 0;
            if (NUCLEUS == "32Ar") ///// ### IF MORE THAN 1 PEAK IN THE WINDOW THE LOWEST IN ENERGY IS TAKEN
            {
                energy = F_SiliconCalibration[(*Silicon)[1].Label]->Eval((*Silicon)[1].Channel / 1000);
                sili_code = (*Silicon)[1].Label;

                for (int peak_number = 0; peak_number < 50; peak_number++)
                {
                    SiPMi = Signal();
                    if (WindowsMap[NUCLEUS][peak_number][sili_code].first == -1 || !WindowsMap[NUCLEUS][peak_number][sili_code].first)
                        continue;

                    if (energy > WindowsMap[NUCLEUS][peak_number][sili_code].first & energy < WindowsMap[NUCLEUS][peak_number][sili_code].second)
                    {
                        peak = peak_number;
                        break;
                    }
                }
            }
            for (int i_group = 0; i_group < (**SiPM_Groups).size(); i_group++)
            {
                if ((**SiPM_Groups)[i_group].size() <= 8)
                    continue;
                for (int i_pair = 0; i_pair < (**SiPM_Groups)[i_group].size(); i_pair++)
                {
                    SiPMi = Signal();
                    if ((**SiPM_Groups)[i_group][i_pair].first.isValid)
                    {
                        // cout << (**SiPM_Groups)[i_group][i_pair].first.Channel << endl;
                        (**SiPM_Groups)[i_group][i_pair].first.Channel = F_MatchingSiPM[GetDetectorChannel((**SiPM_Groups)[i_group][i_pair].first.Label)]->Eval((**SiPM_Groups)[i_group][i_pair].first.Channel);
                        // cout << (**SiPM_Groups)[i_group][i_pair].first.Channel << endl;
                        if ((**SiPM_Groups)[i_group][i_pair].first.Channel < eHighLowLimit)
                        {
                            SiPMi = (**SiPM_Groups)[i_group][i_pair].first;
                        }

                        else if ((**SiPM_Groups)[i_group][i_pair].second.isValid)
                        {
                            SiPMi = (**SiPM_Groups)[i_group][i_pair].second;
                            SiPMi.Channel = F_MatchingSiPM[GetDetectorChannel(SiPMi.Label)]->Eval(SiPMi.Channel);
                            SiPMi.Channel = F_MatchingLowHigh[GetDetectorChannel(SiPMi.Label)]->Eval(SiPMi.Channel);
                        }
                    }

                    if (!SiPMi.isValid)
                        continue;

                    Tree_Peaks[NUCLEUS][peak][GetDetectorChannel(SiPMi.Label)]->Fill();
                }
            }
        }
    }
    WriteTree();
}

void InitTree(string option = "READ")
{
    Info("Init Tree");
    f_tree = MyTFile((DIR_ROOT_DATA_CALIBRATED + "Matching_SiPM_trees.root").c_str(), option);
    NUCLEUS = "32Ar";

    if (option == "RECREATE")
    {
        for (int peak_number = 0; peak_number < 50; peak_number++)
        {
            SiPMi = Signal();
            if (WindowsMap[NUCLEUS][peak_number][11].first == -1 || !WindowsMap[NUCLEUS][peak_number][11].first)
                continue;
            for (int det = 1; det <= 9; det++)
            {
                f_tree->cd();
                Tree_Peaks[NUCLEUS][peak_number][det] = new TTree(("Tree_Peaks_" + NUCLEUS + "_" + to_string(peak_number) + "_" + to_string(det)).c_str(), ("Tree_Peaks_" + NUCLEUS + "_" + to_string(peak_number) + "_" + to_string(det)).c_str());
                Tree_Peaks[NUCLEUS][peak_number][det]->Branch("SiPM", &SiPMi);
            }
        }

        // init trees
        NUCLEUS = "207Bi";
        for (int det = 1; det <= 9; det++)
        {
            f_tree->cd();
            Tree_Peaks["207Bi"][0][det] = new TTree(("Tree_Peaks_" + NUCLEUS + "_" + to_string(0) + "_" + to_string(det)).c_str(), ("Tree_Peaks_" + NUCLEUS + "_" + to_string(0) + "_" + to_string(det)).c_str());
            Tree_Peaks["207Bi"][0][det]->Branch("SiPM", &SiPMi);
        }

        NUCLEUS = "90Sr";
        for (int det = 1; det <= 9; det++)
        {
            f_tree->cd();
            Tree_Peaks["90Sr"][0][det] = new TTree(("Tree_Peaks_" + NUCLEUS + "_" + to_string(0) + "_" + to_string(det)).c_str(), ("Tree_Peaks_" + NUCLEUS + "_" + to_string(0) + "_" + to_string(det)).c_str());
            Tree_Peaks["90Sr"][0][det]->Branch("SiPM", &SiPMi);
        }

        NUCLEUS = "32Ar";
        for (int det = 1; det <= 9; det++)
        {
            f_tree->cd();
            Tree_Peaks["32Ar"][0][det] = new TTree(("Tree_Peaks_" + NUCLEUS + "_" + to_string(0) + "_" + to_string(det)).c_str(), ("Tree_Peaks_" + NUCLEUS + "_" + to_string(0) + "_" + to_string(det)).c_str());
            Tree_Peaks["32Ar"][0][det]->Branch("SiPM", &SiPMi);
        }
        ExtractingTree();
        ReadTree();
    }
    else
    {
        ReadTree();
    }
}


void ReadSimulatedTree()
{
    Info("Read Simulated Tree");
    for (int peak = 0; peak <= 50; peak++)
    {
        if (WindowsMap["32Ar"][peak][11].first == -1 || !WindowsMap["32Ar"][peak][11].first)
                continue;
        Tree_Peaks_Simulated["32Ar"][peak] = (TTree *)f_simulated_tree->Get(("Tree_Peaks_Simulated_32Ar_" + to_string(peak)).c_str());
    }
    Tree_Peaks_Simulated["207Bi"][0] = (TTree *)f_simulated_tree->Get(("Tree_Peaks_Simulated_207Bi_" + to_string(0)).c_str());
    Tree_Peaks_Simulated["90Sr"][0] = (TTree *)f_simulated_tree->Get(("Tree_Peaks_Simulated_90Sr_" + to_string(0)).c_str());
}

void WriteSimulatedTree()
{
    Info("Write Simulated Tree");
    f_simulated_tree->cd();
    for (int peak = 0; peak <= 50; peak++)
    {
        if (WindowsMap["32Ar"][peak][11].first == -1 || !WindowsMap["32Ar"][peak][11].first)
                continue;
        Tree_Peaks_Simulated["32Ar"][peak]->Write();
    }
    Tree_Peaks_Simulated["207Bi"][0]->Write();
    Tree_Peaks_Simulated["90Sr"][0]->Write();
    f_simulated_tree->Close();
    f_simulated_tree = MyTFile((DIR_ROOT_DATA_SIMULATED + "Simulated_SiPM_trees.root").c_str(), "READ");
}

void ExtractingSimulatedTree()
{
    Info("Extracting Simulated Tree");
    for (string NUCLEUS : Nucleis)
    {
        Info("Nucleus: " + NUCLEUS, 1);
        Tree_SIMULATED = (TTree *)SIMULATED_File[NUCLEUS]->Get("PlasticIAS");
        Reader_SIMULATED = new TTreeReader(Tree_SIMULATED);
        Silicon_code = new TTreeReaderValue<int>(*Reader_SIMULATED, "Code");
        Silicon_energy = new TTreeReaderValue<double>(*Reader_SIMULATED, "Energy");
        SiPM_energy = new TTreeReaderValue<double>(*Reader_SIMULATED, "SiPM");

        clock_t start = clock(), Current;
        int Entries = Reader_SIMULATED->GetEntries();
        while (Reader_SIMULATED->Next() && Reader_SIMULATED->GetCurrentEntry() < 2e6)
        {
            ProgressBar(Reader_SIMULATED->GetCurrentEntry(), Entries, start, Current, "Reading Tree");
            int sili_code = **Silicon_code;
            double sili_e = **Silicon_energy;
            double SiPM_e = **SiPM_energy;

            if (NUCLEUS == "207Bi" || NUCLEUS == "90Sr")
            {
                SiPMEnergy = SiPM_e;
                Tree_Peaks_Simulated[NUCLEUS][0]->Fill();
            }
            else
            {
                for (int peak_number = 0; peak_number < 50; peak_number++)
                {
                    if (WindowsMap[NUCLEUS][peak_number][sili_code].first == -1 || !WindowsMap[NUCLEUS][peak_number][sili_code].first)
                        continue;

                    if (sili_e > WindowsMap[NUCLEUS][peak_number][sili_code].first & sili_e < WindowsMap[NUCLEUS][peak_number][sili_code].second)
                    {
                        SiPMEnergy = SiPM_e;
                        Tree_Peaks_Simulated[NUCLEUS][peak_number]->Fill();
                    }
                }
            }
        }
    }
    WriteSimulatedTree();
}

void InitSimulatedTree(string option = "READ")
{
    Info("Init Simulated Tree");
    f_simulated_tree = MyTFile((DIR_ROOT_DATA_SIMULATED + "Simulated_SiPM_trees.root").c_str(), option);
    NUCLEUS = "32Ar";

    if (option == "RECREATE")
    {
        for (int peak_number = 0; peak_number < 50; peak_number++)
        {
            SiPMi = Signal();
            if (WindowsMap[NUCLEUS][peak_number][11].first == -1 || !WindowsMap[NUCLEUS][peak_number][11].first)
                continue;
            for (int det = 1; det <= 9; det++)
            {
                f_simulated_tree->cd();
                Tree_Peaks_Simulated[NUCLEUS][peak_number] = new TTree(("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(peak_number)).c_str(), ("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(peak_number)).c_str());
                Tree_Peaks_Simulated[NUCLEUS][peak_number]->Branch("SiPM", &SiPMEnergy);
            }
        }

        // init trees
        NUCLEUS = "207Bi";
        for (int det = 1; det <= 9; det++)
        {
            f_simulated_tree->cd();
            Tree_Peaks_Simulated["207Bi"][0] = new TTree(("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(0)).c_str(), ("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(0)).c_str());
            Tree_Peaks_Simulated["207Bi"][0]->Branch("SiPM", &SiPMEnergy);
        }

        NUCLEUS = "90Sr";
        for (int det = 1; det <= 9; det++)
        {
            f_simulated_tree->cd();
            Tree_Peaks_Simulated["90Sr"][0] = new TTree(("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(0)).c_str(), ("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(0)).c_str());
            Tree_Peaks_Simulated["90Sr"][0]->Branch("SiPM", &SiPMEnergy);
        }

        NUCLEUS = "32Ar";
        for (int det = 1; det <= 9; det++)
        {
            f_simulated_tree->cd();
            Tree_Peaks_Simulated["32Ar"][0] = new TTree(("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(0)).c_str(), ("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(0)).c_str());
            Tree_Peaks_Simulated["32Ar"][0]->Branch("SiPM", &SiPMEnergy);
        }
        ExtractingSimulatedTree();
        ReadSimulatedTree();
    }
    else
    {
        ReadSimulatedTree();
    }
}

void InitHistograms(int verbose = 0)
{
    SiPM_Range["32Ar"] = make_pair(0, 6000);
    SiPM_Range["207Bi"] = make_pair(0, 2000);
    SiPM_Range["90Sr"] = make_pair(0, 4000);
    Info("Init Histograms");
    for (string Nucleus : Nucleis)
    {
        if (verbose == 1) Info("Nucleus: " + Nucleus, 1);
        NUCLEUS = Nucleus;
        if (CALIBRATED_File != nullptr) dir[NUCLEUS] = CALIBRATED_File->mkdir(NUCLEUS.c_str());       
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorBetaHigh(i))
            {
                if (verbose == 1) Info("Detector: " + detectorName[i], 2);
                for (int peak = 0; peak <= 50; peak++)
                {
                    if ((WindowsMap[NUCLEUS][peak][11].first == -1 || !WindowsMap[NUCLEUS][peak][11].first))
                        continue;

                    if (i == 101)
                    {
                        if (CALIBRATED_File != nullptr) dir_peaks[NUCLEUS][peak] = dir[NUCLEUS]->mkdir(("Peak_" + to_string(peak)).c_str());
                    }

                    H_SiPM[NUCLEUS][peak][GetDetectorChannel(i)] = new TH1D(("H_SiPM_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), ("H_SiPM_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), eHighN / 10, eHighMin, eHighMax);
                    H_SiPM[NUCLEUS][peak][GetDetectorChannel(i)]->GetXaxis()->SetTitle("Channel");
                    H_SiPM[NUCLEUS][peak][GetDetectorChannel(i)]->GetYaxis()->SetTitle("Counts");
                    H_SiPM[NUCLEUS][peak][GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
                    H_SiPM[NUCLEUS][peak][GetDetectorChannel(i)]->GetYaxis()->CenterTitle();

                    H_SiPM_Calibrated[NUCLEUS][peak][GetDetectorChannel(i)] = new TH1D(("H_SiPM_Calibrated_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), ("H_SiPM_Calibrated_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), eSiliN_cal / 10, eSiliMin_cal, eSiliMax_cal);
                    H_SiPM_Calibrated[NUCLEUS][peak][GetDetectorChannel(i)]->GetXaxis()->SetTitle("Energy [keV]");
                    H_SiPM_Calibrated[NUCLEUS][peak][GetDetectorChannel(i)]->GetYaxis()->SetTitle("Counts");
                    H_SiPM_Calibrated[NUCLEUS][peak][GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
                    H_SiPM_Calibrated[NUCLEUS][peak][GetDetectorChannel(i)]->GetYaxis()->CenterTitle();
                }

                if (NUCLEUS == "207Bi" || NUCLEUS == "90Sr")
                {
                    H_SiPM[NUCLEUS][0][GetDetectorChannel(i)] = new TH1D(("H_SiPM_" + NUCLEUS + "_Peak_" + to_string(0) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), ("H_SiPM_" + NUCLEUS + "_Peak_" + to_string(0) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), eHighN / 10, eHighMin, eHighMax);
                    H_SiPM[NUCLEUS][0][GetDetectorChannel(i)]->GetXaxis()->SetTitle("Channel");
                    H_SiPM[NUCLEUS][0][GetDetectorChannel(i)]->GetYaxis()->SetTitle("Counts");
                    H_SiPM[NUCLEUS][0][GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
                    H_SiPM[NUCLEUS][0][GetDetectorChannel(i)]->GetYaxis()->CenterTitle();

                    H_SiPM_Calibrated[NUCLEUS][0][GetDetectorChannel(i)] = new TH1D(("H_SiPM_Calibrated_" + NUCLEUS + "_Peak_" + to_string(0) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), ("H_SiPM_Calibrated_" + NUCLEUS + "_Peak_" + to_string(0) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), eSiliN_cal / 10, eSiliMin_cal, eSiliMax_cal);
                    H_SiPM_Calibrated[NUCLEUS][0][GetDetectorChannel(i)]->GetXaxis()->SetTitle("Energy [keV]");
                    H_SiPM_Calibrated[NUCLEUS][0][GetDetectorChannel(i)]->GetYaxis()->SetTitle("Counts");
                    H_SiPM_Calibrated[NUCLEUS][0][GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
                    H_SiPM_Calibrated[NUCLEUS][0][GetDetectorChannel(i)]->GetYaxis()->CenterTitle();
                }
            }
        }

        for (int peak = 0; peak <= 50; peak++)
        {
            if (WindowsMap[NUCLEUS][peak][11].first == -1 || !WindowsMap[NUCLEUS][peak][11].first)
                continue;

            H_Sim[NUCLEUS][peak] = new TH1D(("H_Sim_" + NUCLEUS + "_Peak_" + to_string(peak)).c_str(), ("H_Sim_" + NUCLEUS + "_Peak_" + to_string(peak)).c_str(), eSiliN_cal / 10, eSiliMin_cal, eSiliMax_cal);
            H_Sim[NUCLEUS][peak]->GetXaxis()->SetTitle("Energy [keV]");
            H_Sim[NUCLEUS][peak]->GetYaxis()->SetTitle("Counts");
            H_Sim[NUCLEUS][peak]->GetXaxis()->CenterTitle();
            H_Sim[NUCLEUS][peak]->GetYaxis()->CenterTitle();
            for (int det = 1; det <= 9; det++)
            {

                H_Sim_Conv[NUCLEUS][peak][det] = new TH1D(("H_Sim_Conv_" + NUCLEUS + "_Peak_" + to_string(peak) + detectorName[100 + det]).c_str(), ("H_Sim_Conv_" + NUCLEUS + "_Peak_" + to_string(peak) + detectorName[100 + det]).c_str(), eSiliN_cal / 10, eSiliMin_cal, eSiliMax_cal);
                H_Sim_Conv[NUCLEUS][peak][det]->GetXaxis()->SetTitle("Energy [keV]");
                H_Sim_Conv[NUCLEUS][peak][det]->GetYaxis()->SetTitle("Counts");
                H_Sim_Conv[NUCLEUS][peak][det]->GetXaxis()->CenterTitle();
                H_Sim_Conv[NUCLEUS][peak][det]->GetYaxis()->CenterTitle();
            }
        }

        if (NUCLEUS == "207Bi" || NUCLEUS == "90Sr")
        {
            H_Sim[NUCLEUS][0] = new TH1D(("H_Sim_" + NUCLEUS + "_Peak_" + to_string(0)).c_str(), ("H_Sim_" + NUCLEUS + "_Peak_" + to_string(0)).c_str(), eSiliN_cal / 10, eSiliMin_cal, eSiliMax_cal);
            H_Sim[NUCLEUS][0]->GetXaxis()->SetTitle("Energy [keV]");
            H_Sim[NUCLEUS][0]->GetYaxis()->SetTitle("Counts");
            H_Sim[NUCLEUS][0]->GetXaxis()->CenterTitle();
            H_Sim[NUCLEUS][0]->GetYaxis()->CenterTitle();
            for (int det = 1; det <= 9; det++)
            {

                H_Sim_Conv[NUCLEUS][0][det] = new TH1D(("H_Sim_Conv_" + NUCLEUS + "_Peak_" + to_string(0) + detectorName[100 + det]).c_str(), ("H_Sim_Conv_" + NUCLEUS + "_Peak_" + to_string(0) + detectorName[100 + det]).c_str(), eSiliN_cal / 10, eSiliMin_cal, eSiliMax_cal);
                H_Sim_Conv[NUCLEUS][0][det]->GetXaxis()->SetTitle("Energy [keV]");
                H_Sim_Conv[NUCLEUS][0][det]->GetYaxis()->SetTitle("Counts");
                H_Sim_Conv[NUCLEUS][0][det]->GetXaxis()->CenterTitle();
                H_Sim_Conv[NUCLEUS][0][det]->GetYaxis()->CenterTitle();
            }
        }
    }
}

double Chi2TreeHist_conv(const double *par)
{
    vector<double> chi2;
    for (string NUCLEUS : Nucleis)
    {
        if (NUCLEUS == "207Bi")
            continue;
        ////////////////////////////////////////////////////////////////
        // Init Histograms
        if (VERBOSE == 1) Info("Init Histograms", 1);
        for (int peak = 0; peak <= 50; peak++)
        {
            if ((WindowsMap[NUCLEUS][peak][11].first == -1 || !WindowsMap[NUCLEUS][peak][11].first))
                continue;
            H_Sim_Conv[NUCLEUS][peak][current_detector]->Reset();
            H_SiPM_Calibrated[NUCLEUS][peak][current_detector]->Reset();
        }

        if (NUCLEUS == "207Bi" || NUCLEUS == "90Sr")
        {
            cout << NUCLEUS << endl;
            H_Sim_Conv[NUCLEUS][0][current_detector]->Reset();
            H_SiPM_Calibrated[NUCLEUS][0][current_detector]->Reset();
        }

        ////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////
        // parameters
        if (VERBOSE == 1) Info("Setting Parameters", 1);
        double Calibration_OffSet = par[0];
        double Calibration = par[1];
        double Resolution_OffSet = par[2];
        double Resolution_SQRT = par[3];
        double Resolution_2 = par[4];
        double Threshold = par[5];
        double Threshold_STD = 0;
        ////////////////////////////////////////////////////////////////

        if (VERBOSE == 1) Info("Calibration of Experimental data", 1);
        int peak_number = 14;
        if (NUCLEUS == "207Bi" || NUCLEUS == "90Sr")
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
            H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector]->Fill(energy);
        }

        if (VERBOSE == 1) Info("Convolution with detector resolution", 1);
        // convoluting resolution on histogram
        // for (int i = 1; i <= H_Sim[NUCLEUS][peak_number]->GetNbinsX(); i++)
        // {
        //     if (H_Sim[NUCLEUS][peak_number]->GetBinContent(i) == 0)
        //         continue;
        //     energy = H_Sim[NUCLEUS][peak_number]->GetBinCenter(i);
        //     double sigma_resolution = sqrt(pow(Resolution_OffSet, 2) + pow(Resolution_SQRT * sqrt(energy), 2) + pow(Resolution_2 * pow(energy, 2), 2));
        //     // double sigma_resolution = sqrt(pow(coefficents[2].second, 2) + pow(coefficents[3].second * sqrt(energy), 2));

        //     gauss = new TF1("gauss", "gaus", 0, eSiliMax_cal);
        //     gauss->SetNpx(eSiliN_cal/10);
        //     gauss->SetParameters(H_Sim[NUCLEUS][peak_number]->GetBinContent(i) / (sigma_resolution * sqrt(2 * M_PI)), energy, sigma_resolution); /// weighted gaussian

        //     H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Add(gauss->GetHistogram());

        //     delete gauss;
        // }
        ///////////////
        // random convoluting resolution on histogram
        TTreeReader *Reader_Sim = new TTreeReader(Tree_Peaks_Simulated[NUCLEUS][peak_number]);
        TTreeReaderValue<double> *Energy = new TTreeReaderValue<double>(*Reader_Sim, "SiPM");

        while (Reader_Sim->Next())
        {
            double sigma_resolution = sqrt(pow(Resolution_OffSet, 2) + pow(Resolution_SQRT * sqrt(**Energy), 2) + pow(Resolution_2 * pow(**Energy, 2), 2));
            // shoot in a gaussian with ramdom engine
            double energy_conv = gRandom->Gaus(**Energy, sigma_resolution);

            if (energy_conv < Threshold)
                continue;

            bool flagUntriggered = false;
            for (int sipm = 1; sipm <= 8; sipm++)
            {
                if (gRandom->Gaus(**Energy, sigma_resolution) < Threshold)
                    flagUntriggered = true;
            }

            if (flagUntriggered)
                continue;

            H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Fill(gRandom->Gaus(**Energy, sigma_resolution));
        }
        ///////////////

        if (VERBOSE == 1) Info("Adding threshold", 1);

        // if (NUCLEUS == "207Bi")
        // {
        //     TF1 *fexp = new TF1("fexp", "expo", 0, 10000);
        //     fexp->SetNpx(10000);

        //     fexp->SetParLimits(0, 0, 200);
        //     fexp->SetParameter(0, 13);
        //     fexp->SetParLimits(1, -0.1, 0);
        //     fexp->SetParameter(1, -0.05);

        //     H_SiPM_Calibrated[NUCLEUS][0][current_detector]->Fit(fexp, "QR", "", 60, 100);

        //     TH1D *h = (TH1D *)fexp->GetHistogram()->Clone();
        //     H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Add(h);
        // }

        ////////// CONVOLUTING THRESHOLD//////////
        // Threshold_f = new TF1("Threshold", "0.5*(1+erf((x-[0])/[1]))", eSiliMin_cal, eSiliMax_cal);
        // Threshold_STD = 1. * sqrt(pow(Resolution_OffSet, 2) + pow(Resolution_SQRT * sqrt(Threshold), 2) + pow(Resolution_2 * pow(Threshold, 2), 2));
        // Threshold_f->SetParameters(Threshold, Threshold_STD);

        // for (int i = 1; i <= H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetNbinsX(); i++)
        // {
        //     H_Sim_Conv[NUCLEUS][peak_number][current_detector]->SetBinContent(i, H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetBinContent(i) * Threshold_f->Eval(H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetBinCenter(i)));
        // }

        // delete Threshold_f;

        //////////////// CHI2 ////////////////

        // NUCLEUS = "207Bi";
        H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector]->GetXaxis()->SetRangeUser(300, -111);
        H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetXaxis()->SetRangeUser(300, -1111);
        chi2.push_back(H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector]->Chi2Test(H_Sim_Conv[NUCLEUS][peak_number][current_detector], "CHI2/NDF"));
        // cout << chi2[chi2.size()-1] << "   " << Calibration_OffSet << "   " << Calibration << "   " << Resolution_OffSet << "   " << Resolution_SQRT << "   " << Resolution_2 << "   " << Threshold << "   " << Threshold_STD << endl;

        H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Scale(H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector]->Integral() / H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Integral());
        H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector]->GetXaxis()->SetRangeUser(0, SiPM_Range[NUCLEUS].second);
        H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetXaxis()->SetRangeUser(0, SiPM_Range[NUCLEUS].second);
    }

    double mean = accumulate(chi2.begin(), chi2.end(), 0.0) / chi2.size();
    Info("Chi2: " + to_string(mean), 1);
    Info("Calibration: " + to_string(par[0]) + "   " + to_string(par[1]) + "   " + to_string(par[2]) + "   " + to_string(par[3]) + "   " + to_string(par[4]) + "   " + to_string(par[5]) + "   " + to_string(par[6]), 1);
    return mean;
}

void CalibrationSiPM()
{
    VERBOSE = 0;
    Info("Calibration", 1);
    res = {0, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5};
    Minimizer *minimizer = Factory::CreateMinimizer("Minuit2", "Migrad");
    ROOT::Math::Functor functor(&Chi2TreeHist_conv, 7);
    vector<double> Par = {30.3, 0.68, 0, 4.0, 1e-05, 50, 15.0801};
    minimizer->SetFunction(functor);
    minimizer->SetLimitedVariable(0, "Calibration_OffSet", Par[0], 10, 0, 100);
    minimizer->SetLimitedVariable(1, "Calibration", Par[1], 0.1, 0.5, 1.2);
    minimizer->SetFixedVariable(2, "Resolution_OffSet", Par[2]);
    minimizer->SetLimitedVariable(3, "Resolution_SQRT", Par[3], 0.1, 2., 10);
    minimizer->SetLimitedVariable(4, "Resolution_2", Par[4], 0.1e-5, 1e-6, 8e-5);
    // minimizer->SetLimitedVariable(5, "Threshold", Par[5], 5, 20, 100);
    // minimizer->SetLimitedVariable(6, "Threshold_STD", Par[6], 10, 2, 20);
    minimizer->SetFixedVariable(5, "Threshold", Par[5]);
    minimizer->SetFixedVariable(6, "Threshold_STD", Par[6]);
    minimizer->SetMaxFunctionCalls(10000000);
    minimizer->SetMaxIterations(100000000);
    minimizer->SetTolerance(1e-3);
    minimizer->SetPrecision(10);
    // const double *bestPar = Par.data();
    minimizer->Minimize();

    const double *bestPar = minimizer->X();
    double chi2 = Chi2TreeHist_conv(bestPar);
}

void WriteHistograms()
{
    gStyle->SetOptStat(0);
    CALIBRATED_File->cd();
    Info("Write Histograms");
    NUCLEUS = "32Ar";
    Info("32Ar", 1);
    dir[NUCLEUS]->cd();

    TCanvas *cExp_Sim_all_32Ar[50];

    for (int peak = 1; peak <= 50; peak++)
    {
        if (WindowsMap[NUCLEUS][peak][11].first == -1 || !WindowsMap[NUCLEUS][peak][11].first)
            continue;
        cExp_Sim_all_32Ar[peak] = new TCanvas(("cExp_Sim1_" + NUCLEUS + "_Peak_" + to_string(peak)).c_str(), ("cExp_Sim1_" + NUCLEUS + "_Peak_" + to_string(peak)).c_str(), 800, 800);
        cExp_Sim_all_32Ar[peak]->Divide(3, 3);
    }

    for (int det = 1; det <= 9; det++)
    {
        for (int peak = 1; peak <= 50; peak++)
        {
            if (WindowsMap[NUCLEUS][peak][11].first == -1 || !WindowsMap[NUCLEUS][peak][11].first)
                continue;

            dir_peaks[NUCLEUS][peak]->cd();
            // dir_SiPM[NUCLEUS][det]->cd();
            TCanvas *cExp_Sim1 = new TCanvas(("cExp_Sim1_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(det)).c_str(), ("cExp_Sim1_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(det)).c_str(), 800, 800);
            H_SiPM_Calibrated[NUCLEUS][peak][det]->Draw("HIST");
            // H_Sim[NUCLEUS][peak]->SetLineColor(kRed);
            // H_Sim[NUCLEUS][peak]->Draw("HIST SAME");
            if (H_Sim_Conv[NUCLEUS][peak][det]->GetEntries() > 0)
            {
                H_Sim_Conv[NUCLEUS][peak][det]->SetLineColor(kRed);
                H_Sim_Conv[NUCLEUS][peak][det]->Draw("HIST SAME");
            }
            cExp_Sim1->Write();

            cExp_Sim_all_32Ar[peak]->cd(det);
            H_SiPM_Calibrated[NUCLEUS][peak][det]->Draw("HIST");
            if (H_Sim_Conv[NUCLEUS][peak][det]->GetEntries() > 0)
            {
                H_Sim_Conv[NUCLEUS][peak][det]->SetLineColor(kRed);
                H_Sim_Conv[NUCLEUS][peak][det]->Draw("HIST SAME");
            }
        }
    }

    for (int peak = 1; peak <= 50; peak++)
    {
        if (WindowsMap[NUCLEUS][peak][11].first == -1 || !WindowsMap[NUCLEUS][peak][11].first)
            continue;
        dir_peaks[NUCLEUS][peak]->cd();
        cExp_Sim_all_32Ar[peak]->Write();
    }

    NUCLEUS = "207Bi";
    Info(NUCLEUS, 1);
    dir[NUCLEUS]->cd();

    TCanvas *cExp_Sim_all = new TCanvas(("cExp_Sim1_" + NUCLEUS + "_Peak_" + to_string(0)).c_str(), ("cExp_Sim1_" + NUCLEUS + "_Peak_" + to_string(0)).c_str(), 800, 800);
    cExp_Sim_all->Divide(3, 3);
    for (int det = 1; det <= 9; det++)
    {
        int peak = 0;
        TCanvas *cExp_Sim1 = new TCanvas(("cExp_Sim1_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(det)).c_str(), ("cExp_Sim1_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(det)).c_str(), 800, 800);
        H_SiPM_Calibrated[NUCLEUS][peak][det]->Draw("HIST");
        // H_Sim[NUCLEUS][peak]->SetLineColor(kRed);
        // H_Sim[NUCLEUS][peak]->Draw("HIST SAME");
        H_Sim_Conv[NUCLEUS][peak][det]->SetLineColor(kRed);
        H_Sim_Conv[NUCLEUS][peak][det]->Draw("HIST SAME");
        cExp_Sim1->Write();

        cExp_Sim_all->cd(det);
        cExp_Sim_all->SetLogy();
        H_SiPM_Calibrated[NUCLEUS][peak][det]->Draw("HIST");
        H_Sim_Conv[NUCLEUS][peak][det]->SetLineColor(kRed);
        H_Sim_Conv[NUCLEUS][peak][det]->Draw("HIST SAME");
    }
    cExp_Sim_all->Write();

    NUCLEUS = "90Sr";
    Info(NUCLEUS, 1);
    dir[NUCLEUS]->cd();

    TCanvas *cExp_Sim_all_90Sr = new TCanvas(("cExp_Sim1_" + NUCLEUS + "_Peak_" + to_string(0)).c_str(), ("cExp_Sim1_" + NUCLEUS + "_Peak_" + to_string(0)).c_str(), 800, 800);
    cExp_Sim_all_90Sr->Divide(3, 3);
    for (int det = 1; det <= 9; det++)
    {
        int peak = 0;
        TCanvas *cExp_Sim1 = new TCanvas(("cExp_Sim1_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(det)).c_str(), ("cExp_Sim1_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(det)).c_str(), 800, 800);
        cExp_Sim1->SetLogy();
        H_SiPM_Calibrated[NUCLEUS][peak][det]->Draw("HIST");
        // H_Sim[NUCLEUS][peak]->SetLineColor(kRed);
        // H_Sim[NUCLEUS][peak]->Draw("HIST SAME");
        H_Sim_Conv[NUCLEUS][peak][det]->SetLineColor(kRed);
        H_Sim_Conv[NUCLEUS][peak][det]->Draw("HIST SAME");
        cExp_Sim1->Write();

        cExp_Sim_all_90Sr->cd(det);
        cExp_Sim_all_90Sr->SetLogy();
        H_SiPM_Calibrated[NUCLEUS][peak][det]->Draw("HIST");
        H_Sim_Conv[NUCLEUS][peak][det]->SetLineColor(kRed);
        H_Sim_Conv[NUCLEUS][peak][det]->Draw("HIST SAME");
    }
    cExp_Sim_all_90Sr->Write();
}

void InitSiliconCalibration(string addpath = "")
{

    string CalibFileName;

    CalibFileName = addpath+"./Config_Files/Calibration.txt";

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