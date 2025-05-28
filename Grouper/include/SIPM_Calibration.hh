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
map<string, double> Qbeta; 
map<string, double[50]> WindowsBetaMap; 
map<string, pair<int, int>> CanvasMap;
double SiPM_Window[10] = {0, 1500e3, 1500e3, 1500e3};
bool FLAG_GAIN_CHANGED = false;

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
map<string, pair<TH1D *, TH1D*>[50][SIGNAL_MAX]> H_SiPM_Calibrated;

//
map<string, TH2D *[SIGNAL_MAX]>H;
//

map<string, TH1D *[50]> H_Sim;
map<string, TH1D *[50][SIGNAL_MAX]> H_Sim_Conv;
/// DIRECTORY ///
map<string, TDirectory *> dir;
map<string, TDirectory *[SIGNAL_MAX]> dir_detector;
map<string, TDirectory *[SIGNAL_MAX]> dir_SiPM;
map<string, TDirectory *[50]> dir_peaks;

// FUCNTIONS //
TF1 *F_SiliconCalibration[SIGNAL_MAX];
map<string, TF1 *[SIGNAL_MAX]> F_MatchingSiPM;
map<string, TF1 *[SIGNAL_MAX]> F_MatchingLowHigh;
TF1 *f_linear = new TF1("f_linear", "[0]*x + [1]", 0, 10000e3);
TF1 *f_erf = new TF1("f_erf", "x < [0] ? [1]*x + [2] : [3]*(erf((x-[4])/[5]))", 0, 10000e3);
pair<double, double> Range_SiPM_LowHigh = make_pair(50e3, 200e3);
TF1 *gauss;
TF1 *Threshold_f;
int current_detector;
pair<TF1 *, TF1*> PileUp_Function[SIGNAL_MAX];

vector<double> Resolution_OffSet_det;
vector<double> Resolution_SQRT_det;
vector<double> Resolution_2_det;
map<string, double[SIGNAL_MAX]> SiPM_Threshold_Ratio;
// DATA //
string Nuclide[3] = {"32Ar", "207Bi", "90Sr"}; //33Ar
double DELTA_PileUp_Range = 200e3;
map<string, vector<int>> Map_Peak;
string NUCLEUS;
double SiliconCalibrationParameter[SIGNAL_MAX][3];
map<string, TF1 *[SIGNAL_MAX]> MatchingLowHigh;
map<string, TF1 *[SIGNAL_MAX]> MatchingLowHigh_erf;
map<string, TF1 *[SIGNAL_MAX]> MatchingSiPM;
map<string, pair<double, double>> SiPM_Range;

bool FLAGSAVINGCHI2 = false;
map<string, double[50][SIGNAL_MAX]> CHI2;

void InitWindows(int verbose = 0, string addpath = "")
{
    Info("Init Windows");

    SiPM_Range["32Ar"] = make_pair(0, 6000);
    SiPM_Range["207Bi"] = make_pair(0, 2000);
    SiPM_Range["90Sr"] = make_pair(0, 4000);

    Qbeta["32Ar"] = 11134-1200;
    Qbeta["207Bi"] = 1827+500;
    WindowsBetaMap["207Bi"][0] = Qbeta["207Bi"];
    Qbeta["90Sr"] = 2300+500;
    WindowsBetaMap["90Sr"][0] = Qbeta["90Sr"];
    
    CanvasMap["32Ar"] = make_pair(7, 5);
    CanvasMap["32Ar_thick"] = make_pair(7, 5);
    CanvasMap["33Ar"] = make_pair(10, 5);
    CanvasMap["18N"] = make_pair(3, 3);
    CanvasMap["207Bi"] = make_pair(1, 1);
    CanvasMap["90Sr"] = make_pair(1, 1);

    Map_Peak["32Ar"] = {6, 7, 12, 14, 15, 18, 23, 25, 27, 30};
    Map_Peak["33Ar"] = {};
    Map_Peak["207Bi"] = {0};
    Map_Peak["90Sr"] = {0};

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

                WindowsBetaMap[nuclei][number] = Qbeta[nuclei] - (energy_low+energy_high)/2;
            }
        }
    }
    Info("Init Windows Done");
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

void InitMatchingSiPM()
{
    Info("Init Matching SiPM");
    TFile *f = MyTFile((DIR_ROOT_DATA_MATCHED + "SiPM_Matching_values.root").c_str(), "READ");

    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBetaHigh(i))
        {
            F_MatchingSiPM["32Ar"][GetDetectorChannel(i)] = InvertingLinear((TF1 *)f->Get(("MatchingSiPM_077_" + detectorName[i]).c_str()));
            F_MatchingSiPM["207Bi"][GetDetectorChannel(i)] = InvertingLinear((TF1 *)f->Get(("MatchingSiPM_137_" + detectorName[i]).c_str()));
            F_MatchingSiPM["90Sr"][GetDetectorChannel(i)] = InvertingLinear((TF1 *)f->Get(("MatchingSiPM_133_" + detectorName[i]).c_str()));
        }

        if (IsDetectorBetaHigh(i))
        {
            F_MatchingLowHigh["32Ar"][GetDetectorChannel(i)] = (TF1 *)f->Get(("MatchingLowHigh_077_" + detectorName[i]).c_str());
            F_MatchingLowHigh["207Bi"][GetDetectorChannel(i)] = (TF1 *)f->Get(("MatchingLowHigh_137_" + detectorName[i]).c_str());
            F_MatchingLowHigh["90Sr"][GetDetectorChannel(i)] = (TF1 *)f->Get(("MatchingLowHigh_133_" + detectorName[i]).c_str());
        }
    }
}

void InitSiliconCalibration(string addpath = "")
{

    string CalibFileName;

    CalibFileName = addpath + "/Config_Files/" + to_string(YEAR) + "/Calibration_" + to_string(YEAR) + ".txt";

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

void InitThresholds()
{
    // // REFERENCE IS 32Ar
    string filename = "./Config_Files/SiPM_Thresholds.txt";
    ifstream file(filename);
    if (!file.is_open())
    {
        Error("Impossible to open the file: " + filename);
    }

    string line;
    double th;
    int det;
    string nuclei;
    while (getline(file, line))
    {
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
        ss >> det >> th;
        SiPM_Threshold_Ratio[nuclei][det] = th;
        if (nuclei == "32Ar")
        {
            continue;
        }
        SiPM_Threshold_Ratio[nuclei][det] = SiPM_Threshold_Ratio[nuclei][det] / SiPM_Threshold_Ratio["32Ar"][det];
        
    }

    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorBetaHigh(det))
        {
            SiPM_Threshold_Ratio["32Ar"][det] = 1.0;
        }
    }
}

void InitPileUp()
{
    Info("Init PileUp");

    // BEFORE GAIN CHANGE //
    TFile *f_pileup1 = MyTFile((DIR_ROOT_DATA_MATCHED + "PileUp_run_077.root").c_str(), "READ");

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

    // AFTER GAIN CHANGE //
    TFile *f_pileup2 = MyTFile((DIR_ROOT_DATA_MATCHED + "PileUp_run_114.root").c_str(), "READ");

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

bool PileUp(Signal SiPM)
{
    double coef = 1;
    if (SiPM.Label > 110)
        coef = 10;
    if (FLAG_GAIN_CHANGED)
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
    for (string NUCLEUS : Nuclide)
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
        int Entry_MAX = 1e6;

        while (Reader->Next() && Reader->GetCurrentEntry() < Entry_MAX)
        {
            ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "Reading Tree");

            ///////////////////////////////////////////////////////////////////////////
            // taking into account the gain chnage of the SiPMs for PileUp rejection //
            if (!(*Silicon)[0].isValid && !(*Silicon)[1].isValid)
            {
                FLAG_GAIN_CHANGED = true;
                continue;
            }
            ///////////////////////////////////////////////////////////////////////////

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
                        // if (PileUp((**SiPM_Groups)[i_group][i_pair].first))
                        //     continue;
                        (**SiPM_Groups)[i_group][i_pair].first.Channel = F_MatchingSiPM[NUCLEUS][GetDetectorChannel((**SiPM_Groups)[i_group][i_pair].first.Label)]->Eval((**SiPM_Groups)[i_group][i_pair].first.Channel);
                        if ((**SiPM_Groups)[i_group][i_pair].first.Channel < eHighLowLimit)
                        {
                            SiPMi = (**SiPM_Groups)[i_group][i_pair].first;
                        }
                    }

                    if ((**SiPM_Groups)[i_group][i_pair].second.isValid && !SiPMi.isValid )
                    {
                        
                        // if (PileUp((**SiPM_Groups)[i_group][i_pair].second))
                        //     continue;
                        
                        SiPMi = (**SiPM_Groups)[i_group][i_pair].second;
                        SiPMi.Channel = F_MatchingLowHigh[NUCLEUS][GetDetectorChannel(SiPMi.Label)]->Eval(SiPMi.Channel);
                        SiPMi.Channel = F_MatchingSiPM[NUCLEUS][GetDetectorChannel(SiPMi.Label)]->Eval(SiPMi.Channel);
                    }

                    if ( (**SiPM_Groups)[i_group][i_pair].second.isValid && (**SiPM_Groups)[i_group][i_pair].first.isValid)
                    {
                        // if (!PileUp((**SiPM_Groups)[i_group][i_pair].first) && !PileUp((**SiPM_Groups)[i_group][i_pair].second))
                        // {
                            H[NUCLEUS][(**SiPM_Groups)[i_group][i_pair].first.Label]->Fill((**SiPM_Groups)[i_group][i_pair].second.Channel, (**SiPM_Groups)[i_group][i_pair].first.Channel);
                        // }
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

void WriteSimulatedTree(vector<string> nucleus_to_recreate)
{
    Info("Write Simulated Tree");
    f_simulated_tree->cd();

    if (find(nucleus_to_recreate.begin(), nucleus_to_recreate.end(), "32Ar") != nucleus_to_recreate.end())
    {
        for (int peak = 0; peak <= 50; peak++)
        {
            if (WindowsMap["32Ar"][peak][11].first == -1 || !WindowsMap["32Ar"][peak][11].first)
                continue;
            Tree_Peaks_Simulated["32Ar"][peak]->Write();
        }
    }
    
    if (find(nucleus_to_recreate.begin(), nucleus_to_recreate.end(), "207Bi") != nucleus_to_recreate.end())
        Tree_Peaks_Simulated["207Bi"][0]->Write();
    if (find(nucleus_to_recreate.begin(), nucleus_to_recreate.end(), "90Sr") != nucleus_to_recreate.end())
        Tree_Peaks_Simulated["90Sr"][0]->Write();
    f_simulated_tree->Close();
    f_simulated_tree = MyTFile((DIR_ROOT_DATA_SIMULATED + "Simulated_SiPM_trees.root").c_str(), "READ");
}

void ExtractingSimulatedTree(vector<string> nucleus_to_recreate)
{
    Info("Extracting Simulated Tree");
    for (string NUCLEUS : nucleus_to_recreate)
    {
        Info("Nucleus: " + NUCLEUS, 1);
        Tree_SIMULATED = (TTree *)SIMULATED_File[NUCLEUS]->Get("PlasticIAS");
        Reader_SIMULATED = new TTreeReader(Tree_SIMULATED);
        Silicon_code = new TTreeReaderValue<int>(*Reader_SIMULATED, "Code");
        Silicon_energy = new TTreeReaderValue<double>(*Reader_SIMULATED, "Energy");
        SiPM_energy = new TTreeReaderValue<double>(*Reader_SIMULATED, "SiPM");

        clock_t start = clock(), Current;
        int Entries = Reader_SIMULATED->GetEntries();
        while (Reader_SIMULATED->Next())
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
    
    WriteSimulatedTree(nucleus_to_recreate);
    
}

void DeletingSimulatedTree(vector<string> nucleus_to_recreate)
{
    Info("Deleting Simulated Tree");
    if (find(nucleus_to_recreate.begin(), nucleus_to_recreate.end(), "32Ar") != nucleus_to_recreate.end())
    {
        Info("Deleting 32Ar");
        for (int peak = 0; peak <= 50; peak++)
        {
            if (WindowsMap["32Ar"][peak][11].first == -1 || !WindowsMap["32Ar"][peak][11].first)
                continue;
            f_simulated_tree->Delete(("Tree_Peaks_Simulated_32Ar_" + to_string(peak)).c_str());
        }
    }
    if (find(nucleus_to_recreate.begin(), nucleus_to_recreate.end(), "207Bi") != nucleus_to_recreate.end())
    {
        Info("Deleting 207Bi");
        f_simulated_tree->Delete(("Tree_Peaks_Simulated_207Bi_" + to_string(0)).c_str());
    }
    if (find(nucleus_to_recreate.begin(), nucleus_to_recreate.end(), "90Sr") != nucleus_to_recreate.end())
    {
        Info("Deleting 90Sr");
        f_simulated_tree->Delete(("Tree_Peaks_Simulated_90Sr_" + to_string(0)).c_str());
    }
}

void InitSimulatedTree(vector<string> nucleus_to_recreate = {})
{
    Info("Init Simulated Tree");
    if (nucleus_to_recreate.empty())
        f_simulated_tree = MyTFile((DIR_ROOT_DATA_SIMULATED + "Simulated_SiPM_trees.root").c_str(), "READ");
    else
    {
        f_simulated_tree = MyTFile((DIR_ROOT_DATA_SIMULATED + "Simulated_SiPM_trees.root").c_str(), "UPDATE");
        DeletingSimulatedTree(nucleus_to_recreate);
    }
    
    

    /// 32Ar ///
    NUCLEUS = "32Ar";
    if (find(nucleus_to_recreate.begin(), nucleus_to_recreate.end(), "32Ar") != nucleus_to_recreate.end())
    {
        for (int peak_number = 0; peak_number < 50; peak_number++)
        {
            if (WindowsMap[NUCLEUS][peak_number][11].first == -1 || !WindowsMap[NUCLEUS][peak_number][11].first)
                continue;
            for (int det = 1; det <= 9; det++)
            {
                f_simulated_tree->cd();
                Tree_Peaks_Simulated[NUCLEUS][peak_number] = new TTree(("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(peak_number)).c_str(), ("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(peak_number)).c_str());
                Tree_Peaks_Simulated[NUCLEUS][peak_number]->Branch("SiPM", &SiPMEnergy);
            }
        }

        for (int det = 1; det <= 9; det++)
        {
            f_simulated_tree->cd();
            Tree_Peaks_Simulated["32Ar"][0] = new TTree(("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(0)).c_str(), ("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(0)).c_str());
            Tree_Peaks_Simulated["32Ar"][0]->Branch("SiPM", &SiPMEnergy);
        }
    }

    /// 207Bi ///
    if (find(nucleus_to_recreate.begin(), nucleus_to_recreate.end(), "207Bi") != nucleus_to_recreate.end())
    {
        NUCLEUS = "207Bi";
        for (int det = 1; det <= 9; det++)
        {
            f_simulated_tree->cd();
            Tree_Peaks_Simulated["207Bi"][0] = new TTree(("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(0)).c_str(), ("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(0)).c_str());
            Tree_Peaks_Simulated["207Bi"][0]->Branch("SiPM", &SiPMEnergy);
        }
    }

    /// 90Sr ///
    if (find(nucleus_to_recreate.begin(), nucleus_to_recreate.end(), "90Sr") != nucleus_to_recreate.end())
    {
        NUCLEUS = "90Sr";
        for (int det = 1; det <= 9; det++)
        {
            f_simulated_tree->cd();
            Tree_Peaks_Simulated["90Sr"][0] = new TTree(("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(0)).c_str(), ("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(0)).c_str());
            Tree_Peaks_Simulated["90Sr"][0]->Branch("SiPM", &SiPMEnergy);
        }
    }

    NUCLEUS = "32Ar";

    ExtractingSimulatedTree(nucleus_to_recreate);
    ReadSimulatedTree();
}

void InitHistograms(int verbose = 0)
{
    Info("Init Histograms");
    for (string Nucleus : Nuclide)
    {
        if (verbose == 1)
            Info("Nucleus: " + Nucleus, 1);
        NUCLEUS = Nucleus;
        if (CALIBRATED_File != nullptr)
            dir[NUCLEUS] = CALIBRATED_File->mkdir(NUCLEUS.c_str());
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorBetaHigh(i))
            {
                H[NUCLEUS][i] = new TH2D(("H_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_" + NUCLEUS + "_" + detectorName[i]).c_str(), eHighN / 2, eHighMin, eHighMax, eLowN / 2, eLowMin, eLowMax);
            }
            if (IsDetectorBetaLow(i))
            {
                H[NUCLEUS][i] = new TH2D(("H_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_" + NUCLEUS + "_" + detectorName[i]).c_str(), eLowN / 100, eLowMin, eLowMax, eLowN / 100, eLowMin, eLowMax);
            }
            if (IsDetectorBetaHigh(i))
            {
                if (verbose == 1)
                    Info("Detector: " + detectorName[i], 2);

                for (int peak = 0; peak <= 50; peak++)
                {
                    if ((WindowsMap[NUCLEUS][peak][11].first == -1 || !WindowsMap[NUCLEUS][peak][11].first))
                        continue;

                    if (i == 101)
                    {
                        dir_peaks[NUCLEUS][peak] = dir[NUCLEUS]->mkdir(("Peak_" + to_string(peak)).c_str());    
                    }

                    if (peak == 14)
                    {
                        dir_SiPM[NUCLEUS][GetDetectorChannel(i)] = dir[NUCLEUS]->mkdir(("SiPM_"+to_string(i)).c_str());
                    }

                    H_SiPM[NUCLEUS][peak][GetDetectorChannel(i)] = new TH1D(("H_SiPM_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), ("H_SiPM_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), eHighN / 10, eHighMin, eHighMax);
                    H_SiPM[NUCLEUS][peak][GetDetectorChannel(i)]->GetXaxis()->SetTitle("Channel");
                    H_SiPM[NUCLEUS][peak][GetDetectorChannel(i)]->GetYaxis()->SetTitle("Counts");
                    H_SiPM[NUCLEUS][peak][GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
                    H_SiPM[NUCLEUS][peak][GetDetectorChannel(i)]->GetYaxis()->CenterTitle();

                    H_SiPM_Calibrated[NUCLEUS][peak][GetDetectorChannel(i)].first = new TH1D(("H_SiPM_Calibrated_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), ("H_SiPM_Calibrated_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), eSiliN_cal / 10, eSiliMin_cal, eSiliMax_cal);
                    H_SiPM_Calibrated[NUCLEUS][peak][GetDetectorChannel(i)].first->GetXaxis()->SetTitle("Energy [keV]");
                    H_SiPM_Calibrated[NUCLEUS][peak][GetDetectorChannel(i)].first->GetYaxis()->SetTitle("Counts");
                    H_SiPM_Calibrated[NUCLEUS][peak][GetDetectorChannel(i)].first->GetXaxis()->CenterTitle();
                    H_SiPM_Calibrated[NUCLEUS][peak][GetDetectorChannel(i)].first->GetYaxis()->CenterTitle();

                    H_SiPM_Calibrated[NUCLEUS][peak][GetDetectorChannel(i)].second = new TH1D(("H_SiPMLow_Calibrated_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), ("H_SiPM_Calibrated_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), eSiliN_cal / 10, eSiliMin_cal, eSiliMax_cal);
                    H_SiPM_Calibrated[NUCLEUS][peak][GetDetectorChannel(i)].second->GetXaxis()->SetTitle("Energy [keV]");
                    H_SiPM_Calibrated[NUCLEUS][peak][GetDetectorChannel(i)].second->GetYaxis()->SetTitle("Counts");
                    H_SiPM_Calibrated[NUCLEUS][peak][GetDetectorChannel(i)].second->GetXaxis()->CenterTitle();
                    H_SiPM_Calibrated[NUCLEUS][peak][GetDetectorChannel(i)].second->GetYaxis()->CenterTitle();
                }

                if (NUCLEUS == "207Bi" || NUCLEUS == "90Sr")
                {
                    dir_peaks[NUCLEUS][0] = NULL;

                    // dir_SiPM[NUCLEUS][GetDetectorChannel(i)] = dir[NUCLEUS]->mkdir(("SiPM_"+to_string(i)).c_str());

                    H_SiPM[NUCLEUS][0][GetDetectorChannel(i)] = new TH1D(("H_SiPM_" + NUCLEUS + "_Peak_" + to_string(0) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), ("H_SiPM_" + NUCLEUS + "_Peak_" + to_string(0) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), eHighN / 10, eHighMin, eHighMax);
                    H_SiPM[NUCLEUS][0][GetDetectorChannel(i)]->GetXaxis()->SetTitle("Channel");
                    H_SiPM[NUCLEUS][0][GetDetectorChannel(i)]->GetYaxis()->SetTitle("Counts");
                    H_SiPM[NUCLEUS][0][GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
                    H_SiPM[NUCLEUS][0][GetDetectorChannel(i)]->GetYaxis()->CenterTitle();

                    H_SiPM_Calibrated[NUCLEUS][0][GetDetectorChannel(i)].first = new TH1D(("H_SiPM_Calibrated_" + NUCLEUS + "_Peak_" + to_string(0) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), ("H_SiPM_Calibrated_" + NUCLEUS + "_Peak_" + to_string(0) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), eSiliN_cal / 10, eSiliMin_cal, eSiliMax_cal);
                    H_SiPM_Calibrated[NUCLEUS][0][GetDetectorChannel(i)].first->GetXaxis()->SetTitle("Energy [keV]");
                    H_SiPM_Calibrated[NUCLEUS][0][GetDetectorChannel(i)].first->GetYaxis()->SetTitle("Counts");
                    H_SiPM_Calibrated[NUCLEUS][0][GetDetectorChannel(i)].first->GetXaxis()->CenterTitle();
                    H_SiPM_Calibrated[NUCLEUS][0][GetDetectorChannel(i)].first->GetYaxis()->CenterTitle();

                    H_SiPM_Calibrated[NUCLEUS][0][GetDetectorChannel(i)].second = new TH1D(("H_SiPMLow_Calibrated_" + NUCLEUS + "_Peak_" + to_string(0) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), ("H_SiPM_Calibrated_" + NUCLEUS + "_Peak_" + to_string(0) + "_SiPM" + to_string(GetDetectorChannel(i))).c_str(), eSiliN_cal / 10, eSiliMin_cal, eSiliMax_cal);
                    H_SiPM_Calibrated[NUCLEUS][0][GetDetectorChannel(i)].second->GetXaxis()->SetTitle("Energy [keV]");
                    H_SiPM_Calibrated[NUCLEUS][0][GetDetectorChannel(i)].second->GetYaxis()->SetTitle("Counts");
                    H_SiPM_Calibrated[NUCLEUS][0][GetDetectorChannel(i)].second->GetXaxis()->CenterTitle();
                    H_SiPM_Calibrated[NUCLEUS][0][GetDetectorChannel(i)].second->GetYaxis()->CenterTitle();
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
    for (string NUCLEUS : Nuclide)
    {
        Info("Nucleus: " + NUCLEUS);
        // if (NUCLEUS == "207Bi")
        //     continue;
        ////////////////////////////////////////////////////////////////
        // Init Histograms
        if (VERBOSE == 1)
            Info("Init Histograms", 1);
        for (int peak_number = 0; peak_number <= 50; peak_number++)
        {
            if ((WindowsMap[NUCLEUS][peak_number][11].first == -1 || !WindowsMap[NUCLEUS][peak_number][11].first) && NUCLEUS == "32Ar")
            {
                continue;
            }
            if (peak_number != 0 && (NUCLEUS == "207Bi" || NUCLEUS == "90Sr"))
            {
                continue;
            }
            // if (peak_number != 14 && NUCLEUS == "32Ar")
            // {
            //     continue;
            // }

            Info("Peak: " + to_string(peak_number), 2);
            
                
            H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Reset();
            H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].first->Reset();
            H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].second->Reset();

            ////////////////////////////////////////////////////////////////

            ////////////////////////////////////////////////////////////////
            // parameters
            if (VERBOSE == 1)
                Info("Setting Parameters", 1);
            double Calibration_OffSet = par[0];
            double Calibration = par[1];
            double Resolution_OffSet = par[2];
            double Resolution_SQRT = par[3];
            double Resolution_2 = par[4];
            double Threshold = SiPM_Threshold_Ratio[NUCLEUS][100+current_detector]*par[5];
            double Threshold_STD = par[6];
            ////////////////////////////////////////////////////////////////

            if (VERBOSE == 1)
                Info("Calibration of Experimental data", 1);

            TTreeReader *Reader_IAS = new TTreeReader(Tree_Peaks[NUCLEUS][peak_number][current_detector]);
            TTreeReaderValue<Signal> *SiPM = new TTreeReaderValue<Signal>(*Reader_IAS, "SiPM");

            Reader_IAS->Restart();
            double energy;
            if (Reader_IAS->GetEntries() == 0) continue;
            while (Reader_IAS->Next())
            {

                energy = Calibration * (**SiPM).Channel / 1000 + Calibration_OffSet;

                if (IsDetectorBetaHigh((**SiPM).Label))
                    H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].first->Fill(energy);
                else
                    H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].second->Fill(energy);
            }

            if (VERBOSE == 1)
                Info("Convolution with detector resolution", 1);
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
            if (Reader_Sim->GetEntries() == 0) continue;

            while (Reader_Sim->Next())
            {
                double th = gRandom->Gaus(Threshold, Threshold_STD);
                double sigma_resolution = sqrt(pow(Resolution_OffSet, 2) + pow(Resolution_SQRT * sqrt(**Energy), 2) + pow(Resolution_2 * pow(**Energy, 2), 2));
                // shoot in a gaussian with ramdom engine
                double energy_conv = gRandom->Gaus(**Energy, sigma_resolution);

                if (energy_conv < th)
                    continue;

                bool flagUntriggered = false;
                for (int sipm = 1; sipm <= 9; sipm++)
                {
                    if (sipm == current_detector)
                        continue;

                    double sigma_resolution_sipm = sqrt(pow(Resolution_OffSet_det[sipm], 2) + pow(Resolution_SQRT_det[sipm] * sqrt(**Energy), 2) + pow(Resolution_2_det[sipm] * pow(**Energy, 2), 2));
                    if (gRandom->Gaus(**Energy, sigma_resolution) < gRandom->Gaus(Threshold, Threshold_STD))
                    {
                        flagUntriggered = true;
                        break;
                    }
                }

                if (flagUntriggered)
                    continue;

                H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Fill(energy_conv);
            }
            ///////////////

            if (VERBOSE == 1)
                Info("Adding threshold", 1);

            ////////// CONVOLUTING THRESHOLD//////////
            // Threshold_f = new TF1("Threshold", "0.5*(1+erf((x-[0])/[1]))", eSiliMin_cal, eSiliMax_cal);
            // Threshold_STD = 1. * sqrt(pow(Resolution_OffSet, 2) + pow(Resolution_SQRT * sqrt(Threshold), 2) + pow(Resolution_2 * pow(Threshold, 2), 2));
            // Threshold_f->SetParameters(Threshold, Threshold_STD);

            // for (int i = 1; i <= H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetNbinsX(); i++)
            // {
            //     H_Sim_Conv[NUCLEUS][peak_number][current_detector]->SetBinContent(i, H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetBinContent(i) * Threshold_f->Eval(H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetBinCenter(i)));
            // }

            // delete Threshold_f;

            //////////////// CHI2. ////////////////
            // Scaling High
            H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].first->GetXaxis()->SetRangeUser(2, 1200);
            H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetXaxis()->SetRangeUser(2, 1200);
            H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].first->Scale(H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Integral() / H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].first->Integral());

            // Scaling Low
            H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].second->GetXaxis()->SetRangeUser(1500, 10000);
            H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetXaxis()->SetRangeUser(1500, 10000);
            H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].second->Scale(H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Integral() / H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].second->Integral());

            // adding Low/High
            H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].first->Add(H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].second);

            H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].first->GetXaxis()->SetRangeUser(0, WindowsBetaMap[NUCLEUS][peak_number]);
            H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetXaxis()->SetRangeUser(0, WindowsBetaMap[NUCLEUS][peak_number]);

            int rebin = Freedman_Diaconis(H_Sim_Conv[NUCLEUS][peak_number][current_detector]);
            H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].first->Rebin(rebin);
            H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Rebin(rebin);

            chi2.push_back(H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].first->Chi2Test(H_Sim_Conv[NUCLEUS][peak_number][current_detector], "CHI2/NDF"));
            // cout << chi2[chi2.size()-1] << "   " << Calibration_OffSet << "   " << Calibration << "   " << Resolution_OffSet << "   " << Resolution_SQRT << "   " << Resolution_2 << "   " << Threshold << "   " << Threshold_STD << endl;

            if (FLAGSAVINGCHI2)
            {
                CHI2[NUCLEUS][peak_number][current_detector] = chi2[chi2.size() - 1];
            }
        }
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
    Minimizer *minimizer = Factory::CreateMinimizer("Minuit2", "Migrad");
    ROOT::Math::Functor functor(&Chi2TreeHist_conv, 7);
    vector<double> Par = {-9.24, 0.7255, 0, Resolution_SQRT_det[current_detector], 2.5e-05, 35, 18};


    for (int i = 0; i <= 10; i++)
    {
        Resolution_OffSet_det.push_back(0);
        Resolution_2_det.push_back(Par[4]);
    }

    minimizer->SetFunction(functor);
    minimizer->SetLimitedVariable(0, "Calibration_OffSet", Par[0], 10, 0, 100);
    minimizer->SetLimitedVariable(1, "Calibration", Par[1], 0.1, 0.5, 1.2);
    minimizer->SetFixedVariable(2, "Resolution_OffSet", Par[2]);
    minimizer->SetLimitedVariable(3, "Resolution_SQRT", Par[3], 0.1, 2., 10);
    minimizer->SetLimitedVariable(4, "Resolution_2", Par[4], 0.1e-5, 1e-6, 8e-5);
    minimizer->SetFixedVariable(5, "Threshold", Par[5]);
    minimizer->SetFixedVariable(6, "Threshold_STD", Par[6]);
    minimizer->SetMaxFunctionCalls(10000000);
    minimizer->SetMaxIterations(100000000);
    minimizer->SetTolerance(1e-3);
    minimizer->SetPrecision(10);
    const double *bestPar = Par.data();
    // minimizer->Minimize();

    // const double *bestPar = minimizer->X();
    FLAGSAVINGCHI2 = true;
    double chi2 = Chi2TreeHist_conv(bestPar);
}

void WriteHistograms()
{
    gStyle->SetOptStat(0);
    CALIBRATED_File->cd();
    Info("Write Histograms");

    // loop on nuclei
    for (string NUCLEUS : Nuclide)
    {
        // Info(NUCLEUS, 1);
        dir[NUCLEUS]->cd();

        TCanvas *cExp_Sim[50];
        TCanvas *cExp_Sim_SiPM_AllPeaks[10];
        for (int peak = 0; peak <= 50; peak++)
        {
            if ((WindowsMap[NUCLEUS][peak][11].first == -1 || !WindowsMap[NUCLEUS][peak][11].first) && NUCLEUS == "32Ar")
                continue;
            
            if (peak != 0 && (NUCLEUS == "207Bi" || NUCLEUS == "90Sr"))
                continue;
            
            cExp_Sim[peak] = new TCanvas(("cExp_Sim_" + NUCLEUS + "_Peak_" + to_string(peak)).c_str(), ("cExp_Sim_" + NUCLEUS + "_Peak_" + to_string(peak)).c_str(), 800, 800);
            cExp_Sim[peak]->Divide(3, 3);
        }

        // loop on SiPM
        for (int det = 1; det <= 9; det++)
        {
            // Info("SiPM: " + to_string(det), 2);
            H[NUCLEUS][100 + det]->Write();
            H[NUCLEUS][110 + det]->Write();
            cExp_Sim_SiPM_AllPeaks[det] = new TCanvas(("cExp_Sim_" + NUCLEUS + "_SiPM" + to_string(det)).c_str(), ("cExp_Sim_" + NUCLEUS + "_SiPM" + to_string(det)).c_str(), 800, 800);
            cExp_Sim_SiPM_AllPeaks[det]->Divide(CanvasMap[NUCLEUS].first, CanvasMap[NUCLEUS].second);
            // loop on peak
            for (int peak = 0; peak <= 50; peak++)
            {

                if ((WindowsMap[NUCLEUS][peak][11].first == -1 || !WindowsMap[NUCLEUS][peak][11].first) && NUCLEUS == "32Ar")
                    continue;

                if (peak != 0 && (NUCLEUS == "207Bi" || NUCLEUS == "90Sr"))
                    continue;

                

                // Info("Peak: " + to_string(peak), 3);
                if (dir_peaks[NUCLEUS][peak] != NULL)
                    dir_peaks[NUCLEUS][peak]->cd();

                H_SiPM_Calibrated[NUCLEUS][peak][current_detector].first->GetXaxis()->SetRangeUser(0, WindowsBetaMap[NUCLEUS][peak]);
                H_Sim_Conv[NUCLEUS][peak][current_detector]->GetXaxis()->SetRangeUser(0, WindowsBetaMap[NUCLEUS][peak]);

                TCanvas *cExp_Sim_det = new TCanvas(("cExp_Sim_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(det)).c_str(), ("cExp_Sim_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(det)).c_str(), 800, 800);
                H_SiPM_Calibrated[NUCLEUS][peak][det].first->Draw("HIST");
                if (H_Sim_Conv[NUCLEUS][peak][det]->GetEntries() > 0)
                {
                    H_Sim_Conv[NUCLEUS][peak][det]->SetLineColor(kRed);
                    H_Sim_Conv[NUCLEUS][peak][det]->Draw("HIST SAME");

                    // write chi2 in txt on the canvas
                    TLatex *latex = new TLatex();
                    latex->SetNDC();
                    latex->SetTextSize(0.04);
                    latex->DrawLatex(0.1, 0.92, ("#chi^{2}_{#nu}: " + to_string(CHI2[NUCLEUS][peak][det])).c_str());
                    latex->Draw("SAME");
                }
                cExp_Sim_det->Write();

                cExp_Sim[peak]->cd(det);
                H_SiPM_Calibrated[NUCLEUS][peak][det].first->Draw("HIST");
                if (H_Sim_Conv[NUCLEUS][peak][det]->GetEntries() > 0)
                {
                    H_Sim_Conv[NUCLEUS][peak][det]->SetLineColor(kRed);
                    H_Sim_Conv[NUCLEUS][peak][det]->Draw("HIST SAME");

                    // write chi2 in txt on the canvas
                    TLatex *latex = new TLatex();
                    latex->SetNDC();
                    latex->SetTextSize(0.04);
                    latex->DrawLatex(0.1, 0.92, ("#chi^{2}_{#nu}: " + to_string(CHI2[NUCLEUS][peak][det])).c_str());
                    latex->Draw("SAME");
                }

                cExp_Sim_SiPM_AllPeaks[det]->cd(peak);
                H_SiPM_Calibrated[NUCLEUS][peak][det].first->Draw("HIST");
                if (H_Sim_Conv[NUCLEUS][peak][det]->GetEntries() > 0)
                {
                    H_Sim_Conv[NUCLEUS][peak][det]->SetLineColor(kRed);
                    H_Sim_Conv[NUCLEUS][peak][det]->Draw("HIST SAME");

                    // write chi2 in txt on the canvas
                    TLatex *latex = new TLatex();
                    latex->SetNDC();
                    latex->SetTextSize(0.04);
                    latex->DrawLatex(0.1, 0.92, ("#chi^{2}_{#nu}: " + to_string(CHI2[NUCLEUS][peak][det])).c_str());
                    latex->Draw("SAME");
                }
            }
        }

        for (int peak = 0; peak <= 50; peak++)
        {
            if ((WindowsMap[NUCLEUS][peak][11].first == -1 || !WindowsMap[NUCLEUS][peak][11].first) && NUCLEUS == "32Ar")
                continue;
            if (peak != 0 && (NUCLEUS == "207Bi" || NUCLEUS == "90Sr"))
                continue;
            if (dir_peaks[NUCLEUS][peak] != NULL)
                dir_peaks[NUCLEUS][peak]->cd();

            cExp_Sim[peak]->Write();
        }

        if (NUCLEUS == "32Ar")
        {
            for (int det = 1; det <= 9; det++)
            {
                dir_SiPM[NUCLEUS][det]->cd();
                cExp_Sim_SiPM_AllPeaks[det]->Write();
            }
        }
    }
}


#endif