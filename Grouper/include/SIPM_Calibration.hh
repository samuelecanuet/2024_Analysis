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
TTreeReaderValue<Signal> *HRS;
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
map<string, pair<double, double>[100][SIGNAL_MAX]> WindowssMap;
// map<string, double> Qbeta; 
// map<string, double[50]> WindowsBetaMap; 
// map<string, pair<int, int>> CanvasMap;
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
map<string, TF1 *[SIGNAL_MAX]> F_MatchingONOFF;
TF1 *f_linear = new TF1("f_linear", "[0]*x + [1]", 0, 10000e3);
TF1 *f_erf = new TF1("f_erf", "x < [0] ? [1]*x + [2] : [3]*(erf((x-[4])/[5]))", 0, 10000e3);
pair<double, double> Range_SiPM_LowHigh = make_pair(50e3, 200e3);
TF1 *gauss;
TF1 *Threshold_f;
int current_detector;
pair<TF1 *, TF1*> PileUp_Function[SIGNAL_MAX];

vector<double> Resolution_OffSet_det;
// vector<double> Resolution_SQRT_det;
vector<double> Resolution_2_det;
map<string, double[SIGNAL_MAX]> SiPM_Threshold_Ratio;
// DATA //
vector<string> Nuclide = {"32Ar", "33Ar"};
double DELTA_PileUp_Range = 200e3;
map<string, vector<int>> Map_Peak;
string NUCLEUS;
double SiliconCalibrationParameter[SIGNAL_MAX][3];
map<string, TF1 *[SIGNAL_MAX]> MatchingLowHigh;
map<string, TF1 *[SIGNAL_MAX]> MatchingLowHigh_erf;
map<string, TF1 *[SIGNAL_MAX]> MatchingSiPM;
map<string, pair<double, double>> SiPM_Range;

vector<double> PreviousPar;
bool FLAGSAVINGCHI2 = false;
map<string, double[50][SIGNAL_MAX]> CHI2;

//try resolution
TGraph *G_Resolution = new TGraph();
TH2D *H_Resolution = new TH2D("H_Resolution", "H_Resolution", 8000, 0, 8000, 2000, 0, 2000);

void InitWindowss(int verbose = 0, string addpath = "", bool ONLYIAS = false)
{
    Info("Init Windows");

    SiPM_Range["32Ar"] = make_pair(0, 6000);
    SiPM_Range["33Ar"] = make_pair(0, 6000);
    SiPM_Range["207Bi"] = make_pair(0, 2000);
    SiPM_Range["90Sr"] = make_pair(0, 4000);

    Qbeta["32Ar"] = 11134-1200;
    Qbeta["33Ar"] = 11134-1200;
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
    Map_Peak["33Ar"] = {21};
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

                if (IAS[nuclei] != number && ONLYIAS)
                    continue;

                // init silicon detector window
                for (int i : Dir2Det(dir, strip))
                {
                    WindowssMap[nuclei][number][i] = make_pair(energy_low, energy_high);
                }

                WindowsBetaMap[nuclei][number] = Qbeta[nuclei] - (energy_low+energy_high)/2;
            }
        }
    }
    Info("Init Windows Done");
}

TF1 *InvertingLinear(TF1 *f)
{
    if (f == nullptr)
        return nullptr;
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
    /* //old
    TFile *f = MyTFile((DIR_ROOT_DATA_MATCHED + "SiPM_Matching_values.root").c_str(), "READ");

    if (f != nullptr)
    {
        if (YEAR == 2024)
        {
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
        else if (YEAR == 2025)
        {
            for (int i = 0; i < SIGNAL_MAX; i++)
            {
                if (IsDetectorBetaHigh(i))
                {
                    F_MatchingSiPM["32Ar"][GetDetectorChannel(i)] = InvertingLinear((TF1 *)f->Get(("MatchingSiPM_066_" + detectorName[i]).c_str()));
                    F_MatchingSiPM["207Bi"][GetDetectorChannel(i)] = InvertingLinear((TF1 *)f->Get(("MatchingSiPM_116_" + detectorName[i]).c_str()));
                    F_MatchingSiPM["90Sr"][GetDetectorChannel(i)] = InvertingLinear((TF1 *)f->Get(("MatchingSiPM_122_" + detectorName[i]).c_str()));
                }

                if (IsDetectorBetaHigh(i))
                {
                    F_MatchingLowHigh["32Ar"][GetDetectorChannel(i)] = (TF1 *)f->Get(("MatchingLowHigh_066_" + detectorName[i]).c_str());
                    F_MatchingLowHigh["207Bi"][GetDetectorChannel(i)] = (TF1 *)f->Get(("MatchingLowHigh_116_" + detectorName[i]).c_str());
                    F_MatchingLowHigh["90Sr"][GetDetectorChannel(i)] = (TF1 *)f->Get(("MatchingLowHigh_122_" + detectorName[i]).c_str());
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorBetaHigh(i))
            {
                F_MatchingSiPM["32Ar"][GetDetectorChannel(i)] = new TF1("dummy","x",0,eHighMax);
                F_MatchingSiPM["207Bi"][GetDetectorChannel(i)] = new TF1("dummy","x",0,eHighMax);
                F_MatchingSiPM["90Sr"][GetDetectorChannel(i)] = new TF1("dummy","x",0,eHighMax);
            }
            if (IsDetectorBetaHigh(i))
            {
                F_MatchingLowHigh["32Ar"][GetDetectorChannel(i)] = new TF1("dummy","10*x",0,eHighMax);
                F_MatchingLowHigh["207Bi"][GetDetectorChannel(i)] = new TF1("dummy","10*x",0,eHighMax);
                F_MatchingLowHigh["90Sr"][GetDetectorChannel(i)] = new TF1("dummy","10*x",0,eHighMax);
            }
        }
    }
    */

   TFile *f = MyTFile((DIR_ROOT_DATA_MATCHED + "SiPM_Matching_Functions.root").c_str(), "READ");

    if (f != nullptr)
    {
        for (string nucleus : Nuclide)
        {
            if (YEAR == 2024 || YEAR == 2025)
            {
                for (int i = 0; i < SIGNAL_MAX; i++)
                {
                    if (IsDetectorBeta(i))
                    {
                        F_MatchingSiPM[nucleus][GetDetectorChannel(i)] = InvertFunction((TF1 *)f->Get(("SiPM_Gain_" + detectorName[i]).c_str()));
                    }

                    if (IsDetectorBetaLow(i))
                    {
                        F_MatchingLowHigh[nucleus][GetDetectorChannel(i)] = (TF1 *)f->Get(("SiPM_HighLow_"+ detectorName[i]).c_str());
                    }

                    if (IsDetectorBeta(i))
                    {
                        F_MatchingONOFF[nucleus][GetDetectorChannel(i)] = ((TF1 *)f->Get(("SiPM_ONOFF_" + detectorName[i]).c_str()));
                    }
                }
            }
        }
    }
    else
    {
        Error("No SiPM Matching file found");
    }


}

void InitSiliconCalibration(string addpath = "")
{

    string CalibFileName;

    CalibFileName = addpath + "Config_Files/" + to_string(YEAR) + "/Calibration_" + to_string(YEAR) + ".txt";

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
    else if (YEAR == 2025)
    {
        // // AFTER GAIN CHANGE //
        TFile *f_pileup2 = MyTFile((DIR_ROOT_DATA_MATCHED + "PileUp_run_122.root").c_str(), "READ");

        if (f_pileup2 != nullptr)
        {
            for (int det = 0; det < SIGNAL_MAX  ; det++)
            {
                if (IsDetectorBeta(det))
                {
                    PileUp_Function[det].first = (TF1 *)f_pileup2->Get((detectorName[det]+"/PileUp_Correction_"+detectorName[det]).c_str());
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

    for (string nucleus : Nuclide)
    {
        for (int peak = 0; peak <= 50; peak++)
        {
            for (int det = 1; det <= 9; det++)
            {
                if (WindowssMap[nucleus][peak][11].first == -1 || !WindowssMap[nucleus][peak][11].first)
                    continue;
                Tree_Peaks[nucleus][peak][det] = (TTree *)f_tree->Get(("Tree_Peaks_" + nucleus + "_" + to_string(peak) + "_" + to_string(det)).c_str());
                if (Tree_Peaks[nucleus][peak][det] == nullptr)
                    Warning("Experimental TTree for " + nucleus + "peak n°" + to_string(peak) + " SIPM" + to_string(det) + " not found");
            }
        }
    }
    // for (int det = 1; det <= 9; det++)
    // {
    //     Tree_Peaks["207Bi"][0][det] = (TTree *)f_tree->Get(("Tree_Peaks_207Bi_" + to_string(0) + "_" + to_string(det)).c_str());
 
    //     Tree_Peaks["90Sr"][0][det] = (TTree *)f_tree->Get(("Tree_Peaks_90Sr_" + to_string(0) + "_" + to_string(det)).c_str());
    // }
    G_Resolution = (TGraph *)f_tree->Get("G_Resolution");
    H_Resolution = (TH2D *)f_tree->Get("H_Resolution");
}

void WriteTree()
{
    Info("Write Tree");
    f_tree->cd();

    for (string nucleus : Nuclide)
    {
        Info("Nucleus: " + nucleus, 1);
        for (int peak = 0; peak <= 50; peak++)
        {
            for (int det = 1; det <= 9; det++)
            {
                if (WindowssMap[nucleus][peak][11].first == -1 || !WindowssMap[nucleus][peak][11].first)
                    continue;
                Tree_Peaks[nucleus][peak][det]->Write();
            }
        }
    }
    // for (int det = 1; det <= 9; det++)
    // {
    //     Tree_Peaks["207Bi"][0][det]->Write();
    //     Tree_Peaks["90Sr"][0][det]->Write();
    // }
    G_Resolution->SetName("G_Resolution");
    G_Resolution->Write();
    H_Resolution->Write();
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

        if (NUCLEUS == "32Ar" || NUCLEUS == "33Ar")
        {
            Tree = (TTree *)GROUPED_File[NUCLEUS]->Get("MERGED_Tree");
            Reader = new TTreeReader(Tree);
            Silicon = new TTreeReaderArray<Signal>(*Reader, "MERGED_Tree_Silicon");
            SiPM_Groups = new TTreeReaderValue<vector<vector<pair<Signal, Signal>>>>(*Reader, "MERGED_Tree_SiPMGroup");
            HRS = new TTreeReaderValue<Signal>(*Reader, "MERGED_Tree_HRS");
        }
        else
        {
            Tree = (TTree *)GROUPED_File[NUCLEUS]->Get("CLEANED_Tree");
            Reader = new TTreeReader(Tree);
            SiPM_Groups = new TTreeReaderValue<vector<vector<pair<Signal, Signal>>>>(*Reader, "CLEANED_Tree_SiPMGroup");
            // HRS = new TTreeReaderValue<Signal>(*Reader, "CLEANED_Tree_HRS");
        }

        clock_t start = clock(), Current;
        int Entries = Reader->GetEntries();
        double Entry_MAX = Entries;

        while (Reader->Next() && Reader->GetCurrentEntry() < Entry_MAX)
        {
            ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "Reading Tree");

            ///////////////////////////////////////////////////////////////////////////
            // taking into account the gain chnage of the SiPMs for PileUp rejection //
            // if (!(*Silicon)[0].isValid && !(*Silicon)[1].isValid)
            // {
            //     FLAG_GAIN_CHANGED = true;
            //     continue;
            // }
            ///////////////////////////////////////////////////////////////////////////

            double energy;
            int sili_code;
            int peak = 0;
            if (NUCLEUS == "32Ar" || NUCLEUS == "33Ar") ///// ### IF MORE THAN 1 PEAK IN THE WINDOW THE LOWEST IN ENERGY IS TAKEN
            {
                
                // cout << **HRS << endl;
                if (IsHRSProton((**HRS).Label))
                    continue;


                energy = F_SiliconCalibration[(*Silicon)[1].Label]->Eval((*Silicon)[1].Channel / 1000.);
                sili_code = (*Silicon)[1].Label;

                for (int peak_number = 0; peak_number < 50; peak_number++)
                {
                    SiPMi = Signal();
                    if (WindowssMap[NUCLEUS][peak_number][sili_code].first == -1 || !WindowssMap[NUCLEUS][peak_number][sili_code].first)
                        continue;

                    if (energy > WindowssMap[NUCLEUS][peak_number][sili_code].first && energy < WindowssMap[NUCLEUS][peak_number][sili_code].second)
                    {
                        peak = peak_number;
                        break;
                    }
                }
            }

            

            // cout << "Peak: " << peak << endl; 
            

            for (int i_group = 0; i_group < (**SiPM_Groups).size(); i_group++)
            {       
                // cout << "Group: " << i_group << " / " << (**SiPM_Groups).size() << endl;

                if ((**SiPM_Groups)[i_group].size() <= 2)
                    continue;
                for (int i_pair = 0; i_pair < (**SiPM_Groups)[i_group].size(); i_pair++)
                {

                    // cout << "     Pair: " << i_pair << endl;

                    SiPMi = Signal();

                    /*
                    ////  ## Selecting High Low with a limit ## ///
                    // taking high is valid
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

                    // - taking low if high not valid
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

                    */

                    ////  ## Collecting All High Low ## ///
                    /////
                    // - taking low if valid
                    // cout << "     Pair: " << i_pair << endl;
                    if ((**SiPM_Groups)[i_group][i_pair].first.isValid)
                    {
                        if (abs((**SiPM_Groups)[i_group][i_pair].first.Time) > 50)
                            continue;
                        // if (PileUp((**SiPM_Groups)[i_group][i_pair].first))
                        //     continue;
                        SiPMi = (**SiPM_Groups)[i_group][i_pair].first;
                        Tree_Peaks[NUCLEUS][peak][GetDetectorChannel(SiPMi.Label)]->Fill();
                    }
                    if ((**SiPM_Groups)[i_group][i_pair].second.isValid)
                    {
                        if (abs((**SiPM_Groups)[i_group][i_pair].second.Time) > 50)
                            continue;
                        // if (PileUp((**SiPM_Groups)[i_group][i_pair].first))
                        //     continue;
                        SiPMi = (**SiPM_Groups)[i_group][i_pair].second;
                        Tree_Peaks[NUCLEUS][peak][GetDetectorChannel(SiPMi.Label)]->Fill();
                    }
                }


                // try to characterise the resolution
                double n = 8;
                if ((**SiPM_Groups)[i_group].size() >= n && peak == 14 && NUCLEUS == "32Ar")
                {
                    vector<double> energies;
                    // in 2025 exluding sipm6
                    for (int i = 0; i < (**SiPM_Groups)[i_group].size(); i++)
                    {
                        if ((**SiPM_Groups)[i_group][i].second.Label == 116) // sipm6
                            continue;

                        // if ((**SiPM_Groups)[i_group][i].first.Label == 106) // sipm6
                        //     continue;
                        // if ((**SiPM_Groups)[i_group][i].first.Label == 107) // sipm7
                        //     continue;
                        // if ((**SiPM_Groups)[i_group][i].first.Label == 108) // sipm8
                        //     continue;
                        
                        if ((**SiPM_Groups)[i_group][i].second.isValid && (**SiPM_Groups)[i_group][i].second.Channel > 0)
                            energies.push_back((**SiPM_Groups)[i_group][i].second.Channel/1000);
                    }

                    if (energies.size() == n)
                    {
                        double mean = 0;
                        for (auto e : energies)
                            mean += e;
                        mean /= energies.size();

                        double variance = 0;
                        for (auto e : energies)
                            variance += (e - mean) * (e - mean);
                        variance /= n;
                        double sigma = sqrt(variance) / sqrt(n);

                        G_Resolution->AddPoint(mean, sigma);
                        H_Resolution->Fill(mean, sigma);
                    }
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
    

    if (option == "RECREATE")
    {

        NUCLEUS = "32Ar";
        for (int peak_number = 0; peak_number < 50; peak_number++)
        {
            SiPMi = Signal();
            if (WindowssMap[NUCLEUS][peak_number][11].first == -1 || !WindowssMap[NUCLEUS][peak_number][11].first)
                continue;
            for (int det = 1; det <= 9; det++)
            {
                f_tree->cd();
                Tree_Peaks[NUCLEUS][peak_number][det] = new TTree(("Tree_Peaks_" + NUCLEUS + "_" + to_string(peak_number) + "_" + to_string(det)).c_str(), ("Tree_Peaks_" + NUCLEUS + "_" + to_string(peak_number) + "_" + to_string(det)).c_str());
                Tree_Peaks[NUCLEUS][peak_number][det]->Branch("SiPM", &SiPMi);
            }
        }
        for (int det = 1; det <= 9; det++)
        {
            f_tree->cd();
            Tree_Peaks[NUCLEUS][0][det] = new TTree(("Tree_Peaks_" + NUCLEUS + "_" + to_string(0) + "_" + to_string(det)).c_str(), ("Tree_Peaks_" + NUCLEUS + "_" + to_string(0) + "_" + to_string(det)).c_str());
            Tree_Peaks[NUCLEUS][0][det]->Branch("SiPM", &SiPMi);
        }

        NUCLEUS = "33Ar";
        for (int peak_number = 0; peak_number < 50; peak_number++)
        {
            SiPMi = Signal();
            if (WindowssMap[NUCLEUS][peak_number][11].first == -1 || !WindowssMap[NUCLEUS][peak_number][11].first)
                continue;
            for (int det = 1; det <= 9; det++)
            {
                f_tree->cd();
                Tree_Peaks[NUCLEUS][peak_number][det] = new TTree(("Tree_Peaks_" + NUCLEUS + "_" + to_string(peak_number) + "_" + to_string(det)).c_str(), ("Tree_Peaks_" + NUCLEUS + "_" + to_string(peak_number) + "_" + to_string(det)).c_str());
                Tree_Peaks[NUCLEUS][peak_number][det]->Branch("SiPM", &SiPMi);
            }
        }
        for (int det = 1; det <= 9; det++)
        {
            f_tree->cd();
            Tree_Peaks[NUCLEUS][0][det] = new TTree(("Tree_Peaks_" + NUCLEUS + "_" + to_string(0) + "_" + to_string(det)).c_str(), ("Tree_Peaks_" + NUCLEUS + "_" + to_string(0) + "_" + to_string(det)).c_str());
            Tree_Peaks[NUCLEUS][0][det]->Branch("SiPM", &SiPMi);
        }

        // init trees
        // NUCLEUS = "207Bi";
        // for (int det = 1; det <= 9; det++)
        // {
        //     f_tree->cd();
        //     Tree_Peaks["207Bi"][0][det] = new TTree(("Tree_Peaks_" + NUCLEUS + "_" + to_string(0) + "_" + to_string(det)).c_str(), ("Tree_Peaks_" + NUCLEUS + "_" + to_string(0) + "_" + to_string(det)).c_str());
        //     if (Tree_Peaks["207Bi"][0][det] == NULL)
        //     {
        //         Error("Tree_Peaks_207Bi_" + to_string(0) + "_" + to_string(det) + " not created");
        //     }
        //     Tree_Peaks["207Bi"][0][det]->Branch("SiPM", &SiPMi);
        // }

        // NUCLEUS = "90Sr";
        // for (int det = 1; det <= 9; det++)
        // {
        //     f_tree->cd();
        //     Tree_Peaks["90Sr"][0][det] = new TTree(("Tree_Peaks_" + NUCLEUS + "_" + to_string(0) + "_" + to_string(det)).c_str(), ("Tree_Peaks_" + NUCLEUS + "_" + to_string(0) + "_" + to_string(det)).c_str());
        //     Tree_Peaks["90Sr"][0][det]->Branch("SiPM", &SiPMi);
        // }
        
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
    TTreeReader *Reader;
    TTreeReaderValue<double> *energy;

    if (find(Nuclide.begin(), Nuclide.end(), "32Ar") != Nuclide.end())
    {
        Info("32Ar", 1);
        for (int peak = 0; peak <= 50; peak++)
        {
            if (WindowssMap["32Ar"][peak][11].first == -1 || !WindowssMap["32Ar"][peak][11].first)
                continue;
            Tree_Peaks_Simulated["32Ar"][peak] = (TTree *)f_simulated_tree->Get(("Tree_Peaks_Simulated_32Ar_" + to_string(peak)).c_str());
            if (Tree_Peaks_Simulated["32Ar"][peak] == nullptr) continue;

            Reader = new TTreeReader(Tree_Peaks_Simulated["32Ar"][peak]);
            energy = new TTreeReaderValue<double>(*Reader, "SiPM");
            while (Reader->Next())
            {
                double e = **energy;
                H_Sim["32Ar"][peak]->Fill(e);
            }
        }
    }

    if (find(Nuclide.begin(), Nuclide.end(), "33Ar") != Nuclide.end())
    {
        Info("33Ar", 1);
        for (int peak = 0; peak <= 50; peak++)
        {
            if (WindowssMap["33Ar"][peak][11].first == -1 || !WindowssMap["33Ar"][peak][11].first)
                continue;
            Tree_Peaks_Simulated["33Ar"][peak] = (TTree *)f_simulated_tree->Get(("Tree_Peaks_Simulated_33Ar_" + to_string(peak)).c_str());
            if (Tree_Peaks_Simulated["33Ar"][peak] == nullptr) continue;

            Reader = new TTreeReader(Tree_Peaks_Simulated["33Ar"][peak]);
            energy = new TTreeReaderValue<double>(*Reader, "SiPM");
            while (Reader->Next())
            {
                double e = **energy;
                H_Sim["33Ar"][peak]->Fill(e);
            }
        }
    }

    // if (find(Nuclide.begin(), Nuclide.end(), "207Bi") != Nuclide.end())
    // {
    //     Tree_Peaks_Simulated["207Bi"][0] = (TTree *)f_simulated_tree->Get(("Tree_Peaks_Simulated_207Bi_" + to_string(0)).c_str());
    //     if (Tree_Peaks_Simulated["207Bi"][0] != nullptr) 
    //     {
    //         Info("Reading 207Bi Simulated Tree");
    //         Reader = new TTreeReader(Tree_Peaks_Simulated["207Bi"][0]);
    //         energy = new TTreeReaderValue<double>(*Reader, "SiPM");
    //         while (Reader->Next())
    //         {
    //             double e = **energy;
    //             H_Sim["207Bi"][0]->Fill(e);
    //         }
    //     }
    // }
    
    // if (find(Nuclide.begin(), Nuclide.end(), "90Sr") != Nuclide.end())
    // {
    //     Tree_Peaks_Simulated["90Sr"][0] = (TTree *)f_simulated_tree->Get(("Tree_Peaks_Simulated_90Sr_" + to_string(0)).c_str());
    //     if (Tree_Peaks_Simulated["90Sr"][0] != nullptr)
    //     {
    //         Info("Reading 90Sr simulated tree");
    //         Reader = new TTreeReader(Tree_Peaks_Simulated["90Sr"][0]);
    //         energy = new TTreeReaderValue<double>(*Reader, "SiPM");
    //         while (Reader->Next())
    //         {
    //             double e = **energy;
    //             H_Sim["90Sr"][0]->Fill(e);
    //         }
    //     }
    // }
    Info("Read Simulated Tree Done");
}

void WriteSimulatedTree(vector<string> nucleus_to_recreate)
{
    Info("Write Simulated Tree");
    f_simulated_tree->cd();

    if (find(nucleus_to_recreate.begin(), nucleus_to_recreate.end(), "32Ar") != nucleus_to_recreate.end())
    {
        for (int peak = 0; peak <= 50; peak++)
        {
            if (WindowssMap["32Ar"][peak][11].first == -1 || !WindowssMap["32Ar"][peak][11].first)
                continue;
            Tree_Peaks_Simulated["32Ar"][peak]->Write();
        }
    }

    if (find(nucleus_to_recreate.begin(), nucleus_to_recreate.end(), "33Ar") != nucleus_to_recreate.end())
    {
        for (int peak = 0; peak <= 50; peak++)
        {
            if (WindowssMap["33Ar"][peak][11].first == -1 || !WindowssMap["33Ar"][peak][11].first)
                continue;
            Tree_Peaks_Simulated["33Ar"][peak]->Write();
        }
    }
    
    // if (find(nucleus_to_recreate.begin(), nucleus_to_recreate.end(), "207Bi") != nucleus_to_recreate.end())
    //     Tree_Peaks_Simulated["207Bi"][0]->Write();
    // if (find(nucleus_to_recreate.begin(), nucleus_to_recreate.end(), "90Sr") != nucleus_to_recreate.end())
    //     Tree_Peaks_Simulated["90Sr"][0]->Write();
    f_simulated_tree->Close();
    f_simulated_tree = MyTFile((DIR_ROOT_DATA_SIMULATED + "Simulated_SiPM_trees.root").c_str(), "READ");
}

void ExtractingSimulatedTree(vector<string> nucleus_to_recreate)
{
    Info("Extracting Simulated Tree");
    for (string NUCLEUS : nucleus_to_recreate)
    {
        if (NUCLEUS == "32Ar" || NUCLEUS == "33Ar")
        {

            Info("Nucleus: " + NUCLEUS, 1);
            Tree_SIMULATED = (TTree *)SIMULATED_File[NUCLEUS]->Get("PlasticCoinc");
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

                for (int peak_number = 0; peak_number < 50; peak_number++)
                {
                    if (WindowssMap[NUCLEUS][peak_number][sili_code].first == -1 || !WindowssMap[NUCLEUS][peak_number][sili_code].first)
                        continue;

                    if (sili_e > WindowssMap[NUCLEUS][peak_number][sili_code].first & sili_e < WindowssMap[NUCLEUS][peak_number][sili_code].second)
                    {
                        SiPMEnergy = SiPM_e;
                        Tree_Peaks_Simulated[NUCLEUS][peak_number]->Fill();
                    }
                }
            }
        }
        else
        {
            Info("Nucleus: " + NUCLEUS, 1);
            Tree_SIMULATED = (TTree *)SIMULATED_File[NUCLEUS]->Get("Plastic");
            Reader_SIMULATED = new TTreeReader(Tree_SIMULATED);
            SiPM_energy = new TTreeReaderValue<double>(*Reader_SIMULATED, "SiPM");

            clock_t start = clock(), Current;
            int Entries = Reader_SIMULATED->GetEntries();
            while (Reader_SIMULATED->Next())
            {
                ProgressBar(Reader_SIMULATED->GetCurrentEntry(), Entries, start, Current, "Reading Tree");

                double SiPM_e = **SiPM_energy;

                SiPMEnergy = SiPM_e;
                Tree_Peaks_Simulated[NUCLEUS][0]->Fill();
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
            if (WindowssMap["32Ar"][peak][11].first == -1 || !WindowssMap["32Ar"][peak][11].first)
                continue;
            f_simulated_tree->Delete(("Tree_Peaks_Simulated_32Ar_" + to_string(peak)).c_str());
        }
    }
    if (find(nucleus_to_recreate.begin(), nucleus_to_recreate.end(), "33Ar") != nucleus_to_recreate.end())
    {
        Info("Deleting 33Ar");
        for (int peak = 0; peak <= 50; peak++)
        {
            if (WindowssMap["33Ar"][peak][11].first == -1 || !WindowssMap["33Ar"][peak][11].first)
                continue;
            f_simulated_tree->Delete(("Tree_Peaks_Simulated_33Ar_" + to_string(peak)).c_str());
        }
    }
    // if (find(nucleus_to_recreate.begin(), nucleus_to_recreate.end(), "207Bi") != nucleus_to_recreate.end())
    // {
    //     Info("Deleting 207Bi");
    //     f_simulated_tree->Delete(("Tree_Peaks_Simulated_207Bi_" + to_string(0)).c_str());
    // }
    // if (find(nucleus_to_recreate.begin(), nucleus_to_recreate.end(), "90Sr") != nucleus_to_recreate.end())
    // {
    //     Info("Deleting 90Sr");
    //     f_simulated_tree->Delete(("Tree_Peaks_Simulated_90Sr_" + to_string(0)).c_str());
    // }
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
            if (WindowssMap[NUCLEUS][peak_number][11].first == -1 || !WindowssMap[NUCLEUS][peak_number][11].first)
                continue;
            for (int det = 1; det <= 9; det++)
            {
                f_simulated_tree->cd();
                Tree_Peaks_Simulated[NUCLEUS][peak_number] = new TTree(("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(peak_number)).c_str(), ("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(peak_number)).c_str());
                if (Tree_Peaks_Simulated[NUCLEUS][peak_number] == nullptr) continue;                    
                Tree_Peaks_Simulated[NUCLEUS][peak_number]->Branch("SiPM", &SiPMEnergy);              
            }
        }

        for (int det = 1; det <= 9; det++)
        {
            f_simulated_tree->cd();
            Tree_Peaks_Simulated["32Ar"][0] = new TTree(("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(0)).c_str(), ("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(0)).c_str());
            if (Tree_Peaks_Simulated["32Ar"][0]) continue;  
            Tree_Peaks_Simulated["32Ar"][0]->Branch("SiPM", &SiPMEnergy);
        }
    }

    NUCLEUS = "33Ar";
    if (find(nucleus_to_recreate.begin(), nucleus_to_recreate.end(), "33Ar") != nucleus_to_recreate.end())
    {
        for (int peak_number = 0; peak_number < 50; peak_number++)
        {
            if (WindowssMap[NUCLEUS][peak_number][11].first == -1 || !WindowssMap[NUCLEUS][peak_number][11].first)
                continue;
            for (int det = 1; det <= 9; det++)
            {
                f_simulated_tree->cd();
                Tree_Peaks_Simulated[NUCLEUS][peak_number] = new TTree(("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(peak_number)).c_str(), ("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(peak_number)).c_str());
                if (Tree_Peaks_Simulated[NUCLEUS][peak_number] == nullptr) continue;                    
                Tree_Peaks_Simulated[NUCLEUS][peak_number]->Branch("SiPM", &SiPMEnergy);              
            }
        }

        for (int det = 1; det <= 9; det++)
        {
            f_simulated_tree->cd();
            Tree_Peaks_Simulated["33Ar"][0] = new TTree(("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(0)).c_str(), ("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(0)).c_str());
            if (Tree_Peaks_Simulated["33Ar"][0]) continue;  
            Tree_Peaks_Simulated["33Ar"][0]->Branch("SiPM", &SiPMEnergy);
        }
    }

    /// 207Bi ///
    // if (find(nucleus_to_recreate.begin(), nucleus_to_recreate.end(), "207Bi") != nucleus_to_recreate.end())
    // {
    //     NUCLEUS = "207Bi";
    //     for (int det = 1; det <= 9; det++)
    //     {
    //         f_simulated_tree->cd();
    //         Tree_Peaks_Simulated["207Bi"][0] = new TTree(("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(0)).c_str(), ("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(0)).c_str());
    //         Tree_Peaks_Simulated["207Bi"][0]->Branch("SiPM", &SiPMEnergy);
    //     }
    // }

    // /// 90Sr ///
    // if (find(nucleus_to_recreate.begin(), nucleus_to_recreate.end(), "90Sr") != nucleus_to_recreate.end())
    // {
    //     NUCLEUS = "90Sr";
    //     for (int det = 1; det <= 9; det++)
    //     {
    //         f_simulated_tree->cd();
    //         Tree_Peaks_Simulated["90Sr"][0] = new TTree(("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(0)).c_str(), ("Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(0)).c_str());
    //         Tree_Peaks_Simulated["90Sr"][0]->Branch("SiPM", &SiPMEnergy);
    //     }
    // }

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
                    if ((WindowssMap[NUCLEUS][peak][11].first == -1 || !WindowssMap[NUCLEUS][peak][11].first))
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
            if (WindowssMap[NUCLEUS][peak][11].first == -1 || !WindowssMap[NUCLEUS][peak][11].first)
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
        // Info("Nucleus: " + NUCLEUS);
        // if (NUCLEUS == "207Bi")
        //     continue;
        ////////////////////////////////////////////////////////////////
        // Init Histograms
        if (VERBOSE == 1)
            Info("Init Histograms", 1);

        // FIRST ias
        vector<double> PeakVector = {IAS[NUCLEUS]};
        for (int i=0;i<IAS[NUCLEUS];++i) PeakVector.push_back(i);
        for (int i=IAS[NUCLEUS]+1;i<=50;++i) PeakVector.push_back(i); 
        //
        for (int peak_number : PeakVector)
        {
            if ((WindowssMap[NUCLEUS][peak_number][11].first == -1 || !WindowssMap[NUCLEUS][peak_number][11].first) && (NUCLEUS == "32Ar" || NUCLEUS == "33Ar"))
            {
                continue;
            }
            if (peak_number != 0 && (NUCLEUS == "207Bi" || NUCLEUS == "90Sr"))
            {
                continue;
            }

            // Info("Peak: " + to_string(peak_number), 2);        

            ////////////////////////////////////////////////////////////////

            ////////////////////////////////////////////////////////////////
            // parameters
            if (VERBOSE == 1)  Info("Setting Parameters", 3);
            double Calibration_OffSet = par[0];
            double Calibration = par[1];
            double Resolution_OffSet = par[2];
            double Resolution_SQRT = par[3];
            double Resolution_2 = par[4];
            // double Threshold = SiPM_Threshold_Ratio[NUCLEUS][100+current_detector]*par[5];
            // double Threshold_STD = par[6];
            ////////////////////////////////////////////////////////////////
            //rebin
            double rebin = 8;
            if (IAS[NUCLEUS] != peak_number)
                rebin = Freedman_Diaconis(H_Sim[NUCLEUS][peak_number]);
            if (H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetNbinsX() == 10000)
            {
                H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].first->Rebin(rebin);
                H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].second->Rebin(rebin);
                H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Rebin(rebin);
            }
            

            // ## CALIBRATION ## //
            if (VERBOSE == 1)
                Info("Calibration of Experimental data", 3);
            if (Tree_Peaks[NUCLEUS][peak_number][current_detector] == nullptr) Error("Experimental Tree is nullptr");
            TTreeReader *Reader_IAS = new TTreeReader(Tree_Peaks[NUCLEUS][peak_number][current_detector]);
            TTreeReaderValue<Signal> *SiPM = new TTreeReaderValue<Signal>(*Reader_IAS, "SiPM");
            Reader_IAS->Restart();
            double energy;
            if (Reader_IAS->GetEntries() == 0) 
            {
                Warning("No entries in Tree_Peaks_" + NUCLEUS + "_" + to_string(peak_number) + "_" + to_string(current_detector));
                continue;
            }

            if (PreviousPar[0] != Calibration_OffSet || PreviousPar[1] != Calibration || H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].first->GetEntries() == 0 || H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].second->GetEntries() == 0)
            {
                PreviousPar[0] = Calibration_OffSet;
                PreviousPar[1] = Calibration;
                H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].first->Reset();
                H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].second->Reset();
                while (Reader_IAS->Next() && Reader_IAS->GetCurrentEntry() < 1e6)
                {
                    energy = Calibration * (**SiPM).Channel / 1000 + Calibration_OffSet;

                    if (IsDetectorBetaHigh((**SiPM).Label))
                        H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].first->Fill(energy);
                    else
                        H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].second->Fill(energy);
                }
            }

            // ## CONVOLUTION ## //
            if (VERBOSE == 1)
                Info("Convolution with detector resolution", 3);
            // convoluting resolution on histogram
            TTreeReader *Reader_Sim = new TTreeReader(Tree_Peaks_Simulated[NUCLEUS][peak_number]);
            TTreeReaderValue<double> *Energy = new TTreeReaderValue<double>(*Reader_Sim, "SiPM");
            if (Reader_Sim->GetEntries() == 0) 
            {
                Warning("No entries in Tree_Peaks_Simulated_" + NUCLEUS + "_" + to_string(peak_number));
                continue;
            }

            if (PreviousPar[2] != Resolution_OffSet || PreviousPar[3] != Resolution_SQRT || PreviousPar[4] != Resolution_2 || H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetEntries() == 0)
            {
                int bin = H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetNbinsX();
                PreviousPar[2] = Resolution_OffSet;
                PreviousPar[3] = Resolution_SQRT;
                PreviousPar[4] = Resolution_2;
                H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Reset();
           
                for (int i = 1; i <= H_Sim[NUCLEUS][peak_number]->GetNbinsX(); i++)
                {
                    if (H_Sim[NUCLEUS][peak_number]->GetBinContent(i) == 0)
                        continue;
                    energy = H_Sim[NUCLEUS][peak_number]->GetBinCenter(i);
                    double sigma_resolution = sqrt(pow(Resolution_OffSet, 2) + pow(Resolution_SQRT * sqrt(energy), 2) + pow(Resolution_2 * pow(energy, 2), 2));
                    // double sigma_resolution = Resolution_OffSet + Resolution_SQRT * sqrt(energy) + Resolution_2 * pow(energy, 2);        

                    // double sigma_resolution = Calibration *(67.20 + 0.115566 * (energy/Calibration) + 3.70031e-6 * pow((energy/Calibration),2)); // FWHM from 32Ar SiPM calibration
                    // sigma_resolution += Resolution_SQRT * sqrt(energy);

                    gauss = new TF1("gauss", "gaus", 0, eSiliMax_cal);
                    gauss->SetNpx(bin);
                    gauss->SetParameters(H_Sim[NUCLEUS][peak_number]->GetBinContent(i) / (sigma_resolution * sqrt(2 * M_PI)), energy, sigma_resolution); /// weighted gaussian

                    H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Add(gauss->GetHistogram());

                    delete gauss;
                }
            }

            ///////////////
            // random convoluting resolution on histogram
            // TTreeReader *Reader_Sim = new TTreeReader(Tree_Peaks_Simulated[NUCLEUS][peak_number]);
            // TTreeReaderValue<double> *Energy = new TTreeReaderValue<double>(*Reader_Sim, "SiPM");
            // if (Reader_Sim->GetEntries() == 0) continue;

            // while (Reader_Sim->Next())
            // {
            //     double th = gRandom->Gaus(Threshold, Threshold_STD);
            //     double sigma_resolution = sqrt(pow(Resolution_OffSet, 2) + pow(Resolution_SQRT * sqrt(**Energy), 2) + pow(Resolution_2 * pow(**Energy, 2), 2));
            //     // shoot in a gaussian with ramdom engine
            //     double energy_conv = gRandom->Gaus(**Energy, sigma_resolution);

            //     if (energy_conv < th)
            //         continue;

            //     bool flagUntriggered = false;
            //     for (int sipm = 1; sipm <= 9; sipm++)
            //     {
            //         if (sipm == current_detector)
            //             continue;

            //         double sigma_resolution_sipm = sqrt(pow(Resolution_OffSet_det[sipm], 2) + pow(Resolution_SQRT_det[sipm] * sqrt(**Energy), 2) + pow(Resolution_2_det[sipm] * pow(**Energy, 2), 2));
            //         if (gRandom->Gaus(**Energy, sigma_resolution) < gRandom->Gaus(Threshold, Threshold_STD))
            //         {
            //             flagUntriggered = true;
            //             break;
            //         }
            //     }

            //     if (flagUntriggered)
            //         continue;

            //     H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Fill(energy_conv);
            // }
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
            double start = 200;
            double end = 6000;
            if (YEAR == 2025)
            {
                // TAKING ONLY THE LOW GAIN (threshold are low for low gain in 2025)
                // H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].first = (TH1D*)H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].second->Clone();

                // full range
                start = 200;
                end = WindowsBetaMap[NUCLEUS][peak_number];

                //scaling sim
                H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].second->GetXaxis()->SetRangeUser(start, end);
                H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetXaxis()->SetRangeUser(start, end);

                if (peak_number == 1 && NUCLEUS == "32Ar")
                {
                    int bin = H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetNbinsX();
                    TH1D *H = (TH1D*)H_Sim_Conv[NUCLEUS][14][current_detector]->Clone();
                    H->Rebin(H->GetNbinsX() / bin);

                    H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Scale(1 / H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Integral());
                    H->Scale(0.7*0.24/0.143 / H->Integral()); // RATIO BR GAMMA TO 2208 OVER BETA TO 2208 (x solid angle gamma for beta gamma pile up)

                    H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Add(H);
                    H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Scale(H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].second->Integral() / H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Integral());
                }   

                else
                {
                    H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Scale(H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].second->Integral() / H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Integral());
                    // H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Scale(H_SiPM_Calibrated[NUCLEUS][(int)IAS[NUCLEUS]][current_detector].second->Integral() / H_Sim_Conv[NUCLEUS][(int)IAS[NUCLEUS]][current_detector]->Integral());
                }
                

                // Chi2
                chi2.push_back(H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].second->Chi2Test(H_Sim_Conv[NUCLEUS][peak_number][current_detector], "CHI2/NDF"));
            }

            else if (YEAR == 2024)
            {
                // TAKING LOW AND HIGH with a limit
                // ranges
                double start_high = 200;
                double end_high = 2000;

                double start_low = 800;
                double end_low = WindowsBetaMap[NUCLEUS][peak_number];

                // Scaling Sim for High
                H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].first->GetXaxis()->SetRangeUser(start_high, end_high);
                H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetXaxis()->SetRangeUser(start_high, end_high);
                H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Scale(H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].first->Integral() / H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Integral());
                //chi2 
                chi2.push_back(H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].first->Chi2Test(H_Sim_Conv[NUCLEUS][peak_number][current_detector], "CHI2"));

                // Scaliong sim for Low
                H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].second->GetXaxis()->SetRangeUser(start_low, end_low);
                H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetXaxis()->SetRangeUser(start_low, end_low);
                H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Scale(H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].second->Integral() / H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Integral());   
                // chi2
                chi2.push_back(H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector].second->Chi2Test(H_Sim_Conv[NUCLEUS][peak_number][current_detector], "CHI2"));
            }

            if (FLAGSAVINGCHI2)
            {
                CHI2[NUCLEUS][peak_number][current_detector] = chi2[chi2.size() - 1];
            }
        }
    }

    
    double mean = accumulate(chi2.begin(), chi2.end(), 0.0) / chi2.size();
    if (FLAGSAVINGCHI2)
    {
        Info("Chi2: " + to_string(mean), 1);
        Info("Calibration: " + to_string(par[0]) + "   " + to_string(par[1]) + "   " + to_string(par[2]) + "   " + to_string(par[3]) + "   " + to_string(par[4] * 1e5) + "e-5   " + to_string(par[5]) + "   " + to_string(par[6]), 1);
    }
    return mean;
}

void CalibrationSiPM()
{
    VERBOSE = 0;
    Info("Calibration", 1);
    Minimizer *minimizer = Factory::CreateMinimizer("Minuit2", "Minuit");
    ROOT::Math::Functor functor(&Chi2TreeHist_conv, 7);
    
    vector<double> Calib;
    vector<double> Resolution_SQRT_det;
    vector<double> Resolution_SQR_det;
    // ## 2024 ## //

    //2024
    // vector<double> Par = {0, 0.64, 0, Resolution_SQRT_det[current_detector], 2.5e-10, 150, 18};
    // vector<double> Par = {0, 0.64, 0, 7.27, 1.7e-5, 150, 18};
    // vector<double> Par = {0, 0.74, 0, 1.92, 1.7e-5, 150, 18};
    // vector<double> Par = {0, 0.64, 0, 7.27, 1.7e-5, 150, 18};

    // ## 2025 ## //
    // Full Range on Low Gain
    if (YEAR == 2025)
    {
        Calib = {0, 1.309, 1.367, 1.359, 1.4099, 1.379, 1.463, 1.338, 1.360, 1.352};
        Resolution_SQRT_det = {0, 11.15, 8.99, 9.40, 5.39, 5.717, 0.01, 10.96, 7.40, 9.632};
        Resolution_SQR_det = {0, 0.178, 2.04, 4.71, 2.97, 1.63, 3.0, 5.047, 2.03, 1.248};

        // 100 - 2000 Range on High Gain
        // vector<double> Calib =              {0,     1.309,    1.367,      1.359,      1.4099,    1.379,    1.463,    1.338,   1.360,    1.352};
        // Resolution_SQRT_det =               {0,     11.15,    8.99,       9.40,       5.39,     5.717,      0.01,    10.96,    7.40,    9.632};
        // vector<double> Resolution_SQR_det = {0,     0.178,    2.04,       4.71,       2.97,     1.63,     3.0,     5.047,     2.03,     1.248};
        // vector<double> Par = {0, Calib[current_detector], 0, Resolution_SQRT_det[current_detector], Resolution_SQR_det[current_detector]*1e-5, 0, 18};
    }
    else if (YEAR == 2024)
    {
        Calib =                 {0, 1.6,    1.617,  1.601,  1.608,  1.620,   1.574,     1.558,  1.634,   1.599};
        Resolution_SQRT_det =   {0, 14.5,   13.71,  11.39,  6.547,  15.74,   9.85,      8.06,   10.96,   13.677};
        Resolution_SQR_det =    {0, 0.5,    0.348,  0.38,   2.17,   1.14,    0.38,      2.26,   1.92,    1.23};
    }
    else Error("YEAR not set correctly");
    

    vector<double> Par = {0, Calib[current_detector], 0, Resolution_SQRT_det[current_detector], Resolution_SQR_det[current_detector] * 1e-5, 0, 18};

    minimizer->SetFunction(functor);
    minimizer->SetFixedVariable(0, "Calibration_OffSet", Par[0]);
    minimizer->SetLimitedVariable(1, "Calibration", Par[1], 0.1, 0.5, 2.0);
    minimizer->SetFixedVariable(2, "Resolution_OffSet", Par[2]);
    minimizer->SetLimitedVariable(3, "Resolution_SQRT", Par[3], 0.5, 0.1, 16.);
    minimizer->SetLimitedVariable(4, "Resolution_2", Par[4], 1e-6, 1e-7, 1e-4);
    minimizer->SetFixedVariable(5, "Threshold", Par[5]);
    minimizer->SetFixedVariable(6, "Threshold_STD", Par[6]);
    minimizer->SetMaxFunctionCalls(10000000);
    minimizer->SetMaxIterations(100000000);
    // minimizer->SetTolerance(1e-5);
    minimizer->SetPrecision(1e-4);
    const double *bestPar = Par.data();
    
    PreviousPar = {Par[0], Par[1], Par[2], Par[3], Par[4], Par[5], Par[6]};
    // minimizer->Minimize();
    // bestPar = minimizer->X();

    FLAGSAVINGCHI2 = true;
    double chi2 = Chi2TreeHist_conv(bestPar);
    FLAGSAVINGCHI2 = false;
    // print results
    minimizer->PrintResults();

    CALIBRATED_File->cd();

    TDirectory *dir_cal = CALIBRATED_File->mkdir(("SiPM_" + to_string(current_detector)).c_str());
    dir_cal->cd();

    TF1 *F_Calibration = new TF1(("F_Calibration_SiPM" + to_string(current_detector)).c_str(), "[0]*x + [1]", 0, 10000);
    F_Calibration->SetParameters(bestPar[1], bestPar[0]);
    F_Calibration->Write();

    TF1 * F_Resolution = new TF1(("F_Resolution_SiPM" + to_string(current_detector)).c_str(), "sqrt(pow([0], 2) + pow([1]*sqrt(x), 2) + pow([2]*x*x, 2))", 0, 10000);
    F_Resolution->SetParameters(bestPar[2], bestPar[3], bestPar[4]);
    F_Resolution->Write();
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

        TCanvas *cExp_Sim_High[50];
        TCanvas *cExp_Sim_Low[50];
        TCanvas *cExp_Sim_SiPM_AllPeaks_High[10];
        TCanvas *cExp_Sim_SiPM_AllPeaks_Low[10];
        for (int peak = 0; peak <= 50; peak++)
        {
            if ((WindowssMap[NUCLEUS][peak][11].first == -1 || !WindowssMap[NUCLEUS][peak][11].first) && (NUCLEUS == "32Ar" || NUCLEUS == "33Ar"))
                continue;
            
            if (peak != 0 && (NUCLEUS == "207Bi" || NUCLEUS == "90Sr"))
                continue;
            
            cExp_Sim_High[peak] = new TCanvas(("cExp_Sim_HIgh_" + NUCLEUS + "_Peak_" + to_string(peak)).c_str(), ("cExp_Sim_HIgh_" + NUCLEUS + "_Peak_" + to_string(peak)).c_str(), 800, 800);
            cExp_Sim_High[peak]->Divide(3, 3);

            cExp_Sim_Low[peak] = new TCanvas(("cExp_Sim_Low_" + NUCLEUS + "_Peak_" + to_string(peak)).c_str(), ("cExp_Sim_Low_" + NUCLEUS + "_Peak_" + to_string(peak)).c_str(), 800, 800);
            cExp_Sim_Low[peak]->Divide(3, 3);
        }

        // loop on SiPM
        for (int det = 1; det <= 9; det++)
        {
            // Info("SiPM: " + to_string(det), 2);
            H[NUCLEUS][100 + det]->Write();
            H[NUCLEUS][110 + det]->Write();
            cExp_Sim_SiPM_AllPeaks_High[det] = new TCanvas(("cExp_Sim_High" + NUCLEUS + "_SiPM" + to_string(det)).c_str(), ("cExp_Sim_High" + NUCLEUS + "_SiPM" + to_string(det)).c_str(), 800, 800);
            cExp_Sim_SiPM_AllPeaks_High[det]->Divide(CanvasMap[NUCLEUS].first, CanvasMap[NUCLEUS].second);
            cExp_Sim_SiPM_AllPeaks_Low[det] = new TCanvas(("cExp_Sim_Low" + NUCLEUS + "_SiPM" + to_string(det)).c_str(), ("cExp_Sim_Low" + NUCLEUS + "_SiPM" + to_string(det)).c_str(), 800, 800);
            cExp_Sim_SiPM_AllPeaks_Low[det]->Divide(CanvasMap[NUCLEUS].first, CanvasMap[NUCLEUS].second);
            // loop on peak
            for (int peak = 0; peak <= 50; peak++)
            {

                if ((WindowssMap[NUCLEUS][peak][11].first == -1 || !WindowssMap[NUCLEUS][peak][11].first) && (NUCLEUS == "32Ar" || NUCLEUS == "33Ar"))
                    continue;

                if (peak != 0 && (NUCLEUS == "207Bi" || NUCLEUS == "90Sr"))
                    continue;
              

                // Info("Peak: " + to_string(peak), 3);
                if (dir_peaks[NUCLEUS][peak] != NULL)
                    dir_peaks[NUCLEUS][peak]->cd();

                // H_SiPM_Calibrated[NUCLEUS][peak][current_detector].first->GetXaxis()->SetRangeUser(200, WindowsBetaMap[NUCLEUS][peak]);
                // H_SiPM_Calibrated[NUCLEUS][peak][current_detector].second->GetXaxis()->SetRangeUser(200, WindowsBetaMap[NUCLEUS][peak]);
                // H_Sim_Conv[NUCLEUS][peak][current_detector]->GetXaxis()->SetRangeUser(200, WindowsBetaMap[NUCLEUS][peak]);
                // H_Sim_Conv[NUCLEUS][peak][current_detector]->Scale(H_SiPM_Calibrated[NUCLEUS][peak][current_detector].first->Integral() / H_Sim_Conv[NUCLEUS][peak][current_detector]->Integral());

                // ##### for each SiPMx
                TCanvas *cExp_Sim_det = new TCanvas(("cExp_Sim_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(det)).c_str(), ("cExp_Sim_" + NUCLEUS + "_Peak_" + to_string(peak) + "_SiPM" + to_string(det)).c_str(), 800, 800);
                cExp_Sim_det->Divide(2, 1);
                // ## - high
                cExp_Sim_det->cd(1);
                // Spectrum
                TPad *padSpectrum_high = new TPad("padSpectrum_high", "padSpectrum_high", 0, 0.3, 1, 1);
                padSpectrum_high->Draw();
                padSpectrum_high->cd();
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
                }
                // Residuals 
                cExp_Sim_det->cd(1);
                TPad *padRes_high = new TPad("padRes_high", "padRes_high", 0, 0, 1, 0.3);
                padRes_high->Draw();
                padRes_high->cd();
                TGraphErrors *gr_res_high = GetResiduals(H_SiPM_Calibrated[NUCLEUS][peak][det].first, H_Sim_Conv[NUCLEUS][peak][det]);
                gr_res_high->GetYaxis()->SetTitle("Residuals (sim-exp)/exp");
                gr_res_high->SetMarkerStyle(20);
                gr_res_high->Draw("AP");
                TLine *line_high = new TLine(gr_res_high->GetXaxis()->GetXmin(), 0, gr_res_high->GetXaxis()->GetXmax(), 0);
                line_high->SetLineColor(kRed);
                line_high->Draw("SAME");
                // ## - low
                cExp_Sim_det->cd(2);
                //Spectrum 
                TPad *padSpectrum_low = new TPad("padSpectrum_low", "padSpectrum_low", 0, 0.3, 1, 1);
                padSpectrum_low->Draw();
                padSpectrum_low->cd();
                H_SiPM_Calibrated[NUCLEUS][peak][det].second->Draw("HIST");
                if (H_Sim_Conv[NUCLEUS][peak][det]->GetEntries() > 0)
                {
                    H_Sim_Conv[NUCLEUS][peak][det]->SetLineColor(kRed);
                    H_Sim_Conv[NUCLEUS][peak][det]->Draw("HIST SAME");
                }
                // Residuals
                cExp_Sim_det->cd(2);
                TPad *padRes_low = new TPad("padRes_low", "padRes_low", 0, 0, 1, 0.3);
                padRes_low->Draw();
                padRes_low->cd();
                TGraphErrors *gr_res_low = GetResiduals(H_SiPM_Calibrated[NUCLEUS][peak][det].second, H_Sim_Conv[NUCLEUS][peak][det]);
                gr_res_low->GetYaxis()->SetTitle("Residuals (sim-exp)/exp");
                gr_res_low->SetMarkerStyle(20);
                gr_res_low->Draw("AP");
                TLine *line_low = new TLine(gr_res_low->GetXaxis()->GetXmin(), 0, gr_res_low->GetXaxis()->GetXmax(), 0);
                line_low->SetLineColor(kRed);
                line_low->Draw("SAME");
                cExp_Sim_det->Write();

                // # for each peak
                // high
                cExp_Sim_High[peak]->cd(det);
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
                    // latex->Draw("SAME");
                }
                //low
                cExp_Sim_Low[peak]->cd(det);
                H_SiPM_Calibrated[NUCLEUS][peak][det].second->Draw("HIST");
                if (H_Sim_Conv[NUCLEUS][peak][det]->GetEntries() > 0)
                {
                    H_Sim_Conv[NUCLEUS][peak][det]->SetLineColor(kRed);
                    H_Sim_Conv[NUCLEUS][peak][det]->Draw("HIST SAME");
                
                    // write chi2 in txt on the canvas
                    TLatex *latex = new TLatex();
                    latex->SetNDC();
                    latex->SetTextSize(0.04);
                    latex->DrawLatex(0.1, 0.92, ("#chi^{2}_{#nu}: " + to_string(CHI2[NUCLEUS][peak][det])).c_str());
                    // latex->Draw("SAME");
                }

                // # for each SiPM all peaks
                //high
                cExp_Sim_SiPM_AllPeaks_High[det]->cd(peak);
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
                    // latex->Draw("SAME");
                }
                //low
                cExp_Sim_SiPM_AllPeaks_Low[det]->cd(peak);
                H_SiPM_Calibrated[NUCLEUS][peak][det].second->Draw("HIST");
                if (H_Sim_Conv[NUCLEUS][peak][det]->GetEntries() > 0)
                {
                    H_Sim_Conv[NUCLEUS][peak][det]->SetLineColor(kRed);
                    H_Sim_Conv[NUCLEUS][peak][det]->Draw("HIST SAME");

                    // write chi2 in txt on the canvas
                    TLatex *latex = new TLatex();
                    latex->SetNDC();
                    latex->SetTextSize(0.04);
                    latex->DrawLatex(0.1, 0.92, ("#chi^{2}_{#nu}: " + to_string(CHI2[NUCLEUS][peak][det])).c_str());
                    // latex->Draw("SAME");
                }
            }
        }

        for (int peak = 0; peak <= 50; peak++)
        {
            if ((WindowssMap[NUCLEUS][peak][11].first == -1 || !WindowssMap[NUCLEUS][peak][11].first) && (NUCLEUS == "32Ar" || NUCLEUS == "33Ar"))
                continue;
            if (peak != 0 && (NUCLEUS == "207Bi" || NUCLEUS == "90Sr"))
                continue;
            if (dir_peaks[NUCLEUS][peak] != NULL)
                dir_peaks[NUCLEUS][peak]->cd();

            cExp_Sim_High[peak]->Write();
            cExp_Sim_Low[peak]->Write();
        }

        if (NUCLEUS == "32Ar" || NUCLEUS == "33Ar")
        {
            for (int det = 1; det <= 9; det++)
            {
                dir_SiPM[NUCLEUS][det]->cd();
                cExp_Sim_SiPM_AllPeaks_High[det]->Write();
                cExp_Sim_SiPM_AllPeaks_Low[det]->Write();
            }
        }
    }
    Info("Write Histograms done");
}


#endif