#ifndef MATCHER_HH
#define MATCHER_HH

#include "../../../lib/SignalDict/Signal.h"
#include "Detectors.hh"

default_random_engine generator;
int VERBOSE = 0;
bool Fitting_Resolution_FLAG = false;
/// FILE ///
string MERGED_filename;
string MERGED_basefilename;
TFile *CALIBRATED_File;
map<string, TFile *> SIMULATED_File;
map<string, TFile *> MERGED_File;
TFile *CRADLE_File;
string NUCLEUS;

map<string, TTree *[SIGNAL_MAX]> MERGED_Tree_Detectors;
TTreeReader *Reader;
TTreeReaderArray<Signal> *Silicon;
double Channel;

// CALIBRATION PEAKS
map<string, vector<double>> CalibrationPeaks;
vector<string> Nuclei;
bool Resolution_applied = false;

// HISTOGRAMS
map<string, TH1D *[SIGNAL_MAX]> H_Exp;
map<string, TH1D *[SIGNAL_MAX]> H_Exp_Channel;
map<string, TH1D *[SIGNAL_MAX]> H_Sim;
map<string, TH1D *[SIGNAL_MAX]> H_Sim_Conv;
map<string, TGraphErrors *[SIGNAL_MAX]> G_Calibration;
TMultiGraph *MG_Global_Calibration[SIGNAL_MAX];
TGraphErrors *G_Calibration_Mean[SIGNAL_MAX];
TGraphErrors *G_Calibration_Alpha[SIGNAL_MAX];
map<string, map<double, TF1 *[SIGNAL_MAX]>> F_Peak_Exp;
map<string, map<double, TF1 *[SIGNAL_MAX]>> F_Peak_Sim;
TGraphErrors *G_Resolution_Peak[SIGNAL_MAX];

TGraphErrors *G_Calibration_FitParameters[3] = {new TGraphErrors(), new TGraphErrors(), new TGraphErrors()};
TGraphErrors *G_Resolution_FitParameter = new TGraphErrors();
TGraphErrors *G_Electronic_Resolution_FitParameter = new TGraphErrors();


map<string, TH1D *[SIGNAL_MAX]> H_Coincidence;
map<string, TH1D *[SIGNAL_MAX]> H_AntiCoincidence;
map<string, TH1D *[SIGNAL_MAX]> H_single;
map<string, TH1D *> H_Exp_All;
map<string, TH1D *> H_Exp_All_Up;
map<string, TH1D *> H_Exp_All_Down;
map<string, TH1D *> H_Sim_All;
map<string, TH1D *> H_Sim_All_Up;
map<string, TH1D *> H_Sim_All_Down;

TGraphErrors *G_Mean[SIGNAL_MAX];
TGraphErrors *G_Chi2[SIGNAL_MAX];
// Function
TF1 *Calibration_Function[SIGNAL_MAX];
TFitResultPtr Calibration_Result[SIGNAL_MAX];
map<string, TF1 *[SIGNAL_MAX]>Background_function;
TF1 *Threshold_function[SIGNAL_MAX];
TF1 *Alpha1_function;
TF1 *Alpha2_function;
TF1 *Alpha3_function;
double CHI2 = 1000000;
bool plotting = false;
int current_detector;
int current_peak;
double alpha1;
double alpha1_error;
double alpha2;
double alpha2_error;

// ARRAY
pair<double, double> Alpha1_Up[6];
pair<double, double> Alpha2_Up[6];
pair<double, double> Alpha1_Down[6];
pair<double, double> Alpha2_Down[6];
map<double, pair<double, double>[SIGNAL_MAX]> ManualCalib;
pair<double, double> ManualCalibLinear[SIGNAL_MAX];
vector<double> ManualCalibFitted[SIGNAL_MAX];
vector<double> ManualCalibFitted_Alpha[SIGNAL_MAX];
pair<double, double> Detector_Resolution[SIGNAL_MAX];
double Detector_ElectronicResolution[SIGNAL_MAX];

map<string, pair<double, double>[SIGNAL_MAX]> PileUp;

TDirectory *dir_detector[SIGNAL_MAX];
TDirectory *dir_peak_detector[SIGNAL_MAX];
map<string, TDirectory *[SIGNAL_MAX]> dir_nuclei_detector;

void InitHistograms(int first, int last)
{
    for (int i = first; i < last; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            if (NUCLEUS == "32Ar")
            {
                dir_detector[i] = CALIBRATED_File->mkdir(detectorName[i].c_str());
                dir_peak_detector[i] = dir_detector[i]->mkdir("Peaks");
                MG_Global_Calibration[i] = new TMultiGraph();
            }

            H_Coincidence[NUCLEUS][i] = new TH1D(("H_Coincidence_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_Coincidence_" + NUCLEUS + "_" + detectorName[i]).c_str(), eSiliN_cal, 0, 10000);
            H_Coincidence[NUCLEUS][i]->GetXaxis()->SetTitle("Energy [keV]");
            H_Coincidence[NUCLEUS][i]->GetYaxis()->SetTitle("Counts / keV");
            H_Coincidence[NUCLEUS][i]->GetXaxis()->CenterTitle();
            H_Coincidence[NUCLEUS][i]->GetYaxis()->CenterTitle();

            H_AntiCoincidence[NUCLEUS][i] = new TH1D(("H_AntiCoincidence_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_AntiCoincidence_" + NUCLEUS + "_" + detectorName[i]).c_str(), eSiliN_cal, 0, 10000);
            H_AntiCoincidence[NUCLEUS][i]->GetXaxis()->SetTitle("Energy [keV]");
            H_AntiCoincidence[NUCLEUS][i]->GetYaxis()->SetTitle("Counts / keV");
            H_AntiCoincidence[NUCLEUS][i]->GetXaxis()->CenterTitle();
            H_AntiCoincidence[NUCLEUS][i]->GetYaxis()->CenterTitle();

            H_single[NUCLEUS][i] = new TH1D(("H_single_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_single_" + NUCLEUS + "_" + detectorName[i]).c_str(), eSiliN_cal, 0, 10000);
            H_single[NUCLEUS][i]->GetXaxis()->SetTitle("Energy [keV]");
            H_single[NUCLEUS][i]->GetYaxis()->SetTitle("Counts / keV");
            H_single[NUCLEUS][i]->GetXaxis()->CenterTitle();
            H_single[NUCLEUS][i]->GetYaxis()->CenterTitle();

            dir_nuclei_detector[NUCLEUS][i] = dir_detector[i]->mkdir(NUCLEUS.c_str());

            H_Exp[NUCLEUS][i] = new TH1D(("H_Exp_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_Exp_" + NUCLEUS + "_" + detectorName[i]).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
            H_Exp[NUCLEUS][i]->GetXaxis()->SetTitle("Energy [keV]");
            H_Exp[NUCLEUS][i]->GetYaxis()->SetTitle("Counts / keV");
            H_Exp[NUCLEUS][i]->GetXaxis()->CenterTitle();
            H_Exp[NUCLEUS][i]->GetYaxis()->CenterTitle();

            H_Exp_Channel[NUCLEUS][i] = new TH1D(("H_Exp_Channel_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_Exp_Channel_" + NUCLEUS + "_" + detectorName[i]).c_str(), eSiliN, 0, eSiliMax / 1000);
            H_Exp_Channel[NUCLEUS][i]->GetXaxis()->SetTitle("Channel");
            H_Exp_Channel[NUCLEUS][i]->GetYaxis()->SetTitle("Counts / keV");
            H_Exp_Channel[NUCLEUS][i]->GetXaxis()->CenterTitle();
            H_Exp_Channel[NUCLEUS][i]->GetYaxis()->CenterTitle();

            G_Calibration[NUCLEUS][i] = new TGraphErrors();
            G_Calibration_Alpha[i] = new TGraphErrors();            
        }
    }

    H_Exp_All[NUCLEUS] = new TH1D(("H_Exp_All_" + NUCLEUS).c_str(), ("H_Exp_All_" + NUCLEUS).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
    H_Exp_All[NUCLEUS]->GetXaxis()->SetTitle("Energy [keV]");
    H_Exp_All[NUCLEUS]->GetYaxis()->SetTitle("Counts / keV");
    H_Exp_All[NUCLEUS]->GetXaxis()->CenterTitle();
    H_Exp_All[NUCLEUS]->GetYaxis()->CenterTitle();

    H_Exp_All_Up[NUCLEUS] = new TH1D(("H_Exp_All_Up_" + NUCLEUS).c_str(), ("H_Exp_All_Up_" + NUCLEUS).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
    H_Exp_All_Up[NUCLEUS]->GetXaxis()->SetTitle("Energy [keV]");
    H_Exp_All_Up[NUCLEUS]->GetYaxis()->SetTitle("Counts / keV");
    H_Exp_All_Up[NUCLEUS]->GetXaxis()->CenterTitle();
    H_Exp_All_Up[NUCLEUS]->GetYaxis()->CenterTitle();

    H_Exp_All_Down[NUCLEUS] = new TH1D(("H_Exp_All_Down_" + NUCLEUS).c_str(), ("H_Exp_All_Down_" + NUCLEUS).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
    H_Exp_All_Down[NUCLEUS]->GetXaxis()->SetTitle("Energy [keV]");
    H_Exp_All_Down[NUCLEUS]->GetYaxis()->SetTitle("Counts / keV");
    H_Exp_All_Down[NUCLEUS]->GetXaxis()->CenterTitle(); 
    H_Exp_All_Down[NUCLEUS]->GetYaxis()->CenterTitle();

    H_Sim_All[NUCLEUS] = new TH1D(("H_Sim_All_" + NUCLEUS).c_str(), ("H_Sim_All_" + NUCLEUS).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
    H_Sim_All[NUCLEUS]->GetXaxis()->SetTitle("Energy [keV]");
    H_Sim_All[NUCLEUS]->GetYaxis()->SetTitle("Counts / keV");
    H_Sim_All[NUCLEUS]->GetXaxis()->CenterTitle();
    H_Sim_All[NUCLEUS]->GetYaxis()->CenterTitle();

    H_Sim_All_Up[NUCLEUS] = new TH1D(("H_Sim_All_Up_" + NUCLEUS).c_str(), ("H_Sim_All_Up_" + NUCLEUS).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
    H_Sim_All_Up[NUCLEUS]->GetXaxis()->SetTitle("Energy [keV]");
    H_Sim_All_Up[NUCLEUS]->GetYaxis()->SetTitle("Counts / keV");
    H_Sim_All_Up[NUCLEUS]->GetXaxis()->CenterTitle();
    H_Sim_All_Up[NUCLEUS]->GetYaxis()->CenterTitle();

}

void InitAlphaPeaks()
{
    /// read txt file
    ifstream file("Config_Files/18N_Up.txt");
    string line;
    double energy1;
    double energy2;
    double energy3;
    double energy4;
    string name;
    int counter = 0;
    while (getline(file, line))
    {
        counter++;

        stringstream ss(line);
        ss >> name >> energy1 >> energy2 >> energy3 >> energy4;
        Alpha1_Up[counter] = make_pair(energy1, energy2);
        Alpha2_Up[counter] = make_pair(energy3, energy4);
    }

    file.close();

    ifstream file2("Config_Files/18N_Down.txt");
    counter = 0;
    while (getline(file2, line))
    {
        counter++;

        stringstream ss(line);
        ss >> name >> energy1 >> energy2 >> energy3 >> energy4;
        Alpha1_Down[counter] = make_pair(energy1, energy2);
        Alpha2_Down[counter] = make_pair(energy3, energy4);
    }
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
            ss >> detname >> channel6 >> channel14 >> channel29 >> energy6 >> energy14 >> energy29;
        if (YEAR == 2025 || YEAR == 2021)
            ss >> detname >> channel6 >> channel14 >> energy6 >> energy14;
        for (int i = 0; i < detectorNum; i++)
        {
            if (IsDetectorSiliStrip(i))
            {
                if (detname == detectorName[i])
                {
                    ManualCalib[6][i] = make_pair(channel6, energy6);
                    ManualCalib[14][i] = make_pair(channel14, energy14);
                    if (YEAR == 2024)
                        ManualCalib[29][i] = make_pair(channel29, energy29);
                    break;
                }
            }
        }
    }

    file.close();
}

void InitElectronicResolution()
{
    ifstream file(("Config_Files/" + to_string(YEAR) + "/Electronic_Resolution_" + to_string(YEAR) + ".txt").c_str());
    if (!file.is_open())
    {
        Error("Impossible to open Electronic_Resolution_" + to_string(YEAR) + ".txt");
    }

    string line;
    int code;
    double resolution;
    double err_resolution;
    int counter = 0;
    while (getline(file, line))
    {
        counter++;
        stringstream ss(line);
        ss >> code >> resolution >> err_resolution;

        Detector_ElectronicResolution[code] = resolution / 1000.;
    }

    file.close();
}

void InitPileUp()
{
    /// setting all to 0 before replacing if given in the file
    for (string nuclei : Nuclei)
    {
        for (int i = 0; i < detectorNum; i++)
        {
            if (IsDetectorSiliStrip(i))
            {
                PileUp[nuclei][i] = make_pair(0, 0);
            }
        }
    }


    ifstream file("Config_Files/" + to_string(YEAR) + "/PileUp_" + to_string(YEAR) + ".txt");
    if (!file.is_open())
        Warning("Impossible to open PileUp_" + to_string(YEAR) + ".txt");

    string line;
    string nuclei;
    string detname;
    double sigma;
    double p;
    string dir;
    int strip;
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
        ss >> detname >> sigma >> p;

        // last character of the string
        strip = atoi(detname.substr(detname.length() - 1).c_str());
        dir = detname.substr(0, detname.length() - 2);

        for (int i : Dir2Det(dir, strip))
        {
            PileUp[nuclei][i] = make_pair(sigma, p);
            // cout << "Nuclei : " << nuclei << " Detector : " << detectorName[i] << " Sigma : " << sigma << " PileUp : " << p << endl;
        }
    }
    file.close();
}

void FillingSimHitograms()
{
    Info("Filling Simulated Histograms");
    /////// HAVE TO BE CHANGE TO TAKE INTO ACCOUNT THE INTERSTRIPS EVENTS //////

    IAS["32Ar"] = 14;
    IAS["32Ar_thick"] = 14;
    IAS["33Ar"] = 21;
    IAS["33Ar_thick"] = 21;
    IAS["18N"] = 1;

    for (auto &pair : SIMULATED_File)
    {
        Info("Nucleus : " + pair.first, 1);
        for (int i = 0; i < detectorNum; i++)
        {
            
            if (IsDetectorSiliStrip(i))
            {
                if (VERBOSE  == 1) Info("Detector : " + detectorName[i], 2);

                string h_name = "Silicon_Detector_Energy_Deposit_" + detectorName[i] + "_All";
                // if (pair.first != "18N")
                // {
                //     h_name = "p/Silicon_Detector/D." + to_string(GetDetector(i)) + "/Silicon_Detector_Energy_Deposit_" + detectorName[i] + "_p";
                //     
                // }
                H_Sim[pair.first][i] = (TH1D *)pair.second->Get((h_name).c_str());
                H_Sim[pair.first][i]->SetName(("H_Sim_" + detectorName[i]).c_str());
                H_Sim[pair.first][i]->SetTitle(("H_Sim_" + detectorName[i]).c_str());
                H_Sim[pair.first][i]->GetXaxis()->SetTitle("Energy [keV]");
                H_Sim[pair.first][i]->GetYaxis()->SetTitle("Counts");
                H_Sim[pair.first][i]->GetXaxis()->CenterTitle();
                H_Sim[pair.first][i]->GetYaxis()->CenterTitle();
                for (int j = 0; j < detectorNum; j++)
                {
                    string h_name = "Silicon_Detector_Energy_Deposit_" + detectorName[i] + "_All";

                    // if (pair.first != "18N")
                    // {
                    //     h_name = "p/Silicon_Detector/D." + to_string(GetDetector(i)) + "/Silicon_Detector_Energy_Deposit_" + detectorName[i] + "_p";
                    //     h_name = "Silicon_Detector_Energy_Deposit_" + detectorName[i] + "_All";
                    // }
                    if ((GetDetector(i) != GetDetector(j)) && (GetDetectorChannel(i) == GetDetectorChannel(j)) && IsDetectorSiliStrip(j))
                    {
                        if (GetDetector(i) <= 4 && GetDetector(j) <= 4)
                            H_Sim[pair.first][i]->Add((TH1D *)pair.second->Get((h_name).c_str()));
                        else if (GetDetector(i) >= 5 && GetDetector(j) >= 5)
                            H_Sim[pair.first][i]->Add((TH1D *)pair.second->Get((h_name).c_str()));
                    }
                }
            }
        }
    }
}

void WriteCalibInFile()
{
    // create and write calibration in txt file
    ofstream file(("Config_Files/Calibration_" + to_string(YEAR) + ".txt").c_str());

    for (int i = 0; i < detectorNum; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            file << i << "\t" << ManualCalibFitted[i][0] << "\t" << ManualCalibFitted[i][1] << "\t" << ManualCalibFitted[i][2] << endl;
        }
    }

    file.close();
}

void WriteResolution()
{
    // create and write resolution in txt file
    ofstream file(("Config_Files/" + to_string(YEAR) + "/Resolution_" + to_string(YEAR) + ".txt").c_str());

    if (!file.is_open())
    {
        Error("Impossible to open Resolution_" + to_string(YEAR) + ".txt");
    }

    for (int i = 0; i < detectorNum; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            file << i << "\t" << Detector_Resolution[i].first << "\t" << Detector_Resolution[i].second << endl;
        }
    }

    file.close();
}

void InitResolution()
{
    if (Fitting_Resolution_FLAG)
    {
        for (int det = 0; det < detectorNum; det++)
        {
            if (IsDetectorSiliStrip(det))
            {
                Detector_Resolution[det] = make_pair(0, 0);
            }
        }
        return;
    } 
    else
    {
        Warning("Resolution will not be fitted");
    }

    ifstream file(("Config_Files/" + to_string(YEAR) + "/Resolution_" + to_string(YEAR) + ".txt").c_str());
    if (!file.is_open())
    {
        Error("Impossible to open Resolution_" + to_string(YEAR) + ".txt");
    }

    string line;
    int code;
    double resolution, err;
    int counter = 0;
    while (getline(file, line))
    {
        counter++;
        stringstream ss(line);
        ss >> code >> resolution >> err;

        Detector_Resolution[code] = make_pair(resolution, err);
    }

    file.close();
    Success("Resolution loaded");
}

TF1 *InvertFunction(TF1 *f)
{
    if (f->GetNpar() == 2)
    {
        // bijective of linear fucntion
        double a = f->GetParameter(1);
        double b = f->GetParameter(0);

        TF1 *f = new TF1("f", Form("(x - %f)/%f", b, a), 0, 10000);

        // same WHTHOUT Form but takng values of a and b

        return f;
    }

    if (f->GetNpar() == 3)
    {
        // bijective of quadratic function
        double a = f->GetParameter(2);
        double b = f->GetParameter(1);
        double c = f->GetParameter(0);

        TF1 *f_plus = new TF1("f_plus", Form("(-%f + sqrt(%f^2 - 4*%f*(%f-x)))/(2*%f)", b, b, a, c, a), 0, 10000);
        TF1 *f_minus = new TF1("f_minus", Form("(-%f - sqrt(%f^2 - 4*%f*(%f-x)))/(2*%f)", b, b, a, c, a), 0, 10000);

        return f_plus;
    }
    else
    {
        Error("Function not implemented : Degree = " + to_string(f->GetNpar() - 1));
        return nullptr;
    }
}

void MakeResidualPlot(string NUCLEUS)
{
    TCanvas *c1 = new TCanvas((detectorName[current_detector] + "_" + NUCLEUS + "_IAS").c_str(), (detectorName[current_detector] + "_" + NUCLEUS + "_IAS").c_str(), 1920, 1080);
    TPad *p1 = new TPad("Spectrum", "Spectrum", 0, 0.3, 1, 1);
    p1->Draw();
    p1->cd();
    p1->SetLogy();
    H_Exp[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][IAS[NUCLEUS]][current_detector].first, WindowsMap[NUCLEUS][IAS[NUCLEUS]][current_detector].second);
    H_Sim_Conv[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][IAS[NUCLEUS]][current_detector].first, WindowsMap[NUCLEUS][IAS[NUCLEUS]][current_detector].second);

    H_Exp[NUCLEUS][current_detector]->Draw("HIST");
    H_Sim_Conv[NUCLEUS][current_detector]->SetLineColor(kRed);
    H_Sim_Conv[NUCLEUS][current_detector]->Draw("HIST SAME");

    c1->cd();
    TPad *p2 = new TPad("Residuals", "Residuals", 0, 0, 1, 0.3);
    p2->Draw();
    p2->cd();
    TH1D *H_Residuals = (TH1D *)H_Exp[NUCLEUS][current_detector]->Clone();
    H_Residuals->Add(H_Sim_Conv[NUCLEUS][current_detector], -1);
    H_Residuals->Divide(H_Sim_Conv[NUCLEUS][current_detector]);

    TGraphErrors *G_Residuals = new TGraphErrors();
    int Entries = H_Exp[NUCLEUS][current_detector]->GetEntries();
    for (int bin = 0; bin < Entries; bin++)
    {
        if (H_Sim_Conv[NUCLEUS][current_detector]->GetBinContent(bin) == 0 || H_Exp[NUCLEUS][current_detector]->GetBinContent(bin) == 0)
            continue;
        G_Residuals->AddPoint(H_Sim_Conv[NUCLEUS][current_detector]->GetBinCenter(bin), H_Residuals->GetBinContent(bin));
        G_Residuals->SetPointError(G_Residuals->GetN() - 1, 0, H_Residuals->GetBinError(bin));
    }
    G_Residuals->SetMarkerStyle(20);
    G_Residuals->SetMarkerSize(0.5);
    G_Residuals->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][IAS[NUCLEUS]][current_detector].first, WindowsMap[NUCLEUS][IAS[NUCLEUS]][current_detector].second);
    G_Residuals->GetYaxis()->SetRangeUser(-1, 1);
    G_Residuals->Draw("AP");

    TLine *line = new TLine(WindowsMap[NUCLEUS][IAS[NUCLEUS]][current_detector].first, 0, WindowsMap[NUCLEUS][IAS[NUCLEUS]][current_detector].second, 0);
    line->SetLineColor(kRed);
    line->Draw("SAME");

    dir_detector[current_detector]->cd();
    c1->Write();
}

TH1D *RemoveBKG(TH1D *H, bool finished, string NUCLEUS)
{
    double MINIMUM = 0;
    double MAXIMUM = 0;
    int N = 0;
    double coef = 1;
    double offset = 0;

    MINIMUM = eSiliMin_cal;
    MAXIMUM = eSiliMax_cal;
    N = eSiliN_cal;
    


    double Threshold = 600;
    TH1D *Selection = (TH1D *)H->Clone("BKG");
    Selection->GetXaxis()->SetRangeUser(MINIMUM, MAXIMUM);
    //////////////////////////////////// with alpha contamination //////////////////////////////////
    if (YEAR == 2021 || YEAR == 2024)
    {
        if (YEAR == 2021)
            Threshold = 300;

        TF1 *BKG_function = new TF1("BKG_function", "( gaus(0) + gaus(3) + [6]*exp([7] * x) ) * 0.5*(1+TMath::Erf((x-[8])/[9]))", MINIMUM, MAXIMUM);
        BKG_function->SetNpx(N);
        double Integral_A1 = 40;
        double Sigma_A1 = 25;
        double Mean_A1_low;
        double Mean_A1_high;
        double Integral_A2 = 40;
        double Mean_A2_low;
        double Mean_A2_high;
        double Sigma_A2 = 25;

        if (GetDetector(current_detector) <= 4)
        {
            Mean_A1_low = (WindowsMap["18N"][1][current_detector].first - offset) / coef;
            Mean_A1_high = (WindowsMap["18N"][1][current_detector].second - offset) / coef;
            Mean_A2_low = (WindowsMap["18N"][2][current_detector].first - offset) / coef;
            Mean_A2_high = (WindowsMap["18N"][2][current_detector].second - offset) / coef;
        }
        else
        {
            Mean_A1_low = (WindowsMap["18N"][1][current_detector].first - offset) / coef;
            Mean_A1_high = (WindowsMap["18N"][1][current_detector].second - offset) / coef;
            Mean_A2_low = (WindowsMap["18N"][2][current_detector].first - offset) / coef;
            Mean_A2_high = (WindowsMap["18N"][2][current_detector].second - offset) / coef;
        }

        // BKG_function->SetParameters(1000, (Mean_A1_low+Mean_A1_high)/2, 25/coef, 500, (Mean_A2_low+Mean_A2_high)/2, 25/coef, 1000, -3e-3*coef, 500., 25.);
        BKG_function->SetParLimits(0, 5, 30000); // ALPHA1 - AMPLITUDE
        BKG_function->SetParameter(0, 1000);
        BKG_function->SetParLimits(1, Mean_A1_low, Mean_A1_high); // ALPHA1 - MEAN
        BKG_function->SetParameter(1, (Mean_A1_low + Mean_A1_high) / 2);
        BKG_function->SetParLimits(2, 5 / coef, 100 / coef); // ALPHA1 - SIGMA
        BKG_function->SetParameter(2, 25 / coef);
        BKG_function->SetParLimits(3, 5, 30000); // ALPHA2 - AMPLITUDE
        BKG_function->SetParameter(3, 500);
        BKG_function->SetParLimits(4, Mean_A2_low, Mean_A2_high); // ALPHA2 - MEAN
        BKG_function->SetParameter(4, (Mean_A2_low + Mean_A2_high) / 2);
        BKG_function->SetParLimits(5, 5 / coef, 75 / coef); // ALPHA2 - SIGMA
        BKG_function->SetParameter(5, 25 / coef);
        BKG_function->SetParLimits(6, 5, 15000); // EXP - AMPLITUDE
        BKG_function->SetParameter(6, 1000);
        BKG_function->SetParLimits(7, -1e-2 * coef, 0); // EXP - SLOPE
        BKG_function->SetParameter(7, -3e-3 * coef);
        BKG_function->SetParLimits(8, Threshold - 200, Threshold + 200); // THRESHOLD - MEAN
        BKG_function->SetParameter(8, Threshold);
        BKG_function->SetParLimits(9, 5 / coef, 100 / coef); // THRESHOLD - SIGMA
        BKG_function->SetParameter(9, 25);
        // BKG_function->SetParLimits(10, 5, 30000); // ALPHA3 - AMPLITUDE
        // BKG_function->SetParameter(10, 2000);
        // BKG_function->SetParLimits(11, 1800, 2500); // ALPHA3 - MEAN
        // BKG_function->SetParameter(11, 1900);
        // BKG_function->SetParLimits(12, 50, 500); // ALPHA3 - SIGMA
        // BKG_function->SetParameter(12, 250);

        Selection->Fit("BKG_function", "QRN", "", (400 - offset) / coef, (1900 - offset) / coef);

        alpha1 = BKG_function->GetParameter(1) * coef + offset;
        alpha2 = BKG_function->GetParameter(4) * coef + offset;

        Alpha1_function = new TF1("Alpha1_function", "gaus(0)", MINIMUM, MAXIMUM);
        Alpha1_function->SetParameters(BKG_function->GetParameter(0), BKG_function->GetParameter(1), BKG_function->GetParameter(2));

        Alpha2_function = new TF1("Alpha2_function", "gaus(0)", MINIMUM, MAXIMUM);
        Alpha2_function->SetParameters(BKG_function->GetParameter(3), BKG_function->GetParameter(4), BKG_function->GetParameter(5));

        // Alpha3_function = new TF1("Alpha3_function", "gaus(0)", MINIMUM, MAXIMUM);
        // Alpha3_function->SetParameters(BKG_function->GetParameter(10), BKG_function->GetParameter(11), BKG_function->GetParameter(12));

        Background_function[NUCLEUS][current_detector] = new TF1(("BKG_function_sub_" + detectorName[current_detector]).c_str(), "  [0]*exp([1] * x) ", MINIMUM, MAXIMUM);
        Background_function[NUCLEUS][current_detector]->SetNpx(N);
        Background_function[NUCLEUS][current_detector]->SetParameters(BKG_function->GetParameter(6), BKG_function->GetParameter(7));
        Threshold_function[current_detector] = new TF1(("Threshold_function_" + detectorName[current_detector]).c_str(), "0.5*(1+TMath::Erf((x-[0])/[1]))", MINIMUM, MAXIMUM);
        Threshold_function[current_detector]->SetNpx(N);
        Threshold_function[current_detector]->SetParameters(BKG_function->GetParameter(8), BKG_function->GetParameter(9));

        /// PLOT//////////////////////////////
        if (finished)
        {
            TCanvas *c = new TCanvas((detectorName[current_detector] + "BKG").c_str(), (detectorName[current_detector] + "BKG").c_str(), 1920, 1080);
            c->Divide(1, 2);
            c->cd(1);
            Selection->GetXaxis()->SetRangeUser((400 - offset) / coef, (1900 - offset) / coef);
            Selection->SetLineColor(kBlack);
            Selection->Draw("HIST");

            //// SUBSTRACTING THE BACKGROUND ////
            BKG_function->SetLineColor(kGreen);
            BKG_function->Draw("SAME");
            TF1 *BKG_function_sub = new TF1("BKG_function_sub", "[0]*exp([1] * x) * 0.5*(1+TMath::Erf((x-[2])/[3]))", MINIMUM, MAXIMUM);
            BKG_function_sub->SetNpx(N);
            BKG_function_sub->SetParameters(BKG_function->GetParameter(6), BKG_function->GetParameter(7), BKG_function->GetParameter(8), BKG_function->GetParameter(9));
            BKG_function_sub->SetLineColor(kRed);
            BKG_function_sub->Draw("SAME");

            c->cd(2);
            TH1D *H_SUB = (TH1D *)BKG_function_sub->GetHistogram()->Clone();
            TH1D *h_sub = (TH1D *)Selection->Clone(("h_sub" + detectorName[current_detector]).c_str());
            h_sub->Add(H_SUB, -1);
            h_sub->SetLineColor(kBlack);
            h_sub->GetYaxis()->SetRangeUser(1, -1111);
            h_sub->GetXaxis()->SetRangeUser((400 - offset) / coef, (1900 - offset) / coef);
            h_sub->Draw("HIST");
            dir_detector[current_detector]->cd();
            c->Write();
            delete BKG_function;
            delete H_SUB;
            delete BKG_function_sub;
            delete Selection;
            return h_sub;
        }
        else
        {
            TF1 *BKG_function_sub = new TF1("BKG_function_sub", "[0]*exp([1] * x) * 0.5*(1+TMath::Erf((x-[2])/[3]))", MINIMUM, MAXIMUM);
            BKG_function_sub->SetParameters(BKG_function->GetParameter(6), BKG_function->GetParameter(7), BKG_function->GetParameter(8), BKG_function->GetParameter(9));
            BKG_function_sub->SetNpx(N);
            TH1D *H_SUB = (TH1D *)BKG_function_sub->GetHistogram()->Clone();
            TH1D *h_sub = (TH1D *)Selection->Clone(("h_sub" + detectorName[current_detector]).c_str());
            h_sub->Add(H_SUB, -1);
            delete BKG_function;
            delete H_SUB;
            delete BKG_function_sub;
            delete Selection;
            return h_sub;
        }
    }
    /////////////////////////////////////
    else
    {
        double max = 1500;
        double min = 400;
        ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);
        TF1 *BKG_function = new TF1("BKG_function", "([0]*exp([1] * x) + [7]*exp([8] * x) ) * 0.5*(1+TMath::Erf((x-[2])/[3])) + gaus(4) + gaus(9)", MINIMUM, MAXIMUM);
        BKG_function->SetNpx(N);

        BKG_function->SetParLimits(0, 0, 3000); // EXP1 - AMPLITUDE
        BKG_function->SetParameter(0, 400);
        BKG_function->SetParLimits(1, -0.005, 0); // EXP1 - SLOPE
        BKG_function->SetParameter(1, -0.003);

        BKG_function->SetParLimits(2, 20, 300); // // THRESHOLD - MEAN
        BKG_function->SetParameter(2, 100);
        BKG_function->SetParLimits(3, 0, 100); // THRESHOLD - SIGMA
        BKG_function->SetParameter(3, 20);

        BKG_function->SetParLimits(4, 0, 1000); // Proton1 - AMPLITUDE
        BKG_function->SetParameter(4, 200);
        BKG_function->SetParLimits(5, 500, 650); // Proton1 - MEAN
        BKG_function->SetParameter(5, 600);
        BKG_function->SetParLimits(6, 1, 50); // Proton1 - SIGMA
        BKG_function->SetParameter(6, 5);

        BKG_function->SetParLimits(7, 5, 30000); // EXP2 - AMPLITUDE
        BKG_function->FixParameter(7, 0); //2500
        BKG_function->SetParLimits(8, -1, -0.006); // EXP2 - SLOPE
        BKG_function->SetParameter(8, -0.015);

        BKG_function->SetParLimits(9, 5, 30000); // Proton2 - AMPLITUDE
        BKG_function->SetParameter(9, 1000);
        // BKG_function->FixParameter(9, 0);
        BKG_function->SetParLimits(10, 1100, 1300); // Proton2 - MEAN
        BKG_function->SetParameter(10, 1200);
        // BKG_function->FixParameter(10, 0);
        BKG_function->SetParLimits(11, 1, 100); // Proton2 - SIGMA
        BKG_function->SetParameter(11, 5);
        // BKG_function->FixParameter(11, 0);


        Selection->Fit("BKG_function", "RQN", "", min, max);

        Background_function[NUCLEUS][current_detector] = new TF1(("BKG_function_sub_" + detectorName[current_detector]).c_str(), "  [0]*exp([1] * x) + [2]*exp([3] * x)", MINIMUM, MAXIMUM);
        Background_function[NUCLEUS][current_detector]->SetNpx(N);
        Background_function[NUCLEUS][current_detector]->SetParameters(BKG_function->GetParameter(0), BKG_function->GetParameter(1), BKG_function->GetParameter(7), BKG_function->GetParameter(8));
        
        Threshold_function[current_detector] = new TF1(("Threshold_function_" + detectorName[current_detector]).c_str(), "0.5*(1+TMath::Erf((x-[0])/[1]))", MINIMUM, MAXIMUM);
        Threshold_function[current_detector]->SetNpx(N);
        Threshold_function[current_detector]->SetParameters(BKG_function->GetParameter(2), BKG_function->GetParameter(3));

        //PLOT
        if (finished)
        {
            TCanvas *c = new TCanvas((detectorName[current_detector] + "BKG").c_str(), (detectorName[current_detector] + "BKG").c_str(), 1920, 1080);
            c->Divide(1, 2);
            c->cd(1);
            Selection->GetXaxis()->SetRangeUser(50, max);
            Selection->SetLineColor(kBlack);
            Selection->Draw("HIST");

            BKG_function->SetLineColor(kGreen);
            BKG_function->Draw("SAME");

            c->cd(2);
            TH1D *H_SUB = (TH1D *)Background_function[NUCLEUS][current_detector]->GetHistogram()->Clone();
            TH1D *h_sub = (TH1D *)Selection->Clone(("h_sub" + detectorName[current_detector]).c_str());
            h_sub->Add(H_SUB, -1);
            h_sub->SetLineColor(kBlack);
            h_sub->GetYaxis()->SetRangeUser(1, -1111);
            h_sub->GetXaxis()->SetRangeUser(50, max);
            h_sub->Draw("HIST");
            dir_detector[current_detector]->cd();
            c->Write();
            delete BKG_function;
            delete H_SUB;
            // delete BKG_function_sub;
            delete Selection;
            return h_sub;
        }
        else
        {
            TH1D *H_SUB = (TH1D *)Background_function[NUCLEUS][current_detector]->GetHistogram()->Clone();
            TH1D *h_sub = (TH1D *)Selection->Clone(("h_sub" + detectorName[current_detector]).c_str());
            h_sub->Add(H_SUB, -1);
            delete BKG_function;
            delete H_SUB;
            delete Selection;
            return h_sub;
        }

    }   
}

void ApplyCalibration(int verbose = 0, double offset_correction = 0, string nuc = "", int indent = 1)
{
    if (verbose == 1) Info("Applying Calibration", indent);

    // LOOP ON NUCLEUS //
    if (nuc == "")
    {
        for (string NUCLEUS : Nuclei)
        {
            if (verbose == 1) Info("NUCLEUS : " + NUCLEUS, indent+1);

            Reader = new TTreeReader(MERGED_Tree_Detectors[NUCLEUS][current_detector]);
            Reader->Restart();
            H_Exp[NUCLEUS][current_detector]->Reset();
            TTreeReaderValue<double> ChannelDet(*Reader, "Channel");
            while (Reader->Next())
            {
                double value = *ChannelDet / 1000;
                H_Exp[NUCLEUS][current_detector]->Fill(Calibration_Function[current_detector]->Eval(value));
            }
        }
    }
    else
    {
        NUCLEUS = nuc;
        if (verbose == 1) Info("NUCLEUS : " + NUCLEUS, indent+1);

        Reader = new TTreeReader(MERGED_Tree_Detectors[NUCLEUS][current_detector]);
        Reader->Restart();
        H_Exp[NUCLEUS][current_detector]->Reset();
        TTreeReaderValue<double> ChannelDet(*Reader, "Channel");
        while (Reader->Next())
        {
            double value = *ChannelDet / 1000;
            H_Exp[NUCLEUS][current_detector]->Fill(Calibration_Function[current_detector]->Eval(value) + offset_correction);
        }
    }
}

TH1D* ApplyResolution(int Verbose, int current_detector, string nucleus, double resolution, int indent=2, bool Hist = false)
{
    // ### CONVOLUTION OF THE SIMULATED SPECTRUM ### //
    if (Verbose == 1)
        Info("CONVOLUTION OF THE SIMULATED SPECTRUM", indent);

    if (Hist)
    {
        TH1D* H_Conv = (TH1D*)SIMULATED_File[nucleus]->Get(("Silicon_Detector_Energy_Deposit_"+detectorName[current_detector]+"_All").c_str())->Clone(("H_Sim_Conv_" + detectorName[current_detector]).c_str());
        H_Conv->Reset();
        for (int det = 1; det <= 4; det++)
        {            
            int detector = det;
            detector += current_detector > 50 ? 4 : 0;
            int detector_number = detector * 10 + GetDetectorChannel(current_detector);

            if (Verbose == 1)
                Info("Detector : " + detectorName[detector_number], indent + 1);

            // loading hist
            TH1D* H = (TH1D*)SIMULATED_File[nucleus]->Get(("Silicon_Detector_Energy_Deposit_"+detectorName[detector_number]+"_All").c_str());

            // init pileup
            double PileUp_Probability = PileUp[NUCLEUS][current_detector].second;
            double PileUp_Sigma = PileUp[NUCLEUS][current_detector].first;
            // cout << H->GetEntries() << endl;
            for (int i = 0; i < H->GetEntries(); i++)
            {   
                double Simulated_Energy = H->GetRandom();
                // ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "Progress: ");
                normal_distribution<double> distribution(Simulated_Energy, resolution);
                double pileup = 0;
                if (gRandom->Uniform() < PileUp_Probability)
                {
                    pileup = abs(gRandom->Gaus(0, PileUp_Sigma));
                }
                H_Conv->Fill(distribution(generator) + pileup);
            }
        }
        return H_Conv;
    }
    
    for (string NUCLEUS : Nuclei)
    {
        if (Verbose == 1)
            Info("NUCLEUS : " + NUCLEUS, indent + 1);

        if (nucleus != NUCLEUS && nucleus != "")
            continue;

        // Resetting histograms
        H_Sim_Conv[NUCLEUS][current_detector] = (TH1D *)H_Sim[NUCLEUS][current_detector]->Clone();
        H_Sim_Conv[NUCLEUS][current_detector]->Reset();

        // Tree Method
        // looping if cylindrical symetry
        for (int det = 1; det <= 4; det++)
        {
            int detector = det;
            // if (current_detector > 50) detector += 4;
            detector += current_detector > 50 ? 4 : 0;
            int detector_number = detector * 10 + GetDetectorChannel(current_detector);

            // loading tree
            TTree *SIMULATED_Tree_Detectors = (TTree *)SIMULATED_File[NUCLEUS]->Get(("Tree_" + detectorName[detector_number]).c_str());
            Reader = new TTreeReader(SIMULATED_Tree_Detectors);
            TTreeReaderValue<double> Simulated_Energy(*Reader, "Energy");
            Reader->Restart();
            int Entries = SIMULATED_Tree_Detectors->GetEntries();
            clock_t start = clock(), Current;

            // init pileup
            double PileUp_Probability = PileUp[NUCLEUS][current_detector].second;
            double PileUp_Sigma = PileUp[NUCLEUS][current_detector].first;
            while (Reader->Next())
            {
                // ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "Progress: ");
                normal_distribution<double> distribution(*Simulated_Energy, resolution);
                double pileup = 0;
                if (gRandom->Uniform() < PileUp_Probability)
                {
                    pileup = abs(gRandom->Gaus(0, PileUp_Sigma));
                }
                H_Sim_Conv[NUCLEUS][current_detector]->Fill(distribution(generator) + pileup);
            }
        }
    }

    return nullptr;
}

double FunctionToMinimize(const double *par)
{
    if (H_Exp_Channel["32Ar"][current_detector]->GetEntries() == 0)
    {
        Warning(detectorName[current_detector] + "   No Entries", 1);
        return 0;
    }

    double chi2 = 0;
    // LOOP ON NUCLEUS
    for (string NUCLEUS : Nuclei)
    {
        Info("NUCLEUS : " + NUCLEUS, 1);
        ApplyResolution(VERBOSE, current_detector, NUCLEUS, par[0]);
        
        H_Exp[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][IAS[NUCLEUS]][current_detector].first, WindowsMap[NUCLEUS][IAS[NUCLEUS]][current_detector].second);
        H_Sim_Conv[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][IAS[NUCLEUS]][current_detector].first, WindowsMap[NUCLEUS][IAS[NUCLEUS]][current_detector].second);


        // ### OFFSET ADDING ON EXP FOR PERFECT COMPARISON OF IAS ### //
        if (VERBOSE == 1) Info("OFFSET ADDING", 2);
        ApplyCalibration(VERBOSE, 0, NUCLEUS);
        double offset = H_Sim_Conv[NUCLEUS][current_detector]->GetMean() - H_Exp[NUCLEUS][current_detector]->GetMean();
        Info("Correction on IAS : " + to_string(offset) + " keV", 2);
        ApplyCalibration(VERBOSE, offset, NUCLEUS);

        // ### SCALING OF THE SIMULATED SPECTRUM ### //
        if (VERBOSE == 1) Info("SCALING OF THE SIMULATED SPECTRUM", 2);
        H_Sim_Conv[NUCLEUS][current_detector]->Scale(H_Exp[NUCLEUS][current_detector]->Integral() / H_Sim_Conv[NUCLEUS][current_detector]->Integral());

        // ### CHI2 CALCULATION ### //
        if (VERBOSE == 1)
            Info("CHI2 CALCULATION", 2);
       
        H_Exp[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][IAS[NUCLEUS]][current_detector].first, WindowsMap[NUCLEUS][IAS[NUCLEUS]][current_detector].second);
        H_Sim_Conv[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][IAS[NUCLEUS]][current_detector].first, WindowsMap[NUCLEUS][IAS[NUCLEUS]][current_detector].second);
        
        gErrorIgnoreLevel = kError;
        double current_chi2 = H_Exp[NUCLEUS][current_detector]->Chi2Test(H_Sim_Conv[NUCLEUS][current_detector], "UW CHI2/NDF");
        gErrorIgnoreLevel = kInfo;

        
        
        Info("χ2 : " + to_string(current_chi2), 2);
        chi2 += current_chi2;

        // ### PLOTTING ### //
        MakeResidualPlot(NUCLEUS);
        ApplyCalibration(VERBOSE, 0, NUCLEUS);
        TH1D *H_Sim_Conv_WOBKG = (TH1D *)H_Sim_Conv[NUCLEUS][current_detector]->Clone(("H_Sim_Conv_WOBKG_" + NUCLEUS + "_" + detectorName[current_detector]).c_str());

        // Adding Beta Background and Applying threshold
        if (NUCLEUS == "32Ar")
        {
            TH1D *H_Exp_Without_Background = RemoveBKG(H_Exp[NUCLEUS][current_detector], true, NUCLEUS);
            H_Sim_Conv[NUCLEUS][current_detector]->Add((TH1D *)Background_function[NUCLEUS][current_detector]->GetHistogram(), 1);
            dir_nuclei_detector[NUCLEUS][current_detector]->cd();
            H_Exp_Without_Background->Write();

            // if (VERBOSE == 1)
            //     Info("Adding 32Cl", 2);
            // TH1D *H_Cl = ApplyResolution(VERBOSE, current_detector, "32Cl", par[0], 2, true);
            // H_Cl->Scale(0.01 / 3 * 32 / 57 * 0.8);
            // H_Sim_Conv["32Cl"][current_detector] = (TH1D *)H_Sim_Conv[NUCLEUS][current_detector]->Clone(("H_Sim_Conv_32Cl_" + detectorName[current_detector]).c_str());
            
            // H_Sim_Conv["32Cl"][current_detector]->Add((TH1D *)Background_function[NUCLEUS][current_detector]->GetHistogram(), 1);
            
            
            // H_Sim_Conv[NUCLEUS][current_detector]->Add(H_Cl);


        }
        for (int i = 0; i < H_Sim_Conv[NUCLEUS][current_detector]->GetNbinsX(); i++)
        {
            H_Sim_Conv[NUCLEUS][current_detector]->SetBinContent(i, H_Sim_Conv[NUCLEUS][current_detector]->GetBinContent(i) * Threshold_function[current_detector]->Eval(H_Sim_Conv[NUCLEUS][current_detector]->GetBinCenter(i)));
        }

        dir_nuclei_detector[NUCLEUS][current_detector]->cd();
        H_Exp[NUCLEUS][current_detector]->Write();
        H_Sim_Conv[NUCLEUS][current_detector]->Write();
        H_Sim_Conv_WOBKG->Write();
    
    }
    
    //////////////////////

    Info("<χ2> : " + to_string(chi2/Nuclei.size()), 1);
    return chi2;
}

void CHI2Minimization()
{

    //////////////////////////////////////////////////////////////////////////////

    if (ManualCalibFitted[current_detector].size() == 0)
        return;
    
    plotting = true;
    const double par2[2] = {Detector_Resolution[current_detector].first, 0};
    FunctionToMinimize(par2);
    plotting = false;
    NUCLEUS = "32Ar";

    ///// ALPHA BACKGROUND ///
    /// 18N
    // as been COMMENTED FOR 2025
    // TH1D *H_Exp_Without_Background = RemoveBKG(H_Exp["32Ar"][current_detector], true, "32Ar");
    // RemoveBKG(H_Exp["32Ar_thick"][current_detector], true, "32Ar_thick");
    // if (SIMULATED_File["18N"])
    // {
    //     double alpha_res = 0.021964;
    //     TH1D *H1 = (TH1D *)H_Sim["18N"][current_detector]->Clone();
    //     H1->Reset();
    //     for (int i = 0; i < H_Sim["18N"][current_detector]->GetEntries(); i++)
    //     {
    //         double value_corrected = (H_Sim["18N"][current_detector]->GetRandom() - ManualCalibFitted_Alpha[current_detector][0]) / ManualCalibFitted_Alpha[current_detector][1];
    //         double value = ManualCalibFitted[current_detector][0] + ManualCalibFitted[current_detector][1] * value_corrected + ManualCalibFitted[current_detector][2] * value_corrected * value_corrected;
    //         normal_distribution<double> distribution(value, 1.4969 + alpha_res * sqrt(value));
    //         H1->Fill(distribution(generator));
    //     }

    //     // FIRST PEAK
    //     TH1D *H11 = (TH1D *)H1->Clone();
    //     H_Exp_Without_Background->GetXaxis()->SetRangeUser(WindowsMap["18N"][1][current_detector].first, WindowsMap["18N"][1][current_detector].second);
    //     H11->GetXaxis()->SetRangeUser(WindowsMap["18N"][1][current_detector].first, WindowsMap["18N"][1][current_detector].second);
    //     H11->Scale(0.016);
    //     for (int bin = 0; bin < H11->GetNbinsX(); bin++)
    //     {
    //         if (H11->GetBinCenter(bin) < WindowsMap["18N"][1][current_detector].first || H11->GetBinCenter(bin) > WindowsMap["18N"][1][current_detector].second)
    //         {
    //             H11->SetBinContent(bin, 0);
    //         }
    //     }

    //     // SECOND PEAK
    //     TH1D *H12 = (TH1D *)H1->Clone();
    //     H_Exp_Without_Background->GetXaxis()->SetRangeUser(WindowsMap["18N"][2][current_detector].first, WindowsMap["18N"][2][current_detector].second);
    //     H12->GetXaxis()->SetRangeUser(WindowsMap["18N"][2][current_detector].first, WindowsMap["18N"][2][current_detector].second);
    //     H12->Scale(0.016);
    //     for (int bin = 0; bin < H12->GetNbinsX(); bin++)
    //     {
    //         if (H12->GetBinCenter(bin) < WindowsMap["18N"][2][current_detector].first || H12->GetBinCenter(bin) > WindowsMap["18N"][2][current_detector].second)
    //         {
    //             H12->SetBinContent(bin, 0);
    //         }
    //     }

        // Other peaks
    //     H1->Reset();
    //     TH1D *sub = (TH1D *)H_Exp[NUCLEUS][current_detector]->Clone();
    //     sub->Add(H_Sim_Conv[NUCLEUS][current_detector], -1);
    //     sub->GetXaxis()->SetRangeUser(-1111, -1111);
    //     TF1 *gauss = new TF1("gaus", "gaus", 1600, 2300);
    //     gauss->SetParameters(1000, 2000, 100);
    //     sub->Fit("gaus", "QRN", "", 1600, 2300);
    //     double first = gauss->GetParameter(1);
    //     for (int i = 0; i < H_Sim["18N"][current_detector]->GetEntries(); i++)
    //     {
    //         double value = H_Sim["18N"][current_detector]->GetRandom();
    //         normal_distribution<double> distribution(value, alpha_res * sqrt(value));
    //         H1->Fill(distribution(generator));
    //     }
    //     H1->GetXaxis()->SetRangeUser(-1111, -1111);
    //     H1->Fit("gaus", "QRN", "", 1600, 2300);
    //     double offset = gauss->GetParameter(1) - first;
    //     H1->Reset();
    //     for (int i = 0; i < H_Sim["18N"][current_detector]->GetEntries(); i++)
    //     {
    //         double value = H_Sim["18N"][current_detector]->GetRandom();
    //         normal_distribution<double> distribution(value, alpha_res * sqrt(value));
    //         H1->Fill(distribution(generator));
    //     }

    //     TH1D *H13 = (TH1D *)H1->Clone();
    //     H13->Scale(0.016);
    //     for (int bin = 0; bin < H13->GetNbinsX(); bin++)
    //     {
    //         if (H13->GetBinCenter(bin) < WindowsMap["18N"][2][current_detector].second)
    //         {
    //             H13->SetBinContent(bin, 0);
    //         }
    //     }
    //     // ADDING
    //     H_Sim_Conv[NUCLEUS][current_detector]->Add(H11, 1);
    //     H_Sim_Conv[NUCLEUS][current_detector]->Add(H12, 1);
    //     H_Sim_Conv[NUCLEUS][current_detector]->Add(H13, 1);
    // }
}


void Manual_Calibration(int Verbose = 0)
{
    for (string Nucleus : Nuclei)
    {
        if (Verbose == 1) Info("NUCLEUS : " + Nucleus, 2);
        Reader = new TTreeReader(MERGED_Tree_Detectors[Nucleus][current_detector]);
        Reader->Restart();
        TTreeReaderValue<double> ChannelDet(*Reader, "Channel");
        H_Exp_Channel[Nucleus][current_detector]->Reset();
        while (Reader->Next())
        {
            H_Exp_Channel[Nucleus][current_detector]->Fill(*ChannelDet / 1000);
        }
        dir_nuclei_detector[Nucleus][current_detector]->cd();
        H_Exp_Channel[Nucleus][current_detector]->Write();
    }

    int counter = 0;
    for (double peak : PeakList["32Ar"])
    {
        if (ManualCalib[peak][current_detector].first == -1 || !ManualCalib[peak][current_detector].first)
            continue;

        if (Verbose == 1) Info("Peak: " + to_string(peak), 2);
        G_Calibration[NUCLEUS][current_detector]->SetPoint(counter, ManualCalib[peak][current_detector].first, ManualCalib[peak][current_detector].second);
        counter++;
    }

    if (G_Calibration[NUCLEUS][current_detector]->GetN() == 0)
    {
        return;
    }

    // FIRST FIT WITH MANUAL
    G_Calibration[NUCLEUS][current_detector]->Fit("pol1", "Q");
    Calibration_Function[current_detector] = G_Calibration[NUCLEUS][current_detector]->GetFunction("pol1");
    ManualCalibLinear[current_detector] = make_pair(G_Calibration[NUCLEUS][current_detector]->GetFunction("pol1")->GetParameter(0), G_Calibration[NUCLEUS][current_detector]->GetFunction("pol1")->GetParameter(1));
    Info(detectorName[current_detector] + "   First Calibration : a = " + to_string(G_Calibration[NUCLEUS][current_detector]->GetFunction("pol1")->GetParameter(0)) + " b = " + 
        to_string(G_Calibration[NUCLEUS][current_detector]->GetFunction("pol1")->GetParameter(1)), 1);
}

void Fitting_Calibration(int Verbose = 0)
{
    if (Verbose == 1) Start("Fitting Calibration", 1);
    // SECOND FIT WITH GAUSSIANS IN WINDOWS
    double coef = ManualCalibLinear[current_detector].second;
    double offset = ManualCalibLinear[current_detector].first;
    TF1 *Calibration_BijectiveFunction = InvertFunction(Calibration_Function[current_detector]);

    if (H_Exp_Channel["32Ar"][current_detector]->GetEntries() == 0)
    {
        Warning(detectorName[current_detector] + "   No Entries", 1);
        return;
    }

    for (string Nucleus : Nuclei)
    {
        if (Verbose == 1) Info("NUCLEUS : " + Nucleus, 1);
        NUCLEUS = Nucleus;
        G_Calibration[NUCLEUS][current_detector] = new TGraphErrors();

        for (double peak : PeakList[Nucleus])
        {

            if (WindowsMap[NUCLEUS][peak][current_detector].first == -1 || !WindowsMap[NUCLEUS][peak][current_detector].first)
                continue;

            if (find(CalibrationPeaks[NUCLEUS].begin(), CalibrationPeaks[NUCLEUS].end(), peak) == CalibrationPeaks[NUCLEUS].end())
                continue;
            

            if (Verbose == 1) Info("Peak: " + to_string(peak), 2);
            

            // #### FITTING THE EXPERIMENTAL SPECTRUM #### //
            if (WindowsMap[NUCLEUS][peak][current_detector].first < 1000 && NUCLEUS.find("32Ar") != string::npos) // taking into account the exponential bkg
            {
                F_Peak_Exp[NUCLEUS][peak][current_detector] = new TF1((detectorName[current_detector] + "_Exp" + to_string(peak) + NUCLEUS).c_str(), "gaus(0) + [3]*exp([4]*x)", 0, 6500);
                F_Peak_Exp[NUCLEUS][peak][current_detector]->SetParLimits(1, Calibration_BijectiveFunction->Eval(WindowsMap[NUCLEUS][peak][current_detector].first), Calibration_BijectiveFunction->Eval(WindowsMap[NUCLEUS][peak][current_detector].second));
                F_Peak_Exp[NUCLEUS][peak][current_detector]->SetParLimits(3, 0, 10000);
                F_Peak_Exp[NUCLEUS][peak][current_detector]->SetParLimits(4, -1, 0);
            }
            else
            {
                F_Peak_Exp[NUCLEUS][peak][current_detector] = new TF1((detectorName[current_detector] + "_Exp" + to_string(peak) + NUCLEUS).c_str(), "gaus", 0, 6500);
                F_Peak_Exp[NUCLEUS][peak][current_detector]->SetParLimits(1, Calibration_BijectiveFunction->Eval(WindowsMap[NUCLEUS][peak][current_detector].first), Calibration_BijectiveFunction->Eval(WindowsMap[NUCLEUS][peak][current_detector].second));
            }

            if (H_Exp_Channel[NUCLEUS][current_detector]->Integral() == 0)
                continue;
            TFitResultPtr r_exp = H_Exp_Channel[NUCLEUS][current_detector]->Fit((detectorName[current_detector] + "_Exp" + to_string(peak) + NUCLEUS).c_str(), "QRNS", "", Calibration_BijectiveFunction->Eval(WindowsMap[NUCLEUS][peak][current_detector].first), Calibration_BijectiveFunction->Eval(WindowsMap[NUCLEUS][peak][current_detector].second));
            if (r_exp != 0)
                Warning("<Exp> Fit failed for peak " + to_string(peak), 2);
            //

            // #### FITTING THE SIMULATED SPECTRUM #### //
            TH1D* H_Simulation;
            if (!Resolution_applied)
                H_Simulation = (TH1D *)H_Sim[NUCLEUS][current_detector]->Clone();
            else
                H_Simulation = (TH1D *)H_Sim_Conv[NUCLEUS][current_detector]->Clone();
            
            H_Simulation->Rebin(10);

            
            F_Peak_Sim[NUCLEUS][peak][current_detector] = new TF1((detectorName[current_detector] + "_Sim" + to_string(peak) + NUCLEUS).c_str(), "gaus(0)", 0, 6500);
            F_Peak_Sim[NUCLEUS][peak][current_detector]->SetParameter(1, (WindowsMap[NUCLEUS][peak][current_detector].first + WindowsMap[NUCLEUS][peak][current_detector].second) / 2);
            F_Peak_Sim[NUCLEUS][peak][current_detector]->SetParLimits(1, WindowsMap[NUCLEUS][peak][current_detector].first, WindowsMap[NUCLEUS][peak][current_detector].second);
            TFitResultPtr r_sim = H_Simulation->Fit((detectorName[current_detector] + "_Sim" + to_string(peak) + NUCLEUS).c_str(), "QRN", "", WindowsMap[NUCLEUS][peak][current_detector].first, WindowsMap[NUCLEUS][peak][current_detector].second);
            if (r_sim != 0)
                Warning("<SIM> Fit failed in " + NUCLEUS + " for peak " + to_string(peak), 2);

            // #### ADDING TO THE CALIBRATION GRAPH #### //
            if (r_exp == 0 && r_sim == 0)
            {
                if (peak == IAS[NUCLEUS])
                {
                    H_Exp_Channel[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(Calibration_BijectiveFunction->Eval(WindowsMap[NUCLEUS][peak][current_detector].first), Calibration_BijectiveFunction->Eval(WindowsMap[NUCLEUS][peak][current_detector].second));
                    H_Simulation->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][peak][current_detector].first, WindowsMap[NUCLEUS][peak][current_detector].second);
                    G_Calibration[NUCLEUS][current_detector]->AddPoint(H_Exp_Channel[NUCLEUS][current_detector]->GetMean(), H_Simulation->GetMean());
                    G_Calibration[NUCLEUS][current_detector]->SetPointError(G_Calibration[NUCLEUS][current_detector]->GetN() - 1, H_Exp_Channel[NUCLEUS][current_detector]->GetMeanError(), sqrt(pow(H_Simulation->GetMeanError(), 2) + pow(PeakData[NUCLEUS][peak][1], 2)));
                    H_Exp_Channel[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(-1111, -1111);
                    H_Simulation->GetXaxis()->SetRangeUser(-1111, -1111);
                }
                else
                {
                    G_Calibration[NUCLEUS][current_detector]->AddPoint(F_Peak_Exp[NUCLEUS][peak][current_detector]->GetParameter(1), F_Peak_Sim[NUCLEUS][peak][current_detector]->GetParameter(1));
                    G_Calibration[NUCLEUS][current_detector]->SetPointError(G_Calibration[NUCLEUS][current_detector]->GetN() - 1, F_Peak_Exp[NUCLEUS][peak][current_detector]->GetParError(1), sqrt(pow(F_Peak_Sim[NUCLEUS][peak][current_detector]->GetParError(1), 2) + pow(PeakData[NUCLEUS][peak][1], 2)));
                }
            }

            // displaying the gaussian fit
            TCanvas *c2 = new TCanvas((detectorName[current_detector] + "_" + NUCLEUS + "_" + to_string(peak)).c_str(), (detectorName[current_detector] + "_" + NUCLEUS + "_" + to_string(peak)).c_str(), 1920, 1080);
            c2->Divide(2, 1);
            c2->cd(1);
            TH1D *H = (TH1D *)H_Exp_Channel[NUCLEUS][current_detector]->Clone();
            H->GetXaxis()->SetRangeUser(Calibration_BijectiveFunction->Eval(WindowsMap[NUCLEUS][peak][current_detector].first), Calibration_BijectiveFunction->Eval(WindowsMap[NUCLEUS][peak][current_detector].second));
            H->Draw("HIST");
            F_Peak_Exp[NUCLEUS][peak][current_detector]->SetNpx(eSiliN_cal);
            F_Peak_Exp[NUCLEUS][peak][current_detector]->SetLineColor(kRed);
            F_Peak_Exp[NUCLEUS][peak][current_detector]->Draw("SAME");
            c2->cd(2);
            TH1D *H1 = (TH1D *)H_Simulation->Clone();
            H1->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][peak][current_detector].first, WindowsMap[NUCLEUS][peak][current_detector].second);
            H1->Draw("HIST");
            F_Peak_Sim[NUCLEUS][peak][current_detector]->SetNpx(eSiliN_cal);
            F_Peak_Sim[NUCLEUS][peak][current_detector]->SetLineColor(kRed);
            F_Peak_Sim[NUCLEUS][peak][current_detector]->Draw("SAME L ");
            dir_peak_detector[current_detector]->cd();
            c2->Write();
        }
    }

    /////// FITTING ALL NUCLEI //////
    MG_Global_Calibration[current_detector] = new TMultiGraph();
    int counter = 0;
    int total_points = 0;
    for (string Nucleus : Nuclei)
    {
        G_Calibration[Nucleus][current_detector]->SetMarkerColor(NucleusColor[Nucleus]);
        G_Calibration[Nucleus][current_detector]->SetMarkerSize(1);
        G_Calibration[Nucleus][current_detector]->SetMarkerStyle(20);
        MG_Global_Calibration[current_detector]->Add(G_Calibration[Nucleus][current_detector]);
        
        counter ++;
        total_points += G_Calibration[Nucleus][current_detector]->GetN();
    }

    if (total_points == 0)
    {
        Warning("No points to fit", 1);
        return;
    }
    

    TCanvas *c1 = new TCanvas((detectorName[current_detector] + "_Calibration").c_str(), (detectorName[current_detector] + "_Calibration").c_str(), 1920, 1080);
    c1->Divide(1, 2);
    c1->cd(1);
    string function = "pol2";
    TF1 *fit = new TF1(function.c_str(), function.c_str(), 0, 6500);
   
    TFitResultPtr r = MG_Global_Calibration[current_detector]->Fit(function.c_str(), "QS");
    Calibration_Result[current_detector] = r;
    if (function == "pol2")
    {
        Info("Calibration : a = " + to_string(MG_Global_Calibration[current_detector]->GetFunction(function.c_str())->GetParameter(0)) + 
            " b = " + to_string(MG_Global_Calibration[current_detector]->GetFunction(function.c_str())->GetParameter(1)) + 
            " c = " + to_string(MG_Global_Calibration[current_detector]->GetFunction(function.c_str())->GetParameter(2)) + 
            "  χ2 = " + to_string(MG_Global_Calibration[current_detector]->GetFunction(function.c_str())->GetChisquare() / MG_Global_Calibration[current_detector]->GetFunction(function.c_str())->GetNDF()), 2);
        
        G_Calibration_FitParameters[0]->AddPoint(current_detector, MG_Global_Calibration[current_detector]->GetFunction(function.c_str())->GetParameter(0));
        G_Calibration_FitParameters[1]->AddPoint(current_detector, MG_Global_Calibration[current_detector]->GetFunction(function.c_str())->GetParameter(1));
        G_Calibration_FitParameters[2]->AddPoint(current_detector, MG_Global_Calibration[current_detector]->GetFunction(function.c_str())->GetParameter(2));
    
    }
    else if (function == "pol1")
    {
        Info("Calibration : a = : a = " + to_string(MG_Global_Calibration[current_detector]->GetFunction(function.c_str())->GetParameter(0)) + 
            " b = " + to_string(MG_Global_Calibration[current_detector]->GetFunction(function.c_str())->GetParameter(1)) + 
            "  χ2 = " + to_string(MG_Global_Calibration[current_detector]->GetFunction(function.c_str())->GetChisquare() / MG_Global_Calibration[current_detector]->GetFunction(function.c_str())->GetNDF()), 2);

        G_Calibration_FitParameters[0]->AddPoint(current_detector, MG_Global_Calibration[current_detector]->GetFunction(function.c_str())->GetParameter(0));
        G_Calibration_FitParameters[1]->AddPoint(current_detector, MG_Global_Calibration[current_detector]->GetFunction(function.c_str())->GetParameter(1));
    }
    else
    {
        Error("Unknown function type for calibration");
    }
    MG_Global_Calibration[current_detector]->GetFunction(function.c_str())->SetLineColor(kRed);
    MG_Global_Calibration[current_detector]->Draw("AP");
    MG_Global_Calibration[current_detector]->GetXaxis()->SetTitle("Channel");
    MG_Global_Calibration[current_detector]->GetYaxis()->SetTitle("Energy [keV]");
    TF1 *linear = new TF1("pol1", "pol1", 0, 6500);
    TMultiGraph *g = (TMultiGraph *)MG_Global_Calibration[current_detector]->Clone();
    g->Fit("pol1", "QN");
    linear->SetLineColor(kGreen);
    linear->Draw("SAME");
    c1->cd(2);
    // draw diff between fit and data
    TMultiGraph *MG_Global_Calibration_diff = new TMultiGraph();

    map<string, TGraphErrors *> G_Calibration_diff;
    for (string Nucleus : Nuclei)
    {
        G_Calibration_diff[Nucleus] = (TGraphErrors *)G_Calibration[Nucleus][current_detector]->Clone();
        for (int i = 0; i < G_Calibration_diff[Nucleus]->GetN(); i++)
        {
            double x, y;
            G_Calibration_diff[Nucleus]->GetPoint(i, x, y);
            G_Calibration_diff[Nucleus]->SetPoint(i, x, y - MG_Global_Calibration[current_detector]->GetFunction(function.c_str())->Eval(x));
        }
        MG_Global_Calibration_diff->Add(G_Calibration_diff[Nucleus]);
    }

    MG_Global_Calibration_diff->Draw("AP");
    MG_Global_Calibration_diff->GetXaxis()->SetTitle("Channel");
    MG_Global_Calibration_diff->GetYaxis()->SetTitle("Energy [keV]");
    TF1 *f = new TF1("f", "0", 0, 6500);
    f->SetLineColor(kRed);
    f->Draw("SAME");

    // TF1 *fit2 = new TF1("fit2", "[1]-[0] + ([3]-[2]) * x - [4] * x * x", 0, 6500);
    // fit2->SetParameters(MG_Global_Calibration[current_detector]->GetFunction(function.c_str())->GetParameter(0), linear->GetParameter(0), MG_Global_Calibration[current_detector]->GetFunction(function.c_str())->GetParameter(1), linear->GetParameter(1), MG_Global_Calibration[current_detector]->GetFunction(function.c_str())->GetParameter(2));
    // fit2->SetLineColor(kGreen);
    // fit2->Draw("SAME");

    // TPAVETEXT for chi2
    TPaveText *pt = new TPaveText(0.1, 0.1, 0.3, 0.3, "NDC");
    pt->SetTextSize(0.05);
    pt->AddText(("#chi^{2}_{#nu} = " + to_string(MG_Global_Calibration[current_detector]->GetFunction(function.c_str())->GetChisquare() / MG_Global_Calibration[current_detector]->GetFunction(function.c_str())->GetNDF())).c_str());
    pt->Draw("SAME");

    TLegend *l = new TLegend(0.1, 0.7, 0.3, 0.9);
    for (string Nucleus : Nuclei)
    {
        l->AddEntry(G_Calibration[Nucleus][current_detector], Nucleus.c_str(), "p");
    }
    l->Draw("SAME");

    dir_detector[current_detector]->cd();
    c1->Write();

    Calibration_Function[current_detector] = MG_Global_Calibration[current_detector]->GetFunction(function.c_str());
    vector<double> vec = {};//{MG_Global_Calibration[current_detector]->GetFunction("pol2")->GetParameter(0), MG_Global_Calibration[current_detector]->GetFunction("pol2")->GetParameter(1), MG_Global_Calibration[current_detector]->GetFunction("pol2")->GetParameter(2)};
    for (int i = 0; i < Calibration_Function[current_detector]->GetNpar(); i++)
        vec.push_back(Calibration_Function[current_detector]->GetParameter(i));
    
    ManualCalibFitted[current_detector] = vec;

    if (Verbose == 1) Info("Background Subtraction", 2);
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    // BACKGROUND ///////////////////////////////////////////////////////////////////////////////////////
    /*
    NUCLEUS = "32Ar";
    counter = 0;
    Reader = new TTreeReader(MERGED_Tree_Detectors[NUCLEUS][current_detector]);
    Reader->Restart();
    H_Exp[NUCLEUS][current_detector]->Reset();
    TTreeReaderValue<double> ChannelDet1(*Reader, "Channel");
    while (Reader->Next())
    {
        double value = *ChannelDet1 / 1000;
        H_Exp[NUCLEUS][current_detector]->Fill(Calibration_Function[current_detector]->Eval(value));
    }

    RemoveBKG(H_Exp[NUCLEUS][current_detector], false, NUCLEUS);

    // NUCLEUS = "32Ar_thick";
    // counter = 0;
    // Reader = new TTreeReader(MERGED_Tree_Detectors[NUCLEUS][current_detector]);
    // Reader->Restart();
    // H_Exp[NUCLEUS][current_detector]->Reset();
    // while (Reader->Next())
    // {
    //     double value = *ChannelDet1 / 1000;
    //     H_Exp[NUCLEUS][current_detector]->Fill(Calibration_Function[current_detector]->Eval(value));
    // }

    // RemoveBKG(H_Exp[NUCLEUS][current_detector], false, NUCLEUS);
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// ALPHA
    for (int peak = 1; peak < CanvasMap["18N"].first * CanvasMap["18N"].second; peak++)
    {
        if (WindowsMap["18N"][peak][current_detector].first == -1 || !WindowsMap["18N"][peak][current_detector].first)
            continue;

        double alpha_low = 0;
        double alpha_high = 0;
        double alpha = 0;
        double alpha_error = 0;
        if (peak == 1)
        {
            if (GetDetector(current_detector) <= 4)
            {
                alpha_low = Alpha1_Up[GetDetectorChannel(current_detector)].first;
                alpha_high = Alpha1_Up[GetDetectorChannel(current_detector)].second;
                alpha = alpha1;
                alpha_error = alpha1_error;
            }
            else
            {
                alpha_low = Alpha1_Down[GetDetectorChannel(current_detector)].first;
                alpha_high = Alpha1_Down[GetDetectorChannel(current_detector)].second;
                alpha = alpha1;
                alpha_error = alpha1_error;
            }
        }
        else if (peak == 2)
        {
            if (GetDetector(current_detector) <= 4)
            {
                alpha_low = Alpha2_Up[GetDetectorChannel(current_detector)].first;
                alpha_high = Alpha2_Up[GetDetectorChannel(current_detector)].second;
                alpha = alpha2;
                alpha_error = alpha2_error;
            }
            else
            {
                alpha_low = Alpha2_Down[GetDetectorChannel(current_detector)].first;
                alpha_high = Alpha2_Down[GetDetectorChannel(current_detector)].second;
                alpha = alpha2;
                alpha_error = alpha2_error;
            }
        }

        // H_Sim["18N"][current_detector]->Write();
        TF1 *gauss2 = new TF1("gauss2", "gaus");
        gauss2->SetParameter(1, (WindowsMap["18N"][peak][current_detector].second - WindowsMap["18N"][peak][current_detector].first) / 2);
        H_Sim["18N"][current_detector]->Fit("gauss2", "QRN", "", alpha_low, alpha_high);

        TF1 *Bijective_Function = InvertFunction(Calibration_Function[current_detector]);
        double x = Bijective_Function->Eval(alpha);
        // double x = (-ManualCalibFitted[current_detector][1] + sqrt(pow(ManualCalibFitted[current_detector][1], 2) - 4 * ManualCalibFitted[current_detector][2] * (ManualCalibFitted[current_detector][0] - alpha))) / (2 * ManualCalibFitted[current_detector][2]);
        G_Calibration_Alpha[current_detector]->SetPoint(counter, x, gauss2->GetParameter(1));
        G_Calibration_Alpha[current_detector]->SetPointError(counter, (alpha_error / coef), gauss2->GetParError(1));
        counter++;
    }

    TCanvas *c2 = new TCanvas((detectorName[current_detector] + "_Calibration_Alpha").c_str(), (detectorName[current_detector] + "_Calibration_Alpha").c_str(), 1920, 1080);
    c2->Divide(1, 2);
    c2->cd(1);
    function = "pol1";
    TF1 *fit_alpha = new TF1(function.c_str(), "[0] + [1] * x", 0, 6500);
    G_Calibration_Alpha[current_detector]->Fit(function.c_str(), "Q");
    G_Calibration_Alpha[current_detector]->GetFunction(function.c_str())->SetLineColor(kRed);
    G_Calibration_Alpha[current_detector]->Draw("*AP");
    G_Calibration_Alpha[current_detector]->GetXaxis()->SetTitle("Channel");
    G_Calibration_Alpha[current_detector]->GetYaxis()->SetTitle("Energy [keV]");
    c2->cd(2);
    // draw diff between fet and data
    TGraphErrors *G_Calibration_diff_alpha = (TGraphErrors *)G_Calibration_Alpha[current_detector]->Clone();
    for (int i = 0; i < G_Calibration_diff_alpha->GetN(); i++)
    {
        double x, y;
        G_Calibration_diff_alpha->GetPoint(i, x, y);
        G_Calibration_diff_alpha->SetPoint(i, x, y - G_Calibration_Alpha[current_detector]->GetFunction(function.c_str())->Eval(x));
    }
    G_Calibration_diff_alpha->Draw("*AP");
    G_Calibration_diff_alpha->GetXaxis()->SetTitle("Channel");
    G_Calibration_diff_alpha->GetYaxis()->SetTitle("Energy [keV]");
    TLine *line_alpha = new TLine(G_Calibration_diff_alpha->GetXaxis()->GetXmin(), 0, G_Calibration_diff_alpha->GetXaxis()->GetXmax(), 0);
    line_alpha->SetLineColor(kRed);
    line_alpha->Draw();

    dir_detector[current_detector]->cd();
    c2->Write();

    vector<double> vec_alpha = {G_Calibration_Alpha[current_detector]->GetFunction("pol1")->GetParameter(0), G_Calibration_Alpha[current_detector]->GetFunction("pol1")->GetParameter(1)};

    // cout << "ALPHA : " << vec_alpha[0] << " " << vec_alpha[1] << endl;
    //  vector<double> vec_alpha = {-5.5938, 7.06253e-2*1000};

    ManualCalibFitted_Alpha[current_detector] = vec_alpha;
    */
}

double CalibrationError(int det, double energy)
{
   TMatrixDSym* covMatrix = new TMatrixDSym(Calibration_Result[det]->GetCovarianceMatrix());
   TVectorD grad(Calibration_Function[det]->GetNpar());
   Calibration_Function[det]->GradientPar(&energy, grad.GetMatrixArray());
   if (!covMatrix) {
      double variance = 0;
      for(int iPar = 0; iPar < Calibration_Function[det]->GetNpar(); iPar++) {
         variance += grad(iPar)*grad(iPar)*Calibration_Function[det]->GetParError(iPar)*Calibration_Function[det]->GetParError(iPar);
      }
      return sqrt(variance);
   }
   return sqrt(covMatrix->Similarity(grad));
}

void FittingResolution(int Verbose = 0)
{

    if (Verbose == 1)
        Start("Fitting Resolution", 1);

    Detector_ElectronicResolution[current_detector] = Calibration_Function[current_detector]->Eval(Detector_ElectronicResolution[current_detector]) - Calibration_Function[current_detector]->GetParameter(0);
    G_Electronic_Resolution_FitParameter->AddPoint(current_detector, Detector_ElectronicResolution[current_detector]);

    if (!Fitting_Resolution_FLAG)
    {
        G_Resolution_FitParameter->AddPoint(current_detector, Detector_Resolution[current_detector].first);
        return;
    }

    NUCLEUS = "33Ar";
    int peak = 21;
    G_Resolution_Peak[current_detector] = new TGraphErrors();
    double start = 2.0;
    double end = 7.0;
    double step = 0.1;  
    int N = (int)((end - start) / step) + 1;
    int counter = 0;
    double minimum_chi2 = MAXFLOAT;
    for (double Res_value = start; Res_value <= end; Res_value += step)
    {
        if (Verbose == 1)
            Info("Try Resolution : " + to_string(Res_value) + " keV", 3);

        ProgressCounter(counter, N, "Resolution Fitting");

        ApplyResolution(VERBOSE, current_detector, NUCLEUS, Res_value, 4);

        // ### OFFSET ADDING ON EXP FOR PERFECT COMPARISON OF IAS ### //
        if (VERBOSE == 1) Info("OFFSET ADDING", 4);
        ApplyCalibration(VERBOSE, 0, NUCLEUS, 4);
        H_Exp[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][peak][current_detector].first, WindowsMap[NUCLEUS][peak][current_detector].second);
        H_Sim_Conv[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][peak][current_detector].first, WindowsMap[NUCLEUS][peak][current_detector].second);
        double offset = H_Sim_Conv[NUCLEUS][current_detector]->GetMean() - H_Exp[NUCLEUS][current_detector]->GetMean();
        ApplyCalibration(VERBOSE, offset, NUCLEUS, 4);
        H_Exp[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][peak][current_detector].first, WindowsMap[NUCLEUS][peak][current_detector].second);
        H_Sim_Conv[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][peak][current_detector].first, WindowsMap[NUCLEUS][peak][current_detector].second);

        // ### SCALING OF THE SIMULATED SPECTRUM ### //
        if (VERBOSE == 1) Info("SCALING OF THE SIMULATED SPECTRUM", 4);
        H_Sim_Conv[NUCLEUS][current_detector]->Scale(H_Exp[NUCLEUS][current_detector]->Integral() / H_Sim_Conv[NUCLEUS][current_detector]->Integral());

        // ### CHI2 CALCULATION ### //
        if (VERBOSE == 1)
            Info("CHI2 CALCULATION", 4);
       
        gErrorIgnoreLevel = kError;
        double current_chi2 = H_Exp[NUCLEUS][current_detector]->Chi2Test(H_Sim_Conv[NUCLEUS][current_detector], "UW CHI2/NDF");
        gErrorIgnoreLevel = kInfo;


        G_Resolution_Peak[current_detector]->AddPoint(Res_value, current_chi2);
        if (minimum_chi2 > current_chi2)
        {
            Detector_Resolution[current_detector].first = Res_value;
            minimum_chi2 = current_chi2;
        }
        if (VERBOSE == 1) Info("Chi2/NDF : " + to_string(current_chi2), 4);
        counter++;
    }



    dir_detector[current_detector]->cd();
    TCanvas *c = new TCanvas((detectorName[current_detector] + "_Resolution").c_str(), (detectorName[current_detector] + "_Resolution").c_str(), 1920, 1080);
    c->cd();
    G_Resolution_Peak[current_detector]->SetMarkerStyle(20);
    G_Resolution_Peak[current_detector]->SetTitle((detectorName[current_detector] + " Resolution").c_str());
    G_Resolution_Peak[current_detector]->GetXaxis()->SetTitle("Resolution [keV]");
    G_Resolution_Peak[current_detector]->GetYaxis()->SetTitle("Chi2/NDF");
    G_Resolution_Peak[current_detector]->Draw("AP");
    c->Write();

    double err = 0;
    ///////// ##### LOOKING FOR THE MINIMUM CHI2 ON THE GRAPH AND CHI2+1 ###### /////////
    // search for the 2 neasrest points to the minimum with minimum chi2 + 1
    // double min_chi2 = minimum_chi2;
    // double min_resolution = Detector_Resolution[current_detector].first;
    // double min_resolution_1 = 0;
    // double delta_chi2_1 = MAXFLOAT;
    // double min_resolution_2 = 0;
    // double delta_chi2_2 = MAXFLOAT;
    // for (int i = 0; i < G_Resolution_Peak[current_detector]->GetN(); i++)
    // {
    //     double res_value, chi2;
    //     G_Resolution_Peak[current_detector]->GetPoint(i, res_value, chi2);
    //     if (abs(min_chi2+1-chi2) < delta_chi2_1 && Detector_Resolution[current_detector].first > res_value)
    //     {
    //         delta_chi2_1 = abs(min_chi2+1-chi2);
    //         min_resolution_1 = res_value;   
    //     }
    //     else if (abs(min_chi2+1-chi2) < delta_chi2_2 && Detector_Resolution[current_detector].first < res_value)
    //     {
    //         delta_chi2_2 = abs(min_chi2+1-chi2);
    //         min_resolution_2 = res_value;
    //     }
    // }
    
    // // keeping the maximal difference resolution between the two points with the minimum ffound
    // Detector_Resolution[current_detector].second = max(abs(min_resolution_1 - Detector_Resolution[current_detector].first), abs(min_resolution_2 - Detector_Resolution[current_detector].first));
    /////////////////////////////////////////////////////////////////////////////////////

    //////// ##### LOOKING FOR THE MINIMUM CHI2 WITH FITTED FUNCTION pol3 ###### /////////
    // searching the minimum with the fucntion
    G_Resolution_Peak[current_detector]->Fit("pol3", "Q");
    TF1 *f = G_Resolution_Peak[current_detector]->GetFunction("pol3");
    Detector_Resolution[current_detector].first = f->GetMinimumX(start, end);


    //searching error bar on found value with chi2+1
    double chi2_at_minimum = f->Eval(Detector_Resolution[current_detector].first);
    TF1 *froot = new TF1("froot", "abs([0] + [1]*x + [2]*x*x + [3]*x*x*x)", start, end);
    froot->SetParameters(f->GetParameter(0) - (chi2_at_minimum+1), f->GetParameter(1), f->GetParameter(2), f->GetParameter(3));

    // seaching the roots in the range
    double x1 = froot->GetMinimumX(start, Detector_Resolution[current_detector].first);
    double x2 = froot->GetMinimumX(Detector_Resolution[current_detector].first, end);

    Detector_Resolution[current_detector].second = max(abs(x1 - Detector_Resolution[current_detector].first), abs(x2 - Detector_Resolution[current_detector].first));
    /////////////////////////////////////////////////////////////////////////////////////
    
    Info("Resolution : " + to_string(Detector_Resolution[current_detector].first) + " ± " + to_string(Detector_Resolution[current_detector].second) + " keV", 2);
    G_Resolution_FitParameter->AddPoint(current_detector, Detector_Resolution[current_detector].first);
    
    VERBOSE = 0;
    ApplyCalibration(VERBOSE, 0, NUCLEUS, 4);
}

void PlottingWindows(int first, int last)
{
    Info("Plotting Windows");
    for (string nucleus : Nuclei)
    {
        for (int i = first; i < last; i++)
        {
            if (IsDetectorSiliStrip(i))
            {
                if (!H_Exp[nucleus][i] || !H_Sim_Conv[nucleus][i])
                    continue;
                
                TCanvas *c = new TCanvas((nucleus + "_" + detectorName[i]).c_str(), (nucleus + "_" + detectorName[i]).c_str(), 1920, 1080);
                c->Divide(CanvasMap[nucleus].first, CanvasMap[nucleus].second);
                int cdcounter = 1;
                for (double peak : PeakList[nucleus])
                {
                    c->cd(cdcounter);

                    if (WindowsMap[nucleus][peak][i].first == -1 || !WindowsMap[nucleus][peak][i].first)
                        continue;
                    TH1D *H_Copy = (TH1D *)H_Exp[nucleus][i]->Clone();
                    gPad->SetLogy();
                    H_Copy->SetTitle((" Peak " + to_string(peak)).c_str());
                    H_Copy->GetXaxis()->SetRangeUser(WindowsMap[nucleus][peak][i].first, WindowsMap[nucleus][peak][i].second);
                    H_Copy->SetStats(false);
                    H_Copy->Draw("HIST");

                    TH1D *H1_Copy = (TH1D *)H_Sim_Conv[nucleus][i]->Clone();
                    H1_Copy->GetXaxis()->SetRangeUser(WindowsMap[nucleus][peak][i].first, WindowsMap[nucleus][peak][i].second);
                    H1_Copy->SetLineColor(kRed);
                    H1_Copy->Draw("HIST SAME");
                    cdcounter++;
                }

                c->cd();
                TPad *FullPad = new TPad("FullPad", "FullPad", 0, 0, 1, 0.2);
                FullPad->Draw();
                FullPad->cd();
                H_Exp[nucleus][i]->GetXaxis()->SetRangeUser(0, 6500);
                H_Exp[nucleus][i]->GetYaxis()->SetRangeUser(1, -1111);
                H_Sim_Conv[nucleus][i]->GetXaxis()->SetRangeUser(0, 6500);
                H_Sim_Conv[nucleus][i]->GetYaxis()->SetRangeUser(1, -1111);
                H_Exp[nucleus][i]->SetStats(false);
                H_Exp[nucleus][i]->Draw("HIST");
                H_Sim_Conv[nucleus][i]->SetLineColor(kRed);
                H_Sim_Conv[nucleus][i]->Draw("HIST SAME");
                FullPad->SetLogy();
                dir_detector[i]->cd();
                c->Write();

                delete c;
                delete H_Sim[nucleus][i];
            }
        }
    }
}

void PlottingSummed(int first, int last)
{
    CALIBRATED_File->cd();
    Info("Summing Histograms");

    for (string Nucleus : Nuclei)
    {
        if (VERBOSE == 1) Info("Nucleus: " + Nucleus);
        for (int i = first; i < last; i++)
        {
            
            if (IsDetectorSiliStrip(i))
            {
                if (VERBOSE == 1) Info(detectorName[i]);
                if (H_Exp_All[Nucleus]->GetEntries() == 0)
                {
                    H_Exp_All[Nucleus] = (TH1D*)H_Exp[Nucleus][i]->Clone((Nucleus + "_Exp_All").c_str());
                    H_Sim_All[Nucleus] = (TH1D*)H_Sim_Conv[Nucleus][i]->Clone((Nucleus + "_Sim_All").c_str());
                }
                else
                {
                    H_Exp_All[Nucleus]->Add(H_Exp[Nucleus][i], 1.);
                    H_Sim_All[Nucleus]->Add(H_Sim_Conv[Nucleus][i], 1.);
                }

                if (GetDetector(i) < 5)
                {
                    if (H_Exp_All_Up[Nucleus]->GetEntries() == 0)
                    {
                        H_Exp_All_Up[Nucleus] = (TH1D*)H_Exp[Nucleus][i]->Clone((Nucleus + "_Exp_All_Up").c_str());
                        H_Sim_All_Up[Nucleus] = (TH1D*)H_Sim_Conv[Nucleus][i]->Clone((Nucleus + "_Sim_All_Up").c_str());
                    }
                    else
                    {
                        H_Exp_All_Up[Nucleus]->Add(H_Exp[Nucleus][i], 1.);
                        H_Sim_All_Up[Nucleus]->Add(H_Sim_Conv[Nucleus][i], 1.);
                    }
                }
                else
                {
                    if (H_Exp_All_Down[Nucleus]->GetEntries() == 0)
                    {
                        H_Exp_All_Down[Nucleus] = (TH1D*)H_Exp[Nucleus][i]->Clone((Nucleus + "_Exp_All_Down").c_str());
                        H_Sim_All_Down[Nucleus] = (TH1D*)H_Sim_Conv[Nucleus][i]->Clone((Nucleus + "_Sim_All_Down").c_str());
                    }
                    else
                    {
                        H_Exp_All_Down[Nucleus]->Add(H_Exp[Nucleus][i], 1.);
                        H_Sim_All_Down[Nucleus]->Add(H_Sim_Conv[Nucleus][i], 1.);
                    }
                }
            }
        }

        TCanvas *c1 = new TCanvas((Nucleus + "_SUM").c_str(), (Nucleus + "_SUM").c_str(), 800, 600);
        c1->cd();
        H_Exp_All[Nucleus]->Draw("HIST");
        H_Sim_All[Nucleus]->SetLineColor(kRed);
        H_Sim_All[Nucleus]->Draw("HIST SAME");
        c1->Write();

        if (first < 50)
        {
            TCanvas *c2 = new TCanvas((Nucleus + "_SUM_Up").c_str(), (Nucleus + "_SUM_Up").c_str(), 800, 600);
            c2->cd();
            H_Exp_All_Up[Nucleus]->Draw("HIST");
            H_Sim_All_Up[Nucleus]->SetLineColor(kRed);
            H_Sim_All_Up[Nucleus]->Draw("HIST SAME");
            c2->Write();
        }


        if (last > 50)
        {
            TCanvas *c3 = new TCanvas((Nucleus + "_SUM_Down").c_str(), (Nucleus + "_SUM_Down").c_str(), 800, 600);
            c3->cd();
            H_Exp_All_Down[Nucleus]->Draw("HIST");
            H_Sim_All_Down[Nucleus]->SetLineColor(kRed);
            H_Sim_All_Down[Nucleus]->Draw("HIST SAME");
            c3->Write();
        }
    }
}

void PlottingFitParameters(int first, int last)
{
    Info("Plotting Fit Parameters");

    /// CALIBRATION PARAMETERS
    TCanvas *c1 = new TCanvas("Calibration_Fit_Parameters", "Calibration_Fit_Parameters", 1920, 1080);
    if (G_Calibration_FitParameters[2]->GetN() == 0)
        c1->Divide(1, 2);
    else
        c1->Divide(1, 3);
    
    
    c1->cd(1);
    G_Calibration_FitParameters[0]->SetMarkerStyle(20);
    G_Calibration_FitParameters[0]->SetTitle("Calibration Fit Parameters");
    G_Calibration_FitParameters[0]->GetXaxis()->SetTitle("Detector");
    G_Calibration_FitParameters[0]->GetYaxis()->SetTitle("a (keV)");
    G_Calibration_FitParameters[0]->Draw("AP");

    c1->cd(2);
    G_Calibration_FitParameters[1]->SetMarkerStyle(20);
    G_Calibration_FitParameters[1]->SetTitle("Calibration Fit Parameters");
    G_Calibration_FitParameters[1]->GetXaxis()->SetTitle("Detector");
    G_Calibration_FitParameters[1]->GetYaxis()->SetTitle("b (keV/CH)");
    G_Calibration_FitParameters[1]->Draw("AP");

    if (G_Calibration_FitParameters[2]->GetN() != 0)
    {
        c1->cd(3);
        G_Calibration_FitParameters[2]->SetMarkerStyle(20);
        G_Calibration_FitParameters[2]->SetTitle("Calibration Fit Parameters");
        G_Calibration_FitParameters[2]->GetXaxis()->SetTitle("Detector");
        G_Calibration_FitParameters[2]->GetYaxis()->SetTitle("c (keV/CH^2)");
        G_Calibration_FitParameters[2]->Draw("AP");
    }

    c1->Write();


    /// RESOLUTION
    TCanvas *c2 = new TCanvas("Resolution_Fit_Parameters", "Resolution_Fit_Parameters", 1920, 1080);
    TMultiGraph *MG_Resolution_FitParameters = new TMultiGraph();
    c2->cd();
    TLegend *legend = new TLegend(0.1, 0.7, 0.3, 0.9);
    
    G_Electronic_Resolution_FitParameter->SetMarkerStyle(20);
    G_Electronic_Resolution_FitParameter->SetMarkerColor(kBlue);
    G_Electronic_Resolution_FitParameter->SetTitle("Electronic Resolution Fit Parameters");
    G_Electronic_Resolution_FitParameter->GetXaxis()->SetTitle("Detector");
    G_Electronic_Resolution_FitParameter->GetYaxis()->SetTitle("Resolution (keV)");
    legend->AddEntry(G_Electronic_Resolution_FitParameter, "Electronic Resolution", "p");
    MG_Resolution_FitParameters->Add(G_Electronic_Resolution_FitParameter);

    G_Resolution_FitParameter->SetMarkerStyle(20);
    G_Resolution_FitParameter->SetMarkerColor(kRed);
    G_Resolution_FitParameter->SetTitle("Resolution Fit Parameters");
    G_Resolution_FitParameter->GetXaxis()->SetTitle("Detector");
    G_Resolution_FitParameter->GetYaxis()->SetTitle("Resolution (keV)");
    legend->AddEntry(G_Resolution_FitParameter, "Resolution", "p");
    MG_Resolution_FitParameters->Add(G_Resolution_FitParameter);

    MG_Resolution_FitParameters->Draw("AP");
    legend->Draw("SAME");
    c2->Write();

    WriteResolution();

    for (int i = first; i < last; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            Calibration_Function[i]->SetName(("Calibration_" + detectorName[i]).c_str());
            Calibration_Function[i]->Write();

            Calibration_Result[i]->SetName(("Calibration_R_" + detectorName[i]).c_str());
            Calibration_Result[i]->Write();
        }
    }

}
#endif