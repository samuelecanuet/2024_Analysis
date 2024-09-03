#ifndef MATCHER_HH
#define MATCHER_HH

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include <ctime>
#include <cstdlib>
#include <algorithm>

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TSpline.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TRandom.h"
#include "TKey.h"
#include "TMinuit.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "TF1Convolution.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"

#include "../../../lib/SignalDict/Signal.h"
#include "Detectors.hh"

using namespace std;
using namespace ROOT::Math;
default_random_engine generator;

/// FILE ///
string GROUPED_filename;
string GROUPED_basefilename;
TFile *MATCHED_File;
map<string, TFile *>SIMULATED_File;
map<string, TFile *>GROUPED_File;

string NUCLEUS;

map<string, TTree *[SIGNAL_MAX]>GROUPED_Tree_Detectors;
TTreeReader *Reader;
TTreeReaderArray<Signal> *Silicon;
double Channel;

// HISTOGRAMS
map<string, TH1D *[SIGNAL_MAX]>H_Exp;
map<string, TH1D *[SIGNAL_MAX]>H_Exp_Channel;
map<string, TH1D *[SIGNAL_MAX]>H_Sim;
map<string, TH1D *[SIGNAL_MAX]>H_Sim_Conv ;
TGraphErrors *G_Calibration[SIGNAL_MAX];
TGraphErrors *G_Calibration_Mean[SIGNAL_MAX];
TGraphErrors *G_Calibration_Alpha[SIGNAL_MAX];
map<string, TF1 *[100][SIGNAL_MAX]>F_Peak_Exp;
map<string, TF1 *[100][SIGNAL_MAX]>F_Peak_Sim;
TGraphErrors *G_Resolution[SIGNAL_MAX];

TGraphErrors *G_Mean[SIGNAL_MAX];
TGraphErrors *G_Chi2[SIGNAL_MAX];
// Function
TF1 *Background_function[SIGNAL_MAX];
TF1 *Threshold_function[SIGNAL_MAX];
TF1 *Alpha1_function;
TF1 *Alpha2_function;
double parameter[SIGNAL_MAX][100][3];
double parameterError[SIGNAL_MAX][100][3];
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
pair<double, double> ManualCalib[100][SIGNAL_MAX];
pair<double, double> ManualCalibLinear[SIGNAL_MAX];
vector<double> ManualCalibFitted[SIGNAL_MAX];
vector<double> ManualCalibFitted_Alpha[SIGNAL_MAX];

map<string, pair<double, double>[100][SIGNAL_MAX]> WindowsMap;
map<string, pair<int, int>> CanvasMap;
map<string, int> ScalerPeak;
map<string, pair<double, double> [100]> EnergyError;

TDirectory *dir_detector[SIGNAL_MAX];
TDirectory *dir_peak_detector[SIGNAL_MAX];
map<string, TDirectory *[SIGNAL_MAX]>dir_nuclei_detector;

void InitHistograms()
{
    for (int i = 0; i < detectorNum; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            if (NUCLEUS == "32Ar")
            {
                dir_detector[i] = MATCHED_File->mkdir(detectorName[i].c_str());
                dir_peak_detector[i] = dir_detector[i]->mkdir("Peaks");
            }

            dir_nuclei_detector[NUCLEUS][i] = dir_detector[i]->mkdir(NUCLEUS.c_str());

            H_Exp[NUCLEUS][i] = new TH1D(("H_Exp_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_Exp_" + NUCLEUS + "_" + detectorName[i]).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
            H_Exp[NUCLEUS][i]->GetXaxis()->SetTitle("Energy [keV]");
            H_Exp[NUCLEUS][i]->GetYaxis()->SetTitle("Counts");
            H_Exp[NUCLEUS][i]->GetXaxis()->CenterTitle();
            H_Exp[NUCLEUS][i]->GetYaxis()->CenterTitle();

            H_Exp_Channel[NUCLEUS][i] = new TH1D(("H_Exp_Channel_" + NUCLEUS + "_" + detectorName[i]).c_str(), ("H_Exp_Channel_" + NUCLEUS + "_" + detectorName[i]).c_str(), eSiliN, 0, eSiliMax/1000);
            H_Exp_Channel[NUCLEUS][i]->GetXaxis()->SetTitle("Channel");
            H_Exp_Channel[NUCLEUS][i]->GetYaxis()->SetTitle("Counts");
            H_Exp_Channel[NUCLEUS][i]->GetXaxis()->CenterTitle();
            H_Exp_Channel[NUCLEUS][i]->GetYaxis()->CenterTitle();

            G_Calibration[i] = new TGraphErrors();
            G_Calibration_Alpha[i] = new TGraphErrors();

            G_Resolution[i] = new TGraphErrors();
        }
    }
}

void InitAlphaPeaks()
{
    /// read txt file
    ifstream file("Config_Files/Alpha_Up.txt");
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

    ifstream file2("Config_Files/Alpha_Down.txt");
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
                    //cout << "Nuclei : " << nuclei << " Number : " << number << " Detector : " << detectorName[i] << " Energy Low : " << energy_low << " Energy High : " << energy_high << endl;
                }
            }
        }
    }
}

void InitManualCalibration()
{
    ifstream file("Config_Files/Manual_Calibration.txt");
    if (!file.is_open())
    {
        Error("Impossible to open Manual_Calibration.txt");
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
        ss >> detname >> channel6 >> channel14 >> channel29 >> energy6 >> energy14 >> energy29;
        for (int i = 0; i < detectorNum; i++)
        {
            if (IsDetectorSiliStrip(i))
            {
                if (detname == detectorName[i])
                {
                    ManualCalib[6][i] = make_pair(channel6, energy6);
                    ManualCalib[14][i] = make_pair(channel14, energy14);
                    ManualCalib[29][i] = make_pair(channel29, energy29);
                    break;
                }
            }
        }
    }
}

void InitEnergyErrors()
{
    ifstream file(("Config_Files/" + NUCLEUS + "_protons.txt").c_str());
    if (!file.is_open())
    {
        Warning("Impossible to open " + NUCLEUS + "_protons.txt");
    }

    string line;
    double error;
    double energy;
    int peak;
    int counter = 0;
    while (getline(file, line))
    {
        counter++;
        stringstream ss(line);
        ss >> peak >> energy >> error;
        EnergyError[NUCLEUS][peak] = make_pair(energy, error);
    }

    file.close();
}

void FillingSimHitograms()
{

    /////// HAVE TO BE CHNAGE TO TAKE INTO ACCOUNT THE INTERSTRIPS EVENTS //////

    CanvasMap["32Ar"] = make_pair(8, 5);
    CanvasMap["33Ar"] = make_pair(10, 5);
    CanvasMap["18N"] = make_pair(3, 3);

    ScalerPeak["32Ar"] = 14;
    ScalerPeak["33Ar"] = 21;
    ScalerPeak["18N"] = 1;

    for (auto &pair : SIMULATED_File)
    {
        for (int i = 0; i < detectorNum; i++)
        {
            if (IsDetectorSiliStrip(i))
            {
                H_Sim[pair.first][i] = (TH1D *)pair.second->Get((detectorName[i] + "_single").c_str());
                H_Sim[pair.first][i]->SetName(("H_Sim_" + detectorName[i]).c_str());
                H_Sim[pair.first][i]->SetTitle(("H_Sim_" + detectorName[i]).c_str());
                H_Sim[pair.first][i]->GetXaxis()->SetTitle("Energy [keV]");
                H_Sim[pair.first][i]->GetYaxis()->SetTitle("Counts");
                H_Sim[pair.first][i]->GetXaxis()->CenterTitle();
                H_Sim[pair.first][i]->GetYaxis()->CenterTitle();
                

                for (int j = 0; j < detectorNum; j++)
                {   
                    if ( (GetDetector(i) != GetDetector(j)) && (GetDetectorChannel(i) == GetDetectorChannel(j)) && IsDetectorSiliStrip(j))
                    {
                        if (GetDetector(i) <=4 && GetDetector(j) <=4 )
                            H_Sim[pair.first][i]->Add((TH1D *)pair.second->Get((detectorName[j] + "_single").c_str()));
                        else if (GetDetector(i) >=5 && GetDetector(j) >=5)
                            H_Sim[pair.first][i]->Add((TH1D *)pair.second->Get((detectorName[j] + "_single").c_str()));
                    }
                }

                H_Sim[pair.first][i]->Rebin(10);
            }
        }
    }
}

TH1D *RemoveBKG(TH1D *H, bool finished)
{
    double MINIMUM = 0;
    double MAXIMUM = 0;
    int N = 0;
    double coef = 1;
    double offset = 0;


    if (offset == 0 && coef == 1)
    {
        MINIMUM = eSiliMin_cal;
        MAXIMUM = eSiliMax_cal;
        N = eSiliN_cal;
    }
    else
    {
        MINIMUM = eSiliMin;
        MAXIMUM = eSiliMax;
        N = eSiliN;
    }

    TH1D *Selection = (TH1D *)H->Clone("test");
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
        Mean_A1_low =  (WindowsMap["18N"][1][current_detector].first-offset)/coef;
        Mean_A1_high = (WindowsMap["18N"][1][current_detector].second-offset)/coef;
        Mean_A2_low =  (WindowsMap["18N"][2][current_detector].first-offset)/coef;
        Mean_A2_high = (WindowsMap["18N"][2][current_detector].second-offset)/coef;
    }
    else
    {
        Mean_A1_low =  (WindowsMap["18N"][1][current_detector].first-offset)/coef;
        Mean_A1_high = (WindowsMap["18N"][1][current_detector].second-offset)/coef;
        Mean_A2_low =  (WindowsMap["18N"][2][current_detector].first-offset)/coef;
        Mean_A2_high = (WindowsMap["18N"][2][current_detector].second-offset)/coef;
    }

    BKG_function->SetParameters(1000, (Mean_A1_low+Mean_A1_high)/2, 25/coef, 500, (Mean_A2_low+Mean_A2_high)/2, 25/coef, 1000, -3e-3*coef, 500, 25);
    BKG_function->SetParLimits(0, 5, 30000);
    BKG_function->SetParLimits(1, Mean_A1_low, Mean_A1_high);
    BKG_function->SetParLimits(2, 5/coef, 50/coef);
    BKG_function->SetParLimits(3, 5, 30000);
    BKG_function->SetParLimits(4, Mean_A2_low, Mean_A2_high);
    BKG_function->SetParLimits(5, 5/coef, 50/coef);
    BKG_function->SetParLimits(6, 5, 15000);
    BKG_function->SetParLimits(7, -1e-2*coef, 0);
    BKG_function->SetParLimits(8, (400-offset)/coef, (Mean_A1_low+Mean_A1_high)/2);
    BKG_function->SetParLimits(9, 5/coef, 100/coef);

    Selection->Fit("BKG_function", "QRN", "", (400-offset)/coef, (1900-offset)/coef);


    alpha1 = BKG_function->GetParameter(1)*coef+offset;
    alpha2 = BKG_function->GetParameter(4)*coef+offset;

    Alpha1_function = new TF1("Alpha1_function", "gaus(0)", MINIMUM, MAXIMUM);
    Alpha1_function->SetParameters(BKG_function->GetParameter(0), BKG_function->GetParameter(1), BKG_function->GetParameter(2));

    Alpha2_function = new TF1("Alpha2_function", "gaus(0)", MINIMUM, MAXIMUM);
    Alpha2_function->SetParameters(BKG_function->GetParameter(3), BKG_function->GetParameter(4), BKG_function->GetParameter(5));

    Background_function[current_detector] = new TF1(("BKG_function_sub_" + detectorName[current_detector]).c_str(), "[0]*exp([1] * x) * 0.5 * erfc((x-[2])/[3])", MINIMUM, MAXIMUM);
    Background_function[current_detector]->SetNpx(N);
    double mean_erfc = ManualCalibFitted[current_detector][0]+ ManualCalibFitted[current_detector][1] * 24.1198*1.02 + ManualCalibFitted[current_detector][2] * 24.1198 * 24.1198*1.02*1.02;
    double sigma_erfc =  ManualCalibFitted[current_detector][0]+ ManualCalibFitted[current_detector][1] * 3.821 + ManualCalibFitted[current_detector][2] * 3.821 * 3.821;
    Background_function[current_detector]->SetParameters(BKG_function->GetParameter(6), BKG_function->GetParameter(7), mean_erfc, sigma_erfc);
    Threshold_function[current_detector] = new TF1(("Threshold_function_" + detectorName[current_detector]).c_str(), "0.5*(1+TMath::Erf((x-[0])/[1]))", MINIMUM, MAXIMUM);
    Threshold_function[current_detector]->SetNpx(N);
    Threshold_function[current_detector]->SetParameters(BKG_function->GetParameter(8), BKG_function->GetParameter(9));

    /// PLOT//////////////////////////////
    if (finished)
    {
        TCanvas *c = new TCanvas((detectorName[current_detector] + "BKG").c_str(), (detectorName[current_detector] + "BKG").c_str(), 1920, 1080);
        c->Divide(1, 2);
        c->cd(1);
        Selection->GetXaxis()->SetRangeUser((400-offset)/coef, (1900-offset)/coef);
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
        h_sub->GetXaxis()->SetRangeUser((400-offset)/coef, (1900-offset)/coef);
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
    /////////////////////////////////////
}

double FunctionToMinimize(const double *par)
{
    /// 32Ar
    NUCLEUS = "32Ar";
    H_Sim_Conv[NUCLEUS][current_detector] = (TH1D *)H_Sim["32Ar"][current_detector]->Clone();
    H_Sim_Conv[NUCLEUS][current_detector]->Reset();
    Reader = new TTreeReader(GROUPED_Tree_Detectors[NUCLEUS][current_detector]);
    Reader->Restart();

    H_Exp[NUCLEUS][current_detector]->Reset();
    TTreeReaderValue<double> ChannelDet(*Reader, "Channel");
    while (Reader->Next())
    {
        double value = *ChannelDet / 1000;
        H_Exp[NUCLEUS][current_detector]->Fill(ManualCalibFitted[current_detector][0]+ ManualCalibFitted[current_detector][1] * value + ManualCalibFitted[current_detector][2] * value * value);
    }

    H_Exp[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(-1111, -1111);
    TH1D* H_Exp_Without_Background = RemoveBKG(H_Exp[NUCLEUS][current_detector], plotting);

    
    TH1D * H = (TH1D *)H_Sim["32Ar"][current_detector]->Clone();
    H->Reset();
    for (int i = 0; i < H_Sim["32Ar"][current_detector]->GetEntries(); i++)
    {
        double value = H_Sim["32Ar"][current_detector]->GetRandom();
        normal_distribution<double> distribution(value, par[0] + par[1]* sqrt(value));
        H->Fill(distribution(generator));
    }

    H_Exp[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(WindowsMap["32Ar"][ScalerPeak["32Ar"]][current_detector].first, WindowsMap["32Ar"][ScalerPeak["32Ar"]][current_detector].second);
    H->GetXaxis()->SetRangeUser(WindowsMap["32Ar"][ScalerPeak["32Ar"]][current_detector].first, WindowsMap["32Ar"][ScalerPeak["32Ar"]][current_detector].second);
    H->Scale(H_Exp[NUCLEUS][current_detector]->GetMaximum() / H->GetMaximum());
    H_Sim_Conv[NUCLEUS][current_detector]->Add(H, 1);
    //Add bkg
    TH1D* h = (TH1D*)Background_function[current_detector]->GetHistogram();
    H_Sim_Conv[NUCLEUS][current_detector]->Add(h, 1);

    /// 18N
    if (SIMULATED_File["18N"])
    {
        double alpha_res = 0.0;
        TH1D * H1 = (TH1D *)H_Sim["18N"][current_detector]->Clone();
        H1->Reset();
        for (int i = 0; i < H_Sim["18N"][current_detector]->GetEntries(); i++)
        {
            double value_corrected = (H_Sim["18N"][current_detector]->GetRandom() - ManualCalibFitted_Alpha[current_detector][0]) / ManualCalibFitted_Alpha[current_detector][1];
            double value = ManualCalibFitted[current_detector][0] + ManualCalibFitted[current_detector][1] * value_corrected + ManualCalibFitted[current_detector][2] * value_corrected * value_corrected;
            normal_distribution<double> distribution(value, alpha_res * sqrt(value));
            H1->Fill(distribution(generator));
        }

        // FIRST PEAK
        TH1D* H11 = (TH1D*)H1->Clone();
        H_Exp_Without_Background->GetXaxis()->SetRangeUser(WindowsMap["18N"][1][current_detector].first, WindowsMap["18N"][1][current_detector].second);
        H11->GetXaxis()->SetRangeUser(WindowsMap["18N"][1][current_detector].first, WindowsMap["18N"][1][current_detector].second);
        H11->Scale(H_Exp_Without_Background->GetMaximum() / H11->GetMaximum());
        for (int bin = 0; bin < H11->GetNbinsX(); bin++)
        {
            if (H11->GetBinCenter(bin) < WindowsMap["18N"][1][current_detector].first || H11->GetBinCenter(bin) > WindowsMap["18N"][1][current_detector].second)
            {
                H11->SetBinContent(bin, 0);
            }
        }
        
        // SECOND PEAK
        TH1D* H12 = (TH1D*)H1->Clone();
        H_Exp_Without_Background->GetXaxis()->SetRangeUser(WindowsMap["18N"][2][current_detector].first, WindowsMap["18N"][2][current_detector].second);
        H12->GetXaxis()->SetRangeUser(WindowsMap["18N"][2][current_detector].first, WindowsMap["18N"][2][current_detector].second);
        H12->Scale(H_Exp_Without_Background->GetMaximum() / H12->GetMaximum());
        for (int bin = 0; bin < H12->GetNbinsX(); bin++)
        {
            if (H12->GetBinCenter(bin) < WindowsMap["18N"][2][current_detector].first || H12->GetBinCenter(bin) > WindowsMap["18N"][2][current_detector].second)
            {
                H12->SetBinContent(bin, 0);
            }
        }

    
        // Other peaks
        H1->Reset();
        TH1D *sub = (TH1D*)H_Exp[NUCLEUS][current_detector]->Clone();
        sub->Add(H_Sim_Conv[NUCLEUS][current_detector], -1);
        sub->GetXaxis()->SetRangeUser(-1111, -1111);
        TF1* gauss = new TF1("gaus", "gaus", 1600, 2300);
        gauss->SetParameters(1000, 2000, 100);
        sub->Fit("gaus", "QRN", "", 1600, 2300);
        double first = gauss->GetParameter(1);
        for (int i = 0; i < H_Sim["18N"][current_detector]->GetEntries(); i++)
        {
            double value = H_Sim["18N"][current_detector]->GetRandom();
            normal_distribution<double> distribution(value, alpha_res * sqrt(value));
            H1->Fill(distribution(generator));
        }
        H1->GetXaxis()->SetRangeUser(-1111, -1111);
        H1->Fit("gaus", "QRN", "", 1600, 2300);
        double offset = gauss->GetParameter(1) - first;
        H1->Reset();
        for (int i = 0; i < H_Sim["18N"][current_detector]->GetEntries(); i++)
        {
            double value = H_Sim["18N"][current_detector]->GetRandom();
            normal_distribution<double> distribution(value, alpha_res * sqrt(value));
            H1->Fill(distribution(generator));
        }

        
        TH1D* H13 = (TH1D*)H1->Clone();
        H13->Scale(0.6);
        for (int bin = 0; bin < H13->GetNbinsX(); bin++)
        {
            if (H13->GetBinCenter(bin) < WindowsMap["18N"][2][current_detector].second)
            {
                H13->SetBinContent(bin, 0);
            }
        }
        //ADDING 
        H_Sim_Conv[NUCLEUS][current_detector]->Add(H11, 1);
        H_Sim_Conv[NUCLEUS][current_detector]->Add(H12, 1);
        H_Sim_Conv[NUCLEUS][current_detector]->Add(H13, 1);
    }

    

    /// apply thrshold  
    for (int bin = 0; bin < H_Exp[NUCLEUS][current_detector]->GetNbinsX(); bin++)
    {
        H_Sim_Conv[NUCLEUS][current_detector]->SetBinContent(bin, H_Sim_Conv[NUCLEUS][current_detector]->GetBinContent(bin) * Threshold_function[current_detector]->Eval(H_Sim_Conv[NUCLEUS][current_detector]->GetBinCenter(bin)));
    }

    /// 33Ar
    NUCLEUS = "33Ar";
    H_Sim_Conv[NUCLEUS][current_detector] = (TH1D *)H_Sim["33Ar"][current_detector]->Clone();
    H_Sim_Conv[NUCLEUS][current_detector]->Reset();
    Reader = new TTreeReader(GROUPED_Tree_Detectors[NUCLEUS][current_detector]);
    Reader->Restart();
    TTreeReaderValue<double> ChannelDet33(*Reader, "Channel");
    H_Exp["33Ar"][current_detector]->Reset();
    while (Reader->Next())
    {
        double value = *ChannelDet33 / 1000;
        H_Exp["33Ar"][current_detector]->Fill(ManualCalibFitted[current_detector][0] + ManualCalibFitted[current_detector][1] * value + ManualCalibFitted[current_detector][2] * value * value);
    }
    
    TH1D* H33Ar = (TH1D *)H_Sim["33Ar"][current_detector]->Clone();
    H33Ar->Reset();
    for (int i = 0; i < H_Sim["33Ar"][current_detector]->GetEntries(); i++)
    {
        double value = H_Sim["33Ar"][current_detector]->GetRandom();
        normal_distribution<double> distribution(value, par[0] + par[1]* sqrt(value));
        H33Ar->Fill(distribution(generator));
    }

    H_Exp[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(WindowsMap["33Ar"][ScalerPeak["33Ar"]][current_detector].first, WindowsMap["33Ar"][ScalerPeak["33Ar"]][current_detector].second);
    H33Ar->GetXaxis()->SetRangeUser(WindowsMap["33Ar"][ScalerPeak["33Ar"]][current_detector].first, WindowsMap["33Ar"][ScalerPeak["33Ar"]][current_detector].second);
    H33Ar->Scale(H_Exp[NUCLEUS][current_detector]->GetMaximum() / H->GetMaximum());
    H_Sim_Conv[NUCLEUS][current_detector]->Add(H33Ar, 1);

    


    // H_Sim_Conv[NUCLEUS][current_detector]->Write();

    // H_Exp[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(1900, 6000);
    // H_Sim[current_detector]->GetXaxis()->SetRangeUser(1900, 6000);
    // double offset = H_Exp[NUCLEUS][current_detector]->GetBinCenter(H_Exp[NUCLEUS][current_detector]->GetMaximumBin()) - H_Sim[current_detector]->GetBinCenter(H_Sim[current_detector]->GetMaximumBin());
    // H_Exp[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(-1111, -1111);
    // H_Sim_Conv[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(-1111, -1111);

    // Reader->Restart();
    // H_Exp[NUCLEUS][current_detector]->Reset();
    // while (Reader->Next())
    // {
    //     H_Exp[NUCLEUS][current_detector]->Fill(-offset + *ChannelDet / 1000 * par[0]);
    // }

    // H_Exp[NUCLEUS][current_detector] = RemoveBKG(H_Exp[NUCLEUS][current_detector], par, plotting);

    // for (int bin = 0; bin < H_Exp[NUCLEUS][current_detector]->GetNbinsX(); bin++)
    // {
    //     if (H_Exp[NUCLEUS][current_detector]->GetBinContent(bin) <= 1)
    //     {
    //         H_Exp[NUCLEUS][current_detector]->SetBinContent(bin, 1);
    //     }

    //     if (H_Sim_Conv[NUCLEUS][current_detector]->GetBinContent(bin) == 0)
    //     {
    //         H_Sim_Conv[NUCLEUS][current_detector]->SetBinContent(bin, 1);
    //     }
    // }

    // H_Exp[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(1900, 6000);
    // H_Sim_Conv[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(1900, 6000);
    // double chi2 = H_Exp[NUCLEUS][current_detector]->Chi2Test(H_Sim_Conv[NUCLEUS][current_detector], "CHI2/NDF");
    // H_Exp[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(-1111, -1111);
    // H_Sim_Conv[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(-1111, -1111);

    NUCLEUS = "32Ar";
    H_Exp[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(400, 6500);
    H_Sim_Conv[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(400, 6500);
    double chi2 = H_Exp[NUCLEUS][current_detector]->Chi2Test(H_Sim_Conv[NUCLEUS][current_detector], "UW CHI2/NDF");
    H_Exp[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(-1111, -1111);
    H_Sim_Conv[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(-1111, -1111);

    cout << "CHI2 : " << chi2;
    cout << setprecision(5) << "    PAR : " << par[0] << " " << par[1];

    if (CHI2 > chi2 && chi2 > 0)
    {
        CHI2 = chi2;
        cout << "   BEST" << endl;
        parameter[current_detector][current_peak][0] = par[0];
        parameter[current_detector][current_peak][1] = par[1];
    }
    else
    {
        cout << endl;
    }

    dir_nuclei_detector["32Ar"][current_detector]->cd();
    H_Exp["32Ar"][current_detector]->Write();
    H_Sim_Conv["32Ar"][current_detector]->Write();
    // dir_nuclei_detector["33Ar"][current_detector]->cd();
    // H_Exp["33Ar"][current_detector]->Write();
    // H_Sim_Conv["33Ar"][current_detector]->Write();

    return abs(chi2);
}

void CHI2Minimization()
{
    //////// GUESS FROM SHOW PEAKS //////
    // INCLUDE ALPHAS
    // REMOVING BKG

    if (ManualCalibFitted[current_detector].size() == 0)
    {
        return;
    }

    CHI2 = 1000000;
    plotting = false;
    Minimizer *minimizer = Factory::CreateMinimizer("Minuit2", "Migrad");
    ROOT::Math::Functor functor(&FunctionToMinimize, 2);
    minimizer->SetFunction(functor);
    minimizer->SetStrategy(2);
    minimizer->SetLimitedVariable(0, "offset-resolution", 5, 5, 0, 10);
    minimizer->SetLimitedVariable(1, "resolution", 0.06, 0.001, 0., 0.1);
    minimizer->SetPrecision(0.001);
    minimizer->SetTolerance(0.001);
    minimizer->SetMaxFunctionCalls(1000000);
    minimizer->SetMaxIterations(1000000);

    // minimizer->Minimize();
    // const double *par = minimizer->X();

    plotting = true;  
    const double par2[2] = {ManualCalibFitted[current_detector][1] * (30/1000) + ManualCalibFitted[current_detector][2] * (30/1000) * (30/1000), sqrt(3.6*1e-3*0.134)};
    // const double par2[2] = {parameter[current_detector][current_peak][0], parameter[current_detector][current_peak][1]};
    FunctionToMinimize(par2);
    plotting = false;

    NUCLEUS="32Ar";
    TCanvas *c1 = new TCanvas((detectorName[current_detector] + "_" + NUCLEUS).c_str(), (detectorName[current_detector] + "_" + NUCLEUS).c_str(), 1920, 1080);
    c1->cd();
    c1->SetLogy();
    H_Exp[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(WindowsMap["32Ar"][14][current_detector].first, WindowsMap["32Ar"][14][current_detector].second);
    H_Sim_Conv[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(WindowsMap["32Ar"][14][current_detector].first, WindowsMap["32Ar"][14][current_detector].second);
    H_Sim_Conv[NUCLEUS][current_detector]->Scale(H_Exp[NUCLEUS][current_detector]->GetMaximum() / H_Sim_Conv[NUCLEUS][current_detector]->GetMaximum());
    H_Exp[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(0, 6500);
    H_Sim_Conv[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(0, 6500);

    H_Exp[NUCLEUS][current_detector]->Draw("HIST");
    H_Sim_Conv[NUCLEUS][current_detector]->SetLineColor(kRed);
    H_Sim_Conv[NUCLEUS][current_detector]->Draw("HIST SAME");

    dir_detector[current_detector]->cd();
    c1->Write();

    NUCLEUS = "33Ar";
    TCanvas *c2 = new TCanvas((detectorName[current_detector] + "_" + NUCLEUS).c_str(), (detectorName[current_detector] + "_" + NUCLEUS).c_str(), 1920, 1080);
    c2->cd();
    c2->SetLogy();
    H_Exp[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(WindowsMap["33Ar"][21][current_detector].first, WindowsMap["33Ar"][21][current_detector].second);
    H_Sim_Conv[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(WindowsMap["33Ar"][21][current_detector].first, WindowsMap["33Ar"][21][current_detector].second);
    H_Sim_Conv[NUCLEUS][current_detector]->Scale(H_Exp[NUCLEUS][current_detector]->GetMaximum() / H_Sim_Conv[NUCLEUS][current_detector]->GetMaximum());
    H_Exp[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(0, 6500);
    H_Sim_Conv[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(0, 6500);

    H_Exp[NUCLEUS][current_detector]->Draw("HIST");
    H_Sim_Conv[NUCLEUS][current_detector]->SetLineColor(kRed);
    H_Sim_Conv[NUCLEUS][current_detector]->Draw("HIST SAME");

    dir_detector[current_detector]->cd();
    c2->Write();
    NUCLEUS = "32Ar";

    dir_nuclei_detector["32Ar"][current_detector]->cd();
    H_Exp["32Ar"][current_detector]->Write();
    H_Sim["32Ar"][current_detector]->Write();
    H_Sim["18N"][current_detector]->Write();
    dir_nuclei_detector["33Ar"][current_detector]->cd();
    H_Sim["33Ar"][current_detector]->Write();
    H_Exp["33Ar"][current_detector]->Write();

}


void Manual_Calibration()
{

    Reader->Restart();
    TTreeReaderValue<double> ChannelDet(*Reader, "Channel");
    while (Reader->Next())
    {
        H_Exp_Channel[NUCLEUS][current_detector]->Fill(*ChannelDet / 1000);
    }
    dir_nuclei_detector[NUCLEUS][current_detector]->cd();
    H_Exp_Channel[NUCLEUS][current_detector]->Write();

    int counter = 0;
    for (int peak = 1; peak < CanvasMap["32Ar"].first * CanvasMap["32Ar"].second; peak++)
    {
        if (ManualCalib[peak][current_detector].first == -1 || !ManualCalib[peak][current_detector].first)
            continue;   

        G_Calibration[current_detector]->SetPoint(counter, ManualCalib[peak][current_detector].first, ManualCalib[peak][current_detector].second);
        counter ++;
    }

    if (G_Calibration[current_detector]->GetN() == 0)
    {
        return;
    }

    //FIRST FIT WITH MANUAL
    G_Calibration[current_detector]->Fit("pol1", "Q");
    ManualCalibLinear[current_detector] = make_pair(G_Calibration[current_detector]->GetFunction("pol1")->GetParameter(0), G_Calibration[current_detector]->GetFunction("pol1")->GetParameter(1));
    Info( detectorName[current_detector] + "   First Calibration : a = " + to_string(G_Calibration[current_detector]->GetFunction("pol1")->GetParameter(0)) + " b = " + to_string(G_Calibration[current_detector]->GetFunction("pol1")->GetParameter(1)));
}

void Fitting_Calibration()
{
    //SECOND FIT WITH GAUSSIANS IN WINDOWS
    double coef = ManualCalibLinear[current_detector].second;
    double offset = ManualCalibLinear[current_detector].first;
    G_Calibration[current_detector]->Set(0);
    int counter = 0;

    NUCLEUS = "32Ar";
    for (int peak = 1; peak < CanvasMap["32Ar"].first * CanvasMap["32Ar"].second; peak++)
    {

        if (WindowsMap["32Ar"][peak][current_detector].first == -1 || !WindowsMap["32Ar"][peak][current_detector].first)
            continue;

        if (peak != 5 && peak != 8 && peak != 9 && peak != 14 && peak != 23 && peak != 25 && peak != 28 && peak != 29 && peak != 30)
        {
            if (peak != 3 && GetDetector(current_detector) <= 4)
            {
                continue;
            }
            if (peak != 3)
            {
                continue;
            }
            continue;
        }

        if (peak != 14 && peak != 8 && peak != 23)
        {
            continue;
        }

        // taking into account the exponential bkg
        
        if (WindowsMap["32Ar"][peak][current_detector].first < 1900)
        {
            F_Peak_Exp[NUCLEUS][peak][current_detector] = new TF1((detectorName[current_detector]  +"_Exp" + to_string(peak) + NUCLEUS).c_str(), "gaus(0) + [3]*exp([4]*x)", 0, 6500);
            F_Peak_Exp[NUCLEUS][peak][current_detector]->SetParLimits(0, 0, 10000);
            F_Peak_Exp[NUCLEUS][peak][current_detector]->SetParameter(1, (WindowsMap["32Ar"][peak][current_detector].second - WindowsMap["32Ar"][peak][current_detector].first) / 2 / coef);
            F_Peak_Exp[NUCLEUS][peak][current_detector]->SetParLimits(1, (WindowsMap["32Ar"][peak][current_detector].first - offset) / coef, (WindowsMap["32Ar"][peak][current_detector].second - offset) / coef);
            F_Peak_Exp[NUCLEUS][peak][current_detector]->SetParLimits(2, 0, 100);
            F_Peak_Exp[NUCLEUS][peak][current_detector]->SetParLimits(3, 0, 10000);
            F_Peak_Exp[NUCLEUS][peak][current_detector]->SetParLimits(4, -1, 0);
            H_Exp_Channel[NUCLEUS][current_detector]->Fit((detectorName[current_detector]  +"_Exp" + to_string(peak) + NUCLEUS).c_str(), "QRN", "", (WindowsMap["32Ar"][peak][current_detector].first - offset) / coef, (WindowsMap["32Ar"][peak][current_detector].second - offset) / coef);
        }
        else
        {
            F_Peak_Exp[NUCLEUS][peak][current_detector] = new TF1((detectorName[current_detector]  +"_Exp" + to_string(peak) + NUCLEUS).c_str(), "gaus", 0, 6500);
            F_Peak_Exp[NUCLEUS][peak][current_detector]->SetParLimits(0, 0, 10000);
            F_Peak_Exp[NUCLEUS][peak][current_detector]->SetParameter(1, (WindowsMap["32Ar"][peak][current_detector].second - WindowsMap["32Ar"][peak][current_detector].first) / 2 / coef);
            F_Peak_Exp[NUCLEUS][peak][current_detector]->SetParLimits(1, (WindowsMap["32Ar"][peak][current_detector].first - offset) / coef, (WindowsMap["32Ar"][peak][current_detector].second - offset) / coef);
            F_Peak_Exp[NUCLEUS][peak][current_detector]->SetParLimits(2, 0, 100);
            H_Exp_Channel[NUCLEUS][current_detector]->Fit((detectorName[current_detector]  +"_Exp" + to_string(peak) + NUCLEUS).c_str(), "QRN", "", (WindowsMap["32Ar"][peak][current_detector].first - offset) / coef, (WindowsMap["32Ar"][peak][current_detector].second - offset) / coef);
        }


        TCanvas *c2 = new TCanvas((detectorName[current_detector] + "_Calibration_" + to_string(peak)).c_str(), (detectorName[current_detector] + "_Calibration_" + to_string(peak)).c_str(), 1920, 1080);
        c2->Divide(2, 1);
        c2->cd(1);
        TH1D *H = (TH1D *)H_Exp_Channel[NUCLEUS][current_detector]->Clone();
        H->GetXaxis()->SetRangeUser((WindowsMap["32Ar"][peak][current_detector].first - offset) / coef, (WindowsMap["32Ar"][peak][current_detector].second - offset) / coef);
        H->Draw("HIST");
        F_Peak_Exp[NUCLEUS][peak][current_detector]->SetNpx(eSiliN_cal*10);
        F_Peak_Exp[NUCLEUS][peak][current_detector]->SetLineColor(kRed);
        F_Peak_Exp[NUCLEUS][peak][current_detector]->Draw("SAME");


        c2->cd(2);
        F_Peak_Sim[NUCLEUS][peak][current_detector] = new TF1((detectorName[current_detector]  +"_Sim" + to_string(peak) + NUCLEUS).c_str(), "gaus(0)", 0, 6500);
        F_Peak_Sim[NUCLEUS][peak][current_detector]->SetParameter(1, (WindowsMap["32Ar"][peak][current_detector].second-WindowsMap["32Ar"][peak][current_detector].first)/2);
        H_Sim["32Ar"][current_detector]->Fit((detectorName[current_detector]  +"_Sim" + to_string(peak) + NUCLEUS).c_str(), "QRN", "", WindowsMap["32Ar"][peak][current_detector].first, WindowsMap["32Ar"][peak][current_detector].second);
        TH1D* H1 = (TH1D*)H_Sim["32Ar"][current_detector]->Clone();
        H1->GetXaxis()->SetRangeUser(WindowsMap["32Ar"][peak][current_detector].first, WindowsMap["32Ar"][peak][current_detector].second);
        H1->Draw("HIST");
        F_Peak_Sim[NUCLEUS][peak][current_detector]->SetNpx(eSiliN_cal*10);
        F_Peak_Sim[NUCLEUS][peak][current_detector]->SetLineColor(kRed);
        F_Peak_Sim[NUCLEUS][peak][current_detector]->Draw("SAME L ");
        
        dir_peak_detector[current_detector]->cd();
        c2->Write();

        

        if (peak == 14)
        {
            H_Exp_Channel[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser((WindowsMap["32Ar"][peak][current_detector].first-offset)/coef, (WindowsMap["32Ar"][peak][current_detector].second-offset)/coef);
            H_Sim[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(WindowsMap["32Ar"][peak][current_detector].first, WindowsMap["32Ar"][peak][current_detector].second);
            G_Calibration[current_detector]->SetPoint(counter, H_Exp_Channel[NUCLEUS][current_detector]->GetMean(),H_Sim[NUCLEUS][current_detector]->GetMean());   
            G_Calibration[current_detector]->SetPointError(counter, H_Exp_Channel[NUCLEUS][current_detector]->GetMeanError(), sqrt(pow(H_Sim[NUCLEUS][current_detector]->GetMeanError(),2)+ pow(EnergyError["32Ar"][peak].second, 2)));
            H_Exp_Channel[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(-1111, -1111);
            H_Sim[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(-1111, -1111);
            counter++;
        }
        else
        {
            G_Calibration[current_detector]->SetPoint(counter, F_Peak_Exp[NUCLEUS][peak][current_detector]->GetParameter(1), F_Peak_Sim[NUCLEUS][peak][current_detector]->GetParameter(1));
            G_Calibration[current_detector]->SetPointError(counter, F_Peak_Exp[NUCLEUS][peak][current_detector]->GetParError(1), sqrt(pow(F_Peak_Sim[NUCLEUS][peak][current_detector]->GetParError(1), 2) + pow(EnergyError["32Ar"][peak].second, 2)));
            counter++;
        }

    }

    NUCLEUS = "33Ar";
    Reader = new TTreeReader(GROUPED_Tree_Detectors[NUCLEUS][current_detector]);
    Reader->Restart();
    TTreeReaderValue<double> ChannelDet33(*Reader, "Channel");
    while (Reader->Next())
    {
        H_Exp_Channel[NUCLEUS][current_detector]->Fill(*ChannelDet33 / 1000);
    }
    dir_nuclei_detector[NUCLEUS][current_detector]->cd();
    H_Exp_Channel[NUCLEUS][current_detector]->Write();

    for (int peak = 1; peak < CanvasMap["33Ar"].first * CanvasMap["33Ar"].second; peak++)
    {
        if (WindowsMap["33Ar"][peak][current_detector].first == -1 || !WindowsMap["33Ar"][peak][current_detector].first)
            continue;
        if (peak != 12 && peak != 21 && peak != 26 && peak != 35 && peak != 39 && peak != 40)
            continue;
        if (peak != 21 && peak != 40)
            continue;

        // taking into account the exponential bkg       
        F_Peak_Exp[NUCLEUS][peak][current_detector] = new TF1((detectorName[current_detector]  +"_Exp" + to_string(peak) + NUCLEUS).c_str(), "gaus", 0, 6500);
        F_Peak_Exp[NUCLEUS][peak][current_detector]->SetParLimits(0, 0, 10000);
        F_Peak_Exp[NUCLEUS][peak][current_detector]->SetParameter(1, (WindowsMap["33Ar"][peak][current_detector].second - WindowsMap["33Ar"][peak][current_detector].first) / 2 / coef);
        F_Peak_Exp[NUCLEUS][peak][current_detector]->SetParLimits(1, (WindowsMap["33Ar"][peak][current_detector].first - offset) / coef, (WindowsMap["33Ar"][peak][current_detector].second - offset) / coef);
        F_Peak_Exp[NUCLEUS][peak][current_detector]->SetParLimits(2, 0, 100);
        H_Exp_Channel[NUCLEUS][current_detector]->Fit((detectorName[current_detector]  +"_Exp" + to_string(peak) + NUCLEUS).c_str(), "QRN", "", (WindowsMap["33Ar"][peak][current_detector].first - offset) / coef, (WindowsMap["33Ar"][peak][current_detector].second - offset) / coef);
        
        TCanvas *c2 = new TCanvas((detectorName[current_detector] + "_Calibration_" + to_string(peak)).c_str(), (detectorName[current_detector] + "_Calibration_" + to_string(peak)).c_str(), 1920, 1080);
        c2->Divide(2, 1);
        c2->cd(1);
        TH1D *H = (TH1D *)H_Exp_Channel[NUCLEUS][current_detector]->Clone();
        H->GetXaxis()->SetRangeUser((WindowsMap["33Ar"][peak][current_detector].first - offset) / coef, (WindowsMap["33Ar"][peak][current_detector].second - offset) / coef);
        H->Draw("HIST");
        F_Peak_Exp[NUCLEUS][peak][current_detector]->SetNpx(eSiliN_cal*10);
        F_Peak_Exp[NUCLEUS][peak][current_detector]->SetLineColor(kRed);
        F_Peak_Exp[NUCLEUS][peak][current_detector]->Draw("SAME");

        c2->cd(2);
        F_Peak_Sim[NUCLEUS][peak][current_detector] = new TF1((detectorName[current_detector]  +"_Sim" + to_string(peak) + NUCLEUS).c_str(), "gaus(0)", 0, 6500);
        F_Peak_Sim[NUCLEUS][peak][current_detector]->SetParameter(1, (WindowsMap["33Ar"][peak][current_detector].second-WindowsMap["33Ar"][peak][current_detector].first)/2);
        H_Sim["33Ar"][current_detector]->Fit((detectorName[current_detector]  +"_Sim" + to_string(peak) + NUCLEUS).c_str(), "QRN", "", WindowsMap["33Ar"][peak][current_detector].first, WindowsMap["33Ar"][peak][current_detector].second);
        TH1D* H1 = (TH1D*)H_Sim["33Ar"][current_detector]->Clone();
        H1->GetXaxis()->SetRangeUser(WindowsMap["33Ar"][peak][current_detector].first, WindowsMap["33Ar"][peak][current_detector].second);
        H1->Draw("HIST");
        F_Peak_Sim[NUCLEUS][peak][current_detector]->SetNpx(eSiliN_cal*10);
        F_Peak_Sim[NUCLEUS][peak][current_detector]->SetLineColor(kRed);
        F_Peak_Sim[NUCLEUS][peak][current_detector]->Draw("SAME");
        
        dir_peak_detector[current_detector]->cd();
        c2->Write();

        if (peak == 21)
        {
            H_Exp_Channel[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser((WindowsMap[NUCLEUS][peak][current_detector].first-offset)/coef, (WindowsMap[NUCLEUS][peak][current_detector].second-offset)/coef);
            H_Sim[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][peak][current_detector].first, WindowsMap[NUCLEUS][peak][current_detector].second);
            G_Calibration[current_detector]->SetPoint(counter, H_Exp_Channel[NUCLEUS][current_detector]->GetMean(),H_Sim[NUCLEUS][current_detector]->GetMean());   
            G_Calibration[current_detector]->SetPointError(counter, H_Exp_Channel[NUCLEUS][current_detector]->GetMeanError(), sqrt(pow(H_Sim[NUCLEUS][current_detector]->GetMeanError(),2)+ pow(EnergyError[NUCLEUS][peak].second, 2)));
            H_Exp_Channel[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(-1111, -1111);
            H_Sim[NUCLEUS][current_detector]->GetXaxis()->SetRangeUser(-1111, -1111);
            counter++;
        }
        else
        {
            G_Calibration[current_detector]->SetPoint(counter, F_Peak_Exp[NUCLEUS][peak][current_detector]->GetParameter(1), F_Peak_Sim[NUCLEUS][peak][current_detector]->GetParameter(1));
            G_Calibration[current_detector]->SetPointError(counter, F_Peak_Exp[NUCLEUS][peak][current_detector]->GetParError(1), sqrt(pow(F_Peak_Sim[NUCLEUS][peak][current_detector]->GetParError(1), 2) + pow(EnergyError[NUCLEUS][peak].second, 2)));
            counter++;
        }

    }
    
    TCanvas *c1 = new TCanvas((detectorName[current_detector] + "_Calibration").c_str(), (detectorName[current_detector] + "_Calibration").c_str(), 1920, 1080);
    c1->Divide(1, 2);
    c1->cd(1);
    string function = "pol2";
    TF1* fit = new TF1(function.c_str(), "[0] + [1] * x + [2] * x * x", 0, 6500);
    fit->SetParLimits(0, -100, 0);
    // fit->FixParameter(2, 0);
    G_Calibration[current_detector]->Fit(function.c_str(), "Q");
    Info( detectorName[current_detector] + "   Second Calibrarion : a = " + to_string(G_Calibration[current_detector]->GetFunction(function.c_str())->GetParameter(0)) + " b = " + to_string(G_Calibration[current_detector]->GetFunction(function.c_str())->GetParameter(1)) + " c = " + to_string(G_Calibration[current_detector]->GetFunction(function.c_str())->GetParameter(2)));
    G_Calibration[current_detector]->GetFunction(function.c_str())->SetLineColor(kRed);
    G_Calibration[current_detector]->Draw("*AP");
    G_Calibration[current_detector]->GetXaxis()->SetTitle("Channel");
    G_Calibration[current_detector]->GetYaxis()->SetTitle("Energy [keV]");
    TF1* linear = new TF1("linear", "[0] + [1] * x", 0, 6500);
    TGraphErrors *g = (TGraphErrors*)G_Calibration[current_detector]->Clone();
    g->Fit("linear", "Q");
    linear->SetLineColor(kGreen);
    linear->Draw("SAME");
    c1->cd(2);
    // draw diff between fet and data
    TGraphErrors *G_Calibration_diff = (TGraphErrors*)G_Calibration[current_detector]->Clone();
    for (int i = 0; i < G_Calibration_diff->GetN(); i++)
    {
        double x, y;
        G_Calibration_diff->GetPoint(i, x, y);
        G_Calibration_diff->SetPoint(i, x, y - G_Calibration[current_detector]->GetFunction(function.c_str())->Eval(x));
    }
    G_Calibration_diff->Draw("*AP");
    G_Calibration_diff->GetXaxis()->SetTitle("Channel");
    G_Calibration_diff->GetYaxis()->SetTitle("Energy [keV]");
    TLine *line = new TLine(0, 0, 100, 0);
    line->SetLineColor(kRed);
    line->Draw("SAME");

    TF1* fit2 = new TF1("fit2", "[1]-[0] + ([3]-[2]) * x - [4] * x * x", 0, 6500);
    fit2->SetParameters(G_Calibration[current_detector]->GetFunction(function.c_str())->GetParameter(0), linear->GetParameter(0), G_Calibration[current_detector]->GetFunction(function.c_str())->GetParameter(1),  linear->GetParameter(1), G_Calibration[current_detector]->GetFunction(function.c_str())->GetParameter(2));
    fit2->SetLineColor(kGreen);
    fit2->Draw("SAME");

    dir_detector[current_detector]->cd();
    c1->Write();

    vector<double> vec = {G_Calibration[current_detector]->GetFunction("pol2")->GetParameter(0), G_Calibration[current_detector]->GetFunction("pol2")->GetParameter(1), G_Calibration[current_detector]->GetFunction("pol2")->GetParameter(2)};
    // vector<double> vec = {G_Calibration[current_detector]->GetFunction("pol2")->GetParameter(0), G_Calibration[current_detector]->GetFunction("pol2")->GetParameter(1)};
    ManualCalibFitted[current_detector] = vec;


    ///ALPHA
    NUCLEUS = "32Ar";
    counter = 0;
    Reader = new TTreeReader(GROUPED_Tree_Detectors[NUCLEUS][current_detector]);
    Reader->Restart();
    H_Exp[NUCLEUS][current_detector]->Reset();
    TTreeReaderValue<double> ChannelDet1(*Reader, "Channel");
    while (Reader->Next())
    {
        double value = *ChannelDet1 / 1000;
        H_Exp[NUCLEUS][current_detector]->Fill(ManualCalibFitted[current_detector][0] + ManualCalibFitted[current_detector][1] * value + ManualCalibFitted[current_detector][2] * value * value);
    }

    RemoveBKG(H_Exp[NUCLEUS][current_detector], false);

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

        //H_Sim["18N"][current_detector]->Write();
        TF1 *gauss2 = new TF1("gauss2", "gaus");
        gauss2->SetParameter(1, (WindowsMap["18N"][peak][current_detector].second - WindowsMap["18N"][peak][current_detector].first) / 2);
        H_Sim["18N"][current_detector]->Fit("gauss2", "QRN", "", alpha_low, alpha_high);
       

        double x = (-ManualCalibFitted[current_detector][1] + sqrt(pow(ManualCalibFitted[current_detector][1], 2) - 4*ManualCalibFitted[current_detector][2] * (ManualCalibFitted[current_detector][0]-alpha) ))/(2*ManualCalibFitted[current_detector][2]);
        G_Calibration_Alpha[current_detector]->SetPoint(counter, x, gauss2->GetParameter(1));
        G_Calibration_Alpha[current_detector]->SetPointError(counter, (alpha_error/coef), gauss2->GetParError(1));
        counter++;
    }

    TCanvas *c2 = new TCanvas((detectorName[current_detector] + "_Calibration_Alpha").c_str(), (detectorName[current_detector] + "_Calibration_Alpha").c_str(), 1920, 1080);
    c2->Divide(1, 2);
    c2->cd(1);
    function = "pol1";
    TF1* fit_alpha = new TF1(function.c_str(), "[0] + [1] * x", 0, 6500);
    G_Calibration_Alpha[current_detector]->Fit(function.c_str(), "Q");
    G_Calibration_Alpha[current_detector]->GetFunction(function.c_str())->SetLineColor(kRed);
    G_Calibration_Alpha[current_detector]->Draw("*AP");
    G_Calibration_Alpha[current_detector]->GetXaxis()->SetTitle("Channel");
    G_Calibration_Alpha[current_detector]->GetYaxis()->SetTitle("Energy [keV]");
    c2->cd(2);
    // draw diff between fet and data
    TGraphErrors *G_Calibration_diff_alpha = (TGraphErrors*)G_Calibration_Alpha[current_detector]->Clone();
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

    //vector<double> vec_alpha = {G_Calibration_Alpha[current_detector]->GetFunction("pol1")->GetParameter(0), G_Calibration_Alpha[current_detector]->GetFunction("pol1")->GetParameter(1)};

    //cout << "ALPHA : " << vec_alpha[0] << " " << vec_alpha[1] << endl;
    vector<double> vec_alpha = {-5.5938, 7.06253e-2*1000};

    ManualCalibFitted_Alpha[current_detector] = vec_alpha;    
}


void PlottingWindows()
    {
        string Nuclei[2] = {"32Ar", "33Ar"};
        for (string nucleus : Nuclei)
        {
            for (int i = 0; i < detectorNum; i++)
            {
                if (IsDetectorSiliStrip(i))
                {
                    if (!H_Exp[nucleus][i] || !H_Sim_Conv[nucleus][i])
                    {
                        continue;
                    }

                    TCanvas *c = new TCanvas((nucleus + "_" + detectorName[i]).c_str(), (nucleus + "_" + detectorName[i]).c_str(), 1920, 1080);
                    c->Divide(CanvasMap[nucleus].first, CanvasMap[nucleus].second);
                    // H_Exp[nucleus][i]->GetXaxis()->SetRangeUser(WindowsMap[nucleus][14][i].first, WindowsMap[nucleus][14][i].second);
                    // H_Sim_Conv[nucleus][i]->GetXaxis()->SetRangeUser(WindowsMap[nucleus][14][i].first, WindowsMap[nucleus][14][i].second);
                    // H_Sim_Conv[nucleus][i]->Scale(H_Exp[nucleus][i]->Integral() / H_Sim_Conv[nucleus][i]->Integral());
                    for (int peak = 1; peak < CanvasMap[nucleus].first * CanvasMap[nucleus].second; peak++)
                    {
                        c->cd(peak);

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
                }
            }
        
    }
}



#endif