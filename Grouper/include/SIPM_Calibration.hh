#ifndef SIPM_CALIBRATION_HH
#define SIPM_CALIBRATION_HH

#include "../../../lib/SignalDict/Signal.h"
#include "Detectors.hh"

/// FILE ///
map<string, TFile *> GROUPED_File;
map<string, TFile *> SIMULATED_File;
TFile *CALIBRATED_File;
TFile *f_tree;

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
string Nucleis[3] = {"32Ar", "207Bi", "90Sr"};
string NUCLEUS;
double SiliconCalibrationParameter[SIGNAL_MAX][3];
map<string, TF1 *[SIGNAL_MAX]> MatchingLowHigh;
map<string, TF1 *[SIGNAL_MAX]> MatchingLowHigh_erf;
map<string, TF1 *[SIGNAL_MAX]> MatchingSiPM;
map<string, pair<double, double>> SiPM_Range;



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

                // init trees
                for (int det = 1; det <= 9; det++)
                {
                    f_tree->cd();
                    Tree_Peaks[NUCLEUS][number][det] = new TTree(("Tree_Peaks_" + NUCLEUS + "_" + to_string(number) + "_" + to_string(det)).c_str(), ("Tree_Peaks_" + NUCLEUS + "_" + to_string(number) + "_" + to_string(det)).c_str());
                    Tree_Peaks[NUCLEUS][number][det]->Branch("SiPM", &SiPMi);
                }
            }
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

void ReadTree()
{
    Info("Read Tree");
    NUCLEUS="32Ar";
    f_tree = new TFile((DIR_ROOT_DATA_CALIBRATED + "Matching_SiPM_trees.root").c_str(), "READ");
    for (int peak = 0; peak <= 50; peak++)
    {
        for (int det = 1; det <= 9; det++)
        {
            if (WindowsMap[NUCLEUS][peak][11].first == -1 || !WindowsMap[NUCLEUS][peak][11].first)
                        continue;
            Tree_Peaks[NUCLEUS][peak][det] = (TTree *)f_tree->Get(("Tree_Peaks_" + NUCLEUS + "_" + to_string(peak) + "_" + to_string(det)).c_str());
        }
    }
    for (int det = 1; det <= 9; det++)
    {
        Tree_Peaks["207Bi"][0][det] = (TTree *)f_tree->Get(("Tree_Peaks_207Bi_" + to_string(0) + "_" + to_string(det)).c_str());
        Tree_Peaks["90Sr"][0][det] = (TTree *)f_tree->Get(("Tree_Peaks_90Sr_" + to_string(0) + "_" + to_string(det)).c_str());
    }

}

void InitHistograms()
{
    SiPM_Range["32Ar"] = make_pair(0, 6000);
    SiPM_Range["207Bi"] = make_pair(0, 2000);
    SiPM_Range["90Sr"] = make_pair(0, 4000);
    Info("Init Histograms");
    for (string Nucleus : Nucleis)
    {
        NUCLEUS = Nucleus;
        dir[NUCLEUS] = CALIBRATED_File->mkdir(NUCLEUS.c_str());
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorBetaHigh(i))
            {
                       
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

                if (NUCLEUS == "207Bi" || NUCLEUS == "90Sr")
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

        if (NUCLEUS == "207Bi" || NUCLEUS == "90Sr")
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

        if (NUCLEUS == "207Bi" || NUCLEUS == "90Sr")
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
            H_SiPM_Calibrated[NUCLEUS][peak_number][GetDetectorChannel((**SiPM).Label)]->Fill(energy);
        }
        
        cout << "convoluting resolution" << endl;
        // convoluting resolution on histogram
        for (int i = 1; i <= H_Sim[NUCLEUS][peak_number]->GetNbinsX(); i++)
        {
            if (H_Sim[NUCLEUS][peak_number]->GetBinContent(i) == 0)
                continue;
            energy = H_Sim[NUCLEUS][peak_number]->GetBinCenter(i);
            double sigma_resolution = sqrt(pow(Resolution_OffSet, 2) + pow(Resolution_SQRT * sqrt(energy), 2) + pow(Resolution_2 * pow(energy, 2), 2));
            // double sigma_resolution = sqrt(pow(coefficents[2].second, 2) + pow(coefficents[3].second * sqrt(energy), 2));

            gauss = new TF1("gauss", "gaus", 0, eSiliMax_cal);
            gauss->SetNpx(eSiliN_cal);
            gauss->SetParameters(H_Sim[NUCLEUS][peak_number]->GetBinContent(i) / (sigma_resolution * sqrt(2 * M_PI)), energy, sigma_resolution); /// weighted gaussian

            H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Add(gauss->GetHistogram());

            delete gauss;
        }
        ///////////////
        // random convoluting resolution on histogram
        // for (int i = 0; i < H_Sim[NUCLEUS][peak_number]->GetEntries() && i < 10000000; i++)
        // {
        //     energy = H_Sim[NUCLEUS][peak_number]->GetRandom();
        //     double sigma_resolution = sqrt(pow(Resolution_OffSet, 2) + pow(Resolution_SQRT * sqrt(energy), 2) + pow(Resolution_2 * pow(energy, 2), 2));
        //     // shoot in a gaussian with ramdom engine
        //     H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Fill(gRandom->Gaus(energy, sigma_resolution));
        // }
        ///////////////

        cout << "convoluting threshold" << endl;

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
        // Threshold_f = new TF1("Threshold", "0.5*(1+erf((x-[0])/[1]))", 0, 8000);
        // Threshold_f->SetParameters(Threshold, Threshold_STD);

        // for (int i = 1; i <= H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetNbinsX(); i++)
        // {
        //     H_Sim_Conv[NUCLEUS][peak_number][current_detector]->SetBinContent(i, H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetBinContent(i) * Threshold_f->Eval(H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetBinCenter(i)));
        // }

        // delete Threshold_f;

        //////////////// CHI2 ////////////////
        
        // NUCLEUS = "207Bi";
        H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector]->GetXaxis()->SetRangeUser(300, 700);
        H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetXaxis()->SetRangeUser(300, 700);
        double chi2 = H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector]->Chi2Test(H_Sim_Conv[NUCLEUS][peak_number][current_detector], "CHI2/NDF");
        cout << chi2 << "   " << Calibration_OffSet << "   " << Calibration << "   " << Resolution_OffSet << "   " << Resolution_SQRT << "   " << Resolution_2 << "   " << Threshold << "   " << Threshold_STD << endl;

        H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Scale(H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector]->Integral()/H_Sim_Conv[NUCLEUS][peak_number][current_detector]->Integral() );
        H_SiPM_Calibrated[NUCLEUS][peak_number][current_detector]->GetXaxis()->SetRangeUser(0, SiPM_Range[NUCLEUS].second);
        H_Sim_Conv[NUCLEUS][peak_number][current_detector]->GetXaxis()->SetRangeUser(0, SiPM_Range[NUCLEUS].second);

        
    }  
    return 0; 
}

void CalibrationSiPM()
{
    res = {0, 8.17, 6.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5}; //#2 8.17
    cout << "Calibration SiPM : " << NUCLEUS << endl;
    Minimizer *minimizer = Factory::CreateMinimizer("Minuit2", "Migrad");
    ROOT::Math::Functor functor(&Chi2TreeHist_conv, 7);
    vector<double> Par = {30.3, 1.7, 0, res[current_detector], 0.e-05, 70, 15.0801};
    minimizer->SetFunction(functor);
    minimizer->SetLimitedVariable(0, "Calibration_OffSet", Par[0], 10, 0, 100);
    minimizer->SetLimitedVariable(1, "Calibration", Par[1], 0.1, 0.5, 4);
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
    const double *bestPar = Par.data();
    // minimizer->Minimize();

    // const double *bestPar = minimizer->X();
    double chi2 = Chi2TreeHist_conv(bestPar);
}

void WriteHistograms()
{
    gStyle->SetOptStat(0);
    Info("Write Histograms");
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
    

    NUCLEUS="207Bi";
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
    dir[NUCLEUS]->cd();

    TCanvas *cExp_Sim_all_90Sr = new TCanvas(("cExp_Sim1_" + NUCLEUS + "_Peak_" + to_string(0)).c_str(), ("cExp_Sim1_" + NUCLEUS + "_Peak_" + to_string(0)).c_str(), 800, 800);
    cExp_Sim_all_90Sr->Divide(3, 3);
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

        cExp_Sim_all_90Sr->cd(det);
        H_SiPM_Calibrated[NUCLEUS][peak][det]->Draw("HIST");
        H_Sim_Conv[NUCLEUS][peak][det]->SetLineColor(kRed);
        H_Sim_Conv[NUCLEUS][peak][det]->Draw("HIST SAME");
    }
    cExp_Sim_all_90Sr->Write();
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