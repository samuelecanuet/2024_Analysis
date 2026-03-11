#include "Detectors.hh"

int VERBOSE = 0;

// TFiles
map<string, TFile*> MERGED_File;
map<string, TFile*> SIMULATED_File;

//Calibrations
TF1* Calibration_Function[SIGNAL_MAX];
TFitResultPtr Calibration_Result[SIGNAL_MAX];

TTreeReader *Reader;
TTreeReaderArray<Signal> *Silicon;
TTreeReaderValue<vector<vector<pair<Signal, Signal>>>> *SiPM_Groups;
TTreeReaderValue<Signal> *HRS;

TFile *FINAL_FILE;

// Silicon Spectrum 
map<string, map<string, TH1D*>> H_Exp;
map<string, map<string, TH1D*>> H_Exp_Coinc;

//Realse Curves
map<string, map<int, map<string, TH1D*>>> H_Release;
map<string, map<int, map<string, TH1D*>>> H_Release_Coinc;

// Beta Spectrum
map<string, map<int, map<string, TH1D*>>> H_Exp_Beta;

// DeltaE
TGraphErrors *G_Sim_DeltaE_Proton = new TGraphErrors();
TGraphErrors *G_Sim_DeltaE_Alpha = new TGraphErrors();
map<string, map<int, TGraphErrors*>> G_DeltaE;

// Utilities
map<string, TDirectory*> dir_nuclei;

vector<string> Directions = {"Up", "Down"};
vector<string> Nuclei;
int start_gate = -20;
int end_gate = 40;

void InitCalibration()
{
    TFile *f = MyTFile((DIR_ROOT_DATA_CALIBRATED + "Calibrated_"+ to_string(YEAR) + "_new.root").c_str(), "READ");
    for (int det = 1;  det < SIGNAL_MAX; det++)
    {
        Calibration_Function[det] = (TF1*)f->Get(("Calibration_" + detectorName[det]).c_str());
        Calibration_Result[det] = (TFitResult*)f->Get(("Calibration_R_" + detectorName[det]).c_str());
    }
    f->Close();
}

void InitExperimentalSpectrum()
{
    TFile *f = MyTFile((DIR_ROOT_DATA_CALIBRATED + "Calibrated_" + to_string(YEAR) + "_new.root").c_str(), "READ");
    for (string Nucleus : Nuclei)
    {
        for (int det = 1; det < SIGNAL_MAX; det++)
        {
            if (IsDetectorSiliStrip(det))
            {
                string dir = det < 50 ? "Up" : "Down";
                if (H_Exp[Nucleus][dir] == nullptr)
                    H_Exp[Nucleus][dir] = (TH1D *)f->Get((detectorName[det] + "/" + Nucleus + "/H_Exp_" + Nucleus + "_" + detectorName[det]).c_str());
                else
                    H_Exp[Nucleus][dir]->Add((TH1D *)f->Get((detectorName[det] + "/" + Nucleus + "/H_Exp_" + Nucleus + "_" + detectorName[det]).c_str()));
            }
        }
    }
}

double gaussBortelsCollaers(double *x, double *p)
{
    double A = p[0];
    double mu = p[1];
    double sigma = p[2]; 
    double eta = p[3];
    double tau1 = p[4];
    double tau2 = p[5];

    double first  = (1 - eta) / tau1    * exp((x[0]-mu)/tau1 + pow(sigma, 2)/(2*pow(tau1, 2))) * erfc(1/sqrt(2) * (x[0] - mu)/sigma + sigma/tau1);
    double second =      eta  / tau2    * exp((x[0]-mu)/tau2 + pow(sigma, 2)/(2*pow(tau2, 2))) * erfc(1/sqrt(2) * (x[0] - mu)/sigma + sigma/tau2);
    // double second = 0;

    return A/2 * (first+second);
}

double gauss(double *x, double *p)
{
    double A = p[0];
    double mu = p[1];
    double sigma = p[2];

    return A * exp(-pow(x[0] - mu, 2) / (2 * pow(sigma, 2))) / (sigma * sqrt(2 * M_PI));
}

double gauss_gaussBortelsCollaers(double *x, double *p)
{
    double A1 = p[0];
    double mu1 = p[1];
    double sigma1 = p[2]; 
    double eta1 = p[3];
    double tau11 = p[4];
    double tau21 = p[5];

    double params[6] = {A1, mu1, sigma1, eta1, tau11, tau21};

    double A2 = p[6];
    double mu2 = p[7] + mu1;
    double sigma2 = p[8]; 

    double params_2[3] = {A2, mu2, sigma2};

    return gaussBortelsCollaers(x, params) + gauss(x, params_2);
}

void InitDirectionDeltaEnergy(bool Analysing_Data)
{
    if (!Analysing_Data)
    {
        TFile *f = MyTFile((DIR_ROOT_DATA_SIMULATED + "DeltaE.root").c_str(), "READ");
        G_Sim_DeltaE_Proton = (TGraphErrors *)f->Get("G_Sim_DeltaE_Proton");
        G_Sim_DeltaE_Alpha = (TGraphErrors *)f->Get("G_Sim_DeltaE_Alpha");
        f->Close();
        if (G_Sim_DeltaE_Proton == nullptr || G_Sim_DeltaE_Alpha == nullptr)
        {
            Error("G_Sim_DeltaE_Proton or G_Sim_DeltaE_Alpha not found in DeltaE.root");
            return;
        }
        else
        {
            Info("ΔE Simulation data loaded");
        }
    }
    else
    {
        TFile *f_res = MyTFile((DIR_ROOT_DATA_SIMULATED + "DeltaE.root").c_str(), "RECREATE");

        TGraphErrors *G_Sim_ELoss_Proton_Up = new TGraphErrors();
        TGraphErrors *G_Sim_ELoss_Proton_Down = new TGraphErrors();
        TGraphErrors *G_Sim_ELoss_Alpha_Up = new TGraphErrors();
        TGraphErrors *G_Sim_ELoss_Alpha_Down = new TGraphErrors();
        TGraphErrors *G_Sim_Eff_Proton_Down = new TGraphErrors();
        TGraphErrors *G_Sim_Eff_Proton_Up = new TGraphErrors();
        TGraphErrors *G_Sim_Eff_Alpha_Down = new TGraphErrors();
        TGraphErrors *G_Sim_Eff_Alpha_Up = new TGraphErrors();
        for (double e = 0.1; e <= 7.0; e+=0.05)
        {
            cout << "Processing energy: " << e << " MeV" << endl;
            // sstream with 1 digit
            stringstream ss;
            ss << fixed << setprecision(2) << e;
            string e_str = ss.str();

            double energy = e * 1000; // Convert MeV to keV

            string path = "/run/media/local1/DATANEX/Samuel-G4/eloss/";//DIR_ROOT_DATA_SIMULATED
            //Protons
            TFile *f = MyTFile((path + "proton_" + e_str + "MeV.root").c_str(), "READ");
            if (f == nullptr)
                continue;
            
            TH1D *hUp = nullptr;
            TH1D *hDown = nullptr;

            for (int det = 1; det <= SIGNAL_MAX; det++)
            {
                if (IsDetectorSiliStrip(det))
                {
                    if (GetDetector(det) < 5)
                    {
                        if (hUp == nullptr)
                        {
                            hUp = (TH1D *)f->Get((detectorName[det] + "_single").c_str());
                        }
                        else
                        {
                            hUp->Add((TH1D *)f->Get((detectorName[det] + "_single").c_str()));
                        }
                    }
                    else if (GetDetector(det) >= 5)
                    {
                        if (hDown == nullptr)
                        {
                            hDown = (TH1D *)f->Get((detectorName[det] + "_single").c_str());
                        }
                        else
                        {
                            hDown->Add((TH1D *)f->Get((detectorName[det] + "_single").c_str()));
                        }
                    }
                }
            }

            hUp->GetXaxis()->SetRangeUser(energy - 25, energy + 25);
            hDown->GetXaxis()->SetRangeUser(energy - 25, energy + 25);

            double peak_position_Down = hDown->GetBinCenter(hDown->GetMaximumBin());
            double peak_position_Up = hUp->GetBinCenter(hUp->GetMaximumBin());

            /// deltaE
            G_Sim_DeltaE_Proton->AddPoint(energy, abs(peak_position_Up - peak_position_Down));
            G_Sim_DeltaE_Proton->SetPointError(G_Sim_DeltaE_Proton->GetN() - 1, 0, sqrt(pow(hUp->GetMeanError(), 2) + pow(hDown->GetMeanError(), 2)));

            // eloss
            G_Sim_ELoss_Proton_Down->AddPoint(peak_position_Down, energy - peak_position_Down);
            G_Sim_ELoss_Proton_Down->SetPointError(G_Sim_ELoss_Proton_Down->GetN() - 1, hDown->GetMeanError(), hDown->GetMeanError());
            G_Sim_ELoss_Proton_Up->AddPoint(peak_position_Up, energy - peak_position_Up);
            G_Sim_ELoss_Proton_Up->SetPointError(G_Sim_ELoss_Proton_Up->GetN() - 1, hUp->GetMeanError(), hUp->GetMeanError());

            

            // efficiency
            double N_event = 10e6;
            G_Sim_Eff_Proton_Down->AddPoint(energy, hDown->Integral() / N_event);
            G_Sim_Eff_Proton_Down->SetPointError(G_Sim_Eff_Proton_Down->GetN() - 1, 0, sqrt(pow(sqrt(hDown->Integral()) / N_event, 2) + pow(hDown->Integral() * sqrt(N_event) / (N_event * N_event), 2)));
            G_Sim_Eff_Proton_Up->AddPoint(energy, hUp->Integral() / N_event);
            G_Sim_Eff_Proton_Up->SetPointError(G_Sim_Eff_Proton_Up->GetN() - 1, 0, sqrt(pow(sqrt(hUp->Integral()) / N_event, 2) + pow(hUp->Integral() * sqrt(N_event) / (N_event * N_event), 2)));


            // // fitting eloss
            // f_res->cd();
            // hDown->GetXaxis()->SetRangeUser(-1111, -1111);
            // TF1 *f_FitPeak = new TF1("f", gauss_gaussBortelsCollaers, 0, 7000, 9);

            // //increasing call limit 
            // ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);

            // f_FitPeak->SetParameter(0, 1);    // Amplitude
            // f_FitPeak->SetParLimits(0, 0, 100000); // Amplitude
            // f_FitPeak->SetParameter(1, hDown->GetBinCenter(hDown->GetMaximumBin())); // Mean
            // f_FitPeak->SetParameter(2, 10);  // Sigma
            // f_FitPeak->SetParLimits(2, 0, 15); // Sigma
            // f_FitPeak->FixParameter(3, 0.);  // Eta
            // // f_FitPeak->SetParLimits(3, 0, 1); // Eta
            // f_FitPeak->SetParameter(4, 100);  // Tau1
            // f_FitPeak->FixParameter(5, 1);  // Tau2
            // f_FitPeak->SetParameter(6, 1);    // Amplitude gauss
            // f_FitPeak->SetParLimits(6, 0, 100000); // Amplitude gauss
            // f_FitPeak->SetParameter(7, -50); // Mean shift gauss
            // f_FitPeak->SetParLimits(7, -1000, 5); // Mean shift gauss
            // f_FitPeak->SetParameter(8, 10);   // Sigma gauss
            // f_FitPeak->SetParLimits(8, 0, 50); // Sigma gauss

            // TFitResultPtr r = hDown->Fit(f_FitPeak, "S MULTITHREAD", "");
            // if (r != 0)
            // {
            //     Warning("Fitting failed for energy: " + to_string(energy) + " keV");
            //     continue;
            // }

            // hDown->SetName(("hDown_" + e_str + "MeV").c_str());
            // TCanvas *c = new TCanvas(("c_" + e_str + "MeV").c_str(), ("c_" + e_str + "MeV").c_str(), 800, 600);
            // TH1D *clone = (TH1D *)hDown->Clone(("clone_" + e_str + "MeV").c_str());
            // clone->SetTitle(("Eloss Down " + e_str + " MeV").c_str());
            // clone->Draw("HIST");    
            // f_FitPeak->SetLineColor(kRed);
            // f_FitPeak->Draw("same");
            // c->Write();


            f->Close();

            //Alphas
            hUp = nullptr;
            hDown = nullptr;
            f = MyTFile((path + "alpha_" + e_str + "MeV.root").c_str(), "READ");
            if (f == nullptr)
                continue;
            
            for (int det = 1; det <= SIGNAL_MAX; det++)
            {
                if (IsDetectorSiliStrip(det))
                {
                    if (GetDetector(det) < 5)
                    {
                        if (hUp == nullptr)
                        {
                            hUp = (TH1D *)f->Get((detectorName[det] + "_single").c_str());
                        }
                        else
                        {
                            hUp->Add((TH1D *)f->Get((detectorName[det] + "_single").c_str()));
                        }
                    }
                    else if (GetDetector(det) >= 5)
                    {
                        if (hDown == nullptr)
                        {
                            hDown = (TH1D *)f->Get((detectorName[det] + "_single").c_str());
                        }
                        else
                        {
                            hDown->Add((TH1D *)f->Get((detectorName[det] + "_single").c_str()));
                        }
                    }
                }

            }

            
            hUp->GetXaxis()->SetRangeUser(hUp->GetBinCenter(hUp->GetMaximumBin()) - 20, hUp->GetBinCenter(hUp->GetMaximumBin()) + 20);
            hDown->GetXaxis()->SetRangeUser(hDown->GetBinCenter(hDown->GetMaximumBin()) - 20, hDown->GetBinCenter(hDown->GetMaximumBin()) + 20);

            G_Sim_DeltaE_Alpha->AddPoint(energy, abs(hUp->GetMean() - hDown->GetMean()));
            G_Sim_DeltaE_Alpha->SetPointError(G_Sim_DeltaE_Alpha->GetN() - 1, 0, sqrt(pow(hUp->GetMeanError(), 2) + pow(hDown->GetMeanError(), 2)));

            G_Sim_ELoss_Alpha_Down->AddPoint(hDown->GetMean(), energy - hDown->GetMean());
            G_Sim_ELoss_Alpha_Down->SetPointError(G_Sim_ELoss_Alpha_Down->GetN() - 1, hDown->GetMeanError(), hDown->GetMeanError());
            G_Sim_ELoss_Alpha_Up->AddPoint(hUp->GetMean(), energy - hUp->GetMean());
            G_Sim_ELoss_Alpha_Up->SetPointError(G_Sim_ELoss_Alpha_Up->GetN() - 1, hUp->GetMeanError(), hUp->GetMeanError());

            f->Close();
        }

        f_res->cd();
        G_Sim_DeltaE_Proton->SetName("G_Sim_DeltaE_Proton");
        G_Sim_DeltaE_Proton->Write();
        G_Sim_DeltaE_Alpha->SetName("G_Sim_DeltaE_Alpha");
        G_Sim_DeltaE_Alpha->Write();
        // Eloss
        G_Sim_ELoss_Proton_Up->SetName("G_Sim_ELoss_Proton_Up");
        G_Sim_ELoss_Proton_Up->Write();
        G_Sim_ELoss_Proton_Down->SetName("G_Sim_ELoss_Proton_Down");
        G_Sim_ELoss_Proton_Down->Write();
        G_Sim_ELoss_Alpha_Up->SetName("G_Sim_ELoss_Alpha_Up");
        G_Sim_ELoss_Alpha_Up->Write();
        G_Sim_ELoss_Alpha_Down->SetName("G_Sim_ELoss_Alpha_Down");
        G_Sim_ELoss_Alpha_Down->Write();
        
        //Efficiency
        G_Sim_Eff_Proton_Down->SetName("G_Sim_Eff_Proton_Down");
        G_Sim_Eff_Proton_Down->Write();
        G_Sim_Eff_Proton_Up->SetName("G_Sim_Eff_Proton_Up");
        G_Sim_Eff_Proton_Up->Write();
        // Ratio Efficiency with 3350keV
        TGraphErrors *G_Sim_Coef_Proton_Down = new TGraphErrors();
        TGraphErrors *G_Sim_Coef_Proton_Up = new TGraphErrors();
        for (int i = 0; i < G_Sim_Eff_Proton_Down->GetN(); i++)
        {
            double x, y, ex, ey;
            G_Sim_Eff_Proton_Down->GetPoint(i, x, y);
            ex = G_Sim_Eff_Proton_Down->GetErrorX(i);
            ey = G_Sim_Eff_Proton_Down->GetErrorY(i);
            G_Sim_Coef_Proton_Down->AddPoint(x, y / G_Sim_Eff_Proton_Down->Eval(3350));
            G_Sim_Coef_Proton_Down->SetPointError(G_Sim_Coef_Proton_Down->GetN() - 1, ex, sqrt(pow(ey / G_Sim_Eff_Proton_Down->Eval(3350), 2) + pow(ey * y / pow(G_Sim_Eff_Proton_Down->Eval(3350), 2), 2)));   
        }
        for (int i = 0; i < G_Sim_Eff_Proton_Up->GetN(); i++)
        {
            double x, y, ex, ey;
            G_Sim_Eff_Proton_Up->GetPoint(i, x, y);
            ex = G_Sim_Eff_Proton_Up->GetErrorX(i);
            ey = G_Sim_Eff_Proton_Up->GetErrorY(i);
            G_Sim_Coef_Proton_Up->AddPoint(x, y / G_Sim_Eff_Proton_Up->Eval(3350));
            G_Sim_Coef_Proton_Up->SetPointError(G_Sim_Coef_Proton_Up->GetN() - 1, ex, sqrt(pow(ey / G_Sim_Eff_Proton_Up->Eval(3350), 2) + pow(ey * y / pow(G_Sim_Eff_Proton_Up->Eval(3350), 2), 2)));   
        }
        G_Sim_Coef_Proton_Down->SetName("G_Sim_Coef_Proton_Down");
        G_Sim_Coef_Proton_Down->Write();
        G_Sim_Coef_Proton_Up->SetName("G_Sim_Coef_Proton_Up");
        G_Sim_Coef_Proton_Up->Write();
        f_res->Close();
        Info("ΔE Simulation data loaded");
    }
}

void InitHistograms(string Nucleus)
{
    Info("Initializing histograms for " + Nucleus, 1);
    for (string dir : Directions)
    {
        H_Release[Nucleus][0][dir] = new TH1D(("H_Release_" + Nucleus + "_" + dir).c_str(), ("H_Release_" + Nucleus + "_" + dir).c_str(), 1200, 0, 1.2);
        H_Release[Nucleus][0][dir]->GetXaxis()->SetTitle("Time (s)");
        H_Release[Nucleus][0][dir]->GetYaxis()->SetTitle("Counts / ms");
        H_Release[Nucleus][0][dir]->GetXaxis()->CenterTitle();
        H_Release[Nucleus][0][dir]->GetYaxis()->CenterTitle();

        H_Release_Coinc[Nucleus][0][dir] = new TH1D(("H_Release_Coinc_" + Nucleus + "_" + dir).c_str(), ("H_Release_Coinc_" + Nucleus + "_" + dir).c_str(), 1200, 0, 1.2);
        H_Release_Coinc[Nucleus][0][dir]->GetXaxis()->SetTitle("Time (s)");
        H_Release_Coinc[Nucleus][0][dir]->GetYaxis()->SetTitle("Counts / ms");
        H_Release_Coinc[Nucleus][0][dir]->GetXaxis()->CenterTitle();
        H_Release_Coinc[Nucleus][0][dir]->GetYaxis()->CenterTitle();

        H_Exp[Nucleus][dir] = new TH1D(("Spectrum_" + Nucleus + "_" + dir).c_str(), ("Spectrum_" + Nucleus + "_" + dir).c_str(), eSiliN_cal/10, eSiliMin_cal, eSiliMax_cal);
        H_Exp[Nucleus][dir]->GetXaxis()->SetTitle("Energy (keV)");
        H_Exp[Nucleus][dir]->GetYaxis()->SetTitle("Counts / keV");
        H_Exp[Nucleus][dir]->GetXaxis()->CenterTitle();
        H_Exp[Nucleus][dir]->GetYaxis()->CenterTitle();

        H_Exp_Coinc[Nucleus][dir] = new TH1D(("Spectrum_Coinc_" + Nucleus + "_" + dir).c_str(), ("Spectrum_Coinc_" + Nucleus + "_" + dir).c_str(), eSiliN_cal/10, eSiliMin_cal, eSiliMax_cal);
        H_Exp_Coinc[Nucleus][dir]->GetXaxis()->SetTitle("Energy (keV)");
        H_Exp_Coinc[Nucleus][dir]->GetYaxis()->SetTitle("Counts / keV");
        H_Exp_Coinc[Nucleus][dir]->GetXaxis()->CenterTitle();
        H_Exp_Coinc[Nucleus][dir]->GetYaxis()->CenterTitle();

        for (int peak = 1; peak <= CanvasMap[Nucleus].first * CanvasMap[Nucleus].second; peak++)
        {
            H_Release[Nucleus][peak][dir] = new TH1D(("H_Release_" + Nucleus + "_" + dir + "_" + to_string(peak)).c_str(), ("H_Release_" + Nucleus + "_" + dir + "_" + to_string(peak)).c_str(), 1200, 0, 1.2);
            H_Release[Nucleus][peak][dir]->GetXaxis()->SetTitle("Time (s)");
            H_Release[Nucleus][peak][dir]->GetYaxis()->SetTitle("Counts / ms");
            H_Release[Nucleus][peak][dir]->GetXaxis()->CenterTitle();
            H_Release[Nucleus][peak][dir]->GetYaxis()->CenterTitle();

            H_Release_Coinc[Nucleus][peak][dir] = new TH1D(("H_Release_Coinc_" + Nucleus + "_" + dir + "_" + to_string(peak)).c_str(), ("H_Release_Coinc_" + Nucleus + "_" + dir + "_" + to_string(peak)).c_str(), 1200, 0, 1.2);
            H_Release_Coinc[Nucleus][peak][dir]->GetXaxis()->SetTitle("Time (s)");
            H_Release_Coinc[Nucleus][peak][dir]->GetYaxis()->SetTitle("Counts / ms");
            H_Release_Coinc[Nucleus][peak][dir]->GetXaxis()->CenterTitle();
            H_Release_Coinc[Nucleus][peak][dir]->GetYaxis()->CenterTitle();

            H_Exp_Beta[Nucleus][peak][dir] = new TH1D(("Spectrum_Beta_" + Nucleus + "_" + dir + "_" + to_string(peak)).c_str(), ("Spectrum_Beta_" + Nucleus + "_" + dir + "_" + to_string(peak)).c_str(), eLowN, eLowMin, eLowMax);
            H_Exp_Beta[Nucleus][peak][dir]->GetXaxis()->SetTitle("Energy (keV)");
            H_Exp_Beta[Nucleus][peak][dir]->GetYaxis()->SetTitle("Counts / 1keV");
            H_Exp_Beta[Nucleus][peak][dir]->GetXaxis()->CenterTitle();
            H_Exp_Beta[Nucleus][peak][dir]->GetYaxis()->CenterTitle();
        }
    }
}

bool IsCoincidence(double P_Time, vector<vector<pair<Signal, Signal>>> SiPM_Groups)
{

  double NEAREST = 1e9;
  int NEAREST_GROUP_INDEX = -1;
  vector<double> nearest_group_time;
  vector<double> mean_group_time = vector<double>(10, 0);
  if ((SiPM_Groups).size() > 0)
  {
    // cout << "Size: " << (SiPM_Groups).size() << endl;
    // Lopping on subgroups
    for (int i_group = 0; i_group < (SiPM_Groups).size(); i_group++)
    {
      nearest_group_time.push_back(1e9);
      int counter_mean_group_time = 0;
      // if ((SiPM_Groups)[i_group].size() > 9)
      // {
      //     for (int i_pair = 0; i_pair < (SiPM_Groups)[i_group].size(); i_pair++)
      //     {
      //         cout << (SiPM_Groups)[i_group][i_pair].first << endl;
      //     }
      // }
      // Looping on pair of the subgroup
      for (int i_pair = 0; i_pair < (SiPM_Groups)[i_group].size(); i_pair++)
      {

        // Taking valid element of pair with High priority
        if ((SiPM_Groups)[i_group][i_pair].first.isValid && !isnan((SiPM_Groups)[i_group][i_pair].first.Time))
        {
          if (abs(nearest_group_time[i_group]) > abs((SiPM_Groups)[i_group][i_pair].first.Time))
          {
            nearest_group_time[i_group] = (SiPM_Groups)[i_group][i_pair].first.Time;
          }

          mean_group_time[i_group] += (SiPM_Groups)[i_group][i_pair].first.Time;
          counter_mean_group_time++;
        }
        else if ((SiPM_Groups)[i_group][i_pair].second.isValid && !isnan((SiPM_Groups)[i_group][i_pair].second.Time))
        {
          if (abs(nearest_group_time[i_group]) > abs((SiPM_Groups)[i_group][i_pair].second.Time))
          {
            nearest_group_time[i_group] = (SiPM_Groups)[i_group][i_pair].second.Time;
          }
          mean_group_time[i_group] += (SiPM_Groups)[i_group][i_pair].first.Time;
          counter_mean_group_time++;
        }
      }
      mean_group_time[i_group] /= counter_mean_group_time;
      // cout << "Mean Time: " << mean_group_time[i_group] << endl;
    }

    // Saving the nearest group in time and index group
    for (int i_group = 0; i_group < nearest_group_time.size(); i_group++)
    {
      if (abs(NEAREST) > abs(nearest_group_time[i_group]))
      {
        NEAREST_GROUP_INDEX = i_group;
        NEAREST = nearest_group_time[i_group];
      }
    }
  }

  // cout << "Nearest: " << NEAREST << "    " << NEAREST_GROUP_INDEX << endl;

  int Multiplicity;
  if (NEAREST_GROUP_INDEX == -1)
    return false;
  else
    Multiplicity = (SiPM_Groups)[NEAREST_GROUP_INDEX].size();
  

  if (Multiplicity >= 3 && (NEAREST > start_gate && NEAREST < end_gate))
  {
    return true;
  }

  return false;
}

void ReadingExperimentalData(string Nucleus)
{
    Info("Reading Experimental Data for " + Nucleus, 1);
    clock_t start = clock(), Current;
    int Entries = Reader->GetEntries();

    double Proton_Pulse = 0;

    while (Reader->Next() && Reader->GetCurrentEntry() < 1e6)
    {
        ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "Reading Tree");

        if ((**HRS).isValid)
        {
            Proton_Pulse = (**HRS).Time * 1e-9;
            continue;
        }

        int Silicon_Label = (*Silicon)[1].Label;
        string dir = Silicon_Label < 50 ? "Up" : "Down";
        double Silicon_Energy = Calibration_Function[Silicon_Label]->Eval((*Silicon)[1].Channel / 1000.);
        double Silicon_Time = (*Silicon)[1].Time * 1e-9;

        // cout << "Proton Pulse: " << Proton_Pulse << "    Silicon Time: " << Silicon_Time << "    Diff: " << Silicon_Time - Proton_Pulse << endl;

        bool Coincidence = IsCoincidence((*Silicon)[1].Time, **SiPM_Groups);

        // All Spectrum
        H_Release[Nucleus][0][dir]->Fill(Silicon_Time - Proton_Pulse);
        H_Exp[Nucleus][dir]->Fill(Silicon_Energy);
        if (Coincidence)
        {
            H_Release_Coinc[Nucleus][0][dir]->Fill(Silicon_Time - Proton_Pulse);
            H_Exp_Coinc[Nucleus][dir]->Fill(Silicon_Energy);
        }

        // Looping through peaks
        for (int peak = 1; peak <= CanvasMap[Nucleus].first * CanvasMap[Nucleus].second; peak++)
        {
            if (WindowsMap[Nucleus][peak][Silicon_Label].first == -1 || !WindowsMap[Nucleus][peak][Silicon_Label].first)
                continue;

            if (WindowsMap[Nucleus][peak][Silicon_Label].first < Silicon_Energy && Silicon_Energy < WindowsMap[Nucleus][peak][Silicon_Label].second)
            {
                H_Release[Nucleus][peak][dir]->Fill(Silicon_Time - Proton_Pulse);
                if (Coincidence)
                {
                    H_Release_Coinc[Nucleus][peak][dir]->Fill(Silicon_Time - Proton_Pulse);
                    if ((**SiPM_Groups)[0].size() == 9)
                    {
                        for (int i = 0; i < (**SiPM_Groups)[0].size(); i++)
                        {
                            if (GetDetectorChannel((**SiPM_Groups)[0][i].first.Label) == 7)
                            {
                                H_Exp_Beta[Nucleus][peak][dir]->Fill((**SiPM_Groups)[0][i].second.Channel);
                            }
                        }
                    }
                }
            }
        }
    }
}

void ComputeDeltaEHist()
{
    Info("Computing ΔE Histograms");
    FINAL_FILE->cd();
    map<string, TGraphErrors*> G_DeltaE_Nuclei;

    TCanvas *canvas_All = new TCanvas("DeltaE", "DeltaE", 800, 600);
    TMultiGraph *mg_All = new TMultiGraph();

    G_Sim_DeltaE_Proton->SetLineColor(kRed);
    G_Sim_DeltaE_Proton->SetLineWidth(2);
    mg_All->Add(G_Sim_DeltaE_Proton, "L");
    G_Sim_DeltaE_Alpha->SetLineColor(kBlue);
    G_Sim_DeltaE_Alpha->SetLineWidth(2);
    mg_All->Add(G_Sim_DeltaE_Alpha, "L");
    mg_All->GetXaxis()->SetTitle("Energy [keV]");
    mg_All->GetYaxis()->SetTitle("#Delta E [keV]");
    
    TLegend *l_all = new TLegend(0.7, 0.7, 0.9, 0.9);

    for (string Nucleus : Nuclei)
    {
        Info(Nucleus, 1);
        dir_nuclei[Nucleus] = FINAL_FILE->mkdir(Nucleus.c_str());
        dir_nuclei[Nucleus]->cd();
        G_DeltaE_Nuclei[Nucleus] = new TGraphErrors();
        for (double peak : PeakList[Nucleus])
        {
            if (WindowsMap[Nucleus][peak][11].first == -1 || !WindowsMap[Nucleus][peak][11].first)
                continue;

            if ( H_Exp[Nucleus]["Up"] == nullptr)
            {
                Error("H_Exp[" + Nucleus + "][Up] is nullptr, skipping peak " + to_string(peak));
            }
            Info("Peak " + to_string(peak), 2);
            H_Exp[Nucleus]["Up"]->GetXaxis()->SetRangeUser(WindowsMap[Nucleus][peak][15].first, WindowsMap[Nucleus][peak][11].second);
            H_Exp[Nucleus]["Down"]->GetXaxis()->SetRangeUser(WindowsMap[Nucleus][peak][55].first, WindowsMap[Nucleus][peak][51].second);
            
            // ## for each peak
            G_DeltaE[Nucleus][peak] = new TGraphErrors();
            G_DeltaE[Nucleus][peak]->AddPoint(H_Exp[Nucleus]["Down"]->GetMean(), abs(H_Exp[Nucleus]["Up"]->GetMean() - H_Exp[Nucleus]["Down"]->GetMean()));
            G_DeltaE[Nucleus][peak]->SetPointError(G_DeltaE[Nucleus][peak]->GetN() - 1, sqrt(pow(H_Exp[Nucleus]["Up"]->GetMeanError(), 2) + pow(H_Exp[Nucleus]["Down"]->GetMeanError(), 2)), sqrt(pow(H_Exp[Nucleus]["Up"]->GetMeanError(), 2) + pow(H_Exp[Nucleus]["Down"]->GetMeanError(), 2)));

            // plotting
            TCanvas *canvas = new TCanvas((Nucleus + "_DeltaE_" + to_string(peak)).c_str(), (Nucleus + "_DeltaE_" + to_string(peak)).c_str(), 800, 600);
            TMultiGraph *mg = new TMultiGraph();
            TLegend *l = new TLegend(0.7, 0.7, 0.9, 0.9);

            mg->Add(G_Sim_DeltaE_Proton, "L");
            mg->Add(G_Sim_DeltaE_Alpha, "L");
            G_DeltaE[Nucleus][peak]->SetMarkerStyle(21);
            G_DeltaE[Nucleus][peak]->SetMarkerColor(NucleusColor[Nucleus]);
            G_DeltaE[Nucleus][peak]->SetDrawOption("AP");
            l->AddEntry(G_DeltaE[Nucleus][peak], "Experimental Peak", "p");
            mg->Add(G_DeltaE[Nucleus][peak], "P");
            mg->SetTitle((Nucleus + " ΔE Peak " + to_string(peak)).c_str());
            mg->GetXaxis()->SetTitle("Energy [keV]");
            mg->GetYaxis()->SetTitle("#Delta E [keV]");
            mg->Draw("A");
            l->Draw("SAME");
            TText *text = new TText(0.8, 0.17, "Proton");
            text->SetNDC();
            text->SetTextColor(kRed);
            text->SetTextSize(0.04);
            text->Draw("SAME");
            TText *text2 = new TText(0.8, 0.35, "Alpha");
            text2->SetNDC();
            text2->SetTextColor(kBlue);
            text2->SetTextSize(0.04);
            text2->Draw("SAME");
            // canvas->Write();

            // ## for all
            G_DeltaE_Nuclei[Nucleus]->AddPoint(H_Exp[Nucleus]["Down"]->GetMean(), abs(H_Exp[Nucleus]["Up"]->GetMean() - H_Exp[Nucleus]["Down"]->GetMean()));
            G_DeltaE_Nuclei[Nucleus]->SetPointError(G_DeltaE_Nuclei[Nucleus]->GetN() - 1, sqrt(pow(H_Exp[Nucleus]["Up"]->GetMeanError(), 2) + pow(H_Exp[Nucleus]["Down"]->GetMeanError(), 2)), sqrt(pow(H_Exp[Nucleus]["Up"]->GetMeanError(), 2) + pow(H_Exp[Nucleus]["Down"]->GetMeanError(), 2)));
        }

        // plotting all
        G_DeltaE_Nuclei[Nucleus]->SetMarkerStyle(21);
        G_DeltaE_Nuclei[Nucleus]->SetMarkerColor(NucleusColor[Nucleus]);
        l_all->AddEntry(G_DeltaE_Nuclei[Nucleus], Nucleus.c_str(), "p");
        // G_DeltaE_Nuclei[Nucleus]->Draw("AP SAME");
        mg_All->Add(G_DeltaE_Nuclei[Nucleus], "P");
    }

    FINAL_FILE->cd();
    canvas_All->cd();

    mg_All->SetTitle("#Delta E All Peaks");
    mg_All->Draw("A");
    l_all->Draw("SAME");
    TText *text = new TText(0.8, 0.17, "Proton");
    text->SetNDC();
    text->SetTextColor(kRed);
    text->SetTextSize(0.04);
    text->Draw("SAME");
    TText *text2 = new TText(0.8, 0.35, "Alpha");
    text2->SetNDC();
    text2->SetTextColor(kBlue);
    text2->SetTextSize(0.04);
    text2->Draw("SAME");
    canvas_All->Write();
}

void PlottingPeak(string Nucleus, double peak)
{
    double E0 = PeakData[Nucleus][peak][0];
    ostringstream oss;
    oss << fixed << setprecision(1) << E0;
    string E0_str = oss.str();
    TCanvas *c = new TCanvas((Nucleus + "_Peak_" + to_string(peak) + "_" + E0_str + "keV").c_str(), (Nucleus + "_Peak_" + to_string(peak) + "_" + E0_str + "keV").c_str(), 1920, 1080);
    c->Divide(2, 2);

    // Up Spectrum Single/Coinc
    c->cd(1);
    TPad *padUp = new TPad("Silicon Spectrum (Up)", "Silicon Spectrum (Up)", 0, 0.5, 1.0, 1.0);
    padUp->Draw();
    padUp->cd();
    string dir = "Up";
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    H_Exp[Nucleus][dir]->GetXaxis()->SetRangeUser(WindowsMap[Nucleus][peak][15].first, WindowsMap[Nucleus][peak][51].second);
    H_Exp[Nucleus][dir]->GetYaxis()->SetRangeUser(1, -1111);
    H_Exp[Nucleus][dir]->SetTitle("Silicon Spectrum (Up)");
    H_Exp[Nucleus][dir]->SetStats(false);
    H_Exp[Nucleus][dir]->Draw("HIST");
    H_Exp_Coinc[Nucleus][dir]->GetXaxis()->SetRangeUser(WindowsMap[Nucleus][peak][15].first, WindowsMap[Nucleus][peak][51].second);
    H_Exp_Coinc[Nucleus][dir]->SetLineColor(kRed);
    H_Exp_Coinc[Nucleus][dir]->SetTitle("Silicon Spectrum (Up)");
    H_Exp_Coinc[Nucleus][dir]->SetStats(false);
    H_Exp_Coinc[Nucleus][dir]->Draw("HIST SAME");
    legend->AddEntry(H_Exp[Nucleus][dir], "Single", "l");
    legend->AddEntry(H_Exp_Coinc[Nucleus][dir], "Coincidence", "l");
    legend->Draw("SAME");

    c->cd(1);
    TPad *padDown = new TPad("Silicon Spectrum (Down)", "Silicon Spectrum (Down)", 0, 0, 1.0, 0.5);
    padDown->Draw();
    padDown->cd();
    dir = "Down";
    H_Exp[Nucleus][dir]->GetXaxis()->SetRangeUser(WindowsMap[Nucleus][peak][15].first, WindowsMap[Nucleus][peak][51].second);
    H_Exp[Nucleus][dir]->GetYaxis()->SetRangeUser(1, -1111);
    H_Exp[Nucleus][dir]->SetTitle("Silicon Spectrum (Down)");
    H_Exp[Nucleus][dir]->SetStats(false);
    H_Exp[Nucleus][dir]->Draw("HIST");
    H_Exp_Coinc[Nucleus][dir]->GetXaxis()->SetRangeUser(WindowsMap[Nucleus][peak][15].first, WindowsMap[Nucleus][peak][51].second);
    H_Exp_Coinc[Nucleus][dir]->SetLineColor(kRed);
    H_Exp_Coinc[Nucleus][dir]->SetTitle("Silicon Spectrum (Down)");
    H_Exp_Coinc[Nucleus][dir]->SetStats(false);
    H_Exp_Coinc[Nucleus][dir]->Draw("HIST SAME");
    legend->Draw("SAME");

    // Particle Identification
    c->cd(2);
    TMultiGraph *mg = new TMultiGraph();
    TLegend *l = new TLegend(0.7, 0.7, 0.9, 0.9);
    mg->Add(G_Sim_DeltaE_Proton, "L");
    mg->Add(G_Sim_DeltaE_Alpha, "L");
    G_DeltaE[Nucleus][peak]->SetMarkerStyle(21);
    G_DeltaE[Nucleus][peak]->SetMarkerColor(NucleusColor[Nucleus]);
    G_DeltaE[Nucleus][peak]->SetDrawOption("AP");
    l->AddEntry(G_DeltaE[Nucleus][peak], "Experimental Peak", "p");
    mg->Add(G_DeltaE[Nucleus][peak], "P");
    mg->SetTitle("Particle identification");
    mg->GetXaxis()->SetTitle("Energy (keV)");
    mg->GetYaxis()->SetTitle("#Delta E (keV)");
    mg->Draw("A");
    l->Draw("SAME");
    TText *text = new TText(0.8, 0.17, "Proton");
    text->SetNDC();
    text->SetTextColor(kRed);
    text->SetTextSize(0.04);
    text->Draw("SAME");
    TText *text2 = new TText(0.8, 0.35, "Alpha");
    text2->SetNDC();
    text2->SetTextColor(kBlue);
    text2->SetTextSize(0.04);
    text2->Draw("SAME");

    // Release curve
    c->cd(3);
    TLegend *legend_release = new TLegend(0.7, 0.7, 0.9, 0.9);
    H_Release[Nucleus][peak]["Down"]->Add(H_Release[Nucleus][peak]["Up"]);
    H_Release_Coinc[Nucleus][peak]["Down"]->Add(H_Release_Coinc[Nucleus][peak]["Up"]);

    H_Release[Nucleus][peak]["Down"]->GetXaxis()->SetRangeUser(0, 1.2);
    H_Release[Nucleus][peak]["Down"]->SetTitle("Release Curve");
    H_Release[Nucleus][peak]["Down"]->SetStats(false);
    H_Release[Nucleus][peak]["Down"]->Draw("HIST");
    legend_release->AddEntry(H_Release[Nucleus][peak]["Down"], "Single", "l");

    H_Release_Coinc[Nucleus][peak]["Down"]->GetXaxis()->SetRangeUser(0, 1.2);
    H_Release_Coinc[Nucleus][peak]["Down"]->SetLineColor(kRed);
    H_Release_Coinc[Nucleus][peak]["Down"]->SetTitle("Release Curve in Coincidence");
    H_Release_Coinc[Nucleus][peak]["Down"]->SetStats(false);
    H_Release_Coinc[Nucleus][peak]["Down"]->Draw("HIST SAME");
    legend_release->AddEntry(H_Release[Nucleus][peak]["Down"], "In Coincidence", "l");

    TH1D *HIAS = (TH1D *)H_Release[Nucleus][IAS[Nucleus]]["Down"]->Clone(("H_Release_IAS_" + Nucleus + "_" + to_string(peak)).c_str());
    if (HIAS->Integral() == 0 || H_Release_Coinc[Nucleus][peak]["Down"]->Integral() == 0)
        return;
    
    HIAS->Scale(H_Release_Coinc[Nucleus][peak]["Down"]->Integral() / HIAS->Integral());
    HIAS->SetLineColor(kBlack);
    HIAS->Draw("HIST SAME");
    legend_release->AddEntry(HIAS, "IAS", "l");
    legend_release->Draw("SAME");

    // Beta Spectrum
    c->cd(4);
    TLegend *legend_beta = new TLegend(0.7, 0.7, 0.9, 0.9);

    for (string dir : Directions)
    {
        // H_Exp_Beta[Nucleus][peak][dir]->GetXaxis()->SetRangeUser(1, WindowsBetaMap[Nucleus][peak]);
        H_Exp_Beta[Nucleus][peak][dir]->GetXaxis()->SetRangeUser(H_Exp_Beta[Nucleus][peak][dir]->GetBinWidth(1), -1111);
        H_Exp_Beta[Nucleus][peak][dir]->SetTitle("#beta Spectrum");
        H_Exp_Beta[Nucleus][peak][dir]->SetStats(false);
        double factor = (double)Freedman_Diaconis(H_Exp_Beta[Nucleus][peak][dir]);
        H_Exp_Beta[Nucleus][peak][dir]->Rebin((int)factor);
        H_Exp_Beta[Nucleus][peak][dir]->GetYaxis()->SetTitle(("Counts / " + to_string(0.1 * factor) + "keV").c_str());
        if (dir == "Up") H_Exp_Beta[Nucleus][peak][dir]->Draw("HIST");
        else H_Exp_Beta[Nucleus][peak][dir]->Draw("HIST SAME");

        legend_beta->AddEntry(H_Exp_Beta[Nucleus][peak][dir], (dir).c_str(), "l");
    }
    legend_beta->Draw("SAME");

    c->Write();
}