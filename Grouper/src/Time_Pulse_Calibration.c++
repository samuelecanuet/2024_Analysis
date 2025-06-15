#include "Time_Pulse_Calibration.hh"

int main()
{
    /////////// INITIALISATION //////////
    FLAG2024 = true;
    // if (YEAR == 2025)
    //     Error("Time_Pulse_Calibration is not needed for 2025");

    InitDetectors("Config_Files/sample.pid");
    InitRuns();

    bool Making_Data_FLAG = true;

    // ReadISOLDE();
    // WriteISOLDE();
    // exit(0);

    if (Making_Data_FLAG)
    {
        fTime = MyTFile((DIR_ROOT_DATA_MATCHED + "Time_Pulse_Calibration_" + to_string(YEAR) + ".root").c_str(), "RECREATE");

        InitCalib();
        InitWindows();

        /////////// NUCLEI ANALYSIS //////////
        for (auto pairs : Map_RunFiles)
        {
            string NUCLEUS = pairs.first;
            vector<string> Runs = pairs.second;
            ////// INITIALISATION //////
            Info("Analysis for " + NUCLEUS);

            for (auto Run : Runs)
            {

                Info("Run " + Run, 1);

                string GROUPED_filename = SearchFiles(DIR_ROOT_DATA_GROUPED, Run);
                GROUPED_File = MyTFile((DIR_ROOT_DATA_GROUPED + GROUPED_filename).c_str(), "READ");

                ReadRunDate(Run, GROUPED_File);
                InitHistograms(Run);
                LoadISOLDE(Run);

                TTree *tree = (TTree *)GROUPED_File->Get("CLEANED_Tree");
                TTreeReader *Reader = new TTreeReader(tree);
                TTreeReaderArray<Signal> *Silicon = new TTreeReaderArray<Signal>(*Reader, "CLEANED_Tree_Silicon");
                TTreeReaderValue<vector<vector<pair<Signal, Signal>>>> *SiPM_Groups = new TTreeReaderValue<vector<vector<pair<Signal, Signal>>>>(*Reader, "CLEANED_Tree_SiPMGroup");
                TTreeReaderValue<Signal> *HRS;
                if (YEAR == 2025)
                    HRS = new TTreeReaderValue<Signal>(*Reader, "CLEANED_Tree_HRS");
                

                /// PULSE PERIOD ANALYSIS ///
                vector<double> Peak_Time;
                vector<double> Peak_Intensity;
                vector<int> Peak_Bin;
                clock_t start = clock(), Current;
                if (VERBOSE == 1)
                    Info("Grouping Pulses", 2);
                for (int i = 0; i < H_Data[Run]->GetNbinsX(); i++)
                {
                    ProgressBar(i, H_Data[Run]->GetNbinsX(), start, Current, "Grouping Pulses : ");
                    if (i == 0)
                    {
                        Peak_Time.push_back(0);
                        Peak_Intensity.push_back(0);
                        Peak_Bin.push_back(0);
                    }

                    if (H_Data[Run]->GetBinContent(i) < Intensity_Threshold)
                    {
                        continue;
                    }

                    // cout << "2" << endl;

                    if (Peak_Time.size() == 1)
                    {
                        Peak_Time.push_back(H_Data[Run]->GetBinCenter(i));
                        Peak_Intensity.push_back(H_Data[Run]->GetBinContent(i));
                        Peak_Bin.push_back(i);

                        // find the next peak for first peak
                        for (int j = i; j < H_Data[Run]->GetNbinsX(); j++)
                        {
                            if (H_Data[Run]->GetBinContent(j) < Intensity_Threshold)
                            {
                                Peak_Time.push_back(H_Data[Run]->GetBinCenter(j));
                                Peak_Intensity.push_back(H_Data[Run]->GetBinContent(j));
                                Peak_Bin.push_back(j);
                                i = j + 1;
                                break;
                            }
                        }
                    }

                    else
                    {
                        Peak_Time[2] = H_Data[Run]->GetBinCenter(i);
                        Peak_Intensity[2] = H_Data[Run]->GetBinContent(i);
                        Peak_Bin[2] = i;
                    }

                    // cout << "3" << endl;

                    double diff10 = Peak_Time[1] - Peak_Time[0];
                    double diff21 = Peak_Time[2] - Peak_Time[1];
                    if (diff10 > Neighboor_A_Threshold && diff21 > Neighboor_C_Threshold)
                    {
                        // H_Data_D[Run]->SetBinContent(Peak_Bin[1], H_Data_D[Run]->GetBinContent(Peak_Bin[1]) + Peak_Intensity[1] * 0.4);
                        H_Data[Run]->SetBinError(Peak_Bin[1], 4);
                    }

                    // cout << "4" << endl;

                    // all left & 3.6s right
                    if (diff10 > Neighboor_A_Threshold && diff21 > Neighboor_B_Threshold && diff21 < Neighboor_C_Threshold)
                    {
                        // H_Data_C[Run]->SetBinContent(Peak_Bin[1], H_Data_D[Run]->GetBinContent(Peak_Bin[1]) + Peak_Intensity[1] * 0.6);
                        H_Data[Run]->SetBinError(Peak_Bin[1], 3);
                    }

                    // all left & 2.4s right
                    if (diff21 > Neighboor_A_Threshold && diff21 < Neighboor_B_Threshold)
                    {
                        // H_Data_B[Run]->SetBinContent(Peak_Bin[1], H_Data_D[Run]->GetBinContent(Peak_Bin[1]) + Peak_Intensity[1] * 0.7);
                        H_Data[Run]->SetBinError(Peak_Bin[1], 2);
                    }

                    // all left & 1.2s right
                    if (diff21 < Neighboor_A_Threshold)
                    {
                        // H_Data_A[Run]->SetBinContent(Peak_Bin[1], H_Data_D[Run]->GetBinContent(Peak_Bin[1]) + Peak_Intensity[1] * 0.8);
                        H_Data[Run]->SetBinError(Peak_Bin[1], 1);
                    }

                    // cout << "5" << endl;

                    Peak_Time[0] = Peak_Time[1];
                    Peak_Time[1] = Peak_Time[2];

                    Peak_Intensity[0] = Peak_Intensity[1];
                    Peak_Intensity[1] = Peak_Intensity[2];

                    Peak_Bin[0] = Peak_Bin[1];
                    Peak_Bin[1] = Peak_Bin[2];

                    // cout << "6" << endl;
                }

                ////// DETECTOR ANALYSIS //////
                if (VERBOSE == 1)
                    Info("Detector Analysis", 2);
                int Entries = Reader->GetEntries();
                start = clock(), Current;
                while (Reader->Next())
                {
                    ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "Reading Tree : ");

                    if (YEAR == 2025)
                    {
                        if ((**HRS).isValid)
                        {
                            continue;
                        }
                    }

                    double Time = (*Silicon)[1].Time * 1e-9;
                    int Label = (*Silicon)[1].Label;
                    double Energy = Calibration[Label]->Eval((*Silicon)[1].Channel / 1000);

                    H_Time_All[Run]->Fill(Time);

                    if (Energy > WindowsMap[NUCLEUS][IAS[NUCLEUS]][Label].first && Energy < WindowsMap[NUCLEUS][IAS[NUCLEUS]][Label].second)
                    {
                        H_Exp[Run][Label]->Fill(Energy);
                        // H_Time[Run][Label]->Fill(Time);
                        // H_Time_All[Run]->Fill(Time);

                        if (IsCoincidence(Time, **SiPM_Groups))
                        {
                            H_Exp_Coinc[Run][Label]->Fill(Energy);
                            H_Time_All_Coinc[Run]->Fill(Time);
                        }
                    }
                }

                // Decay
                if (VERBOSE == 1)
                    Info("Decay Analysis", 2);
                int Peak_Group = 0;
                int Peak_Position = 0;
                Entries = H_Data[Run]->GetNbinsX();
                start = clock();
                for (int i = 0; i < H_Data[Run]->GetNbinsX(); i++)
                {
                    ProgressBar(i, Entries, start, Current, "Reading Data : ", 100000);
                    if (H_Data[Run]->GetBinContent(i) > Intensity_Threshold)
                    {
                        // cout << "Peak found at " << H_Data[Run]->GetBinCenter(i) << " with intensity " << H_Data[Run]->GetBinContent(i) << endl;
                        Peak_Group = Group_Determining(i, Run);
                        // cout << "Peak Group: " << Peak_Group << endl;
                        Peak_Position = i;

                        int j = Peak_Position;
                        while (H_Data[Run]->GetBinCenter(j) < H_Data[Run]->GetBinCenter(Peak_Position) + 4 * Pulse_Period)
                        {
                            // H_Time_All_Group[Run][Peak_Group]->SetBinContent(j, H_Time_All[Run]->GetBinContent(j));
                            // H_Decay[Run][Peak_Group]->SetBinContent(j - Peak_Position, H_Decay[Run][Peak_Group]->GetBinContent(j - Peak_Position) + H_Time_All[Run]->GetBinContent(j));
                            H_Decay_Coinc[Run][Peak_Group]->SetBinContent(j - Peak_Position, H_Decay_Coinc[Run][Peak_Group]->GetBinContent(j - Peak_Position) + H_Time_All_Coinc[Run]->GetBinContent(j));
                            H_Decay_All[Run][Peak_Group]->SetBinContent(j - Peak_Position, H_Decay_All[Run][Peak_Group]->GetBinContent(j - Peak_Position) + H_Time_All[Run]->GetBinContent(j));
                            j++;
                        }
                    }
                }
                WriteHistogram(Run, NUCLEUS);
            }
        }

        fTime->Close();
    }

    fTime = MyTFile((DIR_ROOT_DATA_MATCHED + "Time_Pulse_FitParameters_" + to_string(YEAR) + ".root").c_str(), "RECREATE");
    //// CALIBRATION //////
    Start("Calibration");
    TH1D* H_Decay_Coinc_Shifted_MERGED = new TH1D("H_Decay_Coinc_Shifted_MERGED", "H_Decay_Coinc_Shifted_MERGED", 4*Pulse_Period*1000, 0, 4*Pulse_Period);
    H_Decay_Coinc_Shifted_MERGED->GetXaxis()->SetTitle("Time (s)");
    H_Decay_Coinc_Shifted_MERGED->GetYaxis()->SetTitle("Counts / ms");
    // TH1D *H_Decay_All_Shifted_MERGED = new TH1D("H_Decay_All_Shifted_MERGED", "H_Decay_All_Shifted_MERGED", 4*Pulse_Period*1000, 0, 4*Pulse_Period);
    // H_Decay_All_Shifted_MERGED->GetXaxis()->SetTitle("Time (s)");
    // H_Decay_All_Shifted_MERGED->GetYaxis()->SetTitle("Counts / ms");
    string REFERENCE_RUN_str = REFERENCE_RUN > 100 ? to_string(REFERENCE_RUN) : "0" + to_string(REFERENCE_RUN);
    for (auto pairs : Map_RunFiles)
    {
        string NUCLEUS = pairs.first;
        vector<string> Runs = pairs.second;
        // Info("Analysis for " + NUCLEUS);

        for (auto Run : Runs)
        {
            if (stoi(Run) == REFERENCE_RUN)
                continue;

            Info("Calibration for Run " + Run, 1);
            LoadHistograms(Run);
            fTime->cd();
            // Computing Chi2
            pair<double, double> Best_Chi2 = ComputingChi2(H_Decay_Coinc[Run][1], H_Decay_Coinc[REFERENCE_RUN_str][1], Run);
            int bin_tau = Best_Chi2.first/abs(Best_Chi2.first) * H_Decay_Coinc[Run][1]->FindBin(abs(Best_Chi2.first));

            // building solution for coinc
            TH1D *H_Decay_Coinc_Shifted = (TH1D *)H_Decay_Coinc[Run][1]->Clone(("H_Decay_Coinc_Shifted_" + Run).c_str());
            H_Decay_Coinc_Shifted->Reset();
            for (int i = 1; i <= H_Decay_Coinc[Run][1]->GetNbinsX(); i++)
            {
                double bin_center = H_Decay_Coinc[Run][1]->GetBinCenter(i);
                double bin_content = H_Decay_Coinc[Run][1]->GetBinContent(i);
                double bin_error = H_Decay_Coinc[Run][1]->GetBinError(i);

                H_Decay_Coinc_Shifted->SetBinContent(i + bin_tau, bin_content);
                H_Decay_Coinc_Shifted->SetBinError(i + bin_tau, bin_error);
            }

            double offset = -Best_Chi2.first < 0 ? 0 : -Best_Chi2.first;
            H_Decay_Coinc_Shifted->GetXaxis()->SetRangeUser(offset, 1.2+offset);
            H_Decay_Coinc[Run][1]->GetXaxis()->SetRangeUser(offset, 1.2+offset);
            H_Decay_Coinc[REFERENCE_RUN_str][1]->GetXaxis()->SetRangeUser(offset, 1.2+offset);

            H_Decay_Coinc_Shifted->Scale(H_Decay_Coinc[REFERENCE_RUN_str][1]->Integral() / H_Decay_Coinc_Shifted->Integral());
            H_Decay_Coinc[Run][1]->Scale(H_Decay_Coinc[REFERENCE_RUN_str][1]->Integral() / H_Decay_Coinc[Run][1]->Integral());

            H_Decay_Coinc_Shifted->GetXaxis()->SetRangeUser(0, 1.2);
            H_Decay_Coinc[Run][1]->GetXaxis()->SetRangeUser(0, 1.2);
            H_Decay_Coinc[REFERENCE_RUN_str][1]->GetXaxis()->SetRangeUser(0, 1.2);

            // building solution for all    
            // TH1D *H_Decay_All_Shifted = (TH1D *)H_Decay_All[Run][1]->Clone(("H_Decay_All_Shifted_" + Run).c_str());
            // H_Decay_All_Shifted->Reset();
            // for (int i = 1; i <= H_Decay_All[Run][1]->GetNbinsX(); i++)
            // {
            //     double bin_center = H_Decay_All[Run][1]->GetBinCenter(i);
            //     double bin_content = H_Decay_All[Run][1]->GetBinContent(i);
            //     double bin_error = H_Decay_All[Run][1]->GetBinError(i);

            //     H_Decay_All_Shifted->SetBinContent(i + bin_tau, bin_content);
            //     H_Decay_All_Shifted->SetBinError(i + bin_tau, bin_error);
            // }

            // H_Decay_All_Shifted->GetXaxis()->SetRangeUser(offset, 1.2+offset);
            // H_Decay_All[Run][1]->GetXaxis()->SetRangeUser(offset, 1.2+offset);
            // H_Decay_All[REFERENCE_RUN_str][1]->GetXaxis()->SetRangeUser(offset, 1.2+offset);

            // H_Decay_All_Shifted->Scale(H_Decay_All[REFERENCE_RUN_str][1]->Integral() / H_Decay_All_Shifted->Integral());
            // H_Decay_All[Run][1]->Scale(H_Decay_All[REFERENCE_RUN_str][1]->Integral() / H_Decay_All[Run][1]->Integral());

            // H_Decay_All_Shifted->GetXaxis()->SetRangeUser(0, 1.2);
            // H_Decay_All[Run][1]->GetXaxis()->SetRangeUser(0, 1.2);
            // H_Decay_All[REFERENCE_RUN_str][1]->GetXaxis()->SetRangeUser(0, 1.2);

            // plotting
            // dir[Run]->cd();
            TCanvas *c = new TCanvas(("Calibration_" + Run).c_str(), ("Calibration_" + Run).c_str(), 1920, 1080);
            c->Divide(2, 1);
            c->cd(1);
            G_TCalibration[Run]->SetMarkerStyle(20);
            G_TCalibration[Run]->SetTitle(("Time Calibration for " + Run).c_str());
            G_TCalibration[Run]->GetXaxis()->SetTitle("#tau (s)");
            G_TCalibration[Run]->GetYaxis()->SetTitle("#chi^{2}_{#nu}");
            G_TCalibration[Run]->Draw("AP");

            // Determining error bar : Get Y Minimim of TGraphErrors
            double min_y = 1e10;
            double min_x = -1;
            for (int i = 0; i < G_TCalibration[Run]->GetN(); i++)
            {
                double x, y;
                G_TCalibration[Run]->GetPoint(i, x, y);
                if (y < min_y)
                {
                    min_y = y;
                    min_x = x;
                }
            }

            // find roots of the fucntion pol2-(chi2+1)
            gStyle->SetOptStat(0);
            G_TCalibration[Run]->Fit("pol2", "QR", "", min_x - 10 * StepCalibration, min_x + 10 * StepCalibration);
            TF1 *f = G_TCalibration[Run]->GetFunction("pol2");
            TF1 *froot = new TF1("froot", "abs([0] + [1]*x + [2]*x*x)", min_x - 10 * StepCalibration, min_x + 10 * StepCalibration);
            froot->SetParameters(f->GetParameter(0) - (min_y + 1), f->GetParameter(1), f->GetParameter(2));
            double x1 = froot->GetMinimumX(min_x - 10 * StepCalibration, min_x);
            double x2 = froot->GetMinimumX(min_x, min_x + 10 * StepCalibration);
            double tau_error = max(abs(x1 - min_x), abs(x2 - min_x));

            G_TimeShift->AddPoint(stoi(Run), Best_Chi2.first);
            G_TimeShift->SetPointError(G_TimeShift->GetN() - 1, 0, tau_error);

            c->cd(2);
            TLegend *legend = new TLegend(0.5, 0.8, 0.9, 0.9);
            H_Decay_Coinc[REFERENCE_RUN_str][1]->SetTitle(("Coincidence realease curve for run : " + Run).c_str());
            H_Decay_Coinc[REFERENCE_RUN_str][1]->GetXaxis()->SetTitle("Time (s)");
            H_Decay_Coinc[REFERENCE_RUN_str][1]->GetYaxis()->SetTitle("Counts / ms");
            H_Decay_Coinc[REFERENCE_RUN_str][1]->Draw("HIST");
            legend->AddEntry(H_Decay_Coinc[REFERENCE_RUN_str][1], "Reference Run", "l");
            H_Decay_Coinc[REFERENCE_RUN_str][1]->SetLineColor(kRed);
            H_Decay_Coinc[Run][1]->SetLineColor(kBlue);
            H_Decay_Coinc[Run][1]->Draw("HIST SAME");
            legend->AddEntry(H_Decay_Coinc[Run][1], "Before calibration", "l");
            H_Decay_Coinc_Shifted->SetLineColor(kGreen);
            H_Decay_Coinc_Shifted->Draw("HIST SAME");

            ostringstream osschi2;
            osschi2 << fixed << setprecision(2) << Best_Chi2.second;
            string chi2_str = osschi2.str();
            ostringstream osstau;
            osstau << fixed << setprecision(2) << Best_Chi2.first;
            string tau_str = osstau.str();
            ostringstream osstau_error;
            osstau_error << fixed << setprecision(2) << tau_error;
            string tau_error_str = osstau_error.str();
            legend->AddEntry(H_Decay_Coinc_Shifted, ("After calibration (#tau = " + tau_str + " #pm " + tau_error_str + "s)").c_str(), "l");
            legend->Draw("SAME");

            Info("τ = " + tau_str + " ± " + tau_error_str + " s | χ² = " + chi2_str, 2);

            TLatex *latex = new TLatex();
            latex->SetNDC();
            latex->SetTextSize(0.04);
            latex->DrawLatex(0.1, 0.92, ("#chi^{2}_{#nu}: " + chi2_str).c_str());
            latex->Draw("SAME");

            c->Write();
            
            if (NUCLEUS == "32Ar")
            {
                H_Decay_Coinc_Shifted_MERGED->Add(H_Decay_Coinc_Shifted);
                // H_Decay_All_Shifted_MERGED->Add(H_Decay_All_Shifted);
            }
        }
    }
    
    fTime->cd();
    G_TimeShift->SetName("G_Time_Pulse_Calibration");
    G_TimeShift->SetTitle("Time_Pulse_Calibration");
    G_TimeShift->GetXaxis()->SetTitle("Run");
    G_TimeShift->GetYaxis()->SetTitle("Tau (s)");
    G_TimeShift->Write("G_Time_Pulse_Calibration");

    H_Decay_Coinc_Shifted_MERGED->SetTitle("H_Decay_Coinc_Shifted_MERGED");
    H_Decay_Coinc_Shifted_MERGED->GetXaxis()->SetTitle("Time (s)");
    H_Decay_Coinc_Shifted_MERGED->GetYaxis()->SetTitle("Counts / ms");
    H_Decay_Coinc_Shifted_MERGED->Write();

    // H_Decay_All_Shifted_MERGED->SetTitle("H_Decay_All_Shifted_MERGED");
    // H_Decay_All_Shifted_MERGED->GetXaxis()->SetTitle("Time (s)");
    // H_Decay_All_Shifted_MERGED->GetYaxis()->SetTitle("Counts / ms");
    // H_Decay_All_Shifted_MERGED->Write();

    fTime->Close();
}