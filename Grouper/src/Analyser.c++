#include "Analyser.hh"

int main(int argc, char **argv)
{
    
    string filename = "32Ar";
    int peak = 14;
    string filename_for_nucleus;
    string addyear = "";

    if (argc == 1)
    {
        Warning("No Run Number Given doing 32Ar in 2024");
    }
    else if (argc >= 2) // 2024 default with run number given
    {
        FLAG2021 = false;
        // Run number
        Run = atoi(argv[1]);
        if (Run < 10)
            Run_string = "00" + to_string(Run);
        else if (Run < 100)
            Run_string = "0" + to_string(Run);
        else
            Run_string = to_string(Run);

        if (argc == 3)
        {
            if (string(argv[2]) == "2021")
            {
                FLAG2021 = true;
                addyear = "/2021/";
            }
            else
            {
                FLAG2021 = false;
            }
        }

        
        filename_for_nucleus = SearchFiles(DIR_ROOT_DATA_GROUPED+addyear, Run_string, "Q");

        if (filename_for_nucleus.find("32Ar"))
        {
            filename = "32Ar";
        }
        else if (filename_for_nucleus.find("33Ar"))
        {
            filename = "33Ar";
            peak = 21;
        }
        else
        {
            Error("Impossible to extract nucleus name in : " + filename_for_nucleus);
        }
    }

    // init detectors with year
    if (FLAG2021)
    {
        InitDetectors("Config_Files/sample.pid", 2021);
    }
    else
    {
        InitDetectors("Config_Files/sample.pid");
    }

    ///////// READING FILE //////////
    TTreeReader *Reader;
    TTreeReaderArray<Signal> *Silicon;
    TTreeReaderValue<vector<vector<pair<Signal, Signal>>>> *SiPM_Groups;
    if (Run == 0)
    {
        MERGED_File = MyTFile((DIR_ROOT_DATA_MERGED + filename + "_merged.root").c_str(), "READ");
        Reader = new TTreeReader((TTree *)MERGED_File->Get("MERGED_Tree"));
        Silicon = new TTreeReaderArray<Signal>(*Reader, "MERGED_Tree_Silicon");
        SiPM_Groups = new TTreeReaderValue<vector<vector<pair<Signal, Signal>>>>(*Reader, "MERGED_Tree_SiPMGroup");
    }
    else
    {
        MERGED_File = MyTFile((DIR_ROOT_DATA_GROUPED + addyear + filename_for_nucleus).c_str(), "READ");
        Reader = new TTreeReader((TTree *)MERGED_File->Get("CLEANED_Tree"));
        Silicon = new TTreeReaderArray<Signal>(*Reader, "CLEANED_Tree_Silicon");
        SiPM_Groups = new TTreeReaderValue<vector<vector<pair<Signal, Signal>>>>(*Reader, "CLEANED_Tree_SiPMGroup");
    }

    CALIBRATED_File = MyTFile((DIR_ROOT_DATA_CALIBRATED + "Calibrated.root").c_str(), "READ");
    InitCalib();

    ///////// NEW FILE //////////
    string f = MERGED_File->GetName();
    string ff = f.substr(0, f.rfind('_'));
    ANALYSED_File = MyTFile((DIR_ROOT_DATA_ANALYSED + addyear + ff.substr(ff.rfind("/") + 1) + "_analysed.root").c_str(), "RECREATE");
    TTree *ANALYSED_Tree = new TTree("ANALYSED_Tree", "ANALYSED_Tree");

    double TH[10] = {0, 300, 300, 300, 400, 200, 400, 400, 400, 300};

    InitWindows();
    InitHistograms();
    Info("Start Fake Coincidence Correction");
    clock_t start = clock(), Current;
    int Entries = Reader->GetEntries();
    int limit = Entries;
    while (Reader->Next() && Reader->GetCurrentEntry() < limit)
    {
        ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "Reading Tree");
        int Strip_Label = (*Silicon)[1].Label;
        double energy = Calibration[Strip_Label]->Eval((*Silicon)[1].Channel / 1000);

        // cout << "Energy: " << energy << "    " << Strip_Label << endl;
        double NEAREST = 1e9;

        if (WindowsMap[filename][peak][Strip_Label].first > energy || energy > WindowsMap[filename][peak][Strip_Label].second) // only for ias
        {
            continue;
        }
        H_Single[Strip_Label]->Fill(energy);
        int NEAREST_GROUP_INDEX = -1;
        vector<double> nearest_group_time;
        vector<double> mean_group_time = vector<double>(10, 0);
        if ((**SiPM_Groups).size() > 0)
        {
            // cout << "Size: " << (**SiPM_Groups).size() << endl;
            // Lopping on subgroups
            for (int i_group = 0; i_group < (**SiPM_Groups).size(); i_group++)
            {
                nearest_group_time.push_back(1e9);
                int counter_mean_group_time = 0;
                // if ((**SiPM_Groups)[i_group].size() > 9)
                // {
                //     for (int i_pair = 0; i_pair < (**SiPM_Groups)[i_group].size(); i_pair++)
                //     {
                //         cout << (**SiPM_Groups)[i_group][i_pair].first << endl;
                //     }
                // }
                // Looping on pair of the subgroup
                for (int i_pair = 0; i_pair < (**SiPM_Groups)[i_group].size(); i_pair++)
                {

                    // Taking valid element of pair with High priority
                    if ((**SiPM_Groups)[i_group][i_pair].first.isValid && !isnan((**SiPM_Groups)[i_group][i_pair].first.Time))
                    {
                        if (abs(nearest_group_time[i_group]) > abs((**SiPM_Groups)[i_group][i_pair].first.Time))
                        {
                            nearest_group_time[i_group] = (**SiPM_Groups)[i_group][i_pair].first.Time;
                        }

                        mean_group_time[i_group] += (**SiPM_Groups)[i_group][i_pair].first.Time;
                        counter_mean_group_time++;
                    }
                    else if ((**SiPM_Groups)[i_group][i_pair].second.isValid && !isnan((**SiPM_Groups)[i_group][i_pair].second.Time))
                    {
                        if (abs(nearest_group_time[i_group]) > abs((**SiPM_Groups)[i_group][i_pair].second.Time))
                        {
                            nearest_group_time[i_group] = (**SiPM_Groups)[i_group][i_pair].second.Time;
                        }
                        mean_group_time[i_group] += (**SiPM_Groups)[i_group][i_pair].first.Time;
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

                for (int mul = 1; mul <= BETA_SIZE; mul++)
                {
                    if (mul <= (**SiPM_Groups)[i_group].size())
                    {
                        H_SiPM_Time[Strip_Label][mul]->Fill(nearest_group_time[i_group]);
                        if (mul == 9)
                        {
                            H_Time_Channel[Strip_Label]->Fill(mean_group_time[i_group], energy);
                        }
                    }
                }
            }
        }

        // cout << "Nearest: " << NEAREST << "    " << NEAREST_GROUP_INDEX << endl;

        int Multiplicity;
        if (NEAREST_GROUP_INDEX == -1)
            Multiplicity = 0;
        else
        {
            Multiplicity = (**SiPM_Groups)[NEAREST_GROUP_INDEX].size();
        }

        // Filling Hist with repect to Multiplicity
        for (int mul = 1; mul <= BETA_SIZE; mul++)
        {
            if (mul <= Multiplicity                              // Mulitplicity sorting condition
                && (NEAREST > start_gate && NEAREST < end_gate)) // Time gate condition
            {

                H_Coinc[Strip_Label][mul]->Fill(energy);
                H_SiPM_Time_Coinc[Strip_Label][mul]->Fill(NEAREST);
                
                if (mul == Multiplicity) // for equal Multiplicity
                    H_Coinc_Mulitplicity[Strip_Label][mul]->Fill(energy);

                for (int i = 0; i < Multiplicity; i++)
                {
                    if (!(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.isValid)
                        continue;
                    H_SiPM_Channel_M_coinc[(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.Label][mul]->Fill((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.Channel);
                }
                for (int i = 0; i < Multiplicity; i++)
                {
                    if (!(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.isValid)
                        continue;
                    H_SiPM_Channel_M_coinc[(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.Label][mul]->Fill((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.Channel);
                }
            }
            else
            {
                H_NoCoinc[Strip_Label][mul]->Fill(energy);

                // for (int i = 0; i < Multiplicity; i++)
                // {
                //     if (!(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.isValid)
                //         continue;
                //     H_SiPM_Channel_M_nocoinc[(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.Label][mul]->Fill((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.Channel);
                // }
                // for (int i = 0; i < Multiplicity; i++)
                // {
                //     if (!(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.isValid)
                //         continue;
                //     H_SiPM_Channel_M_nocoinc[(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.Label][mul]->Fill((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.Channel);
                // }
            }
        }

        for (int i_group = 0; i_group < (**SiPM_Groups).size(); i_group++)
        {
            for (int mul = 1; mul <= BETA_SIZE; mul++)
            {
                if (mul <= (**SiPM_Groups)[i_group].size()                                                                         // Mulitplicity sorting condition
                    && ((**SiPM_Groups)[i_group][0].first.Time < start_gate || (**SiPM_Groups)[i_group][0].first.Time > end_gate)) // Time gate condition
                {
                    for (int i = 0; i < (**SiPM_Groups)[i_group].size(); i++)
                    {
                        if ((**SiPM_Groups)[i_group][i].first.isValid)
                            H_SiPM_Channel_M_nocoinc[(**SiPM_Groups)[i_group][i].first.Label][mul]->Fill((**SiPM_Groups)[i_group][i].first.Channel);
                        if ((**SiPM_Groups)[i_group][i].second.isValid)
                            H_SiPM_Channel_M_nocoinc[(**SiPM_Groups)[i_group][i].second.Label][mul]->Fill((**SiPM_Groups)[i_group][i].second.Channel);
                    }
                }
            }
        }
    }

    ///////// COMPUTE COINCIDENCE CORRECTION //////////
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            for (int mul = 1; mul <= BETA_SIZE; mul++)
            {
                H_SiPM_Time[i][mul]->GetXaxis()->SetRangeUser(start_gate_fake, end_gate_fake);
                HFake[i][mul] = H_SiPM_Time[i][mul]->Integral() / (abs(start_gate_fake - end_gate_fake));
                H_SiPM_Time[i][mul]->GetXaxis()->SetRangeUser(-1111, -1111);
                NFake[i][mul] = HFake[i][mul] * abs(end_gate - start_gate);
            }
        }
    }

    ///////// APPLY CORRECTION OF FAKE COINCIDENCES //////////

    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            for (int mul = 1; mul <= BETA_SIZE; mul++)
            {
                H_NoCoinc_Corrected[i][mul] = (TH1D *)H_NoCoinc[i][mul]->Clone(("H_NoCoinc_Corrected" + detectorName[i] + "_" + to_string(mul)).c_str());
                H_NoCoinc_Corrected[i][mul]->SetTitle(("H_NoCoinc_Corrected" + detectorName[i] + "_" + to_string(mul)).c_str());
                H_Coinc_Corrected[i][mul] = (TH1D *)H_Coinc[i][mul]->Clone(("H_Coinc_Corrected" + detectorName[i] + "_" + to_string(mul)).c_str());
                H_Coinc_Corrected[i][mul]->SetTitle(("H_Coinc_Corrected" + detectorName[i] + "_" + to_string(mul)).c_str());

                // perfect
                //  H_NoCoinc[i][mul]->GetXaxis()->SetRangeUser(WindowsMap[filename][peak][i].first, WindowsMap[filename][peak][i].second );
                //  H_Coinc[i][mul]->GetXaxis()->SetRangeUser(WindowsMap[filename][peak][i].first, WindowsMap[filename][peak][i].second );
                //  double integral_nocoinc = H_NoCoinc[i][mul]->Integral();
                //  double integral_coinc = H_Coinc[i][mul]->Integral();
                //  double ratio = abs(integral_coinc - integral_nocoinc) / 2 ;

                // compute
                H_NoCoinc[i][mul]->GetXaxis()->SetRangeUser(WindowsMap[filename][peak][i].first, WindowsMap[filename][peak][i].second);
                double integral = H_NoCoinc[i][mul]->Integral();
                H_NoCoinc[i][mul]->Scale(1. / integral);
                H_NoCoinc[i][mul]->GetXaxis()->SetRangeUser(-1111, -1111);

                // apply
                H_Coinc_Corrected[i][mul]->Add(H_NoCoinc[i][mul], -NFake[i][mul]);
                H_NoCoinc_Corrected[i][mul]->Add(H_NoCoinc[i][mul], NFake[i][mul]);

                H_Fake[i][mul] = (TH1D *)H_NoCoinc[i][mul]->Clone(("H_Fake" + detectorName[i] + "_" + to_string(mul)).c_str());
                H_Fake[i][mul]->Scale(NFake[i][mul]);
                H_NoCoinc[i][mul]->Scale(integral);

                // for display
                G_NFake[mul]->AddPoint(i, NFake[i][mul]);
            }
        }
    }

    TCanvas *cRatioCoinc_NoCoinc[BETA_SIZE];
    TGraph *G_RatioCoinc_NoCoinc[BETA_SIZE];
    TCanvas *cRatioCoinc_NoCoinc_Corrected[BETA_SIZE];
    TGraph *G_RatioCoinc_NoCoinc_Corrected[BETA_SIZE];
    for (int mul = 1; mul <= BETA_SIZE; mul++)
    {
        cRatioCoinc_NoCoinc[mul] = new TCanvas(("RatioCoinc_NoCoinc_" + to_string(mul)).c_str(), ("RatioCoinc_NoCoinc_" + to_string(mul)).c_str(), 1920, 1080);
        G_RatioCoinc_NoCoinc[mul] = new TGraph();
        cRatioCoinc_NoCoinc_Corrected[mul] = new TCanvas(("RatioCoinc_NoCoinc_Corrected_" + to_string(mul)).c_str(), ("RatioCoinc_NoCoinc_Corrected_" + to_string(mul)).c_str(), 1920, 1080);
        G_RatioCoinc_NoCoinc_Corrected[mul] = new TGraph();
    }

    /// DRAWING
    ANALYSED_File->cd();
    double all_coinc = 0;
    double all_nocoinc = 0;
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            if (H_Single[i]->GetEntries() == 0)
                continue;
            cout << endl;
            Info(detectorName[i]);
            for (int mul = 1; mul <= BETA_SIZE; mul++)
            {
                // Info("Multiplicity: " + to_string(mul), 1);
                dir_FakeCorrection_Strip[i]->cd();
                TCanvas *c = new TCanvas(("RAW_" + detectorName[i] + "_" + to_string(mul)).c_str(), ("RAW_" + detectorName[i] + "_" + to_string(mul)).c_str(), 1920, 1080);
                H_Single[i]->SetLineColor(kBlack);
                H_Single[i]->Draw("HIST");
                H_Coinc[i][mul]->SetLineColor(kRed);
                H_Coinc[i][mul]->GetXaxis()->SetRangeUser(WindowsMap[filename][peak][i].first, WindowsMap[filename][peak][i].second);
                H_Coinc[i][mul]->Draw("HIST SAME");
                H_NoCoinc[i][mul]->SetLineColor(kBlue);
                H_NoCoinc[i][mul]->GetXaxis()->SetRangeUser(WindowsMap[filename][peak][i].first, WindowsMap[filename][peak][i].second);
                H_NoCoinc[i][mul]->Draw("HIST SAME");
                H_Fake[i][mul]->SetLineColor(kGreen);
                H_Fake[i][mul]->GetXaxis()->SetRangeUser(WindowsMap[filename][peak][i].first, WindowsMap[filename][peak][i].second);
                H_Fake[i][mul]->Draw("HIST SAME");
                c->Write();

                TCanvas *c_Corrected = new TCanvas(("Corrected_" + detectorName[i] + "_" + to_string(mul)).c_str(), ("Corrected_" + detectorName[i] + "_" + to_string(mul)).c_str(), 1920, 1080);
                H_Single[i]->SetLineColor(kBlack);
                H_Single[i]->Draw("HIST");
                H_Coinc_Corrected[i][mul]->SetLineColor(kRed);
                H_Coinc_Corrected[i][mul]->GetXaxis()->SetRangeUser(WindowsMap[filename][peak][i].first, WindowsMap[filename][peak][i].second);
                H_Coinc_Corrected[i][mul]->Draw("HIST SAME");
                H_NoCoinc_Corrected[i][mul]->SetLineColor(kBlue);
                H_NoCoinc_Corrected[i][mul]->GetXaxis()->SetRangeUser(WindowsMap[filename][peak][i].first, WindowsMap[filename][peak][i].second);
                H_NoCoinc_Corrected[i][mul]->Draw("HIST SAME");
                c_Corrected->Write();

                TCanvas *c_Time = new TCanvas(("Time_" + detectorName[i] + "_" + to_string(mul)).c_str(), ("Time_" + detectorName[i] + "_" + to_string(mul)).c_str(), 1920, 1080);
                H_SiPM_Time[i][mul]->Draw("HIST");
                H_SiPM_Time_Coinc[i][mul]->SetLineColor(kRed);
                H_SiPM_Time_Coinc[i][mul]->Draw("HIST SAME");
                c_Time->Write();

                // printing //
                
                // Info("## RAW ##", 2);
                int nocoinc = H_NoCoinc[i][mul]->Integral();
                int coinc = H_Coinc[i][mul]->Integral();
                int single = H_Single[i]->Integral();
                if (mul == 9)
                {
                Info("Single: " + to_string(single), 3);
                }
                // Info("No Coinc: " + to_string(nocoinc), 3);
                // Info("Coinc: " + to_string(coinc), 3);
                // Info("Coinc/NoCoinc: " + to_string(coinc / (double)nocoinc), 3);
                G_RatioCoinc_NoCoinc[mul]->AddPoint(i, coinc / (double)nocoinc);
                // Info("## CORRECTED ##", 2);
                // Info("## Fake: " + to_string(NFake[i][mul]), 3);
                nocoinc = H_NoCoinc_Corrected[i][mul]->Integral();
                coinc = H_Coinc_Corrected[i][mul]->Integral();
                // Info("No Coinc: " + to_string(nocoinc), 3);
                // Info("Coinc: " + to_string(coinc), 3);
                // Info("Coinc/NoCoinc: " + to_string(coinc / (double)nocoinc), 3);
                G_RatioCoinc_NoCoinc_Corrected[mul]->AddPoint(i, coinc / (double)nocoinc);
                // Info("## SHIFT CORRECTED ##", 2);
                // Info("No Coinc: " + to_string(H_NoCoinc_Corrected[i][mul]->GetMean()), 3);
                // Info("Coinc: " + to_string(H_Coinc_Corrected[i][mul]->GetMean()), 3);
                // Info("E: " + to_string(0.5 * abs(H_NoCoinc_Corrected[i][mul]->GetMean() - H_Coinc_Corrected[i][mul]->GetMean())), 3);

                dir_FakeCorrection_Strip_Write[i]->cd();
                H_NoCoinc[i][mul]->Write();
                H_Coinc[i][mul]->Write();
                H_NoCoinc_Corrected[i][mul]->Write();
                H_Coinc_Corrected[i][mul]->Write();
                H_SiPM_Time[i][mul]->Write();
                H_SiPM_Time_Coinc[i][mul]->Write();

                // all_coinc+= H_Coinc_Corrected[i][mul]->Integral();
                // all_nocoinc += H_NoCoinc_Corrected[i][mul]->Integral();
            }

            // Compare coinc with multiplicity
            dir_FakeCorrection_Strip[i]->cd();
            TCanvas *c = new TCanvas(("H_Coinc_AtLeastMulitplicity_" + detectorName[i]).c_str(), ("H_Coinc_AtLeastMulitplicity_" + detectorName[i]).c_str(), 1920, 1080);
            TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);
            H_Single[i]->SetLineColor(kBlack);
            H_Single[i]->Draw("HIST");
            for (int mul = 1; mul <= BETA_SIZE; mul++)
            {
                H_Coinc[i][mul]->SetLineColor(mul);
                H_Coinc[i][mul]->Draw("HIST SAME");
                legend->AddEntry(H_Coinc[i][mul], ("Multiplicity " + to_string(mul)).c_str(), "l");
            }
            legend->Draw("SAME");
            c->Write();

            TCanvas *c1 = new TCanvas(("H_Coinc_EqualMulitplicity_" + detectorName[i]).c_str(), ("H_Coinc_EqualMulitplicity_" + detectorName[i]).c_str(), 1920, 1080);
            TLegend *legend1 = new TLegend(0.1, 0.7, 0.48, 0.9);
            H_Single[i]->SetLineColor(kBlack);
            H_Single[i]->Draw("HIST");
            for (int mul = 1; mul <= BETA_SIZE; mul++)
            {
                H_Coinc_Mulitplicity[i][mul]->SetLineColor(mul);
                H_Coinc_Mulitplicity[i][mul]->Draw("HIST SAME");
                legend1->AddEntry(H_Coinc_Mulitplicity[i][mul], ("Multiplicity " + to_string(mul)).c_str(), "l");
            }
            legend1->Draw("SAME");
            c1->Write();

            H_Time_Channel[i]->Write();
        }
    }

    cout << "All coinc: " << all_coinc << endl;
    cout << "All nocoinc: " << all_nocoinc << endl;
    cout << "ratio: " << all_coinc / all_nocoinc << endl;

    ANALYSED_File->cd();
    // Compare NFake with multiplicity
    for (int mul = 1; mul <= BETA_SIZE; mul++)
    {
        TCanvas *c = new TCanvas(("H_NFake_" + to_string(mul)).c_str(), ("H_NFake_" + to_string(mul)).c_str(), 1920, 1080);
        G_NFake[mul]->SetMarkerStyle(20);
        G_NFake[mul]->SetMarkerSize(1);
        G_NFake[mul]->GetXaxis()->SetTitle("Strip");
        G_NFake[mul]->GetYaxis()->SetTitle("Fake coincidences");
        G_NFake[mul]->Draw("AP");
        c->Write();
    }

    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBeta(i))
        {
            TCanvas *c = new TCanvas(("H_SiPM_Channel_nocoinc" + detectorName[i]).c_str(), ("H_SiPM_Channel_nocoinc" + detectorName[i]).c_str(), 1920, 1080);

            TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);
            c->cd();
            for (int m = 1; m <= MAX_MULTIPLICTY; m++)
            {
                H_SiPM_Channel_M_nocoinc[i][m]->SetLineColor(m);
                if (m == 1)
                    H_SiPM_Channel_M_nocoinc[i][m]->Draw("HIST");
                else
                    H_SiPM_Channel_M_nocoinc[i][m]->Draw("HIST SAME");
                legend->AddEntry(H_SiPM_Channel_M_nocoinc[i][m], ("Multiplicity " + to_string(m)).c_str(), "l");
            }
            legend->Draw("SAME");
            c->Write();

            TCanvas *c1 = new TCanvas(("H_SiPM_Channel_coinc" + detectorName[i]).c_str(), ("H_SiPM_Channel_coinc" + detectorName[i]).c_str(), 1920, 1080);
            TLegend *legend1 = new TLegend(0.1, 0.7, 0.48, 0.9);
            c1->cd();
            for (int m = 1; m <= MAX_MULTIPLICTY; m++)
            {
                H_SiPM_Channel_M_coinc[i][m]->SetLineColor(m);
                if (m == 1)
                    H_SiPM_Channel_M_coinc[i][m]->Draw("HIST");
                else
                    H_SiPM_Channel_M_coinc[i][m]->Draw("HIST SAME");
                legend1->AddEntry(H_SiPM_Channel_M_coinc[i][m], ("Multiplicity " + to_string(m)).c_str(), "l");
            }
            legend1->Draw("SAME");
            c1->Write();
        }
    }

    // for each multiplciity
    for (int m = 1; m <= MAX_MULTIPLICTY; m++)
    {
        TCanvas *cH = new TCanvas(("H_SiPMHigh_Channel_coinc_M" + to_string(m)).c_str(), ("H_SiPMHigh_Channel_coinc_M" + to_string(m)).c_str(), 1920, 1080);
        TCanvas *cL = new TCanvas(("H_SiPMLow_Channel_coinc_M" + to_string(m)).c_str(), ("H_SiPMLow_Channel_coinc_M" + to_string(m)).c_str(), 1920, 1080);
        TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorBetaHigh(i))
            {
                cH->cd();
                H_SiPM_Channel_M_coinc[i][m]->SetLineColor(GetDetectorChannel(i));
                if (i == 0)
                    H_SiPM_Channel_M_coinc[i][m]->Draw("HIST");
                else
                    H_SiPM_Channel_M_coinc[i][m]->Draw("HIST SAME");

                legend->AddEntry(H_SiPM_Channel_M_coinc[i][m], ("SiPM" + to_string(GetDetectorChannel(i))).c_str(), "l");
            }

            if (IsDetectorBetaLow(i))
            {
                cL->cd();
                H_SiPM_Channel_M_coinc[i][m]->SetLineColor(GetDetectorChannel(i));
                if (i == 0)
                    H_SiPM_Channel_M_coinc[i][m]->Draw("HIST");
                else
                    H_SiPM_Channel_M_coinc[i][m]->Draw("HIST SAME");
            }
        }
        cH->cd();
        legend->Draw("SAME");
        cH->Write();
        cL->cd();
        legend->Draw("SAME");
        cL->Write();
    }

    int counter = 0;

    for (int mul = 1; mul <= BETA_SIZE; mul++)
    {
        cRatioCoinc_NoCoinc[mul]->cd();
        G_RatioCoinc_NoCoinc[mul]->SetMarkerStyle(20);
        G_RatioCoinc_NoCoinc[mul]->SetMarkerSize(1);
        G_RatioCoinc_NoCoinc[mul]->GetXaxis()->SetTitle("Strip");
        G_RatioCoinc_NoCoinc[mul]->GetYaxis()->SetTitle("Ratio Coinc/NoCoinc");
        G_RatioCoinc_NoCoinc[mul]->Draw("AP");
        cRatioCoinc_NoCoinc[mul]->Write();

        cRatioCoinc_NoCoinc_Corrected[mul]->cd();
        G_RatioCoinc_NoCoinc_Corrected[mul]->SetMarkerStyle(20);
        G_RatioCoinc_NoCoinc_Corrected[mul]->SetMarkerSize(1);
        G_RatioCoinc_NoCoinc_Corrected[mul]->GetXaxis()->SetTitle("Strip");
        G_RatioCoinc_NoCoinc_Corrected[mul]->GetYaxis()->SetTitle("Ratio Coinc/NoCoinc Corrected");
        G_RatioCoinc_NoCoinc_Corrected[mul]->Draw("AP");
        cRatioCoinc_NoCoinc_Corrected[mul]->Write();
    }

    ANALYSED_File->Close();

    return 0;
}