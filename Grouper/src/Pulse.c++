#include "Pulse.hh"

int main()
{
    /////////// INITIALISATION //////////
    InitDetectors("Config_Files/sample.pid");
    ReadAllRunsDate();

    // ReadISOLDE();
    // WriteISOLDE();

    

    fTime = new TFile((DIR_ROOT_DATA_GROUPED + "../Pulses.root").c_str(), "RECREATE");
    Init();
    InitCalib();
    InitWindows();

    /////////// NUCLEI ANALYSIS //////////
    for (auto NUCLEUS : NUCLEI_LIST)
    {
        ////// INITIALISATION //////
        Info("Analysis for " + NUCLEUS);

        for (auto Run : MAP_Run_Number[NUCLEUS])
        {
            Info("Run " + to_string(Run));

            LoadISOLDE(64);

            string GROUPED_filename = SearchFiles(DIR_ROOT_DATA_GROUPED, "0"+to_string(Run));
            GROUPED_File = MyTFile((DIR_ROOT_DATA_GROUPED + GROUPED_filename).c_str(), "READ");

            TTree *tree = (TTree *)GROUPED_File->Get("CLEANED_Tree");
            TTreeReader *Reader = new TTreeReader(tree);
            TTreeReaderArray<Signal> *Silicon = new TTreeReaderArray<Signal>(*Reader, "CLEANED_Tree_Silicon");
            TTreeReaderValue<vector<vector<pair<Signal, Signal>>>> *SiPM = new TTreeReaderValue<vector<vector<pair<Signal, Signal>>>>(*Reader, "CLEANED_Tree_SiPMGroup");

            /// PULSE PERIOD ANALYSIS ///

            vector<double> Peak_Time;
            vector<double> Peak_Intensity;
            vector<int> Peak_Bin;
            for (int i = 0; i < H_Data[Run]->GetNbinsX(); i++)
            {
                ProgressBar(i, H_Data[Run]->GetNbinsX(), 0, 0, "Grouping Pulses : ");
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

                double diff10 = Peak_Time[1] - Peak_Time[0];
                double diff21 = Peak_Time[2] - Peak_Time[1];
                if (diff10 > Neighboor_A_Threshold && diff21 > Neighboor_C_Threshold)
                {
                    H_Data_D[NUCLEUS]->SetBinContent(Peak_Bin[1], H_Data_D[NUCLEUS]->GetBinContent(Peak_Bin[1]) + Peak_Intensity[1] * 0.4);
                }

                // all left & 3.6s right
                if (diff10 > Neighboor_A_Threshold && diff21 > Neighboor_B_Threshold && diff21 < Neighboor_C_Threshold)
                {
                    H_Data_C[NUCLEUS]->SetBinContent(Peak_Bin[1], H_Data_D[NUCLEUS]->GetBinContent(Peak_Bin[1]) + Peak_Intensity[1] * 0.6);
                }

                // all left & 2.4s right
                if (diff21 > Neighboor_A_Threshold && diff21 < Neighboor_B_Threshold)
                {
                    H_Data_B[NUCLEUS]->SetBinContent(Peak_Bin[1], H_Data_D[NUCLEUS]->GetBinContent(Peak_Bin[1]) + Peak_Intensity[1] * 0.7);
                }

                // all left & 1.2s right
                if (diff21 < Neighboor_A_Threshold)
                {
                    H_Data_A[NUCLEUS]->SetBinContent(Peak_Bin[1], H_Data_D[NUCLEUS]->GetBinContent(Peak_Bin[1]) + Peak_Intensity[1] * 0.8);
                }

                Peak_Time[0] = Peak_Time[1];
                Peak_Time[1] = Peak_Time[2];

                Peak_Intensity[0] = Peak_Intensity[1];
                Peak_Intensity[1] = Peak_Intensity[2];

                Peak_Bin[0] = Peak_Bin[1];
                Peak_Bin[1] = Peak_Bin[2];
            }

            ////// DETECTOR ANALYSIS //////
            if (NUCLEUS == "18N")
            {
                while (Reader->Next())
                {
                    ProgressBar(Reader->GetCurrentEntry(), Reader->GetEntries(), 0, 0, "Reading Tree : ");
                    double Time = (*Silicon)[1].Time * 1e-9;
                    int Label = (*Silicon)[1].Label;
                    double Energy = Calibration[Label]->Eval((*Silicon)[1].Channel);

                    if ((**SiPM).size() < 9)
                        H_Channel18N[Label]->Fill(Energy);
                    else
                        H_Channel18N_coinc[Label]->Fill(Energy);

                    if ((Energy > peaks_window_F[Label][0].first && Energy < peaks_window_F[Label][0].second) || (Energy > peaks_window_F[Label][1].first && Energy < peaks_window_F[Label][1].second))
                    {
                        if ((**SiPM).size() < 9)
                            continue;
                        H_Exp[NUCLEUS][Label]->Fill(Energy);
                        H_Time[NUCLEUS][Label]->Fill(Time);
                        H_Time_All[NUCLEUS]->Fill(Time);
                    }

                    if ((*Silicon).GetSize() == 3 && Label < 50)
                    {
                        H_Channel18N_3a[Label]->Fill(Energy);
                        H_Time_All_3a->Fill(Time);
                    }
                }
            }
            else
            {   
                int peak;
                if (NUCLEUS == "32Ar")
                {peak=14;}
                else if (NUCLEUS == "33Ar")
                {peak=21;}
                else
                {peak=0;}
                while (Reader->Next())
                {
                    ProgressBar(Reader->GetCurrentEntry(), Reader->GetEntries(), 0, 0, "Reading Tree : ");
                    double Time = (*Silicon)[1].Time * 1e-9;
                    int Label = (*Silicon)[1].Label;
                    double Energy = Calibration[Label]->Eval((*Silicon)[1].Channel/1000);
                    
                    if (Energy > WindowsMap[NUCLEUS][peak][Label].first && Energy < WindowsMap[NUCLEUS][peak][Label].second)
                    {
                        H_Exp[NUCLEUS][Label]->Fill(Energy);
                        H_Time[NUCLEUS][Label]->Fill(Time);
                        H_Time_All[NUCLEUS]->Fill(Time);
                        
                        if (IsCoincidence((**SiPM)))
                        {
                            H_Exp_Coinc[NUCLEUS][Label]->Fill(Energy);
                            H_Time_All_Coinc[NUCLEUS]->Fill(Time);
                        }
                    }
                }
            }

            // Decay
            Info("Decay Analysis");
            int Peak_Group = 0;
            int Peak_Position = 0;
            int Entries = H_Data[Run]->GetNbinsX();
            for (int i = 0; i < H_Data[Run]->GetNbinsX(); i++)
            {
                ProgressCounter(i, Entries, "Reading Data", 100000);
                if (H_Data[Run]->GetBinContent(i) > Intensity_Threshold)
                {
                    Peak_Group = Group_Determining(i, NUCLEUS);
                    Peak_Position = i + H_Data[Run]->FindBin(MAP_Delta_Time[NUCLEUS][0]);

                    int j = Peak_Position;
                    while (H_Data[Run]->GetBinCenter(j) < H_Data[Run]->GetBinCenter(Peak_Position) + 4 * Pulse_Period)
                    {
                        H_Time_All_Group[NUCLEUS][Peak_Group]->SetBinContent(j - H_Data[Run]->FindBin(MAP_Delta_Time[NUCLEUS][0]), H_Time_All[NUCLEUS]->GetBinContent(j - H_Data[Run]->FindBin(MAP_Delta_Time[NUCLEUS][0])));
                        H_Decay[NUCLEUS][Peak_Group]->SetBinContent(j - Peak_Position, H_Decay[NUCLEUS][Peak_Group]->GetBinContent(j - Peak_Position) + H_Time_All[NUCLEUS]->GetBinContent(j));
                        H_Decay_Coinc[NUCLEUS][Peak_Group]->SetBinContent(j - Peak_Position, H_Decay_Coinc[NUCLEUS][Peak_Group]->GetBinContent(j - Peak_Position) + H_Time_All_Coinc[NUCLEUS]->GetBinContent(j));
                        j++;
                    }
                }
            }
        }
        WriteHistogram(NUCLEUS);
    }

    fTime->Close();
}