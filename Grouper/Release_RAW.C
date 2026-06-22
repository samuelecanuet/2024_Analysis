#include "include/Detectors.hh"
#include "include/Utilities.hh"
#include <iostream>
#include <cmath>

TF1 *F_Calibration[SIGNAL_MAX];

void InitCalibration()
{
    TFile *f = MyTFile((DIR_ROOT_DATA_CALIBRATED + "Calibrated_"+ to_string(YEAR) + ".root").c_str(), "READ");
    for (int det = 1;  det < SIGNAL_MAX; det++)
    {
        F_Calibration[det] = (TF1*)f->Get(("Calibration_" + detectorName[det]).c_str());
    }
    f->Close();
}


int Release_RAW()
{

    FLAG2025 = true;
    InitDetectors("Config_Files/sample.pid");
    InitWindows();
    InitCalibration();

    TFile *fout = MyTFile("/mnt/hgfs/shared-2/2025_DATA/DETECTOR_DATA/ANALYSED/Release_RAW_pulsecut.root", "RECREATE");

    map<string, TH1D*[SIGNAL_MAX]> H_Energy_RAW;
    map<string, TH1D*> H_Release_RAW;
    map<string, TH1D*> H_Release_RAW_all;
    map<string, TH1D*> H_Release_N;
    double maxi = 3.6;
    double mini = 0;
    H_Release_RAW["32Ar"] = new TH1D("H_Release_RAW_32Ar", "H_Release_RAW_32Ar", maxi*1000, mini, maxi);
    H_Release_RAW["33Ar"] = new TH1D("H_Release_RAW_33Ar", "H_Release_RAW_33Ar", maxi*1000, mini, maxi);
    H_Release_RAW_all["32Ar"] = new TH1D("H_Release_RAW_all_32Ar", "H_Release_RAW_all_32Ar", maxi*1000, mini, maxi);
    H_Release_RAW_all["33Ar"] = new TH1D("H_Release_RAW_all_33Ar", "H_Release_RAW_all_33Ar", maxi*1000, mini, maxi);
    H_Release_N["32Ar"] = new TH1D("H_Release_N_32Ar", "H_Release_N_32Ar", maxi*1000, mini, maxi);
    H_Release_N["33Ar"] = new TH1D("H_Release_N_33Ar", "H_Release_N_33Ar", maxi*1000, mini, maxi);

    map<string, TH1D*> H_HRS_precision;
    H_HRS_precision["32Ar"] = new TH1D("H_HRS_precision_32Ar", "H_HRS_precision", 1000, -0.0002, 0.0002);
    H_HRS_precision["33Ar"] = new TH1D("H_HRS_precision_33Ar", "H_HRS_precision", 1000, -0.0002, 0.0002);

    TH1D *H_Current_Cycle = new TH1D("", "", maxi*1000, mini, maxi);
    TH1D *H_Current_Cycle_all = new TH1D("", "", maxi*1000, mini, maxi);
    map<string, string> NucleiFiles = {
        {"32Ar", "run_000_data_32Ar.root"},
        {"33Ar", "run_066_data_33Ar.root"},
    };

    // map<string, string> NucleiFiles = {
    //     {"32Ar", "../GROUPED/run_000_data_32Ar_grouped_full.root"},
    //     {"33Ar", "../GROUPED/run_066_data_33Ar_grouped_full.root"},
    // };

    for (auto const& [Nucleus, filename] : NucleiFiles)
    {
        H_Current_Cycle->Reset();
        Info("Nucleus : " + Nucleus, 1);
        TFile *ROOT_File = MyTFile((DIR_ROOT_DATA + filename).c_str(), "READ");
        TTree *Tree = (TTree *)ROOT_File->Get("Tree_Group");
        TTreeReader *Reader;
        TTreeReaderArray<Signal> *signals = nullptr;
        TTreeReaderValue<Signal> *HRS = nullptr;
        if (Tree != nullptr)
        {
            Reader = new TTreeReader(Tree);
            signals = new TTreeReaderArray<Signal>(*Reader, "Signal");
        }
        else
        {
            TTree *Tree = (TTree *)ROOT_File->Get("CLEANED_Tree");
            if (Tree != nullptr)
            {
                Reader = new TTreeReader(Tree);
                signals = new TTreeReaderArray<Signal>(*Reader, "CLEANED_Tree_Silicon");
                HRS = new TTreeReaderValue<Signal>(*Reader, "CLEANED_Tree_HRS");
            }
        }
        

        // Read data and fill release histograms
        clock_t start = clock(), Current;   
        double Pulse_Time = DBL_MAX;
        int Entries = Reader->GetEntries();
        while (Reader->Next())
        {
            ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "Filling Release Histograms : ");
            
            for (int index = 0; index < signals->GetSize(); index++)
            {
                int current_label = (*signals)[index].Label;
                if (IsDetectorSiliStrip(current_label))
                {
                    // H_Channel_RAW[current_label]->Fill(signals[index].Channel);
                    (*signals)[index].Channel = (*signals)[index].Channel - 0.5 + gRandom->Rndm();

                    (*signals)[index].Time = (*signals)[index].Time + (-0.5 + gRandom->Rndm()) * 8;
                    double Energy = F_Calibration[current_label]->Eval((*signals)[index].Channel/1000);

                    if (Energy > WindowsMap[Nucleus][IAS[Nucleus]][current_label].first && Energy < WindowsMap[Nucleus][IAS[Nucleus]][current_label].second)
                    {
                        // H_Energy_RAW[Nucleus][current_label]->Fill(Energy);
                        // H_Release_RAW[Nucleus]->Fill((signals[index].Time - Pulse_Time)*1e-9);
                        H_Current_Cycle->Fill(((*signals)[index].Time - Pulse_Time)*1e-9);
                    }
                    H_Current_Cycle_all->Fill(((*signals)[index].Time - Pulse_Time)*1e-9);
                }

                else if (IsHRSProton(current_label))
                {
                    if (H_Current_Cycle != nullptr)
                    {
                        if (((*signals)[index].Time - Pulse_Time)*1e-9 > 2.6)
                        {                                
                            H_Release_RAW[Nucleus]->Add(H_Current_Cycle);
                            H_Release_RAW_all[Nucleus]->Add(H_Current_Cycle_all);   
                        }
                        H_HRS_precision[Nucleus]->Fill((double)std::remainder(((*signals)[index].Time - Pulse_Time)*1e-9, 1.20));
                        H_Current_Cycle->Reset();
                        H_Current_Cycle_all->Reset();
                    }
                    else
                    {
                        H_Current_Cycle = new TH1D(("H_Current_Cycle_" + Nucleus).c_str(), ("H_Current_Cycle_" + Nucleus).c_str(), maxi*1000, mini, maxi);
                    }

                    Pulse_Time = (*signals)[index].Time;
                    H_Release_N[Nucleus]->Fill(((*signals)[index].Time - Pulse_Time)*1e-9);
                }
            }

            if (HRS != nullptr)
                {
                    if ((*HRS)->isValid)
                    {
                        if (H_Current_Cycle != nullptr)
                        {
                            if (((*HRS)->Time - Pulse_Time)*1e-9 > 2.6)
                            {                                
                                H_Release_RAW[Nucleus]->Add(H_Current_Cycle);
                                H_Release_RAW_all[Nucleus]->Add(H_Current_Cycle_all);   
                            }
                            H_HRS_precision[Nucleus]->Fill((double)std::remainder(((*HRS)->Time - Pulse_Time)*1e-9, 1.20));
                            H_Current_Cycle->Reset();
                            H_Current_Cycle_all->Reset();
                        }
                        else
                        {
                            H_Current_Cycle = new TH1D(("H_Current_Cycle_" + Nucleus).c_str(), ("H_Current_Cycle_" + Nucleus).c_str(), maxi*1000, mini, maxi);
                            H_Current_Cycle_all = new TH1D(("H_Current_Cycle_all_" + Nucleus).c_str(), ("H_Current_Cycle_all_" + Nucleus).c_str(), maxi*1000, mini, maxi);
                        }

                        Pulse_Time = (*HRS)->Time;
                        H_Release_N[Nucleus]->Fill(((*HRS)->Time - Pulse_Time)*1e-9);
                    }
                }
        }

        ROOT_File->Close();

        fout->cd();
        // for (int det = 1; det < SIGNAL_MAX; det++)
        // {
        //     if (IsDetectorSiliStrip(det))
        //     {
        //         H_Energy_RAW[Nucleus][det]->Write();
        //     }
        // }
        H_Release_RAW[Nucleus]->Write();
        H_Release_RAW_all[Nucleus]->Write();
        H_HRS_precision[Nucleus]->Write();
        H_Release_N[Nucleus]->Write();
    }
    fout->Close();

    return 0;
}