#include "include/Detectors.hh"
#include "include/Utilities.hh"
#include <iostream>
#include <cmath>

TF1 *Calibration[SIGNAL_MAX];

void InitCalib(int sigma_calib = 0)
{
    TFile *CALIBRATED_File = MyTFile((DIR_ROOT_DATA_CALIBRATED + "Calibrated_" + to_string(YEAR) + ".root").c_str(), "READ");
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            Calibration[i] = (TF1 *)CALIBRATED_File->Get(("Calibration_" + detectorName[i]).c_str());
            

            if (Calibration[i] == NULL)
            {
                Error("No calibration found for " + detectorName[i]);
            }

            // error progagtion
            Calibration[i]->SetParameter(1, Calibration[i]->GetParameter(1) + sigma_calib * Calibration[i]->GetParError(1));
            Calibration[i]->SetParameter(0, Calibration[i]->GetParameter(0) + sigma_calib * Calibration[i]->GetParError(0));            
        }
    }
    Success("Silicon Calibration loaded");
}

int PSB_Synchro()
{
    FLAG2025 = true;
    InitDetectors("Config_Files/sample.pid");
    InitCalib();
    InitWindows();

    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliStrip(det))
        {
            WindowsMap["18N"][0][det] = make_pair(0, 10e3);
        }
    }
    
    int N = 4;

    vector<string> Nuclei = {"32Ar", "33Ar"};
    map<string, TH1D *[10]> H_Release;
    map<string, TH1D *[10]> H_Release_Full;
    map<string, TH1D*> H_Single;

    TFile *RESULT_File = MyTFile((DIR_ROOT_DATA_ANALYSED + "../PSB_Synchro.root").c_str(), "RECREATE");

    for (const auto &Nucleus : Nuclei)
    {
        Info("Reading " + Nucleus);
        TFile *f = MyTFile((DIR_ROOT_DATA_MERGED + Nucleus + "_merged.root").c_str(), "READ");
        TTree *Tree = (TTree*)f->Get("MERGED_Tree");
        if (Tree == NULL)        {
            Warning("No Tree found in " + Nucleus + "_merged.root");
            continue;
        }

        TTreeReader *Reader = new TTreeReader(Tree);
        TTreeReaderArray<Signal> *Silicon = new TTreeReaderArray<Signal>(*Reader, "MERGED_Tree_Silicon");
        TTreeReaderValue<vector<vector<pair<Signal, Signal>>>> *SiPM_Groups = new TTreeReaderValue<vector<vector<pair<Signal, Signal>>>>(*Reader, "MERGED_Tree_SiPMGroup"); 
        TTreeReaderValue<Signal> *HRS = new TTreeReaderValue<Signal>(*Reader, "MERGED_Tree_HRS");

        // InitHIst
        for (int i = 0; i < 10; i++)
        {
                H_Release[Nucleus][i] = new TH1D(("H_Release_T" + to_string(i)).c_str(), ("H_Release_T" + to_string(i)).c_str(), 120*N, 0, 1.20*N);
                H_Release[Nucleus][i]->GetXaxis()->SetTitle("Time (s)");
                H_Release[Nucleus][i]->GetYaxis()->SetTitle("Counts / 10 ms");
                H_Release[Nucleus][i]->GetXaxis()->CenterTitle();
                H_Release[Nucleus][i]->GetYaxis()->CenterTitle();

                H_Release_Full[Nucleus][i] = new TH1D(("H_Release_Full_T" + to_string(i)).c_str(), ("H_Release_Full_T" + to_string(i)).c_str(), 120*(N+2), 0, 1.20*(N+2));
                H_Release_Full[Nucleus][i]->GetXaxis()->SetTitle("Time (s)");
                H_Release_Full[Nucleus][i]->GetYaxis()->SetTitle("Counts / 10 ms");
                H_Release_Full[Nucleus][i]->GetXaxis()->CenterTitle();
                H_Release_Full[Nucleus][i]->GetYaxis()->CenterTitle();
        }

        H_Single[Nucleus] = new TH1D(("H_Single_" + Nucleus).c_str(), ("H_Single_" + Nucleus).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
        H_Single[Nucleus]->GetXaxis()->SetTitle("Energy (keV)");
        H_Single[Nucleus]->GetYaxis()->SetTitle("Counts");
        H_Single[Nucleus]->GetXaxis()->CenterTitle();
        H_Single[Nucleus]->GetYaxis()->CenterTitle();
    
        clock_t start = clock(), Current;   
        int Entries = Reader->GetEntries();
        double TimeHRS = 0;
        double LastTimeHRS = 0;
        double FullTimeHRS = 0;
        while (Reader->Next() && Reader->GetCurrentEntry() < Entries)
        {
            ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "Merging : ");         
            // Info("Current Entry : " + to_string(Reader->GetCurrentEntry()) + " / " + to_string(Entries), 2);
            if ((**HRS).isValid)
            {
                LastTimeHRS = TimeHRS;
                if (TimeHRS > (**HRS).Time*1e-9 - 1.3)
                    H_Release[Nucleus][1]->Add(H_Release[Nucleus][0], 1);
                else if (TimeHRS > (**HRS).Time*1e-9 - 2.5)
                    H_Release[Nucleus][2]->Add(H_Release[Nucleus][0], 1);
                else if (TimeHRS > (**HRS).Time*1e-9 - 3.7)
                    H_Release[Nucleus][3]->Add(H_Release[Nucleus][0], 1);
                else if (TimeHRS > (**HRS).Time*1e-9 - 4.9)
                    H_Release[Nucleus][4]->Add(H_Release[Nucleus][0], 1);
                else if (TimeHRS > (**HRS).Time*1e-9 - 6.2)
                    H_Release[Nucleus][5]->Add(H_Release[Nucleus][0], 1);
                TimeHRS = (**HRS).Time*1e-9;               
                H_Release[Nucleus][0]->Reset();

                if (TimeHRS-FullTimeHRS > 6.2)
                {
                    FullTimeHRS = TimeHRS;
                }
                continue;
            }
            else 
            {
                // Info("Silicon event", 2);
                int det = (*Silicon)[1].Label;
                if (det > 50) continue; // rear strips only
                double energy = Calibration[det]->Eval((*Silicon)[1].Channel / 1000.);
                double time = (*Silicon)[1].Time*1e-9;
                
                H_Single[Nucleus]->Fill(energy);
                // Info("Energy : " + to_string(energy) + " keV, Time : " + to_string(time - TimeHRS) + " ns", 3);

                // energy gate 
                if (energy > WindowsMap[Nucleus][IAS[Nucleus]][det].first && energy < WindowsMap[Nucleus][IAS[Nucleus]][det].second)
                {
                    // cout << TimeHRS-LastTimeHRS << endl;
                    // Info("Accepted event", 3);
                    H_Release[Nucleus][0]->Fill((time - TimeHRS));
                    
                    if (TimeHRS-LastTimeHRS < 1.3)
                        H_Release_Full[Nucleus][1]->Fill((time - FullTimeHRS));
                    else if (TimeHRS-LastTimeHRS < 2.5)
                        H_Release_Full[Nucleus][2]->Fill((time - FullTimeHRS));
                    else if (TimeHRS-LastTimeHRS < 3.7)
                        H_Release_Full[Nucleus][3]->Fill((time - FullTimeHRS));
                    else if (TimeHRS-LastTimeHRS < 4.9)
                        H_Release_Full[Nucleus][4]->Fill((time - FullTimeHRS));
                    else if (TimeHRS-LastTimeHRS < 6.2)
                        H_Release_Full[Nucleus][5]->Fill((time - FullTimeHRS));
                }
            }
        }
    }

    RESULT_File->cd();
    // Energy

    TCanvas *c = new TCanvas("c_Energy", "c_Energy", 800, 600);
    int count = 0;
    for (const auto &Nucleus : Nuclei)
    {
        H_Single[Nucleus]->Scale(1./H_Single[Nucleus]->Integral());
        H_Single[Nucleus]->SetLineColor(kRed);
        H_Single[Nucleus]->Draw(count==0 ? "HIST" : "HIST SAME");
        count++;
    }
    c->BuildLegend();
    c->Write();

    // Time
    for (const auto &Nucleus : Nuclei)
    {
        TCanvas *c_Time = new TCanvas(Form("c_Time_%s", Nucleus.c_str()), Form("c_Time_%s", Nucleus.c_str()), 800, 600);
        for (int i = 1; i < 6; i++)
        {
            // H_Release[Nucleus][i]->Scale(1./H_Release[Nucleus][i]->Integral());
            H_Release[Nucleus][i]->SetLineColor(i+1);
            H_Release[Nucleus][i]->Draw(i==0 ? "HIST" : "HIST SAME");
        }
        c_Time->BuildLegend();
        c_Time->Write();

        TCanvas *c_Time_Full = new TCanvas(Form("c_Time_Full_%s", Nucleus.c_str()), Form("c_Time_Full_%s", Nucleus.c_str()), 800, 600);
        for (int i = 1; i < 6; i++)        {
            // H_Release_Full[Nucleus][i]->Scale(1./H_Release_Full[Nucleus][i]->Integral());
            H_Release_Full[Nucleus][i]->SetLineColor(i+1);
            H_Release_Full[Nucleus][i]->Draw(i==0 ? "HIST" : "HIST SAME");
         }
        c_Time_Full->BuildLegend();
        c_Time_Full->Write();
    }
    


    RESULT_File->Close();
    return 0;
}