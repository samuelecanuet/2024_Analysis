#include "Analyser.hh"  

int main()
{
    InitDetectors("Config_Files/sample.pid");
    
    ///////// READING FILE //////////
    string filename = "32Ar";
    int peak = 14;
    MERGED_File = new TFile((DIR_ROOT_DATA_MERGED + filename + "_merged.root").c_str(), "READ");
    // MERGED_File = new TFile((DIR_ROOT_DATA_GROUPED + "run_078_multifast_33Ar_grouped.root").c_str(), "READ");
    TTreeReader *Reader = new TTreeReader((TTree*)MERGED_File->Get("MERGED_Tree"));
    TTreeReaderArray<Signal> *Silicon = new TTreeReaderArray<Signal>(*Reader, "MERGED_Tree_Silicon");
    TTreeReaderArray<Signal> *SiPM_High = new TTreeReaderArray<Signal>(*Reader, "MERGED_Tree_SiPMHigh");
    TTreeReaderArray<Signal> *SiPM_Low = new TTreeReaderArray<Signal>(*Reader, "MERGED_Tree_SiPMLow");
    
    
    CALIBRATED_File = new TFile((DIR_ROOT_DATA_CALIBRATED + "Calibrated.root").c_str(), "READ");    
    InitCalib();

    ///////// NEW FILE //////////
    ANALYSED_File = new TFile((DIR_ROOT_DATA_ANALYSED + filename + "_analysed.root").c_str(), "RECREATE");
    TTree *ANALYSED_Tree = new TTree("ANALYSED_Tree", "ANALYSED_Tree");
    
    double TH[10] = {0, 300, 300, 300, 400, 200, 400, 400, 400, 300};

    InitWindows();
    InitHistograms();
    Info("Start Fake Coincidence Correction");
    clock_t start = clock(), Current;
    int Entries = Reader->GetEntries();
    int limit =5e6;//Entries;
    while (Reader->Next() && Reader->GetCurrentEntry() < limit)
    {
        ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "Reading Tree");
        int Strip_Label = (*Silicon)[1].Label;
        double energy = Calibration[Strip_Label]->Eval((*Silicon)[1].Channel / 1000);
        H_Single[Strip_Label]->Fill(energy);
        // cout << "Energy: " << energy << "    " << Strip_Label << endl;
        double nearest = 1e9;

        // if (WindowsMap[filename][peak][Strip_Label].first < energy && energy < WindowsMap[filename][peak][Strip_Label].second ) // only for ias
        // {
            if ((*SiPM_High).GetSize() > 0)
            {
                // cout << "SiPM_Low: " << (*SiPM_Low).GetSize() << endl;
                for (int i = 0; i < (*SiPM_High).GetSize(); i++)
                {
                    if (abs(nearest) > abs((*SiPM_High)[i].Time))
                    {
                        nearest = (*SiPM_High)[i].Time;
                    }
                }
            }

            for (int mul = 1; mul <= BETA_SIZE; mul++)
            {
                if (mul <= (*SiPM_High).GetSize() && (nearest > start_gate && nearest < end_gate) && TH[GetDetectorChannel((*SiPM_High)[0].Label)] * 1e3 < (*SiPM_High)[0].Channel) // coinc for at least mul
                {

                    H_Coinc[Strip_Label][mul]->Fill(energy);
                    H_SiPM_Time_Coinc[Strip_Label][mul]->Fill(nearest);

                    for (int i = 0; i < (*SiPM_High).GetSize(); i++)
                    {
                        H_SiPM_Channel_M_coinc[(*SiPM_High)[i].Label][mul]->Fill((*SiPM_High)[i].Channel);
                    }
                    for (int i = 0; i < (*SiPM_Low).GetSize(); i++)
                    {
                        H_SiPM_Channel_M_coinc[(*SiPM_Low)[i].Label][mul]->Fill((*SiPM_Low)[i].Channel);
                    }
                }
                else
                {
                    H_NoCoinc[Strip_Label][mul]->Fill(energy);

                    for (int i = 0; i < (*SiPM_High).GetSize(); i++)
                    {
                        H_SiPM_Channel_M_nocoinc[(*SiPM_High)[i].Label][mul]->Fill((*SiPM_High)[i].Channel);
                    }
                    for (int i = 0; i < (*SiPM_Low).GetSize(); i++)
                    {
                        H_SiPM_Channel_M_nocoinc[(*SiPM_Low)[i].Label][mul]->Fill((*SiPM_Low)[i].Channel);
                    }
                }

                if (mul <= (*SiPM_High).GetSize())
                    H_SiPM_Time[Strip_Label][mul]->Fill(nearest);

            }

            if ((*SiPM_High).GetSize() > 0)
            {                   
                if (nearest > start_gate && nearest < end_gate) //coinc
                {
                    H_Coinc_Mulitplicity[Strip_Label][(*SiPM_High).GetSize()]->Fill(energy);
                }
            }
        // }
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
                NFake[i][mul] = 2*HFake[i][mul] * abs(end_gate - start_gate); 
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
                H_NoCoinc_Corrected[i][mul] = (TH1D*)H_NoCoinc[i][mul]->Clone(("H_NoCoinc_Corrected" + detectorName[i] + "_" + to_string(mul)).c_str());
                H_NoCoinc_Corrected[i][mul]->SetTitle(("H_NoCoinc_Corrected" + detectorName[i] + "_" + to_string(mul)).c_str());
                H_Coinc_Corrected[i][mul] = (TH1D*)H_Coinc[i][mul]->Clone(("H_Coinc_Corrected" + detectorName[i] + "_" + to_string(mul)).c_str());
                H_Coinc_Corrected[i][mul]->SetTitle(("H_Coinc_Corrected" + detectorName[i] + "_" + to_string(mul)).c_str());

                //perfect
                // H_NoCoinc[i][mul]->GetXaxis()->SetRangeUser(WindowsMap[filename][peak][i].first, WindowsMap[filename][peak][i].second );
                // H_Coinc[i][mul]->GetXaxis()->SetRangeUser(WindowsMap[filename][peak][i].first, WindowsMap[filename][peak][i].second );
                // double integral_nocoinc = H_NoCoinc[i][mul]->Integral();
                // double integral_coinc = H_Coinc[i][mul]->Integral();
                // double ratio = abs(integral_coinc - integral_nocoinc) / 2 ;

                //compute 
                H_NoCoinc[i][mul]->GetXaxis()->SetRangeUser(WindowsMap[filename][peak][i].first, WindowsMap[filename][peak][i].second );
                double integral = H_NoCoinc[i][mul]->Integral();
                H_NoCoinc[i][mul]->Scale(1. / integral);
                H_NoCoinc[i][mul]->GetXaxis()->SetRangeUser(-1111, -1111);

                

                //apply              
                H_Coinc_Corrected[i][mul]->Add(H_NoCoinc[i][mul], -2*NFake[i][mul]);
                H_NoCoinc_Corrected[i][mul]->Add(H_NoCoinc[i][mul], 2*NFake[i][mul]);

                H_Fake[i][mul] = (TH1D*)H_NoCoinc[i][mul]->Clone(("H_Fake" + detectorName[i] + "_" + to_string(mul)).c_str());
                H_Fake[i][mul]->Scale(NFake[i][mul]);
                H_NoCoinc[i][mul]->Scale(integral);

                //for display
                G_NFake[mul]->AddPoint(i, NFake[i][mul]);
            }
        }
    }

    /// DRAWING
    ANALYSED_File->cd();
    double all_coinc = 0;
    double all_nocoinc = 0;
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            cout << endl;
            Info(detectorName[i]);
            for (int mul = 1; mul <= BETA_SIZE; mul++)
            {
                cout << "Multiplicity: " << mul << endl;
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
                cout << "## RAW ##" << endl;
                int nocoinc = H_NoCoinc[i][mul]->Integral();
                int coinc = H_Coinc[i][mul]->Integral();
                cout << "No Coinc: " << nocoinc << endl;
                cout << "Coinc: " << coinc << endl;
                cout << "Coinc/NoCoinc: " << coinc / (double)nocoinc << endl;
                cout << "## CORRECTED ##" << endl;
                cout << "## Fake: " <<  NFake[i][mul] << endl;
                nocoinc = H_NoCoinc_Corrected[i][mul]->Integral();
                coinc = H_Coinc_Corrected[i][mul]->Integral();
                cout << "No Coinc: " << nocoinc << endl;
                cout << "Coinc: " << coinc << endl;
                cout << "Coinc/NoCoinc: " << coinc / (double)nocoinc << endl;
                cout << "## SHIFT CORRECTED ##" << endl;
                cout << "No Coinc: " << H_NoCoinc_Corrected[i][mul]->GetMean() << endl;
                cout << "Coinc: " << H_Coinc_Corrected[i][mul]->GetMean() << endl;
                cout << "E: " << 0.5*abs(H_NoCoinc_Corrected[i][mul]->GetMean() -  H_Coinc_Corrected[i][mul]->GetMean())<< endl;

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

            //Compare coinc with multiplicity
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
        }
    }

    cout << "All coinc: " << all_coinc << endl;
    cout << "All nocoinc: " << all_nocoinc << endl;
    cout << "ratio: " << all_coinc / all_nocoinc << endl;


    ANALYSED_File->cd();
    //Compare NFake with multiplicity
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
            c->cd();
            for (int m = 1; m < MAX_MULTIPLICTY; m++)
            {
                if (m == 1)
                    H_SiPM_Channel_M_nocoinc[i][m]->Draw("HIST");
                else
                    H_SiPM_Channel_M_nocoinc[i][m]->Draw("HIST SAME");
            }
            c->Write();

            TCanvas *c1 = new TCanvas(("H_SiPM_Channel_coinc" + detectorName[i]).c_str(), ("H_SiPM_Channel_coinc" + detectorName[i]).c_str(), 1920, 1080);
            c1->cd();
            for (int m = 1; m < MAX_MULTIPLICTY; m++)
            {
                if (m == 1)
                    H_SiPM_Channel_M_coinc[i][m]->Draw("HIST");
                else
                    H_SiPM_Channel_M_coinc[i][m]->Draw("HIST SAME");
            }
            c1->Write();
        }
    }
    

    ANALYSED_File->Close();

    return 0;
}