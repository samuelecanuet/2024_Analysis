#include "LifeTime.hh"

int main()
{
    /////////// INITIALISATION //////////
    InitDetectors("Config_Files/sample.pid");
    ReadAllRunsDate();
    ReadISOLDE();

    fTime = new TFile((DIR_ROOT_DATA_GROUPED + "../Time.root").c_str(), "RECREATE");
    Init();

    /////////// NUCLEI ANALYSIS //////////
    for (auto NUCLEUS : NUCLEI_LIST)
    {

        ////// INITIALISATION //////
        cout << endl;
        Info("Analysis for " + NUCLEUS);
        InitPeakWindow(NUCLEUS);

        for (int run_i = 0; run_i < MAP_Run_Number[NUCLEUS].size(); run_i++)
        {
            Info("Run " + to_string(MAP_Run_Number[NUCLEUS][run_i]));
            run_number = MAP_Run_Number[NUCLEUS][run_i];

            if (run_number < 100)
                run_number_str = "0" + to_string(run_number);
            else
                run_number_str = to_string(run_number);

            GROUPED_File = new TFile((DIR_ROOT_DATA_GROUPED + MAP_Filename[NUCLEUS][run_i]).c_str(), "READ");

            if (!GROUPED_File->IsOpen())
                Error("Impossible to open ROOT file");

            TTree *tree = (TTree *)GROUPED_File->Get("CLEANED_Tree");
            TTreeReader *Reader = new TTreeReader(tree);
            TTreeReaderArray<Signal> *Silicon = new TTreeReaderArray<Signal>(*Reader, "CLEANED_Tree_Silicon");
            TTreeReaderArray<Signal> *SiPMHigh = new TTreeReaderArray<Signal>(*Reader, "CLEANED_Tree_SiPMHigh");

            /// PULSE PERIOD ANALYSIS ///

            vector<double> Peak_Time;
            vector<double> Peak_Intensity;
            vector<int> Peak_Bin;
            for (int i = 0; i < H_Data[run_number]->GetNbinsX(); i++)
            {
                ProgressBar(i, H_Data[run_number]->GetNbinsX(), 0, 0, "Grouping Pulses : ");
                if (i == 0)
                {
                    Peak_Time.push_back(0);
                    Peak_Intensity.push_back(0);
                    Peak_Bin.push_back(0);
                }

                if (H_Data[run_number]->GetBinContent(i) > Intensity_Threshold)
                {
                    if (Peak_Time.size() == 1)
                    {
                        Peak_Time.push_back(H_Data[run_number]->GetBinCenter(i));
                        Peak_Intensity.push_back(H_Data[run_number]->GetBinContent(i));
                        Peak_Bin.push_back(i);

                        // find the next peak for first peak
                        for (int j = i; j < H_Data[run_number]->GetNbinsX(); j++)
                        {
                            if (H_Data[run_number]->GetBinContent(j) < Intensity_Threshold)
                            {
                                Peak_Time.push_back(H_Data[run_number]->GetBinCenter(j));
                                Peak_Intensity.push_back(H_Data[run_number]->GetBinContent(j));
                                Peak_Bin.push_back(j);
                                i = j + 1;
                                break;
                            }
                        }
                    }

                    else
                    {
                        Peak_Time[2] = H_Data[run_number]->GetBinCenter(i);
                        Peak_Intensity[2] = H_Data[run_number]->GetBinContent(i);
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
            }

            ////// DETECTOR ANALYSIS //////
            if (NUCLEUS == "18N")
            {
                while (Reader->Next())
                {
                    ProgressBar(Reader->GetCurrentEntry(), Reader->GetEntries(), 0, 0, "Reading Tree : ");
                    double Time = (*Silicon)[1].Time * 1e-9;
                    double Channel = (*Silicon)[1].Channel;
                    int Label = (*Silicon)[1].Label;

                    if ((*SiPMHigh).GetSize() < 9)
                        H_Channel18N[Label]->Fill(Channel);  
                    else
                        H_Channel18N_coinc[Label]->Fill(Channel);

                    if ((Channel > peaks_window_F[Label][0].first && Channel < peaks_window_F[Label][0].second) || (Channel > peaks_window_F[Label][1].first && Channel < peaks_window_F[Label][1].second))
                    {
                        if ((*SiPMHigh).GetSize() < 9)
                            continue;
                        H_Exp[NUCLEUS][Label]->Fill(Channel);
                        H_Time[NUCLEUS][Label]->Fill(Time);
                        H_Time_All[NUCLEUS]->Fill(Time);
                    }

                    if ((*Silicon).GetSize() == 3 && Label < 50)
                    {
                        H_Channel18N_3a[Label]->Fill(Channel);  
                        H_Time_All_3a->Fill(Time);
                    }

                                      
                }
            }
            else
            {
                while (Reader->Next())
                {
                    ProgressBar(Reader->GetCurrentEntry(), Reader->GetEntries(), 0, 0, "Reading Tree : ");
                    double Time = (*Silicon)[1].Time * 1e-9;
                    double Channel = (*Silicon)[1].Channel;
                    int Label = (*Silicon)[1].Label;

                    if (Channel > peaks_window_F[Label][0].first && Channel < peaks_window_F[Label][0].second)
                    {
                        H_Exp[NUCLEUS][Label]->Fill(Channel);
                        H_Time[NUCLEUS][Label]->Fill(Time);
                        H_Time_All[NUCLEUS]->Fill(Time);
                    }
                }
            }

            // Decay
            int Peak_Group = 0;
            int Peak_Position = 0;
            for (int i = 0; i < H_Data[run_number]->GetNbinsX(); i++)
            {
                if (H_Data[run_number]->GetBinContent(i) > Intensity_Threshold)
                {
                    Peak_Group = Group_Determining(i, NUCLEUS);
                    Peak_Position = i + H_Data[run_number]->FindBin(MAP_Delta_Time[NUCLEUS][run_i]);

                    int j = Peak_Position;
                    while (H_Data[run_number]->GetBinCenter(j) < H_Data[run_number]->GetBinCenter(Peak_Position) + 4 * Pulse_Period)
                    {
                        H_Time_All_Group[NUCLEUS][Peak_Group]->SetBinContent(j - H_Data[run_number]->FindBin(MAP_Delta_Time[NUCLEUS][run_i]), H_Time_All[NUCLEUS]->GetBinContent(j - H_Data[run_number]->FindBin(MAP_Delta_Time[NUCLEUS][run_i])));
                        H_Decay[NUCLEUS][Peak_Group]->SetBinContent(j - Peak_Position, H_Decay[NUCLEUS][Peak_Group]->GetBinContent(j - Peak_Position) + H_Time_All[NUCLEUS]->GetBinContent(j));
                        j++;
                    }
                }
            }
        }

        // Fitting Decay only for case C
        int group = 3;
        if (NUCLEUS == "18N")
        {
            group = 4;
            H_Decay[NUCLEUS][group]->Rebin(20);
        }
        else
        {
            H_Decay[NUCLEUS][group]->Rebin(5);
        }
            

        ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
        fReleaseCurve[NUCLEUS] = new TF1(("fReleaseCurve_" + NUCLEUS).c_str(), ReleaseCurve, 0.01, group * Pulse_Period, 6);
        fReleaseCurve[NUCLEUS]->SetNpx(1000);
        fReleaseCurve[NUCLEUS]->SetParLimits(0, 0, 1);
        fReleaseCurve[NUCLEUS]->FixParameter(1, 0);
        fReleaseCurve[NUCLEUS]->SetParLimits(2, 0, 100);
        fReleaseCurve[NUCLEUS]->FixParameter(3, 0);
        fReleaseCurve[NUCLEUS]->SetParLimits(4, 0, 1);
        fReleaseCurve[NUCLEUS]->SetParameter(4, MAP_LifeTime[NUCLEUS]);
        fReleaseCurve[NUCLEUS]->SetParLimits(5, -0.5, 0.5);

        if ( NUCLEUS == "33Ar")
        {
            fReleaseCurve[NUCLEUS]->SetParLimits(0, 0, 1);
            fReleaseCurve[NUCLEUS]->SetParameter(0, 6.19e-8);
            fReleaseCurve[NUCLEUS]->FixParameter(1, 0);
            fReleaseCurve[NUCLEUS]->SetParLimits(2, 0, 100);
            fReleaseCurve[NUCLEUS]->SetParameter(2, 1e-2);
            fReleaseCurve[NUCLEUS]->FixParameter(3, 0);
            fReleaseCurve[NUCLEUS]->SetParLimits(4, 0, 1);
            fReleaseCurve[NUCLEUS]->SetParameter(4, MAP_LifeTime[NUCLEUS]);
            fReleaseCurve[NUCLEUS]->SetParLimits(5, -0.5, 0.5);
            fReleaseCurve[NUCLEUS]->SetParameter(5, 5.6e-2);
        }
        TFitResultPtr resf = H_Decay[NUCLEUS][group]->Fit(("fReleaseCurve_" + NUCLEUS).c_str(), "RES");

        // Printing results
        Success("Results for " + NUCLEUS);
        Success("Realese Time  : " + formatValueWithError(fReleaseCurve[NUCLEUS]->GetParameter(2), fReleaseCurve[NUCLEUS]->GetParError(2), "nolatex") + "s");
        Success("Halflife : " + formatValueWithError(fReleaseCurve[NUCLEUS]->GetParameter(4) * 1000, fReleaseCurve[NUCLEUS]->GetParError(4) * 1000, "nolatex") + " ms");
        Success("Chi2 : " + to_string(fReleaseCurve[NUCLEUS]->GetChisquare() / fReleaseCurve[NUCLEUS]->GetNDF()));

        fTime->cd();
        dir[NUCLEUS]->cd();
        TCanvas *cFit = new TCanvas(("cFit_" + NUCLEUS).c_str(), ("cFit_" + NUCLEUS).c_str(), 800, 600);
        H_Decay[NUCLEUS][group]->Draw("HIST");
        H_Decay[NUCLEUS][group]->GetXaxis()->SetTitle("Time [s]");
        H_Decay[NUCLEUS][group]->GetYaxis()->SetTitle("Counts");
        fReleaseCurve[NUCLEUS]->Draw("SAME");
        TPaveText *pt = new TPaveText(0.7, 0.7, 0.9, 0.9, "NDC");
        pt->AddText(("Halflife : " + formatValueWithError(fReleaseCurve[NUCLEUS]->GetParameter(4) * 1000, fReleaseCurve[NUCLEUS]->GetParError(4) * 1000, "nolatex") + " ms").c_str());
        pt->AddText(("Chi2 : " + to_string(fReleaseCurve[NUCLEUS]->GetChisquare() / fReleaseCurve[NUCLEUS]->GetNDF())).c_str());
        pt->Draw("SAME");
        cFit->Write();


        /////// trying gaussan x exp //////
        TF1 *f = new TF1("f", ConvRelease, 0, 4.6, 4);
        f->SetParLimits(0, 0, 10000);
        f->SetParLimits(1, 0, 4.6);
        f->SetParLimits(2, 0, 5);
        f->SetParameter(3, log(2)/MAP_LifeTime[NUCLEUS]);
        H_Decay[NUCLEUS][group]->Fit(f, "MULTITHREAD R", "", 0, 4.6);
        
        TCanvas *cFitConv = new TCanvas(("cFitConv_" + NUCLEUS).c_str(), ("cFitConv_" + NUCLEUS).c_str(), 800, 600);
        H_Decay[NUCLEUS][group]->Draw("HIST");
        H_Decay[NUCLEUS][group]->GetXaxis()->SetTitle("Time [s]");
        H_Decay[NUCLEUS][group]->GetYaxis()->SetTitle("Counts");
        f->Draw("SAME");
        cFitConv->Write();
    }

    fTime->cd();

    ////#################### 18N THIRD PEAK ###################
    //### Decay
    
    // H_Time_All_3a->Write();
    // run_number = 114;
    // int run_i = 0;
    // string NUCLEUS = "18N";
    // int Peak_Group = 0;
    // int Peak_Position = 0;
    // for (int i = 0; i < H_Data[run_number]->GetNbinsX(); i++)
    // {
    //     if (H_Data[run_number]->GetBinContent(i) > Intensity_Threshold)
    //     {
    //         Peak_Group = Group_Determining(i, NUCLEUS);
    //         if (Peak_Group != 4)
    //             { continue; }
    //         Peak_Position = i + H_Data[run_number]->FindBin(MAP_Delta_Time[NUCLEUS][run_i]);
    //         int j = Peak_Position;
    //         while (H_Data[run_number]->GetBinCenter(j) < H_Data[run_number]->GetBinCenter(Peak_Position) + 4 * Pulse_Period)
    //         {
    //             H_Decay_3a->SetBinContent(j - Peak_Position, H_Decay_3a->GetBinContent(j - Peak_Position) + H_Time_All_3a->GetBinContent(j));
    //             j++;
    //         }
    //     }
    // }

    // H_Decay_3a->Rebin(20);
    // ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
    // TF1 *fReleaseCurve3a = new TF1(("fReleaseCurve_3a" + NUCLEUS).c_str(), ReleaseCurve, 0.1, 4 * Pulse_Period, 6);
    // fReleaseCurve3a->SetNpx(1000);
    // fReleaseCurve3a->SetParLimits(0, 0, 1);
    // fReleaseCurve3a->FixParameter(1, 0);
    // fReleaseCurve3a->SetParLimits(2, 0, 100);
    // fReleaseCurve3a->FixParameter(3, 0);
    // fReleaseCurve3a->SetParLimits(4, 0, 1);
    // fReleaseCurve3a->SetParameter(4, MAP_LifeTime[NUCLEUS]);
    // fReleaseCurve3a->SetParLimits(5, -0.5, 0.5);

    // TFitResultPtr resf = H_Decay_3a->Fit(("fReleaseCurve_3a" + NUCLEUS).c_str(), "RES");

    // // Printing results
    // Success("Results for third alpha peak of " + NUCLEUS);
    // Success("Realese Time  : " + formatValueWithError(fReleaseCurve[NUCLEUS]->GetParameter(2), fReleaseCurve[NUCLEUS]->GetParError(2), "nolatex") + "s");
    // Success("Halflife : " + formatValueWithError(fReleaseCurve3a->GetParameter(4) * 1000, fReleaseCurve3a->GetParError(4) * 1000, "nolatex") + " ms");
    // Success("Chi2 : " + to_string(fReleaseCurve3a->GetChisquare() / fReleaseCurve3a->GetNDF()));

    // fTime->cd();
    // dir[NUCLEUS]->cd();
    // TCanvas *cFit = new TCanvas(("cFit_3alpha_" + NUCLEUS).c_str(), ("cFit_3alpha" + NUCLEUS).c_str(), 800, 600);
    // H_Decay_3a->Draw("HIST");
    // H_Decay_3a->GetXaxis()->SetTitle("Time [s]");
    // H_Decay_3a->GetYaxis()->SetTitle("Counts");
    // fReleaseCurve3a->Draw("SAME");
    // TPaveText *pt = new TPaveText(0.7, 0.7, 0.9, 0.9, "NDC");
    // pt->AddText(("Halflife : " + formatValueWithError(fReleaseCurve3a->GetParameter(4) * 1000, fReleaseCurve3a->GetParError(4) * 1000, "nolatex") + " ms").c_str());
    // pt->AddText(("Chi2 : " + to_string(fReleaseCurve3a->GetChisquare() / fReleaseCurve3a->GetNDF())).c_str());
    // pt->Draw("SAME");
    // cFit->Write();


    ////////////
    cout << "ok" << endl;
    TH1D* h = (TH1D*)H_Decay["32Ar"][2]->Clone("h");
    h->Rebin(10);
    TH1D* hh = (TH1D*)H_Decay["32Ar"][2]->Clone("hh");
    hh->Reset();
    hh->Rebin(10);
    TH1D* h_exp;
    TF1* f = new TF1("f", "[0]*exp(-0.69*(x-[2])/[1])", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    f->SetParameter(1, MAP_LifeTime["32Ar"]);
    f->SetNpx(480);
    TH1D* h_res = (TH1D*)h->Clone("h_res");
    h_res->Reset();

    for (int bin = 0; bin < h->GetNbinsX()-1; bin++)
    {
        if ((h->GetBinContent(bin-1) < h->GetBinContent(bin)) && h->GetBinCenter(bin) > 0.4)
            h->SetBinContent(bin, h->GetBinContent(bin-1));
    }

    TCanvas *cFitConv = new TCanvas("cFitConv", "cFitConv", 800, 600);
    cFitConv->cd();
    for (int bin = 0; bin < h->GetNbinsX(); bin++)
    {
        double value = h->GetBinContent(bin);
        double x = h->GetBinCenter(bin);

        if (value <= 0)
            continue;

        f->SetParameter(0, value);
        f->SetParameter(2, x);
        h_exp = (TH1D *)f->GetHistogram()->Clone("h_exp");
        for (int bin_exp = 0; bin_exp < h_exp->GetNbinsX(); bin_exp++)
        {
            if (h_exp->GetBinCenter(bin_exp) < x)
                h_exp->SetBinContent(bin_exp, 0);
        }

        h_exp->Draw("SAME");
        h_res->Add(h_exp);
        h->Add(h_exp, -1);
        hh->SetBinContent(bin, value);

    }
    cFitConv->Write();

    h_res->Write();
    TCanvas *cFitConv2 = new TCanvas("cFitConv2", "cFitConv2", 800, 600);
    h_res->SetLineColor(kRed);
    h_res->Draw("HIST");
    H_Decay["32Ar"][2]->Draw("HIST SAME");
    cFitConv2->Write();
    hh->Write();
    h->Write();
    ///////////

    cout << "Writing histograms" << endl;
    WriteHistograms();
    fTime->Close();
}