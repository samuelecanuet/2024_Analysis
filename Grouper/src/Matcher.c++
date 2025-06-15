#include "Matcher.hh"

int main(int argc, char *argv[])
{
    string Run_string;

    if (argc == 1)
    {
        Error("No Run Number Given");
    }
    else
    {
        if (strcmp(argv[1], "all") == 0)
        {
            Info("All Runs");
        }
        else
        {
            Runs.clear();
            int Run = atoi(argv[1]);
            if (Run < 100)
                Runs.push_back("0" + to_string(Run));
            else
                Runs.push_back(to_string(Run));
        }
    }

    FLAG2021 = true;    
    InitDetectors("Config_Files/sample.pid");
    InitRuns();
    InitManualCalibration();
    InitWindows();

    ///////////////////////////////////  REFERENCE  //////////////////////////////
    REFERENCE_filename = SearchFiles(DIR_ROOT_DATA_GROUPED, "0"+to_string(MATCHING_RUN));
    REFERENCE_File = MyTFile((DIR_ROOT_DATA_GROUPED + REFERENCE_filename).c_str(), "READ");
    MATCHED_File = MyTFile((DIR_ROOT_DATA_MATCHED + "matched.root").c_str(), "RECREATE");

    ///////////////////////////////////  INIT HIST  //////////////////////////////
    Info("Init Histograms");
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            C_Det[i] = new TCanvas((detectorName[i]).c_str(), (detectorName[i]).c_str(), 800, 600);
            H_Run_Ref[i] = (TH1D *)REFERENCE_File->Get(("Strip/Strip_Channel/" + detectorName[i] + "/Channel_Cleaned_" + detectorName[i]).c_str());
            H_Run_Ref[i]->SetLineColor(kBlack);
            H_Run_Ref[i]->Draw("HIST");
            C_Det_corr[i] = new TCanvas((detectorName[i] + "_corr").c_str(), (detectorName[i] + "_corr").c_str(), 800, 600);
            H_Run_Ref[i]->Draw("HIST");

            C_Det_vs[i] = new TCanvas((detectorName[i] + "_vs").c_str(), (detectorName[i] + "_vs").c_str(), 800, 600);
            H_Run_BEFORE[i] = (TH1D *)H_Run_Ref[i]->Clone();
            H_Run_BEFORE[i]->SetTitle("Before Correction");
            H_Run_BEFORE[i]->GetXaxis()->SetTitle("Channel");
            H_Run_BEFORE[i]->GetYaxis()->SetTitle("Counts");
            H_Run_BEFORE[i]->GetXaxis()->CenterTitle();
            H_Run_AFTER[i] = (TH1D *)H_Run_Ref[i]->Clone();
            H_Run_AFTER[i]->SetTitle("After Correction");

            for (int peak : Peaks)
            {
                G_Mean_Run[peak][i] = new TGraphErrors();
                G_Mean_Run_corr[peak][i] = new TGraphErrors();
            }

            G_resPar0[i] = new TGraphErrors();
            G_resPar1[i] = new TGraphErrors();
        }
    }

    ///////////////////////////////////////////////////////////////////////////////

    int counter_R[SIGNAL_MAX] = {0};
    for (auto Run_string : Map_RunFiles["32Ar"])
    {
        if (atoi(Run_string.c_str()) == MATCHING_RUN)
            continue;
        Info("Current Run : " + Run_string);
        InitHistograms(Run_string);
        Run = atoi(Run_string.c_str());

        ///////////////////////////////////  INPUT ///////////////////////////////////
        GROUPED_filename = SearchFiles(DIR_ROOT_DATA_GROUPED, Run_string);
        GROUPED_basefilename = GROUPED_filename.substr(0, GROUPED_filename.find(".root"));
        GROUPED_File = MyTFile((DIR_ROOT_DATA_GROUPED + GROUPED_filename).c_str(), "READ");

        ///////////////////////////////////  OUTPUT ///////////////////////////////////
        MATCHED_File->cd();

        ///////////////////////////////////////////////////////////////////////////////
        clock_t start = clock(), Current;

        TTree *Tree = (TTree *)GROUPED_File->Get("CLEANED_Tree");
        if (Tree == NULL)
        {
            Warning("No Tree found in " + GROUPED_filename);
            continue;
        }

        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorSiliStrip(i))
            {
                GROUPED_Tree_Detectors[i] = (TTree *)GROUPED_File->Get(("CLEANED_Tree_" + detectorName[i]).c_str());
            }
        }   

        
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorSiliStrip(i))
            {
                // Info(detectorName[i], 1);
                Reader = new TTreeReader(GROUPED_Tree_Detectors[i]);
                CHI2Minimization(i);

                H_Run_BEFORE[i]->Add(H_Run[Run][i]);
                H_Run_AFTER[i]->Add(H_Run_corr[Run][i]);

                G_resPar0[i]->SetPoint(counter_R[i], Run, fpol1[Run][i]->GetParameter(0));
                G_resPar0[i]->SetPointError(counter_R[i], 0, fpol1[Run][i]->GetParError(0));
                G_resPar1[i]->SetPoint(counter_R[i], Run, fpol1[Run][i]->GetParameter(1));
                G_resPar1[i]->SetPointError(counter_R[i], 0, fpol1[Run][i]->GetParError(1));
                counter_R[i]++;
                fpol1[Run][i]->Write();
            }
        }
        GROUPED_File->Close();
    }

    MATCHED_File->cd();
    TCanvas *c_std = new TCanvas("STD", "STD", 800, 600);
    TCanvas *c_mean = new TCanvas("MEAN", "MEAN", 800, 600);
    TGraphErrors *G_std = new TGraphErrors();
    TGraphErrors *G_mean = new TGraphErrors();
    TCanvas *c_std_e = new TCanvas("STD_e", "STD_e", 800, 600);
    TCanvas *c_mean_e = new TCanvas("MEAN_e", "MEAN_e", 800, 600);
    TGraphErrors *G_std_e = new TGraphErrors();
    TGraphErrors *G_mean_e = new TGraphErrors();
    int counter = 0;
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            C_Det[i]->Write();
            C_Det_corr[i]->Write();

            TCanvas *c = new TCanvas((detectorName[i] + "_Mean").c_str(), (detectorName[i] + "_Mean").c_str(), 800, 600);
            c->Divide(2, 2);

            int counter_c = 1;
            for (double peak : Peaks)
            {
                c->cd(counter_c);
                H_Run_Ref[i]->GetXaxis()->SetRangeUser(WindowsMap["32Ar"][peak][i].first, WindowsMap["32Ar"][peak][i].second);

                G_Mean_Run[peak][i]->SetMarkerStyle(20);
                G_Mean_Run[peak][i]->SetMarkerColor(kGray);
                G_Mean_Run[peak][i]->SetLineColor(kGray);
                G_Mean_Run[peak][i]->SetTitle(("Peak " + to_string(peak)).c_str());
                G_Mean_Run[peak][i]->GetXaxis()->SetTitle("Run");
                G_Mean_Run[peak][i]->GetYaxis()->SetTitle("Mean (Channel)");
                G_Mean_Run[peak][i]->GetXaxis()->CenterTitle();
                G_Mean_Run[peak][i]->GetYaxis()->CenterTitle();
                G_Mean_Run[peak][i]->Draw("AP");

                TF1 *fpol0 = new TF1("fpol0", "[0]", 0, 10000);
                fpol0->FixParameter(0, H_Run_Ref[i]->GetMean());
                fpol0->SetLineColor(kRed);
                G_Mean_Run[peak][i]->Fit(fpol0, "Q");
                fpol0->Draw("SAME");

                G_Mean_Run_corr[peak][i]->SetMarkerStyle(20);
                G_Mean_Run_corr[peak][i]->SetMarkerColor(kBlack);
                G_Mean_Run_corr[peak][i]->SetLineColor(kBlack);
                G_Mean_Run_corr[peak][i]->Draw("P SAME");

                TF1 *fpol0_corr = new TF1("fpol0_corr", "[0]", 0, 10000);
                fpol0_corr->FixParameter(0, H_Run_Ref[i]->GetMean());
                fpol0_corr->SetLineColor(kRed);
                G_Mean_Run_corr[peak][i]->Fit(fpol0_corr, "Q");
                fpol0_corr->Draw("SAME");

                TPaveText *pt = new TPaveText(0.7, 0.7, 0.9, 0.9, "NDC");
                pt->AddText(TString::Format("#chi^{2}_{Before} : %.2f", fpol0->GetChisquare() / fpol0->GetNDF()));
                pt->AddText(TString::Format("#chi^{2}_{After} : %.2f", fpol0_corr->GetChisquare() / fpol0_corr->GetNDF()));
                pt->Draw("SAME");

                TGraph *G_Ref = new TGraph();
                G_Ref->SetPoint(0, MATCHING_RUN, H_Run_Ref[i]->GetMean());
                G_Ref->SetMarkerStyle(20);
                G_Ref->SetMarkerColor(kRed);
                G_Ref->SetLineColor(kRed);
                G_Ref->Draw("P SAME");
                counter_c++;
            }
            c->Write();

            C_Det_vs[i]->cd();
            H_Run_BEFORE[i]->Draw("HIST");
            H_Run_AFTER[i]->SetLineColor(kRed);
            H_Run_AFTER[i]->Draw("HIST SAME");

            TPaveText *pt = new TPaveText(0.7, 0.7, 0.9, 0.9, "NDC");
            H_Run_BEFORE[i]->GetXaxis()->SetRangeUser(WindowsMap["32Ar"][IAS["32Ar"]][i].first, WindowsMap["32Ar"][IAS["32Ar"]][i].second);
            H_Run_AFTER[i]->GetXaxis()->SetRangeUser(WindowsMap["32Ar"][IAS["32Ar"]][i].first, WindowsMap["32Ar"][IAS["32Ar"]][i].second);
            pt->AddText(("STD Before : " + to_string(H_Run_BEFORE[i]->GetRMS())).c_str());
            pt->AddText(("STD After : " + to_string(H_Run_AFTER[i]->GetRMS())).c_str());
            pt->Draw("SAME");
            C_Det_vs[i]->Write();

            // STD&MEAN FOR EACH DETECTOR IN CHANNEL
            G_std->SetPoint(counter, i, H_Run_BEFORE[i]->GetRMS() - H_Run_AFTER[i]->GetRMS());
            G_std->SetPointError(counter, 0, sqrt(pow(H_Run_BEFORE[i]->GetRMSError(), 2) + pow(H_Run_AFTER[i]->GetRMSError(), 2)));
            G_mean->SetPoint(counter, i, abs(H_Run_BEFORE[i]->GetMean() - H_Run_AFTER[i]->GetMean()));
            G_mean->SetPointError(counter, 0, sqrt(pow(H_Run_BEFORE[i]->GetMeanError(), 2) + pow(H_Run_AFTER[i]->GetMeanError(), 2)));

            // STD&MEAN FOR EACH DETECTOR IN ENERGY
            G_std_e->SetPoint(counter, i, (H_Run_BEFORE[i]->GetRMS() - H_Run_AFTER[i]->GetRMS()) * 70 / 1000);
            G_std_e->SetPointError(counter, 0, (sqrt(pow(H_Run_BEFORE[i]->GetRMSError(), 2) + pow(H_Run_AFTER[i]->GetRMSError(), 2))) * 70 / 1000);
            G_mean_e->SetPoint(counter, i, (abs(H_Run_BEFORE[i]->GetMean() - H_Run_AFTER[i]->GetMean())) * 70 / 1000);
            G_mean_e->SetPointError(counter, 0, (sqrt(pow(H_Run_BEFORE[i]->GetMeanError(), 2) + pow(H_Run_AFTER[i]->GetMeanError(), 2))) * 70 / 1000);
            counter++;
        }
    }

    c_std->cd();
    G_std->SetTitle("STD Before - STD After");
    G_std->GetXaxis()->SetTitle("Detector");
    G_std->GetYaxis()->SetTitle("#Delta STD (Channel)");
    G_std->SetMarkerStyle(20);
    G_std->Draw("AP");
    c_std->Write();

    c_mean->cd();
    G_mean->SetTitle("Mean STD Before - STD After");
    G_mean->GetXaxis()->SetTitle("Detector");
    G_mean->GetYaxis()->SetTitle("#Delta Mean (Channel)");
    G_mean->SetMarkerStyle(20);
    G_mean->Draw("AP");
    c_mean->Write();

    c_std_e->cd();
    G_std_e->SetTitle("STD Before - STD After");
    G_std_e->GetXaxis()->SetTitle("Detector");
    G_std_e->GetYaxis()->SetTitle("#Delta STD (keV)");
    G_std_e->SetMarkerStyle(20);
    G_std_e->Draw("AP");
    c_std_e->Write();

    c_mean_e->cd();
    G_mean_e->SetTitle("Mean STD Before - STD After");
    G_mean_e->GetXaxis()->SetTitle("Detector");
    G_mean_e->GetYaxis()->SetTitle("#Delta Mean (keV)");
    G_mean_e->SetMarkerStyle(20);
    G_mean_e->Draw("AP");
    c_mean_e->Write();

    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            // result
            G_resPar0[i]->SetMarkerStyle(20);
            G_resPar0[i]->SetTitle(detectorName[i].c_str());
            G_resPar0[i]->SetName(("Par0_" + detectorName[i]).c_str());
            G_resPar0[i]->SetMarkerColor(kBlack);
            G_resPar0[i]->GetXaxis()->SetTitle("Run");
            G_resPar0[i]->GetYaxis()->SetTitle("Offset");
            G_resPar0[i]->Write();
            G_resPar1[i]->SetMarkerStyle(20);
            G_resPar1[i]->SetTitle(detectorName[i].c_str());
            G_resPar1[i]->SetName(("Par1_" + detectorName[i]).c_str());
            G_resPar1[i]->SetMarkerColor(kBlack);
            G_resPar1[i]->GetXaxis()->SetTitle("Run");
            G_resPar1[i]->GetYaxis()->SetTitle("Proportionnality");
            G_resPar1[i]->Write();
        }
    }

    MATCHED_File->Close();
    REFERENCE_File->Close();

    return 0;
}
