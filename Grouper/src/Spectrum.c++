#include "Spectrum.hh"



int main()
{

    FLAG2025 = true;
    
    InitDetectors("Config_Files/sample.pid");
    VERBOSE = 0;

    //////////////////////////////  FILES PER YEAR /////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  2021 //////////////////////////////////
    if (YEAR == 2021)
    {
        Nuclei = {"32Ar", "33Ar"};
        Start("DATA Files");
        MERGED_File["32Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
        MERGED_File["33Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "33Ar_merged.root").c_str(), "READ");
        Start("SIMULATED Files");
        SIMULATED_File["32Ar"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/04-01/32Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["33Ar"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/03-26/33Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["18N"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/04-01/18N_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");

    }
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  2024 //////////////////////////////////
    else if (YEAR == 2024)
    {
        Nuclei = {"32Ar", "33Ar"};
        Start("DATA Files");
        MERGED_File["32Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
        MERGED_File["33Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "33Ar_merged.root").c_str(), "READ");
        Start("SIMULATED Files");
        SIMULATED_File["32Ar"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/04-01/32Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["33Ar"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/03-26/33Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["18N"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/04-01/18N_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
    }
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  2025 //////////////////////////////////
    else if (YEAR == 2025)
    {
        Nuclei = {"32Ar", "33Ar"};  
        Start("DATA Files");
        // MERGED_File["32Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
        MERGED_File["32Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_merged_around33.root").c_str(), "READ");
        // MERGED_File["32Cl"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
        MERGED_File["33Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "33Ar_merged.root").c_str(), "READ");
        Start("SIMULATED Files");
        SIMULATED_File["32Ar"] = MyTFile(DIR_DATA_HDD + "../SIMULATED_DATA/32Ar_ENSDF_2025_Default_analysed.root", "READ");
        SIMULATED_File["33Ar"] = MyTFile(DIR_DATA_HDD + "../SIMULATED_DATA/33Ar_ENSDF_2025_Default_analysed.root", "READ");
    }

    
    FINAL_FILE = MyTFile(DIR_ROOT_DATA_ANALYSED + "Spectrum.root", "RECREATE");
    InitWindows();
    for (string Nucleus : Nuclei)
        InitPeakData(Nucleus);

    // #All
    InitCalib();

    InitDirectionDeltaEnergy(false);

    //// READING DATA FOR IDENTIFICATION AND SPECTROCOPY ////
    Start("Reading Experimental Data");
    for (string Nucleus : Nuclei)
    {
        InitHistograms(Nucleus);
        Reader = new TTreeReader((TTree*)MERGED_File[Nucleus]->Get("MERGED_Tree"));
        Silicon = new TTreeReaderArray<Signal>(*Reader, "MERGED_Tree_Silicon");
        SiPM_Groups = new TTreeReaderValue<vector<vector<pair<Signal, Signal>>>>(*Reader, "MERGED_Tree_SiPMGroup");
        HRS = new TTreeReaderValue<Signal>(*Reader, "MERGED_Tree_HRS");
        
        ReadingExperimentalData(Nucleus);
    }    

    ComputeDeltaEHist();
    
    for (string Nucleus : Nuclei)
    {
        Info("Nucleus : " + Nucleus, 1);
        dir_nuclei[Nucleus]->cd();

        for (double peak : PeakList[Nucleus])
        {
            if (WindowsMap[Nucleus][peak][11].first == -1 || !WindowsMap[Nucleus][peak][11].first)
                continue;

            PlottingPeak(Nucleus, peak);

            //Doing spectrocopic fit
                // testing fullfunction
        }
    }
    

    // Plotting spectrum with scaling on release curves
    PlottingReleaseScaling("32Ar");
    //

    for (string Nucleus : Nuclei)
    {
        Info("Nucleus : " + Nucleus, 1);
        dir_nuclei[Nucleus]->cd();

        TCanvas *c_EpEb = new TCanvas(Form("BetaProton_Correlation_%s", Nucleus.c_str()), Form("BetaProton_Correlation_%s", Nucleus.c_str()), 1920, 1080);
        c_EpEb->Divide(2, 1);
        c_EpEb->cd(1);
        H_EpEb[Nucleus]["Up"]->Draw("COLZ");
        c_EpEb->cd(2);
        H_EpEb[Nucleus]["Down"]->Draw("COLZ");

        double Qbeta = 11030; // keV, to be updated with real Qbeta
        double Sp = 1580; // keV, to be updated with real Sp
        TF1 *f_GS = new TF1("f_GS", "[0]-[1] - x", 0, 10000);
        f_GS->SetParameters(Qbeta, Sp);
        f_GS->SetLineColor(kRed);
        TF1 *f_Ex1 = new TF1("f_Ex1", "[0]-[1]-2230 - x", 0, 10000);
        f_Ex1->SetParameters(Qbeta, Sp);
        f_Ex1->SetLineColor(kRed);
        c_EpEb->cd(1);
        f_GS->Draw("SAME");
        f_Ex1->Draw("SAME");
        c_EpEb->cd(2);
        f_GS->Draw("SAME");
        f_Ex1->Draw("SAME");
        c_EpEb->Write();
        

        // SIngle /coinc rescale
        TCanvas *c_Spectrum_SingleCoinc = new TCanvas(Form("Spectrum_SingleCoinc_%s", Nucleus.c_str()), Form("Spectrum_SingleCoinc_%s", Nucleus.c_str()), 1920, 1080);
        string dir = "Down";
        double peak = 14.;
        TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        H_Exp_Coinc[Nucleus][dir]->GetXaxis()->SetRangeUser(WindowsMap[Nucleus][peak][51].first, WindowsMap[Nucleus][peak][55].second);
        H_Exp[Nucleus][dir]->GetXaxis()->SetRangeUser(WindowsMap[Nucleus][peak][51].first, WindowsMap[Nucleus][peak][55].second);
        double scale = (double)H_Exp_Coinc[Nucleus][dir]->Integral() / (double)H_Exp[Nucleus][dir]->Integral();
        // H_Exp[Nucleus][dir]->Scale(scale);
        H_Exp_Coinc[Nucleus][dir]->GetXaxis()->SetRangeUser(-1111, -1111);
        H_Exp[Nucleus][dir]->GetXaxis()->SetRangeUser(-1111, -1111);
        H_Exp[Nucleus][dir]->SetTitle("Silicon Spectrum (Down)");
        H_Exp[Nucleus][dir]->SetStats(false);
        H_Exp[Nucleus][dir]->Draw("HIST");
        legend->AddEntry(H_Exp[Nucleus][dir], "Single", "l");
        H_Exp_Coinc[Nucleus][dir]->SetLineColor(kRed);
        H_Exp_Coinc[Nucleus][dir]->SetTitle("Silicon Spectrum (Down)");
        H_Exp_Coinc[Nucleus][dir]->SetStats(false);
        H_Exp_Coinc[Nucleus][dir]->Draw("HIST SAME");
        legend->AddEntry(H_Exp_Coinc[Nucleus][dir], "Coincidence", "l");
        legend->Draw("SAME");
        c_Spectrum_SingleCoinc->Write();

        // Single coinc rescale Litt gs 1st 2nd
        TCanvas *c_Spectrum_SingleCoinc_Litt = new TCanvas(Form("Spectrum_SingleCoinc_Litt_%s_Ex", Nucleus.c_str()), Form("Spectrum_SingleCoinc_Litt_%s_Ex", Nucleus.c_str()), 1920, 1080);
        c_Spectrum_SingleCoinc_Litt->Divide(1, 3);
        string dir_Litt = "Down";
        double peak_Litt = 14.;
        TLegend *legend_Litt = new TLegend(0.7, 0.7, 0.9, 0.9);
        H_Exp_Coinc_Litt[Nucleus][dir_Litt]->GetXaxis()->SetRangeUser(WindowsMap[Nucleus][peak_Litt][51].first, WindowsMap[Nucleus][peak_Litt][55].second);
        H_Exp_Litt[Nucleus][dir_Litt]->GetXaxis()->SetRangeUser(WindowsMap[Nucleus][peak_Litt][51].first, WindowsMap[Nucleus][peak_Litt][55].second);
        double scale_Litt = (double)H_Exp_Coinc_Litt[Nucleus][dir_Litt]->Integral() / (double)H_Exp_Litt[Nucleus][dir_Litt]->Integral();
        // H_Exp_Litt[Nucleus][dir_Litt]->Scale(scale_Litt);
        H_Exp_Coinc_Litt[Nucleus][dir_Litt]->GetXaxis()->SetRangeUser(0, 7000);
        H_Exp_Litt[Nucleus][dir_Litt]->GetXaxis()->SetRangeUser(0, 7000);
        H_Exp_Litt[Nucleus][dir_Litt]->SetTitle("Silicon Spectrum (Down)");
        H_Exp_Litt[Nucleus][dir_Litt]->SetStats(false);
        c_Spectrum_SingleCoinc_Litt->cd(1);
        // GS 
        H_Exp_Litt[Nucleus][dir_Litt]->Draw("HIST");
        legend_Litt->AddEntry(H_Exp_Litt[Nucleus][dir_Litt], "Single", "l");
        // 1st Excited
        c_Spectrum_SingleCoinc_Litt->cd(2);
        TH1D *H_Exp_Litt_1st = (TH1D*)H_Exp_Litt[Nucleus][dir_Litt]->Clone(Form("%s_1st", H_Exp_Litt[Nucleus][dir_Litt]->GetName()));
        // ScaleXaxis(H_Exp_Litt_1st, Scale, 1., -1248.87);
        H_Exp_Litt_1st->GetXaxis()->SetRangeUser(-1248.87, 7000-1248.87);
        H_Exp_Litt_1st->GetXaxis()->SetTitle("E_p + E_{#gamma_{1}} (keV)");
        H_Exp_Litt_1st->Draw("HIST");
        gPad->Modified();
        gPad->Update();
        // 2nd Excited
        c_Spectrum_SingleCoinc_Litt->cd(3);
        TH1D *H_Exp_Litt_2nd = (TH1D*)H_Exp_Litt[Nucleus][dir_Litt]->Clone(Form("%s_2nd", H_Exp_Litt[Nucleus][dir_Litt]->GetName()));
        // ScaleXaxis(H_Exp_Litt_2nd, Scale, 1., -2230.16);
        H_Exp_Litt_2nd->GetXaxis()->SetRangeUser(-2230.16, 7000-2230.16);
        H_Exp_Litt_2nd->GetXaxis()->SetTitle("E_p + E_{#gamma_{2}} (keV)");
        H_Exp_Litt_2nd->Draw("HIST");
        legend_Litt->AddEntry(H_Exp_Coinc_Litt[Nucleus][dir_Litt], "Coincidence", "l");
        legend_Litt->Draw("SAME");
        gPad->Modified();
        gPad->Update();
        
        c_Spectrum_SingleCoinc_Litt->cd(1);
        // new upper axis for pad 1
        TGaxis *axis = new TGaxis(0, 0, 7000+1581.1, 0, 0, 7000+1581.1, 510, "+L");
        axis->SetTitle("E_x (keV)");  
        axis->SetTitleOffset(1.2);
        axis->Draw("SAME");


        c_Spectrum_SingleCoinc_Litt->Write();

    }

    FINAL_FILE->Close();    
    
}