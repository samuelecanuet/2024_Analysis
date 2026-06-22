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
        Nuclei = {"32Ar"};  
        Start("DATA Files");
        MERGED_File["32Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
        // MERGED_File["32Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_merged_around33.root").c_str(), "READ");
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
    SiPM_used = 7;

    // #
    InitBetaSpectrumFit();
    InitResolution_SiPM();

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

            if (peak != 14.0)
                continue;

            PlottingPeak(Nucleus, peak);
        }
    }
    

    // Plotting spectrum with scaling on release curves
    PlottingReleaseScaling("32Ar");
    //

    for (string Nucleus : Nuclei)
    {
        Info("Nucleus : " + Nucleus, 1);
        dir_nuclei[Nucleus]->cd();        

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
    }

    FINAL_FILE->Close();    
    
}