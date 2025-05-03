#include "Calibration.hh"

int main(int argc, char *argv[])
{
    FLAG2024 = true;
    
    InitDetectors("Config_Files/sample.pid");
    VERBOSE = 0;

    //////////////////////////////  FILES PER YEAR /////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  2021 //////////////////////////////////
    if (YEAR == 2021)
    {
        Nuclei = {"32Ar", "33Ar", "32Ar_thick"};
        Start("DATA Files");
        MERGED_File["32Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
        MERGED_File["32Ar_thick"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_thick_merged.root").c_str(), "READ");
        MERGED_File["33Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "33Ar_merged.root").c_str(), "READ");
        Start("SIMULATED Files");
        SIMULATED_File["32Ar"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/04-01/32Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["32Ar_thick"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/03-26/32Ar_full_thick_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["33Ar"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/03-26/33Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["18N"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/04-01/18N_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");

        CalibrationPeaks["32Ar"] = {5, 8, 9, 14, 23};
        CalibrationPeaks["32Ar_thick"] = {5, 8, 9, 14, 23};
        CalibrationPeaks["33Ar"] = {2, 12, 21, 26, 35, 40};
    }
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  2024 //////////////////////////////////
    else if (YEAR == 2024)
    {
        Nuclei = {"32Ar", "33Ar", "32Ar_thick"};
        Start("DATA Files");
        MERGED_File["32Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
        MERGED_File["32Ar_thick"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_thick_merged.root").c_str(), "READ");
        MERGED_File["33Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "33Ar_merged.root").c_str(), "READ");
        Start("SIMULATED Files");
        SIMULATED_File["32Ar"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/04-01/32Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["32Ar_thick"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/03-26/32Ar_full_thick_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["33Ar"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/03-26/33Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["18N"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/04-01/18N_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        // SIMULATED_File["18N_thick"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/04-01/18N_thick_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");

        CalibrationPeaks["32Ar"] = {5, 8, 9, 14, 23, 25, 28, 29, 30};
        CalibrationPeaks["32Ar_thick"] = {5, 8, 14};
        CalibrationPeaks["33Ar"] = {2, 12, 21, 26, 35, 37, 40};
    }
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  2025 //////////////////////////////////
    else if (YEAR == 2025)
    {
        Nuclei = {"32Ar", "33Ar", "33Ar_thick"};
        Start("DATA Files");
        MERGED_File["32Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
        MERGED_File["33Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "33Ar_merged.root").c_str(), "READ");
        MERGED_File["33Ar_thick"] = MyTFile((DIR_ROOT_DATA_MERGED + "33Ar_thick_merged.root").c_str(), "READ");
        Start("SIMULATED Files");
        SIMULATED_File["32Ar"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/04-01/32Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["33Ar"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/03-26/33Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["33Ar_thick"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/03-26/33Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["18N"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/04-01/18N_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");

        CalibrationPeaks["32Ar"] = {5, 8, 9, 14, 23, 25, 28, 29, 30};
        CalibrationPeaks["33Ar"] = {2, 12, 21, 26, 35};
        CalibrationPeaks["33Ar_thick"] = {2, 12, 21, 26, 35};
    }
    ////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////  OUTPUT ///////////////////////////////////
    Start("OUTPUT Files");
    CALIBRATED_File = MyTFile((DIR_ROOT_DATA_CALIBRATED + "Calibrated_"+ to_string(YEAR) + ".root").c_str(), "RECREATE");
    CALIBRATED_File->cd();

    ///////////////////////////////////////////////////////////////////////////////
    int first = 51;
    int last = 86;

    FillingSimHitograms();
    InitAlphaPeaks();
    InitWindows();
    InitManualCalibration();
    InitElectronicResolution();
    InitPileUp();


    // SET ALL EXPERIMENTAL TREE FOR EACH DETECTOR
    for (string Nucleus : Nuclei)
    {
        NUCLEUS = Nucleus;
        Info("Current Nucleus : " + NUCLEUS);

        InitEnergyErrors();
        InitHistograms();

        // taking exp tree for each detector
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorSiliStrip(i))
            {
                MERGED_Tree_Detectors[NUCLEUS][i] = (TTree *)MERGED_File[NUCLEUS]->Get(("MERGED_Tree_" + detectorName[i]).c_str());
                if (MERGED_Tree_Detectors[NUCLEUS][i] == NULL)
                {
                    MERGED_Tree_Detectors[NUCLEUS][i] = (TTree *)MERGED_File[NUCLEUS]->Get(("CLEANED_Tree_" + detectorName[i]).c_str());
                }
            }
        }
    }

    /// LINEAR MANUAL CALIB ON 32Ar
    NUCLEUS = "32Ar";
    Info("Manual Calibration");
    for (int i = first; i < last; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            Reader = new TTreeReader(MERGED_Tree_Detectors[NUCLEUS][i]);
            current_detector = i;
            Manual_Calibration(VERBOSE);
            ApplyCalibration(VERBOSE);
        }
    }


    // Fitting all peak with quadratic function
    Info("Fitting Calibration Peaks");
    for (int i = first; i < last; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            current_detector = i;
            Fitting_Calibration(VERBOSE);
            ApplyCalibration(VERBOSE);
        }
    }

    // WriteCalibInFile();

    NUCLEUS = "32Ar";

    // Apply Calibration & Resolution minimization
    Info("Apply Calibration & Resolution minimization");
    for (int i = first; i < last; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            Info(detectorName[i]);
            current_detector = i;
            CHI2Minimization();
            Resolution_applied = true;
            // Fitting_Calibration();
            // ApplyCalibration(VERBOSE);
            // CHI2Minimization();
        }
    }


    PlottingWindows();

    // for (string Nucleus : Nuclei)
    // {
    //     TTreeReader *Reader = new TTreeReader((TTree *)MERGED_File[Nucleus]->Get("MERGED_Tree"));
    //     TTreeReaderArray<Signal> *Silicon_Detector = new TTreeReaderArray<Signal>(*Reader, "MERGED_Tree_Silicon");
    //     TTreeReaderArray<Signal> *SIPM_High = new TTreeReaderArray<Signal>(*Reader, "MERGED_Tree_SiPMHigh");

    //     clock_t start = clock(), Current ;
    //     int entries = Reader->GetEntries();
    //     while (Reader->Next())
    //     {
    //         int silicon_code = (*Silicon_Detector)[1].Label;
    //         double silicon_channel = (*Silicon_Detector)[1].Channel;
    //         if (IsDetectorSiliStrip(silicon_code) && silicon_code >= first && silicon_code < last)
    //         {
    //             if (ManualCalibFitted[silicon_code][0] == 0)
    //                 continue;
                
    //             ProgressBar(entries, Reader->GetCurrentEntry(), start, Current);
    //             double silicon_channel_calibrated = Calibration_Function[silicon_code]->Eval(silicon_channel / 1000);

                
    //             if ((*SIPM_High).GetSize() >= 6)
    //             {
    //                 H_Coincidence[Nucleus][silicon_code]->Fill(silicon_channel_calibrated);
    //             }
    //             else
    //             {
    //                 H_AntiCoincidence[Nucleus][silicon_code]->Fill(silicon_channel_calibrated);
    //             }
    //             H_single[Nucleus][silicon_code]->Fill(silicon_channel_calibrated);
    //         }
    //     }

    //     CALIBRATED_File->cd();
    //     for (int i = 0; i < SIGNAL_MAX; i++)
    //     {
    //         if (IsDetectorSiliStrip(i) && i >= first && i < last)
    //         {
    //             dir_detector[i]->cd();
    //             TCanvas *c1 = new TCanvas((Nucleus + "_" + detectorName[i] + "_Coincidence").c_str(), (Nucleus + "_" + detectorName[i] + "_Coincidence").c_str(), 800, 600);
    //             c1->cd();
    //             H_single[Nucleus][i]->SetLineColor(kBlack);
    //             H_single[Nucleus][i]->Draw("HIST");
    //             H_Coincidence[Nucleus][i]->SetLineColor(kRed);
    //             H_Coincidence[Nucleus][i]->Draw("SAME");
    //             H_AntiCoincidence[Nucleus][i]->SetLineColor(kBlue);
    //             H_AntiCoincidence[Nucleus][i]->Draw("SAME");
    //             TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);
    //             legend->AddEntry(H_single[Nucleus][i], "Single proton", "l");
    //             legend->AddEntry(H_Coincidence[Nucleus][i], "#beta-p Coincidence", "l");
    //             legend->AddEntry(H_AntiCoincidence[Nucleus][i], "#beta-p Anti-Coincidence", "l");
    //             legend->Draw("SAME");
    //             c1->Write();
    //         }
    //     }
    // }

    CALIBRATED_File->cd();
    Info("Summing Histograms");

    for (string Nucleus : Nuclei)
    {
        for (int i = first; i < last; i++)
        {
            if (IsDetectorSiliStrip(i) && GetDetectorChannel(i) == 1)
            {
                current_detector = i;
                H_Exp_All[Nucleus]->Add(H_Exp[Nucleus][i], 1.);
                H_Sim_All[Nucleus]->Add(H_Sim_Conv[Nucleus][i], 1.);
            }
        }

        TCanvas *c1 = new TCanvas((Nucleus + "_SUM").c_str(), (Nucleus + "_SUM").c_str(), 800, 600);
        c1->cd();
        H_Exp_All[Nucleus]->Draw("HIST");
        H_Sim_All[Nucleus]->SetLineColor(kRed);
        H_Sim_All[Nucleus]->Draw("HIST SAME");
        c1->Write();

    }

    for (int i = first; i < last; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            Calibration_Function[i]->SetName(("Calibration_" + detectorName[i]).c_str());
            Calibration_Function[i]->Write();
        }
    }
    

    CALIBRATED_File->Close();
    MERGED_File["32Ar"]->Close();
    MERGED_File["32Ar_thick"]->Close();
    MERGED_File["33Ar"]->Close();
    return 0;
}
