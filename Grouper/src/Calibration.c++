#include "Calibration.hh"
#include "CLHEP/Vector/ThreeVector.h"

int main(int argc, char *argv[])
{
    InitDetectors("Config_Files/sample.pid");

    ///////////////////////////////////  FILES //////////////////////////////////
    GROUPED_File["32Ar"] = new TFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
    GROUPED_File["32Ar_thick"] = new TFile((DIR_ROOT_DATA_MERGED + "32Ar_thick_merged.root").c_str(), "READ");
    GROUPED_File["33Ar"] = new TFile((DIR_ROOT_DATA_MERGED + "33Ar_merged.root").c_str(), "READ");

    SIMULATED_File["32Ar"] = new TFile((DIR_ROOT_DATA_SIMULATED + "/17-09/32Ar_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
    SIMULATED_File["32Ar_thick"] = new TFile((DIR_ROOT_DATA_SIMULATED + "/17-09/32Ar_thick_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
    SIMULATED_File["33Ar"] = new TFile((DIR_ROOT_DATA_SIMULATED + "/30-09/33Ar_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
    SIMULATED_File["18N"] = new TFile((DIR_ROOT_DATA_SIMULATED + "18N__CS0_CSP0_CV1_CVP1.root").c_str(), "READ");
    // SIMULATED_File["18N"] = new TFile("../../WISArD/SAMPLE.root", "READ");

    CRADLE_File = new TFile((DIR_ROOT_DATA_SIMULATED + "/new/32ArRMATRIX__CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");

    ///////////////////////////////////  OUTPUT ///////////////////////////////////
    CALIBRATED_File = new TFile((DIR_ROOT_DATA_CALIBRATED + "Calibrated.root").c_str(), "RECREATE");
    CALIBRATED_File->cd();
    // WriteTime(GROUPED_File[], CALIBRATED_File);

    ///////////////////////////////////////////////////////////////////////////////
    int first = 11;
    int last = 16;

    FillingSimHitograms();
    InitAlphaPeaks();
    InitWindows();
    InitManualCalibration();
    InitElectronicResolution();

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
                GROUPED_Tree_Detectors[NUCLEUS][i] = (TTree *)GROUPED_File[NUCLEUS]->Get(("MERGED_Tree_" + detectorName[i]).c_str());
            }
        }
    }

    /// LINEAR MANUAL CALIB ON 32Ar
    NUCLEUS = "32Ar";
    Info("Manuel Calibration");
    for (int i = first; i < last; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            Reader = new TTreeReader(GROUPED_Tree_Detectors[NUCLEUS][i]);
            current_detector = i;
            Manual_Calibration();
            ApplyCalibration();
        }
    }


    // Fitting all peak with quadratic function
    Info("Fitting Calibration Peaks");
    for (int i = first; i < last; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            current_detector = i;
            Fitting_Calibration();
            ApplyCalibration();
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
        }
    }


    PlottingWindows();

    // for (string Nucleus : Nuclei)
    // {
    //     TTreeReader *Reader = new TTreeReader((TTree *)GROUPED_File[Nucleus]->Get("MERGED_Tree"));
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
            if (IsDetectorSiliStrip(i))
            {
                current_detector = i;
                // Fitting_Calibration_E0();
                // ApplyCalibration();
                // CHI2Minimization();
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
    

    CALIBRATED_File->Close();
    GROUPED_File["32Ar"]->Close();
    GROUPED_File["32Ar_thick"]->Close();
    GROUPED_File["33Ar"]->Close();
    return 0;
}
