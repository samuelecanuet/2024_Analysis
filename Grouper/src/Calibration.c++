#include "Calibration.hh"


int main(int argc, char *argv[])
{
    InitDetectors("Config_Files/sample.pid");
   
    ///////////////////////////////////  FILES //////////////////////////////////
    GROUPED_File["32Ar"] = new TFile((DIR_ROOT_DATA_GROUPED+"merge3.root").c_str(), "READ");
    GROUPED_File["33Ar"] = new TFile((DIR_ROOT_DATA_GROUPED+"run_078_multifast_33Ar_grouped.root").c_str(), "READ");

    SIMULATED_File["32Ar"] = new TFile((DIR_ROOT_DATA_SIMULATED+"32ArRMATRIX__CS0_CSP0_CV1_CVP1.root").c_str(), "READ");
    SIMULATED_File["33Ar"] = new TFile((DIR_ROOT_DATA_SIMULATED+"33ArRMATRIX__CS0_CSP0_CV1_CVP1.root").c_str(), "READ");
    SIMULATED_File["18N"] = new TFile((DIR_ROOT_DATA_SIMULATED+"18N__CS0_CSP0_CV1_CVP1.root").c_str(), "READ");
    //SIMULATED_File["18N"] = new TFile("../../WISArD/SAMPLE.root", "READ");
    
    ///////////////////////////////////  OUTPUT ///////////////////////////////////
    MATCHED_File = new TFile((DIR_ROOT_DATA_MATCHED+"Calibrated.root").c_str(), "RECREATE");
    MATCHED_File->cd();
    // WriteTime(GROUPED_File[], MATCHED_File);

    ///////////////////////////////////////////////////////////////////////////////
    string Nuclei[2] = {"32Ar", "33Ar"};
    int first = 11;
    int last = 86;

    FillingSimHitograms();
    InitAlphaPeaks();
    InitWindows();
    InitManualCalibration();

    /// LINEAR MANUAL CALIB ON 32Ar
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
                GROUPED_Tree_Detectors[NUCLEUS][i] = (TTree *)GROUPED_File[NUCLEUS]->Get(("CLEANED_Tree_" + detectorName[i]).c_str());
            }
        }

        if (NUCLEUS == "32Ar")
        {
            Info("Calibration Peaks for each detector");
            for (int i = first; i < last; i++)
            {
                if (IsDetectorSiliStrip(i))
                {
                    Reader = new TTreeReader(GROUPED_Tree_Detectors[NUCLEUS][i]);
                    current_detector = i;
                    Manual_Calibration();
                }
            }
        }
    }

    Info("Fitting Calibration Peaks");
    for (int i = first; i < last; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            current_detector = i;
            Fitting_Calibration();
            // GaussiannFit(i, H_Exp);
            //  WindowFit(i);
        }
    }

    NUCLEUS="32Ar";

    Info("Apply Calibration & Resolution minimization");
    for (int i = first; i < last; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            Info(detectorName[i]);
            Reader = new TTreeReader(GROUPED_Tree_Detectors[NUCLEUS][i]);
            current_detector = i;
            CHI2Minimization();
        }
    }


    PlottingWindows();
    
    MATCHED_File->Close();
    // GROUPED_File["32Ar"]->Close();
    // GROUPED_File["33Ar"]->Close();
    return 0;
}
