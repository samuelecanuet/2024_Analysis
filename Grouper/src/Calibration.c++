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
        Nuclei = {"32Ar", "33Ar"};
        Start("DATA Files");
        MERGED_File["32Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
        // MERGED_File["32Ar_thick"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_thick_merged.root").c_str(), "READ");
        MERGED_File["33Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "33Ar_merged.root").c_str(), "READ");
        Start("SIMULATED Files");
        SIMULATED_File["32Ar"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/04-01/32Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        // SIMULATED_File["32Ar_thick"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/03-26/32Ar_full_thick_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["33Ar"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/03-26/33Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["18N"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/04-01/18N_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");

        CalibrationPeaks["32Ar"] = {5, 8, 9, 14, 23};
        // CalibrationPeaks["32Ar_thick"] = {5, 8, 9, 14, 23};
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
        SIMULATED_File["32Ar"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/06-01/32Ar_ENSDF_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["32Ar_thick"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/03-26/32Ar_full_thick_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["33Ar"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/06-09/33Ar_ENSDF_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["18N"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/04-01/18N_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        // SIMULATED_File["18N_thick"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/04-01/18N_thick_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");

        CalibrationPeaks["32Ar"] = {5, 8, 9, 14, 23, 25, 28, 29, 30};
        // CalibrationPeaks["32Ar_thick"] = {5, 8, 14};
        CalibrationPeaks["33Ar"] = {2, 12, 21, 26, 35, 37, 40};
    }
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  2025 //////////////////////////////////
    else if (YEAR == 2025)
    {
        Nuclei = {"32Ar", "33Ar"};
        Start("DATA Files");
        MERGED_File["32Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
        MERGED_File["33Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "33Ar_merged.root").c_str(), "READ");
        MERGED_File["33Ar_thick"] = MyTFile((DIR_ROOT_DATA_MERGED + "33Ar_thick_merged.root").c_str(), "READ");
        Start("SIMULATED Files");
        SIMULATED_File["32Ar"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/06-01/32Ar_ENSDF_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["33Ar"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/06-09/33Ar_ENSDF_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["33Ar_thick"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/03-26/33Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["32Cl"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/06-01/32Cl_ENSDF_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");   
        SIMULATED_File["18N"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/04-01/18N_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");

        CalibrationPeaks["32Ar"] = {5, 8, 9, 14, 23, 25, 28, 29, 30};
        CalibrationPeaks["33Ar"] = {2, 12, 21, 26, 35};
        // CalibrationPeaks["33Ar_thick"] = {2, 12, 21, 26, 35};
    }
    ////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////  OUTPUT ///////////////////////////////////
    Start("OUTPUT Files");
    CALIBRATED_File = MyTFile((DIR_ROOT_DATA_CALIBRATED + "Calibrated_"+ to_string(YEAR) + "test.root").c_str(), "RECREATE");
    CALIBRATED_File->cd();

    ///////////////////////////////////////////////////////////////////////////////
    int first = 11;
    int last = 86;
    Fitting_Resolution_FLAG = false;
    // Fitting_Calibration_FLAG = false;

    FillingSimHitograms();
    InitAlphaPeaks();
    InitWindows();
    InitManualCalibration();
    InitElectronicResolution();
    InitResolution();
    InitPileUp();
    
    
    Info("Loading Experimental Trees");
    // SET ALL EXPERIMENTAL TREE FOR EACH DETECTOR
    for (string Nucleus : Nuclei)
    {
        NUCLEUS = Nucleus;
        Info(NUCLEUS, 1);

        InitPeakData(NUCLEUS);
        InitHistograms(first, last);

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
    Start("Manual Calibration");
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
    Start("Fitting Calibration & Resolution");   
    for (int i = first; i < last; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            Info(detectorName[i], 1);
            current_detector = i;
            Fitting_Calibration(VERBOSE);
            ApplyCalibration(VERBOSE);
            FittingResolution(VERBOSE);
        }
    }

    // WriteCalibInFile();

    NUCLEUS = "32Ar";
    // Apply Calibration & Resolution minimization
    Start("Applying Calibration & Resolution");
    for (int i = first; i < last; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            Info(detectorName[i]);
            current_detector = i;
            CHI2Minimization();
            Resolution_applied = true;
        }
    }

    Start("Writting Plots");
    PlottingWindows(first, last);
    PlottingSummed(first, last);
    PlottingFitParameters(first, last);

    for (string Nucleus : Nuclei)
    {
        MERGED_File[Nucleus]->Close();
    }
    CALIBRATED_File->Close();

    
    return 0;
}
