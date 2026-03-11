#include "SIPM_Calibration.hh"

int main()
{
    // INITIALIZATION
    FLAG2025 = true;
    InitDetectors("Config_Files/sample.pid");
    VERBOSE = 1;


    ///////////////////////////////////  FILES //////////////////////////////////
    if (YEAR == 2024)
    {
        Start("DATA Files");
        GROUPED_File["32Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
        // GROUPED_File["33Ar"] = MyTFile((DIR_ROOT_DATA_GROUPED + "run_078_multifast_33Ar_grouped.root").c_str(), "READ");
        GROUPED_File["33Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "33Ar_merged.root").c_str(), "READ");
        // GROUPED_File["207Bi"] = MyTFile((DIR_ROOT_DATA_GROUPED + "run_137_multifast_207Bi_grouped.root").c_str(), "READ");
        // GROUPED_File["90Sr"] = MyTFile((DIR_ROOT_DATA_GROUPED + "run_133_data_90Sr_grouped.root").c_str(), "READ");
        Start("SIMULATED Files");
        // SIMULATED_File["32Ar"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/03-26/32Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["32Ar"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/12-16/32Ar_ENSDF_2024_Default_analysed.root").c_str(), "READ");
        // SIMULATED_File["33Ar"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/03-26/33Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["33Ar"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/12-16/33Ar_ENSDF_2024_Default_analysed.root").c_str(), "READ");
        // SIMULATED_File["207Bi"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/03-17/207Bi_100um_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        // SIMULATED_File["207Bi"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/11-24/207Bi_2024_7p5um_analysed.root").c_str(), "READ");
        // SIMULATED_File["90Sr"] = MyTFile((DIR_ROOT_DATA_SIMULATED + "02-24/90Sr__analysed.root").c_str(), "READ");
        // SIMULATED_File["90Sr"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/11-24/90Sr_2024_7p5um_analysed.root").c_str(), "READ");
    }
    if (YEAR == 2025)
    {
        Start("DATA Files");
        GROUPED_File["32Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
        GROUPED_File["33Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "33Ar_merged.root").c_str(), "READ");
        // GROUPED_File["207Bi"] = MyTFile((DIR_ROOT_DATA_GROUPED + "run_116_data_207Bi_4T_grouped.root").c_str(), "READ");
        // GROUPED_File["90Sr"] = MyTFile((DIR_ROOT_DATA_GROUPED + "run_122_data_90Sr_4T_grouped.root").c_str(), "READ");
        Start("SIMULATED Files");
        // SIMULATED_File["32Ar"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/32Ar_IAS_2025_Default_analysed.root").c_str(), "READ");
        SIMULATED_File["32Ar"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/32Ar_ENSDF_2025_Default_analysed.root").c_str(), "READ");
        SIMULATED_File["33Ar"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/33Ar_ENSDF_2025_Default_analysed.root").c_str(), "READ");
        // SIMULATED_File["207Bi"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/207Bi_2025_7p5um_analysed.root").c_str(), "READ");
        // SIMULATED_File["90Sr"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/90Sr_2025_7p5um_analysed.root").c_str(), "READ");
    }

    ///////////////////////////////////  OUTPUT ///////////////////////////////////
    Start("OUTPUT Files");
    CALIBRATED_File = MyTFile((DIR_ROOT_DATA_CALIBRATED + "SiPM_Calibrated_" + to_string(YEAR) + "_res.root").c_str(), "RECREATE");
    CALIBRATED_File->cd();
    
    ///////////////////////////////  INITIALISATION ///////////////////////////////
    InitWindowss(0, "");
    InitSiliconCalibration();
    InitHistograms();
    InitThresholds();
    InitPileUp();

    ////////////////////////  EXPERIMENTAL TREE FOR EACH PEAK //////////////////////
    InitTree("READ");    

    /////////////////////////  SIMULATION TREE FOR EACH PEAK //////////////////////
    InitSimulatedTree();

    
    for (int det = 5; det <= 5; det++)
    {
        Start("Detector: SiPM " + to_string(det));
        current_detector = det;
        CalibrationSiPM();
    }
   
    
    WriteHistograms();
    CALIBRATED_File->Close();
    SIMULATED_File["32Ar"]->Close();
    SIMULATED_File["33Ar"]->Close();
    // if (SIMULATED_File["207Bi"] != nullptr) SIMULATED_File["207Bi"]->Close();
    // SIMULATED_File["90Sr"]->Close();
    GROUPED_File["32Ar"]->Close();
    GROUPED_File["33Ar"]->Close();
    // GROUPED_File["207Bi"]->Close();
    // GROUPED_File["90Sr"]->Close();
    return 0;
}