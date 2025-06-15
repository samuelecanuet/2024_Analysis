#include "SIPM_Calibration.hh"


int main()
{
    
    // INITIALIZATION

    InitDetectors("Config_Files/sample.pid");


    ///////////////////////////////////  FILES //////////////////////////////////
    Start("DATA Files");
    GROUPED_File["32Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
    GROUPED_File["33AR"] = MyTFile((DIR_ROOT_DATA_GROUPED + "run_078_multifast_33Ar_grouped.root").c_str(), "READ");
    GROUPED_File["207Bi"] = MyTFile((DIR_ROOT_DATA_GROUPED + "run_137_multifast_207Bi_grouped.root").c_str(), "READ");
    GROUPED_File["90Sr"] = MyTFile((DIR_ROOT_DATA_GROUPED + "run_133_data_90Sr_grouped.root").c_str(), "READ");
    Start("SIMULATED Files");
    SIMULATED_File["32Ar"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/03-26/32Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
    SIMULATED_File["33Ar"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/03-26/33Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
    SIMULATED_File["207Bi"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/03-17/207Bi_100um_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
    SIMULATED_File["90Sr"] = MyTFile((DIR_ROOT_DATA_SIMULATED + "02-24/90Sr__analysed.root").c_str(), "READ");

    ///////////////////////////////////  OUTPUT ///////////////////////////////////
    Start("OUTPUT Files");
    CALIBRATED_File = MyTFile((DIR_ROOT_DATA_CALIBRATED + "SiPM_Calibrated.root").c_str(), "RECREATE");
    CALIBRATED_File->cd();
    
    ///////////////////////////////  INITIALISATION ///////////////////////////////
    InitWindowss();
    InitSiliconCalibration();
    InitHistograms();
    InitThresholds();
    InitPileUp();

    ////////////////////////  EXPERIMETAL TREE FOR EACH PEAK //////////////////////
    InitTree("RECREATE");    

    /////////////////////////  SIMULATION TREE FOR EACH PEAK //////////////////////
    InitSimulatedTree();

    Resolution_SQRT_det = {0, 5.8, 5.4, 4.3987, 4.3987, 5.8, 4.8, 4.3987, 4.3987, 5.6};
    
    for (int det = 7; det <= 7; det++)
    {
        Start("Detector: SiPM " + to_string(det));
        current_detector = det;
        CalibrationSiPM();
    }
   
    
    WriteHistograms();
    CALIBRATED_File->Close();
    SIMULATED_File["32Ar"]->Close();
    SIMULATED_File["33Ar"]->Close();
    SIMULATED_File["207Bi"]->Close();
    SIMULATED_File["90Sr"]->Close();
    GROUPED_File["32Ar"]->Close();
    GROUPED_File["33Ar"]->Close();
    GROUPED_File["207Bi"]->Close();
    GROUPED_File["90Sr"]->Close();
    return 0;
}