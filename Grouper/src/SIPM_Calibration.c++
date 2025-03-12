#include "SIPM_Calibration.hh"


int main()
{
    // INITIALIZATION
    InitDetectors("Config_Files/sample.pid");


    ///////////////////////////////////  FILES //////////////////////////////////
    GROUPED_File["32Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
    GROUPED_File["207Bi"] = MyTFile((DIR_ROOT_DATA_GROUPED + "run_137_multifast_207Bi_grouped.root").c_str(), "READ");
    GROUPED_File["90Sr"] = MyTFile((DIR_ROOT_DATA_GROUPED + "run_133_data_90Sr_grouped.root").c_str(), "READ");
    SIMULATED_File["32Ar"] = MyTFile((DIR_ROOT_DATA_SIMULATED + "s02-24/32Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
    SIMULATED_File["207Bi"] = MyTFile((DIR_ROOT_DATA_SIMULATED + "02-24/207Bi_thin_analysed.root").c_str(), "READ");
    SIMULATED_File["90Sr"] = MyTFile((DIR_ROOT_DATA_SIMULATED + "02-24/90Sr__analysed.root").c_str(), "READ");

    ///////////////////////////////////  OUTPUT ///////////////////////////////////
    CALIBRATED_File = MyTFile((DIR_ROOT_DATA_CALIBRATED + "SiPM_Calibrated.root").c_str(), "RECREATE");
    CALIBRATED_File->cd();
    
    ///////////////////////////////  INITIALISATION ///////////////////////////////
    InitWindows();
    InitSiliconCalibration();
    InitHistograms();

    ////////////////////////  EXPERIMETAL TREE FOR EACH PEAK //////////////////////
    InitTree("READ");    

    ////////////////////////  SIMULATION TREE FOR EACH PEAK //////////////////////
    InitSimulatedTree("READ");

    Resolution_SQRT_det = {0, 5.8, 5.4, 4.3987, 4.3987, 5.8, 4.8, 4.3987, 4.3987, 5.6};
    
    for (int det  = 1; det <= 9; det++)
    {
        Start("Detector: SiPM " + to_string(det));
        current_detector = det;
        CalibrationSiPM();
    }
   
    
    WriteHistograms();
    CALIBRATED_File->Close();
    SIMULATED_File["32Ar"]->Close();
    GROUPED_File["32Ar"]->Close();
    return 0;
}