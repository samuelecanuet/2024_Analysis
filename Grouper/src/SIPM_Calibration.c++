#include "SIPM_Calibration.hh"


int main()
{
    // INITIALIZATION
    InitDetectors("Config_Files/sample.pid");


    ///////////////////////////////////  FILES //////////////////////////////////
    GROUPED_File["32Ar"] = new TFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
    GROUPED_File["32Ar_thick"] = new TFile((DIR_ROOT_DATA_MERGED + "32Ar_thick_merged.root").c_str(), "READ");
    GROUPED_File["33Ar"] = new TFile((DIR_ROOT_DATA_MERGED + "33Ar_merged.root").c_str(), "READ");

    ///////////////////////////////////  OUTPUT ///////////////////////////////////
    CALIBRATED_File = new TFile((DIR_ROOT_DATA_CALIBRATED + "SiPM_Calibrated.root").c_str(), "RECREATE");
    CALIBRATED_File->cd();


    InitWindows();
    InitSiliconCalibration();
    InitHistograms();


    // Matching SiPM LOW / HIGH
    NUCLEUS = "32Ar";

    Tree = (TTree *)GROUPED_File[NUCLEUS]->Get("MERGED_Tree");
    Reader = new TTreeReader(Tree);
    Silicon = new TTreeReaderArray<Signal>(*Reader, "MERGED_Tree_Silicon");
    SiPM_High = new TTreeReaderArray<Signal>(*Reader, "MERGED_Tree_SiPMHigh");
    SiPM_Low = new TTreeReaderArray<Signal>(*Reader, "MERGED_Tree_SiPMLow");


    FittingLowHigh();
    FittingSiPMs();
    MergingSiPMs();

    
    
    WriteHistograms();
    CALIBRATED_File->Close();
    GROUPED_File["32Ar"]->Close();
    GROUPED_File["32Ar_thick"]->Close();
    GROUPED_File["33Ar"]->Close();

    return 0;
}