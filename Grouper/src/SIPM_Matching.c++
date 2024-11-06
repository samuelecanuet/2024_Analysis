#include "SIPM_Matching.hh"


int main()
{
    // INITIALIZATION
    InitDetectors("Config_Files/sample.pid");


    ///////////////////////////////////  FILES //////////////////////////////////
    GROUPED_File["32Ar"] = new TFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
    GROUPED_File["207Bi"] = new TFile((DIR_ROOT_DATA + "run_137_multifast_207Bi.root").c_str(), "READ");
    GROUPED_File["90Sr"] = new TFile((DIR_ROOT_DATA + "run_133_data_90Sr.root").c_str(), "READ");

    ///////////////////////////////////  OUTPUT ///////////////////////////////////
    CALIBRATED_File = new TFile((DIR_ROOT_DATA_CALIBRATED + "SiPM_Matching.root").c_str(), "RECREATE");
    CALIBRATED_File->cd();

    NUCLEUS = "32Ar";
    f_tree = new TFile((DIR_ROOT_DATA_CALIBRATED + "Matching_SiPM_trees.root").c_str(), "RECREATE");

    InitWindows();
    InitSiliconCalibration();
    InitHistograms();

    NUCLEUS = "32Ar";
    // Reader for grouped
    Tree = (TTree *)GROUPED_File[NUCLEUS]->Get("MERGED_Tree");
    Reader = new TTreeReader(Tree);
    Silicon = new TTreeReaderArray<Signal>(*Reader, "MERGED_Tree_Silicon");
    SiPM_High = new TTreeReaderArray<Signal>(*Reader, "MERGED_Tree_SiPMHigh");
    SiPM_Low = new TTreeReaderArray<Signal>(*Reader, "MERGED_Tree_SiPMLow");
    
    Entry_MAX = 100000000;

    NUCLEUS = "32Ar";

    FittingLowHigh();
    FittingSiPMs();
    WriteValues();

    // Read values
    // ReadValues();
    MergingSiPMs();

    // Bi207
    NUCLEUS = "207Bi";
    FittingLowHighBi207();
    FittingSiPMsBi207();
    MergingSiPMsBi207();

    NUCLEUS = "90Sr";
    FittingLowHighBi207();
    FittingSiPMsBi207();
    MergingSiPMsBi207();

    WriteTree();
  
    
    WriteHistograms();
    CALIBRATED_File->Close();
    GROUPED_File["32Ar"]->Close();
    return 0;
}