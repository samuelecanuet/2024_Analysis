#include "SIPM_Matching.hh"


int main()
{
    // INITIALIZATION
    InitDetectors("Config_Files/sample.pid");


    ///////////////////////////////////  FILES //////////////////////////////////
    GROUPED_File["32Ar"] = MyTFile(DIR_ROOT_DATA_GROUPED + "run_074_multifast_32Ar_grouped.root", "READ");
    GROUPED_File["207Bi"] = MyTFile(DIR_ROOT_DATA_GROUPED + "run_137_multifast_207Bi_grouped.root", "READ");
    GROUPED_File["90Sr"] = MyTFile(DIR_ROOT_DATA + "run_133_data_90Sr.root", "READ");

    ///////////////////////////////////  OUTPUT ///////////////////////////////////
    CALIBRATED_File = MyTFile(DIR_ROOT_DATA_CALIBRATED + "SiPM_Matching.root", "RECREATE");
    CALIBRATED_File->cd();

    NUCLEUS = "32Ar";
    f_tree = MyTFile(DIR_ROOT_DATA_CALIBRATED + "Matching_SiPM_trees.root", "RECREATE");

    InitWindows();
    InitSiliconCalibration();
    InitHistograms(1);

    Entry_MAX = 7e7;

    NUCLEUS = "32Ar";

    if (NUCLEUS == "32Ar")
    {
        TYPE = "CLEANED";
    }
    else
    {
        TYPE = "CLEANED";
    }


    // LOADING TREE
    Tree = (TTree *)GROUPED_File[NUCLEUS]->Get((TYPE+"_Tree").c_str());
    Reader = new TTreeReader(Tree);
    if (NUCLEUS == "32Ar")
    {
        Silicon = new TTreeReaderArray<Signal>(*Reader, (TYPE+"_Tree_Silicon").c_str());
    }

    SiPM_Groups = new TTreeReaderValue<vector<vector<pair<Signal, Signal>>>>(*Reader, (TYPE+"_Tree_SiPMGroup").c_str());
    

    // Fitting Low-High SiPMs

    FittingLowHigh();
    // FittingSiPMs();
    // WriteValues();

    // // Read values
    // // ReadValues();
    // MergingSiPMs();

    // Bi207
    // NUCLEUS = "207Bi";
    // FittingLowHighBi207();
    // FittingSiPMsBi207();
    // MergingSiPMsBi207();

    // NUCLEUS = "90Sr";
    // FittingLowHighBi207();
    // FittingSiPMsBi207();
    // MergingSiPMsBi207();

    // WriteTree();
  
    
    WriteHistograms(0);
    CALIBRATED_File->Close();
    GROUPED_File["32Ar"]->Close();
    return 0;
}