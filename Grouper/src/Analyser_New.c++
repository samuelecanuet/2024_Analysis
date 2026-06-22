#include "Analyser_New.hh"

int main(int argc, char **argv)
{
    
    NUCLEUS = "32Ar";
    FLAG2025 = true;
    VERBOSE = 0;

    InitDetectors("Config_Files/sample.pid");

    PEAK = IAS[NUCLEUS];
    MULTIPLICITY = 3; // at least
    THRESHOLD = 100; //keV

    if (argc == 3)  
    {
        MULTIPLICITY = atoi(argv[1]);
        THRESHOLD = atoi(argv[2]);
        Info("MULTIPLICITY: " + to_string(MULTIPLICITY), 1);
        Info("THRESHOLD: " + to_string(THRESHOLD) + " keV", 1);
    }
    else if (argc < 3)
    {
        Warning("-- DEFAULT PARAMETERS --");
        Info("MULTIPLICITY: " + to_string(MULTIPLICITY), 1);
        Info("THRESHOLD: " + to_string(THRESHOLD) + " keV", 1);
    }
    else
    {
        Error("Too many arguments. Usage: ./Analyser_new [MULTIPLICITY] [THRESHOLD]");
        return 1;
    }
        
    ///////// READING FILE //////////
    TTreeReader *Reader;
    TTreeReaderArray<Signal> *Silicon;
    TTreeReaderValue<vector<vector<pair<Signal, Signal>>>> *SiPM_Groups;

    MERGED_File = MyTFile((DIR_ROOT_DATA_MERGED + NUCLEUS + "_merged.root").c_str(), "READ");
   
    
    ///////// NEW FILE //////////
    ANALYSED_File = MyTFile(Form("%s%s_%d_analysed_new_N%.1f_M%d_T%.0f.root", DIR_ROOT_DATA_ANALYSED.c_str(), NUCLEUS.c_str(), YEAR, PEAK, MULTIPLICITY, THRESHOLD), "RECREATE");
    TTree *ANALYSED_Tree = new TTree("ANALYSED_Tree", "ANALYSED_Tree");
    InitCalib(0);
    InitWindows();
    InitHistograms();

    ///////
    LoadFitParameters();
    //////

    ReaderData();

    RandomCorrection();

    WriteHistograms();

    ANALYSED_File->Close();

    return 0;
}