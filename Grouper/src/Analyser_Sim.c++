#include "Analyser_Sim.hh"

int main(int argc, char **argv)
{
    
    NUCLEUS = "32Ar";
    FLAG2025 = true;
    VERBOSE = 0;


    InitDetectors("Config_Files/sample.pid");
        
    ///////// READING FILE //////////
    TTreeReader *Reader;
    TTreeReaderArray<Signal> *Silicon;
    TTreeReaderValue<vector<vector<pair<Signal, Signal>>>> *SiPM_Groups;

    // SIM_File = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA32Ar_ENSDF_2025_Default_analysed_new.root").c_str(), "READ");
    // SIM_File = MyTFile("/run/media/local1/DATANEX/Samuel-G4/32Ar_ENSDF_"+to_string(YEAR)+"_Default_analysed_new.root", "READ");
    SIM_File = MyTFile(argv[1], "READ");
    // SIM_File = MyTFile("/run/media/local1/DATANEX/Samuel-G4/32Ar_IAS_2025_x-0.5000_y1.5000_sx0.3364_sy0.3933_catcher5deg_analysed_new.root", "READ");
   
   
    
    ///////// NEW FILE //////////
    string name_to_handle_at_end = string(argv[1]);
    size_t lastindex = name_to_handle_at_end.find_last_of("_analysed_new.root");
    string rawname = name_to_handle_at_end.substr(0, lastindex - 17);
    ANALYSED_File = MyTFile((rawname + "_result.root").c_str(), "READ");
    if (ANALYSED_File != nullptr)
        Error("File " + (rawname + "_result.root") + " already exists!");
    ANALYSED_File = MyTFile((rawname + "_result.root").c_str(), "RECREATE");
    // ANALYSED_File = MyTFile((DIR_ROOT_DATA_SIMULATED + "Result/" + NUCLEUS + "_" + to_string(YEAR) + "_result_1deg.root").c_str(), "RECREATE");

    InitResolution_SiPM();
    // InitResolution_Silicon();
    InitWindows();
    InitHistograms();

    ReaderData();
    
    WriteHistograms();
    ANALYSED_File->Close();
    SIM_File->Close();

    return 0;
}