#include "Systematic.hh"


int main()
{
    FLAG2025 = true;
    NUCLEUS = "32Ar";
    InitDetectors("../Grouper/Config_Files/sample.pid");
    InitWindows();

    // string sys = "DL";
    // string sys = "Bfield";
    // string sys = "Catcher_thickness";
    string sys = "Silicon_Position_z";

    // limit_a
    a_min = 0.90;
    a_max = 1.1;

    string syspath = "/" + sys + "/";
    if (sys.find("Silicon_Position") != std::string::npos)
        syspath = "/Silicon_Position/";
    
    FINAL_File = MyTFile(("/run/media/local1/DATANEX/Samuel-G4/Systematics/" + syspath + "/Systematics_" + sys + "_" + to_string(YEAR) + ".root").c_str(), "RECREATE");


    // Get the calibration curve a - Eshift
    Start("Initializing a calibration");
    InitAnalysis();
    InitSampling();

    ComputingParameter(sys);
    AnalyseParameter(sys);

    WriteHistograms(sys);
    FINAL_File->Close();

    return 0;
} 