#include "Systematic_MC.hh"


int main()
{
    FLAG2024 = true;
    NUCLEUS = "32Ar";
    InitDetectors("../Grouper/Config_Files/sample.pid");
    InitWindows();

    string sys = "beam_pos";

    FINAL_File = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/Sytematic_" + sys + "_" + to_string(YEAR) + ".root").c_str(), "RECREATE");


    InitAnalysis();

    vector<string> par = {"y"};

    ComputingParameter(par);
    AnalyseParameter(par);

    WriteHistograms(par);
    FINAL_File->Close();

    return 0;
}