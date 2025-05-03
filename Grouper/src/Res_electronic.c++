#include "Res_electronic.hh"



int main()
{

    FLAG2025 = true;
    InitDetectors("Config_Files/sample.pid");

    /////////////////////////////  FILES PER YEAR /////////////////////////////
    if (YEAR == 2021)
    {
        file = MyTFile(DIR_ROOT_DATA + "run_047_PulserTest.root", "READ");
        FLAG_PEAK_FINDER = true;
    }
    else if (YEAR == 2024)
    {
        file = MyTFile(DIR_ROOT_DATA + "010_SiPMs_SiDet_Pulser_Po208_-15C_4T.root", "READ");
        FLAG_PEAK_FINDER = true;
    }
    else if (YEAR == 2025)
    {
        file = MyTFile(DIR_ROOT_DATA + "/faster_grouped/run_017_208Po_test_pulser.root", "READ");
        Mean = 38000;
    }
    ////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////  OUTPUT ///////////////////////////////////
    FINAL_File = MyTFile(DIR_ROOT_DATA_CALIBRATED + "Electronic_Resolution_" + to_string(YEAR) + ".root", "RECREATE");
    std::ofstream outfile(("Config_Files/" + to_string(YEAR) + "/Electronic_Resolution_" + to_string(YEAR) + ".txt").c_str());
    ////////////////////////////////////////////////////////////////////////////

    BuildHistograms();

    FittingHistograms();

    WriteHistograms();

    return 0;
}