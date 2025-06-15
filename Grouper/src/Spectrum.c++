#include "Spectrum.hh"



int main()
{

    FLAG2024 = true;
    
    InitDetectors("Config_Files/sample.pid");
    VERBOSE = 0;

    //////////////////////////////  FILES PER YEAR /////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  2021 //////////////////////////////////
    if (YEAR == 2021)
    {
        Nuclei = {"32Ar", "33Ar"};
        Start("DATA Files");
        MERGED_File["32Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
        MERGED_File["33Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "33Ar_merged.root").c_str(), "READ");
        Start("SIMULATED Files");
        SIMULATED_File["32Ar"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/04-01/32Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["33Ar"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/03-26/33Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["18N"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/04-01/18N_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");

    }
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  2024 //////////////////////////////////
    else if (YEAR == 2024)
    {
        Nuclei = {"32Ar", "33Ar"};
        Start("DATA Files");
        MERGED_File["32Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
        MERGED_File["33Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "33Ar_merged.root").c_str(), "READ");
        Start("SIMULATED Files");
        SIMULATED_File["32Ar"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/04-01/32Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["33Ar"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/03-26/33Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["18N"] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/04-01/18N_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
    }
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  2025 //////////////////////////////////
    else if (YEAR == 2025)
    {
        Nuclei = {"32Ar", "33Ar", "32Cl"};  
        Start("DATA Files");
        MERGED_File["32Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
        MERGED_File["32Cl"] = MyTFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
        MERGED_File["33Ar"] = MyTFile((DIR_ROOT_DATA_MERGED + "33Ar_merged.root").c_str(), "READ");
        Start("SIMULATED Files");
        SIMULATED_File["32Ar"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/06-01/32Ar_ENSDF_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
        SIMULATED_File["33Ar"] = MyTFile((DIR_DATA_HDD + "../../2024_DATA/SIMULATED_DATA/03-26/33Ar_full_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
    }

    
    FINAL_FILE = MyTFile(DIR_ROOT_DATA_ANALYSED + "Spectrum.root", "RECREATE");
    InitWindows();
    for (string Nucleus : Nuclei)
        InitPeakData(Nucleus);
    
    // #SILICON
    InitCalibration();
    // InitResolution();
    // #SiPM
    // InitSiPM_Calibration();
    // InitSiPM_Resolution();
    InitExperimentalSpectrum();
    InitDirectionDeltaEnergy(false);

    //// READING DATA FOR IDENTIFICATION AND SPECTROCOPY ////
    Start("Reading Experimental Data");
    for (string Nucleus : Nuclei)
    {
        InitHistograms(Nucleus);
        Reader = new TTreeReader((TTree*)MERGED_File[Nucleus]->Get("MERGED_Tree"));
        Silicon = new TTreeReaderArray<Signal>(*Reader, "MERGED_Tree_Silicon");
        SiPM_Groups = new TTreeReaderValue<vector<vector<pair<Signal, Signal>>>>(*Reader, "MERGED_Tree_SiPMGroup");
        HRS = new TTreeReaderValue<Signal>(*Reader, "MERGED_Tree_HRS");
        
        ReadingExperimentalData(Nucleus);
    }    

    ComputeDeltaEHist();
    
    for (string Nucleus : Nuclei)
    {
        Info("Nucleus : " + Nucleus, 1);
        dir_nuclei[Nucleus]->cd();

        for (double peak : PeakList[Nucleus])
        {
            if (WindowsMap[Nucleus][peak][11].first == -1 || !WindowsMap[Nucleus][peak][11].first)
                continue;

            PlottingPeak(Nucleus, peak);

            //Doing spectrocopic fit
                // testing fullfunction
        }
    }

    FINAL_FILE->Close();    
}