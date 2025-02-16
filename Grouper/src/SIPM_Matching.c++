#include "SIPM_Matching.hh"


int main(int argc, char **argv)
{
    if (argc == 2 && string(argv[1]) == "FULL")
    {
        FULL = true;
        Info("FULL option USED");
    }
    else
    {
        Info("Genereate only histograms");
    }

    // INITIALIZATION
    InitDetectors("Config_Files/sample.pid");


    ///////////////////////////////////  FILES //////////////////////////////////
    
    ///////////////////////////////////  OUTPUT ///////////////////////////////////
    MATCHED_File = MyTFile(DIR_ROOT_DATA_MATCHED + "SiPM_Matching.root", "RECREATE");
    MATCHED_File->cd();

    init();
    InitWindows();
    InitSiliconCalibration();

    bool first = true;


    Entry_MAX = 7e7;

    // FITTING LOW HIGH ON THE PRINCIPAL LINE 
    // # if FULL option do it and saved in SiPM_Matching_Values.root
    // # else loaded from it
    for (auto &pairr : Map_RunFiles)
    {
        string type = "multifast";
        NUCLEI = pairr.first;
        if (pairr.first == "90Sr") type = "data";
        for (string run : pairr.second)
        {
            Info("Run: " + run);
            RUN = run;
           
            // LOADING FILE
            GROUPED_File[run] = MyTFile(DIR_ROOT_DATA_GROUPED + "run_" + run + "_"+type+"_" + pairr.first + "_grouped.root", "READ");
            if (GROUPED_File[run] == NULL)
                continue;
                
            InitHistograms(RUN, 0);

            // LOADING TREE
            Tree = (TTree *)GROUPED_File[run]->Get("CLEANED_Tree");
            Reader = new TTreeReader(Tree);
            Silicon = NULL;
            if (pairr.first == "32Ar" || pairr.first == "32Ar"|| pairr.first == "33Ar")
            {
                Silicon = new TTreeReaderArray<Signal>(*Reader, "CLEANED_Tree_Silicon");
            }

            SiPM_Groups = new TTreeReaderValue<vector<vector<pair<Signal, Signal>>>>(*Reader, "CLEANED_Tree_SiPMGroup");

            // IF FULL OPTION
            if (FULL)
            {
                // Read Data fill histograms
                ReadData();
                // Fitting Low-High SiPMs
                FittingLowHigh(1);
                FittingSiPM(1);
                // Saving fit in file
                WriteValues(first);
                // Write Histograms betfore correction applyed (where fits comes from)
                WriteHistogram(1);
                
               
            }
            else
            {   
                // Read Data fill histograms
                ReadData();
                // Load fit from file
                int res = LoadValues();
                if (res == 1)
                {
                    FittingLowHigh(1);
                    FittingSiPM(1);
                    // Write Histograms betfore correction applyed (where fits comes from)
                    WriteValues(first);
                }
                // Write Histograms betfore correction applyed (where fits comes from)
                WriteHistogram(1);

                if (first)
                    first = false;
            }      
        }
    }

    WriteHistograms(1);

    Info("APPLYING MATCHING");
    // Applying corrections
    for (auto &pairr : Map_RunFiles)
    {
        string type = "multifast";
        NUCLEI = pairr.first;
        if (pairr.first == "90Sr") type = "data";
        for (string run : pairr.second)
        {
            Info("Run: " + run);
            RUN = run;
           
            // LOADING FILE
            GROUPED_File[run] = MyTFile(DIR_ROOT_DATA_GROUPED + "run_" + run + "_"+type+"_" + pairr.first + "_grouped.root", "READ");
            if (GROUPED_File[run] == NULL)
                continue;
                
            // LOADING TREE
            Tree = (TTree *)GROUPED_File[run]->Get("CLEANED_Tree");
            Reader = new TTreeReader(Tree);
            Silicon = NULL;
            if (pairr.first == "32Ar" || pairr.first == "32Ar"|| pairr.first == "33Ar")
            {
                Silicon = new TTreeReaderArray<Signal>(*Reader, "CLEANED_Tree_Silicon");
            }

            SiPM_Groups = new TTreeReaderValue<vector<vector<pair<Signal, Signal>>>>(*Reader, "CLEANED_Tree_SiPMGroup");

            InitHistograms(RUN, 0, true);
            ReadDataWithCorrections();
            WriteHistogramAfterCorrection(1);      
        }
    }

    MATCHED_File->Close();
    return 0;
}