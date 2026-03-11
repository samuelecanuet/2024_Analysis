#include "SIPM_Matching.hh"


int main(int argc, char **argv)
{

    FLAG2025 = true;
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
    MATCHED_Filename = "SiPM_Matching_test.root";
    MATCHED_File = MyTFile(DIR_ROOT_DATA_MATCHED + MATCHED_Filename, "RECREATE");
    MATCHED_File->cd();

    init();
    InitWindows();
    InitSiliconCalibration();

    InitPileUp();

    bool first = true;

    Entry_MAX = 1e7;

    // NEW ROUTINE TO MATCH DETECTORS BETWEEN OFF/ON
    Info("FITTING MATCHING PARAMETERS");
    if (FULL && YEAR == 2025)
    {
        InitHistogramsONOFF();
        for (auto &pairr : Map_RunFiles)
        {
            NUCLEI = pairr.first;
            for (string run : pairr.second)
            {
                Info("Run: " + run);
                RUN = run;

                // init hist per run
                MATCHED_File->cd();
                for (int i = 0; i < SIGNAL_MAX; i++)
                {
                    if (IsDetectorBetaHigh(i))
                    {
                        H_SiPM_ONOFF_RUN[i][stoi(run)]= new TH1D(("H_SiPM_ONOFF_RUN_" + detectorName[i] + "_run_" + run).c_str(), ("H_SiPM_ONOFF_RUN_" + detectorName[i] + "_run_" + run).c_str(), eHighN, eHighMin, eHighMax);
                        H_SiPM_ONOFF_RUN[i][stoi(run)]->GetXaxis()->SetTitle("Channel");
                        H_SiPM_ONOFF_RUN[i][stoi(run)]->GetYaxis()->SetTitle("Counts");
                        H_SiPM_ONOFF_RUN[i][stoi(run)]->GetXaxis()->CenterTitle();
                        H_SiPM_ONOFF_RUN[i][stoi(run)]->GetYaxis()->CenterTitle();

                        Tree_SiPM_ONOFF_RUN[i][stoi(run)] = new TTree(("Tree_SiPM_ONOFF_before_" + detectorName[i] + "_run_" + run).c_str(), ("Tree_SiPM_ONOFF_before_" + detectorName[i] + "_run_" + run).c_str());
                        Tree_SiPM_ONOFF_RUN[i][stoi(run)]->Branch("Channel", &Channel);

                    }
                    if (IsDetectorBetaLow(i))
                    {
                        H_SiPM_ONOFF_RUN[i][stoi(run)]= new TH1D(("H_SiPM_ONOFF_RUN_" + detectorName[i] + "_run_" + run).c_str(), ("H_SiPM_ONOFF_RUN_" + detectorName[i] + "_run_" + run).c_str(), eLowN, eLowMin, eLowMax);
                        H_SiPM_ONOFF_RUN[i][stoi(run)]->GetXaxis()->SetTitle("Channel");
                        H_SiPM_ONOFF_RUN[i][stoi(run)]->GetYaxis()->SetTitle("Counts");
                        H_SiPM_ONOFF_RUN[i][stoi(run)]->GetXaxis()->CenterTitle();
                        H_SiPM_ONOFF_RUN[i][stoi(run)]->GetYaxis()->CenterTitle();

                        Tree_SiPM_ONOFF_RUN[i][stoi(run)] = new TTree(("Tree_SiPM_ONOFF_before_" + detectorName[i] + "_run_" + run).c_str(), ("Tree_SiPM_ONOFF_before_" + detectorName[i] + "_run_" + run).c_str());
                        Tree_SiPM_ONOFF_RUN[i][stoi(run)]->Branch("Channel", &Channel);

                    }
                }
            
                // LOADING FILE
                string filename = SearchFiles(DIR_ROOT_DATA_GROUPED, run);
                GROUPED_File[run] = MyTFile(DIR_ROOT_DATA_GROUPED + filename, "READ");
                if (GROUPED_File[run] == NULL)
                    continue;
                    
                    
                // LOADING TREE
                Tree = (TTree *)GROUPED_File[run]->Get("CLEANED_Tree");
                Reader = new TTreeReader(Tree);
                Silicon = nullptr;
                if (pairr.first == "32Ar" || pairr.first == "33Ar")
                {
                    Silicon = new TTreeReaderArray<Signal>(*Reader, "CLEANED_Tree_Silicon");
                    HRS = new TTreeReaderValue<Signal>(*Reader, "CLEANED_Tree_HRS");
                }

                SiPM_Groups = new TTreeReaderValue<vector<vector<pair<Signal, Signal>>>>(*Reader, "CLEANED_Tree_SiPMGroup");
                ReadDataONOFF(stoi(run));
            }
        }
        FittingONOFF();
        // FittingONOFF_Runs();
        WriteONOFFValues(first);
        first = false;
    }
    else if (YEAR != 2025)
    {
        WriteONOFFValues(first, true);
        first = false;
    }

    LoadONOFFValues();
    


    // FITTING LOW HIGH ON THE PRINCIPAL LINE 
    // # if FULL option do it and saved in SiPM_Matching_Values.root
    // # else loaded from it
    for (auto &pairr : Map_RunFiles)
    {
        NUCLEI = pairr.first;
        for (string run : pairr.second)
        {
            Info("Run: " + run);
            RUN = run;
           
            // LOADING FILE
            string filename = SearchFiles(DIR_ROOT_DATA_GROUPED, run);
            GROUPED_File[run] = MyTFile(DIR_ROOT_DATA_GROUPED + filename, "READ");
            if (GROUPED_File[run] == NULL)
                continue;
                
            InitHistograms(RUN, 1);

            // LOADING TREE
            Tree = (TTree *)GROUPED_File[run]->Get("CLEANED_Tree");
            Reader = new TTreeReader(Tree);
            Silicon = nullptr;
            if (pairr.first == "32Ar" || pairr.first == "33Ar")
            {
                Silicon = new TTreeReaderArray<Signal>(*Reader, "CLEANED_Tree_Silicon");
                HRS = new TTreeReaderValue<Signal>(*Reader, "CLEANED_Tree_HRS");
            }

            SiPM_Groups = new TTreeReaderValue<vector<vector<pair<Signal, Signal>>>>(*Reader, "CLEANED_Tree_SiPMGroup");
            

            // IF FULL OPTION
            if (FULL)
            {
                // Read Data fill histograms
                ReadData();
                // ReadData();
                FittingLowHigh(1);
                FittingSiPM(1);
                // Saving fit in file
                WriteValues(first);
                // Write Histograms betfore correction applyed (where fits comes from)
                WriteHistogram(1);
                if (first) first = false;
            }
            else
            {   
                // Read Data fill histograms
                ReadData();
                // Load fit from file
                int res = LoadValues(false);
                if (res == 1)
                {
                    FittingLowHigh(1);
                    FittingSiPM(1);
                    // Write Histograms betfore correction applyed (where fits comes from)
                    WriteValues(first);
                }
                // Write Histograms betore correction applyed (where fits comes from)
                WriteHistogram(1);

                if (first) first = false;
            }      
        }
    }

    WriteHistograms(1);

    Info("APPLYING MATCHING");
    // Applying corrections
    for (auto &pairr : Map_RunFiles)
    {
        NUCLEI = pairr.first;
        for (string run : pairr.second)
        {
            Info("Run: " + run);
            RUN = run;
           
            // LOADING FILE
            string filename = SearchFiles(DIR_ROOT_DATA_GROUPED, run);
            GROUPED_File[run] = MyTFile(DIR_ROOT_DATA_GROUPED + filename, "READ");
            if (GROUPED_File[run] == NULL)
                continue;
                
            // LOADING TREE
            Tree = (TTree *)GROUPED_File[run]->Get("CLEANED_Tree");
            Reader = new TTreeReader(Tree);
            Silicon = NULL;
            if (pairr.first == "32Ar" || pairr.first == "33Ar")
            {
                Silicon = new TTreeReaderArray<Signal>(*Reader, "CLEANED_Tree_Silicon");
                HRS = new TTreeReaderValue<Signal>(*Reader, "CLEANED_Tree_HRS");
            }

            SiPM_Groups = new TTreeReaderValue<vector<vector<pair<Signal, Signal>>>>(*Reader, "CLEANED_Tree_SiPMGroup");

            InitHistograms(RUN, 0, true);
            LoadValues();
            ReadDataWithCorrections();
            WriteHistogramAfterCorrection(1);      
        }
    }

    MATCHED_File->Close();

    WriteHistogramsAfterCorrection(1);

    
    return 0;
}