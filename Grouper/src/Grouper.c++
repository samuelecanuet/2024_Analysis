#include "Grouper.hh"
#include <ctime>

int main(int argc, char *argv[])
{
    
    string Run_string;

    if (argc == 1)
    {
        Error("No Run Number Given");
    }
    else
    {
        Run = atoi(argv[1]);
        if (Run < 10)
            Run_string = "00" + to_string(Run);
        else if (Run < 100)
            Run_string = "0" + to_string(Run);
        else
            Run_string = to_string(Run);
    }

    Info("Current Run : " + Run_string);
    ///////////////////////////////////  INPUT ///////////////////////////////////
    ROOT_filename = SearchFiles(DIR_ROOT_DATA, Run_string);
    ROOT_basefilename = ROOT_filename.substr(0, ROOT_filename.find(".root"));

    TFile *ROOT_File = new TFile((DIR_ROOT_DATA+ROOT_filename).c_str(), "READ");
    TTree *Tree = (TTree *)ROOT_File->Get("Tree_Group");
    TTreeReader *Reader = new TTreeReader(Tree);
    TTreeReaderArray<Signal> signals(*Reader, "Signal");

    ///////////////////////////////////  Grouped ///////////////////////////////////
    GROUPED_File = new TFile((DIR_ROOT_DATA_GROUPED+ROOT_basefilename+"_grouped.root").c_str(), "RECREATE");
    GROUPED_File->cd();
    WriteTime(ROOT_File, GROUPED_File);
    ///////////////////////////////////  INITIALISATION ///////////////////////////////////
    InitDetectors("Config_Files/sample.pid");
    InitHistograms_Grouped();
    InitCalibration();

    /////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  SELECTING TREE /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    clock_t start = clock(), Current;
    GROUPED_Tree = new TTree("GROUPED_Tree", "GROUPED_Tree");
    GROUPED_Tree->Branch("GROUPED_Tree_Silicon", &GROUPED_Tree_Silicon);
    GROUPED_Tree->Branch("GROUPED_Tree_SiPMHigh", &GROUPED_Tree_SiPMHigh);
    GROUPED_Tree->Branch("GROUPED_Tree_SiPMLow", &GROUPED_Tree_SiPMLow);
    ULong64_t TotalEntries = Reader->GetEntries();
    while (Reader->Next())
    {
        ProgressBar(Reader->GetCurrentEntry(), TotalEntries, start, Current, "Selecting Groups : ");
        
        // for (int i = 0; i < signals.GetSize(); i++)
        //     Verbose(signals[i], VERBOSE, 2);

        SearchForCoincidence(signals);    
    }

    WriteHistograms_Grouped();
    /////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  CUTTING GROUPS /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    start = clock(), Current;

    CUTTED_Tree = new TTree("CUTTED_Tree", "CUTTED_Tree");
    CUTTED_Tree->Branch("CUTTED_Tree_Silicon", &CUTTED_Tree_Silicon);
    CUTTED_Tree->Branch("CUTTED_Tree_SiPMHigh", &CUTTED_Tree_SiPMHigh);
    CUTTED_Tree->Branch("CUTTED_Tree_SiPMLow", &CUTTED_Tree_SiPMLow);
    
    Reader = new TTreeReader(GROUPED_Tree);
    Silicon = new TTreeReaderArray<Signal>(*Reader, "GROUPED_Tree_Silicon");
    SiPM_High = new TTreeReaderArray<Signal>(*Reader, "GROUPED_Tree_SiPMHigh");
    SiPM_Low = new TTreeReaderArray<Signal>(*Reader, "GROUPED_Tree_SiPMLow");

    Reader->Restart();
    TotalEntries = Reader->GetEntries();
    while (Reader->Next())
    {
        ProgressBar(Reader->GetCurrentEntry(), TotalEntries, start, Current, "Cutting Groups : ");
        
        // for (int i = 0; i < signals.GetSize(); i++)
        //     Verbose(signals[i], VERBOSE, 2);

        CuttingGroups();    
    }

    WriteHistograms_Cutted();
    /////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////// SILICON WALK GROUPS /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    start = clock(), Current;

    SILICON_WALK_Tree = new TTree("SILICON_WALK_Tree", "SILICON_WALK_Tree");
    SILICON_WALK_Tree->Branch("SILICON_WALK_Tree_Silicon", &SILICON_WALK_Tree_Silicon);
    SILICON_WALK_Tree->Branch("SILICON_WALK_Tree_SiPMHigh", &SILICON_WALK_Tree_SiPMHigh);
    SILICON_WALK_Tree->Branch("SILICON_WALK_Tree_SiPMLow", &SILICON_WALK_Tree_SiPMLow);

    Reader = new TTreeReader(CUTTED_Tree);
    Silicon = new TTreeReaderArray<Signal>(*Reader, "CUTTED_Tree_Silicon");
    SiPM_High = new TTreeReaderArray<Signal>(*Reader, "CUTTED_Tree_SiPMHigh");
    SiPM_Low = new TTreeReaderArray<Signal>(*Reader, "CUTTED_Tree_SiPMLow");
       
    Reader->Restart();
    TotalEntries = Reader->GetEntries();
    while (Reader->Next())
    {
        ProgressBar(Reader->GetCurrentEntry(), TotalEntries, start, Current, "Removing Silicon walk : ");
        
        // for (int i = 0; i < signals.GetSize(); i++)
        //     Verbose(signals[i], VERBOSE, 2);

        SiliconWalkCorrection();
    }

    WriteHistograms_SiliconWalk();
    /////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////// SiPM WALK GROUPS ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    start = clock(), Current;
    Reader->Restart();
    while (Reader->Next())
    {
        ProgressBar(Reader->GetCurrentEntry(), TotalEntries, start, Current, "Removing SiPM walk : ");

        // for (int i = 0; i < signals.GetSize(); i++)
        //     Verbose(signals[i], VERBOSE, 2);

        SiPMWalkCorrection();
    }

    WriteHistograms_SiPMWalk();
    /////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  FINAL CLEANING /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    start = clock(), Current;
    CLEANED_Tree = new TTree("CLEANED_Tree", "CLEANED_Tree");
    CLEANED_Tree->Branch("CLEANED_Tree_Silicon", &CLEANED_Tree_Silicon);
    CLEANED_Tree->Branch("CLEANED_Tree_SiPMHigh", &CLEANED_Tree_SiPMHigh);
    CLEANED_Tree->Branch("CLEANED_Tree_SiPMLow", &CLEANED_Tree_SiPMLow);

    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            CLEANED_Tree_detector[i] = new TTree(("CLEANED_Tree_" + detectorName[i]).c_str(), ("CLEANED_Tree_" + detectorName[i]).c_str());
            CLEANED_Tree_detector[i]->Branch("Channel", &Tree_Channel_detector);
        }
    }

    Reader->Restart();
    while (Reader->Next())
    {
        ProgressBar(Reader->GetCurrentEntry(), TotalEntries, start, Current, "Final Cleaning : ");

        // for (int i = 0; i < signals.GetSize(); i++)
        //     Verbose(signals[i], VERBOSE, 2);

        CleaningGroups();
    }

    WriteHistograms_Cleaned();

    //////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Counting IAS losses /////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    WriteIASLosses();


    delete GROUPED_Tree;
    delete CUTTED_Tree;
    delete SILICON_WALK_Tree;
    WriteTree_Grouped();
    
    GROUPED_File->Close();
    ROOT_File->Close();
    Success("Grouped File Created");
}
