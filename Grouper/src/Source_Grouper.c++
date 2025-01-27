#include "Source_Grouper.hh"
#include <ctime>

int main(int argc, char *argv[])
{

    string Run_string;

    if (argc < 2)
    {
        Error("No Run Number Given");
    }
    else
    {
        // option
        if (argc == 3)
        {
            if (string(argv[2]) == "full")
            {
                FULL = true;
            }
        }


        // Run number
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

    TFile *ROOT_File = new TFile((DIR_ROOT_DATA + ROOT_filename).c_str(), "READ");
    TTree *Tree = (TTree *)ROOT_File->Get("Tree_Group");
    
    GROUPED_File = new TFile((DIR_ROOT_DATA_GROUPED + ROOT_basefilename + "_grouped.root").c_str(), "RECREATE");
    GROUPED_File->cd();
    WriteTime(ROOT_File, GROUPED_File);
    ///////////////////////////////////  INITIALISATION ///////////////////////////////////
    InitDetectors("Config_Files/sample.pid");
    InitHistograms_Grouped();
    InitCalibration();

    clock_t start = clock(), Current;
    LoadFitParameters();
    
    ///////////////////////////////////////////////////////////
    start = clock(), Current;
    CLEANED_Tree = new TTree("CLEANED_Tree", "CLEANED_Tree");
    CLEANED_Tree->Branch("CLEANED_Tree_SiPMGroup", &CLEANED_Tree_SiPMGroup);

    TTreeReader *Reader = new TTreeReader(Tree);
    TTreeReaderArray<Signal> signals(*Reader, "Signal");
    Reader->Restart();
    int TotalEntries = 10e6;//Reader->GetEntries();
    while (Reader->Next() && Reader->GetCurrentEntry() < TotalEntries)
    {
        ProgressBar(Reader->GetCurrentEntry(), TotalEntries, start, Current, "Final Cleaning : ");

        // for (int i = 0; i < signals.GetSize(); i++)
        //     cout << signals[i] << endl;

        CleaningGroups(signals);
    }

    WriteHistograms_Cleaned();

    //////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Counting IAS losses /////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    WriteTree_Grouped();

    GROUPED_File->Close();
    ROOT_File->Close();
    Success("Grouped File Created");
}
