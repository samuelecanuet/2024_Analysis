#include "BasicAnalyser.hh"



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
    string ROOT_filename = SearchFiles(DIR_ROOT_DATA, Run_string);
    string ROOT_basefilename = ROOT_filename.substr(0, ROOT_filename.find(".root"));

    TFile *ROOT_File = new TFile((DIR_ROOT_DATA + ROOT_filename).c_str(), "READ");
    TTree *Tree = (TTree *)ROOT_File->Get("Tree_Group");
    TTreeReader *Reader = new TTreeReader(Tree);
    TTreeReaderArray<Signal> signals(*Reader, "Signal");

    InitDetectors("Config_Files/sample.pid");
    

    GROUPED_File = new TFile((DIR_ROOT_DATA + ROOT_basefilename + "_basic.root").c_str(), "RECREATE");
    InitHistograms();

    clock_t start = clock(), Current;
    CLEANED_Tree = new TTree("CLEANED_Tree", "CLEANED_Tree");
    CLEANED_Tree->Branch("CLEANED_Tree_Silicon", &CLEANED_Tree_Silicon);
    CLEANED_Tree->Branch("CLEANED_Tree_SiPMGroup", &CLEANED_Tree_SiPMGroup);

    int TotalEntries = Reader->GetEntries();
    while (Reader->Next())
    {
        ProgressBar(Reader->GetCurrentEntry(), TotalEntries, start, Current, "Final Cleaning : ");

        // for (int i = 0; i < signals.GetSize(); i++)
        //     cout << signals[i] << endl;

        CleaningGroups(signals, 0);
    }

    WriteHistograms();
    GROUPED_File->Close();
    ROOT_File->Close();

}