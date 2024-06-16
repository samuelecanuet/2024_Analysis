#include "Grouper.hh"
#include <ctime>

int main(int argc, char *argv[])
{
    clock_t start = clock(), Current;
    string Run_string;

    if (argc == 1)
    {
        Error("No Run Number Given");
        exit(0);
    }
    else
    {
        Run = atoi(argv[1]);
        if (Run < 100)
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

    // ///////////////////////////////////  DELETE OLD FILE ///////////////////////////////////
    // int status = remove((dirNameGrouped+baseFileName+"_grouped.root").c_str());
    // status = remove((dirNameCleaned+baseFileName+"_cleaned.root").c_str());

    ///////////////////////////////////  Grouped ///////////////////////////////////
    GROUPED_File = new TFile((DIR_ROOT_DATA_GROUPED+ROOT_basefilename+"_grouped.root").c_str(), "RECREATE");
    
    ///////////////////////////////////  INITIALISATION ///////////////////////////////////
    InitDetectors("Config_Files/sample.pid");

    InitHistograms_Grouped();
    InitTree_Grouped();
    
    ///////////////////////////////////  PROCESSING TREE ///////////////////////////////////
    Verbose = 0;
    ULong64_t TotalEntries = Tree->GetEntries();
    while (Reader->Next())
    {
        ProgressBar(Reader->GetCurrentEntry(), TotalEntries, start, Current, "Selecting Groups : ");
        
        if (Verbose > 0)
        {
            for (int i = 0; i < signals.GetSize(); i++)
            {
                cout << signals[i] << endl;
            }
        }

        SearchForCoincidence(signals);    
    }

    Info("Writting Histograms");

    ProcessCuttingGroups();
    WriteHistograms_Grouped();


    Reader->Restart();
    while (Reader->Next())
    {
        ProgressBar(Reader->GetCurrentEntry(), TotalEntries, start, Current, "Cutting Groups : ");
        
        if (Verbose > 0)
        {
            for (int i = 0; i < signals.GetSize(); i++)
                cout << signals[i] << endl;
        }
        CuttingGroups(signals);    
    }

    Info("Writting Histograms");
    WriteHistograms_Cutted();
    

    // WriteTree_Grouped();
    // WriteTime(ROOT_File, GROUPED_File);
    GROUPED_File->Close();
    ROOT_File->Close();
    cout << "Grouped File Created" << endl;


}
