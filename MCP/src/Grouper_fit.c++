#include "Grouper_fit.hh"

int main(int argc, char *argv[])
{
    FLAG2025 = true;
    string year = to_string(YEAR);
    InitDetectors("Config_Files/sample.pid");

    if (argc != 2)
    {
        Error("Usage: ./Grouper <Run>");
    }

    ///////////////////////////////////  SWIPPING ///////////////////////////////////
    string Run_string = argv[1];
    ROOT_filename = SearchFiles(DIR_ROOT_DATA_MCP, Run_string);
    ROOT_basefilename = ROOT_filename.substr(0, ROOT_filename.find(".root"));
    TFile *ROOT_File = MyTFile((DIR_ROOT_DATA_MCP + ROOT_filename).c_str(), "READ");


    ///////////////////////////////////  OUTPUT ///////////////////////////////////
    OUTPUT_File = MyTFile((DIR_ROOT_DATA_MCP + ROOT_basefilename + "_grouped_fit.root").c_str(), "RECREATE");
    
    TTree *Tree = (TTree *)ROOT_File->Get("Tree_Group");
    TTreeReader *Reader = new TTreeReader(Tree);
    TTreeReaderArray<Signal> signals(*Reader, "Signal");


    InitHistograms();

    clock_t start = clock(), Current;
    int Entries = Reader->GetEntries();
    while (Reader->Next())
    {
        ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "Final Cleaning : ");

        // cout << "#### " << Reader->GetCurrentEntry() << " ####" << endl;
        ReadTree_MCP(signals);
        // cout << endl;
    }   

        FitTry();


    WriteHistograms();

    ROOT_File->Close();
}