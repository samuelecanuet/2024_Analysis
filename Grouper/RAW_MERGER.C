#include "include/Detectors.hh"

int RAW_MERGER()
{

    FLAG2024 = true;
    InitDetectors("Config_Files/sample.pid");
    InitRuns();

    TFile *fout = MyTFile(DIR_ROOT_DATA + "run_000_data_32Ar.root", "RECREATE");

    vector<string> filenames;
    for (int i = 0; i < Map_RunFiles["32Ar"].size(); i++)
    {
        string Run = Map_RunFiles["32Ar"][i];
        int Run_int = atoi(Run.c_str());

        string filename = SearchFiles(DIR_ROOT_DATA, Run);

        if (filename.empty())
        {
            Warning("No file found for Run " + Run);
            continue;
        }
        filenames.push_back(filename);
    }

    // TChain merging all the files of the same nucleus
    TChain *chain = new TChain("Tree_Group");
    for (string filename : filenames)
    {
        chain->Add((DIR_ROOT_DATA + filename).c_str());
    }

    TTree *merged_tree = chain->CloneTree(0);
    merged_tree->CopyEntries(chain);
    fout->cd();
    merged_tree->Write("Tree_Group");
    fout->Close();

    return 1;
}