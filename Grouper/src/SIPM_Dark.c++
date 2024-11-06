#include "SIPM_Dark.hh"

int main()
{

    InitDetectors("Config_Files/sample.pid");
    InitHistograms();
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-20);

    for (string run : Runs)
    {
        ReadThreshold(run);
    }
    // processing data for each run and sipms
    // for (string run : Runs)
    // {
    //     Info("Reading run " + run);
    //     string ROOTFILE_Name = DIR_DATA_HDD + "run_" + run + "_lossless_NoBeam.fast/run_" + run + "_lossless_NoBeam.root";

    //     TFile *file = new TFile(ROOTFILE_Name.c_str(), "READ");
    //     if (file->IsZombie())
    //     {
    //         Error("Could not open file " + ROOTFILE_Name);
    //     }

    //     TTree *tree = (TTree *)file->Get("Tree");
    //     TTreeReader *Reader = new TTreeReader(tree);
    //     TTreeReaderValue<Signal> SignalReader(*Reader, "Signal");

    //     int Entries = Reader->GetEntries();
    //     clock_t start = clock(), Current;
    //     while (Reader->Next() && Reader->GetCurrentEntry() < 1e8)
    //     {
    //         ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "");
    //         if (IsDetectorBetaHigh((*SignalReader).Label))
    //             H_High[run][GetDetectorChannel((*SignalReader).Label)]->Fill((*SignalReader).Channel);
    //     }

    //     file->Close();
    // }

    // WriteInFile();
    ReadInFile();

    WriteHistograms();
    return 0;
}