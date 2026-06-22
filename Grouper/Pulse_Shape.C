#include "include/Detectors.hh"
#include "include/Utilities.hh"
#include <iostream>
#include <cmath>

TF1 *F_Calibration[SIGNAL_MAX];

void InitCalibration()
{
    TFile *f = MyTFile((DIR_ROOT_DATA_CALIBRATED + "Calibrated_"+ to_string(YEAR) + ".root").c_str(), "READ");
    for (int det = 1;  det < SIGNAL_MAX; det++)
    {
        F_Calibration[det] = (TF1*)f->Get(("Calibration_" + detectorName[det]).c_str());
    }
    f->Close();
}

int Pulse_Shape()
{

    FLAG2025 = true;
    InitDetectors("Config_Files/sample.pid");
    InitWindows();
    InitCalibration();

    TFile *fout = MyTFile("/mnt/hgfs/shared-2/2025_DATA/DETECTOR_DATA/GROUPED/Pulse_Shape.root", "RECREATE");

    map<string, string> NucleiFiles = {
        {"32Ar", "run_000_data_32Ar_grouped_full.root"},
    };
    
    map<string, TH2D*> H_Energy_DT;
    H_Energy_DT["32Ar"] = new TH2D("H_Energy_DT_32Ar", "H_Energy_DT_32Ar", eSiliN_cal, eSiliMin_cal, eSiliMax_cal, 500, 6000, 10000);

    for (auto const& [Nucleus, filename] : NucleiFiles)
    {
        Info("Current Run : " + Nucleus);
        TFile *ROOT_File = MyTFile((DIR_ROOT_DATA_GROUPED + filename).c_str(), "READ");
        TTree *Tree = (TTree *)ROOT_File->Get("CLEANED_Tree");
        TTreeReader *Reader = new TTreeReader(Tree);
        TTreeReaderArray<Signal> signals(*Reader, "CLEANED_Tree_Silicon");

        clock_t start = clock(), Current;   
        int Entries = Reader->GetEntries();
        while (Reader->Next())
        {
            ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "Filling Release Histograms : ");
            
            for (int index = 0; index < signals.GetSize(); index++)
            {
                int current_label = (signals)[index].Label;
                if (IsDetectorSiliStrip(current_label) && current_label > 50)
                {
                    double Energy = F_Calibration[current_label]->Eval((signals)[index].Channel/1000);
                    double DT = (signals)[index].dt;
                    H_Energy_DT["32Ar"]->Fill(Energy, DT);
                }
            }
        }
    }

    fout->cd();
    H_Energy_DT["32Ar"]->Write();
    fout->Close();  

    return 0;
}