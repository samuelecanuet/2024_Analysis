#include "SIPM_PileUp.hh"

int main()
{

    InitDetectors("Config_Files/sample.pid");

    string RUN = "077";
     ///////////////////////////////////  FILES //////////////////////////////////
    MERGED_File = MyTFile((DIR_ROOT_DATA + "run_"+RUN+"_multifast_32Ar.root").c_str(), "READ");
    

    ///////////////////////////////////  OUTPUT ///////////////////////////////////
    FINAL_File = MyTFile((DIR_ROOT_DATA_MATCHED + "PileUp_run_"+RUN+".root").c_str(), "RECREATE");

    TTree *Tree = (TTree *)MERGED_File->Get("Tree_Group");
    Reader = new TTreeReader(Tree);
    signals = new TTreeReaderArray<Signal>(*Reader, "Signal");

    InitHistograms();

    clock_t start = clock(), Current;
    int Entries = Reader->GetEntries();
    while (Reader->Next())
    {
        ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "Reading : ");
        for (int i = 0; i < (*signals).GetSize(); i++)
        {
            if (IsDetectorBetaHigh((*signals)[i].Label))
            {
                H_Q1[(*signals)[i].Label]->Fill((*signals)[i].Channel);
                H_Q2[(*signals)[i].Label]->Fill((*signals)[i].Pileup);
                H_Q1_Q2[(*signals)[i].Label]->Fill((*signals)[i].Channel, (*signals)[i].Pileup);
            }

            if (IsDetectorBetaLow((*signals)[i].Label))
            {
                H_Q1[(*signals)[i].Label]->Fill((*signals)[i].Channel);
                H_Q2[(*signals)[i].Label]->Fill((*signals)[i].Pileup);
                H_Q1_Q2[(*signals)[i].Label]->Fill((*signals)[i].Channel, (*signals)[i].Pileup);
                P_Q1_Q2[(*signals)[i].Label]->Fill((*signals)[i].Channel, (*signals)[i].Pileup);
            }
        }
    }

    FittingQ1Q2();

    Reader->Restart();

    while (Reader->Next())
    {
        ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "Reading : ");
        for (int i = 0; i < (*signals).GetSize(); i++)
        {
            if (IsDetectorBetaHigh((*signals)[i].Label))
            {
                if (!PileUp((*signals)[i], 1))
                {
                    H_Q1_Q2_Cleaned[(*signals)[i].Label]->Fill((*signals)[i].Channel, (*signals)[i].Pileup);
                    H_Q1_Cleaned[(*signals)[i].Label]->Fill((*signals)[i].Channel);
                }

                H_Q1_Q2_Modified[(*signals)[i].Label]->Fill((*signals)[i].Pileup - Correction[(*signals)[i].Label]->Eval((*signals)[i].Channel), (*signals)[i].Channel);
            }

            if (IsDetectorBetaLow((*signals)[i].Label))
            {
                if (!PileUp((*signals)[i], 10))
                {
                    H_Q1_Q2_Cleaned[(*signals)[i].Label]->Fill((*signals)[i].Channel, (*signals)[i].Pileup);
                    H_Q1_Cleaned[(*signals)[i].Label]->Fill((*signals)[i].Channel);
                }

                H_Q1_Q2_Modified[(*signals)[i].Label]->Fill((*signals)[i].Pileup - Correction[(*signals)[i].Label]->Eval((*signals)[i].Channel), (*signals)[i].Channel);
            }
        }
    }

    WriteHistograms();


    return 0;
}