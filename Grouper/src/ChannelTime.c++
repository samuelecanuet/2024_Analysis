
#include "ChannelTime.hh"

int main()
{

    TFile *f = new TFile("../../../../../mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/GROUPED/run_114_multifast_32Ar_grouped.root", "READ");
    TTree *Tree = (TTree *)f->Get("Tree_Group");

    //verify if file is open
    if (!f->IsOpen())
    {
        cout << "Error: File not opened" << endl;
        return 0;
    }

    // check if tree exists
    if (!f->GetListOfKeys()->Contains("CLEANED_Tree"))
    {
        cout << "Error: Tree not found" << endl;
        return 0;
    }
     
    TTreeReader *Reader = new TTreeReader(Tree);
    Signal_Silicon = new TTreeReaderArray<Signal>(*Reader, "Signal");
    int counter = 0;
    double t0;
    TH2D* H = new TH2D("H", "H", 7200, 0, 7200, 10000, 0, 100000);
    int Entries = Reader->GetEntries();
    clock_t start = clock(), Current;
    while (Reader->Next())
    {
        ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "");
        // cout << (*Signal_Silicon)[1] << endl;
        // if ((*Signal_Silicon)[1].Label == 11)
        // {
        //     if (counter == 0)
        //         t0 = (*Signal_Silicon)[1].Time*1e-9;
        //     H->Fill((*Signal_Silicon)[1].Time*1e-9-t0, (*Signal_Silicon)[1].Channel);
        // }   
        // counter++;

        for (int i = 0; Signal_Silicon->GetSize(); i++)
        {
            if ((*Signal_Silicon)[i].Pileup == 1)
            {
                cout << "PileUp" << endl;
            }
        }
    }
    
    // TFile *ff = new TFile ("CHANNELTime.root", "RECREATE");   
    // ff->cd();
    // H->Write();
    // ff->Close();
}