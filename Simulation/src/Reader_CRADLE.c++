#include "Reader_CRADLE.hh"

int main()
{
    string name = "../../../../../../../mnt/hgfs/shared-2/2024_DATA/SIMULATED_DATA/fe/fe/test";
    // name = "241Am_700nm_width";
    TFile *SIMULATED_File = new TFile((DIR_ROOT_DATA_SIMULATED + name + ".root").c_str(), "READ");
    ANALYSIS_File = new TFile((DIR_ROOT_DATA_SIMULATED + name + "_analysed.root").c_str(), "RECREATE");
    TTree *tree = (TTree *)SIMULATED_File->Get("ParticleTree");

    Reader = new TTreeReader(tree);
    Tree_PDG = new TTreeReaderArray<int>(*Reader, "code");
    Tree_E0 = new TTreeReaderArray<double>(*Reader, "energy");
    Tree_time = new TTreeReaderArray<double>(*Reader, "time");
    Tree_Ex = new TTreeReaderArray<double>(*Reader, "excitation_energy");
    Tree_px = new TTreeReaderArray<double>(*Reader, "px");
    Tree_py = new TTreeReaderArray<double>(*Reader, "py");
    Tree_pz = new TTreeReaderArray<double>(*Reader, "pz");


    int Entries = tree->GetEntries();
    clock_t start = clock(), Current;

    InitHistograms();

    Info("Starting Loop");
    int Event_MAX = 1e6;

    while (Reader->Next() && Reader->GetCurrentEntry() < Event_MAX)
    {
        ProgressBar(Reader->GetCurrentEntry(), Event_MAX, start, Current, "Reading Tree");

        // cout << "Event: " << **Tree_Event << endl;
        for (int part_i = 0; part_i < Tree_PDG->GetSize(); part_i++)
        {
            VerifyPDG(Tree_PDG->At(part_i));
            int PDG = Tree_PDG->At(part_i);
            //cout << "Particle: " << **Tree_PDG << " Energy: " << **Tree_E0 << " Ex: " << **Tree_Ex << " Time: " << **Tree_time << " px: " << **Tree_px << " py: " << **Tree_py << " pz: " << **Tree_pz << endl;
            H_E0[PDG]->Fill(Tree_E0->At(part_i));
            H_Ex[PDG]->Fill(Tree_Ex->At(part_i));
            H_Time[PDG]->Fill(Tree_time->At(part_i));
            H_px[PDG]->Fill(Tree_px->At(part_i));
            H_py[PDG]->Fill(Tree_py->At(part_i));
            H_pz[PDG]->Fill(Tree_pz->At(part_i));
        }        
    }

    WriteHistograms();

    SIMULATED_File->Close();
    ANALYSIS_File->Close();
}