#include "Reader_CRADLE.hh"

int main()
{
    string name = "../../../../../../../mnt/hgfs/shared-2/2024_DATA/SIMULATED_DATA/fe/fe/test";
    // name = "241Am_700nm_width";
    TFile *SIMULATED_File = new TFile((DIR_ROOT_DATA_SIMULATED + name + ".root").c_str(), "READ");
    ANALYSIS_File = new TFile((DIR_ROOT_DATA_SIMULATED + name + "_analysed.root").c_str(), "RECREATE");
    TTree *tree = (TTree *)SIMULATED_File->Get("ParticleTree");

    Reader = new TTreeReader(tree);
    Tree_Event = new TTreeReaderValue<int>(*Reader, "event");
    Tree_PDG = new TTreeReaderValue<int>(*Reader, "code");
    Tree_E0 = new TTreeReaderValue<double>(*Reader, "energy");
    Tree_time = new TTreeReaderValue<double>(*Reader, "time");
    Tree_Ex = new TTreeReaderValue<double>(*Reader, "excitation_energy");
    Tree_px = new TTreeReaderValue<double>(*Reader, "px");
    Tree_py = new TTreeReaderValue<double>(*Reader, "py");
    Tree_pz = new TTreeReaderValue<double>(*Reader, "pz");


    int Entries = tree->GetEntries();
    clock_t start = clock(), Current;

    InitHistograms();

    Info("Starting Loop");
    int Event_MAX = 1e6;
    int Event = -1;

    vector<double> pe = {0, 0, 0};
    vector<double> e_nu = {0, 0, 0};
    while (Event < Event_MAX)
    {
        ProgressBar(Event, Event_MAX, start, Current, "Reading Tree");
        if (!Reader->Next())
            break;
        Event++;

        // cout << "Event: " << **Tree_Event << endl;
        while (Event == **Tree_Event)
        {
            VerifyPDG(**Tree_PDG);
            //cout << "Particle: " << **Tree_PDG << " Energy: " << **Tree_E0 << " Ex: " << **Tree_Ex << " Time: " << **Tree_time << " px: " << **Tree_px << " py: " << **Tree_py << " pz: " << **Tree_pz << endl;
            H_E0[**Tree_PDG]->Fill(**Tree_E0);
            H_Ex[**Tree_PDG]->Fill(**Tree_Ex);
            H_Time[**Tree_PDG]->Fill(**Tree_time);
            H_px[**Tree_PDG]->Fill(**Tree_px);
            H_py[**Tree_PDG]->Fill(**Tree_py);
            H_pz[**Tree_PDG]->Fill(**Tree_pz);


            Reader->Next();
            if (!Reader->Next())
            break;

        }
    }

    WriteHistograms();

    SIMULATED_File->Close();
    ANALYSIS_File->Close();
}