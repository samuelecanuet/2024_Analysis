#include "Sensitivity.hh"


int main()
{
    FLAG2025 = true;
    NUCLEUS = "32Ar";
    InitDetectors("../Grouper/Config_Files/sample.pid");
    InitWindows();



    TFile *fout = MyTFile(("/run/media/local1/DATANEX/Samuel-G4/Systematics/Sensibility_B/Sensitivity_B_" + to_string(YEAR) + "_a0.9.root").c_str(),"RECREATE");

    
    
    vector<string> filenames;
    Search_files("/run/media/local1/DATANEX/Samuel-G4/Systematics/Sensibility_B/", {"a0.9", "_result.root"}, filenames);

    map<double, TGraphErrors*> G_Eshift;


    for (auto &filename : filenames)
    {
        TFile *file = MyTFile(filename.c_str(), "READ");
        if (!file || file->IsZombie())
        {
            Warning("Could not open file: " + filename);
            continue;
        }   

        cout << filename.substr(filename.find("5_B") + 3, filename.find("T") - (filename.find("5_B") + 3)) << endl;

        double Bfield = stod(filename.substr(filename.find("5_B") + 3, filename.find("T") - (filename.find("5_B") + 3)));

        for (int det = 1; det < SIGNAL_MAX; det++)
        {
            if (IsDetectorSiliStrip(det))
            {
                pair<double, double> Eshift = ComputeEshift(file, det);
                if (G_Eshift[Bfield] == nullptr)
                {
                    G_Eshift[Bfield] = new TGraphErrors();
                    G_Eshift[Bfield]->SetName(("G_Eshift_Bfield_" + to_string(Bfield) + "_" + detectorName[det]).c_str());
                }
                // adding data to graph for each det
                G_Eshift[Bfield]->AddPoint(det, Eshift.first);
                G_Eshift[Bfield]->SetPointError(G_Eshift[Bfield]->GetN() - 1, 0, Eshift.second);
            }
        }

        file->Close();
    }

    fout ->cd();
    TCanvas *c = new TCanvas("c", "c", 800, 600);
    TGraphErrors *g_mean = new TGraphErrors();
    for (auto &entry : G_Eshift)
    {
        if (G_Eshift[entry.first] != nullptr)
        {
            G_Eshift[entry.first]->Fit("pol0", "QE");
            g_mean->AddPoint(entry.first, G_Eshift[entry.first]->GetFunction("pol0")->GetParameter(0));
            g_mean->SetPointError(g_mean->GetN() - 1, 0, G_Eshift[entry.first]->GetFunction("pol0")->GetParError(0));

            G_Eshift[entry.first]->SetTitle(("Eshift for Bfield = " + to_string(entry.first) + " T").c_str());
            G_Eshift[entry.first]->GetXaxis()->SetTitle("Detector");
            G_Eshift[entry.first]->GetYaxis()->SetTitle("Eshift [keV]");
            G_Eshift[entry.first]->Write(); 
        }
    }
    g_mean->Draw("AP");
    c->Write();
    fout->Close();  



    return 0;
}