#include "Reader.hh"

int main()
{
    string name = "30-09/32Ar_100_6000_100_CS0_CSP0_CV1_CVP1";
    // name = "241Am_700nm_width";
    TFile *SIMULATED_File = new TFile((DIR_ROOT_DATA_SIMULATED + name + ".root").c_str(), "READ");
    ANALYSIS_File = new TFile((DIR_ROOT_DATA_SIMULATED + name + "_analysed.root").c_str(), "RECREATE");
    TTree *tree = (TTree *)SIMULATED_File->Get("Tree");

    Reader = new TTreeReader(tree);
    Tree_PDG = new TTreeReaderArray<int>(*Reader, "Particle_PDG");
    Tree_E0 = new TTreeReaderArray<double>(*Reader, "Kinetic_Energy");
    Tree_x = new TTreeReaderArray<double>(*Reader, "x");
    Tree_y = new TTreeReaderArray<double>(*Reader, "y");
    Tree_z = new TTreeReaderArray<double>(*Reader, "z");
    Tree_Catcher_Central_Energy_Deposit = new TTreeReaderArray<double>(*Reader, "Catcher_Central_Energy_Deposit");
    Tree_Catcher_Side_Energy_Deposit = new TTreeReaderArray<double>(*Reader, "Catcher_Side_Energy_Deposit");
    Tree_PlasticScintillator_Energy_Deposit = new TTreeReaderArray<double>(*Reader, "PlasticScintillator_Energy_Deposit");
    Tree_PlasticScintillator_Hit_Position = new TTreeReaderArray<Hep3Vector>(*Reader, "PlasticScintillator_Hit_Position");
    Tree_PlasticScintillator_Hit_Angle = new TTreeReaderArray<double>(*Reader, "PlasticScintillator_Hit_Angle");
    Tree_Silicon_Detector_Energy_Deposit = new TTreeReaderArray<vector<double>>(*Reader, "Silicon_Detector_Energy_Deposit");
    Tree_Silicon_Detector_Hit_Position = new TTreeReaderArray<vector<Hep3Vector>>(*Reader, "Silicon_Detector_Hit_Position");
    Tree_Silicon_Detector_Hit_Angle = new TTreeReaderArray<vector<double>>(*Reader, "Silicon_Detector_Hit_Angle");
    Tree_Silicon_Detector_Code = new TTreeReaderArray<vector<int>>(*Reader, "Silicon_Detector_Code");
    Tree_Silicon_Detector_DL_Energy_Deposit = new TTreeReaderArray<vector<double>>(*Reader, "Silicon_Detector_DL_Energy_Deposit");
    Tree_Silicon_Detector_InterStrip_Code = new TTreeReaderArray<vector<int>>(*Reader, "Silicon_Detector_InterStrip_Code");
    tree->SetBranchAddress("Silicon_Detector_InterStrip_Energy_Deposit", &Tree_Silicon_Detector_InterStrip_Energy_Deposit);
    tree->SetBranchAddress("Silicon_Detector_InterStrip_Hit_Position", &Tree_Silicon_Detector_InterStrip_Hit_Position);

    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            OutTree[i] = new TTree(("Tree_" + detectorName[i]).c_str(), ("Tree_" + detectorName[i]).c_str());
            OutTree[i]->Branch("Energy", &e);
        }
    }

    int Entries = tree->GetEntries();
    clock_t start = clock(), Current;

    Init();
    InitDetectors("../Grouper/Config_Files/sample.pid");
    InitHistograms();

    Info("Starting Loop");
    while (Reader->Next())
    {
        ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "Reading Tree");

        //// LOOP ON PARTICLE ////
        for (int part_i = 0; part_i < Tree_PDG->GetSize(); part_i++)
        {
            int PDG = Tree_PDG->At(part_i);
            double E0 = Tree_E0->At(part_i);
            double x = Tree_x->At(part_i);
            double y = Tree_y->At(part_i);
            double z = Tree_z->At(part_i);
            double Catcher_Central_Energy_Deposit = Tree_Catcher_Central_Energy_Deposit->At(part_i);
            double Catcher_Side_Energy_Deposit = Tree_Catcher_Side_Energy_Deposit->At(part_i);

            int index = PDGtoIndex[PDG];

            H_E0[index]->Fill(E0);
            H_x[index]->Fill(x);
            H_y[index]->Fill(y);
            H_z[index]->Fill(z);
            H_Catcher_Central_Energy_Deposit[index]->Fill(Catcher_Central_Energy_Deposit);
            H_Catcher_Side_Energy_Deposit[index]->Fill(Catcher_Side_Energy_Deposit);
            if (Tree_PlasticScintillator_Energy_Deposit->At(part_i) != 0)
            {
                H_PlasticScintillator_Energy_Deposit[index]->Fill(Tree_PlasticScintillator_Energy_Deposit->At(part_i));
                H_PlasticScintillator_Hit_Position_x[index]->Fill(Tree_PlasticScintillator_Hit_Position->At(part_i).x());
                H_PlasticScintillator_Hit_Position_y[index]->Fill(Tree_PlasticScintillator_Hit_Position->At(part_i).y());
                H_PlasticScintillator_Hit_Position_z[index]->Fill(Tree_PlasticScintillator_Hit_Position->At(part_i).z());
                H_PlasticScintillator_Hit_Position_xy[index]->Fill(Tree_PlasticScintillator_Hit_Position->At(part_i).x(), Tree_PlasticScintillator_Hit_Position->At(part_i).y());
                H_PlasticScintillator_Hit_Angle[index]->Fill(Tree_PlasticScintillator_Hit_Angle->At(part_i));
            }


            //// LOOP ON HITTED SILICON DETECTORS ////
            for (int sili_i = 0; sili_i < Tree_Silicon_Detector_Energy_Deposit->At(part_i).size(); sili_i++)
            {
                double X = Tree_Silicon_Detector_Hit_Position->At(part_i).at(sili_i).x();
                double Y = Tree_Silicon_Detector_Hit_Position->At(part_i).at(sili_i).y();
                double Z = Tree_Silicon_Detector_Hit_Position->At(part_i).at(sili_i).z();    
                double THETA = Tree_Silicon_Detector_Hit_Angle->At(part_i).at(sili_i);
                int Silicon_Detector_Code = Tree_Silicon_Detector_Code->At(part_i).at(sili_i);
                H_Silicon_Detector_Hit_Position_x[index]->Fill(X);
                H_Silicon_Detector_Hit_Position_y[index]->Fill(Y);
                H_Silicon_Detector_Hit_Position_z[index]->Fill(Z);
                H_Silicon_Detector_Hit_Position_xy[index]->Fill(X, Y);
                H_Silicon_Detector_Hit_Position_xz[index]->Fill(X, Z);
                H_Silicon_Detector_Hit_Position_yz[index]->Fill(Y, Z);
                H_Silicon_Detector_Hit_Angle[index]->Fill(THETA);
                H_Silicon_Detector_Hit_Anglez[index]->Fill(THETA, Z);
                H_Silicon_Detector_Hit_Anglexy[index]->Fill(THETA, X);
                H_Silicon_Detector_Code[index]->Fill(Silicon_Detector_Code);

                H_Silicon_Detector_Energy_Deposit_Det[index][Silicon_Detector_Code]->Fill(Tree_Silicon_Detector_Energy_Deposit->At(part_i).at(sili_i));
                Full_energy[Silicon_Detector_Code] += Tree_Silicon_Detector_Energy_Deposit->At(part_i).at(sili_i);
                H_Silicon_Detector_Hit_Angle_Det[index][Silicon_Detector_Code]->Fill(THETA);
            }

            for (int sili_inter_i = 0; sili_inter_i < Tree_Silicon_Detector_InterStrip_Code->At(part_i).size(); sili_inter_i++)
            {
                int Silicon_Detector_InterStrip_Code = Tree_Silicon_Detector_InterStrip_Code->At(part_i).at(sili_inter_i);
                double energy_interstrip = 0;
                int code_strip_for = (Silicon_Detector_InterStrip_Code+50)/100;

                tree->GetEntry(Reader->GetCurrentEntry());
                for (int step = 0; step < Tree_Silicon_Detector_InterStrip_Energy_Deposit->at(part_i).at(sili_inter_i).size(); step++)
                {
                    double X = Tree_Silicon_Detector_InterStrip_Hit_Position->at(part_i).at(sili_inter_i).at(step).x();
                    double Y = Tree_Silicon_Detector_InterStrip_Hit_Position->at(part_i).at(sili_inter_i).at(step).y();
                    double Z = Tree_Silicon_Detector_InterStrip_Hit_Position->at(part_i).at(sili_inter_i).at(step).z();
                    H_Silicon_Detector_InterStrip_xy[index]->Fill(X, Y);
                    H_Silicon_Detector_InterStrip_z[index]->Fill(Z);
                    energy_interstrip += Tree_Silicon_Detector_InterStrip_Energy_Deposit->at(part_i).at(sili_inter_i).at(step);

                    // giving energy to strips

                    // if (Z + (300 * 1e-3 ) / 2 > 300 * 1e-3 - 800 * 1e-6)
                        // continue;
                    // new 
                    if (Z + (300 * 1e-3 + 800 * 1e-6) / 2 > 300 * 1e-3)
                        continue;
                    double ratio = 1;//-abs(1/(9000*x)) + 1;
                    // if (ratio < 0.5)
                        // ratio = 0.5;
                    
                    if ( Y > 0)
                    {
                        InterStrip[code_strip_for] += ratio * Tree_Silicon_Detector_InterStrip_Energy_Deposit->at(part_i).at(sili_inter_i).at(step);
                        InterStrip[code_strip_for-1] += (1-ratio) * Tree_Silicon_Detector_InterStrip_Energy_Deposit->at(part_i).at(sili_inter_i).at(step);
                        Full_energy[code_strip_for] += ratio * Tree_Silicon_Detector_InterStrip_Energy_Deposit->at(part_i).at(sili_inter_i).at(step);
                        Full_energy[code_strip_for-1] += (1-ratio) * Tree_Silicon_Detector_InterStrip_Energy_Deposit->at(part_i).at(sili_inter_i).at(step);
                    }
                    else
                    {
                        InterStrip[code_strip_for] += (1-ratio) * Tree_Silicon_Detector_InterStrip_Energy_Deposit->at(part_i).at(sili_inter_i).at(step);
                        InterStrip[code_strip_for-1] += ratio * Tree_Silicon_Detector_InterStrip_Energy_Deposit->at(part_i).at(sili_inter_i).at(step);
                        Full_energy[code_strip_for] += (1-ratio) * Tree_Silicon_Detector_InterStrip_Energy_Deposit->at(part_i).at(sili_inter_i).at(step);
                        Full_energy[code_strip_for-1] += ratio * Tree_Silicon_Detector_InterStrip_Energy_Deposit->at(part_i).at(sili_inter_i).at(step);
                    }
                    // double ratio = Y_real / size_interstrip;
                    
                    // cout << "Code: " << Silicon_Detector_InterStrip_Code << " Strip: " << code_strip_for << " Ratio: " << ratio << endl;
                    
                    
                }
                // cout << energy << endl;

                //energy deposit in the interstrip for each particle
                if (energy_interstrip != 0)
                    H_Silicon_Detector_InterStrip_Energy_Deposit_Det[index][Silicon_Detector_InterStrip_Code]->Fill(energy_interstrip);
            }
        }

        // energy despoit in the interstrip for all the particle
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorSiliStrip(i))
            {
                if (InterStrip[i] != 0)
                {
                    H_Silicon_Detector_Energy_Deposit_Det[Particle_Size - 2][i]->Fill(InterStrip[i]);
                    InterStrip[i] = 0;
                }
            }
        }


        // counting if x strips triggered
        int counter_strip_trigger = 0;
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorSiliStrip(i))
            {
                if (Full_energy[i] > 400)
                {
                    counter_strip_trigger++;
                }
            }
        }

        // filling energy deposit by all particle 
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorSiliStrip(i) && counter_strip_trigger == 1)
            {
                if (Full_energy[i] > 400)
                {
                    H_Silicon_Detector_Energy_Deposit_Det[Particle_Size-1][i]->Fill(Full_energy[i]);    
                    e = Full_energy[i];
                    OutTree[i]->Fill();
                }
            }
            Full_energy[i] = 0;
                e = 0;
        }
        

        //// LOOP ON PARTICLE ////
        if (Verbosee == 1)
        {
            Info("Event: " + to_string((int)Reader->GetCurrentEntry()));
            for (int i = 0; i < Tree_PDG->GetSize(); i++)
            {
                cout << " -| Particle: " << SearchName(Tree_PDG->At(i)) << endl;
                cout << " ----| E0: " << Tree_E0->At(i) << " keV" << endl;
                cout << " ----| x: " << Tree_x->At(i) << " um  y: " << Tree_y->At(i) << " um  z: " << Tree_z->At(i) << " nm" << endl;
                cout << " ----| Catcher: " << Tree_Catcher_Central_Energy_Deposit->At(i) << " keV" << endl;

                if (Tree_PlasticScintillator_Energy_Deposit->At(i) != 0)
                {
                    cout << " ----| Plastic Scintillator: " << Tree_PlasticScintillator_Energy_Deposit->At(i) << " keV" << endl;
                    cout << " ----| Plastic Scintillator Hit Angle: " << Tree_PlasticScintillator_Hit_Angle->At(i) << " deg" << endl;
                }
            }
        }
    }

    WriteHistograms();
    ANALYSIS_File->cd();
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            OutTree[i]->Write();
        }
    }
    SIMULATED_File->Close();
    ANALYSIS_File->Close();
}