#include "ReaderNew.hh"
default_random_engine generator;
int main()
{
    FLAG2025=true;
    NUCLEUS="32Ar";
    InitDetectors("../Grouper/Config_Files/sample.pid");


    // string name = "/run/media/local1/DATANEX/Samuel-G4/32Ar_ION_2025_CAD";
    string name = "/run/media/local1/DATANEX/Samuel-G4/32Ar_IAS_2025_Default";
    
    TFile *SIMULATED_File = MyTFile((name + ".root").c_str(), "READ");
    ANALYSIS_File = MyTFile((name + "_analysed_new.root").c_str(), "RECREATE");
    TTree *tree = (TTree *)SIMULATED_File->Get("Tree");

    Reader = new TTreeReader(tree);
    Tree_PDG = new TTreeReaderArray<int>(*Reader, "Particle_PDG");
    Tree_E0 = new TTreeReaderArray<double>(*Reader, "Kinetic_Energy");
    Tree_x = new TTreeReaderArray<double>(*Reader, "x");
    Tree_y = new TTreeReaderArray<double>(*Reader, "y");
    Tree_z = new TTreeReaderArray<double>(*Reader, "z");
    Tree_px = new TTreeReaderArray<double>(*Reader, "px");
    Tree_py = new TTreeReaderArray<double>(*Reader, "py");
    Tree_pz = new TTreeReaderArray<double>(*Reader, "pz");
    Tree_T0 = new TTreeReaderArray<double>(*Reader, "T0");
    Tree_Catcher_Central_Energy_Deposit = new TTreeReaderArray<double>(*Reader, "Catcher_Central_Energy_Deposit");
    Tree_Catcher_Side_Energy_Deposit = new TTreeReaderArray<double>(*Reader, "Catcher_Side_Energy_Deposit");
    Tree_PlasticScintillator_Energy_Deposit = new TTreeReaderArray<double>(*Reader, "PlasticScintillator_Energy_Deposit");
    Tree_PlasticScintillator_Visible_Energy_Deposit = new TTreeReaderArray<double>(*Reader, "PlasticScintillator_Visible_Energy_Deposit");
    Tree_PlasticScintillator_Hit_Position = new TTreeReaderArray<Hep3Vector>(*Reader, "PlasticScintillator_Hit_Position");
    Tree_PlasticScintillator_Hit_Angle = new TTreeReaderArray<double>(*Reader, "PlasticScintillator_Hit_Angle");
    Tree_PlasticScintillator_Hit_Time = new TTreeReaderArray<double>(*Reader, "PlasticScintillator_Hit_Time");
    Tree_Silicon_Detector_Energy_Deposit = new TTreeReaderArray<vector<double>>(*Reader, "Silicon_Detector_Energy_Deposit");
    Tree_Silicon_Detector_Hit_Position = new TTreeReaderArray<vector<Hep3Vector>>(*Reader, "Silicon_Detector_Hit_Position");
    Tree_Silicon_Detector_Hit_Angle = new TTreeReaderArray<vector<double>>(*Reader, "Silicon_Detector_Hit_Angle");
    Tree_Silicon_Detector_Hit_Time = new TTreeReaderArray<vector<double>>(*Reader, "Silicon_Detector_Hit_Time");
    Tree_Silicon_Detector_Code = new TTreeReaderArray<vector<int>>(*Reader, "Silicon_Detector_Code");
    Tree_Silicon_Detector_DL_Energy_Deposit = new TTreeReaderArray<vector<double>>(*Reader, "Silicon_Detector_DL_Energy_Deposit");


    // Init Tree for each silicon detector
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            OutTree[i] = new TTree(("Tree_" + detectorName[i]).c_str(), ("Tree_" + detectorName[i]).c_str());
            OutTree[i]->Branch("Energy", &e);
        }
    }

    // Init Tree for Plastic Scintillator
    PlasticIASTree = new TTree("PlasticIAS", "PlasticIAS");
    PlasticIASTree->Branch("Code", &sili_code);
    PlasticIASTree->Branch("Energy", &sili_e);
    PlasticIASTree->Branch("SiPM", &SiPM_e);

    // Init final output Tree
    FINAL_OutTree = new TTree("Tree", "Tree");
    FINAL_OutTree->Branch("Silicon", &FINAL_OutTree_Silicon);
    FINAL_OutTree->Branch("SiPM", &FINAL_OutTree_SiPM);

    int Entries = tree->GetEntries();
    clock_t start = clock(), Current;

    Verbosee = 0;

    Init();
    
    InitHistograms(0);
    InitWindows();

    bool FLAGPRINTINGEVENT = false;


    double Time_begining = 0;
   

    Info("Starting Loop");
    while (Reader->Next() && Reader->GetCurrentEntry() < Entries)
    {
        ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "Reading Tree");

        vector<Signal> DetectorSignals_RAW[SIGNAL_MAX];

        if (FLAGPRINTINGEVENT)
        {
            Verbosee = 1;
        }

        if (Verbosee == 1) cout << endl;
        Time_begining = 0;

        //// LOOP ON PARTICLE ////
        for (int part_i = 0; part_i < Tree_PDG->GetSize(); part_i++)
        {
            // if (part_i == 1) Time_begining = Tree_T0->At(part_i);
            // int PDG = Tree_PDG->At(part_i);

            // if (Verbosee == 1) 
            // {
            //     Info("Particle: " + SearchName(PDG), 1);
            //     Info("E0: " + to_string(Tree_E0->At(part_i)) + " keV", 2);
            //     Info("T0: " + to_string((Tree_T0->At(part_i)-Time_begining)/1e9) + " s", 2);
            // }
            // if (PDG > 1000000000) // removing ion excited states information
            //     PDG = (PDG / 10) * 10;
            
            // if (H_E0.find(PDG) == H_E0.end())
            //     InitHistograms(PDG);

            //// INITIAL ////
            // if (Verbosee == 1) Info("Initial state", 2);
            // double E0 = Tree_E0->At(part_i);
            // double T0 = Tree_T0->At(part_i)/1e9;
            // double x = Tree_x->At(part_i)*1e3;
            // double y = Tree_y->At(part_i)*1e3;
            // double z = Tree_z->At(part_i)*1e6;
            // double px = Tree_px->At(part_i);
            // double py = Tree_py->At(part_i);
            // double pz = Tree_pz->At(part_i);
            // H_E0[PDG]->Fill(E0);
            // H_T0[PDG]->Fill(T0);
            // H_x[PDG]->Fill(x);
            // H_y[PDG]->Fill(y);
            // H_z[PDG]->Fill(z);
            // H_px[PDG]->Fill(px);
            // H_py[PDG]->Fill(py);
            // H_pz[PDG]->Fill(pz);

            //// CATCHER ////
            // if (Verbosee == 1) Info("Catcher", 2);
            // double Catcher_Central_Energy_Deposit = Tree_Catcher_Central_Energy_Deposit->At(part_i);
            // double Catcher_Side_Energy_Deposit = Tree_Catcher_Side_Energy_Deposit->At(part_i);
            // H_Catcher_Central_Energy_Deposit[PDG]->Fill(Catcher_Central_Energy_Deposit);
            // H_Catcher_Side_Energy_Deposit[PDG]->Fill(Catcher_Side_Energy_Deposit);
            

            //// PLASTIC SCINTILLATOR ////
            // if (Verbosee == 1) Info("Plastic Scintillator", 2);
            if (Tree_PlasticScintillator_Energy_Deposit->At(part_i) != 0)
            {
                // H_PlasticScintillator_Energy_Deposit[PDG]->Fill(Tree_PlasticScintillator_Energy_Deposit->At(part_i));
                // H_PlasticScintillator_Hit_Position_x[PDG]->Fill(Tree_PlasticScintillator_Hit_Position->At(part_i).x());
                // H_PlasticScintillator_Hit_Position_y[PDG]->Fill(Tree_PlasticScintillator_Hit_Position->At(part_i).y());
                // H_PlasticScintillator_Hit_Position_z[PDG]->Fill(Tree_PlasticScintillator_Hit_Position->At(part_i).z());
                // H_PlasticScintillator_Hit_Position_xy[PDG]->Fill(Tree_PlasticScintillator_Hit_Position->At(part_i).x(), Tree_PlasticScintillator_Hit_Position->At(part_i).y());
                // H_PlasticScintillator_Hit_Angle[PDG]->Fill(Tree_PlasticScintillator_Hit_Angle->At(part_i));
                // H_PlasticScintillator_Hit_Time[PDG]->Fill(Tree_PlasticScintillator_Hit_Time->At(part_i));
                // time_e_IAS = Tree_PlasticScintillator_Hit_Time->At(part_i);
                // Plastic_Full_energy += Tree_PlasticScintillator_Energy_Deposit->At(part_i);
                // Plastic_Full_energy_vec.push_back(Tree_PlasticScintillator_Energy_Deposit->At(part_i));
                // H_PlasticScintillator_Energy_Deposit[0]->Fill(Tree_PlasticScintillator_Energy_Deposit->At(part_i));

                DetectorSignals_RAW[100].push_back(Signal(100, Tree_PlasticScintillator_Hit_Time->At(part_i), Tree_PlasticScintillator_Energy_Deposit->At(part_i), 0, 0.));
            }
        
            //// LOOP ON HITTED SILICON DETECTORS ////
            // if (Verbosee == 1) Info("Silicon Detector", 2);
            for (int sili_i = 0; sili_i < Tree_Silicon_Detector_Energy_Deposit->At(part_i).size(); sili_i++)
            {
                int Silicon_Detector_Code = Tree_Silicon_Detector_Code->At(part_i).at(sili_i);
                double TrigTime = Tree_Silicon_Detector_Hit_Time->At(part_i).at(sili_i);
                double TotalEnergyDeposit = Tree_Silicon_Detector_Energy_Deposit->At(part_i).at(sili_i);
                // H_Silicon_Detector_Hit_Position_x[PDG]->Fill(X);
                // H_Silicon_Detector_Hit_Position_y[PDG]->Fill(Y);
                // H_Silicon_Detector_Hit_Position_z[PDG]->Fill(Z);
                // H_Silicon_Detector_Hit_Position_xy[PDG]->Fill(X, Y);
                // H_Silicon_Detector_Hit_Position_xz[PDG]->Fill(X, Z);
                // H_Silicon_Detector_Hit_Position_yz[PDG]->Fill(Y, Z);
                // H_Silicon_Detector_Hit_Angle[PDG]->Fill(THETA);
                // H_Silicon_Detector_Hit_Time[PDG]->Fill(Tree_Silicon_Detector_Hit_Time->At(part_i).at(sili_i));
                // H_Silicon_Detector_Hit_Anglez[PDG]->Fill(THETA, Z);
                // H_Silicon_Detector_Hit_Anglexy[PDG]->Fill(THETA, X);
                // H_Silicon_Detector_Code[PDG]->Fill(Silicon_Detector_Code);
            
                // H_Silicon_Detector_Energy_Deposit_Det[PDG][Silicon_Detector_Code]->Fill(TotalEnergyDeposit);

                // Full_energy[Silicon_Detector_Code] += TotalEnergyDeposit;
                // Full_energy_without_interstrip[Silicon_Detector_Code] += TotalEnergyDeposit;
                // H_Silicon_Detector_Hit_Angle_Det[PDG][Silicon_Detector_Code]->Fill(THETA);
                
                DetectorSignals_RAW[Silicon_Detector_Code].push_back(Signal(Silicon_Detector_Code, TrigTime, TotalEnergyDeposit, 0, 0.));
            }
        }

        // Printer for the all detector signals 
        if (Verbosee == 1) 
        {   
            Info("All Detector signals", 2); 
            for (int det = 0; det < SIGNAL_MAX; det++)
            {
                if (DetectorSignals_RAW[det].empty())
                    continue;
                for (Signal s : DetectorSignals_RAW[det])
                {
                   Info("Detector: " + detectorName[s.Label] + " Time: " + to_string(s.Time) + " Energy: " + to_string(s.Channel), 3);
                }
            }
        }

        // Adding energyDeposit per detector considering Time Integration of detectors
        vector<Signal> DetectorSignals = TimeIntegration(DetectorSignals_RAW);
        if (Verbosee == 1) 
        {
            Info("Time Integration", 2);
            for (Signal s : DetectorSignals)
            {
                Info("Detector: " + detectorName[s.Label] + " Time: " + to_string(s.Time) + " Energy: " + to_string(s.Channel), 3);
            }
            if (DetectorSignals.empty()) 
                continue;
        }
        // erasing low energy signals for strips
        int size = DetectorSignals.size();
        for (int int_s = DetectorSignals.size() - 1; int_s >= 0; int_s--)
        {
            Signal s = DetectorSignals[int_s];
            if (IsDetectorSiliStrip(s.Label) && s.Channel < 100.)
            {
                DetectorSignals.erase(DetectorSignals.begin() + int_s);
            }
        }

        // Grouping the signals as faster
        vector<vector<Signal>> DetectorSignals_Group = LossLessToGroup(DetectorSignals);
        if (Verbosee == 1) 
        {
            Info("Grouping signals", 2);
            for (int igroup = 0; igroup < DetectorSignals_Group.size(); igroup++)
            {
                vector<Signal> vecs = DetectorSignals_Group[igroup];
                Info("Group: " + to_string(igroup), 3);
                for (Signal s : vecs)
                {
                    Info("Detector: " + detectorName[s.Label] + " Time: " + to_string(s.Time) + " Energy: " + to_string(s.Channel), 4);
                }
            }
        }

        // Saving //
        // Selcting only 1 strip per subgroup
        for (auto Group : DetectorSignals_Group)
        {

            // mul strip == 1
            int multiplicity = 0;
            for (auto s : Group)
            {
                if (IsDetectorSiliStrip(s.Label))
                {
                    multiplicity++;
                }
            }
            if (multiplicity != 1)
                continue;
            
            // filling final tree
            for (auto s : Group)
            {
                if (IsDetectorSiliStrip(s.Label))
                {
                    FINAL_OutTree_Silicon = s;
                }
                else
                {
                    FINAL_OutTree_SiPM.push_back(s);
                }
            }
            FINAL_OutTree->Fill();
            FINAL_OutTree_Silicon = Signal();
            FINAL_OutTree_SiPM.clear();
        }
        
        if (Verbosee == 1) Info("Silicon detectors filling", 2);       
        // Pre-analysis // 
        for (int i_group = 0; i_group < DetectorSignals_Group.size(); i_group++)
        {
            vector<Signal> vecs = DetectorSignals_Group[i_group];
            for (Signal s : vecs)
            {
                if (IsDetectorSiliStrip(s.Label))
                {
                    if (s.Channel > 0)
                    {                       
                        H_Silicon_Detector_Energy_Deposit_Det[0][s.Label]->Fill(s.Channel);
                        e = s.Channel;
                        OutTree[s.Label]->Fill();
                        H_Silicon_Detector_Energy_Deposit_SINGLE[0][s.Label]->Fill(s.Channel);

                        bool coinc = false;
                        for (Signal ss : vecs)
                        {
                            if (ss.Label == 100)
                            {
                                if (abs(ss.Time - s.Time) < 300)
                                {
                                    H_Silicon_Detector_Energy_Deposit_COINC[0][s.Label]->Fill(s.Channel);
                                    H_DeltaCoincTime->Fill(ss.Time - s.Time);
                                    coinc = true;
                                }
                                break;
                            }
                        }
                        if (!coinc)
                        {
                            H_Silicon_Detector_Energy_Deposit_NOCOINC[0][s.Label]->Fill(s.Channel);
                        }
                    }
                }
                Full_energy[s.Label] = 0;
                Full_energy_without_interstrip[s.Label] = 0;
                e = 0;
            }
        }

        sili_e = 0;
        sili_code = 0;
        SiPM_e = 0;
        Plastic_Full_energy = 0;
        Plastic_Full_energy_vec.clear();

        if (Verbosee == 1) Info("End", 2);
    }

    WriteHistograms();
    ANALYSIS_File->cd();
    // for complet analysis
    FINAL_OutTree->Write();
    // for sipm calib
    PlasticIASTree->Write();
    // for silicon calib
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            OutTree[i]->Write();
        }
    }
    SIMULATED_File->Close();
    ANALYSIS_File->Close();

    return 0;
}