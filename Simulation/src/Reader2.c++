#include "Reader2.hh"
#include <chrono>
#include <thread>
default_random_engine generator;
int main(int argc, char *argv[])
{
    FLAG2024 = true;
    InitDetectors("../Grouper/Config_Files/sample.pid");
    // string name = "../../../../../../mnt/hgfs/shared-2/2024_DATA/SIMULATED_DATA/24-02/207Bi_thin";
    // string name = "/run/media/local1/Disque_Dur/2024_DATA/SIMULATED_DATA/03-17/207Bi_100um_CS0_CSP0_CV1_CVP1"; 
    // string name = "../../../../../../mnt/hgfs/shared-2/2024_DATA/Time_test_2";
    // string name = DIR_DATA_HDD + "../SIMULATED_DATA/05-18/32Arx0.0_y4.0_z3.0_CS0_CSP0_CV1_CVP1";

    // y and z from argv
    // string y = argv[1];
    // string z = argv[2];
    // string theta = argv[3];

    // string e_str = argv[1];

    string path = DIR_DATA_HDD + "../SIMULATED_DATA/06-09/";
    // string path = "/run/media/local1/DATANEX/Samuel-G4/06-03/";
    // string name = "proton_" + e_str + "MeV";
    string name = "33Ar_ENSDF_CS0_CSP0_CV1_CVP1";

    
    // name = "241Am_700nm_width";
    TFile *SIMULATED_File = MyTFile((path+name + ".root").c_str(), "READ");
    if (SIMULATED_File == nullptr)
    {
        Error("Impossible to open " + path + name + ".root");
    }

    // if _analaysed already exists exit(O)
    ANALYSIS_File = MyTFile((path + name + "_analysed.root").c_str(), "READ", "Q");
    if (ANALYSIS_File != nullptr)
    {
        Error("File already analysed");
    }


    ANALYSIS_File = MyTFile((path + name + "_analysed.root").c_str(), "RECREATE");
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
    Tree_Catcher_Central_Energy_Deposit = new TTreeReaderArray<double>(*Reader, "Catcher_Central_Energy_Deposit");
    Tree_Catcher_Side_Energy_Deposit = new TTreeReaderArray<double>(*Reader, "Catcher_Side_Energy_Deposit");
    Tree_PlasticScintillator_Energy_Deposit = new TTreeReaderArray<double>(*Reader, "PlasticScintillator_Energy_Deposit");
    Tree_PlasticScintillator_Hit_Position = new TTreeReaderArray<Hep3Vector>(*Reader, "PlasticScintillator_Hit_Position");
    Tree_PlasticScintillator_Hit_Angle = new TTreeReaderArray<double>(*Reader, "PlasticScintillator_Hit_Angle");
    Tree_PlasticScintillator_Hit_Time = new TTreeReaderArray<double>(*Reader, "PlasticScintillator_Hit_Time");
    Tree_Silicon_Detector_Energy_Deposit = new TTreeReaderArray<vector<double>>(*Reader, "Silicon_Detector_Energy_Deposit");
    Tree_Silicon_Detector_Hit_Position = new TTreeReaderArray<vector<Hep3Vector>>(*Reader, "Silicon_Detector_Hit_Position");
    Tree_Silicon_Detector_Hit_Angle = new TTreeReaderArray<vector<double>>(*Reader, "Silicon_Detector_Hit_Angle");
    Tree_Silicon_Detector_Hit_Time = new TTreeReaderArray<vector<double>>(*Reader, "Silicon_Detector_Hit_Time");
    Tree_Silicon_Detector_Code = new TTreeReaderArray<vector<int>>(*Reader, "Silicon_Detector_Code");
    Tree_Silicon_Detector_DL_Energy_Deposit = new TTreeReaderArray<vector<double>>(*Reader, "Silicon_Detector_DL_Energy_Deposit");
    // Tree_Silicon_Detector_InterStrip_Code = new TTreeReaderArray<vector<int>>(*Reader, "Silicon_Detector_InterStrip_Code");
    // tree->SetBranchAddress("Silicon_Detector_InterStrip_Energy_Deposit", &Tree_Silicon_Detector_InterStrip_Energy_Deposit);
    // tree->SetBranchAddress("Silicon_Detector_InterStrip_Hit_Position", &Tree_Silicon_Detector_InterStrip_Hit_Position);

    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            OutTree[i] = new TTree(("Tree_" + detectorName[i]).c_str(), ("Tree_" + detectorName[i]).c_str());
            OutTree[i]->Branch("Energy", &e);
        }
    }
    PlasticIASTree = new TTree("PlasticIAS", "PlasticIAS");
    PlasticIASTree->Branch("Code", &sili_code);
    PlasticIASTree->Branch("Energy", &sili_e);
    PlasticIASTree->Branch("SiPM", &SiPM_e);

    OutTree_Total = new TTree("Tree", "Tree");
    OutTree_Total->Branch("Silicon_Energy", &sili_e);
    OutTree_Total->Branch("Silicon_Code", &sili_code);
    OutTree_Total->Branch("SiPM_Energy", &SiPM_e);

    int Entries = tree->GetEntries();
    clock_t start = clock(), Current;

    Verbosee = 0;

    Init();
    InitHistograms(0);

    Info("Starting Loop");
    while (Reader->Next())
    {
        ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "Reading Tree");

        //// LOOP ON PARTICLE ////
        for (int part_i = 0; part_i < Tree_PDG->GetSize(); part_i++)
        {
            if (Verbosee == 1) Info("Particle: " + to_string(part_i), 1);
            int PDG = Tree_PDG->At(part_i);
            if (PDG > 1000000000)
            {
                PDG = (PDG / 10) * 10;
            }
            if (H_E0.find(PDG) == H_E0.end())
                InitHistograms(PDG);

            //// INITIAL ////
            if (Verbosee == 1) Info("Initial state", 2);
            double E0 = Tree_E0->At(part_i);
            double x = Tree_x->At(part_i)*1e3;
            double y = Tree_y->At(part_i)*1e3;
            double z = Tree_z->At(part_i)*1e6;
            double px = Tree_px->At(part_i);
            double py = Tree_py->At(part_i);
            double pz = Tree_pz->At(part_i);
            double Catcher_Central_Energy_Deposit = Tree_Catcher_Central_Energy_Deposit->At(part_i);
            double Catcher_Side_Energy_Deposit = Tree_Catcher_Side_Energy_Deposit->At(part_i);

            H_E0[PDG]->Fill(E0);
            H_x[PDG]->Fill(x);
            H_y[PDG]->Fill(y);
            H_z[PDG]->Fill(z);
            H_px[PDG]->Fill(px);
            H_py[PDG]->Fill(py);
            H_pz[PDG]->Fill(pz);

            //// CATCHER ////
            if (Verbosee == 1) Info("Catcher", 2);
            H_Catcher_Central_Energy_Deposit[PDG]->Fill(Catcher_Central_Energy_Deposit);
            H_Catcher_Side_Energy_Deposit[PDG]->Fill(Catcher_Side_Energy_Deposit);

            //// PLASTIC SCINTILLATOR ////
            if (Verbosee == 1) Info("Plastic Scintillator", 2);
            if (Tree_PlasticScintillator_Energy_Deposit->At(part_i) != 0)
            {
                H_PlasticScintillator_Energy_Deposit[PDG]->Fill(Tree_PlasticScintillator_Energy_Deposit->At(part_i));
                H_PlasticScintillator_Hit_Position_x[PDG]->Fill(Tree_PlasticScintillator_Hit_Position->At(part_i).x());
                H_PlasticScintillator_Hit_Position_y[PDG]->Fill(Tree_PlasticScintillator_Hit_Position->At(part_i).y());
                H_PlasticScintillator_Hit_Position_z[PDG]->Fill(Tree_PlasticScintillator_Hit_Position->At(part_i).z());
                H_PlasticScintillator_Hit_Position_xy[PDG]->Fill(Tree_PlasticScintillator_Hit_Position->At(part_i).x(), Tree_PlasticScintillator_Hit_Position->At(part_i).y());
                H_PlasticScintillator_Hit_Angle[PDG]->Fill(Tree_PlasticScintillator_Hit_Angle->At(part_i));
                // H_PlasticScintillator_Hit_Time[PDG]->Fill(Tree_PlasticScintillator_Hit_Time->At(part_i));
                time_e_IAS = Tree_PlasticScintillator_Hit_Time->At(part_i);
                Plastic_Full_energy += Tree_PlasticScintillator_Energy_Deposit->At(part_i);
                Plastic_Full_energy_vec.push_back(Tree_PlasticScintillator_Energy_Deposit->At(part_i));
                H_PlasticScintillator_Energy_Deposit[0]->Fill(Tree_PlasticScintillator_Energy_Deposit->At(part_i));
                // SiPM_e = Tree_PlasticScintillator_Energy_Deposit->At(part_i);
                
            }
        
            //// LOOP ON HITTED SILICON DETECTORS ////
            if (Verbosee == 1) Info("Silicon Detector", 2);
            for (int sili_i = 0; sili_i < Tree_Silicon_Detector_Energy_Deposit->At(part_i).size(); sili_i++)
            {
                double X = Tree_Silicon_Detector_Hit_Position->At(part_i).at(sili_i).x();
                double Y = Tree_Silicon_Detector_Hit_Position->At(part_i).at(sili_i).y();
                double Z = Tree_Silicon_Detector_Hit_Position->At(part_i).at(sili_i).z();    
                double THETA = Tree_Silicon_Detector_Hit_Angle->At(part_i).at(sili_i);
                int Silicon_Detector_Code = Tree_Silicon_Detector_Code->At(part_i).at(sili_i);
                H_Silicon_Detector_Hit_Position_x[PDG]->Fill(X);
                H_Silicon_Detector_Hit_Position_y[PDG]->Fill(Y);
                H_Silicon_Detector_Hit_Position_z[PDG]->Fill(Z);
                H_Silicon_Detector_Hit_Position_xy[PDG]->Fill(X, Y);
                H_Silicon_Detector_Hit_Position_xz[PDG]->Fill(X, Z);
                H_Silicon_Detector_Hit_Position_yz[PDG]->Fill(Y, Z);
                H_Silicon_Detector_Hit_Angle[PDG]->Fill(THETA);
                H_Silicon_Detector_Hit_Time[PDG]->Fill(Tree_Silicon_Detector_Hit_Time->At(part_i).at(sili_i));
                H_Silicon_Detector_Hit_Anglez[PDG]->Fill(THETA, Z);
                H_Silicon_Detector_Hit_Anglexy[PDG]->Fill(THETA, X);
                H_Silicon_Detector_Code[PDG]->Fill(Silicon_Detector_Code);
            
                H_Silicon_Detector_Energy_Deposit_Det[PDG][Silicon_Detector_Code]->Fill(Tree_Silicon_Detector_Energy_Deposit->At(part_i).at(sili_i));

                Full_energy[Silicon_Detector_Code] += Tree_Silicon_Detector_Energy_Deposit->At(part_i).at(sili_i);
                Full_energy_without_interstrip[Silicon_Detector_Code] += Tree_Silicon_Detector_Energy_Deposit->At(part_i).at(sili_i);
                H_Silicon_Detector_Hit_Angle_Det[PDG][Silicon_Detector_Code]->Fill(THETA);
            }

            // for (int sili_inter_i = 0; sili_inter_i < Tree_Silicon_Detector_InterStrip_Code->At(part_i).size(); sili_inter_i++)
            // {
            //     int Silicon_Detector_InterStrip_Code = Tree_Silicon_Detector_InterStrip_Code->At(part_i).at(sili_inter_i);
            //     double energy_interstrip = 0;
            //     int code_strip_for = (Silicon_Detector_InterStrip_Code+50)/100;

            //     tree->GetEntry(Reader->GetCurrentEntry());
                // for (int step = 0; step < Tree_Silicon_Detector_InterStrip_Energy_Deposit->at(part_i).at(sili_inter_i).size(); step++)
                // {
                //     double X = Tree_Silicon_Detector_InterStrip_Hit_Position->at(part_i).at(sili_inter_i).at(step).x();
                //     double Y = Tree_Silicon_Detector_InterStrip_Hit_Position->at(part_i).at(sili_inter_i).at(step).y();
                //     double Z = Tree_Silicon_Detector_InterStrip_Hit_Position->at(part_i).at(sili_inter_i).at(step).z();
                //     H_Silicon_Detector_InterStrip_xy[PDG]->Fill(X, Y);
                //     H_Silicon_Detector_InterStrip_z[PDG]->Fill(Z);
                //     energy_interstrip += Tree_Silicon_Detector_InterStrip_Energy_Deposit->at(part_i).at(sili_inter_i).at(step);

                //     // giving energy to strips
                //     // new 
                //     if (Z + (300 * 1e-3 + 800 * 1e-6) / 2 > 300 * 1e-3)
                //         continue;
                //     double ratio = (Y+size_interstrip)/size_interstrip;
                    
                //     if ( Y > 0)
                //     {
                //         InterStrip[code_strip_for] += ratio * Tree_Silicon_Detector_InterStrip_Energy_Deposit->at(part_i).at(sili_inter_i).at(step);
                //         InterStrip[code_strip_for-1] += (1-ratio) * Tree_Silicon_Detector_InterStrip_Energy_Deposit->at(part_i).at(sili_inter_i).at(step);
                //         Full_energy[code_strip_for] += ratio * Tree_Silicon_Detector_InterStrip_Energy_Deposit->at(part_i).at(sili_inter_i).at(step);
                //         Full_energy[code_strip_for-1] += (1-ratio) * Tree_Silicon_Detector_InterStrip_Energy_Deposit->at(part_i).at(sili_inter_i).at(step);
                //     }
                //     else
                //     {
                //         InterStrip[code_strip_for] += (1-ratio) * Tree_Silicon_Detector_InterStrip_Energy_Deposit->at(part_i).at(sili_inter_i).at(step);
                //         InterStrip[code_strip_for-1] += ratio * Tree_Silicon_Detector_InterStrip_Energy_Deposit->at(part_i).at(sili_inter_i).at(step);
                //         Full_energy[code_strip_for] += (1-ratio) * Tree_Silicon_Detector_InterStrip_Energy_Deposit->at(part_i).at(sili_inter_i).at(step);
                //         Full_energy[code_strip_for-1] += ratio * Tree_Silicon_Detector_InterStrip_Energy_Deposit->at(part_i).at(sili_inter_i).at(step);
                //     }
                   
                //     // cout << "Code: " << Silicon_Detector_InterStrip_Code << " Strip: " << code_strip_for << " Ratio: " << ratio << endl;
                // }

                //energy deposit in the interstrip for each particle
                // if (energy_interstrip != 0)
                    // H_Silicon_Detector_InterStrip_Energy_Deposit_Det[PDG][Silicon_Detector_InterStrip_Code]->Fill(energy_interstrip);
            // }
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

        // calculating "rear" energy
        
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorSiliStrip(i))
            {
                RearEnergy[GetDetector(i)] += Full_energy[i];
            }
        }

        if (Verbosee == 1) Info("Silicon detectors filling", 2);
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorSiliStrip(i))
            {
                if (Full_energy[i] > 0)
                {
                    H_Silicon_Detector_Energy_Deposit_Det[0][i]->Fill(Full_energy[i]); 
                    H_Silicon_Detector_Energy_Deposit_Det_without_interstrip[0][i]->Fill(Full_energy_without_interstrip[i]);   
                    H_Silicon_Detector_Energy_Deposit_Det_Rear[0][i]->Fill(RearEnergy[GetDetector(i)], Full_energy[i]);
                    e = Full_energy[i];
                    if (counter_strip_trigger == 1)
                        OutTree[i]->Fill();
                    

                    //plastic scintillator tree
                    
                    sili_code = i;
                    sili_e = Full_energy[i];

                    for (int j = 0; j < Plastic_Full_energy_vec.size(); j++)
                    {
                        SiPM_e = Plastic_Full_energy_vec[j];
                        PlasticIASTree->Fill();
                    }
                    
                    // double offset = Calibration_Function[i]->GetParameter(1) * Detector_Resolution[i] + Calibration_Function[i]->GetParameter(2) * pow(Detector_Resolution[i], 2);
                    // normal_distribution<double> distribution(Full_energy[i], offset + sqrt(3.6 * 1e-3 * 0.134) * sqrt(Full_energy[i]));
                    // Full_energy[i] = distribution(generator);
                    H_Silicon_Detector_Energy_Deposit_SINGLE[0][i]->Fill(Full_energy[i]);
                    if (Plastic_Full_energy > 100)
                    {
                        H_Silicon_Detector_Energy_Deposit_COINC[0][i]->Fill(Full_energy[i]);
                    }
                    else
                    {
                        H_Silicon_Detector_Energy_Deposit_NOCOINC[0][i]->Fill(Full_energy[i]);
                    }
                }
            }
            Full_energy[i] = 0;
            Full_energy_without_interstrip[i] = 0;
            e = 0;

            if (IsDetectorSiliBack(i))
            {
                RearEnergy[GetDetector(i)] = 0;
            }
        }

        //
        if (Verbosee == 1) Info("Platic trees", 2);
        if (sili_e > 3200 && sili_e < 3400)
        {
            // H_PlasticScintillator_Hit_Time[-11]->Fill(time_e_IAS);
            H->Fill(Plastic_Full_energy);
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
    PlasticIASTree->Write();
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
