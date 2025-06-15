#include "Efficiency.hh"
default_random_engine generator;
#include <regex>
#include <filesystem>

int VERBOSE = 0;



int main()
{
    FLAG2024 = true;
    InitDetectors("../Grouper/Config_Files/sample.pid");
    InitWindows();

    if (YEAR == 2025)
    {
        DIR_DATA_HDD += "../../2024_DATA/DETECTOR_DATA/";
    }

    bool display_ref = true;

    TFile *output = MyTFile(("Efficiency_" + to_string(YEAR) + ".root").c_str(), "RECREATE");


    map<string, vector<vector<int>>> Peaks;

    Peaks["32Ar"] = {{1}, {3}, {5, 6}, {25}, {29}, {30}};
    Peaks["33Ar"] = {{2}, {3, 4, 5, 6, 7}, {12}, {26}};


    ////////////// ######################## EXPERIMENTAL ######################## /////////////
    TFile *EXP_FILE = MyTFile((DIR_ROOT_DATA_CALIBRATED + "Calibrated_" + to_string(YEAR) + ".root").c_str(), "READ");
    if (EXP_FILE == nullptr)
        Error("Impossible to open " + DIR_ROOT_DATA_CALIBRATED + "Calibrated_" + to_string(YEAR) + "test.root");

    map<string, TH1D*[SIGNAL_MAX]> H_Exp;
    map<string, map<int, TGraphErrors *>> G_Efficiency_Exp;
    map<int, TGraphErrors *> G_Counts_Exp_Det;

    for (string Nucleus : {"32Ar", "33Ar"})
    {
        InitPeakData(Nucleus);
        for (int det = 0; det < SIGNAL_MAX; det++)
        {
            if (IsDetectorSiliStrip(det))
            {
                Info("Loading Experimental Efficiency for " + Nucleus + " in " + detectorName[det], 1);
                H_Exp[Nucleus][det] = (TH1D *)EXP_FILE->Get((detectorName[det] + "/" + Nucleus + "/h_sub" + detectorName[det]).c_str());
                if (H_Exp[Nucleus][det] == nullptr)
                {
                    H_Exp[Nucleus][det] = (TH1D *)EXP_FILE->Get((detectorName[det] + "/" + Nucleus + "/H_Exp_" + Nucleus + "_" + detectorName[det]).c_str());
                }

                G_Efficiency_Exp[Nucleus][det] = new TGraphErrors();

                for (int i_peak = 0; i_peak < Peaks[Nucleus].size(); i_peak++)
                {
                    // int peak = Peaks[Nucleus][i_peak];

                    // if (peak == IAS[Nucleus])
                    // {
                    //     continue; // Skip IAS peak
                    // }

                    // double energy = PeakData[Nucleus][peak].at(0);

                    // H_Exp[Nucleus][det]->GetXaxis()->SetRangeUser(WindowsMap[Nucleus][peak][det].first, WindowsMap[Nucleus][peak][det].second);
                    // double integral1 = H_Exp[Nucleus][det]->Integral();
                    // H_Exp[Nucleus][det]->GetXaxis()->SetRangeUser(WindowsMap[Nucleus][IAS[Nucleus]][det].first, WindowsMap[Nucleus][IAS[Nucleus]][det].second);
                    // double integralIAS = H_Exp[Nucleus][det]->Integral();

                    // double efficiency = integral1 / integralIAS / (PeakData[Nucleus][peak].at(2) / PeakData[Nucleus][IAS[Nucleus]].at(2));    

                    // double efficiency_err_stat = sqrt(pow(sqrt(integral1)/integralIAS, 2) + pow((sqrt(integralIAS)*integral1)/pow(integralIAS, 2), 2)) / (PeakData[Nucleus][peak].at(2) / PeakData[Nucleus][IAS[Nucleus]].at(2));
                    // // double efficiency_err_syst = sqrt(pow(efficiency, 2) * (pow(PeakData[Nucleus][peak][1], 2) + pow(PeakData[Nucleus][IAS[Nucleus]][1], 2))) / (PeakData[Nucleus][peak].at(2) / PeakData[Nucleus][IAS[Nucleus]].at(2));


                    // G_Efficiency_Exp[Nucleus][det]->AddPoint(energy, efficiency);
                    // G_Efficiency_Exp[Nucleus][det]->SetPointError(G_Efficiency_Exp[Nucleus][det]->GetN() - 1, PeakData[Nucleus][peak][1], efficiency_err_stat);


                    double energy = 0; 
                    double Br = 0;
                    for (int subpeak : Peaks[Nucleus][i_peak])
                    {
                        energy += PeakData[Nucleus][subpeak].at(0) * PeakData[Nucleus][subpeak].at(2);
                        Br += PeakData[Nucleus][subpeak].at(2);
                    }

                    energy /= Br;

                    H_Exp[Nucleus][det]->GetXaxis()->SetRangeUser(WindowsMap[Nucleus][Peaks[Nucleus][i_peak][0]][det].first, WindowsMap[Nucleus][Peaks[Nucleus][i_peak].back()][det].second);
                    double integral1 = H_Exp[Nucleus][det]->Integral();
                    H_Exp[Nucleus][det]->GetXaxis()->SetRangeUser(WindowsMap[Nucleus][IAS[Nucleus]][det].first, WindowsMap[Nucleus][IAS[Nucleus]][det].second);
                    double integralIAS = H_Exp[Nucleus][det]->Integral();

                    double efficiency = integral1 / integralIAS / (Br / PeakData[Nucleus][IAS[Nucleus]].at(2));

                    double efficiency_err_stat = sqrt(pow(sqrt(integral1)/integralIAS, 2) + pow((sqrt(integralIAS)*integral1)/pow(integralIAS, 2), 2)) / (Br / PeakData[Nucleus][IAS[Nucleus]].at(2));


                    G_Efficiency_Exp[Nucleus][det]->AddPoint(energy, efficiency);
                    G_Efficiency_Exp[Nucleus][det]->SetPointError(G_Efficiency_Exp[Nucleus][det]->GetN() - 1, PeakData[Nucleus][Peaks[Nucleus][i_peak][0]][1], efficiency_err_stat);            
                
                }
            }
        }
    }

    output->cd();
    for (int detector= 1; detector <= SILI_NUM; detector++)
    {
        Info("Drawing Efficiency for Detector " + to_string(detector), 1);
        TCanvas *canvas = new TCanvas(("Efficiency_D" + to_string(detector)).c_str(), ("Efficiency_D" + to_string(detector)).c_str(), 800, 600);
        TMultiGraph *mg = new TMultiGraph();
        TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);

        for (string Nucleus: {"32Ar", "33Ar"})
        {
            if (G_Efficiency_Exp[Nucleus][detector * 10 + 1] == nullptr)
            {
                Warning("No efficiency data for " + Nucleus + " in detector " + to_string(detector));
                continue;
            }

            for (int strip = 1; strip <= 5; strip++)
            {
                if (G_Efficiency_Exp[Nucleus][detector * 10 + strip] == nullptr)
                {
                    Warning("No efficiency data for " + Nucleus + " in detector " + to_string(detector) + " strip " + to_string(strip));
                    continue;
                }

                G_Efficiency_Exp[Nucleus][detector * 10 + strip]->SetMarkerStyle(20);
                G_Efficiency_Exp[Nucleus][detector * 10 + strip]->SetMarkerSize(0.5*strip);
                G_Efficiency_Exp[Nucleus][detector * 10 + strip]->SetLineColor(NucleusColor[Nucleus]);
                G_Efficiency_Exp[Nucleus][detector * 10 + strip]->SetMarkerColor(NucleusColor[Nucleus]);
                G_Efficiency_Exp[Nucleus][detector * 10 + strip]->SetLineWidth(1);

                legend->AddEntry(G_Efficiency_Exp[Nucleus][detector * 10 + strip], (Nucleus + "  Strip " + to_string(strip)).c_str(), "P");

            mg->Add(G_Efficiency_Exp[Nucleus][detector * 10 + strip], "PL");
            }

            
        }

        mg->GetXaxis()->SetTitle("Energy [keV]");
        mg->GetYaxis()->SetTitle("Efficiency compared to IAS");
        mg->GetXaxis()->CenterTitle();
        mg->GetYaxis()->CenterTitle();
        mg->SetTitle(("Efficiency of Detector " + to_string(detector)).c_str());
        mg->Draw();

        legend->Draw("SAME");

        canvas->Write();
    }

//     output->Close();

//     EXP_FILE->Close();

//     exit(0);


//   /*
    
    ////////////// ######################## SIMULATION ######################## /////////////
    vector<string> FileName;
    string path = DIR_DATA_HDD+"../SIMULATED_DATA/05-24/";
    DIR* dir = opendir(path.c_str());
    if (dir == nullptr) {
        Error("Unable to open directory: " + path);
        return 1;
    }


    //analysed
    // struct dirent* entry;
    // while ((entry = readdir(dir)) != nullptr) {
    //     string file_name = entry->d_name;

    //     if (file_name == "." || file_name == "..") continue;

    //     if (file_name.find("_analysed.root") != string::npos) {
    //         TFile *file = new TFile((path + file_name).c_str(), "READ");
    //         if (file != nullptr)
    //         {
    //             FileName.push_back(path+file_name);
    //             file->Close();
    //         }
    //     }
    // }
    // closedir(dir);
    

    //other
    // path = "/run/media/local1/DATANEX/Samuel-G4/new/new/";
    // dir = opendir(path.c_str());
    // if (dir == nullptr) {
    //     Error("Unable to open directory: " + path);
    //     return 1;
    // }
    // while ((entry = readdir(dir)) != nullptr) {
    //     string file_name = entry->d_name;

    //     if (file_name == "." || file_name == "..") continue;

    //     if (file_name.find("_CVP1.root") != string::npos) {
    //         TFile *file = new TFile((path + file_name).c_str(), "READ");
    //         if (file != nullptr)
    //         {
    //             FileName.push_back(path+file_name);
    //             file->Close();
    //         }
    //     }
    // }
    // closedir(dir);

    FileName.push_back("/run/media/local1/Disque_Dur/2024_DATA/SIMULATED_DATA/06-01/32Ar_ENSDF_CS0_CSP0_CV1_CVP1_analysed.root");
    FileName.push_back("/run/media/local1/Disque_Dur/2024_DATA/SIMULATED_DATA/06-01/33Ar_ENSDF_CS0_CSP0_CV1_CVP1_analysed.root");
    map<int, map<string, map<int, TGraphErrors*>>> G_Efficiency_Sim;
    // map position nucleus det
    map<int, Position> Position_Map;
    int counter = 0;    
    for (const auto &file_name : FileName)
    {
        // extracting values from filename
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        double theta = 0.0;

        regex pattern(R"x(x(-?\d+\.?\d*)_y(-?\d+\.?\d*)_z(-?\d+\.?\d*)(?:_theta(-?\d+\.?\d*))?_CS)x");
        smatch matches;

        Position pos;

        if (regex_search(file_name, matches, pattern))
        {
            pos.x = stod(matches[1]);
            pos.y = stod(matches[2]);
            pos.z = stod(matches[3]);
            if (matches[4].matched)
                pos.theta = stod(matches[4]);
            else
                pos.theta = 0.0;
        }
        else
        {
            Warning("Filename does not match expected pattern: " + file_name);
            pos.x = 0.0;
            pos.y = 0.0;
            pos.z = 0.0;
            pos.theta = 0.0;
            // continue;
        }

        if (!Accepting_Position(pos))
            continue;

        TFile *SIM_FILE = MyTFile((file_name).c_str(), "READ");
        if (SIM_FILE == nullptr)
            continue;

        string Nucleus;
        if (file_name.find("32Ar") != string::npos)
            Nucleus = "32Ar";
        else if (file_name.find("33Ar") != string::npos)
            Nucleus = "33Ar";
        else
        {
            Error("Nucleus not recognized in file: " + file_name);
            continue;
        }
        
        //Verify if exist 
        if (Position_Map.find(counter) == Position_Map.end())
        {
            // position don't exist, add it
            counter++;
            Position_Map[counter] = pos;
        }
       
        TH1D *H_Sim[2][SILI_SIZE + 1] = {nullptr};

        // adding histogram from file in hemisphere cynlindric symetry
        for (int det = 0; det < SIGNAL_MAX; det++)
        {
            if (IsDetectorSiliStrip(det))
            {
                if (VERBOSE == 1) Info("Strip " + to_string(det), 1);

                int dir = GetDetector(det) <= 4 ? 0 : 1;

                if (GetDetector(det) == 1 || GetDetector(det) == 5)
                {
                    TH1D* h = (TH1D *)SIM_FILE->Get(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_All").c_str());
                    if (h != nullptr)
                        H_Sim[dir][GetDetectorChannel(det)] = (TH1D *)SIM_FILE->Get(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_All").c_str());
                    else
                        H_Sim[dir][GetDetectorChannel(det)] = (TH1D *)SIM_FILE->Get((detectorName[det] + "_single").c_str());
                }
                else
                {
                    TH1D* h = (TH1D *)SIM_FILE->Get(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_All").c_str());
                    if (h != nullptr)
                        H_Sim[dir][GetDetectorChannel(det)]->Add((TH1D *)SIM_FILE->Get(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_All").c_str()));
                    else
                        H_Sim[dir][GetDetectorChannel(det)]->Add((TH1D *)SIM_FILE->Get((detectorName[det] + "_single").c_str()));
                }
            }
        }

        Info("File loaded", 1);
        for (int detector = 1; detector <= SILI_NUM; detector++)
        {
            if (VERBOSE == 1)
                Info("D " + to_string(detector), 1);

            int dir = detector <= 4 ? 0 : 1;

            for (int strip = 1; strip < SILI_SIZE; strip++)
            {
                 if (VERBOSE == 1)
                        Info("Strip " + to_string(strip), 3);

                int det = detector*10 + strip;

                G_Efficiency_Sim[counter][Nucleus][det] = new TGraphErrors();

                for (int i_peak = 0; i_peak < Peaks[Nucleus].size(); i_peak++)
                {
                    vector<int> peaks = Peaks[Nucleus][i_peak];
                    double energy = 0;
                    double br = 0;
                    double br_err = 0;

                    for (int peak : peaks)
                    {
                        energy += PeakData[Nucleus][peak][0] * PeakData[Nucleus][peak][2];
                        br += PeakData[Nucleus][peak][2];
                        br_err += pow(PeakData[Nucleus][peak][3], 2);
                    }

                    energy /= br;


                    // int peak = Peaks[Nucleus][i_peak];
                    // double energy = PeakData[Nucleus][peak][0];

                    if (VERBOSE == 1)
                        Info("Peak " + to_string(peaks[0]) + " Energy: " + to_string(energy), 3);   

                    H_Sim[dir][strip]->GetXaxis()->SetRangeUser(WindowsMap[Nucleus][peaks[0]][det].first, WindowsMap[Nucleus][peaks.back()][det].second);
                    double integral1 = H_Sim[dir][strip]->Integral();

                    H_Sim[dir][strip]->GetXaxis()->SetRangeUser(WindowsMap[Nucleus][IAS[Nucleus]][det].first, WindowsMap[Nucleus][IAS[Nucleus]][det].second);
                    double integralIASexp = H_Sim[dir][strip]->Integral();

                    double efficiency = integral1 / integralIASexp / (br / PeakData[Nucleus][IAS[Nucleus]][2]);

                    double efficiency_err_stat = sqrt(pow(sqrt(integral1) / integralIASexp, 2) + pow(integral1 * sqrt(integralIASexp) / (integralIASexp * integralIASexp), 2)) / (br / PeakData[Nucleus][IAS[Nucleus]][2]);
                    // double efficiency_err_syst = 

                    // double efficiency_err = sqrt(pow(efficiency_err_stat, 2) + pow(efficiency_err_syst, 2));

                    G_Efficiency_Sim[counter][Nucleus][det]->AddPoint(energy, efficiency);
                    G_Efficiency_Sim[counter][Nucleus][det]->SetPointError(G_Efficiency_Sim[counter][Nucleus][det]->GetN() - 1, 0, efficiency_err_stat);
                }
            }
        }

        SIM_FILE->Close();
    }



    ////////////// ######################## DISPLAY ######################## //////////////
    Info("Displaying Efficiency", 1);
    output->cd();
    VERBOSE =1;
    TH3D *H_yztheta[SILI_NUM + 1] = {nullptr};

    map<int, map<int, vector<int>>> Best_Position_Det;

    TMultiGraph *MG_Efficiency_Exp[SILI_NUM + 1];
    map<int, TMultiGraph *[SILI_NUM + 1]> MG_Efficiency_Sim;

    map<double, map<int, int>> Chi2_Value;


    // Building the MULTIGRAPH for each detector and computing chi2
    for (int det = 1; det <= SILI_NUM; det++)
    {
        if (VERBOSE == 1)
            Info("Detector " + to_string(det), 1);

        int dir = det <= 4 ? 0 : 1;

        MG_Efficiency_Exp[det] = new TMultiGraph();

        for (auto pair : Position_Map)
        {

            int counter_pos = pair.first;
            Position pos = pair.second;

            MG_Efficiency_Sim[counter_pos][det] = new TMultiGraph();

            if (VERBOSE == 1)
                Info("Position " + to_string(counter_pos), 2);

            for (int s = 1; s < SILI_SIZE; s++)
            {
                if (VERBOSE == 1)
                    Info("Strip " + to_string(s), 3);

                for (string Nucleus : {"32Ar", "33Ar"})
                {
                    if (VERBOSE == 1)
                        Info("Nucleus " + Nucleus, 4);

                    /// EXPERIMENTAL
                    G_Efficiency_Exp[Nucleus][det * 10 + s]->SetMarkerColor(NucleusColor[Nucleus]);
                    MG_Efficiency_Exp[det]->Add(G_Efficiency_Exp[Nucleus][det * 10 + s], "P");

                    G_Efficiency_Sim[counter_pos][Nucleus][det * 10 + s]->SetMarkerColor(NucleusColor[Nucleus]);
                    G_Efficiency_Sim[counter_pos][Nucleus][det * 10 + s]->SetLineColor(NucleusColor[Nucleus]);
                    

                    
                    MG_Efficiency_Sim[counter_pos][det]->Add(G_Efficiency_Sim[counter_pos][Nucleus][det * 10 + s], "L");
                }
            }

            Chi2_Value[det][Chi2(MG_Efficiency_Exp[det], MG_Efficiency_Sim[counter_pos][det])] = counter_pos;
        }
    }



    // Displaying the efficiency detector
    int Maximum_Position = 10;
    counter = 0;
    for (int det = 1; det <= SILI_NUM; det++)
    {
        TMultiGraph *MG_Efficiency = new TMultiGraph();
        if (VERBOSE == 1)
            Info("Detector " + to_string(det), 1);

        TCanvas *c = new TCanvas(("Efficiency_Detector_" + to_string(det)).c_str(), ("Efficiency_Detector_" + to_string(det)).c_str(), 1920, 1080);
        TLegend *leg = new TLegend(0.7, 0.0, 1.0, 1.0);
        c->cd();
        // MG_Efficiency_Exp[det]->Draw();
        MG_Efficiency->Add(MG_Efficiency_Exp[det], "P");
        MG_Efficiency_Exp[det]->GetXaxis()->SetTitle("Energy (keV)");
        MG_Efficiency_Exp[det]->GetYaxis()->SetTitle("Efficiency");
        MG_Efficiency_Exp[det]->SetTitle(("Efficiency of Detector " + to_string(det)).c_str());

        //loop on Chi2_Value[det]
        for (auto &pair : Chi2_Value[det])
        {
            int counter_pos = pair.second;
            double chi2_value = pair.first;

            if (VERBOSE == 1)
                Info("Chi2 value " + to_string(chi2_value) + " for position " + to_string(counter_pos), 4);

            Position pos = Position_Map[counter_pos];

            ostringstream ossx;
            ossx << fixed << setprecision(1) << pos.x;
            string str_value_x = ossx.str();
            ostringstream ossy;
            ossy << fixed << setprecision(1) << pos.y;
            string str_value_y = ossy.str();
            ostringstream ossz;
            ossz << fixed << setprecision(1) << pos.z;
            string str_value_z = ossz.str();
            ostringstream osstheta;
            osstheta << fixed << setprecision(1) << pos.theta;
            string str_value_theta = osstheta.str();
            ostringstream oss_chi2;
            oss_chi2 << fixed << setprecision(2) << chi2_value;
            string str_value_chi2 = oss_chi2.str();


            // MG_Efficiency_Sim[counter_pos][det]->Draw("SAME");
            MG_Efficiency->Add(MG_Efficiency_Sim[counter_pos][det], "L");
            MG_Efficiency_Sim[counter_pos][det]->GetXaxis()->SetTitle("Energy (keV)");
            MG_Efficiency_Sim[counter_pos][det]->GetYaxis()->SetTitle("Efficiency");
            MG_Efficiency_Sim[counter_pos][det]->SetTitle(("Efficiency of Detector " + to_string(det) + " at position " + to_string(counter_pos)).c_str());

            leg->AddEntry(MG_Efficiency_Sim[counter_pos][det], ("x = " + str_value_x + " y = " + str_value_y + " z = " + str_value_z + " #theta = " + str_value_theta + " #Chi^{2} = " + str_value_chi2).c_str(), "L");
        
        
            if (counter > Maximum_Position)
                break;
            counter++;
        }
        MG_Efficiency->Draw("SAME");
        c->Write();
    }

    output->Close();

}



