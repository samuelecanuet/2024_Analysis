#include "Efficiency2.hh"
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

    strip_ref = 1;
    detector_ref = 4;
    bool display_ref = true;

    string NUCLEUS = "32Ar";

    vector<int> Peaks = {1, 3, 14};

    TTree *Result_Tree[SILI_NUM + 1];
    double y_tree, z_tree, theta_tree, chi2_tree;
    for (int det = 1; det <= SILI_NUM; det++)
    {
        Result_Tree[det] = new TTree(("Efficiency_Ratio_" + to_string(det)).c_str(), ("Efficiency Ratio for " + detectorName[det]).c_str());
        Result_Tree[det]->Branch("y", &y_tree, "y/D");
        Result_Tree[det]->Branch("z", &z_tree, "z/D");
        Result_Tree[det]->Branch("theta", &theta_tree, "theta/D");
        Result_Tree[det]->Branch("chi2", &chi2_tree, "chi2/D");  
    }

    TNtuple *Result_Ntuple[SILI_NUM + 1];
    for (int det = 1; det <= SILI_NUM; det++)
    {
        Result_Ntuple[det] = new TNtuple(("Efficiency_Ratio_Ntuple_" + to_string(det)).c_str(), ("Efficiency Ratio Ntuple for " + detectorName[det]).c_str(), "y:z:theta:chi2");
    }

    
    ////////////// ######################## SIMULATION ######################## /////////////
    vector<string> FileName;
    string path = DIR_DATA_HDD+"../SIMULATED_DATA/05-24/";
    DIR* dir = opendir(path.c_str());
    if (dir == nullptr) {
        Error("Unable to open directory: " + path);
        return 1;
    }


    //analysed
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        string file_name = entry->d_name;

        if (file_name == "." || file_name == "..") continue;

        if (file_name.find("_analysed.root") != string::npos) {
            TFile *file = new TFile((path + file_name).c_str(), "READ");
            if (file != nullptr)
            {
                FileName.push_back(path+file_name);
                file->Close();
            }
        }
    }
    closedir(dir);
    

    //other
    path = "/run/media/local1/DATANEX/Samuel-G4/new/new/";
    dir = opendir(path.c_str());
    if (dir == nullptr) {
        Error("Unable to open directory: " + path);
        return 1;
    }
    while ((entry = readdir(dir)) != nullptr) {
        string file_name = entry->d_name;

        if (file_name == "." || file_name == "..") continue;

        if (file_name.find("_CVP1.root") != string::npos) {
            TFile *file = new TFile((path + file_name).c_str(), "READ");
            if (file != nullptr)
            {
                FileName.push_back(path+file_name);
                file->Close();
            }
        }
    }
    closedir(dir);


    map<int, map<int, map<int, TGraphErrors*>>> G_EfficiencyRatio_Sim;
    map<int, map<int, map<int, double>>> G_Counts_Sim_Det;
    map<int, Position> Position_Map;
    int counter = 0;    
    for (const auto &file_name : FileName)
    {
        counter++;
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
            continue;
        }

        if (!Accepting_Position(pos))
            continue;

        TFile *SIM_FILE = MyTFile((file_name).c_str(), "READ");
        if (SIM_FILE == nullptr)
            continue;
        
       
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

        

        if (VERBOSE == 1) Info("File loaded", 1);
        for (int dir = 0; dir <= 1; dir++)
        {
            if (VERBOSE == 1) Info("Direction " + to_string(dir), 1);

            for (int i_peak = 0; i_peak < Peaks.size(); i_peak++)
            {
                int peak = Peaks[i_peak];
                if (VERBOSE == 1) Info("Peak " + to_string(peak), 3);
                
                G_EfficiencyRatio_Sim[counter][dir][peak] = new TGraphErrors();
                G_Counts_Sim_Det[counter][dir][peak] = 0.;
                Position_Map[counter] = pos;
                for (int strip = 1; strip < SILI_SIZE; strip++)
                {
                    if (VERBOSE == 1) Info("Strip " + to_string(strip), 3);

                    int d = dir == 0 ? 50 : 10; // for windows

                    H_Sim[dir][strip]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][peak][d + strip].first, WindowsMap[NUCLEUS][peak][d + strip].second);
                    double integral1 = H_Sim[dir][strip]->Integral();
                    H_Sim[dir][strip_ref]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][peak][d + strip_ref].first, WindowsMap[NUCLEUS][peak][d + strip_ref].second);
                    double integral2 = H_Sim[dir][strip_ref]->Integral();

                    G_EfficiencyRatio_Sim[counter][dir][peak]->AddPoint(strip, integral1 / integral2);
                    G_EfficiencyRatio_Sim[counter][dir][peak]->SetPointError(G_EfficiencyRatio_Sim[counter][dir][peak]->GetN() - 1, 0, sqrt(pow(sqrt(integral1) / integral2, 2) + pow(integral1 * sqrt(integral2) / (integral2 * integral2), 2)));
                
                    G_Counts_Sim_Det[counter][dir][peak] += integral1;
                }
            }
        }

        SIM_FILE->Close();
    }



    ////////////// ######################## EXPERIMENTAL ######################## /////////////
    TFile *EXP_FILE = MyTFile((DIR_ROOT_DATA_CALIBRATED + "Calibrated_" + to_string(YEAR) + ".root").c_str(), "READ");
    if (EXP_FILE == nullptr)
        Error("Impossible to open " + DIR_ROOT_DATA_CALIBRATED + "Calibrated_" + to_string(YEAR) + ".root");

    TH1D* H_Exp[SIGNAL_MAX] = {nullptr};
    map<int, map<int, TGraphErrors*>> G_EfficiencyRatio_Exp;
    map<int, TGraphErrors*> G_Counts_Exp_Det;


    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliStrip(det))
        {
            H_Exp[det] = (TH1D *)EXP_FILE->Get((detectorName[det] + "/32Ar/h_sub" + detectorName[det]).c_str());
        }
    }

    for (int i_peak = 0; i_peak < Peaks.size(); i_peak++)
    {
        int peak = Peaks[i_peak];
        double integral_detector = 0;

        G_Counts_Exp_Det[peak] = new TGraphErrors();

        for (int det = 0; det < SIGNAL_MAX; det++)
        {
            if (IsDetectorSiliStrip(det))
            {
                H_Exp[det]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][peak][det].first, WindowsMap[NUCLEUS][peak][det].second);
                double integral1 = H_Exp[det]->Integral();
                H_Exp[GetDetector(det) * 10 + strip_ref]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][peak][GetDetector(det) * 10 + strip_ref].first, WindowsMap[NUCLEUS][peak][GetDetector(det) * 10 + strip_ref].second);
                double integral2 = H_Exp[GetDetector(det) * 10 + strip_ref]->Integral();

                if (GetDetectorChannel(det) == 1)
                    G_EfficiencyRatio_Exp[GetDetector(det)][peak] = new TGraphErrors();
                G_EfficiencyRatio_Exp[GetDetector(det)][peak]->AddPoint(GetDetectorChannel(det), integral1 / integral2);
                G_EfficiencyRatio_Exp[GetDetector(det)][peak]->SetPointError(G_EfficiencyRatio_Exp[GetDetector(det)][peak]->GetN() - 1, 0, sqrt(integral1) / integral2 + integral1 * sqrt(integral2) / (integral2 * integral2));
                
                integral_detector += integral1;
                if (GetDetectorChannel(det) == 5)
                {
                    G_Counts_Exp_Det[peak]->AddPoint(GetDetector(det), integral_detector);
                    G_Counts_Exp_Det[peak]->SetPointError(G_Counts_Exp_Det[peak]->GetN() - 1, 0, sqrt(integral_detector));
                    integral_detector = 0; // reset for next detector
                }
            }
        }
    }

    ////////////// ######################## DISPLAY ######################## /////////////
    TFile *output = MyTFile(("EfficiencyRatio_" + to_string(YEAR) + ".root").c_str(), "RECREATE");

    TH2D *Chi2_2D_xy = new TH2D("Chi2_2D_xy", "Chi2_2D_xy", 10, -5.0, 5.0, 10, -5.0, 5.0);
    TH2D *Chi2_2D_yz = new TH2D("Chi2_2D_xz", "Chi2_2D_xz", 10, -5.0, 5.0, 10, -5.0, 5.0);
    TH3D *H_yztheta[SILI_NUM + 1] = {nullptr};

    map<int, map<int, vector<int>>> Best_Position_Det;

    
    for (int i_peak = 0; i_peak < Peaks.size(); i_peak++)
    {
        if (VERBOSE == 1)
            Info("Peak " + to_string(Peaks[i_peak]), 1);
        int peak = Peaks[i_peak];

        TGraphErrors *G_EfficiencyRatio_Exp_Det = new TGraphErrors();
        TGraphErrors *G_EfficiencyRatio_Sim_Det = new TGraphErrors();
        TGraphErrors *G_EfficiencyRatio_SimExp_Det = new TGraphErrors();

        for (int detector = 1; detector <= 8; detector++)
        {
            if (VERBOSE == 1) Info("Detector " + to_string(detector), 2); 

            int dir = detector <= 4 ? 0 : 1;
            TCanvas *c = new TCanvas(("EfficiencyRatio_D"+to_string(detector) + "_" + to_string(peak)).c_str(), ("EfficiencyRatio_D"+to_string(detector) + "_" + to_string(peak)).c_str(), 1920, 1080);
            TPad *p1 = new TPad("p1", "p1", 0.0, 0.0, 0.7, 1.0);
            p1->Draw();
            p1->cd();
            TMultiGraph *mg_det = new TMultiGraph();
            TLegend *leg_det = new TLegend(0.7, 0.0, 1.0, 1.0);
            int counter = 0;

            // Simulation
            if (VERBOSE == 1) Info("Simulation", 3);
            map<double , int> chi2_map;
            for (auto pair_ : Position_Map)
            {
                int pos_i = pair_.first;
                Position pos = pair_.second;
                
                TGraphErrors *ge = G_EfficiencyRatio_Sim[pos_i][dir][peak];
                if (ge == nullptr)
                    continue;

                
                
                double chi2 = Chi2(ge, G_EfficiencyRatio_Exp[detector][peak]);
                chi2_map[chi2] = pos_i;
                if (peak == 14)
                {
                    Chi2_2D_xy->Fill(pos.x, pos.y, chi2);
                    Chi2_2D_yz->Fill(pos.y, pos.z, chi2);

                    if (H_yztheta[detector] == nullptr)
                    {
                        H_yztheta[detector] = new TH3D(("H_yztheta_" + to_string(detector)).c_str(), ("H_yztheta_" + to_string(detector)).c_str(),
                                                       21, -5.0, 5.0,
                                                       21, -5.0, 5.0,
                                                       21, -5.0, 5.0);
                    }
                    H_yztheta[detector]->Fill(pos.y, pos.z, pos.theta, chi2);

                    y_tree = pos.y;
                    z_tree = pos.z;
                    theta_tree = pos.theta;
                    chi2_tree = chi2;
                    Result_Tree[detector]->Fill();
                }

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
                oss_chi2 << fixed << setprecision(2) << chi2;
                string str_value_chi2 = oss_chi2.str();
                // leg_det->AddEntry(G_EfficiencyRatio_Sim[pos_i][dir][peak], ("x = " + str_value_x + " y = " + str_value_y + " z = " + str_value_z + " #theta = " + str_value_theta + " #Chi^{2} = " + str_value_chi2).c_str(), "p");
                
                // cout << "Chi2 : " << Chi2(ge, G_EfficiencyRatio_Exp[detector][peak]) << endl;
                counter++;
            }

            //Experimental
            if (VERBOSE == 1) Info("Experimental", 3);
            G_EfficiencyRatio_Exp[detector][peak]->SetMarkerColor(kRed);
            G_EfficiencyRatio_Exp[detector][peak]->SetMarkerStyle(20);
            G_EfficiencyRatio_Exp[detector][peak]->SetMarkerSize(2);
            G_EfficiencyRatio_Exp[detector][peak]->SetLineColor(kRed);
            G_EfficiencyRatio_Exp[detector][peak]->SetLineWidth(3);
            
            



            // #######################  DISPAYING BEST FITTING POSITIONS IN ORDER ####################### //
            //Legend
            int counter_graph = 0;
            int chi_best = chi2_map.begin()->first;
            for (auto it = chi2_map.begin(); it != chi2_map.end(); ++it)
            {
                int pos_i = it->second;
                double chi2 = it->first;

                if (abs(chi_best - chi2) <= 1.)
                    Best_Position_Det[detector][peak].push_back(pos_i); 
                
                if (counter_graph > 10)
                    continue;

                cout << "Detector " << detector << "   " << "χ2 = " << chi2 << endl;

                Position pos = Position_Map[pos_i];
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
                oss_chi2 << fixed << setprecision(2) << chi2;
                string str_value_chi2 = oss_chi2.str();
                leg_det->AddEntry(G_EfficiencyRatio_Sim[pos_i][dir][peak], ("x = " + str_value_x + " y = " + str_value_y + " z = " + str_value_z + " #theta = " + str_value_theta + " #Chi^{2} = " + str_value_chi2).c_str(), "p");
                
                if (counter_graph == 0)
                {
                    G_EfficiencyRatio_Sim[pos_i][dir][peak]->SetMarkerSize(2);
                    G_EfficiencyRatio_Sim[pos_i][dir][peak]->SetLineWidth(2);
                    G_EfficiencyRatio_Sim[pos_i][dir][peak]->SetMarkerStyle(21);

                    if (peak == 14)
                    {
                        cout << "Detector " << detector << "   " << "χ2 = " << str_value_chi2 << endl;
                        cout << "Best position : " << "x = " << str_value_x << "    y = " << str_value_y << "   z = " << str_value_z << " 	θ = " << str_value_theta << endl;                
                    } 
                }
                G_EfficiencyRatio_Sim[pos_i][dir][peak]->SetMarkerColor(counter_graph + 1);
                G_EfficiencyRatio_Sim[pos_i][dir][peak]->SetLineColor(counter_graph + 1);
                G_EfficiencyRatio_Sim[pos_i][dir][peak]->SetMarkerStyle(20);
                G_EfficiencyRatio_Sim[pos_i][dir][peak]->SetDrawOption("ALP");
                mg_det->Add(G_EfficiencyRatio_Sim[pos_i][dir][peak]);

                counter_graph++;
            }

            mg_det->Add(G_EfficiencyRatio_Exp[detector][peak]);

            p1->cd();
            mg_det->Draw("ALP");

            c->cd();
            leg_det->Draw("SAME");
            c->Write();
        }

        ////////// DETECTOR RATIO ///////////

        for (int detector = 1; detector <= SILI_NUM; detector++)
        {
            int best_pos_i = Best_Position_Det[detector][peak][0];
            int dir = detector <= 4 ? 0 : 1;

            int best_pos_i_ref = Best_Position_Det[detector_ref][peak][0];
            int dir_ref = detector_ref <= 4 ? 0 : 1;

            // cout << "Best position for detector " << detector << " : " << best_pos_i << endl;
            // cout << "Counts Sim : " << G_Counts_Sim_Det[best_pos_i][dir][peak] << endl;
            // cout << "Counts Exp: " << G_Counts_Exp_Det[peak]->GetY()[detector - 1] << endl;
        

            if (VERBOSE == 1)
                Info("Detector " + to_string(detector), 3);


            // SIMULATION
            double ratio = G_Counts_Sim_Det[best_pos_i][dir][peak] / G_Counts_Sim_Det[best_pos_i_ref][dir_ref][peak];
            double ratio_err = sqrt(pow(sqrt(G_Counts_Sim_Det[best_pos_i][dir][peak]) / G_Counts_Sim_Det[best_pos_i_ref][dir_ref][peak], 2) + pow(G_Counts_Sim_Det[best_pos_i][dir][peak] * sqrt(G_Counts_Sim_Det[best_pos_i_ref][dir_ref][peak]) / (G_Counts_Sim_Det[best_pos_i_ref][dir_ref][peak] * G_Counts_Sim_Det[best_pos_i_ref][dir_ref][peak]), 2));
            G_EfficiencyRatio_Sim_Det->AddPoint(detector, ratio);
            G_EfficiencyRatio_Sim_Det->SetPointError(G_EfficiencyRatio_Sim_Det->GetN() - 1, 0, ratio_err);


            // EXPERIMENTAL
            double ratio_exp = G_Counts_Exp_Det[peak]->GetY()[detector - 1] / G_Counts_Exp_Det[peak]->GetY()[detector_ref - 1];
            double ratio_exp_err = sqrt(pow(sqrt(G_Counts_Exp_Det[peak]->GetY()[detector - 1]) / G_Counts_Exp_Det[peak]->GetY()[detector_ref - 1], 2) + pow(G_Counts_Exp_Det[peak]->GetY()[detector - 1] * sqrt(G_Counts_Exp_Det[peak]->GetY()[detector_ref - 1]) / (G_Counts_Exp_Det[peak]->GetY()[detector_ref - 1] * G_Counts_Exp_Det[peak]->GetY()[detector_ref - 1]), 2));
            G_EfficiencyRatio_Exp_Det->AddPoint(detector, ratio_exp);
            G_EfficiencyRatio_Exp_Det->SetPointError(G_EfficiencyRatio_Exp_Det->GetN() - 1, 0, ratio_exp_err);

            // SIMULATION / EXPERIMENTAL RATIO

            for (int best_pos_i : Best_Position_Det[detector][peak])
            {
                cout << "Best position for detector " << detector << " : " << best_pos_i << endl;
                double ratio_sim_exp = G_Counts_Sim_Det[best_pos_i][dir][peak] / G_Counts_Exp_Det[peak]->GetY()[detector - 1];
                double ratio_sim_exp_err = sqrt(pow(sqrt(G_Counts_Sim_Det[best_pos_i][dir][peak]) / G_Counts_Exp_Det[peak]->GetY()[detector - 1], 2) + pow(G_Counts_Sim_Det[best_pos_i][dir][peak] * sqrt(G_Counts_Exp_Det[peak]->GetY()[detector - 1]) / (G_Counts_Exp_Det[peak]->GetY()[detector - 1] * G_Counts_Exp_Det[peak]->GetY()[detector - 1]), 2));
                G_EfficiencyRatio_SimExp_Det->AddPoint(detector, ratio_sim_exp);
                G_EfficiencyRatio_SimExp_Det->SetPointError(G_EfficiencyRatio_SimExp_Det->GetN() - 1, 0, ratio_sim_exp_err);
            }
            
        }

        TCanvas *c2 = new TCanvas(("EfficiencyRatio_Detector_Ratio_" + to_string(peak)).c_str(), ("EfficiencyRatio_Detector_Ratio_" + to_string(peak)).c_str(), 1920, 1080);
        c2->cd();
        TMultiGraph *mg_det_ratio = new TMultiGraph();
        TLegend *leg_det_ratio = new TLegend(0.7, 0.0, 1.0, 1.0);
        G_EfficiencyRatio_Sim_Det->SetMarkerColor(kRed);
        G_EfficiencyRatio_Sim_Det->SetMarkerStyle(20);
        G_EfficiencyRatio_Sim_Det->SetMarkerSize(2);
        G_EfficiencyRatio_Sim_Det->SetLineColor(kRed);
        G_EfficiencyRatio_Sim_Det->SetLineWidth(3);
        G_EfficiencyRatio_Exp_Det->SetMarkerColor(kBlack);
        G_EfficiencyRatio_Exp_Det->SetMarkerStyle(20);
        G_EfficiencyRatio_Exp_Det->SetMarkerSize(2);
        G_EfficiencyRatio_Exp_Det->SetLineColor(kBlack);
        G_EfficiencyRatio_Exp_Det->SetLineWidth(3);
        mg_det_ratio->Add(G_EfficiencyRatio_Sim_Det);
        mg_det_ratio->Add(G_EfficiencyRatio_Exp_Det);
        mg_det_ratio->Draw("ALP");
        c2->Write();

        TCanvas *c3 = new TCanvas(("EfficiencyRatio_SimExp_Detector_Ratio_" + to_string(peak)).c_str(), ("EfficiencyRatio_SimExp_Detector_Ratio_" + to_string(peak)).c_str(), 1920, 1080);
        c3->cd();
        G_EfficiencyRatio_SimExp_Det->Draw("AP");
        c3->Write();

    }

    TCanvas *c = new TCanvas("Chi2_2D", "Chi2_2D", 1920, 1080);
    Chi2_2D_yz->Draw("COLZ");
    Chi2_2D_yz->GetXaxis()->SetTitle("y");
    Chi2_2D_yz->GetYaxis()->SetTitle("z");
    Chi2_2D_yz->Write();
    
    TCanvas *c2 = new TCanvas("Chi2_2D_xy", "Chi2_2D_xy", 1920, 1080);
    Chi2_2D_xy->Draw("COLZ");
    Chi2_2D_xy->GetXaxis()->SetTitle("x");
    Chi2_2D_xy->GetYaxis()->SetTitle("y");
    Chi2_2D_xy->Write();
    
    for (int detector = 1; detector <= 8; detector++)
    {
        TCanvas *c3 = new TCanvas(("H_yztheta_" + to_string(detector)).c_str(), ("H_yztheta_" + to_string(detector)).c_str(), 1920, 1080);
        c3->cd();
        H_yztheta[detector]->Draw("BOX2Z");
        H_yztheta[detector]->GetXaxis()->SetTitle("y");
        H_yztheta[detector]->GetYaxis()->SetTitle("z");
        H_yztheta[detector]->GetZaxis()->SetTitle("theta");
        c3->Write();
    }
    
    

    output->Close();
    EXP_FILE->Close();  

    return 0;
}





















