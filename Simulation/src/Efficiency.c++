#include "Efficiency.hh"
default_random_engine generator;


int main()
{
    FLAG2024 = true;
    InitDetectors("../Grouper/Config_Files/sample.pid");
    InitWindows();

    if (YEAR == 2025)
    {
        DIR_DATA_HDD += "../../2024_DATA/DETECTOR_DATA/";
    }

    int strip_ref = 1;
    bool display_ref = true;

    const double miny = .0;
    const double maxy = 5.0;
    const double minx = -3.0;
    const double maxx = 5.0;
    const double step = 1.0;
    
    map<double, TGraphErrors *[3]>G_EfficiencyRatio_Up;
    map<double, TGraphErrors *[3]>G_EfficiencyRatio_Down;

    map<int, map<double, TGraphErrors *[3]>> G_EfficiencyRatio_Up_PerStrip;
    map<int, map<double, TGraphErrors *[3]>> G_EfficiencyRatio_Down_PerStrip; 

    string NUCLEUS = "32Ar";

    vector<int> Peaks = {1, 3, 14};

    for (double value_y = miny; value_y <= maxy; value_y += step)
    {
        
        G_EfficiencyRatio_Up[value_y][0] = nullptr;
        G_EfficiencyRatio_Up[value_y][1] = nullptr;
        G_EfficiencyRatio_Up[value_y][2] = nullptr;
        G_EfficiencyRatio_Down[value_y][0] = nullptr;
        G_EfficiencyRatio_Down[value_y][1] = nullptr;
        G_EfficiencyRatio_Down[value_y][2] = nullptr;

        for (double value_x = minx; value_x <= maxx; value_x += step)
        {
            if (value_x == -0.5)
                continue;

            // if ((value_x != 3.0 && value_y != 3.0) || (value_y != 4.0 && value_x != 3.0))
            //     continue;

            // if (value_x != value_y)
            //     continue;
            // only 2 digit
            ostringstream oss_x;
            oss_x << fixed << setprecision(1) << value_x;
            string str_value_x = oss_x.str();

            ostringstream oss_y;
            oss_y << fixed << setprecision(1) << value_y;
            string str_value_y = oss_y.str();

            
            // TFile *SIM_FILE = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/05-18/32Arx"+str_value_x+"_y"+str_value_y+"_z0_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");

            TFile *SIM_FILE = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/05-18/32Arx0.0_y" + str_value_y + "_z" + str_value_x + "_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
            if (str_value_x == "0.0" && SIM_FILE == nullptr)
            {
                SIM_FILE = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/05-18/32Arx0.0_y" + str_value_y + "_z0_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
            }

            if (value_y == 4.0 && value_x == 3.0)
            {
                SIM_FILE = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/05-18/new/32Arx0.0_y3.0_z3.0_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
            }

            if (SIM_FILE == nullptr)
            {
                continue;
            }

            // int Scale = 10e6/(int)SIM_FILE->Get("Tree")->GetEntries(); // not used because ratio

            
            TH1D *H_UP[SILI_SIZE+1] = {nullptr};
            TH1D *H_DOWN[SILI_SIZE+1] = {nullptr};


            // adding histogram from file in hemisphere cynlindric symetry
            for (int det = 0; det < SIGNAL_MAX; det++)
            {
                if (IsDetectorSiliStrip(det))
                {
                    if (GetDetector(det) <= 4)
                    {
                        if (H_UP[GetDetectorChannel(det)] == nullptr)
                            H_UP[GetDetectorChannel(det)] = (TH1D *)SIM_FILE->Get(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_All").c_str());
                        else
                            H_UP[GetDetectorChannel(det)]->Add((TH1D *)SIM_FILE->Get(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_All").c_str()));
                    }
                    else if (GetDetector(det) >= 5)
                    {
                        if (H_DOWN[GetDetectorChannel(det)] == nullptr)
                            H_DOWN[GetDetectorChannel(det)] = (TH1D *)SIM_FILE->Get(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_All").c_str());
                        else
                            H_DOWN[GetDetectorChannel(det)]->Add((TH1D *)SIM_FILE->Get(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_All").c_str()));
                    }
                }
            }

            // compute ratio fo each peak between strips
            for (int strip = 1; strip < SILI_SIZE; strip++)
            {
                // Info("Strip " + to_string(strip), 2);
                for (int i_peak = 0; i_peak < Peaks.size(); i_peak++)
                {
                    // Info("Peak " + to_string(Peaks[i_peak]), 2);
                    int peak = Peaks[i_peak];

                    H_UP[strip]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][peak][10 + strip].first, WindowsMap[NUCLEUS][peak][10 + strip].second);
                    double integral1 = H_UP[strip]->Integral();
                    H_UP[strip_ref]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][peak][10+strip_ref].first, WindowsMap[NUCLEUS][peak][10+strip_ref].second);
                    double integral2 = H_UP[strip_ref]->Integral();

                    if (G_EfficiencyRatio_Up[value_y][i_peak] == nullptr)
                        G_EfficiencyRatio_Up[value_y][i_peak] = new TGraphErrors();
                
                    G_EfficiencyRatio_Up[value_y][i_peak]->AddPoint(value_x, integral1 / integral2);
                    G_EfficiencyRatio_Up[value_y][i_peak]->SetPointError(G_EfficiencyRatio_Up[value_y][i_peak]->GetN() - 1, 0, sqrt(integral1) / integral2 + integral1 * sqrt(integral2) / (integral2 * integral2));  

                    if (G_EfficiencyRatio_Up_PerStrip[strip][value_y][i_peak] == nullptr)
                        G_EfficiencyRatio_Up_PerStrip[strip][value_y][i_peak] = new TGraphErrors();

                    G_EfficiencyRatio_Up_PerStrip[strip][value_y][i_peak]->AddPoint(value_x, integral1 / integral2);
                    G_EfficiencyRatio_Up_PerStrip[strip][value_y][i_peak]->SetPointError(G_EfficiencyRatio_Up_PerStrip[strip][value_y][i_peak]->GetN() - 1, 0, sqrt(integral1) / integral2 + integral1 * sqrt(integral2) / (integral2 * integral2));

                    H_DOWN[strip]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][peak][50 + strip].first, WindowsMap[NUCLEUS][peak][50 + strip].second);
                    double integral3 = H_DOWN[strip]->Integral();
                    H_DOWN[strip_ref]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][peak][50+strip_ref].first, WindowsMap[NUCLEUS][peak][50+strip_ref].second);
                    double integral4 = H_DOWN[strip_ref]->Integral();


                    if (G_EfficiencyRatio_Down[value_y][i_peak] == nullptr)
                        G_EfficiencyRatio_Down[value_y][i_peak] = new TGraphErrors();

                    G_EfficiencyRatio_Down[value_y][i_peak]->AddPoint(value_x, integral3 / integral4);
                    G_EfficiencyRatio_Down[value_y][i_peak]->SetPointError(G_EfficiencyRatio_Down[value_y][i_peak]->GetN() - 1, 0, sqrt(integral3) / integral4 + integral3 * sqrt(integral4) / (integral4 * integral4));

                    if (G_EfficiencyRatio_Down_PerStrip[strip][value_y][i_peak] == nullptr)
                        G_EfficiencyRatio_Down_PerStrip[strip][value_y][i_peak] = new TGraphErrors();
                    
                    G_EfficiencyRatio_Down_PerStrip[strip][value_y][i_peak]->AddPoint(value_x, integral3 / integral4);
                    G_EfficiencyRatio_Down_PerStrip[strip][value_y][i_peak]->SetPointError(G_EfficiencyRatio_Down_PerStrip[strip][value_y][i_peak]->GetN() - 1, 0, sqrt(integral3) / integral4 + integral3 * sqrt(integral4) / (integral4 * integral4));
                }
            }
        }
    }

    // SAME FOR EXPERIMENTAL // 
    TFile *EXP_FILE = MyTFile((DIR_ROOT_DATA_CALIBRATED + "Calibrated_"+ to_string(YEAR) + ".root").c_str(), "READ");
    if (EXP_FILE == nullptr)
        return 0;
    
    TH1D *H_EXP_UP[SIGNAL_MAX] = {nullptr};
    TH1D *H_EXP_DOWN[SIGNAL_MAX] = {nullptr};

    pair<double, double> Value_Up[SIGNAL_MAX][3];
    pair<double, double> Value_Down[SIGNAL_MAX][3];

    // adding histogram from file in hemisphere cylindric symetry
    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliStrip(det))
        {
            if (GetDetector(det) <= 4)
            {
                H_EXP_UP[det] = (TH1D *)EXP_FILE->Get((detectorName[det] + "/32Ar/h_sub" + detectorName[det]).c_str());
            }
            else if (GetDetector(det) >= 5)
            {
                H_EXP_DOWN[det] = (TH1D *)EXP_FILE->Get((detectorName[det] + "/32Ar/h_sub" + detectorName[det]).c_str());
            }
        }
    }

    for (int det = 1; det < SIGNAL_MAX; det++)
    {
        if (!IsDetectorSiliStrip(det))
            continue;
        // Info("Strip " + to_string(det), 2);
        for (int i_peak = 0; i_peak < Peaks.size(); i_peak++)
        {
            int peak = Peaks[i_peak];

            if (GetDetector(det) <= 4)
            {
                H_EXP_UP[det]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][peak][det].first, WindowsMap[NUCLEUS][peak][det].second);
                double integral1 = H_EXP_UP[det]->Integral();
                H_EXP_UP[GetDetector(det)*10 + strip_ref]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][peak][GetDetector(det)*10 + strip_ref].first, WindowsMap[NUCLEUS][peak][GetDetector(det)*10 + strip_ref].second);
                double integral2 = H_EXP_UP[GetDetector(det)*10 + strip_ref]->Integral();

                Value_Up[det][i_peak].first = integral1 / integral2;
                Value_Up[det][i_peak].second = sqrt(integral1) / integral2 + integral1 * sqrt(integral2) / (integral2 * integral2);
            }
            else if (GetDetector(det) >= 5)
            {
                H_EXP_DOWN[det]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][peak][det].first, WindowsMap[NUCLEUS][peak][det].second);
                double integral3 = H_EXP_DOWN[det]->Integral();
                H_EXP_DOWN[GetDetector(det)*10 + strip_ref]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][peak][GetDetector(det)*10 + strip_ref].first, WindowsMap[NUCLEUS][peak][GetDetector(det)*10 + strip_ref].second);
                double integral4 = H_EXP_DOWN[GetDetector(det)*10 + strip_ref]->Integral();

                Value_Down[det][i_peak].first = integral3 / integral4;
                Value_Down[det][i_peak].second = sqrt(integral3) / integral4 + integral3 * sqrt(integral4) / (integral4 * integral4);
            }
        }
    }

    TFile *output = MyTFile(("EfficiencyRatio_" + to_string(YEAR) + ".root").c_str(), "RECREATE");
    


    // for (int i_peak = 0; i_peak < Peaks.size(); i_peak++)
    // {

    //     TCanvas *c_Up = new TCanvas(("EfficiencyRatio_Up_" + to_string(Peaks[i_peak])).c_str(), ("EfficiencyRatio_Up_" + to_string(Peaks[i_peak])).c_str(), 800, 600);
    //     TMultiGraph *mg_Up = new TMultiGraph();
    //     TCanvas *c_Down = new TCanvas(("EfficiencyRatio_Down_" + to_string(Peaks[i_peak])).c_str(), ("EfficiencyRatio_Down_" + to_string(Peaks[i_peak])).c_str(), 800, 600);
    //     TMultiGraph *mg_Down = new TMultiGraph();
    //     TLegend *leg_Up = new TLegend(0.7, 0.7, 0.9, 0.9);
    //     TLegend *leg_Down = new TLegend(0.7, 0.7, 0.9, 0.9);
    //     int counter = 0; 

        
    //     //DISPLAY TGRAPHERROR FROM SIMS
    //     for (double value_y = miny; value_y <= maxy; value_y += step)
    //     {
            
    //         if (G_EfficiencyRatio_Up[value_y][i_peak] != nullptr)
    //         {
    //             c_Up->cd();
    //             G_EfficiencyRatio_Up[value_y][i_peak]->SetMarkerStyle(20);
    //             G_EfficiencyRatio_Up[value_y][i_peak]->SetMarkerSize(1);
    //             G_EfficiencyRatio_Up[value_y][i_peak]->SetMarkerColor(counter+1);
    //             G_EfficiencyRatio_Up[value_y][i_peak]->GetXaxis()->SetTitle("x");
    //             G_EfficiencyRatio_Up[value_y][i_peak]->GetYaxis()->SetTitle("Efficiency Ratio");
    //             mg_Up->Add(G_EfficiencyRatio_Up[value_y][i_peak]);
    //             leg_Up->AddEntry(G_EfficiencyRatio_Up[value_y][i_peak], ("y = " + to_string(value_y)).c_str(), "p");
                
    //             c_Down->cd();
    //             G_EfficiencyRatio_Down[value_y][i_peak]->SetMarkerStyle(20);
    //             G_EfficiencyRatio_Down[value_y][i_peak]->SetMarkerSize(1);
    //             G_EfficiencyRatio_Down[value_y][i_peak]->SetMarkerColor(counter + 1);
    //             G_EfficiencyRatio_Down[value_y][i_peak]->GetXaxis()->SetTitle("x");
    //             G_EfficiencyRatio_Down[value_y][i_peak]->GetYaxis()->SetTitle("Efficiency Ratio");
    //             mg_Down->Add(G_EfficiencyRatio_Down[value_y][i_peak]);
    //             leg_Down->AddEntry(G_EfficiencyRatio_Down[value_y][i_peak], ("y = " + to_string(value_y)).c_str(), "p");

    //             counter++;
    //         }
    //     }

    //     c_Up->cd();
    //     mg_Up->Draw("AP");
    //     c_Down->cd();
    //     mg_Down->Draw("AP");

    //     for (int det = 1; det < SIGNAL_MAX; det++)
    //     {
    //         if (!IsDetectorSiliStrip(det))
    //             continue;
    //         // DISPLAY EXPERIMENTAL

    //         if (GetDetector(det) <= 4)
    //         {

    //             c_Up->cd();
    //             TGraphErrors *line_Exp_Up = new TGraphErrors();
    //             line_Exp_Up->SetPoint(0, minx, Value_Up[det][i_peak].first);
    //             line_Exp_Up->SetPointError(0, 0, Value_Up[det][i_peak].second);
    //             line_Exp_Up->SetPoint(1, maxx, Value_Up[det][i_peak].first);
    //             line_Exp_Up->SetPointError(1, 0, Value_Up[det][i_peak].second);
    //             line_Exp_Up->SetLineColor(kRed);
    //             line_Exp_Up->SetLineWidth(2);
    //             line_Exp_Up->SetFillStyle(3005);
    //             line_Exp_Up->SetFillColor(kRed);
    //             line_Exp_Up->Draw("E3 SAME");
    //         }
    //         else if (GetDetector(det) >= 5)
    //         {

    //             c_Down->cd();
    //             TGraphErrors *line_Exp_Down = new TGraphErrors();
    //             line_Exp_Down->SetPoint(0, minx, Value_Down[det][i_peak].first);
    //             line_Exp_Down->SetPointError(0, 0, Value_Down[det][i_peak].second);
    //             line_Exp_Down->SetPoint(1, maxx, Value_Down[det][i_peak].first);
    //             line_Exp_Down->SetPointError(1, 0, Value_Down[det][i_peak].second);
    //             line_Exp_Down->SetLineColor(kRed);
    //             line_Exp_Down->SetLineWidth(2);
    //             line_Exp_Down->SetFillStyle(3005);
    //             line_Exp_Down->SetFillColor(kRed);
    //             line_Exp_Down->Draw("E3 SAME");
    //         }
    //     }

    //     // leg_Up->AddEntry(line_Exp_Up, "Experimental", "l");
    //     c_Up->cd();
    //     leg_Up->Draw("SAME");

    //     // leg_Down->AddEntry(line_Exp_Down, "Experimental", "l");
    //     c_Down->cd();
    //     leg_Down->Draw("SAME");

    //     c_Up->Write();
    //     c_Down->Write();
    // }





    // // ############ FOR EACH STRIP

    // TDirectory *dir_detector[SILI_SIZE+1] = {nullptr};

    // for (int i_peak = 0; i_peak < Peaks.size(); i_peak++)
    // {
    //     for (int strip = 2; strip < SILI_SIZE; strip++)
    //     {
    //         // if (dir_detector[strip] == nullptr) dir_detector[strip] = output->mkdir(("Up_Strip_" + to_string(strip)).c_str());
            

    //         TCanvas *c_Up = new TCanvas(("EfficiencyRatio_Up_Strip_" + to_string(strip) + "_" + to_string(Peaks[i_peak])).c_str(), ("EfficiencyRatio_Up_Strip_" + to_string(strip) + "_" + to_string(Peaks[i_peak])).c_str(), 800, 600);
    //         TMultiGraph *mg_Up = new TMultiGraph();
    //         TCanvas *c_Down = new TCanvas(("EfficiencyRatio_Down_Strip_" + to_string(strip) + "_" + to_string(Peaks[i_peak])).c_str(), ("EfficiencyRatio_Down_Strip_" + to_string(strip) + "_" + to_string(Peaks[i_peak])).c_str(), 800, 600);
    //         TMultiGraph *mg_Down = new TMultiGraph();
    //         TLegend *leg_Up = new TLegend(0.7, 0.7, 0.9, 0.9);
    //         TLegend *leg_Down = new TLegend(0.7, 0.7, 0.9, 0.9);
    //         int counter = 0;

    //         // DISPLAY TGRAPHERROR FROM SIMS
    //         for (double value_y = miny; value_y <= maxy; value_y += step)
    //         {
    //             G_EfficiencyRatio_Up_PerStrip[strip][value_y][i_peak]->SetMarkerStyle(20);
    //             G_EfficiencyRatio_Up_PerStrip[strip][value_y][i_peak]->SetMarkerSize(1);
    //             G_EfficiencyRatio_Up_PerStrip[strip][value_y][i_peak]->SetMarkerColor(counter + 1);
    //             G_EfficiencyRatio_Up_PerStrip[strip][value_y][i_peak]->GetXaxis()->SetTitle("x");
    //             G_EfficiencyRatio_Up_PerStrip[strip][value_y][i_peak]->GetYaxis()->SetTitle("Efficiency Ratio");
    //             mg_Up->Add(G_EfficiencyRatio_Up_PerStrip[strip][value_y][i_peak]);
    //             leg_Up->AddEntry(G_EfficiencyRatio_Up_PerStrip[strip][value_y][i_peak], ("y = " + to_string(value_y)).c_str(), "p");

    //             G_EfficiencyRatio_Down_PerStrip[strip][value_y][i_peak]->SetMarkerStyle(20);
    //             G_EfficiencyRatio_Down_PerStrip[strip][value_y][i_peak]->SetMarkerSize(1);
    //             G_EfficiencyRatio_Down_PerStrip[strip][value_y][i_peak]->SetMarkerColor(counter + 1);
    //             G_EfficiencyRatio_Down_PerStrip[strip][value_y][i_peak]->GetXaxis()->SetTitle("x");
    //             G_EfficiencyRatio_Down_PerStrip[strip][value_y][i_peak]->GetYaxis()->SetTitle("Efficiency Ratio");
    //             mg_Down->Add(G_EfficiencyRatio_Down_PerStrip[strip][value_y][i_peak]);
    //             leg_Down->AddEntry(G_EfficiencyRatio_Down_PerStrip[strip][value_y][i_peak], ("y = " + to_string(value_y)).c_str(), "p");

    //             counter++;
    //         }

    //         c_Up->cd();
    //         mg_Up->GetXaxis()->SetRangeUser(minx, maxx);
    //         mg_Up->GetYaxis()->SetRangeUser(0.5, 1.1);
    //         mg_Up->Draw("AP");

    //         c_Down->cd();
    //         mg_Down->GetXaxis()->SetRangeUser(minx, maxx);
    //         mg_Down->GetYaxis()->SetRangeUser(0.5, 1.1);    
    //         mg_Down->Draw("AP");

    //         // DISPLAY EXPERIMENTAL
    //         for (int detector = 1; detector <= 8; detector++)
    //         {
    //             int det = 10 * detector + strip;

    //             if (GetDetector(det) <= 4)
    //             {
    //                 c_Up->cd();
    //                 TGraphErrors *line_Exp_Up = new TGraphErrors();
    //                 line_Exp_Up->SetPoint(0, minx, Value_Up[det][i_peak].first);
    //                 line_Exp_Up->SetPointError(0, 0, Value_Up[det][i_peak].second);
    //                 line_Exp_Up->SetPoint(1, maxx, Value_Up[det][i_peak].first);
    //                 line_Exp_Up->SetPointError(1, 0, Value_Up[det][i_peak].second);
    //                 line_Exp_Up->SetLineColor(kRed);
    //                 line_Exp_Up->SetLineWidth(2);
    //                 line_Exp_Up->SetFillStyle(3005);
    //                 line_Exp_Up->SetFillColor(kRed);
    //                 line_Exp_Up->Draw("E3 SAME");
                    
    //             }
    //             else if (GetDetector(det) >= 5)
    //             {

    //                 c_Down->cd();
    //                 TGraphErrors *line_Exp_Down = new TGraphErrors();
    //                 line_Exp_Down->SetPoint(0, minx, Value_Down[det][i_peak].first);
    //                 line_Exp_Down->SetPointError(0, 0, Value_Down[det][i_peak].second);
    //                 line_Exp_Down->SetPoint(1, maxx, Value_Down[det][i_peak].first);
    //                 line_Exp_Down->SetPointError(1, 0, Value_Down[det][i_peak].second);
    //                 line_Exp_Down->SetLineColor(kRed);
    //                 line_Exp_Down->SetLineWidth(2);
    //                 line_Exp_Down->SetFillStyle(3005);
    //                 line_Exp_Down->SetFillColor(kRed);
    //                 line_Exp_Down->Draw("E3 SAME");
    //             }
    //         }

    //         c_Up->cd();
    //         leg_Up->Draw("SAME");
    //         c_Up->Write();

    //         c_Down->cd();
    //         leg_Down->Draw("SAME");
    //         c_Down->Write();
    //     }
    // }

    // ############ FOR EACH DETECTOR

    for (int i_peak = 0; i_peak < Peaks.size(); i_peak++)
    {
        for (int detector = 1; detector <= 8; detector++)
        {
            // Info("Detector " + to_string(detector), 2); 
            TCanvas *c = new TCanvas(("EfficiencyRatio_D"+to_string(detector) + "_" + to_string(Peaks[i_peak])).c_str(), ("EfficiencyRatio_D"+to_string(detector) + "_" + to_string(Peaks[i_peak])).c_str(), 1920, 1080);
            TPad *p1 = new TPad("p1", "p1", 0.0, 0.0, 0.7, 1.0);
            p1->Draw();
            p1->cd();
            TMultiGraph *mg_det = new TMultiGraph();
            TLegend *leg_det = new TLegend(0.7, 0.0, 1.0, 1.0);
            int counter = 0;


            map<double, map<double, TGraphErrors *>> ge_new;
            for (int strip = 1; strip < SILI_SIZE; strip++)
            {
                // if (strip == strip_ref) 
                //     continue;
                // Info("Strip " + to_string(strip), 3);
                for (double value_y = miny; value_y <= maxy; value_y += step)
                {
                    TGraphErrors *ge;
                    if (GetDetector(detector * 10 + strip) <= 4)
                        ge = G_EfficiencyRatio_Up_PerStrip[strip][value_y][i_peak];
                    else
                        ge = G_EfficiencyRatio_Down_PerStrip[strip][value_y][i_peak];

                    if (ge == nullptr)
                        continue;

                    for (int ix = 0; ix < ge->GetN(); ix++)
                    {
                        if (strip == 5) counter++;
                        double x_value, Ratio;
                        ge->GetPoint(ix, x_value, Ratio);
                        c->cd();
                        if (strip == 1) ge_new[x_value][value_y] = new TGraphErrors();
                        if (!display_ref && strip == strip_ref) continue;
                        ge_new[x_value][value_y]->AddPoint(strip, Ratio);
                        ge_new[x_value][value_y]->SetPointError(ge_new[x_value][value_y]->GetN()-1, 0, ge->GetErrorY(ix));
                        if (strip == 5)
                        {
                            ge_new[x_value][value_y]->SetMarkerStyle(20);
                            ge_new[x_value][value_y]->SetMarkerSize(1);
                            ge_new[x_value][value_y]->SetDrawOption("ALP");
                            ge_new[x_value][value_y]->SetMarkerColor(counter + 1);
                            ge_new[x_value][value_y]->SetLineColor(counter + 1);
                            mg_det->Add(ge_new[x_value][value_y]);
                        }
                    }
                }
                // Info("Strip " + to_string(strip), 3);
            }

            // c->cd();
            TGraphErrors *line_Exp = new TGraphErrors();
            for (int strip = 1; strip < SILI_SIZE; strip++)
            {
                if (!display_ref && strip == strip_ref) continue;

                if (GetDetector(detector * 10 + strip) <= 4)
                {
                    line_Exp->AddPoint(strip, Value_Up[detector * 10 + strip][i_peak].first);
                    line_Exp->SetPointError(line_Exp->GetN()-1, 0, Value_Up[detector * 10 + strip][i_peak].second);
                }
                else if (GetDetector(detector * 10 + strip) >= 5)
                {
                    line_Exp->AddPoint(strip, Value_Down[detector * 10 + strip][i_peak].first);
                    line_Exp->SetPointError(line_Exp->GetN()-1, 0, Value_Down[detector * 10 + strip][i_peak].second);
                }
            }
            line_Exp->SetMarkerColor(kRed);
            line_Exp->SetMarkerStyle(20);
            line_Exp->SetMarkerSize(2);
            line_Exp->SetLineColor(kRed);
            line_Exp->SetLineWidth(2);
            // line_Exp->Draw("AP");
            mg_det->Add(line_Exp);

            // chi2 and legend
            // ordered map of chi2 and string legend
            map<double, string> chi2_map;
            map<double, pair<double, double>> chi2_map_value;
            for (double value_y = miny; value_y <= maxy; value_y += step)
            {
                for (double value_x = minx; value_x <= maxx; value_x += step)
                {
                    TGraphErrors *ge = ge_new[value_x][value_y];
                    if (ge != nullptr)
                    {
                        ostringstream oss;
                        oss << fixed << setprecision(1) << value_y;
                        string str_value_y = oss.str();
                        ostringstream oss_x;
                        oss_x << fixed << setprecision(1) << value_x;
                        string str_value_x = oss_x.str();
                        double chi2 = Chi2(ge, line_Exp);
                        ostringstream oss_chi2;
                        oss_chi2 << fixed << setprecision(2) << chi2;
                        string str_value_chi2 = oss_chi2.str();
                        // chi2_map[chi2] = ("y = " + str_value_y + " x = " + str_value_x + "  (Chi2 = " + str_value_chi2 + ")");
                        chi2_map[chi2] = ("y = " + str_value_y + " z = " + str_value_x + "  (Chi2 = " + str_value_chi2 + ")");
                        chi2_map_value[chi2] = make_pair(value_x, value_y);
                        // leg_det->AddEntry(ge, chi2_map[chi2].c_str(), "p");
                    }
                }
            }

            for (auto it = chi2_map.begin(); it != chi2_map.end(); ++it)
            {
                leg_det->AddEntry(ge_new[chi2_map_value[it->first].first][chi2_map_value[it->first].second], it->second.c_str(), "p");
                if (it->first == chi2_map.begin()->first)
                {
                    ge_new[chi2_map_value[it->first].first][chi2_map_value[it->first].second]->SetMarkerSize(2);
                    ge_new[chi2_map_value[it->first].first][chi2_map_value[it->first].second]->SetLineWidth(2);
                    ge_new[chi2_map_value[it->first].first][chi2_map_value[it->first].second]->SetMarkerStyle(21);
                }
            }

            p1->cd();
            mg_det->Draw("ALP");

            c->cd();
            leg_det->Draw("SAME");
            c->Write();
        }
    }

    output->Close();

    return 0;
}
