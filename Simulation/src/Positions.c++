#include "Positions.hh"

// Purpose
// - taking simulated data with different positions of the catcher 
// - Computing the Eshift, coinc/single and proton efficiency for each position
// - Comparing to experimental data 

// Parameteres
// - x y z of the catcher (angle?)
// - Threshold and resolution of SiPMs (fixed)
// - SiPM multiplicity (fixed M>=3)


int main(int argc, char *argv[])
{
    FLAG2025 = true;
    InitDetectors("../Grouper/Config_Files/sample.pid");
    InitWindows();
    NUCLEUS = "32Ar";

    if (argc > 1)
        TYPE = argv[1];
    else
        TYPE = "Catcher";

    //  ################################################# //
    // OUTPUT FILE ///
    OUTPUT_FILE = MyTFile((DIR_ROOT_DATA_ANALYSED + "Positions_" + TYPE + "_" + to_string(YEAR) + ".root").c_str(), "RECREATE");
    //  ################################################# //


    //  ################################################# //
    // LOADING EXPERIMENTAL DATA ///
    // - M3
    // - threshold at 100keV
    // EXP
    Exp_FILE = MyTFile((DIR_ROOT_DATA_ANALYSED + "32Ar_" + to_string(YEAR) + "_analysed.root").c_str(), "READ");   

    // Eshift 
    TCanvas *cEshiftCompare_Corrected = (TCanvas *)Exp_FILE->Get("Eshift/Eshift_M3/Eshift_3");
    G_Exp_Eshift = nullptr;
    for (auto key : *cEshiftCompare_Corrected->GetListOfPrimitives())
    {
        if (string(key->ClassName()).find("TGraphErrors") != string::npos)
        {
            G_Exp_Eshift = (TGraphErrors *)key;
        }
    }
    if (G_Exp_Eshift == nullptr) Error("Eshift graph not found in experimental data");
    Info("Eshift graph loaded from experimental data", 1);

    // Coinc/Single
    TCanvas *cCoincSingle = (TCanvas *)Exp_FILE->Get("RatioCoinc_NoCoinc_Corrected_3");
    G_Exp_CoincSingle = nullptr;
    for (auto key : *cCoincSingle->GetListOfPrimitives())
    {
        if (string(key->ClassName()).find("TGraphErrors") != string::npos)
        {
            G_Exp_CoincSingle = (TGraphErrors *)key;
        }
    }
    if (G_Exp_CoincSingle == nullptr) Error("Coinc/Single graph not found in experimental data");
    Info("Coinc/Single graph loaded from experimental data", 1);

    // Proton efficiency (without cuts)
    Exp_FILE_grouped = MyTFile((DIR_ROOT_DATA_GROUPED + "merged_" + to_string(YEAR) + "_grouped.root").c_str(), "READ");
    LoadCalibrations();
    OUTPUT_FILE->cd();
    TCanvas *c_ProtonCounts = new TCanvas("c_ProtonCounts", "c_ProtonCounts", 1200, 800);
    c_ProtonCounts->Divide(5, 8);
    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliStrip(det))
        {
            TH1D* H_Exp = (TH1D *)Exp_FILE_grouped->Get(("Strip/Strip_Channel/" + detectorName[det] + "/Channel_Cleaned_" + detectorName[det]).c_str());
            if (H_Exp == nullptr) Error("Channel_Cleaned histogram not found for " + detectorName[det] + " in experimental grouped file");
            H_Exp->GetXaxis()->SetRangeUser(Silicon_Calibration[det]->Eval(WindowsMap[NUCLEUS][IAS[NUCLEUS]][det].first)*1e3, Silicon_Calibration[det]->Eval(WindowsMap[NUCLEUS][IAS[NUCLEUS]][det].second)*1e3);
            double integral_exp = H_Exp->Integral();
            G_Exp_ProtonCounts->AddPoint(det, integral_exp);
            G_Exp_ProtonCounts->SetPointError(G_Exp_ProtonCounts->GetN()-1, 0, sqrt(integral_exp));

            c_ProtonCounts->cd(ConvertBase5ToBase10(det) - 5);
            H_Exp->Draw("HIST");
        }
    }
    c_ProtonCounts->Write();
    Info("Proton counts graph loaded from experimental data", 1);
    //  ################################################# //

    PlottingExperimental();

    //  ################################################# //
    // LOADING SIMULATED DATA ///
    // - M3
    // - Threshold 100keV

    // path to files 
    // string path = "/run/media/local1/DATANEX/Samuel-G4/CatcherRelativePosition/";
    string path = "/run/media/local1/DATANEX/Samuel-G4/Positions/" + TYPE + "/";

    // Finding all interesting files in the directory
    vector<string> filenames = FindAllSimulatedFiles(path, to_string(YEAR));

    int counter = 0;
    // Looping over files and filing maps TGraphErrors
    for (auto f : filenames)
    {
        counter++;
        ProgressCounter(counter, filenames.size(), "File ");
        double x, y, z;
        string par_string = AnalyseSim(f, x, y, z);
        if (par_string == "")
            continue;

        // if (x != -2 && x != -1 && x != 0 && x != 1 && x != 2)
        //     continue;
        // if (y != 0 && y != 1 && y != 2 && y != 3 && y != 4)
        //     continue;

        if (z == 0.7 && y == 1.7)
            continue;
        
        PlottingSimulation(par_string);   
    
        // COMPUTING CHI2   
        // 
        ComputeRelativeEshift(par_string);
        Chi2_Eshift[par_string] = Chi2(G_Exp_Eshift, G_Sim_Eshift[par_string]);
        // Chi2_Eshift[par_string] = Sum(G_Sim_Eshift_Relative[par_string])*10;
        Chi2_CoincSingle[par_string] = Chi2(G_Exp_CoincSingle, G_Sim_CoincSingle[par_string]);
        // Chi2_CoincSingle[par_string] = Chi2(G_Exp_CoincSingle, G_Sim_CoincSingle[par_string]);
        Chi2_ProtonCounts[par_string] = Chi2(G_Exp_ProtonCounts, G_Sim_ProtonCounts[par_string]);
        
        ComputingRelativeProtonCounts(par_string);
        // delete G_Sim_ProtonCounts[par_string];
        // Chi2_ProtonCounts[par_string] = Sum(G_Sim_ProtonCounts_Relative[par_string]);
        Chi2_Total[par_string] = Chi2_Eshift[par_string] + 4*Chi2_CoincSingle[par_string] + 2*Chi2_ProtonCounts[par_string];

        

        if (x == 0.0)
        {
            G2_CoincSingle_YZ->AddPoint(y, z, Chi2_CoincSingle[par_string]);
            // G2_ProtonCounts_YZ->AddPoint(y, z, Chi2_ProtonCounts[par_string]);
            pair<double, double> value = MeanDistance(G_Sim_ProtonCounts_Relative[par_string]);
            G2_ProtonCounts_YZ->AddPoint(y, z, value.first);
            G2_ProtonCounts_YZ->SetPointError(G2_ProtonCounts_YZ->GetN() - 1, 0, 0, value.second);
            G2_Eshift_YZ->AddPoint(y, z, Chi2_Eshift[par_string]);
            // pair<double, double> value_eshift = MeanDistance(G_Sim_Eshift_Relative[par_string]);
            // G2_Eshift_YZ->AddPoint(y, z, value_eshift.first);
            // G2_Eshift_YZ->SetPointError(G2_Eshift_YZ->GetN() - 1, 0, 0, value_eshift.second);
            G2_Total_YZ->AddPoint(y, z, Chi2_Total[par_string]);
        }
        if ((z == 0.0 && TYPE == "Catcher") || (z == 0.0 && TYPE == "Detectors"))
        {
            G2_CoincSingle_XY->AddPoint(x, y, Chi2_CoincSingle[par_string]);
            // G2_ProtonCounts_XY->AddPoint(x, y, Chi2_ProtonCounts[par_string]);
            pair<double, double> value = MeanDistance(G_Sim_ProtonCounts_Relative[par_string]);
            G2_ProtonCounts_XY->AddPoint(x, y, value.first);
            G2_ProtonCounts_XY->SetPointError(G2_ProtonCounts_XY->GetN() - 1, 0, 0, value.second);
            G2_Eshift_XY->AddPoint(x, y, Chi2_Eshift[par_string]);
            pair<double, double> value_eshift = MeanDistance(G_Sim_Eshift_Relative[par_string]);
            G2_Eshift_XY->AddPoint(x, y, value_eshift.first);
            G2_Eshift_XY->SetPointError(G2_Eshift_XY->GetN() - 1, 0, 0, value_eshift.second);
            G2_Total_XY->AddPoint(x, y, Chi2_Total[par_string]);

            // for catcher and vakleur entirer de x and y 
            if (TYPE == "Catcher")
            {
                pair<double, double> value = MeanDiff(G_Exp_CoincSingle, G_Sim_CoincSingle[par_string]);
                G2_MeanDiff_CoincSingle_XY->AddPoint(x, y, abs(value.first));
                G2_MeanDiff_CoincSingle_XY->SetPointError(G2_MeanDiff_CoincSingle_XY->GetN() - 1, 0, 0, value.second);
            }
            
        }
        // PLOTTING
        PlottingExpvsSim(par_string);

        // Info(Form("File: %s, Eshift Chi2: %.2f, Coinc/Single Chi2: %.2f, Proton Counts Chi2: %.2f", par_string.c_str(), Chi2_Eshift[par_string], Chi2_CoincSingle[par_string], Chi2_ProtonCounts[par_string]));

    }

    

    
    Info("Ranking:");
    int top = 20;   
    // - best cases (sorting)
    // sorting map by value
    vector<pair<string, double>> Chi2_Eshift_vec(Chi2_Eshift.begin(), Chi2_Eshift.end());
    sort(Chi2_Eshift_vec.begin(), Chi2_Eshift_vec.end(), [](const pair<string, double>& a, const pair<string, double>& b) {
        return a.second < b.second;
    });
    
    PlottingExpvsSim(Chi2_Eshift_vec[0].first, "Eshift");
    Info(Form("Top %d Eshift Chi2 Results:", top), 1);
    for (int j = 0; j < top && j < Chi2_Eshift_vec.size(); j++) 
    {
        string par_string = Chi2_Eshift_vec[j].first;
        //regex to extract position from par_string
        std::regex rgx("x=([-+]?[0-9]*\\.?[0-9]+) y=([-+]?[0-9]*\\.?[0-9]+) z=([-+]?[0-9]*\\.?[0-9]+)");
        std::smatch match;
        if (std::regex_search(par_string, match, rgx) && match.size() == 4)
        {
            string x_pos = to_string(stod(match[1]));
            string y_pos = to_string(stod(match[2]));
            string z_pos = to_string(stod(match[3]));
            Info(Form("%d) X: %.2f mm, Y: %.2f mm, Z: %.2f mm, Chi2: %.6f", j+1, stod(x_pos), stod(y_pos), stod(z_pos), Chi2_Eshift_vec[j].second), 2);
        }       
    }

    vector<pair<string, double>> Chi2_CoincSingle_vec(Chi2_CoincSingle.begin(), Chi2_CoincSingle.end());
    sort(Chi2_CoincSingle_vec.begin(), Chi2_CoincSingle_vec.end(), [](const pair<string, double>& a, const pair<string, double>& b) {
        return a.second < b.second;
    });

    PlottingExpvsSim(Chi2_CoincSingle_vec[0].first, "CoincSingle");
    Info(Form("Top %d Coinc/Single Chi2 Results:", top), 1);
    for (int j = 0; j < top && j < Chi2_CoincSingle_vec.size(); j++) 
    {
        string par_string = Chi2_CoincSingle_vec[j].first;
        //regex to extract position from par_string
        std::regex rgx("x=([-+]?[0-9]*\\.?[0-9]+) y=([-+]?[0-9]*\\.?[0-9]+) z=([-+]?[0-9]*\\.?[0-9]+)");
        std::smatch match;
        if (std::regex_search(par_string, match, rgx) && match.size() == 4)
        {
            string x_pos = to_string(stod(match[1]));
            string y_pos = to_string(stod(match[2]));
            string z_pos = to_string(stod(match[3]));
            Info(Form("%d) X: %.2f mm, Y: %.2f mm, Z: %.2f mm, Chi2: %.6f", j+1, stod(x_pos), stod(y_pos), stod(z_pos), Chi2_CoincSingle_vec[j].second), 2);
        }       
    }


    vector<pair<string, double>> Chi2_ProtonCounts_vec(Chi2_ProtonCounts.begin(), Chi2_ProtonCounts.end());
    sort(Chi2_ProtonCounts_vec.begin(), Chi2_ProtonCounts_vec.end(), [](const pair<string, double>& a, const pair<string, double>& b) {
        return a.second < b.second;
    });

    PlottingExpvsSim(Chi2_ProtonCounts_vec[0].first, "ProtonCounts");
    Info(Form("Top %d Proton Counts Chi2 Results:", top), 1);
    for (int j = 0; j < top && j < Chi2_ProtonCounts_vec.size(); j++) 
    {
        string par_string = Chi2_ProtonCounts_vec[j].first;
        //regex to extract position from par_string
        std::regex rgx("x=([-+]?[0-9]*\\.?[0-9]+) y=([-+]?[0-9]*\\.?[0-9]+) z=([-+]?[0-9]*\\.?[0-9]+)");
        std::smatch match;
        if (std::regex_search(par_string, match, rgx) && match.size() == 4)
        {
            string x_pos = to_string(stod(match[1]));
            string y_pos = to_string(stod(match[2]));
            string z_pos = to_string(stod(match[3]));
            Info(Form("%d) X: %.2f mm, Y: %.2f mm, Z: %.2f mm, Chi2: %.6f", j+1, stod(x_pos), stod(y_pos), stod(z_pos), Chi2_ProtonCounts_vec[j].second), 2);
        }       
    }



    vector<pair<string, double>> Chi2_Total_vec(Chi2_Total.begin(), Chi2_Total.end());
    sort(Chi2_Total_vec.begin(), Chi2_Total_vec.end(), [](const pair<string, double>& a, const pair<string, double>& b) {
        return a.second < b.second;
    });
    PlottingExpvsSim(Chi2_Total_vec[0].first, "Total");
    Info(Form("Top %d Total Chi2 Results:", top), 1);
    for (int j = 0; j < top && j < Chi2_Total_vec.size(); j++) 
    {
        string par_string = Chi2_Total_vec[j].first;
        //regex to extract position from par_string
        std::regex rgx("x=([-+]?[0-9]*\\.?[0-9]+) y=([-+]?[0-9]*\\.?[0-9]+) z=([-+]?[0-9]*\\.?[0-9]+)");
        std::smatch match;
        if (std::regex_search(par_string, match, rgx) && match.size() == 4)
        {
            string x_pos = to_string(stod(match[1]));
            string y_pos = to_string(stod(match[2]));
            string z_pos = to_string(stod(match[3]));
            Info(Form("%d) X: %.2f mm, Y: %.2f mm, Z: %.2f mm, Chi2: %.6f", j+1, stod(x_pos), stod(y_pos), stod(z_pos), Chi2_Total_vec[j].second), 2);
        }       
    }


    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
    PlottingScans();
    PlottingStackedScan();
    // Plotting3D();

    

    OUTPUT_FILE->Close();
    return 0;


}