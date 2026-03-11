#include "Systematic_MC.hh"


vector<string> GetFilenames(const char* dir, const char* pattern)
{
    vector<string> filenames;
    DIR* dp;
    struct dirent* dirp;
    if ((dp = opendir(dir)) == NULL)
    {
        Error("GetFilenames - Couldn't open the directory");
    }

    while ((dirp = readdir(dp)) != NULL)
    {
        string filename = string(dirp->d_name);
        if (filename.find(pattern) != string::npos)
        {
            filenames.push_back(string(dir) + filename);
        }
    }
    closedir(dp);
    return filenames;
}

double *GetParameter(const string& filename, vector<string> par)
{
    // fin par element in fielname and take the 6 character in the string for its value
    double *param = new double[par.size()];
    for (size_t i = 0; i < par.size(); i++)
    {
        size_t pos = filename.find(par[i]);
        if (pos != string::npos)
        {
            param[i] = stod(filename.substr(pos + par[i].size(), 6));
        }
        else
        {
            Error("GetParameter - Parameter not found in filename");
        }
    }

    return param;
}

TGraphErrors *ComputeParameter(string filename)
{
    TGraphErrors *G_Eshift;
    TFile *file = MyTFile(filename.c_str(), "READ");
    if (!file || file->IsZombie())
    {
        Warning("Could not open file: " + filename);
    }
    bool changed = false;
    // looping on all detectors
    for (int det = 1; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliStrip(det) && G_Calibration_a[det] != nullptr)
        {

            TCanvas *c = (TCanvas *)file->Get(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_Coinc").c_str());
            if (c == nullptr)
            {
                Warning("Canvas not found for detector: " + detectorName[det]);
                continue;
            }

            TH1D *h_single = (TH1D *)c->GetPrimitive(("Silicon_Detector_Energy_Deposit_SINGLE_" + detectorName[det] + "_All").c_str());
            TH1D *h_coinc = (TH1D *)c->GetPrimitive(("Silicon_Detector_Energy_Deposit_COINC_" + detectorName[det] + "_All").c_str());
            TH1D *h_nocoinc = (TH1D *)c->GetPrimitive(("Silicon_Detector_Energy_Deposit_NOCOINC_" + detectorName[det] + "_All").c_str());
            if (h_single == nullptr || h_coinc == nullptr || h_nocoinc == nullptr)
            {
                Warning("Histograms not found for detector: " + detectorName[det]);
                continue;
            }

            pair<double, double> Eshift = ComputeEshift(det, h_single, h_coinc, h_nocoinc);
            
            // adding data to the graph for each det
            G_Eshift->AddPoint(det, Eshift.first);
            G_Eshift->SetPointError(G_Eshift->GetN() - 1, 0, Eshift.second);
        }
    }
    file->Close();
    return G_Eshift;
}

int main()
{
    FLAG2024 = true;
    NUCLEUS = "32Ar";
    InitDetectors("../Grouper/Config_Files/sample.pid");
    InitWindows();

    string sys = "beam_pos";
    vector<string> par = {"x", "y", "sx", "sy"};    
    map<string, TH1D*> H_Parameter;
    H_Parameter["x"] = new TH1D("H_x", "H_x", 300, -3, 3);
    H_Parameter["y"] = new TH1D("H_y", "H_y", 300, -3, 3);
    H_Parameter["sx"] = new TH1D("H_sx", "H_sx", 200, 0, 2);
    H_Parameter["sy"] = new TH1D("H_sy", "H_sy", 200, 0, 2);

    TH1D *H_a = new TH1D("H_a", "H_a", 400, 0.8, 1.2);

    FINAL_File = MyTFile(("Sytematic_MC_" + sys + "_" + to_string(YEAR) + ".root").c_str(), "RECREATE");

    // Computing the linear fucntion Eshift(a) with the 3 points
    InitAnalysis();

    // Creating az vector<string> of the filename with "*analysed.root" in the directory dir
    vector<string> filenames = GetFilenames("/data333/lecanuet/Result/Sampling/", "analysed.root");
    vector<double*> Paramater_vec;
    for (auto filename : filenames)
    {
        TFile *f = MyTFile(filename.c_str(), "READ");
        Paramater_vec.push_back(GetParameter(filename, par));
        for (size_t i = 0; i < par.size(); i++)
        {
            H_Parameter[par[i]]->Fill(Paramater_vec.back()[i]);
        }
        TGraphErrors *G_Eshift = ComputeParameter(filename);

        TGraphErrors *G_a = new TGraphErrors();
        for (int i = 0; i < G_Eshift->GetN(); i++)
        {
            double det, Eshift, Eshift_err;
            G_Eshift->GetPoint(i, det, Eshift);
            Eshift_err = G_Eshift->GetErrorY(i);

            double a = G_Calibration_a[(int)det]->GetFunction("pol1")->Eval(Eshift);
            double err_low_a = abs(G_Calibration_a[(int)det]->GetFunction("pol1")->Eval(Eshift - Eshift_err) - a);
            double err_high_a = abs(G_Calibration_a[(int)det]->GetFunction("pol1")->Eval(Eshift + Eshift_err) - a);
            double err_a = err_low_a < err_high_a ? err_high_a : err_low_a;

            G_a->AddPoint(det, a);
            G_a->SetPointError(G_a->GetN() - 1, 0, err_a);
        }
        G_a->Fit("pol0", "QE");

        double a_global = G_a->GetFunction("pol0")->GetParameter(0);
        double a_global_err = G_a->GetFunction("pol0")->GetParError(0);

        H_a->Fill(a_global);
    }

    // Writting
    FINAL_File->cd();
    for (auto &pair : H_Parameter)
    {
        pair.second->Write();
    }
    H_a->Write();

    FINAL_File->Close();
    return 0;
}