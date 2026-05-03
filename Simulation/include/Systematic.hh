#include "../../Grouper/include/Detectors.hh"
#include "../../Grouper/include/Analyser.hh"
#include "TKey.h"

#include <dirent.h>
#include <sys/stat.h>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>

TGraphErrors *G_Calibration_a[SIGNAL_MAX];
map<double, TGraphErrors *[SIGNAL_MAX]> G_Parameter;
TGraphErrors *G_Eshift[SIGNAL_MAX];
map<double, TGraphErrors *> G_Eshift_per_ParValue;
map<double, TGraphErrors *> G_DiffEshift_per_ParValue;
map<double, TGraphErrors *> G_Result_per_ParValue;
TGraphErrors *G_Result_per_Par = nullptr;
TFile *FINAL_File;
TH1D* MC_Stat;

double MC_STATISTICAL_ERROR;

double a_min;
double a_max;

map<string, double> Default = {
    {"Catcher_thickness", 610.},
    {"Bfield", 4.},
    {"DL", 92.},
    {"Silicon_Position_x", -0.30},
    {"Silicon_Position_y", 1.83},
    {"Silicon_Position_z", -0.66}
    };

map<string, double> Sigma = {
    {"Catcher_thickness", 25},
    {"Bfield", 1e-3},
    {"DL", 20.},
    {"Silicon_Position_x", 0.1},
    {"Silicon_Position_y", 0.1},
    {"Silicon_Position_z", 0.1}
    };

map<string, string> Name = {
    {"Catcher_thickness", "#varepsilon_{catcher} (nm)"},
    {"Bfield", "B (T)"},
    {"DL", "#varepsilon_{DL} (nm)"},
    {"Silicon_Position_x", "x (mm)"},
    {"Silicon_Position_y", "y (mm)"},
    {"Silicon_Position_z", "z (mm)"}
};

string GetParameterNameMacro(string par)
{
    if (par == "sx")
        return "/Beam/Sigma_X";
    else if (par == "sy")
        return "/Beam/Sigma_Y";
    else if (par == "y")
        return "/Beam/Y";
    else if (par == "x")
        return "/Beam/X";
    else
    {
        Warning("Unknown parameter: " + par);
        return par;
    }
}



double GetParameterValue(string par, string filename)
{
    Info("Getting parameter value for " + par + " from file: " + filename, 2);
    // oppening tfile same as filname without _analysed
    TFile *f = MyTFile(filename.substr(0, filename.find("_analysed")) + ".root", "READ", "Q");

    if (f == nullptr || f->IsZombie())
    {
        Warning("Could not open root file with the logs, Trying to get the paramater value in the filename");
        // if the file is not found, try to get the parameter value from the filename
        // take value of the folling 3 characters after the parameter name
        size_t pos = filename.find("_" + par);
        if (pos != string::npos)
        {
            // Extract the substring after the parameter name
            string value_str = filename.substr(pos + par.length() + 1);
            // Find the next underscore or end of string
            size_t end_pos = value_str.find('_');
            if (end_pos == string::npos)
            {
                end_pos = value_str.length();
            }
            // Convert to double
            return stod(value_str.substr(0, end_pos));
        }
    }

    // looping on all element of file TString
    TIter next(f->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)next()))
    {
        TString name = key->GetName();
        if (name.Contains(GetParameterNameMacro(par).c_str()))
        {
            return stod(name(name.Index(" ") + 1, name.Length() - name.Index(" ") - 1).Data());
        }
    }

    f->Close();
    Error("Parameter " + par + " not found in file: " + filename);
    return -1;
}

void Search_files(std::string &dir_search, std::string filename_base, std::vector<std::string> &filenames)
{
    if (filename_base.find("Catcher_thickness") != std::string::npos)
        filename_base = "Al";
    if (filename_base.find("Silicon_Position") != std::string::npos)
    {
        dir_search = "/run/media/local1/DATANEX/Samuel-G4/Systematics/Silicon_Position/";
        filename_base = "x";
    }
    if (filename_base.find("Bfield") != std::string::npos)
        filename_base = "B";
    Info("Searching for files in directory: " + dir_search + " with base name: " + filename_base, 1);
    DIR *dir = opendir(dir_search.c_str());
    if (!dir)
    {
        Error("Error opening directory: " + dir_search);
        return;
    }

    struct dirent *entry;
    while ((entry = readdir(dir)) != nullptr)
    {
        std::string name = entry->d_name;

        // cout << "Checking file: " << name << endl;

        // Skip '.' and '..'
        // if (name == "." || name == "..") continue;

        std::string full_path = dir_search + "/" + name;

        // cout << "Full path: " << full_path << endl;

        struct stat file_stat;
        if (stat(full_path.c_str(), &file_stat) == 0 && S_ISREG(file_stat.st_mode))
        {
            // if (name.find(filename_base) != std::string::npos && name.find("_analysed") != std::string::npos)
            if (name.find(filename_base) != std::string::npos && name.find("_result") != std::string::npos)
            {
                filenames.push_back(full_path);
            }
        }
    }

    closedir(dir);
}

pair<double, double> ComputeEshift(TFile *file, int det)
{
    TCanvas *c = nullptr;
    TH1D *h_single = nullptr, *h_coinc = nullptr, *h_nocoinc = nullptr;
    if (string(file->GetName()).find("analysed_new") != std::string::npos)
    {
        c = (TCanvas *)file->Get(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_Coinc").c_str());
        h_single = (TH1D *)c->GetPrimitive(("Silicon_Detector_Energy_Deposit_SINGLE_" + detectorName[det] + "_All").c_str());
        h_coinc = (TH1D *)c->GetPrimitive(("Silicon_Detector_Energy_Deposit_COINC_" + detectorName[det] + "_All").c_str());
        h_nocoinc = (TH1D *)c->GetPrimitive(("Silicon_Detector_Energy_Deposit_NOCOINC_" + detectorName[det] + "_All").c_str());
    }
    else if (string(file->GetName()).find("result") != std::string::npos)
    {
        c = (TCanvas *)file->Get(("H_Single_" + detectorName[det]).c_str());
        h_single = (TH1D *)c->GetPrimitive(("H_Single_" + detectorName[det]).c_str());
        h_coinc = (TH1D *)c->GetPrimitive(("H_Coinc_" + detectorName[det]).c_str());
        h_nocoinc = (TH1D *)c->GetPrimitive(("H_NoCoinc_" + detectorName[det]).c_str());
    }
    else
    {
        Error("Unknown file format for file: " + string(file->GetName()));
    }

    if (h_single == nullptr || h_coinc == nullptr || h_nocoinc == nullptr)
    {
        Error("Histograms not found for detector: " + detectorName[det]);
    }

    return ComputeEshift(det, h_single, h_coinc, h_nocoinc);
}

void InitAnalysis()
{
    // COMPUTING THE LINEAR FUNCTION ESHIFT vs a
    map<double, TFile *> files;

    // Find files
    string in = ".root";
    string path = "/run/media/local1/DATANEX/Samuel-G4/Systematics/a_Calibration/";
    vector<string> filenames;
    Search_files(path, in, filenames);

    for (auto &filename : filenames)
    {
        // pos last /
        size_t pos = filename.find("_a");
        if (pos != string::npos)
        {
            // Extract the substring after "a_"
            string value_str = filename.substr(pos + 2);
            // Find the next underscore or end of string
            size_t end_pos = value_str.find('_');
            if (end_pos == string::npos)
            {
                end_pos = value_str.length();
            }
            // Convert to double
            double a_value = stod(value_str.substr(0, end_pos));
            if (a_value < a_min || a_value > a_max)
                continue;

            files[a_value] = MyTFile(filename, "READ");
        }
    }

    // Open and fill TGraphs
    for (auto &entry : files)
    {
        double a = entry.first;
        TFile *file = entry.second;
        if (file == nullptr || file->IsZombie())
        {
            Warning("Could not open file for a = " + to_string(a));
            continue;
        }

        for (int det = 1; det < SIGNAL_MAX; det++)
        {
            if (IsDetectorSiliStrip(det))
            {
                pair<double, double> Eshift = ComputeEshift(file, det);
                if (G_Calibration_a[det] == nullptr)
                {
                    G_Calibration_a[det] = new TGraphErrors();
                    G_Calibration_a[det]->SetName(("G_Calibration_a_" + detectorName[det]).c_str());
                }

                G_Calibration_a[det]->AddPoint(Eshift.first, a);
                G_Calibration_a[det]->SetPointError(G_Calibration_a[det]->GetN() - 1, Eshift.second, 1e-10);
            }
        }
    }

    // fitting each detector curve
    // for (int det = 1; det < SIGNAL_MAX; det++)
    // {
    //     if (IsDetectorSiliStrip(det) && G_Calibration_a[det] != nullptr)
    //     {
    //         TFitResultPtr r = G_Calibration_a[det]->Fit("pol1", "SE");
    //         if (r->IsValid())
    //             Warning("Fit failed for detector: " + detectorName[det]);
    //     }
    // }
}

void ComputingParameter(string sys)
{
    Info("Computing paramater : " + sys, 1);
    string dir_search = "/run/media/local1/DATANEX/Samuel-G4/Systematics/" + sys + "/";
    // creatoing vector of filename containing filenmaebase
    vector<string> filenames;
    Search_files(dir_search, sys, filenames);

    filenames.push_back(dir_search + "../a_Calibration/32Ar_IAS_" + to_string(YEAR) + "_a1.0000_b0.0000_result.root");

    // looping on all files
    for (string filename : filenames)
    {
        // get par value
        // double par_value = GetParameterValue(sys, filename);
        // from sys to next _ or .root
        double par_value;

        try
        {
            if (filename == filenames[filenames.size() - 1])
            {
                cout << "Filename: " << filename << endl;
                if (sys != "DL")
                    par_value = Default[sys];
                else 
                    par_value = Default[sys];
            }
            else
            {
               // par_value = stod(filename.substr(filename.find("_" + sys) + sys.length() + 1, filename.find(".root") - (filename.find("_" + sys) + sys.length() + 1)));
                // par_value = stod(filename.substr(filename.find("_" + sys) + sys.length() + 1, filename.find("_") - (filename.find("_" + sys) + sys.length() + 1)));
                if (sys == "Bfield")
                {
                    par_value = stod(filename.substr(filename.find("_B") + 2, filename.find("T") - (filename.find("_B") + 2)));
                    // excluding 
                    if (abs(par_value - Default[sys]) < 0.01)
                        continue;   
                }
                else if (sys == "Catcher_thickness")
                {
                    par_value = stod(filename.substr(filename.find("_Al") + 3, filename.find("_Mylar") - (filename.find("_Al") + 3))) + stod(filename.substr(filename.find("_Mylar") + 6, filename.find("_") - (filename.find("_Mylar") + 6)));
                }
                else if (sys.find("Silicon_Position") != std::string::npos)
                {
                    double x = stod(filename.substr(filename.find("_x") + 2, filename.find("_y") - (filename.find("_x") + 2)));
                    double y = stod(filename.substr(filename.find("_y") + 2, filename.find("_z") - (filename.find("_y") + 2)));
                    double z = stod(filename.substr(filename.find("_z") + 2, filename.find("_Rx") - (filename.find("_z") + 2)));
                    if (x != Default["Silicon_Position_x"] && y == Default["Silicon_Position_y"] && z == Default["Silicon_Position_z"] && sys == "Silicon_Position_x")
                        par_value = x;
                    else if (x == Default["Silicon_Position_x"] && y != Default["Silicon_Position_y"] && z == Default["Silicon_Position_z"] && sys == "Silicon_Position_y")
                        par_value = y;
                    else if (x == Default["Silicon_Position_x"] && y == Default["Silicon_Position_y"] && z != Default["Silicon_Position_z"] && sys == "Silicon_Position_z")
                        par_value = z;
                    else 
                    {
                        Warning("Could not determine which silicon position parameter is varied in file: " + filename);
                        continue;
                    }                    
                }
                else 
                {
                    par_value = stod(filename.substr(filename.find("_" + sys) + sys.length() + 1, filename.find("_") - (filename.find("_" + sys) + sys.length() + 1)));
                }
            }

            Info("Parameter value for " + sys + ": " + to_string(par_value), 2);
        }
        catch (const std::exception &e)
        {
            Warning("Could not convert parameter value for file: " + filename);
            continue;
        }

        Info("Parameter value for " + sys + ": " + to_string(par_value), 2);

        // openning file
        TFile *file = MyTFile(filename.c_str(), "READ");
        if (!file || file->IsZombie())
        {
            Warning("Could not open file: " + filename);
            continue;
        }
        // looping on all detectors
        for (int det = 1; det < SIGNAL_MAX; det++)
        {
            if (IsDetectorSiliStrip(det) && G_Calibration_a[det] != nullptr)
            {
                pair<double, double> Eshift = ComputeEshift(file, det);
                // cout << "Eshift for detector " << detectorName[det] << ": " << Eshift.first << " +/- " << Eshift.second << "      MC stat error: " << MC_STATISTICAL_ERROR[det] << endl;
                // Eshift.second = sqrt(Eshift.second * Eshift.second - MC_STATISTICAL_ERROR[det] * MC_STATISTICAL_ERROR[det]);

                if (G_Parameter[par_value][det] == nullptr)
                {
                    G_Parameter[par_value][det] = new TGraphErrors();
                    G_Parameter[par_value][det]->SetName(("G_Parameter_" + sys + "_" + to_string(par_value) + "_" + detectorName[det]).c_str());
                }
                // adding data to graph for each det
                G_Parameter[par_value][det]->AddPoint(Eshift.first, Eshift.first);
                G_Parameter[par_value][det]->SetPointError(G_Parameter[par_value][det]->GetN() - 1, Eshift.second, Eshift.second);

                if (G_Eshift[det] == nullptr)
                {
                    G_Eshift[det] = new TGraphErrors();
                    G_Eshift[det]->SetName(("G_Eshift_" + detectorName[det]).c_str());
                }
                // adding data to graph for each det
                G_Eshift[det]->AddPoint(par_value, Eshift.first);
                G_Eshift[det]->SetPointError(G_Eshift[det]->GetN() - 1, 0, Eshift.second);
            }
        }
        file->Close();
    }
}



void InitSampling()
{
    Start("Sampling", 1);
    string dir_search = "/run/media/local1/DATANEX/Samuel-G4/Systematics/Sampling/";
    // creatoing vector of filename containing filenmaebase
    vector<string> filenames;
    Search_files(dir_search, "Default", filenames);

    vector<pair<double, double>> value[SIGNAL_MAX];

    for (string filename : filenames)
    {
        TFile *file = MyTFile(filename.c_str(), "READ", "Q");
        if (!file || file->IsZombie())
        {
            Warning("Could not open file: " + filename);
            continue;
        }
        for (int det = 1; det < SIGNAL_MAX; det++)
        {
            if (IsDetectorSiliStrip(det))
            {
                pair<double, double> Eshift = ComputeEshift(file, det);
                // Adding to calib also
                if (G_Calibration_a[det] == nullptr)
                {
                    G_Calibration_a[det] = new TGraphErrors();
                    G_Calibration_a[det]->SetName(("G_Calibration_a_" + detectorName[det]).c_str());
                }

                // G_Calibration_a[det]->AddPoint(Eshift.first, 1.0);
                // G_Calibration_a[det]->SetPointError(G_Calibration_a[det]->GetN() - 1, Eshift.second, 1e-10);              
            }
        }
        file->Close();
    }

    // fit calibration curve with sampling points 
    for (int det = 1; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliStrip(det) && G_Calibration_a[det] != nullptr)
        {
            TFitResultPtr r = G_Calibration_a[det]->Fit("pol1", "SE");
            if (r->IsValid())
                Warning("Fit failed for detector: " + detectorName[det]);
        }
    }

    // for each a=1.0 get a_fit from Eshift
    for (string filename : filenames)
    {
        TFile *file = MyTFile(filename.c_str(), "READ", "Q");
        if (!file || file->IsZombie())
        {
            Warning("Could not open file: " + filename);
            continue;
        }
        TGraphErrors *g = new TGraphErrors();
        for (int det = 1; det < SIGNAL_MAX; det++)
        {
            if (IsDetectorSiliStrip(det))
            {
                pair<double, double> Eshift = ComputeEshift(file, det);
                g->AddPoint(det, G_Calibration_a[det]->GetFunction("pol1")->Eval(Eshift.first));
                g->SetPointError(g->GetN() - 1, 0, abs(G_Calibration_a[det]->GetFunction("pol1")->Derivative(Eshift.first)) * Eshift.second);
            }
        }
        file->Close();

        g->Fit("pol0", "QE");
        if (MC_Stat == nullptr)        {
            MC_Stat = new TH1D("MC_Stat", "MC_Stat", 100, 0.99, 1.01);
        }
        MC_Stat->Fill(g->GetFunction("pol0")->GetParameter(0));
    }

    MC_STATISTICAL_ERROR = MC_Stat->GetStdDev();
    return;
}

void AnalyseParameter(string sys)
{
    for (int det = 1; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliStrip(det) && G_Calibration_a[det] != nullptr)
        {
            Info("Analysing parameter: " + sys + " for detector: " + detectorName[det], 1);
            // TF1 *fit = new TF1("pol1", "[0]*x + (1. - [0] * [1])");
            // // find x for y = 1 in G_Calibration_a[det] and fix parameter 2 to this value
            // for (int i = 0; i < G_Calibration_a[det]->GetN(); i++)
            // {
            //     double x, y;
            //     G_Calibration_a[det]->GetPoint(i, x, y);
            //     if (y == 1.)
            //     {
            //         fit->FixParameter(1, x);
            //         Info("Fixing parameter 1 to " + to_string(x) + " for detector: " + detectorName[det], 2);
            //         break;
            //     }
            // }
            for (auto &pair : G_Parameter)
            {
                double par_value = pair.first;
                TGraphErrors *graph = pair.second[det];
                if (graph != nullptr)
                {
                    double Eshift, dummy, Eshift_err;
                    graph->GetPoint(0, Eshift, dummy);
                    Eshift_err = graph->GetErrorY(0);
                    double a = G_Calibration_a[det]->GetFunction("pol1")->Eval(Eshift);
                    // double err_low_a = G_Calibration_a[det]->GetFunction("pol1")->Eval(Eshift - Eshift_err) - a;
                    // double err_high_a = G_Calibration_a[det]->GetFunction("pol1")->Eval(Eshift + Eshift_err) - a;
                    double err_a = abs(G_Calibration_a[det]->GetFunction("pol1")->Derivative(Eshift)) * Eshift_err;

                    // cout << "Parameter: " << sys << ", Value: " << par_value << ", Detector: " << detectorName[det] << endl;
                    // cout << "    Eshift: " << Eshift << ", a: " << a << endl;
                    // cout << "    Eshift error: " << Eshift_err << ", a error: " << err_a << endl;
                    // creating graph for each param and par_value

                    if (G_Result_per_ParValue[par_value] == nullptr)
                    {
                        G_Result_per_ParValue[par_value] = new TGraphErrors();
                        G_Result_per_ParValue[par_value]->SetName(("G_Result_" + sys + "_" + detectorName[det]).c_str());
                    }
                    // adding data to graph for each det
                    G_Result_per_ParValue[par_value]->AddPoint(det, a);
                    G_Result_per_ParValue[par_value]->SetPointError(G_Result_per_ParValue[par_value]->GetN() - 1, 0, err_a);                
                }
            }
        }
    }

    // Computing a for each paramater and paramater value in G_Result_per_Par
    for (auto &pair : G_Parameter)
    {
        double par_value = pair.first;
        G_Result_per_ParValue[par_value]->Fit("pol0", "QE");
        double a = G_Result_per_ParValue[par_value]->GetFunction("pol0")->GetParameter(0);
        double err_a = G_Result_per_ParValue[par_value]->GetFunction("pol0")->GetParError(0);
        err_a = MC_STATISTICAL_ERROR;// sqrt(err_a * err_a - MC_STATISTICAL_ERROR * MC_STATISTICAL_ERROR);

        if (G_Result_per_Par == nullptr)
        {
            G_Result_per_Par = new TGraphErrors();
            G_Result_per_Par->SetName(("G_Result_" + sys).c_str());
        }
        // adding data to graph for each param
        G_Result_per_Par->AddPoint(par_value, a);
        G_Result_per_Par->SetPointError(G_Result_per_Par->GetN() - 1, 0, err_a);
    }
}

void WriteHistograms(string sys)
{
    FINAL_File->cd();
    for (int det = 1; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliStrip(det) && G_Calibration_a[det] != nullptr)
        {
            TCanvas *canvas = new TCanvas(("Calibration_a_" + detectorName[det]).c_str(), ("Calibration_a_" + detectorName[det]).c_str(), 800, 600);
            G_Calibration_a[det]->SetTitle(("Calibration a for " + detectorName[det]).c_str());
            G_Calibration_a[det]->GetYaxis()->SetTitle("a");
            G_Calibration_a[det]->GetXaxis()->SetTitle("Eshift [keV]");
            G_Calibration_a[det]->SetMarkerStyle(20);
            G_Calibration_a[det]->SetMarkerSize(1);
            G_Calibration_a[det]->Draw("AP");
            TLatex *latex = new TLatex(0.15, 0.8, Form("#chi^{2}_{#nu} = %.2f", G_Calibration_a[det]->GetFunction("pol1")->GetChisquare() / G_Calibration_a[det]->GetFunction("pol1")->GetNDF()));
            latex->SetNDC();
            latex->SetTextSize(0.04);
            latex->Draw("SAME");
            canvas->Write();
        }
    }

    for (auto &pair : G_Parameter)
    {
        TCanvas *canvas_param = new TCanvas(("Parameter_" + sys + "_" + to_string(pair.first)).c_str(), ("Parameter_" + sys + "_" + to_string(pair.first)).c_str(), 800, 600);
        TGraphErrors *graph = G_Result_per_ParValue[pair.first];
        if (graph != nullptr)
        {
            graph->SetTitle(("Parameter " + sys + " for value " + to_string(pair.first)).c_str());
            graph->GetYaxis()->SetTitle("a");
            graph->GetXaxis()->SetTitle("Detector");
            graph->SetMarkerStyle(20);
            graph->SetMarkerSize(1);
            graph->Draw("AP");
            TLatex *latex = new TLatex(0.15, 0.8, ("Fitted a: " + to_string(graph->GetFunction("pol0")->GetParameter(0)) + " #pm " + to_string(graph->GetFunction("pol0")->GetParError(0))).c_str());
            latex->SetNDC();
            latex->SetTextSize(0.04);
            latex->Draw("SAME");
            canvas_param->Write();
        }
    }

    TCanvas *canvas_result = new TCanvas(("Result_" + sys).c_str(), ("Result_" + sys).c_str(), 800, 600);
    TGraphErrors *result_graph = G_Result_per_Par;

    // adding delta a from default value to all the points of the graph
    // double shift;
    // for (int i = 0; i < result_graph->GetN(); i++)
    // {
    //     double x, y;
    //     result_graph->GetPoint(i, x, y);
    //     if (x == Default[sys])
    //     {
    //         shift = y - 1.;
    //         break;
    //     }
    // }
    // 


    if (result_graph != nullptr)
    {
        result_graph->SetTitle(("Result for parameter " + sys).c_str());
        result_graph->GetYaxis()->SetTitle("a");
        result_graph->GetXaxis()->SetTitle(Name[sys].c_str());
        result_graph->SetMarkerStyle(20);
        result_graph->SetMarkerSize(1);
        
        double maxmax = result_graph->GetYaxis()->GetXmax();
        double minlow = result_graph->GetYaxis()->GetXmin();
        double maxi = std::max(abs(maxmax-1.), abs(1.-minlow)); 
        result_graph->GetYaxis()->SetRangeUser(1. - maxi*1.2, 1. + maxi*1.2);
        result_graph->Draw("AP");

        TF1 *fit = new TF1("pol1", "pol1");
        TFitResultPtr r = result_graph->Fit(fit, "QSE");
        if (r->IsValid())
            Warning("Fit failed for result graph of parameter: " + sys);
        TLatex *latex = new TLatex(0.15, 0.8, Form("#chi^{2}_{#nu} = %.2f", r.Get()->Chi2() / r.Get()->Ndf()));
        latex->SetNDC();
        latex->SetTextSize(0.04);
        latex->Draw("SAME");

        // shifting to have default a= 1.
        double shift = 1 - fit->Eval(Default[sys]);
        for (int i = 0; i < result_graph->GetN(); i++)
        {
            double x, y;
            result_graph->GetPoint(i, x, y);
            result_graph->SetPoint(i, x, y + shift);
        }
        result_graph->Fit(fit, "QE");

        double sigma = abs(fit->Derivative(Default[sys])) * Sigma[sys];
        TLatex *latex_fit = new TLatex(0.15, 0.75, Form("#Delta a = %.6f ", sigma));
        latex_fit->SetNDC();
        latex_fit->SetTextSize(0.04);
        latex_fit->Draw("SAME");
        TLine *hline = new TLine(Default[sys], result_graph->GetYaxis()->GetXmin(), Default[sys], 1.);
        hline->SetLineColor(kRed);
        hline->SetLineStyle(2);
        hline->Draw("SAME");
        TLine *vline = new TLine(result_graph->GetXaxis()->GetXmin(), 1., Default[sys], 1.);
        vline->SetLineColor(kRed);
        vline->SetLineStyle(2);
        vline->Draw("SAME");
        TLine *vline_psigma = new TLine(Default[sys]+Sigma[sys], result_graph->GetYaxis()->GetXmin(), Default[sys]+Sigma[sys], 1. + sigma);
        vline_psigma->SetLineColor(kRed);
        vline_psigma->SetLineStyle(2);
        vline_psigma->Draw("SAME");
        TLine *vline_msigma = new TLine(Default[sys]-Sigma[sys], result_graph->GetYaxis()->GetXmin(), Default[sys]-Sigma[sys], 1. - sigma);
        vline_msigma->SetLineColor(kRed);
        vline_msigma->SetLineStyle(2);
        vline_msigma->Draw("SAME");
        TLine *hline_psigma = new TLine(result_graph->GetXaxis()->GetXmin(), 1. + sigma, Default[sys]+Sigma[sys], 1. + sigma);
        hline_psigma->SetLineColor(kRed);
        hline_psigma->SetLineStyle(2);
        hline_psigma->Draw("SAME");
        TLine *hline_msigma = new TLine(result_graph->GetXaxis()->GetXmin(), 1. - sigma, Default[sys]-Sigma[sys], 1. - sigma);     
        hline_msigma->SetLineColor(kRed);
        hline_msigma->SetLineStyle(2);
        hline_msigma->Draw("SAME");
        canvas_result->Write();
    }

    TCanvas *MC_Stat_Canvas = new TCanvas("MC_Stat", "MC_Stat", 800, 600);
    MC_Stat->SetTitle("MC Statistical error");
    MC_Stat->GetXaxis()->SetTitle("a");
    MC_Stat->GetYaxis()->SetTitle("Counts");
    MC_Stat->Draw("HIST");
    MC_Stat_Canvas->Write();

}