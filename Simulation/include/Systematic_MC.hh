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
unordered_map<string, map<double ,TGraphErrors *[SIGNAL_MAX]>>G_Parameter;
TGraphErrors *G_Eshift[SIGNAL_MAX];
map<string, map<double, TGraphErrors *>> G_Eshift_per_ParValue;
map<string, map<double, TGraphErrors *>> G_DiffEshift_per_ParValue;
map<string, map<double, TGraphErrors *>> G_Result_per_ParValue;
map<string, TGraphErrors *> G_Result_per_Par;
TFile* FINAL_File;

void InitAnalysis()
{
    // COMPUTING THE LINEAR FUNCTION ESHIFT vs a

    vector<double> a_values = {1, 0, -1}; // a = 1, 0, -1
    map<double, TFile*> files;

    // a = 1
    files[1] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/07-17/32Ar_sx0.5_sy0.5_y0.0_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");

    // a = 0
    files[0] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/07-17/32Ar_sx0.5_sy0.5_y0.0_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");

    // a = -1
    files[-1] = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/07-17/32Ar_sx0.5_sy0.5_y0.0_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");


    for (int det = 1; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliStrip(det))
        {

            for (double a : a_values)
            {
                // getting data in the file
                TCanvas *c = (TCanvas *)files[a]->Get(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_Coinc").c_str());
                if (c == nullptr)
                {
                    Warning("Canvas not found for detector: " + detectorName[det] + " with a = " + to_string(a));
                    continue;
                }

                TH1D *h_single = (TH1D *)c->GetPrimitive(("Silicon_Detector_Energy_Deposit_SINGLE_" + detectorName[det] + "_All").c_str());
                TH1D *h_coinc = (TH1D *)c->GetPrimitive(("Silicon_Detector_Energy_Deposit_COINC_" + detectorName[det] + "_All").c_str());
                TH1D *h_nocoinc = (TH1D *)c->GetPrimitive(("Silicon_Detector_Energy_Deposit_NOCOINC_" + detectorName[det] + "_All").c_str());
                if (h_single == nullptr || h_coinc == nullptr || h_nocoinc == nullptr)
                {
                    Warning("Histograms not found for detector: " + detectorName[det] + " with a = " + to_string(a));
                    continue;
                }

                // computing eshift
                pair<double, double> Eshift = ComputeEshift(det, h_single, h_coinc, h_nocoinc);
                if (G_Calibration_a[det] == nullptr)
                {
                    G_Calibration_a[det] = new TGraphErrors();
                    G_Calibration_a[det]->SetName(("G_Calibration_a_" + detectorName[det]).c_str());
                }

                // adding data to grpah for each det
                G_Calibration_a[det]->AddPoint(Eshift.first + a - 1, a);
                G_Calibration_a[det]->SetPointError(G_Calibration_a[det]->GetN() - 1, Eshift.second, 0);
            }
            G_Calibration_a[det]->Fit("pol1", "QE");
        }
    }
}

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
        Warning ("Unknown parameter: " + par);
        return par;
    }
}

double GetParameterValue(string par, string filename)
{
    Info("Getting parameter value for " + par + " from file: " + filename, 2);
    //oppening tfile same as filname without _analysed
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
            string value_str = filename.substr(pos + par.length()+1);
            // Find the next underscore or end of string
            size_t end_pos = value_str.find('_');
            if (end_pos == string::npos) {
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

void Search_files(const std::string& dir_search, const std::string& filename_base, std::vector<std::string>& filenames)
{
    Info("Searching for files in directory: " + dir_search + " with base name: " + filename_base, 1);
    DIR* dir = opendir(dir_search.c_str());
    if (!dir) {
        Error("Error opening directory: " + dir_search);
        return;
    }

    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        std::string name = entry->d_name;

        cout << "Checking file: " << name << endl;

        // Skip '.' and '..'
        if (name == "." || name == "..") continue;

        std::string full_path = dir_search + "/" + name;

        cout << "Full path: " << full_path << endl;

        struct stat file_stat;
        if (stat(full_path.c_str(), &file_stat) == 0 && S_ISREG(file_stat.st_mode)) {
            if (name.find(filename_base) != std::string::npos && name.find("_analysed") != std::string::npos && name.find("CS0_CSP0_CV1_CVP1") != std::string::npos) {
                filenames.push_back(full_path);
            }
        }
    }

    closedir(dir);
}

void ComputingParameter(vector<string> parameters)
{
    for (const string &param : parameters)
    {
        Info("Computing paramater : " + param, 1);
        string filename_base = "32Ar_sx0.5_sy0.5_y";
        string dir_search = DIR_DATA_HDD + "../SIMULATED_DATA/07-17/";
        // creatoing vector of filename containing filenmaebase
        vector<string> filenames;
        Search_files(dir_search, filename_base, filenames);

        // looping on all files
        for (string filename : filenames)
        {

            // get par value 
            double par_value = GetParameterValue(param, filename);

            if (par_value == 3.0)
            {
                filename = DIR_DATA_HDD + "../SIMULATED_DATA/07-17/32Ar_sx0.5_sy0.5_y3.0_CS0_CSP0_CV1_CVP1_1um_analysed.root";
            }

            // if (par_value == 0.0)
            // {
            //     filename = DIR_DATA_HDD + "../SIMULATED_DATA/07-17/32Ar_sx0.5_sy0.5_y0.0_CS0_CSP0_CV1_CVP1_analysed.root";
            // }

            // if (par_value == 0.0)
            // {
            //     filename = DIR_DATA_HDD + "../SIMULATED_DATA/07-17/32Ar_sx0.5_sy0.5_y0.0_CS0_CSP0_CV1_CVP1_proton0_analysed.root";
            // }

        
            Info("Parameter value for " + param + ": " + to_string(par_value), 2);

            // openning file 
            TFile *file = MyTFile(filename.c_str(), "READ");
            if (!file || file->IsZombie())
            {
                Warning("Could not open file: " + filename);
                continue;
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

                    if (G_Parameter[param][par_value][det] == nullptr)
                    {
                        G_Parameter[param][par_value][det] = new TGraphErrors();
                        G_Parameter[param][par_value][det]->SetName(("G_Parameter_" + param + "_" + to_string(par_value) + "_" + detectorName[det]).c_str());
                    }
                    // adding data to graph for each det
                    G_Parameter[param][par_value][det]->AddPoint(Eshift.first, Eshift.first);
                    G_Parameter[param][par_value][det]->SetPointError(G_Parameter[param][par_value][det]->GetN() - 1, Eshift.second, Eshift.second);
                
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
}

void AnalyseParameter(vector<string> parameters)
{

    for (int det = 1; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliStrip(det) && G_Calibration_a[det] != nullptr)
        {
            G_Calibration_a[det]->Fit("pol1", "QE");
            for (const string &param : parameters)
            {           
                for (auto &pair : G_Parameter[param])
                {
                    double par_value = pair.first;
                    TGraphErrors *graph = pair.second[det];
                    if (graph != nullptr)
                    {
                        double Eshift, dummy, Eshift_err;
                        graph->GetPoint(0, Eshift, dummy);
                        Eshift_err = graph->GetErrorY(0);
                        double a = G_Calibration_a[det]->GetFunction("pol1")->Eval(Eshift);
                        double err_low_a = G_Calibration_a[det]->GetFunction("pol1")->Eval(Eshift - Eshift_err) - a;
                        double err_high_a = G_Calibration_a[det]->GetFunction("pol1")->Eval(Eshift + Eshift_err) - a;

                        double err_a = err_low_a < err_high_a ? err_high_a : err_low_a;

                        // cout << "Parameter: " << param << ", Value: " << par_value << ", Detector: " << detectorName[det] <<endl;
                        // cout << "    Eshift: " << Eshift << ", a: " << a << endl;
                        // cout << "    Eshift error: " << Eshift_err << ", a error: " << err_a << endl;
                        // creating graph for each param and par_value

                        if (G_Result_per_ParValue[param][par_value] == nullptr)
                        {
                            G_Result_per_ParValue[param][par_value] = new TGraphErrors();
                            G_Result_per_ParValue[param][par_value]->SetName(("G_Result_" + param + "_" + detectorName[det]).c_str());
                        }
                        // adding data to graph for each det
                        G_Result_per_ParValue[param][par_value]->AddPoint(det, a);
                        G_Result_per_ParValue[param][par_value]->SetPointError(G_Result_per_ParValue[param][par_value]->GetN() - 1, 0, err_a);

                        if (G_Eshift_per_ParValue[param][par_value] == nullptr)
                        {
                            G_Eshift_per_ParValue[param][par_value] = new TGraphErrors();
                            G_Eshift_per_ParValue[param][par_value]->SetName(("G_Eshift_" + param + "_" + detectorName[det]).c_str());
                            
                            G_DiffEshift_per_ParValue[param][par_value] = new TGraphErrors();
                            G_DiffEshift_per_ParValue[param][par_value]->SetName(("G_DiffEshift_" + param + "_" + detectorName[det]).c_str());
                        }
                        // adding data to graph for each det
                        G_Eshift_per_ParValue[param][par_value]->AddPoint(det, Eshift);
                        G_Eshift_per_ParValue[param][par_value]->SetPointError(G_Eshift_per_ParValue[param][par_value]->GetN() - 1, 0, Eshift_err);

                        if (par_value != 0)
                        {
                            // adding data to graph for each det
                            cout << Eshift_err << endl;
                            G_DiffEshift_per_ParValue[param][par_value]->AddPoint(det, Eshift - G_Eshift_per_ParValue[param][0.0]->GetPointY(G_DiffEshift_per_ParValue[param][par_value]->GetN()));
                            G_DiffEshift_per_ParValue[param][par_value]->SetPointError(G_DiffEshift_per_ParValue[param][par_value]->GetN() - 1, 0, Eshift_err + G_Eshift_per_ParValue[param][0.0]->GetErrorY(G_DiffEshift_per_ParValue[param][par_value]->GetN()-1));
                        }
                    }
                }
            }
        }
    }

    // Computing a for each paramater and paramater value in G_Result_per_Par
    for (const string &param : parameters)
    {
        for (auto &pair : G_Parameter[param])
        {
            double par_value = pair.first;
            G_Result_per_ParValue[param][par_value]->Fit("pol0", "QE");
            double a = G_Result_per_ParValue[param][par_value]->GetFunction("pol0")->GetParameter(0);
            double err_a = G_Result_per_ParValue[param][par_value]->GetFunction("pol0")->GetParError(0);

            if (G_Result_per_Par[param] == nullptr)
            {
                G_Result_per_Par[param] = new TGraphErrors();
                G_Result_per_Par[param]->SetName(("G_Result_" + param).c_str());
            }
            // adding data to graph for each param
            G_Result_per_Par[param]->AddPoint(par_value, a);
            G_Result_per_Par[param]->SetPointError(G_Result_per_Par[param]->GetN() - 1, 0, err_a);
        }
    }
}

void WriteHistograms(vector<string> parameters)
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
            G_Calibration_a[det]->Write();

            TCanvas *canvas_Eshift = new TCanvas(("Eshift_" + detectorName[det]).c_str(), ("Eshift_" + detectorName[det]).c_str(), 800, 600);
            G_Eshift[det]->SetTitle(("Eshift for " + detectorName[det]).c_str());
            G_Eshift[det]->GetYaxis()->SetTitle("Eshift [keV]");
            G_Eshift[det]->GetXaxis()->SetTitle("Parameter value");
            G_Eshift[det]->Write();
        }
    }

    for (const string &param : parameters)
    {
        for (auto &pair : G_Parameter[param])
        {
            TCanvas *canvas_param = new TCanvas(("Parameter_" + param + "_" + to_string(pair.first)).c_str(), ("Parameter_" + param + "_" + to_string(pair.first)).c_str(), 800, 600);
            TGraphErrors *graph = G_Result_per_ParValue[param][pair.first];
            if (graph != nullptr)
            {
                graph->SetTitle(("Parameter " + param + " for value " + to_string(pair.first)).c_str());
                graph->GetYaxis()->SetTitle("a");
                graph->GetXaxis()->SetTitle("Detector");
                graph->SetMarkerStyle(20);
                graph->SetMarkerSize(1);
                graph->Draw("AP");
                TLatex *latex = new TLatex(0.1, 0.9, ("Fitted a: " + to_string(graph->GetFunction("pol0")->GetParameter(0)) + " #pm " + to_string(graph->GetFunction("pol0")->GetParError(0))).c_str());
                latex->SetNDC();
                latex->SetTextSize(0.04);
                latex->Draw("SAME");               
                canvas_param->Write();
            }

            TCanvas *canvas_Eshift_param = new TCanvas(("Eshift_" + param + "_" + to_string(pair.first)).c_str(), ("Eshift_" + param + "_" + to_string(pair.first)).c_str(), 800, 600);
            if (G_Eshift_per_ParValue[param][pair.first] != nullptr)
            {
                G_Eshift_per_ParValue[param][pair.first]->SetTitle(("Eshift for parameter " + param + " with value " + to_string(pair.first)).c_str());
                G_Eshift_per_ParValue[param][pair.first]->GetYaxis()->SetTitle("Eshift [keV]");
                G_Eshift_per_ParValue[param][pair.first]->GetXaxis()->SetTitle("Detector");
                G_Eshift_per_ParValue[param][pair.first]->SetMarkerStyle(20);
                G_Eshift_per_ParValue[param][pair.first]->SetMarkerSize(1);
                G_Eshift_per_ParValue[param][pair.first]->Draw("AP");
                canvas_Eshift_param->Write();
            }

            TCanvas *canvas_DiffEshift_param = new TCanvas(("DiffEshift_" + param + "_" + to_string(pair.first)).c_str(), ("DiffEshift_" + param + "_" + to_string(pair.first)).c_str(), 800, 600);
            if (G_DiffEshift_per_ParValue[param][pair.first] != nullptr)
            {
                G_DiffEshift_per_ParValue[param][pair.first]->SetTitle(("Difference Eshift for parameter " + param + " with value " + to_string(pair.first)).c_str());
                G_DiffEshift_per_ParValue[param][pair.first]->GetYaxis()->SetTitle("Difference Eshift [keV]");
                G_DiffEshift_per_ParValue[param][pair.first]->GetXaxis()->SetTitle("Detector");
                G_DiffEshift_per_ParValue[param][pair.first]->SetMarkerStyle(20);
                G_DiffEshift_per_ParValue[param][pair.first]->SetMarkerSize(1);
                G_DiffEshift_per_ParValue[param][pair.first]->Draw("AP");
                canvas_DiffEshift_param->Write();
            }
        }

        TCanvas *canvas_result = new TCanvas(("Result_" + param).c_str(), ("Result_" + param).c_str(), 800, 600);
        TGraphErrors *result_graph = G_Result_per_Par[param];
        if (result_graph != nullptr)
        {
            result_graph->SetTitle(("Result for parameter " + param).c_str());
            result_graph->GetYaxis()->SetTitle("a");
            result_graph->GetXaxis()->SetTitle("Parameter value");
            result_graph->SetMarkerStyle(20);
            result_graph->SetMarkerSize(1);
            result_graph->Draw("AP");
            canvas_result->Write();
        }
    }    
}