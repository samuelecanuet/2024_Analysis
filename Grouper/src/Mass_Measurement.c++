#include "Mass_Measurement.hh"
#include "TROOT.h"

double ExtractParValue(const string &filename, const string &par_name)
{
    size_t pos = filename.find(par_name);
    if (pos == string::npos)
    {
        Error("Parameter name " + par_name + " not found in filename: " + filename);
        return 0;
    }

    size_t start = pos + par_name.length();
    size_t end = filename.find("_", start);
    if (end == string::npos)
    {
        end = filename.find(".", start);
        if (end == string::npos)
        {
            Error("Could not find end of parameter value in filename: " + filename);
            return 0;
        }
    }

    string value_str = filename.substr(start, end - start);
    try
    {
        return stod(value_str);
    }
    catch (const invalid_argument &e)
    {
        Error("Invalid parameter value in filename: " + filename);
        return 0;
    }
}

TGraphErrors *G_Eshift[SIGNAL_MAX];

int main(int argc, char **argv)
{
    
    NUCLEUS = "32Ar";
    FLAG2025 = true;
    VERBOSE = 0;


    InitDetectors("Config_Files/sample.pid");
    InitWindows();
        
    ///////// READING FILE //////////
    TFile *MASS_File = MyTFile((DIR_ROOT_DATA_ANALYSED + NUCLEUS + "_" + to_string(YEAR) + "_mass.root").c_str(), "RECREATE");

    TFile *ANALYSED_File_CENTRAL = MyTFile((DIR_ROOT_DATA_ANALYSED + NUCLEUS + "_" + to_string(YEAR) + "_analysed_new_N14.0_M3_T100.root").c_str(), "READ");
    TCanvas *c_CENTRAL = (TCanvas *)ANALYSED_File_CENTRAL->Get("Eshift/Eshift_M3/Eshift_M3");
    TGraphErrors *g_CENTRAL = (TGraphErrors *)c_CENTRAL->GetListOfPrimitives()->FindObject("");

    vector<string> filenames = SearchFilesIn("/run/media/local1/DATANEX/Samuel-G4/Systematics/a_Calibration/", "result");

    TFile *CALIBRATION_SHIFT = MyTFile("/run/media/local1/DATANEX/Samuel-G4/Systematics/DL/Systematics_DL_2025.root", "READ");

    // # plotting deltaM vs Eshift
    TGraphErrors *GMass = new TGraphErrors();

    for (auto filename : filenames)
    {
        TFile *ANALYSED_File_SIMULATED = MyTFile(("/run/media/local1/DATANEX/Samuel-G4/Systematics/a_Calibration/" + filename).c_str(), "READ");

        double par_value = ExtractParValue(filename, "a");
        
        // looping on all detectors
        for (int det = 1; det < SIGNAL_MAX; det++)
        {
            if (IsDetectorSiliStrip(det))// && G_Calibration_a[det] != nullptr)
            {
                // pair<double, double> Eshift = ComputeEshift(ANALYSED_File_SIMULATED, det);
                pair<double, double> Eshift;
                TCanvas *c_Eshift = (TCanvas *)ANALYSED_File_SIMULATED->Get("Eshift");
                TGraphErrors *G_Eshift_sim = (TGraphErrors*)c_Eshift->GetListOfPrimitives()->At(0);
                for (int i = 0; i < G_Eshift_sim->GetN(); i++)
                {
                    double x, y;
                    G_Eshift_sim->GetPoint(i, x, y);
                    if ((int)x == det)                    {
                        Eshift = make_pair(y, G_Eshift_sim->GetErrorY(i));
                        break;
                    }
                }
                if (G_Eshift[det] == nullptr)
                {
                    G_Eshift[det] = new TGraphErrors();
                    G_Eshift[det]->SetName(("G_Eshift_" + detectorName[det]).c_str());
                }
                // adding data to graph for each det
                G_Eshift[det]->AddPoint(Eshift.first, par_value);
                G_Eshift[det]->SetPointError(G_Eshift[det]->GetN() - 1, Eshift.second, 0.00001);
            }
        }
    }

    MASS_File->cd();
    for (int det = 1; det < SIGNAL_MAX; det++)
    {
        if (G_Eshift[det] != nullptr)
        {
            TCanvas *canvas = new TCanvas(("Eshift_vs_Mass_" + detectorName[det]).c_str(), ("Eshift_vs_Mass_" + detectorName[det]).c_str(), 800, 600);

            TF1 *fit = new TF1("fit", "[0]*x + [1]", -1e6, 1e6);            
            G_Eshift[det]->Fit(fit, "SE");

            double Eshiftexp = 0;
            double Eshiftexp_err = 0;
            for (int i = 0; i < g_CENTRAL->GetN(); i++)
            {
                double x, y;
                g_CENTRAL->GetPoint(i, x, y);
                if (x == det)
                {
                    Eshiftexp = y;
                    Eshiftexp_err = g_CENTRAL->GetErrorY(i);
                    break;
                }
            }

            G_Eshift[det]->SetTitle(("Eshift vs Mass for " + detectorName[det]).c_str());
            G_Eshift[det]->GetYaxis()->SetTitle("Mass [keV]");
            G_Eshift[det]->GetXaxis()->SetTitle("Eshift [keV]");
            G_Eshift[det]->SetMarkerStyle(20);
            G_Eshift[det]->SetMarkerSize(1);
            G_Eshift[det]->Draw("AP");

            TLine *line = new TLine(Eshiftexp, canvas->GetUymin(), Eshiftexp, canvas->GetUymax());
            line->SetLineColor(kRed);
            line->SetLineStyle(2);
            line->Draw("SAME");

            fit->Draw("SAME");
            
            canvas->Write();

            GMass->AddPoint(det, G_Eshift[det]->GetFunction("fit")->Eval(Eshiftexp));
            GMass->SetPointError(GMass->GetN() - 1, 0, sqrt(pow(G_Eshift[det]->GetFunction("fit")->Eval(Eshiftexp + Eshiftexp_err) - G_Eshift[det]->GetFunction("fit")->Eval(Eshiftexp), 2) + 
        pow(G_Eshift[det]->GetFunction("fit")->GetParError(0) * Eshiftexp, 2)));
        }
    }

    TCanvas *canvas_result = new TCanvas("Mass_vs_Eshift", "Mass_vs_Eshift", 800, 600);
    GMass->SetTitle("Mass vs Eshift");
    GMass->GetYaxis()->SetTitle("Mass [keV]");
    GMass->GetXaxis()->SetTitle("Detector");
    GMass->SetMarkerStyle(20);
    GMass->SetMarkerSize(1);
    GMass->Draw("AP");
    GMass->Fit("pol0", "SR", "", 1, 90);
    cout << "Mass fit result: " << GMass->GetFunction("pol0")->GetParameter(0) << " +/- " << GMass->GetFunction("pol0")->GetParError(0) << endl;
    GMass->Draw("AP");
    canvas_result->Write();
    MASS_File->Close();

    return 0;
}