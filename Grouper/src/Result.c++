#include "Result.hh"
#include "TROOT.h"
int main(int argc, char **argv)
{
    
    NUCLEUS = "32Ar";
    FLAG2025 = true;
    VERBOSE = 0;


    InitDetectors("Config_Files/sample.pid");
        
    ///////// READING FILE //////////
    TFile *ANALYSED_File = MyTFile((DIR_ROOT_DATA_RESULTS + NUCLEUS + "_" + to_string(YEAR) + "_result.root").c_str(), "RECREATE");

    TFile *ANALYSED_File_CENTRAL = MyTFile((DIR_ROOT_DATA_ANALYSED + NUCLEUS + "_" + to_string(YEAR) + "_analysed.root").c_str(), "READ");
    TCanvas *c_CENTRAL = (TCanvas *)ANALYSED_File_CENTRAL->Get("Eshift/Eshift_M3/Eshift_3");
    TGraphErrors *g_CENTRAL = (TGraphErrors *)c_CENTRAL->GetListOfPrimitives()->FindObject("");

    TFile *ANALYSED_File_p1 = MyTFile((DIR_ROOT_DATA_ANALYSED + NUCLEUS + "_" + to_string(YEAR) + "_analysed_p1.root").c_str(), "READ");
    TCanvas *c_p1 = (TCanvas *)ANALYSED_File_p1->Get("Eshift/Eshift_M3/Eshift_3");
    TGraphErrors *g_p1 = (TGraphErrors *)c_p1->GetListOfPrimitives()->FindObject("");

    TFile *CALIBRATION_SHIFT = MyTFile("/run/media/local1/DATANEX/Samuel-G4/Systematics/DL/Systematics_DL_2025.root", "READ");

    TGraphErrors *G = new TGraphErrors();

    vector<double> values;
    for (int i = 0; i < g_CENTRAL->GetN(); i++)
    {
        double x, y, ex, ey;
        g_CENTRAL->GetPoint(i, x, y);
        ex = g_CENTRAL->GetErrorX(i);
        ey = g_CENTRAL->GetErrorY(i);

        double xp, yp, exp, eyp;
        g_p1->GetPoint(i, xp, yp);
        exp = g_p1->GetErrorX(i);
        eyp = g_p1->GetErrorY(i);
        gROOT->SetBatch(kTRUE);
        // Info( + ": " + to_string(y-yp) , 1);

        // calib
        TCanvas *c_calib = (TCanvas *)CALIBRATION_SHIFT->Get(("Calibration_a_" + detectorName[(int)x]).c_str());
        TGraphErrors *g_calib = (TGraphErrors *)c_calib->GetListOfPrimitives()->FindObject(Form("G_Calibration_a_%s", detectorName[(int)x].c_str()));
        TF1 *f_calib = g_calib->GetFunction("pol1");

        Info(detectorName[(int)x] + "   delata: " + to_string(f_calib->Eval(y)-f_calib->Eval(yp)), 1);

        values.push_back(f_calib->Eval(y)-f_calib->Eval(yp));

        G->AddPoint(x, f_calib->Eval(y)-f_calib->Eval(yp));
        G->SetPointError(G->GetN() - 1, 0, abs(f_calib->Eval(y)-f_calib->Eval(y+ey)));
    }

    cout << "Average: " << accumulate(values.begin(), values.end(), 0.0) / values.size() << endl;

    ANALYSED_File->cd();
    G->Fit("pol0");
    cout << "Fit result: " << G->GetFunction("pol0")->GetParameter(0) << " +/- " << G->GetFunction("pol0")->GetParError(0) << endl;
    G->Write("AP");
    ANALYSED_File->Close();

    return 0;
}