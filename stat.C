#include "Grouper/include/Detectors.hh"
#include "Grouper/include/Utilities.hh"

int stat()
{
    FLAG2025 = true;
    InitDetectors("Grouper/Config_Files/sample.pid");

    TFile *Exp = MyTFile((DIR_ROOT_DATA_ANALYSED + "32Ar_" + to_string(YEAR) + "_analysed_all.root").c_str(), "READ");    
    TCanvas *cEshiftCompare_Corrected = (TCanvas *)Exp->Get("Eshift/Eshift_M3/Eshift_3");
    TGraphErrors *G_Exp = nullptr;
    for (auto key : *cEshiftCompare_Corrected->GetListOfPrimitives())
    {
        if (string(key->ClassName()).find("TGraphErrors") != string::npos)
        {
            G_Exp = (TGraphErrors *)key;
        }
    }

    double sum = 0.0;

    for (int i = 0; i < G_Exp->GetN(); i++)
    {
        double x_err, y_err, x, y;
        G_Exp->GetPoint(i, x, y);
        y_err = G_Exp->GetErrorY(i);
        
        sum += y_err;
    }

    cout << "Average error: " << 3./4. * sum / G_Exp->GetN() / sqrt(G_Exp->GetN()) << endl;


    return 0;
}