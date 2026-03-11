#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include <vector>


using namespace std;


int MacroResolution_temp()
{


    TFile *fin = new TFile("/mnt/hgfs/shared-2/2025_DATA/MCP_DATA/CALIBRATED/run_001_StableBeamScan_secondstep_calibrated.root", "READ");

    vector<int> holes = {16, 17, 18, 23, 24, 25, 30, 31, 32};

    TGraphErrors *Res_X = new TGraphErrors();
    TGraphErrors *Res_Y = new TGraphErrors();

    TH1D *H_X = new TH1D("H_X", "H_X", 1000, 0, 0.3);
    TH1D *H_Y = new TH1D("H_Y", "H_Y", 1000, 0, 0.3);
    for (int hole : holes)
    {
        TCanvas *c = (TCanvas*)fin->Get(Form("Cell_fit_%.0f", (double)hole));

        for (int i = 0; i < c->GetListOfPrimitives()->GetEntries(); i++)
        {
            TObject *obj = c->GetListOfPrimitives()->At(i);
            if (obj->InheritsFrom("TPad") && string(obj->GetName()) == Form("Cell_fit_%.0f_3", (double)hole))
            {
                for (int j = 0; j < ((TPad*)obj)->GetListOfPrimitives()->GetEntries(); j++)
                {
                    TObject *obj2 = ((TPad*)obj)->GetListOfPrimitives()->At(j);
                    if (obj2->InheritsFrom("TF1"))
                    {
                        cout << "Hole " << hole << " X Res: " << ((TF1*)obj2)->GetParameter(2) << " +/- " << ((TF1*)obj2)->GetParError(2) << endl;
                        Res_X->AddPoint(Res_X->GetN(), ((TF1*)obj2)->GetParameter(2));
                        Res_X->SetPointError(Res_X->GetN()-1, 0, ((TF1*)obj2)->GetParError(2));
                        cout << "Res_X N: " << Res_X->GetN() << endl;

                        H_X->Fill(((TF1*)obj2)->GetParameter(2));
                    }
                }
            }

            if (obj->InheritsFrom("TPad") && string(obj->GetName()) == Form("Cell_fit_%.0f_4", (double)hole))
            {
                for (int j = 0; j < ((TPad*)obj)->GetListOfPrimitives()->GetEntries(); j++)
                {
                    TObject *obj2 = ((TPad*)obj)->GetListOfPrimitives()->At(j);
                    if (obj2->InheritsFrom("TF1"))
                    {
                        Res_Y->AddPoint(Res_Y->GetN(), ((TF1*)obj2)->GetParameter(2));
                        Res_Y->SetPointError(Res_Y->GetN()-1, 0, ((TF1*)obj2)->GetParError(2));

                        H_Y->Fill(((TF1*)obj2)->GetParameter(2));
                    }
                }
            }
        }
    }

    TFile *fout = new TFile("test.root", "RECREATE");
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    Res_X->SetTitle("MCP Resolution X;Hole Number;Resolution (mm)");
    Res_X->SetMarkerStyle(20);
    Res_X->Draw("AP");
    Res_X->Fit("pol0");
    c1->Write();

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    Res_Y->SetTitle("MCP Resolution Y;Hole Number;Resolution (mm)");
    Res_Y->SetMarkerStyle(20);
    Res_Y->Draw("AP");
    Res_Y->Fit("pol0");
    c2->Write();

    cout << "Mean X Resolution: " << H_X->GetMean() << " +/- " << H_X->GetMeanError() << endl;
    cout << "Mean Y Resolution: " << H_Y->GetMean() << " +/- " << H_Y->GetRMS()/sqrt(H_Y->GetEntries()) << endl;


    fout->Close();


    return 1;
}