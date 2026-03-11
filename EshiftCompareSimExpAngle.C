#include "Grouper/include/Detectors.hh"
#include "Grouper/include/Utilities.hh"


int EshiftCompareSimExpAngle()
{
    FLAG2025 = true;
    InitDetectors("Grouper/Config_Files/sample.pid");

    map<int, TGraphErrors*> SimGraphs;

    TCanvas *cEshiftCompare = new TCanvas("cEshiftCompare", "cEshiftCompare", 1200, 800);
    TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);


    for (int i = -4; i <= 5; i++) 


    {

        // if (i != 0 && i != -4 && i != 1) 
        //     continue;

        string name = Form("32Ar_%d_result_%ddeg.root", YEAR, i);
        TFile *Sim = MyTFile((DIR_ROOT_DATA_SIMULATED + "Result/" + name.c_str()), "READ");

        //sim
        SimGraphs[i] = nullptr;
        TCanvas *c_Sim = (TCanvas*)Sim->Get("Eshift");
        for (auto key: *c_Sim->GetListOfPrimitives()) {
            if (string(key->ClassName()).find("TGraphErrors") != string::npos) {
                SimGraphs[i] = (TGraphErrors*)key;
            }
        }


        cEshiftCompare->cd();
        SimGraphs[i]->SetName(Form("G_Sim_%ddeg", i));
        SimGraphs[i]->SetTitle(Form("Eshift Simulated %d deg", i));
        SimGraphs[i]->SetMarkerColor(i + 5 == 10 ? kOrange : i + 6);
        SimGraphs[i]->SetLineColor(i + 5 == 10 ? kOrange : i + 6);
        SimGraphs[i]->SetMarkerStyle(21);
        if (i == -4)
            SimGraphs[i]->Draw("AP");
        else 
            SimGraphs[i]->Draw("P SAME");   
        leg->AddEntry(SimGraphs[i], Form("Simulated %d deg", i), "lp");
    }

    TFile *Exp = MyTFile((DIR_ROOT_DATA_ANALYSED + "32Ar_" + to_string(YEAR) + "_analysed_1.root").c_str(), "READ");

    //exp
    TCanvas *cEshiftCompare_Corrected = (TCanvas*)Exp->Get("Eshift/Eshift_M1/Eshift_1");
    TGraphErrors *G_Exp = nullptr;
    for (auto key: *cEshiftCompare_Corrected->GetListOfPrimitives()) {
        if (string(key->ClassName()).find("TGraphErrors") != string::npos) {
            G_Exp = (TGraphErrors*)key;
        }
    } 

    cEshiftCompare->cd();
    G_Exp->SetMarkerColor(kBlack);
    G_Exp->SetLineColor(kBlack);
    G_Exp->SetMarkerStyle(20);
    G_Exp->SetMarkerSize(2);
    G_Exp->SetTitle("Eshift Comparison Simulated vs Experimental;Code;Eshift (keV)");
    G_Exp->Draw("P SAME");
    // G_Sim->SetMarkerColor(kRed);
    // G_Sim->SetLineColor(kRed);
    // G_Sim->SetMarkerStyle(21);
    // G_Sim->Draw("P SAME");
    // TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
    leg->AddEntry(G_Exp, "Experimental", "lp");
    // leg->AddEntry(G_Sim, "Simulated", "lp");
    leg->Draw("SAME");
    cEshiftCompare->Draw();



    return 0;
}