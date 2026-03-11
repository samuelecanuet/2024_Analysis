#include "Grouper/include/Detectors.hh"
#include "Grouper/include/Utilities.hh"


int EshiftCompareSimExp()
{
    FLAG2025 = true;
    InitDetectors("Grouper/Config_Files/sample.pid");
    // TFile *Sim = MyTFile((DIR_ROOT_DATA_SIMULATED + "Result/" + "32Ar_" + to_string(YEAR) + "_result.root").c_str(), "READ");
    TFile *Sim = MyTFile(("/run/media/local1/DATANEX/Samuel-G4/32Ar_IAS_" + to_string(YEAR) + "_Default_result.root").c_str(), "READ");
    // TFile *Exp = MyTFile((DIR_ROOT_DATA_SIMULATED + "Result/" + "32Ar_" + to_string(YEAR) + "_result_5deg.root").c_str(), "READ");
    TFile *Exp = MyTFile((DIR_ROOT_DATA_ANALYSED + "32Ar_" + to_string(YEAR) + "_analysed_all.root").c_str(), "READ");

    TCanvas *cEshiftCompare = new TCanvas("cEshiftCompare", "cEshiftCompare", 1200, 800);

    //exp
    TCanvas *cEshiftCompare_Corrected = (TCanvas*)Exp->Get("Eshift/Eshift_M1/Eshift_1");
    TGraphErrors *G_Exp = nullptr;
    for (auto key: *cEshiftCompare_Corrected->GetListOfPrimitives()) {
        if (string(key->ClassName()).find("TGraphErrors") != string::npos) {
            G_Exp = (TGraphErrors*)key;
        }
    } 
    // TGraphErrors *G_Exp = nullptr;
    // TCanvas *cEshiftCompare_Corrected = (TCanvas*)Exp->Get("Eshift");
    // for (auto key: *cEshiftCompare_Corrected->GetListOfPrimitives()) {
    //     if (string(key->ClassName()).find("TGraphErrors") != string::npos) {
    //         G_Exp = (TGraphErrors*)key;
    //     }
    // }

    //sim
    TGraphErrors *G_Sim = nullptr;
    TCanvas *c_Sim = (TCanvas*)Sim->Get("Eshift");
    for (auto key: *c_Sim->GetListOfPrimitives()) {
        if (string(key->ClassName()).find("TGraphErrors") != string::npos) {
            G_Sim = (TGraphErrors*)key;
        }
    }


    cEshiftCompare->cd();
    G_Exp->SetMarkerColor(kBlack);
    G_Exp->SetLineColor(kBlack);
    G_Exp->SetMarkerStyle(20);
    G_Exp->SetTitle("Eshift Comparison Simulated vs Experimental;Code;Eshift (keV)");
    G_Exp->Draw("AP");
    G_Sim->SetMarkerColor(kRed);
    G_Sim->SetLineColor(kRed);
    G_Sim->SetMarkerStyle(21);
    G_Sim->Draw("P SAME");
    TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
    leg->AddEntry(G_Exp, "Experimental", "lp");
    leg->AddEntry(G_Sim, "Simulated", "lp");
    leg->Draw("SAME");
    cEshiftCompare->Draw();



    return 0;
}