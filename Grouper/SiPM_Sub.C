#include "include/Detectors.hh"


void SiPM_Sub()
{
    TFile *f137 = new TFile((DIR_ROOT_DATA_CALIBRATED + "SiPM_Matching_137.root").c_str(), "READ");
    TFile *f134 = new TFile((DIR_ROOT_DATA_CALIBRATED + "SiPM_Matching_134.root").c_str(), "READ");

    double time137 = 15*60+12;
    double time134 = 58*60+13;


    TH1D* H137[9];
    TH1D* H134[9];

    TH1D* Sub[9];

    for (int det = 1; det <= 9; det++)
    {
        H137[det] = (TH1D *)f137->Get(("207Bi/SiPM_"+to_string(det)+"/H_SiPM_High_207Bi_BetaHi"+to_string(det)).c_str());
        H134[det] = (TH1D *)f134->Get(("207Bi/SiPM_"+to_string(det)+"/H_SiPM_High_207Bi_BetaHi"+to_string(det)).c_str());

        Sub[det] = (TH1D *)H137[det]->Clone();
        Sub[det]->Add(H134[det], -time137/time134);

        H134[det]->Scale(time137/time134);
    }

    TCanvas *c = new TCanvas("c", "c", 800, 800);
    c->Divide(3, 3);
    for (int det = 1; det <= 9; det++)
    {
        c->cd(det);
        H137[det]->SetLineColor(kBlack);
        H137[det]->Draw("HIST");
        H134[det]->SetLineColor(kRed);
        H134[det]->Draw("SAME");
        Sub[det]->Draw("SAME");
    }
    c->Draw();

    // f137->Close();  
    // f134->Close();
}