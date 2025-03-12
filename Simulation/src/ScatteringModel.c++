#include "ScatteringModel.hh"

int main()
{

    TFile *fModel = MyTFile(DIR_ROOT_DATA_SIMULATED + "models.root", "RECREATE");

    TFile *fopt4 = MyTFile(DIR_ROOT_DATA_SIMULATED + "24-02/32Ar_opt4_a1_b0_analysed.root", "READ");
    TFile *fliv = MyTFile(DIR_ROOT_DATA_SIMULATED + "24-02/32Ar_liv_a1_b0_analysed.root", "READ");
    TFile *fgs = MyTFile(DIR_ROOT_DATA_SIMULATED + "24-02/32Ar_gs_a1_b0_analysed.root", "READ");
    TFile *fpen = MyTFile(DIR_ROOT_DATA_SIMULATED + "24-02/32Ar_pen_a1_b0_analysed.root", "READ");
    TFile *fwvi = MyTFile(DIR_ROOT_DATA_SIMULATED + "24-02/32Ar_wvi_a1_b0_analysed.root", "READ");

    TH1D *hopt4 = (TH1D *)fopt4->Get("e+/PlasticScintillator/PlasticScintillator_Energy_Deposit_e+");
    TH1D *hliv = (TH1D *)fliv->Get("e+/PlasticScintillator/PlasticScintillator_Energy_Deposit_e+");
    TH1D *hgs = (TH1D *)fgs->Get("e+/PlasticScintillator/PlasticScintillator_Energy_Deposit_e+");
    TH1D *hpen = (TH1D *)fpen->Get("e+/PlasticScintillator/PlasticScintillator_Energy_Deposit_e+");
    TH1D *hwvi = (TH1D *)fwvi->Get("e+/PlasticScintillator/PlasticScintillator_Energy_Deposit_e+");

    fModel->cd();
    

    TCanvas *call = new TCanvas("Scattering Model", "Scattering Model", 800, 600);
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    hopt4->SetLineColor(kRed);
    leg->AddEntry(hopt4, "Opt4", "l");
    hopt4->Draw("HIST");
    hliv->SetLineColor(kBlue);
    leg->AddEntry(hliv, "Livermore", "l");
    hliv->Draw("HIST SAME");
    hgs->SetLineColor(kGreen);
    leg->AddEntry(hgs, "GS", "l");
    hgs->Draw("HIST SAME");
    leg->Draw("SAME");
    hpen->SetLineColor(kBlack);
    leg->AddEntry(hpen, "Penelope", "l");
    hpen->Draw("HIST SAME");
    hwvi->SetLineColor(kMagenta);
    leg->AddEntry(hwvi, "WVI", "l");
    hwvi->Draw("HIST SAME");
    call->Write();

    cout << "Writing" << endl;

    // OPT4 vs GS
    TCanvas *cgs = new TCanvas("Scattering Model GS", "Scattering Model GS", 800, 600);
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
    pad1->Draw();
    pad1->cd();
    hopt4->SetLineColor(kRed);
    hopt4->Draw("HIST");
    hgs->SetLineColor(kGreen);
    hgs->Draw("HIST SAME");
    leg->Draw("SAME");
    cgs->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
    pad2->Draw();
    pad2->cd();
    TH1D *hgs_ratio = (TH1D *)hopt4->Clone("hgs_ratio");
    hgs_ratio->Add(hgs, -1);
    hgs_ratio->Divide(hopt4);
    hgs_ratio->Draw();
    cgs->Write();

    // OPT4 vs Livermore
    TCanvas *cliv = new TCanvas("Scattering Model Livermore", "Scattering Model Livermore", 800, 600);
    TPad *pad3 = new TPad("pad3", "pad3", 0, 0.3, 1, 1);
    pad3->Draw();
    pad3->cd();
    hopt4->SetLineColor(kRed);
    hopt4->Draw("HIST");
    hliv->SetLineColor(kBlue);
    hliv->Draw("HIST SAME");
    leg->Draw("SAME");
    cliv->cd();
    TPad *pad4 = new TPad("pad4", "pad4", 0, 0, 1, 0.3);
    pad4->Draw();
    pad4->cd();
    TH1D *hliv_ratio = (TH1D *)hopt4->Clone("hliv_ratio");
    hliv_ratio->Add(hliv, -1);
    hliv_ratio->Divide(hopt4);
    hliv_ratio->Draw();
    cliv->Write();

    // OPT4 vs Penelope
    TCanvas *cpen = new TCanvas("Scattering Model Penelope", "Scattering Model Penelope", 800, 600);
    TPad *pad5 = new TPad("pad5", "pad5", 0, 0.3, 1, 1);
    pad5->Draw();
    pad5->cd();
    hopt4->SetLineColor(kRed);
    hopt4->Draw("HIST");
    hpen->SetLineColor(kBlack);
    hpen->Draw("HIST SAME");
    leg->Draw("SAME");
    cpen->cd();
    TPad *pad6 = new TPad("pad6", "pad6", 0, 0, 1, 0.3);
    pad6->Draw();
    pad6->cd();
    TH1D *hpen_ratio = (TH1D *)hopt4->Clone("hpen_ratio");
    hpen_ratio->Add(hpen, -1);
    hpen_ratio->Divide(hopt4);
    hpen_ratio->Draw();
    cpen->Write();

    // OPT4 vs WVI
    TCanvas *cwvi = new TCanvas("Scattering Model WVI", "Scattering Model WVI", 800, 600);
    TPad *pad7 = new TPad("pad7", "pad7", 0, 0.3, 1, 1);
    pad7->Draw();
    pad7->cd();
    hopt4->SetLineColor(kRed);
    hopt4->Draw("HIST");
    hwvi->SetLineColor(kMagenta);
    hwvi->Draw("HIST SAME");
    leg->Draw("SAME");
    cwvi->cd();
    TPad *pad8 = new TPad("pad8", "pad8", 0, 0, 1, 0.3);
    pad8->Draw();
    pad8->cd();
    TH1D *hwvi_ratio = (TH1D *)hopt4->Clone("hwvi_ratio");
    hwvi_ratio->Add(hwvi, -1);
    hwvi_ratio->Divide(hopt4);
    hwvi_ratio->Draw();
    cwvi->Write();
    
    fModel->Close();

    return 0;
}