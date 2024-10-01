
#include "TEST.hh"

int main()
{

    MERGED_File = new TFile("../../../../../mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/CALIBRATED/Calibrated.root", "READ");
    SIMULATED_File = new TFile("../../../../../mnt/hgfs/shared-2/2024_DATA/SIMULATED_DATA/32ArRMATRIX__CS0_CSP0_CV1_CVP1.root", "READ");
    TEST_File = new TFile("test.root", "RECREATE");

    H_Exp = (TH1D*)MERGED_File->Get("D1.1/32Ar/H_Exp_32Ar_D1.1");
    H_Sim = (TH1D*)SIMULATED_File->Get("D1.1_single");
    H_Sim_Conv = (TH1D*)H_Sim->Clone("H_Sim_Conv");
    H_Sim_Conv->Reset();
    H_Sim_Conv2 = (TH1D*)H_Sim->Clone("H_Sim_Conv2");
    H_Sim_Conv2->Reset();

    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
    TCanvas *c2 = new TCanvas("c2", "c2", 1920, 1080);
    fBC = new TF1("fBC", gaussBortelsCollaers, 0, 10000, 6);
    fBC->SetNpx(10000);
    H_Exp->GetXaxis()->SetRangeUser(3200, 3500);
    fBC->SetParameter(0, 2e6);
    fBC->SetParLimits(1, 3200, 3500);
    fBC->SetParameter(1, 3352);
    fBC->SetParLimits(2, 0, 100);
    fBC->SetParameter(2, 76);
    fBC->SetParLimits(3, 0, 1);
    fBC->SetParameter(3, 0.002);
    fBC->SetParLimits(4, 0, 100);
    fBC->SetParameter(4, 3);
    fBC->SetParLimits(5, 0, 100);
    fBC->SetParameter(5, 30);
    fBC->SetLineColor(kGreen);
    H_Exp->Fit("fBC", "RN");
    H_Exp->Draw("HIST");
    fBC->Draw("SAME");
    c2->Draw();

    H_Sim_Conv->Rebin(10);
    H_Sim_Conv2->Rebin(10);
    H_Sim->Rebin(10);
    fBC->SetParameter(2, 0.1);
    fBC->SetParameter(3, 0.008);
    fBC->SetParameter(4, 4);
    fBC->SetParameter(5, 60);

    Chi2Mini();


    TCanvas *c = new TCanvas("c", "c", 1920, 1080);
    c->SetLogy();
    // H_Sim_Conv->Rebin(10);
    H_Exp->GetXaxis()->SetRangeUser(3300, 3400);
    H_Sim_Conv->GetXaxis()->SetRangeUser(3300, 3400);
    H_Sim->GetXaxis()->SetRangeUser(3300, 3400);
    H_Sim_Conv->Scale(H_Exp->Integral() / H_Sim_Conv->Integral());
    H_Sim->Scale(H_Exp->Integral() / H_Sim->Integral());
    


    H_Exp->GetXaxis()->SetRangeUser(3200, 3500);
    H_Sim_Conv->GetXaxis()->SetRangeUser(3200, 3500);
    H_Sim_Conv->SetLineColor(kRed);
    H_Exp->Draw("HIST");
    H_Sim->SetLineColor(kBlue);
    H_Sim->Draw("HIST SAME");
    H_Sim_Conv->Draw("HIST SAME");
    c->Write();

    

    H_Sim_Conv->GetXaxis()->SetRangeUser(3275, 3375);
    H_Sim->GetXaxis()->SetRangeUser(3275, 3375);

    cout << H_Sim_Conv->GetMean() << " " << H_Sim_Conv->GetRMS() << endl;
    cout << H_Sim->GetMean() << " " << H_Sim->GetRMS() << endl;

    TEST_File->Close();

    // cout << "Chi2: " << H_Exp->Chi2Test(H_Sim_Conv, "UW CHI2/NDF") << endl;


}

