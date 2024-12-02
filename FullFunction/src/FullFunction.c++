#include "FullFunction.hh"




int main()
{

    // Input parameters
    double A = 32.0;
    double Z = 18;
    double E_p = 3356.2;                       // keV
    double Gamma = 20e-3;                       // keV
    double W0 = sqrt(5546 * 5546 + 511 * 511); // keV
    double a = 1;

    double Phi_max = 0;
    double Phi_min = M_PI;

    double M = A * m;

    // Generate T values
    const int N = 2000;
    double T[N];
    double Wt[N];
    double k = K(E_p, M);

    double par[9] = {2700, E_p, W0, k, a, A, Z, Phi_min, Phi_max};

    for (int i = 0; i < N; ++i)
    {
        T[i] = -20 + i * 40 / N + E_p;
        Wt[i] = fTotalWeightKinematicShift(&T[i], par);
    }


    TF1 *KinematicShift = new TF1("KinematicShift", fTotalWeightKinematicShift, 3300, 3400, 9);
    KinematicShift->SetParameter(0, 1e2);
    KinematicShift->SetParLimits(0, 0, 1e5);
    KinematicShift->FixParameter(1, E_p);
    KinematicShift->FixParameter(2, W0);
    KinematicShift->FixParameter(3, k);
    KinematicShift->FixParameter(4, a);
    KinematicShift->FixParameter(5, A);
    KinematicShift->FixParameter(6, Z);
    KinematicShift->FixParameter(7, Phi_min);
    KinematicShift->FixParameter(8, Phi_max);

    TF1 *NuclearBroadering = new TF1("NuclearBroadering", fNuclearBroadering, 3300, 3400, 2);
    NuclearBroadering->FixParameter(0, E_p);
    NuclearBroadering->FixParameter(1, Gamma);

    TF1Convolution *conv = new TF1Convolution("KinematicShift", "NuclearBroadering", 3300, 3400, true);
    conv->SetRange(3300, 3400);
    conv->SetNofPointsFFT(1000);

    TF1 *convolution = new TF1("convolution", *conv, 3300, 3400, conv->GetNpar());
    convolution->SetParameter(0, 1e2);
    convolution->SetParLimits(0, 0, 1e5);
    convolution->FixParameter(1, E_p);
    convolution->FixParameter(2, W0);
    convolution->FixParameter(3, k);
    convolution->FixParameter(4, a);
    convolution->FixParameter(5, A);
    convolution->FixParameter(6, Z);
    convolution->FixParameter(7, Phi_min);
    convolution->FixParameter(8, Phi_max);
    convolution->FixParameter(9, Gamma);
    convolution->FixParameter(10, E_p);
  


    // Plotting with ROOT
    TFile *file = TFile::Open("../../../../../../../mnt/hgfs/shared-2/2024_DATA/SIMULATED_DATA/fe/fe/test_analysed.root", "READ");
    TH1 *h = dynamic_cast<TH1 *>(file->Get("p/H_E0_2212"));

    TCanvas *c1 = new TCanvas("c1", "Recoil Broadening", 800, 600);
    h->GetXaxis()->SetRangeUser(3300, 3400);
    h->Draw();

    TGraph *graph = new TGraph(N);
    for (int i = 0; i < N; ++i)
    {
        graph->SetPoint(i, T[i] , Wt[i]);
    }
    graph->SetLineColor(kRed);
    graph->Draw("L SAME");

    c1->SaveAs("recoil_broadening.png");

    TCanvas *c2 = new TCanvas("c2", "Recoil Broadening", 800, 600);
    h->Draw();
    h->Fit("KinematicShift");
    KinematicShift->SetLineColor(kRed);
    KinematicShift->Draw("SAME");

    c2->SaveAs("recoil_broadening_fit.png");

    TCanvas *c3 = new TCanvas("c3", "Recoil Broadening", 800, 600);
    h->Draw();
    h->Fit("convolution");
    convolution->SetLineColor(kRed);
    convolution->Draw("SAME");

    c3->SaveAs("recoil_broadening_convolution.png");


    return 0;
}
