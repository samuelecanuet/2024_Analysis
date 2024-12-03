#include "FullFunction.hh"
#include "Math/Minimizer.h"
#include <random>

int main()
{

    default_random_engine gen;

    // Input parameters
    double A = 32.0;
    double Z = 18;
    double E_p = 3356.6;                       // keV
    double Gamma = 20e-3;                       // keV
    double Q = 5046.6;                          // keV
    double W0 = sqrt((Q+511) * (Q+511) + 511 * 511); // keV
    double a = 1;

    double Phi_max = 0;
    double Phi_min = M_PI;

    double M = A * m;

    double E_min = 3300;
    double E_max = 3400;


    double k = K(E_p, M);


    Peak *IAS = new Peak(E_p, 1, Gamma, Q, A, Z, a, Phi_max, Phi_min, E_min, E_max);

    MyTF1 * convolution1 = IAS->GetFunction();


    // TF1 *KinematicShift = new TF1("KinematicShift", fTotalWeightKinematicShift, E_min, E_max, 9);
    // KinematicShift->FixParameter(0, 1);
    // KinematicShift->FixParameter(1, E_p);
    // KinematicShift->FixParameter(2, W0);
    // KinematicShift->FixParameter(3, k);
    // KinematicShift->FixParameter(4, a);
    // KinematicShift->FixParameter(5, A);
    // KinematicShift->FixParameter(6, Z);
    // KinematicShift->FixParameter(7, Phi_min);
    // KinematicShift->FixParameter(8, Phi_max);

    // TF1 *NuclearBroadering = new TF1("NuclearBroadering", fNuclearBroadering, E_min, E_max, 2);
    // NuclearBroadering->SetParameter(0, 1e5);
    // NuclearBroadering->SetParLimits(0, 1e2, 1e7);
    // NuclearBroadering->FixParameter(1, Gamma);


    // // Create a Convolution object
    // Convolution conv1(KinematicShift, NuclearBroadering, E_min, E_max);

    // Define the convolution function
    // TF1 *convolution1 = new TF1("convolution1", [&conv1](double *x, double *p) {
    //     return conv1.Evaluate(x, p);
    // }, E_min, E_max, 11);
    // convolution1->FixParameter(0, 1);
    // convolution1->FixParameter(1, E_p);
    // convolution1->FixParameter(2, W0);
    // convolution1->FixParameter(3, k);
    // convolution1->FixParameter(4, a);
    // convolution1->FixParameter(5, A);
    // convolution1->FixParameter(6, Z);
    // convolution1->FixParameter(7, Phi_min);
    // convolution1->FixParameter(8, Phi_max);
    // convolution1->SetParameter(9, 1e5);
    // convolution1->SetParLimits(9, 1e2, 1e7);
    // convolution1->FixParameter(10, Gamma);

    TF1* Gaussian = new TF1("Gaussian", fGauss, E_min, E_max, 3);
    Gaussian->SetParameter(0, 0);
    Gaussian->SetParameter(1, 5);

    Convolution conv2(convolution1, Gaussian, E_min, E_max);

    // Define the convolution function
    TF1 *convolution2 = new TF1("convolution2", [&conv2](double *x, double *p) {
        return conv2.Evaluate(x, p);
    }, E_min, E_max, 13);
    convolution2->FixParameter(0, 1);
    convolution2->FixParameter(1, E_p);
    convolution2->FixParameter(2, W0);
    convolution2->FixParameter(3, k);
    convolution2->FixParameter(4, a);
    convolution2->FixParameter(5, A);
    convolution2->FixParameter(6, Z);
    convolution2->FixParameter(7, Phi_min);
    convolution2->FixParameter(8, Phi_max);
    convolution2->SetParameter(9, 1e5);
    convolution2->SetParLimits(9, 0, 1e7);
    convolution2->SetParameter(10, Gamma);
    convolution2->SetParLimits(10, 0, 10);  
    convolution2->FixParameter(11, 0);
    convolution2->SetParameter(12, 2);
    convolution2->SetParLimits(12, 0, 10);
    



    // Plotting with ROOT
    TFile *file = TFile::Open("../../../../../../../mnt/hgfs/shared-2/2024_DATA/SIMULATED_DATA/fe/fe/test_analysed.root", "READ");
    TH1 *hh = dynamic_cast<TH1 *>(file->Get("p/H_E0_2212"));
    
    TH1D* h = (TH1D*)hh->Clone("h");
    h->Reset();

    std::cauchy_distribution<double> distribution(0, 0.5/2);
    std::normal_distribution<double> distribution2(0, 2);
    for (int i = 0; i < hh->GetEntries(); ++i)
    {
        h->Fill(hh->GetRandom() + distribution(gen) + distribution2(gen));
    }

    TCanvas *c1 = new TCanvas("c1", "Recoil Broadening", 800, 600);
    h->GetXaxis()->SetRangeUser(E_min, E_max);
    h->Draw();

    // TCanvas *c2 = new TCanvas("c2", "Recoil Broadening", 800, 600);
    // h->Draw();
    // h->Fit("KinematicShift", "N");
    // KinematicShift->SetLineColor(kRed);
    // KinematicShift->Draw("SAME");

    // c2->SaveAs("recoil_broadening_fit.png");

    TCanvas *c3 = new TCanvas("c3", "Recoil Broadening", 800, 600);
    h->Draw();
    clock_t start = clock();
    h->Fit("convolution2", "MULTITHREAD RNEL", "", E_min, E_max);
    convolution1->SetLineColor(kRed);
    convolution1->Draw("SAME");
    cout << "Time: " << (clock() - start) / (double) CLOCKS_PER_SEC / 8 << endl;

    c3->SetLogy();
    c3->SaveAs("recoil_broadening_convolution.png");


    return 0;
}
