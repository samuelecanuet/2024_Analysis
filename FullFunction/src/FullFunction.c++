#include "FullFunction.hh"
#include "Math/Minimizer.h"
#include <random>

int main()
{
    

    default_random_engine gen;

    // Spectrum *spectrum = new Spectrum(18, 33, "proton");
    // spectrum->GenerateSpectrum();
    // TF1* fTotal = spectrum->GetFunction();
    // vector <TF1*> fPeaks = spectrum->GetTF1Peaks();
    // vector <Peak*> Peaks = spectrum->GetPeaks();

    // Input parameters
    double step = 0.1; //keV
    double A = 33.;
    double Z = 18.;
    double Phi_max = 0;
    double Phi_min = M_PI;
    

    // ////////////// SPECTRUM //////////////
    double E_p = 3356.6;                       
    double Gamma = 20e-3;                       
    double Q = 5046.6;                          
    double W0 = sqrt((Q+511) * (Q+511) + 511 * 511);
    double a = 1;
    double E_min = (int)E_p - 50;
    double E_max = (int)E_p + 50;

    Peak *IAS = new Peak(E_p, 1, Gamma, W0, A, Z, a, Phi_max, Phi_min, E_min, E_max, step);
    TF1 *fIAS = (TF1*)IAS->GetFunction();

    // GT1
    E_p = 2190.6*31/32;
    Gamma = 3;//SecondTokeV(1.52e-19);
    Q = 7362.32-1022;
    W0 = sqrt((Q+511) * (Q+511) + 511 * 511);
    a = -1/3.;
    E_min = (int)E_p - 500;
    E_max = (int)E_p + 500;

    Peak *GT1 = new Peak(E_p, 1, Gamma, W0, A, Z, a, Phi_max, Phi_min, E_min, E_max, step);
    TF1 *fGT1 = (TF1*)GT1->GetFunction();

    // // GT2
    E_p = 2499.0*31/32;
    Gamma = 3;//SecondTokeV(2.68e-20);
    Q = 7052.72-1022;
    W0 = sqrt((Q+511) * (Q+511) + 511 * 511);
    a = -1/3.;
    E_min = (int)E_p - 500;
    E_max = (int)E_p + 500;

    Peak *GT2 = new Peak(E_p, 1, Gamma, W0, A, Z, a, Phi_max, Phi_min, E_min, E_max, step);
    TF1 *fGT2 = (TF1*)GT2->GetFunction();

    //GT3
    E_p = 2593.5*31/32;
    Gamma = 3;//SecondTokeV(1.3e-20);
    Q = 5708.92;                          
    W0 = sqrt((Q+511) * (Q+511) + 511 * 511);
    a = -1/3.;
    E_min = (int)E_p - 500;
    E_max = (int)E_p + 500;

    Peak *GT3 = new Peak(E_p, 1, Gamma, W0, A, Z, a, Phi_max, Phi_min, E_min, E_max, step);
    TF1 *fGT3 = (TF1*)GT3->GetFunction();


    // 33Ar one peak to test the lozentizn 
    // double E_p = 5103*32/33;
    // double Gamma = 0;
    // double Q = 7540.0 - 1022; // Q value in keV
    // double W0 = sqrt((Q + 511) * (Q + 511) + 511 * 511);
    // double a = -1/3;
    // double E_min = E_p - 500;
    // double E_max = E_p + 500;
    // Peak *p = new Peak(E_p, 1, Gamma, W0, A, Z, a, Phi_max, Phi_min, E_min, E_max, step);
    // TF1 *fp = (TF1*)p->GetFunction();



    // vector fPeaks = {fp};
    // vector <Peak*> Peaks = {p};


    // // TOTAL Spectrum

    vector <TF1*> fPeaks = {fGT1, fGT2, fGT3};
    vector <Peak*> Peaks = {GT1, GT2, GT3};
    Addition Total(fPeaks);
    TF1 *fTotal = new TF1("fTotal", [&Total](double *x, double *p) {
        return Total.Evaluate(x, p);
    }, Total.GetXmin(), Total.GetXmax(), Total.GetNpar());
    ////////////////////////////////////////

    ////////////// RESOLUTION //////////////
    TF1* Gaussian = new TF1("Gaussian", fGauss, -500, 500, 2);
    ////////////////////////////////////////

    ////////////// CONVOLUTION Spectrum x Resolution //////////////
    Convolution conv2(fTotal, Gaussian, fTotal->GetXmin(), fTotal->GetXmax(), step);
    TF1 *convolution2 = new TF1("convolution2", [&conv2](double *x, double *p) {
        return conv2.Evaluate(x, p);
    }, fTotal->GetXmin(), fTotal->GetXmax(), fTotal->GetNpar()+2);


    for (int i = 0; i < Peaks.size(); ++i)
    {
        Peaks[i]->ExternalFixAllParameters(convolution2, i, {9});
        Peaks[i]->ExternalSetParameter(convolution2, i, 9, 1e2);
        Peaks[i]->ExternalSetParLimits(convolution2, i, 9, 1e-7, 1e7);
        // Peaks[i]->ExternalSetParameter(convolution2, i, 0, 10);
        // Peaks[i]->ExternalSetParLimits(convolution2, i, 0, 0, 20);

        // Peaks[i]->ExternalSetParameter(convolution2, i, 1, fPeaks[i]->GetParameter(1));
        // Peaks[i]->ExternalSetParLimits(convolution2, i, 1, fPeaks[i]->GetParameter(1) - 10, fPeaks[i]->GetParameter(1) + 10);
    }

    convolution2->FixParameter(fTotal->GetNpar(), 1);
    // convolution2->SetParLimits(11*IAS->GetNPeak(), 0, 1e7);
    convolution2->SetParameter(fTotal->GetNpar()+1, 5);
    convolution2->SetParLimits(fTotal->GetNpar()+1, 0, 20);
    
    ////////////////////////////////////////////////////////////////
    
    ////////////// CREATING DATA //////////////
    TFile *file = TFile::Open("/mnt/hgfs/shared-2/2024_DATA/SIMULATED_DATA/06-01/32Ar_ENSDF_CS0_CSP0_CV1_CVP1_analysed.root", "READ");
    TH1D *hh = (TH1D*)(file->Get("p/Initial/E0_p"));
    
    TH1D* h = (TH1D*)hh->Clone("h");
    h->Reset();

    std::cauchy_distribution<double> distribution(0, 10/2);
    std::normal_distribution<double> distribution2(0, 4);
    for (int i = 0; i < hh->GetEntries(); ++i)
    {
        h->Fill(hh->GetRandom() + distribution2(gen));
    }
    //////////////////////////////////////////



    TCanvas *c1 = new TCanvas("c1", "Recoil Broadening", 800, 600);
    h->GetXaxis()->SetRangeUser(fTotal->GetXmin(), fTotal->GetXmax());
    h->Draw();
    c1->SaveAs("recoil_broadening.png");

    TCanvas *c3 = new TCanvas("c3", "Recoil Broadening", 800, 600);
    h->Draw();
    clock_t start = clock();
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(2); // 2 for verbose
    h->Fit(convolution2, "MULTITHREAD RNE");
    convolution2->SetNpx((3400-3300)*10);
    convolution2->SetLineColor(kRed);
    convolution2->Draw("SAME");
    cout << "Time: " << (clock() - start) / (double) CLOCKS_PER_SEC / 8 << endl;

    c3->SetLogy();
    c3->SaveAs("recoil_broadening_convolution.png");

    


    return 0;
}
