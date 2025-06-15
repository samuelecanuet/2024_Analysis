#include "FullFunction.hh"
#include "Math/Minimizer.h"
#include <random>

int main()
{

    // TODO
         // INTEGRATE DATA FROM TXT FILE
         // USING IT IN SPECTRUM.C++
         // REMOVE USELESS PARAMETERS IN TF1
         // GetIntegral FUCNTION AT THIS END 
         // CREATE SAME FOR ALPHAS
         // VERIFY BEHAVIOR WITH RESOLUTION


    

    default_random_engine gen;

    // Spectrum *spectrum = new Spectrum(18, 33, "proton");
    // spectrum->GenerateSpectrum();
    // TF1* fTotal = spectrum->GetFunction();
    // vector <TF1*> fPeaks = spectrum->GetTF1Peaks();
    // vector <Peak*> Peaks = spectrum->GetPeaks();

    // Input parameters
    double step = 0.1; //keV
    double A = 32.;
    double Z = 18.;
    double Phi_max = 0;
    double Phi_min = M_PI;
    

    // ////////////// SPECTRUM //////////////
    // double E_p = 3356.2;                       
    // double Gamma = 20e-3;                       
    // double Q = 5046.6;                          
    // double W0 = sqrt((Q+511) * (Q+511) + 511 * 511);
    // double a = 1;
    // double E_min = (int)E_p - 50;
    // double E_max = (int)E_p + 50;

    // Peak *IAS = new Peak(E_p, 1, Gamma, W0, A, Z, a, Phi_max, Phi_min, E_min, E_max, step);
    // TF1 *fIAS = (TF1*)IAS->GetFunction();

    // n°2
    double E_p = 1209.034025;
    double Gamma = 17; //keV
    double Q = 7256.3 + 1022; // Q value in keV
    double W0 = sqrt((Q + 511) * (Q + 511) + 511 * 511);
    double a = -1./3.;
    double E_min = (int)E_p - 80;
    double E_max = (int)E_p + 80;

    Peak *IAS = new Peak(E_p, 1, Gamma, W0, A, Z, a, Phi_max, Phi_min, E_min, E_max, step);
    TF1 *fIAS = (TF1*)IAS->GetFunction();

    E_p = 1183.034025;
    Gamma = 0; //keV
    Q = 9982 + 1022; // Q value in keV
    W0 = sqrt((Q + 511) * (Q + 511) + 511 * 511);
    a = -1./3.;
    E_min = (int)E_p - 80;
    E_max = (int)E_p + 80;

    Peak *IAS1 = new Peak(E_p, 1, Gamma, W0, A, 17, a, Phi_max, Phi_min, E_min, E_max, step);
    TF1 *fIAS1 = (TF1*)IAS->GetFunction();


    // GT1
    // E_p = 2190.6*31/32;
    // Gamma = 3;//SecondTokeV(1.52e-19);
    // Q = 7362.32-1022;
    // W0 = sqrt((Q+511) * (Q+511) + 511 * 511);
    // a = -1/3.;
    // E_min = (int)E_p - 500;
    // E_max = (int)E_p + 500;

    // Peak *GT1 = new Peak(E_p, 1, Gamma, W0, A, Z, a, Phi_max, Phi_min, E_min, E_max, step);
    // TF1 *fGT1 = (TF1*)GT1->GetFunction();

    // // // GT2
    // E_p = 2499.0*31/32;
    // Gamma = 3;//SecondTokeV(2.68e-20);
    // Q = 7052.72-1022;
    // W0 = sqrt((Q+511) * (Q+511) + 511 * 511);
    // a = -1/3.;
    // E_min = (int)E_p - 500;
    // E_max = (int)E_p + 500;

    // Peak *GT2 = new Peak(E_p, 1, Gamma, W0, A, Z, a, Phi_max, Phi_min, E_min, E_max, step);
    // TF1 *fGT2 = (TF1*)GT2->GetFunction();

    // //GT3
    // E_p = 2593.5*31/32;
    // Gamma = 3;//SecondTokeV(1.3e-20);
    // Q = 5708.92;                          
    // W0 = sqrt((Q+511) * (Q+511) + 511 * 511);
    // a = -1/3.;
    // E_min = (int)E_p - 500;
    // E_max = (int)E_p + 500;

    // Peak *GT3 = new Peak(E_p, 1, Gamma, W0, A, Z, a, Phi_max, Phi_min, E_min, E_max, step);
    // TF1 *fGT3 = (TF1*)GT3->GetFunction();


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


    ////////////// ADDITION Spectrum //////////////
    vector <TF1*> fPeaks = {fIAS, fIAS1};
    vector <Peak*> Peaks = {IAS, IAS1};
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
    TF1 *convolution1 = new TF1("convolution1", [&conv2](double *x, double *p) {
        return conv2.Evaluate(x, p);
    }, fTotal->GetXmin(), fTotal->GetXmax(), fTotal->GetNpar()+Gaussian->GetNpar());
    ///////////////////////////////////////////////////////////////
    
    ////////////// ADDITION (Spectrum x Resolution) + BKG //////////////
    TF1 *bkg = new TF1("bkg", "[0]+x*[1]", fTotal->GetXmin(), fTotal->GetXmax(), 2);
    Addition addition({convolution1, bkg});
    TF1 *convolution2 = new TF1("convolution2", [&addition](double *x, double *p) {
        return addition.Evaluate(x, p);
    }, fTotal->GetXmin(), fTotal->GetXmax(), addition.GetNpar());
    ////////////////////////////////////////////////////////////////////


    for (int i = 0; i < Peaks.size(); ++i)
    {
        Peaks[i]->ExternalFixAllParameters(convolution2, i, {});
        
        // Peaks[i]->ExternalSetParameter(convolution2, i, 9, 0.3);
        // Peaks[i]->ExternalSetParLimits(convolution2, i, 9, 0, 1000);
        // Peaks[i]->ExternalSetParameter(convolution2, i, 1, Peaks[i]->fEnergy);
        // Peaks[i]->ExternalSetParLimits(convolution2, i, 1, Peaks[i]->fEnergy-5, Peaks[i]->fEnergy+5);
        // Peaks[i]->ExternalSetParameter(convolution2, i, 10, 17);
        // Peaks[i]->ExternalSetParLimits(convolution2, i, 10, 10, 25);

        // Peaks[i]->ExternalSetParameter(convolution2, i, 1, fPeaks[i]->GetParameter(1));
        // Peaks[i]->ExternalSetParLimits(convolution2, i, 1, fPeaks[i]->GetParameter(1) - 10, fPeaks[i]->GetParameter(1) + 10);
    }

    convolution2->FixParameter(9, 7.66); // amplitude

    // convolution2->SetParameter(convolution2->GetNpar()-2, 50); // amplitude of the total function
    // convolution2->SetParLimits(convolution2->GetNpar()-2, 0, 1000); // amplitude of the total function
    convolution2->FixParameter(convolution2->GetNpar()-2, 69.4); // amplitude of the total function
    // convolution2->SetParameter(convolution2->GetNpar()-1, -0.03); // amplitude of the total function
    // convolution2->SetParLimits(convolution2->GetNpar()-1, -10, 0); // amplitude of the total function
    convolution2->FixParameter(convolution2->GetNpar()-1, -0.0426); // amplitude of the total function

    convolution2->FixParameter(fTotal->GetNpar(), 1); ///amplitude gaussian resolution
    convolution2->SetParName(fTotal->GetNpar(), "Amp_Resolution");
    convolution2->FixParameter(fTotal->GetNpar()+1, 4.28);
    // convolution2->SetParLimits(fTotal->GetNpar()+1, 0, 10);  // resolution sigma
    convolution2->SetParName(fTotal->GetNpar()+1, "Sigma_Resolution");

    ////////////////////////////////////////////////////////////////

    

    ////////////// CREATING DATA //////////////$
    // from CRADLE
    // TFile *file = TFile::Open("/run/media/local1/Disque_Dur/2024_DATA/SIMULATED_DATA/06-01/32Ar_ENSDF_CS0_CSP0_CV1_CVP1_analysed.root", "READ");
    // TH1D *hh = (TH1D*)(file->Get("Silicon_Detector_Energy_Deposit_D5.1_All"));

    TFile *file = TFile::Open("/mnt/hgfs/shared-2/2025_DATA/DETECTOR_DATA/CALIBRATED/Calibrated_2025.root", "READ");
    TCanvas *c = (TCanvas *)(file->Get("32Ar_SUM_Down"));
    TH1D *hh = (TH1D *)c->GetPrimitive("32Ar_Exp_All_Down");

    TH1D *h = (TH1D *)hh->Clone("h");
    // h->Reset();
    //////////////////////////////////////////



    TFile *result = new TFile("FullFunction_Result.root", "RECREATE");
    result->cd();   
    
    clock_t start = clock();
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(3); // 2 for verbose
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
    h->Fit(convolution2, "MULTITHREAD RN", "" , fTotal->GetXmin(), fTotal->GetXmax());


    //Plot
    TCanvas *c3 = new TCanvas("c3", "Recoil Broadening", 800, 600);
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    h->GetXaxis()->SetRangeUser(fTotal->GetXmin(), fTotal->GetXmax());
    h->Draw("HIST");
    legend->AddEntry(h, "Experimental Data", "l");
    //Complete function
    convolution2->SetNpx(10000);
    convolution2->SetLineColor(kRed);
    convolution2->Draw("SAME");
    double chi2 = convolution2->GetChisquare() / (convolution2->GetNDF() - 1);

    if (Peaks.size() > 1)
    {
        for (int ipeak = 0; ipeak < Peaks.size(); ipeak++)
        {
            Convolution convolution(fPeaks[ipeak], Gaussian, fTotal->GetXmin(), fTotal->GetXmax(), step);
            TF1 *fPeakconv = new TF1("fPeakconv", [&convolution](double *x, double *p) {
                return convolution.Evaluate(x, p);
            }, fTotal->GetXmin(), fTotal->GetXmax(), fPeaks[ipeak]->GetNpar() + Gaussian->GetNpar());
            Addition additionPeak({fPeakconv, bkg});
            TF1 *PeakAndBkg = new TF1("PeakAndBkg", [additionPeak](double *x, double *p) {
                return additionPeak.Evaluate(x, p);
            }, fTotal->GetXmin(), fTotal->GetXmax(), fPeakconv->GetNpar() + bkg->GetNpar());

            PeakAndBkg->SetParameter(0, convolution2->GetParameter(ipeak * fPeaks[ipeak]->GetNpar() + 0));
            PeakAndBkg->SetParameter(1, convolution2->GetParameter(ipeak * fPeaks[ipeak]->GetNpar() + 1));
            PeakAndBkg->SetParameter(2, convolution2->GetParameter(ipeak * fPeaks[ipeak]->GetNpar() + 2));
            PeakAndBkg->SetParameter(3, convolution2->GetParameter(ipeak * fPeaks[ipeak]->GetNpar() + 3));
            PeakAndBkg->SetParameter(4, convolution2->GetParameter(ipeak * fPeaks[ipeak]->GetNpar() + 4));
            PeakAndBkg->SetParameter(5, convolution2->GetParameter(ipeak * fPeaks[ipeak]->GetNpar() + 5));
            PeakAndBkg->SetParameter(6, convolution2->GetParameter(ipeak * fPeaks[ipeak]->GetNpar() + 6));
            PeakAndBkg->SetParameter(7, convolution2->GetParameter(ipeak * fPeaks[ipeak]->GetNpar() + 7));
            PeakAndBkg->SetParameter(8, convolution2->GetParameter(ipeak * fPeaks[ipeak]->GetNpar() + 8));
            PeakAndBkg->SetParameter(9, convolution2->GetParameter(ipeak * fPeaks[ipeak]->GetNpar() + 9));
            PeakAndBkg->SetParameter(10, convolution2->GetParameter(ipeak * fPeaks[ipeak]->GetNpar() + 10));
            PeakAndBkg->SetParameter(11, convolution2->GetParameter(fTotal->GetNpar()));
            PeakAndBkg->SetParameter(12, convolution2->GetParameter(fTotal->GetNpar() + 1));
            PeakAndBkg->SetParameter(13, convolution2->GetParameter(convolution2->GetNpar() - 2));
            PeakAndBkg->SetParameter(14, convolution2->GetParameter(convolution2->GetNpar() - 1));
            PeakAndBkg->SetLineColor(kBlue);
            PeakAndBkg->Draw("SAME");

            legend->AddEntry(PeakAndBkg, ("Peak " + to_string(ipeak + 1)).c_str(), "l");
        }
    }
    bkg->SetParameter(0, convolution2->GetParameter(convolution2->GetNpar() - 2));
    bkg->SetParameter(1, convolution2->GetParameter(convolution2->GetNpar() - 1));
    bkg->SetLineColor(kGreen);
    bkg->Draw("SAME");
    legend->AddEntry(bkg, "Background", "l");
    legend->Draw("SAME");
    c3->SetLogy();
    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.04);
    ostringstream oss;
    oss << setprecision(2) << fixed << chi2;
    string chi2_str = oss.str();
    latex->DrawLatex(0.1, 0.92, ("#chi^{2}_{#nu}: " + chi2_str).c_str());
    latex->Draw("SAME");

    c3->Write();

    cout << "χ2: " << chi2 << endl;
    cout << "Time: " << (clock() - start) / (double) CLOCKS_PER_SEC / 8 << endl;
    convolution2->Write();
    result->Close();


    return 0;
}
