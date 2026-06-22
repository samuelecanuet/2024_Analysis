#include "TFile.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "include/Utilities.hh"


// TH1D *ConvolveWithExponential(TH1D *hist, double t_half)
// {
//     double lambda = log(2.) / t_half; // decay constant
//     int nbins = hist->GetNbinsX();
//     double bin_width = hist->GetBinWidth(1);

//     // Prepare the convolved histogram
//     TH1D *convolved_hist = (TH1D *)hist->Clone((string(hist->GetName()) + "_Conv").c_str());
//     convolved_hist->Reset();

//     // Perform discrete convolution
//     for (int i = 1; i <= nbins; i++)
//     {
//         double sum = 0;
//         double t_i = hist->GetBinCenter(i);

//         // Convolve over all previous bins
//         for (int j = 1; j <= i; j++)
//         {
//             double t_j = hist->GetBinCenter(j);
//             double decay = exp(-lambda * (t_i - t_j));
//             sum += hist->GetBinContent(j) * decay * bin_width;
//         }

//         convolved_hist->SetBinContent(i, sum);
//     }
//     return convolved_hist;
// }


TH1D *ConvolveWithExponential(const TH1D *hist, double t_half, bool normalize_kernel = false)
{
    if (!hist || t_half <= 0)
        return nullptr;

    const double lambda = 0.69314718056 / t_half;

    const int nbins = hist->GetNbinsX();

    TH1D *convolved_hist = (TH1D *)hist->Clone(
        (std::string(hist->GetName()) + "_Conv").c_str()
    );

    convolved_hist->Reset("ICES");

    for (int i = 1; i <= nbins; i++)
    {
        const double t_i = hist->GetBinCenter(i);
        double sum = 0.0;
        double err2 = 0.0;

        for (int j = 1; j <= i; j++)
        {
            const double t_j = hist->GetBinCenter(j);
            const double dt = t_i - t_j;

            double kernel = std::exp(-lambda * dt);

            if (normalize_kernel)
                kernel *= lambda;

            sum += hist->GetBinContent(j) * kernel * hist->GetBinWidth(j);

            double w = kernel * hist->GetBinWidth(j);
            err2 += std::pow(hist->GetBinError(j) * w, 2);
        }

        convolved_hist->SetBinContent(i, sum);
        convolved_hist->SetBinError(i, std::sqrt(err2));
    }

    return convolved_hist;
}

void MinimizationChi2(const TH1D *H32Ar, const TH1D *H33Ar, double t_half_32Ar_min, double t_half_33Ar_min, double t_half_range, double t_half_step, TH2D *H2D, double &t_half_32Ar_result, double &t_half_33Ar_result, TH1D *&H32Ar_Conv_final, TH1D *&H33Ar_Conv_final)
{
    double chi2_min = std::numeric_limits<double>::max();
    int counter = 0;
    int TOTAL = 2*t_half_range/t_half_step*2*t_half_range/t_half_step;
    for (double t_half32 = t_half_32Ar_min - t_half_range; t_half32 < t_half_32Ar_min + t_half_range; t_half32 += t_half_step)
    {
        TH1D *H33Ar_Conv = ConvolveWithExponential(H33Ar, t_half32);
        H33Ar_Conv->Scale(1./H33Ar_Conv->Integral());
        H33Ar_Conv->GetXaxis()->SetRangeUser(0.1, 1.2);
        for (double t_half33 = t_half_33Ar_min - t_half_range; t_half33 < t_half_33Ar_min + t_half_range; t_half33 += t_half_step)
        {
            cout << counter << " / " << TOTAL << "\r";
            counter++;
            TH1D *H32Ar_Conv = ConvolveWithExponential(H32Ar, t_half33);
            H32Ar_Conv->Scale(1./H32Ar_Conv->Integral());
                        
            H32Ar_Conv->GetXaxis()->SetRangeUser(0.1, 1.2);
            double chi2 = H32Ar_Conv->Chi2Test(H33Ar_Conv, "CHI2");

            H2D->Fill(t_half32, t_half33, chi2);

            if (chi2 < chi2_min)
            {
                chi2_min = chi2;
                t_half_32Ar_result = t_half32;
                t_half_33Ar_result = t_half33;
                H32Ar_Conv_final = (TH1D *)H32Ar_Conv->Clone("H32Ar_Conv_final");
                H33Ar_Conv_final = (TH1D *)H33Ar_Conv->Clone("H33Ar_Conv_final");
            }
        }
    }
}

int LifeTime(string OPTION="")
{
    if (OPTION != "raw" && OPTION != "cut")
    {
        Error("Option not recognized : " + OPTION);
        return -1;
    }


    TFile *fout;
    if (OPTION == "raw")
    {
        fout = new TFile("/mnt/hgfs/shared-2/Lifetime_Raw.root", "RECREATE");
    }
    else if (OPTION == "cut")
    {
        fout = new TFile("/mnt/hgfs/shared-2/Lifetime_Cut.root", "RECREATE");
    }

    TH1D *H32Ar = nullptr;
    TH1D *H33Ar = nullptr;

    TFile *f;
    if (OPTION == "raw")
    {
        f = MyTFile("/mnt/hgfs/shared-2/2025_DATA/DETECTOR_DATA/ANALYSED/Release_RAW.root", "READ");
        H32Ar = (TH1D *)f->Get("H_Release_RAW_32Ar");
        H33Ar = (TH1D *)f->Get("H_Release_RAW_33Ar");
    }
    else if (OPTION == "cut")
    {
        f = MyTFile("/mnt/hgfs/shared-2/2025_DATA/DETECTOR_DATA/ANALYSED/Spectrum (copy).root", "READ");

        // 32Ar
        TDirectory *dir = f->GetDirectory("32Ar");
        TIter next(dir->GetListOfKeys());
        TKey *key;
        
        while ((key = (TKey *)next()))
        {
            if (string(key->GetName()).find("_14_") != string::npos || string(key->GetName()).find("_14.000000_") != string::npos)
            {
                TCanvas *c = (TCanvas *)key->ReadObj();
                for (int i = 0; i < c->GetListOfPrimitives()->GetEntries(); i++)
                {
                    if (string(c->GetListOfPrimitives()->At(i)->GetName()).find("V_3") != string::npos)
                    {
                        TPad *pad = (TPad *)c->GetListOfPrimitives()->At(i);
                        H32Ar = (TH1D *)pad->GetPrimitive("H_Release_32Ar_Down_14");
                        H32Ar->Add((TH1D *)pad->GetPrimitive("H_Release_32Ar_Up_14"));
                    }
                }
            }
        }

        // 33Ar
        dir = f->GetDirectory("33Ar");
        next = TIter(dir->GetListOfKeys());
        
        while ((key = (TKey *)next()))
        {
            if (string(key->GetName()).find("_21_") != string::npos || string(key->GetName()).find("_21.000000_") != string::npos)
            {
                TCanvas *c = (TCanvas *)key->ReadObj();
                for (int i = 0; i < c->GetListOfPrimitives()->GetEntries(); i++)
                {
                    if (string(c->GetListOfPrimitives()->At(i)->GetName()).find("V_3") != string::npos)
                    {
                        TPad *pad = (TPad *)c->GetListOfPrimitives()->At(i);
                        H33Ar = (TH1D *)pad->GetPrimitive("H_Release_33Ar_Down_21");
                        H33Ar->Add((TH1D *)pad->GetPrimitive("H_Release_33Ar_Up_21"));
                    }
                }
            }
        }
    }    

    // doing 2D
    H33Ar->Rebin(5);
    H32Ar->Rebin(5);
    /// ROUGHT SCAN TO FIND MINIMUM
    double t_half_32Ar = 0.070;
    double t_half_33Ar = 0.170;
    double t_half_range = 0.04;
    double t_half_step = 0.005;

    double t_half_32Ar_min = 0;
    double t_half_33Ar_min = 0;

    fout->cd();
    gErrorIgnoreLevel = kWarning;

    TH1D *H32Ar_Conv_final = nullptr;
    TH1D *H33Ar_Conv_final = nullptr;
    TH2D *H2D = new TH2D("H2D", "H2D", 2.*t_half_range/t_half_step, t_half_32Ar-t_half_range, t_half_32Ar+t_half_range, 2.*t_half_range/t_half_step, t_half_33Ar-t_half_range, t_half_33Ar+t_half_range);

    cout << "Starting rought minimization..." << endl;
    MinimizationChi2(H32Ar, H33Ar, t_half_32Ar, t_half_33Ar, t_half_range, t_half_step, H2D, t_half_32Ar_min, t_half_33Ar_min, H32Ar_Conv_final, H33Ar_Conv_final);

    cout << "Rought minimization done, minimum found at t_half_32Ar = " << t_half_32Ar_min << " s and t_half_33Ar = " << t_half_33Ar_min << " s." << endl;
    /// FINE SCAN AROUND MINIMUM
    t_half_32Ar_min = (int)(t_half_32Ar_min*1000)/1000.;
    t_half_33Ar_min = (int)(t_half_33Ar_min*1000)/1000.;
    cout << "Starting fine minimization around minimum..." << endl;
    t_half_range = 0.02;
    t_half_step = 0.0005;
    TH2D *H2D_fine = new TH2D("H2D_fine", "H2D_fine", 2.*t_half_range/t_half_step, t_half_32Ar_min-t_half_range, t_half_32Ar_min+t_half_range, 2.*t_half_range/t_half_step, t_half_33Ar_min-t_half_range, t_half_33Ar_min+t_half_range);
    double t_half_32Ar_result = 0;
    double t_half_33Ar_result = 0;
    H32Ar_Conv_final = nullptr;
    H33Ar_Conv_final = nullptr;
    MinimizationChi2(H32Ar, H33Ar, t_half_32Ar_min, t_half_33Ar_min, t_half_range, t_half_step, H2D_fine, t_half_32Ar_result, t_half_33Ar_result, H32Ar_Conv_final, H33Ar_Conv_final);

    cout << "Fine minimization done, minimum found at t_half_32Ar = " << t_half_32Ar_result << " s and t_half_33Ar = " << t_half_33Ar_result << " s." << endl;

    gErrorIgnoreLevel = kInfo;

    TCanvas *c = new TCanvas("c", "c", 800, 600);
    H2D->Draw("COLZ");
    c->Write();
    TCanvas *c_fine = new TCanvas("c_fine", "c_fine", 800, 600);
    H2D_fine->Draw("COLZ");
    TLine *line32 = new TLine(0.098, t_half_33Ar_min - t_half_range, 0.098, t_half_33Ar_min + t_half_range);
    line32->SetLineColor(kRed);
    line32->Draw("SAME");

    TLine *line33 = new TLine(t_half_32Ar_min - t_half_range, 0.174, t_half_32Ar_min + t_half_range, 0.174);
    line33->SetLineColor(kRed);
    line33->Draw("SAME");
    c_fine->Write();
    TCanvas *c_all = new TCanvas("c_all", "c_all", 800, 600);
    H32Ar->SetLineColor(kRed);
    H33Ar->SetLineColor(kBlue);
    H32Ar->Scale(1./H32Ar->Integral());
    H33Ar->Scale(1./H33Ar->Integral());
    H32Ar->Draw("HIST");
    H33Ar->Draw("HIST SAME");
    H32Ar_Conv_final->SetLineColor(kRed);
    H33Ar_Conv_final->SetLineColor(kBlue);
    H32Ar_Conv_final->Draw("HIST SAME");
    H33Ar_Conv_final->Draw("HIST SAME");
    c_all->Write();
    H32Ar->Write();
    H33Ar->Write();
    H32Ar_Conv_final->Write();
    H33Ar_Conv_final->Write();

    fout->Close();

    return 0;
}


