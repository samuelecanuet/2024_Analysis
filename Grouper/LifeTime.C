#include "TFile.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"


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

int LifeTime()
{
    TFile *f = new TFile("/mnt/hgfs/shared-2/2025_DATA/DETECTOR_DATA/ANALYSED/Spectrum.root", "READ");
    
    TFile *fout = new TFile("/mnt/hgfs/shared-2/Lifetime_2D.root", "RECREATE");

    // 32Ar
    TDirectory *dir = f->GetDirectory("32Ar");
    TIter next(dir->GetListOfKeys());
    TKey *key;
    TH1D *H32Ar = nullptr;
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
    TH1D *H33Ar = nullptr;
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

    // doing 2D
    double t_half_32Ar = 0.080;
    double t_half_33Ar = 0.160;
    double t_half_range = 0.04;
    double t_half_step = 0.002;

    fout->cd();
    gErrorIgnoreLevel = kWarning;

    double chi2_min = std::numeric_limits<double>::max();

    int TOTAL = 2*t_half_range/t_half_step*2*t_half_range/t_half_step;
    int counter = 0;
    TH2D *H2D = new TH2D("H2D", "H2D", 2.*t_half_range/t_half_step, t_half_32Ar-t_half_range, t_half_32Ar+t_half_range, 2.*t_half_range/t_half_step, t_half_33Ar-t_half_range, t_half_33Ar+t_half_range);
    for (double t_half32 = -t_half_range; t_half32 < t_half_range; t_half32 += t_half_step)
    {
        TH1D *H32Ar_Conv = ConvolveWithExponential(H33Ar, t_half_32Ar+t_half32);
        H32Ar_Conv->Scale(1./H32Ar_Conv->Integral());
        for (double t_half33 = -t_half_range; t_half33 < t_half_range; t_half33 += t_half_step)
        {
            cout << counter << " / " << TOTAL << "\r";
            counter++;
            TH1D *H33Ar_Conv = ConvolveWithExponential(H32Ar, t_half_33Ar+t_half33);
            H33Ar_Conv->Scale(1./H33Ar_Conv->Integral());
            
            H33Ar_Conv->GetXaxis()->SetRangeUser(0.4, 1.2);
            H32Ar_Conv->GetXaxis()->SetRangeUser(0.4, 1.2);
            double chi2 = H32Ar_Conv->Chi2Test(H33Ar_Conv, "WW CHI2/NDF");

            chi2_min = std::min(chi2_min, chi2);


            TCanvas *c = new TCanvas(Form("c_%f_%f", t_half_32Ar+t_half32, t_half_33Ar+t_half33), Form("c_%f_%f", t_half_32Ar+t_half32, t_half_33Ar+t_half33), 800, 600);
            H32Ar_Conv->SetLineColor(kRed);
            H32Ar_Conv->Draw("HIST");
            H33Ar_Conv->SetLineColor(kBlue);
            H33Ar_Conv->Draw("HIST SAME");
            c->Write();

            H2D->Fill(t_half_32Ar+t_half32, t_half_33Ar+t_half33, chi2);
        }
    }

    gErrorIgnoreLevel = kInfo;

    TCanvas *c = new TCanvas("c", "c", 800, 600);
    H2D->Draw("COLZ");
    c->Write();

    fout->Close();

    return 0;
}


