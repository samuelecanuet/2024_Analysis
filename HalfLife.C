
#include <cmath>
#include <vector>
#include <iostream>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraph.h>
#include "TF1.h"

// Reconstruit g(t) à partir de h(t) et d’un T1/2 supposé
TH1D* ReconstructG(const TH1D* h, double T12) {
    const int N = h->GetNbinsX();
    const double dt = h->GetXaxis()->GetBinWidth(1);
    const double lambda = 0.6930 / T12;

    TH1D* g = (TH1D*)h->Clone("g");
    g->Reset();

    double ref = 0.0;

    for (int n = 1; n < N; ++n) 
    {
        double tn = h->GetXaxis()->GetBinCenter(n);

        if (tn < 0.02) {
            continue; 
        }

        // first bin
        double Height = h->GetBinContent(n);
        TF1 *f = new TF1("f", "x < [1] ? 0.0 : exp(-[0]*(x-[1]))", 0, 1.2);
        f->SetNpx(N);
        f->SetParameter(0, lambda);
        f->SetParameter(1, tn);
        TH1D* hexpo = (TH1D*)f->GetHistogram();
        g->Add(hexpo, Height-g->GetBinContent(n) > 0.0 ? Height-g->GetBinContent(n) : 0.0);

        g->SetBinError(n, sqrt(g->GetBinContent(n)));
    }


    return g;
}

TH1D* ReconstructGG(const TH1D* h, double T12) {
    const int N = h->GetNbinsX();
    const double dt = h->GetXaxis()->GetBinWidth(1);
    const double lambda = 0.6930 / T12;

    TH1D* gg = (TH1D*)h->Clone("g");
    gg->Reset();

    TH1D* g = (TH1D*)h->Clone("g");
    g->Reset();


    double ref = 0.0;

    for (int n = 1; n < N; ++n) 
    {
        double tn = h->GetXaxis()->GetBinCenter(n);

        if (tn < 0.02) {
            continue; 
        }

        // first bin
        double Height = h->GetBinContent(n);
        TF1 *f = new TF1("f", "x < [1] ? 0.0 : exp(-[0]*(x-[1]))", 0, 1.2);
        f->SetNpx(N);
        f->SetParameter(0, lambda);
        f->SetParameter(1, tn);
        TH1D* hexpo = (TH1D*)f->GetHistogram();
        double height_release = Height-g->GetBinContent(n);
        g->Add(hexpo, height_release);
        g->SetBinError(n, sqrt(g->GetBinContent(n)));

        gg->SetBinContent(n,  height_release);
    }


    return gg;
}
TH1D* ReconstructH(const TH1D* g, double T12) {
    const int N = g->GetNbinsX();
    const double dt = g->GetXaxis()->GetBinWidth(1);
    const double lambda = 0.6930 / T12;

    TH1D* h = (TH1D*)g->Clone("g");
    h->Reset();

    for (int n = 1; n < N; ++n) 
    {
        double tn = h->GetXaxis()->GetBinCenter(n);

        if (tn < 0.02) {
            continue; 
        }

        // first bin
        double Height = g->GetBinContent(n);
        TF1 *f = new TF1("f", "x < [1] ? 0.0 : exp(-[0]*(x-[1]))", 0, 1.2);
        f->SetNpx(N);
        f->SetParameter(0, lambda);
        f->SetParameter(1, tn);
        TH1D* hexpo = (TH1D*)f->GetHistogram();
        h->Add(hexpo, Height);
        h->SetBinError(n, sqrt(h->GetBinContent(n)));
    }


    return h;
}

// Critère de qualité : somme des carrés négatifs de g(t)
double CostFunction(TH1D* h, double T12, TH1D* h32) {
    TH1D* g = ReconstructG(h, T12);
    TH1D* gg = ReconstructGG(h, T12);
    TH1D* hh = ReconstructH(gg, 0.1743);

    TCanvas *c = new TCanvas("c", "Reconstructed g(t)", 800, 600);
    c->cd();
    

    hh->Scale(h32->Integral() / hh->Integral());
    hh->SetLineColor(kRed);
    hh->Draw("HIST");

    h32->SetLineColor(kBlack);
    h32->Draw("SAME");

    c->SaveAs("Reconstructed_g.root");


    double chi2 = h32->Chi2Test(hh, "SCALE CHI2/NDF");
    cout << "Chi2 for T1/2 = " << T12 << " s: " << chi2 << endl;
    return chi2;
}

int HalfLife()
{


    TFile *f =  new TFile("/mnt/hgfs/shared-2/2025_DATA/DETECTOR_DATA/ANALYSED/Spectrum (copy).root", "READ");

    TCanvas *c = (TCanvas *)f->Get("32Ar/32Ar_Peak_14_3356.200196keV");
    if (!c) {
        std::cerr << "Canvas not found!" << std::endl;
        return 1;
    }

    TH1D* h;
    for (int i = 0; i < c->GetListOfPrimitives()->GetEntries(); ++i) {
        TPad *obj = (TPad *)c->GetListOfPrimitives()->At(i);
        string name = obj->GetName();
        if (name == "32Ar_Peak_14_3356.200196keV_3") {
            for (int j = 0; j < obj->GetListOfPrimitives()->GetEntries(); ++j) {
                auto subObj = obj->GetListOfPrimitives()->At(j);
                if (subObj->InheritsFrom("TH1D")) {
                    h = (TH1D *)subObj;
                    break;
                }
            }
        }
    }

    if(!h) {
        std::cerr << "Histogram not found!" << std::endl;
        return 1;
    }


    TCanvas *c33 = (TCanvas *)f->Get("33Ar/33Ar_Peak_21_3171.859004keV");
    if (!c33) {
        std::cerr << "Canvas not found!" << std::endl;
        return 1;
    }

    TH1D* h33;
    for (int i = 0; i < c33->GetListOfPrimitives()->GetEntries(); ++i) {
        TPad *obj = (TPad *)c33->GetListOfPrimitives()->At(i);
        string name = obj->GetName();
        if (name == "33Ar_Peak_21_3171.859004keV_3") {
            for (int j = 0; j < obj->GetListOfPrimitives()->GetEntries(); ++j) {
                auto subObj = obj->GetListOfPrimitives()->At(j);
                if (subObj->InheritsFrom("TH1D")) {
                    h33 = (TH1D *)subObj;
                    break;
                }
            }
        }
    }

    if(!h33) {
        std::cerr << "Histogram not found!" << std::endl;
        return 1;
    }


    // Balayage de T1/2 entre 10 ms et 500 ms
    const int Nscan = 100;
    std::vector<double> T12_list, cost_list;

    h33->Rebin(20); // Rebinning pour améliorer la résolution
    h->Rebin(20); // Rebinning pour améliorer la résolution

    double cost = CostFunction(h, 0.098, h33);

    // for (int i = 0; i < Nscan; ++i) {
    //     cout << "T1/2 = " << 0.01 + i * (0.2 - 0.01) / (Nscan - 1) << " s" << endl;
    //     double T12 = 0.01 + i * (0.2 - 0.01) / (Nscan - 1); // en secondes
    //     double cost = CostFunction(h, T12);
    //     T12_list.push_back(T12 * 1000.0);     // en ms
    //     cost_list.push_back(cost); // log-échelle
    // }

    // Tracé du critère de coût
    TGraph* gr = new TGraph(Nscan, T12_list.data(), cost_list.data());
    gr->SetTitle("Scan T_{1/2};T_{1/2} (ms);log_{10}(coût)");
    gr->SetLineColor(kBlue);
    gr->SetLineWidth(2);

    TCanvas* c1 = new TCanvas("c1", "Scan", 800, 600);
    gr->Draw("AL");
    return 1;
}

