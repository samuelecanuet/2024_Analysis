#ifndef SIPM_PILEUP_HH
#define SIPM_PILEUP_HH

#include "../../../lib/SignalDict/Signal.h"
#include "Detectors.hh"

TFile *MERGED_File;
TFile *FINAL_File;

TTreeReader *Reader;
TTreeReaderArray<Signal> *signals;

TH1D *H_Q1[SIGNAL_MAX];
TH1D *H_Q2[SIGNAL_MAX];
TH2D *H_Q1_Q2[SIGNAL_MAX];
TProfile *P_Q1_Q2[SIGNAL_MAX];
TH2D *H_Q1_Q2_Cleaned[SIGNAL_MAX];
TH2D *H_Q1_Q2_Modified[SIGNAL_MAX];
TH1D *H_Q1_Cleaned[SIGNAL_MAX];
TGraphErrors *graph[SIGNAL_MAX];
TF1 *Correction[SIGNAL_MAX];
TDirectory *dir[SIGNAL_MAX];

double DELTA = 200e3;

void InitHistograms()
{
    Info("Init Histograms");
    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorBetaHigh(det))
        {
            Info(detectorName[det], 1);

            dir[det] = FINAL_File->mkdir(detectorName[det].c_str());

            H_Q1[det] = new TH1D(("H_Q1_" + detectorName[det]).c_str(), ("H_Q1_" + detectorName[det]).c_str(), eHighN, eHighMin, eHighMax);
            H_Q1[det]->GetXaxis()->SetTitle("Q1 [Channel]");
            H_Q1[det]->GetYaxis()->SetTitle("Counts");
            H_Q1[det]->GetXaxis()->CenterTitle();
            H_Q1[det]->GetYaxis()->CenterTitle();

            H_Q2[det] = new TH1D(("H_Q2_" + detectorName[det]).c_str(), ("H_Q2_" + detectorName[det]).c_str(), eHighN, eHighMin, eHighMax);
            H_Q2[det]->GetXaxis()->SetTitle("Q2 [Channel]");
            H_Q2[det]->GetYaxis()->SetTitle("Counts");
            H_Q2[det]->GetXaxis()->CenterTitle();
            H_Q2[det]->GetYaxis()->CenterTitle();

            H_Q1_Q2[det] = new TH2D(("H_Q1_Q2_" + detectorName[det]).c_str(), ("H_Q1_Q2_" + detectorName[det]).c_str(), eHighN/10, eHighMin, eHighMax, eHighN/10, eHighMin, eHighMax);
            H_Q1_Q2[det]->GetXaxis()->SetTitle("Q1 [Channel]");
            H_Q1_Q2[det]->GetYaxis()->SetTitle("Q2 [Channel]");
            H_Q1_Q2[det]->GetXaxis()->CenterTitle();
            H_Q1_Q2[det]->GetYaxis()->CenterTitle();

            H_Q1_Q2_Cleaned[det] = new TH2D(("H_Q1_Q2_Cleaned_" + detectorName[det]).c_str(), ("H_Q1_Q2_Cleaned_" + detectorName[det]).c_str(), eHighN/10, eHighMin, eHighMax, eHighN/10, eHighMin, eHighMax);
            H_Q1_Q2_Cleaned[det]->GetXaxis()->SetTitle("Q1 [Channel]");
            H_Q1_Q2_Cleaned[det]->GetYaxis()->SetTitle("Q2 [Channel]");
            H_Q1_Q2_Cleaned[det]->GetXaxis()->CenterTitle();
            H_Q1_Q2_Cleaned[det]->GetYaxis()->CenterTitle();

            H_Q1_Q2_Modified[det] = new TH2D(("H_Q1_Q2_Modified_" + detectorName[det]).c_str(), ("H_Q1_Q2_Modified_" + detectorName[det]).c_str(), eHighN/10, -1e6, 1e6, eHighN/10, eHighMin, eHighMax);
            H_Q1_Q2_Modified[det]->GetXaxis()->SetTitle("Q1 [Channel]");
            H_Q1_Q2_Modified[det]->GetYaxis()->SetTitle("Q2 [Channel]");
            H_Q1_Q2_Modified[det]->GetXaxis()->CenterTitle();
            H_Q1_Q2_Modified[det]->GetYaxis()->CenterTitle();

            H_Q1_Cleaned[det] = new TH1D(("H_Q1_Cleaned_" + detectorName[det]).c_str(), ("H_Q1_Cleaned_" + detectorName[det]).c_str(), eHighN, eHighMin, eHighMax);
            H_Q1_Cleaned[det]->GetXaxis()->SetTitle("Q1 [Channel]");
            H_Q1_Cleaned[det]->GetYaxis()->SetTitle("Counts");
            H_Q1_Cleaned[det]->GetXaxis()->CenterTitle();
            H_Q1_Cleaned[det]->GetYaxis()->CenterTitle();

            graph[det] = new TGraphErrors();
            
        }

        if (IsDetectorBetaLow(det))
        {
            Info(detectorName[det], 1);

            dir[det] = FINAL_File->mkdir(detectorName[det].c_str());
            H_Q1[det] = new TH1D(("H_Q1_" + detectorName[det]).c_str(), ("H_Q1_" + detectorName[det]).c_str(), eLowN, eLowMin, eLowMax);
            H_Q1[det]->GetXaxis()->SetTitle("Q1 [Channel]");
            H_Q1[det]->GetYaxis()->SetTitle("Counts");
            H_Q1[det]->GetXaxis()->CenterTitle();
            H_Q1[det]->GetYaxis()->CenterTitle();

            H_Q2[det] = new TH1D(("H_Q2_" + detectorName[det]).c_str(), ("H_Q2_" + detectorName[det]).c_str(), eLowN, eLowMin, eLowMax);
            H_Q2[det]->GetXaxis()->SetTitle("Q1 [Channel]");
            H_Q2[det]->GetYaxis()->SetTitle("Q2 [Channel]");
            H_Q2[det]->GetXaxis()->CenterTitle();
            H_Q2[det]->GetYaxis()->CenterTitle();

            H_Q1_Q2[det] = new TH2D(("H_Q1_Q2_" + detectorName[det]).c_str(), ("H_Q1_Q2_" + detectorName[det]).c_str(), eLowN/10, eLowMin, eLowMax, eLowN/10, eLowMin, eLowMax);
            H_Q1_Q2[det]->GetXaxis()->SetTitle("Q1 [Channel]");
            H_Q1_Q2[det]->GetYaxis()->SetTitle("Q2 [Channel]");
            H_Q1_Q2[det]->GetXaxis()->CenterTitle();
            H_Q1_Q2[det]->GetYaxis()->CenterTitle();

            P_Q1_Q2[det] = new TProfile(("P_Q1_Q2_" + detectorName[det]).c_str(), ("P_Q1_Q2_" + detectorName[det]).c_str(), eLowN/10, eLowMin, eLowMax, eLowMin, eLowMax);
            P_Q1_Q2[det]->GetXaxis()->SetTitle("Q1 [Channel]");
            P_Q1_Q2[det]->GetYaxis()->SetTitle("Q2 [Channel]");
            P_Q1_Q2[det]->GetXaxis()->CenterTitle();
            P_Q1_Q2[det]->GetYaxis()->CenterTitle();

            H_Q1_Q2_Cleaned[det] = new TH2D(("H_Q1_Q2_Cleaned_" + detectorName[det]).c_str(), ("H_Q1_Q2_Cleaned_" + detectorName[det]).c_str(), eLowN/10, eLowMin, eLowMax, eLowN/10, eLowMin, eLowMax);
            H_Q1_Q2_Cleaned[det]->GetXaxis()->SetTitle("Q1 [Channel]");
            H_Q1_Q2_Cleaned[det]->GetYaxis()->SetTitle("Q2 [Channel]");
            H_Q1_Q2_Cleaned[det]->GetXaxis()->CenterTitle();
            H_Q1_Q2_Cleaned[det]->GetYaxis()->CenterTitle();

            H_Q1_Q2_Modified[det] = new TH2D(("H_Q1_Q2_Modified_" + detectorName[det]).c_str(), ("H_Q1_Q2_Modified_" + detectorName[det]).c_str(), eLowN/10, -1e6, 1e6, eLowN/10, eLowMin, eLowMax);
            H_Q1_Q2_Modified[det]->GetXaxis()->SetTitle("Q1 [Channel]");
            H_Q1_Q2_Modified[det]->GetYaxis()->SetTitle("Q2 [Channel]");
            H_Q1_Q2_Modified[det]->GetXaxis()->CenterTitle();
            H_Q1_Q2_Modified[det]->GetYaxis()->CenterTitle();

            H_Q1_Cleaned[det] = new TH1D(("H_Q1_Cleaned_" + detectorName[det]).c_str(), ("H_Q1_Cleaned_" + detectorName[det]).c_str(), eLowN, eLowMin, eLowMax);
            H_Q1_Cleaned[det]->GetXaxis()->SetTitle("Q1 [Channel]");
            H_Q1_Cleaned[det]->GetYaxis()->SetTitle("Counts");
            H_Q1_Cleaned[det]->GetXaxis()->CenterTitle();
            H_Q1_Cleaned[det]->GetYaxis()->CenterTitle();



            graph[det] = new TGraphErrors();

        }
    }
}

double FittingFunction(double *x, double *par)
{
    return par[0] * 0.5 * (1-erf( (x[0] - par[1] )/par[2]));
}

double GraphFit(double *x, double *par)
{
    double A = par[0];
    double mu = par[1];
    double sigma = par[2];
    double offset = par[3];
    double L = par[4];

    double delta = par[5];

    double a = A/sqrt(M_PI)/sigma * exp(-pow((L - mu)/sigma, 2));
    double b = par[0] * 0.5 * (offset+erf( (L - par[1] )/par[2])) - a*L;

    if (x[0] < L )
    {
        return a*x[0] + b + delta;
    }
    else
    {
        return par[0] * 0.5 * (offset+erf( (x[0] - par[1] )/par[2])) + delta;
    }
}

double POL2(double *x, double *par)
{
    return par[0] + par[1]*x[0];
}

void FittingQ1Q2()
{   
    TF1 *fit_High = new TF1("fit", FittingFunction, eHighMin, eHighMax, 3); 
    TF1 *fit_Low = new TF1("fit", "gaus", eLowMin, eLowMax);
    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorBetaHigh(det))
        {
            Info(detectorName[det], 1);

            TH1D *H_ProjectionX = H_Q1_Q2[det]->ProjectionX("H_Projection", 10, 10);
            for (int bin = 0; bin < H_Q1_Q2[det]->GetNbinsX(); bin++)
            {
                
                TH1D *H_Projection = H_Q1_Q2[det]->ProjectionY("H_Projection", bin, bin);
                if (H_Projection->GetEntries() < 100) continue;
                H_Projection->GetXaxis()->SetRangeUser(H_Projection->GetBinCenter(H_Projection->GetMaximumBin()), eHighMax);
                fit_High->SetParameter(0, H_Projection->GetMaximum());
                fit_High->SetParameter(1, H_Projection->GetBinCenter(H_Projection->GetMaximumBin()));
                fit_High->SetParLimits(1, H_Projection->GetBinCenter(H_Projection->GetMaximumBin()) - 10e3, H_Projection->GetBinCenter(H_Projection->GetMaximumBin()) + 10e3);
                fit_High->SetParameter(2, 45e3);
                TFitResultPtr r = H_Projection->Fit(fit_High, "QRS");

                if (r->Status() == 0)
                {
                    graph[det]->AddPoint(H_ProjectionX->GetBinCenter(bin), fit_High->GetParameter(1));
                    // graph[det]->AddPointError(bin, H_ProjectionX->GetBinWidth(bin)/2, fit->GetParError(1));
                }
            }

            Correction[det] = new TF1(("PileUp_Correction_"+detectorName[det]).c_str(), GraphFit, eHighMin, eHighMax, 6);
            Correction[det]->SetParameter(0, 2e6);
            Correction[det]->SetParameter(2, 1e5);
            // Correction[det]->SetParLimits(2, 0, 8e5);
            Correction[det]->SetParameter(2, 1e6);
            Correction[det]->SetParameter(4, 3e6);
            Correction[det]->SetParLimits(4, 2.0e6, 3.2e6);
            Correction[det]->FixParameter(5, 0);
            graph[det]->Fit(Correction[det], "R");
        }

        if (IsDetectorBetaLow(det))
        {
            Info(detectorName[det], 1);
            // Correction[det] = new TF1("fit", "pol2", eLowMin, eLowMax);
            // P_Q1_Q2[det]->Fit(Correction[det], "R", "", 60e3, eLowMax);

            TH1D *H_ProjectionX = H_Q1_Q2[det]->ProjectionX("H_Projection", 10, 10);
            for (int bin = 0; bin < H_Q1_Q2[det]->GetNbinsX(); bin++)
            {
                
                TH1D *H_Projection = H_Q1_Q2[det]->ProjectionY("H_Projection", bin, bin);
                double value = H_ProjectionX->GetBinCenter(bin);
                if (H_Projection->GetEntries() < 10) continue;
                fit_Low->SetParameter(1, H_Projection->GetBinCenter(H_Projection->GetMaximumBin()));
                fit_Low->SetParLimits(1, H_Projection->GetBinCenter(H_Projection->GetMaximumBin())-20e3, H_Projection->GetBinCenter(H_Projection->GetMaximumBin())+20e3);
                fit_Low->SetParameter(2, 8e3);
                fit_Low->SetParLimits(2, 1e3, 20e3);
                H_Projection->Fit(fit_Low, "QRS", "", eLowMin, eLowMax);

                // if (r->Status() == 0)
                // {
                    graph[det]->AddPoint(H_ProjectionX->GetBinCenter(bin), fit_Low->GetParameter(1));
                    // graph[det]->SetPointError(graph[det]->GetN(), H_ProjectionX->GetBinWidth(bin)/2, fit_Low->GetParError(1));
                // }
            }

            Correction[det] = new TF1(("PileUp_Correction_"+detectorName[det]).c_str(), POL2, eLowMin, eLowMax, 2);
            Correction[det]->SetParLimits(1, 0, 1);
            Correction[det]->SetParameter(1, 0.5);
            graph[det]->Fit(Correction[det], "R", "", 60e3, eLowMax);          

        }
    }
}

bool PileUp(Signal signal, double coef)
{
    if (abs(Correction[signal.Label]->Eval(signal.Channel) - signal.Pileup) < DELTA/coef)
    {
        return false;
    }
    return true;
}

void WriteHistograms()
{
    TCanvas *c_AllHigh = new TCanvas("AllHigh", "AllHigh", 800, 800);
    c_AllHigh->Divide(3, 3);
    TCanvas *c_AllLow = new TCanvas("AllLow", "AllLow", 800, 800);
    c_AllLow->Divide(3, 3);

    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        gStyle->SetOptStat(0);
        if (IsDetectorBetaHigh(det))
        {
            dir[det]->cd();
            H_Q1[det]->Write();
            H_Q2[det]->Write();
            H_Q1_Q2[det]->Write();
            Correction[det]->Write();

            TCanvas *c = new TCanvas(detectorName[det].c_str(), detectorName[det].c_str(), 800, 800);   
            H_Q1_Q2[det]->Draw("COLZ");
            graph[det]->Draw("P SAME");
            Correction[det]->SetParameter(5, DELTA);
            TGraph *g = new TGraph();
            for (int bin = 0; bin < eHighN; bin++)
            {
                g->AddPoint(H_Q1_Q2[det]->GetXaxis()->GetBinCenter(bin), Correction[det]->Eval(H_Q1_Q2[det]->GetXaxis()->GetBinCenter(bin)));
            }
            g->SetLineStyle(2);
            g->SetLineWidth(2);
            g->SetLineColor(kRed);
            g->Draw("SAME");
            Correction[det]->SetParameter(5, -DELTA);
            TGraph *g2 = new TGraph();
            for (int bin = 0; bin < eHighN; bin++)
            {
                g2->AddPoint(H_Q1_Q2[det]->GetXaxis()->GetBinCenter(bin), Correction[det]->Eval(H_Q1_Q2[det]->GetXaxis()->GetBinCenter(bin)));
            }
            g2->SetLineStyle(2);
            g2->SetLineWidth(2);
            g2->SetLineColor(kRed);
            g2->Draw("SAME");            


            c->Write();

            H_Q1_Q2_Cleaned[det]->Write();
            H_Q1_Q2_Modified[det]->Write();

            TCanvas *c2 = new TCanvas(("Corrected_"+detectorName[det]).c_str(), ("Corrected_"+detectorName[det]).c_str(), 800, 800);
            H_Q1[det]->Draw("HIST");
            H_Q1_Cleaned[det]->SetLineColor(kRed);
            H_Q1_Cleaned[det]->Draw("HIST SAME");
            c2->Write();

            c_AllHigh->cd(GetDetectorChannel(det));
            TPad *pad = (TPad*)c_AllHigh->GetPad(GetDetectorChannel(det));
            pad->SetLogz();
            H_Q1_Q2[det]->GetXaxis()->SetRangeUser(0, 5000e3);
            H_Q1_Q2[det]->GetYaxis()->SetRangeUser(0, 2800e3);
            H_Q1_Q2[det]->Draw("COLZ");
            g->Draw("SAME");
            g2->Draw("SAME");
        }

        if (IsDetectorBetaLow(det))
        {
            dir[det]->cd();
            H_Q1[det]->Write();
            H_Q2[det]->Write();
            H_Q1_Q2[det]->Write();
            P_Q1_Q2[det]->Write();
            Correction[det]->Write();

            TCanvas *c = new TCanvas(detectorName[det].c_str(), detectorName[det].c_str(), 800, 800);
            H_Q1_Q2[det]->Draw("COLZ");
            graph[det]->Draw("P SAME");
            Correction[det]->SetParameter(0, Correction[det]->GetParameter(0) + DELTA/10);
            TGraph *g = new TGraph();
            for (int bin = 0; bin < eLowN; bin++)
            {
                g->AddPoint(H_Q1_Q2[det]->GetXaxis()->GetBinCenter(bin), Correction[det]->Eval(H_Q1_Q2[det]->GetXaxis()->GetBinCenter(bin)));
            }
            g->SetLineStyle(2);
            g->SetLineWidth(2);
            g->SetLineColor(kRed);
            g->Draw("SAME");
            Correction[det]->SetParameter(0, Correction[det]->GetParameter(0) - 2*DELTA/10);
            TGraph *g2 = new TGraph();
            for (int bin = 0; bin < eLowN; bin++)
            {
                g2->AddPoint(H_Q1_Q2[det]->GetXaxis()->GetBinCenter(bin), Correction[det]->Eval(H_Q1_Q2[det]->GetXaxis()->GetBinCenter(bin)));
            }
            g2->SetLineStyle(2);
            g2->SetLineWidth(2);
            g2->SetLineColor(kRed);
            g2->Draw("SAME");

            c->Write();

            H_Q1_Q2_Cleaned[det]->Write();
            H_Q1_Q2_Modified[det]->Write();

            TCanvas *c2 = new TCanvas(("Corrected_"+detectorName[det]).c_str(), ("Corrected_"+detectorName[det]).c_str(), 800, 800);
            H_Q1[det]->Draw("HIST");
            H_Q1_Cleaned[det]->SetLineColor(kRed);
            H_Q1_Cleaned[det]->Draw("HIST SAME");
            c2->Write();

            c_AllLow->cd(GetDetectorChannel(det));
            TPad *pad = (TPad*)c_AllLow->GetPad(GetDetectorChannel(det));
            pad->SetLogz();
            H_Q1_Q2[det]->GetXaxis()->SetRangeUser(0, 1000e3);
            H_Q1_Q2[det]->GetYaxis()->SetRangeUser(0, 800e3);
            H_Q1_Q2[det]->Draw("COLZ");
            g->Draw("SAME");
            g2->Draw("SAME");
        }
    }

    FINAL_File->cd();   
    c_AllLow->Write();
    c_AllHigh->Write();
    

    FINAL_File->Close();
}

#endif