#include <TCanvas.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TMath.h>

#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <iostream>
#include <cmath>

// 2D Gaussian function
double Gaussian2D(double x, double y, double mu_x, double mu_y, double sigma_x, double sigma_y) {
    double norm = 1.0 / (2 * M_PI * sigma_x * sigma_y);
    double exponent = -0.5 * (pow((x - mu_x) / sigma_x, 2) + pow((y - mu_y) / sigma_y, 2));
    return norm * exp(exponent);
}

// Convolution of histogram h(x, y) with Gaussian, evaluated at y0
TH1D* ConvolveAndProject(TH2D* histo, double sigma_x, double sigma_y, double y0) {
    int nbinsX = histo->GetNbinsX();
    int nbinsY = histo->GetNbinsY();

    double x_min = histo->GetXaxis()->GetXmin();
    double x_max = histo->GetXaxis()->GetXmax();
    TH1D* projectedHist = new TH1D("proj", "Projected Convolution; X; Convolved Counts", nbinsX, x_min, x_max);

    double bin_width_x = histo->GetXaxis()->GetBinWidth(1);
    double bin_width_y = histo->GetYaxis()->GetBinWidth(1);

    // Loop over x bins to compute f(x, y0)
    for (int i = 1; i <= nbinsX; i++) {
        double x = histo->GetXaxis()->GetBinCenter(i);
        double sum = 0.0;

        // Perform convolution over all (x', y') in histogram
        for (int j = 1; j <= nbinsX; j++) {
            double x_prime = histo->GetXaxis()->GetBinCenter(j);
            for (int k = 1; k <= nbinsY; k++) {
                double y_prime = histo->GetYaxis()->GetBinCenter(k);
                double h_val = histo->GetBinContent(j, k);
                double gauss_weight = Gaussian2D(x, y0, x_prime, y_prime, sigma_x, sigma_y);
                sum += h_val * gauss_weight * bin_width_x * bin_width_y; // Discrete convolution
            }
        }
        projectedHist->SetBinContent(i, sum);
    }

    return projectedHist;
}

void Test2() {
    // Example: Create a 2D histogram
    TH2D* histo = new TH2D("histo", "Example Histogram; X; Y", 50, -5, 5, 50, -5, 5);
    
    // Fill histogram with some values (e.g., Gaussian blob centered at (0,0))
    for (int i = 0; i < histo->GetNbinsX(); i++) {
        histo->SetBinContent(i, i, 1);
    }

    // Define Gaussian smoothing parameters
    double sigma_x = 1.0;
    double sigma_y = 1.0;
    double y0 = 0.0;

    // Compute the convolved projection
    TH1D* result = ConvolveAndProject(histo, sigma_x, sigma_y, y0);

    // Draw results
    TCanvas* c = new TCanvas("c", "Convolution Result", 800, 600);
    result->Draw();

    TCanvas *c2 = new TCanvas("c2", "c2", 1920, 1080);
    histo->Draw("COLZ");
    c2->Write();
}
