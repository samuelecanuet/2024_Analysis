#ifndef MATCHER_HH
#define MATCHER_HH

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TCut.h"
#include "TCutG.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TRandom.h"
#include "TKey.h"
#include "TMinuit.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "TF1Convolution.h"
#include "TSpectrum.h"

#include <gsl/gsl_statistics.h>

#include "../../../lib/SignalDict/Signal.h"
#include "Detectors.hh"

using namespace std;
using namespace ROOT::Math;
default_random_engine generator;

/// FILE ///

TFile *ROOT_file;
TFile *FINAL_file;
TTreeReader *Reader;
TTree *Tree;
TTree *Cleaned_Tree;

TH1D *H_RAW[5];
TH2D *H_RAW_2D;
TH2D *H;
TH2D *H_precorrrected;
TH2D *H_corrected;
TH2D *H_reconstruction;
TH2D *H_measurement;
TH2D *H_measurement_reconstruted;

int A = 1;
int B = 3;
int C = 4;
int D = 2;

// double rho_real = 2;
// double l_real = 1.2;
double rho = 2.0;
double l = 1.2;
double l2 = l / 2;
double cercle = 7.5;
int n = 8;
// double limit = limit_real * rho / rho_real;
int MAXI = 10;
int DegX = 2;
int DegY = 2;
TF2 *FittedFunction;
TF2 *FinalFunction;
TF2 *MeasurementFunction;
TF2 *f_X;
TF2 *f_Y;

double init_par[5] = {0, 1, 1, 1, 1};
double data[5] = {0, 0, 0, 0, 0};

bool M8[8][8] = {
    {false, false, false, false, false, false, false, false},
    {false, false, false, false, false, false, false, false},
    {false, false, true, true, true, true, false, false},
    {false, false, true, true, true, true, false, false},
    {false, false, true, true, true, true, false, false},
    {false, false, true, true, true, true, false, false},
    {false, false, false, false, false, false, false, false},
    {false, false, false, false, false, false, false, false}};

    // 6x6 without corners
    // {false, false, false, false, false, false, false, false},
    // {false, false, true, true, true, true, false, false},
    // {false, true, true, true, true, true, true, false},
    // {false, true, true, true, true, true, true, false},
    // {false, true, true, true, true, true, true, false},
    // {false, true, true, true, true, true, true, false},
    // {false, false, true, true, true, true, false, false},
    // {false, false, false, false, false, false, false, false}};

    // FULL
    // {false, false, true, true, true, true, false, false},
    // {false, true, true, true, true, true, true, false},
    // {true, true, true, true, true, true, true, true},
    // {true, true, true, true, true, true, true, true},
    // {true, true, true, true, true, true, true, true},
    // {true, true, true, true, true, true, true, true},
    // {false, true, true, true, true, true, true, false},
    // {false, false, true, true, true, true, false, false}};

bool M_final[8][8] = {
    {false, false, false, false, false, false, false, false},
    {false, false, false, false, false, false, false, false},
    {false, false, true, true, true, true, false, false},
    {false, false, true, true, true, true, false, false},
    {false, false, true, true, true, true, false, false},
    {false, false, true, true, true, true, false, false},
    {false, false, false, false, false, false, false, false},
    {false, false, false, false, false, false, false, false}};

bool M_measurement[8][8] = {
    {false, false, false, false, false, false, false, false},
    {false, false, false, false, false, false, false, false},
    {false, false, false, false, false, false, false, false},
    {false, false, false, true, true, false, false, false},
    {false, false, false, true, true, false, false, false},
    {false, false, false, false, false, false, false, false},
    {false, false, false, false, false, false, false, false},
    {false, false, false, false, false, false, false, false}};

pair<double, double> M8_center[64];

/////////////////////////////////////////////////////////////////////
//////////// NEW 2D FULL FITTING ////////////////////////////////////
/////////////////////////////////////////////////////////////////////

///////////////// Fitting Points /////////////////////////////////////
void LoadCenters()
{
    ifstream file("guess_center.txt");
    string line;
    while (getline(file, line))
    {
        istringstream iss(line);
        int cell;
        double x, y;
        iss >> cell >> x >> y;
        M8_center[cell] = make_pair(x, y);
    }
    file.close();
}

double Function(double *x, double *par)
{
    double res = 0;
    for (int i = 0; i <= DegX; i++)
    {
        for (int j = 0; j <= DegY; j++)
        {
            res += par[i * (DegX + 1) + j] * pow(x[0], i) * pow(x[1], j);
        }
    }

    return res;

    // return 1 * (par[1] * x[0] + par[2] + pow(x[0], 2) + par[3] * pow(x[0], 3) + par[6] * x[1] * x[1] * x[1] * x[0] + par[0] * pow(x[0], 4) + par[7] * pow(x[0], 5));
}

void FittingReconstruction()
{
    int counter = 0;
    TGraphErrors *G_coord_fitgrid = new TGraphErrors();
    std::ifstream file("guess_center.txt");
    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        double x, y;
        int N;
        if (!(iss >> N >> x >> y))
        {
            break;
        }
        if (!M8[N / n][N % n])
        {
            continue;
        }
        G_coord_fitgrid->SetPoint(counter, x, y);
        counter++;
    }
    TCanvas *c_coord_fitgrid = new TCanvas("c_coord_fitgrid", "c_coord_fitgrid", 800, 800);
    G_coord_fitgrid->SetMarkerStyle(20);
    G_coord_fitgrid->GetXaxis()->SetTitle("x [mm]");
    G_coord_fitgrid->GetYaxis()->SetTitle("y [mm]");
    G_coord_fitgrid->Draw("AP");
    c_coord_fitgrid->Write();

    TGraph2D *G2D_X = new TGraph2D();
    TGraph2D *G2D_Y = new TGraph2D();

    counter = 0;
    int rho = 2.;
    for (int i = 0; i < n * n; ++i)
    {
        if (M8[i / n][i % n])
        {
            double X_real = i % 8 * rho - rho * (8 - 1) / 2;
            double Y_real = i / 8 * rho - rho * (8 - 1) / 2;

            double x, y;
            G_coord_fitgrid->GetPoint(counter, x, y);
            G2D_X->SetPoint(counter, x, y, X_real);
            G2D_Y->SetPoint(counter, y, x, Y_real);
            counter++;
        }
    }

    // fitting X
    TCanvas *c_fit_X = new TCanvas("c_fit_X", "c_fit_X", 800, 800);
    G2D_X->SetMarkerStyle(20);
    f_X = new TF2("f_X", Function, G2D_X->GetXmin(), G2D_X->GetXmax(), G2D_X->GetYmin(), G2D_X->GetYmax(), (DegX + 1) * (DegY + 1));
    G2D_X->Fit("f_X", "MULTITHREAD");
    f_X->Draw("SURF2");
    G2D_X->Draw("AP SAME");
    c_fit_X->Write();

    // fitting Y
    TCanvas *c_fit_Y = new TCanvas("c_fit_Y", "c_fit_Y", 800, 800);
    G2D_Y->SetMarkerStyle(20);

    f_Y = new TF2("f_Y", Function, G2D_Y->GetXmin(), G2D_Y->GetXmax(), G2D_Y->GetYmin(), G2D_Y->GetYmax(), (DegX + 1) * (DegY + 1));
    G2D_Y->Fit("f_Y", "MULTITHREAD");
    f_Y->Draw("SURF2");
    G2D_Y->Draw("AP SAME");
    c_fit_Y->Write();

    // Plotting residus from fit X and Y
    TCanvas *c_residus_XY = new TCanvas("c_residus_XY", "c_residus_XY", 800, 800);
    TGraph *G_residus_XY = new TGraph();
    TGraph *G_real_XY = new TGraph();
    for (int i = 0; i < G2D_X->GetN(); ++i)
    {
        double x = G2D_X->GetZ()[i];
        double y = G2D_Y->GetZ()[i];
        double x_fit = f_X->Eval(G2D_X->GetX()[i], G2D_X->GetY()[i]);
        double y_fit = f_Y->Eval(G2D_Y->GetX()[i], G2D_Y->GetY()[i]);

        G_residus_XY->SetPoint(i, x_fit, y_fit);
        G_real_XY->SetPoint(i, x, y);
    }

    G_residus_XY->SetMarkerStyle(20);
    G_residus_XY->SetMarkerColor(kRed);
    G_residus_XY->Draw("AP");
    G_real_XY->SetMarkerStyle(20);
    G_real_XY->Draw("P SAME");
    c_residus_XY->Write();

    // Reconstruction of the grid
    TCanvas *c_reconstruction = new TCanvas("c_reconstruction_from_guess", "c_reconstruction_from_guess", 800, 800);
    H_reconstruction = new TH2D("H_reconstruction", "H_reconstruction", 160, -4, 4, 160, -4, 4);
    for (int i = 0; i < 1000000; ++i)
    {
        double x, y;
        H_precorrrected->GetRandom2(x, y);
        double x_fit = f_X->Eval(x, y);
        double y_fit = f_Y->Eval(y, x);
        H_reconstruction->Fill(x_fit, y_fit);
    }
    H_reconstruction->Draw("COLZ");
    c_reconstruction->Write();
}

///////////////// FIITING DEFORMED IMAGE /////////////////////////////
double MyFullFittedFunction2D(double *x, double *par)
{
    double yx[2] = {x[1], x[0]};
    double parX[(DegX + 1) * (DegY + 1)];
    double parY[(DegX + 1) * (DegY + 1)];

    for (int i = 0; i < (DegX + 1) * (DegY + 1); i++)
    {
        parX[i] = par[8 + n * n + i];
        parY[i] = par[8 * n * n + (DegX + 1) * (DegY + 1) + i];
    }

    double X = Function(x, parX);
    double Y = Function(yx, parY);

    double sigma_x = par[0];
    double sigma_y = par[1];
    double A_g = par[2];
    double mu_gx = par[3];
    double mu_gy = par[4];
    double sigma_gx = par[5];
    double sigma_gy = par[6];
    double bkg = par[7];
    double A[n][n];

    double result_sum = 0.0;
    double grid_sum = 0.0;

    double r = sqrt(x[0] * x[0] + x[1] * x[1]);

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {

            // if (!M8[i][j])
            //     continue;

            A[i][j] = par[8 + i * n + j];

            double x_center = i  * rho - rho * (n - 1) / 2;
            double y_center = j  * rho - rho * (n - 1) / 2;

            // if (abs(sqrt(x_center * x_center + y_center * y_center) - R) > 2 * rho)
            // {
            //     continue;
            // }

            double erf_x_pos = erf((X - (x_center - l2)) / (sqrt(2) * sigma_x));
            double erf_x_neg = erf((X - (x_center + l2)) / (sqrt(2) * sigma_x));
            double erf_y_pos = erf((Y - (y_center - l2)) / (sqrt(2) * sigma_y));
            double erf_y_neg = erf((Y - (y_center + l2)) / (sqrt(2) * sigma_y));

            // result_sum += (erf_x_pos - erf_x_neg) * (erf_y_pos - erf_y_neg) ;
            grid_sum += A[0][0] * (erf_x_pos - erf_x_neg) * (erf_y_pos - erf_y_neg);
        }
    }

    // RAW EXPRESSION  OF 2D GAUSSIAN
    double gauss = 0; // A_g * exp(-0.5 * ((x[0] - mu_gx) * (x[0] - mu_gx) / (sigma_gx * sigma_gx) + (x[1] - mu_gy) * (x[1] - mu_gy) / (sigma_gy * sigma_gy)));

    // result_sum *= gauss;

    double result = result_sum + grid_sum;

    // EDGE
    // double R = sqrt(X * X + Y * Y);
    // if (R >= cercle)
    // {
    //     double theta = atan2(Y, X);
    //     double sigma_r = 1 / sqrt(2) * sqrt(pow(cos(theta) * sigma_x, 2) + pow(sin(theta) * sigma_y, 2));
    //     double edge = 0.5 * erfc(((R - cercle) / (sqrt(2) * sigma_r)));
    //     result *= edge;
    // }

    // if (result > MAXI)
    // {
    //     return MAXI;
    // }

    return result;
}

void FullFunctionToMinimize2D()
{
    double chi2 = 0.0;
    int N = n * n;
    FittedFunction = new TF2("FittedFunction", MyFullFittedFunction2D, -0.25, 0.18, -0.17, 0.24, N + 8 + 2 * (DegX + 1) * (DegY + 1));
    FittedFunction->SetNpx(500);
    FittedFunction->SetNpy(500);

    // for (int i = 0; i < H_precorrrected->GetNbinsX(); i++)
    // {
    //     for (int j = 0; j < H_precorrrected->GetNbinsY(); j++)
    //     {
    //         if (H_precorrrected->GetBinContent(i, j) > MAXI)
    //         {
    //             H_precorrrected->SetBinContent(i, j, MAXI);
    //         }
    //     }
    // }

    // sigma x
    FittedFunction->SetParLimits(0, 0., 2);
    FittedFunction->SetParameter(0, 0.5);
    // FittedFunction->FixParameter(1, 0.0);

    // sigma y
    FittedFunction->SetParLimits(1, 0., 2);
    FittedFunction->SetParameter(1, 0.5);
    // FittedFunction->FixParameter(2, 0.0);

    // amplitude gaus
    FittedFunction->SetParLimits(2, 10, 5000);
    FittedFunction->SetParameter(2, 100);
    FittedFunction->FixParameter(2, 0);

    // mu ga3s x
    FittedFunction->SetParLimits(3, -2, 2);
    FittedFunction->SetParameter(3, -0.2);
    FittedFunction->FixParameter(3, 0);

    // mu ga4s y
    FittedFunction->SetParLimits(4, -2, 2);
    FittedFunction->SetParameter(4, -0.2);
    FittedFunction->FixParameter(4, 0);

    // sigma5gaus x
    FittedFunction->SetParLimits(5, 0.3, 2);
    FittedFunction->SetParameter(5, 0.8);
    FittedFunction->FixParameter(5, 0);

    // sigma6gaus y
    FittedFunction->SetParLimits(6, 0.2, 2);
    FittedFunction->SetParameter(6, 0.6);
    FittedFunction->FixParameter(6, 0);

    // bkg
    FittedFunction->SetParLimits(7, 0., 20);
    FittedFunction->SetParameter(7, 2);
    FittedFunction->FixParameter(7, 0);

    for (int i = 0; i < N; ++i)
    {
        if (i ==0)
            {
                FittedFunction->SetParLimits(8 + i, 0, 200);
                FittedFunction->SetParameter(8 + i, 20);
            }
        else
            {
                FittedFunction->FixParameter(8 + i, 0);
            }
        // if (M8[i % n][i / n])
        // {
        //     // amp ij
            
            
            
        // }
    }

    for (int i = 0; i <= (DegX + 1) * (DegY + 1); i++)
    {
        FittedFunction->SetParameter(8 + n * n + i, f_X->GetParameter(i));
        FittedFunction->SetParameter(8 + n * n + (DegX + 1) * (DegY + 1) + i, f_Y->GetParameter(i));
    }

    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0);
    // H_precorrrected->Rebin2D(2, 2);
    // H_precorrrected->GetXaxis()->SetRangeUser(-0.4, 0.4);
    // H_precorrrected->GetYaxis()->SetRangeUser(-0.4, 0.4);

    H_precorrrected->SetBinContent(182, 178, 0);
    H_precorrrected->Fit(FittedFunction, "MULTITHREAD RN", "", -0.4, 0.4);

    cout << "chi2 = " << FittedFunction->GetChisquare() / FittedFunction->GetNDF() << endl;

    TCanvas *c = new TCanvas("c", "c", 800, 800);
    H_precorrrected->Draw("COLZ");
    FittedFunction->Draw("SAME");
    c->Write();

    //// writing center of grid cell in a txt file with cell number x an y
    ofstream file;
    file.open("out_centers.txt");
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            // double x_center = i * rho - rho * (n-1)/2-0.01;
            // double y_center = j * rho - rho * (n-1)/2+0.017;
            file << i * n + j << " " << FittedFunction->GetParameter(10 + i * n + j) << " " << FittedFunction->GetParameter(N + 10 + i * n + j) << endl;
        }
    }
    file.close();

    // TH2D* hf = (TH2D*)FittedFunction->GetHistogram("hf");
    // TCanvas *cprojx = new TCanvas("cprojx", "cprojx", 800, 800);
    // TH1D* Hprojx = H_precorrrected->ProjectionX("Hprojx");
    // Hprojx->Draw("HIST");
    // TH1D* hfx = hf->ProjectionX("hfx");
    // hfx->SetLineColor(kRed);
    // hfx->Draw("SAME");
    // cprojx->Write();

    // TCanvas *cprojy = new TCanvas("cprojy", "cprojy", 800, 800);
    // TH1D* Hprojy = H_precorrrected->ProjectionY("Hprojy");
    // Hprojy->Draw("HIST");
    // TH1D* hfy = hf->ProjectionY("hfy");
    // hfy->SetLineColor(kRed);
    // hfy->Draw("SAME");
    // cprojy->Write();

    TCanvas *c_center = new TCanvas("c_center", "c_center", 800, 800);
    H_precorrrected->Draw("COLZ");
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (!M8[i][j])
                continue;
            double x_center = M8_center[i * n + j].first;
            double y_center = M8_center[i * n + j].second;

            TGraph *c_center = new TGraph(1);
            c_center->SetPoint(0, FittedFunction->GetParameter(10 + i * n + j), FittedFunction->GetParameter(N + 10 + i * n + j));
            c_center->SetMarkerStyle(20);
            c_center->SetMarkerSize(2);
            c_center->SetMarkerColor(kRed);
            c_center->Draw("P SAME");
            TGraph *c_center_old = new TGraph(1);
            c_center_old->SetPoint(0, x_center, y_center);
            c_center_old->SetMarkerStyle(20);
            c_center_old->SetMarkerSize(2);
            c_center_old->Draw("P SAME");
        }
    }
    c_center->Write();

    TCanvas *c_corner = new TCanvas("c_corner", "c_corner", 800, 800);
    H_precorrrected->Draw("COLZ");
    l = FittedFunction->GetParameter(N + 8);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            // double x_center = i * rho - rho * (n-1)/2-0.01;
            // double y_center = j * rho - rho * (n-1)/2+0.017;

            TGraph *c_corner = new TGraph(4);
            c_corner->SetPoint(0, FittedFunction->GetParameter(10 + i * n + j) + l / 2, FittedFunction->GetParameter(N + 10 + i * n + j) + l / 2);
            c_corner->SetPoint(1, FittedFunction->GetParameter(10 + i * n + j) - l / 2, FittedFunction->GetParameter(N + 10 + i * n + j) + l / 2);
            c_corner->SetPoint(2, FittedFunction->GetParameter(10 + i * n + j) + l / 2, FittedFunction->GetParameter(N + 10 + i * n + j) - l / 2);
            c_corner->SetPoint(3, FittedFunction->GetParameter(10 + i * n + j) - l / 2, FittedFunction->GetParameter(N + 10 + i * n + j) - l / 2);

            c_corner->SetMarkerStyle(20);
            c_corner->SetMarkerSize(2);
            c_corner->Draw("P SAME");
        }
    }
    c_corner->Write();

    FittedFunction->Write();

    // Reconstruction of the grid
    TCanvas *c_reconstruction = new TCanvas("c_reconstruction", "c_reconstruction", 800, 800);
    H_reconstruction = new TH2D("H_reconstruction1", "H_reconstruction1", 160, -4, 4, 160, -4, 4);

    for (int i = 0; i <= (DegX + 1) * (DegY + 1); i++)
    {
        f_X->SetParameter(i, FittedFunction->GetParameter(8 + n * n + i));
        f_Y->SetParameter(i, FittedFunction->GetParameter(8 + n * n + (DegX + 1) * (DegY + 1) + i));
    }

    for (int i = 0; i < 1000000; ++i)
    {
        double x, y;
        H_precorrrected->GetRandom2(x, y);
        double x_fit = f_X->Eval(x, y);
        double y_fit = f_Y->Eval(y, x);
        H_reconstruction->Fill(x_fit, y_fit);
    }
    H_reconstruction->Draw("COLZ");
    c_reconstruction->Write();
}

#endif