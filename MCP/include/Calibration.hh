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

double rho_real = 2;
double l_real = 1.2;
double rho = 0.1;
double l = 0.08;
double limit_real = 7.5;
int n = 8;
double limit = limit_real * rho / rho_real;
int MAXI = 10;
int DegX = 4;
int DegY = 4;
TF2 *FittedFunction;
TF2 *FinalFunction;
TF2 *MeasurementFunction;
TF2 *f_X;
TF2 *f_Y;

pair<double, double> SIGMA_X;
pair<double, double> SIGMA_Y;
pair<double, double> SIGMA_BEAM_X;
pair<double, double> SIGMA_BEAM_Y;
pair<double, double> MEAN_BEAM_X;
pair<double, double> MEAN_BEAM_Y;
pair<double, double> AMPLITUDE_BEAM;

double init_par[5] = {0, 1, 1, 1, 1};
double data[5] = {0, 0, 0, 0, 0};

const int NparX = 3;
const int NparY = 3;
const int Nparam = (NparX + 1) * (NparY + 1);

double param_X[NparX];
double param_Y[NparY];

bool M8[8][8] = {
    {false, false, false, false, false, false, false, false},
    {false, false, true, true, true, true, false, false},
    {false, true, true, true, true, true, true, false},
    {false, true, true, true, true, true, true, false},
    {false, true, true, true, true, true, true, false},
    {false, true, true, true, true, true, true, false},
    {false, false, true, true, true, true, false, false},
    {false, false, false, false, false, false, false, false}};

double x_final_max = 4;
double x_final_min = -4;
double y_final_max = 4;
double y_final_min = -4;
bool M_final[8][8] = {
    {false, false, false, false, false, false, false, false},
    {false, false, false, false, false, false, false, false},
    {false, false, true, true, true, true, false, false},
    {false, false, true, true, true, true, false, false},
    {false, false, true, true, true, true, false, false},
    {false, false, true, true, true, true, false, false},
    {false, false, false, false, false, false, false, false},
    {false, false, false, false, false, false, false, false}};

double x_measurement_min = -2.5;
double x_measurement_max = 2.5;
double y_measurement_min = -2.5;
double y_measurement_max = 2.5;

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

Double_t fit_points(pair<double, double> x, Double_t *par)
{
    int i, j, ij;
    double res = 0;
    for (int i = 0; i <= NparX; i++)
        for (int j = 0; j <= NparY; j++)
        {
            ij = i * (NparY + 1) + j;
            res = res + par[ij] * pow(x.first, i) * pow(x.second, j);
        }

    return res;
}

void readfit()
{
    ifstream infile("fit_params.txt");
    // double xmax, xmin, ymax, ymin;
    // infile >> xmax >> xmin >> ymax >> ymin;

    std::string line;
    std::vector<std::vector<double>> data;
    while (std::getline(infile, line))
    {
        std::vector<double> values;
        std::istringstream iss(line);
        double value;
        while (iss >> value)
        {
            cout << value << endl;
            values.push_back(value);
        }
        data.push_back(values);
    }

    std::copy(data[0].begin(), data[0].end(), param_X);
    std::copy(data[1].begin(), data[1].end(), param_Y);

    infile.close();
}

pair<double, double> LogFunction(double *x)
{
    double Va = x[A];
    double Vb = x[B];
    double Vc = x[C];
    double Vd = x[D];

    double norm = 1 / log((Va * Vb * Vc * Vd) / pow(Va + Vb + Vc + Vd, 4));
    double X = -log(Vb * Vc / (Va * Vd)) * norm;
    double Y = -log(Vb * Va / (Vc * Vd)) * norm;

    // double X = (Vd + Vb - Va - Vc) / (Va + Vb + Vc + Vd);
    // double Y = (Vd - Vb + Va - Vc) / (Va + Vb + Vc + Vd);

    // cout << "X: " << X << " Y: " << Y << endl;
    return make_pair(X, Y);
}

double MyFittedFunction1D(double *x, double *par)
{
    double result = 0.0;
    double A[n];
    double sigma_x[n];
    double sigma_y[n];

    l = par[n + 1];

    for (int i = 0; i < n; ++i)
    {
        // A[i] = par[i];
        sigma_x[i] = par[i];
    }

    for (int i = 0; i < n; ++i)
    {
        double x_center = i * rho - rho * (n - 1) / 2 + 0.1;

        double erf_x_pos = erf((x[0] - x_center - par[n + 6] + l / 2) / (sqrt(2) * sigma_x[i]));
        double erf_x_neg = erf((x[0] - x_center - par[n + 6] - l / 2) / (sqrt(2) * sigma_x[i]));

        result += (erf_x_pos - erf_x_neg);
    }

    double grid = result;
    grid *= par[n + 2];

    // raw forula gauss
    double gauss = par[n + 3] * exp(-0.5 * ((x[0] - par[n + 4]) * (x[0] - par[n + 4]) / (par[n + 5] * par[n + 5])));
    result *= gauss;

    return par[n] + result + grid;
}

double FunctionToMinimize1x()
{
    double chi2 = 0.0;
    TH1D *H1D = H_precorrrected->ProjectionX("H1Dx");

    int N = n;
    TF1 *FittedFunction = new TF1("FittedFunction", MyFittedFunction1D, -8, 8, N + 7);

    FittedFunction->SetNpx(150);

    for (int i = 0; i < n; ++i)
    {
        FittedFunction->SetParLimits(i, 0., 0.2);
        FittedFunction->SetParameter(i, 0.1);
    }

    FittedFunction->SetParLimits(N, 400, 800);
    FittedFunction->SetParameter(N, 500);

    FittedFunction->SetParLimits(N + 1, 1.0, 1.4);
    FittedFunction->SetParameter(N + 1, 1.2);

    FittedFunction->SetParLimits(N + 2, 1, 5000);
    FittedFunction->SetParameter(N + 2, 10);

    FittedFunction->SetParLimits(N + 3, 0, 10000);
    FittedFunction->SetParameter(N + 3, 10);

    FittedFunction->SetParLimits(N + 4, -0.4, 0.4);
    FittedFunction->SetParameter(N + 4, 0);

    FittedFunction->SetParLimits(N + 5, 0, 2);
    FittedFunction->SetParameter(N + 5, 0.5);

    FittedFunction->SetParLimits(N + 6, 0., 0.);
    FittedFunction->SetParameter(N + 6, 0.);

    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);
    H1D->Fit(FittedFunction, "MULTITHREAD R", "", -2, 2);
    H1D->Write();
    chi2 = FittedFunction->GetChisquare() / FittedFunction->GetNDF();

    return chi2;
}

double FunctionToMinimize1y()
{
    double chi2 = 0.0;
    TH1D *H1D = H_precorrrected->ProjectionY("H1Dy");

    int N = n;
    TF1 *FittedFunction = new TF1("FittedFunction", MyFittedFunction1D, -8, 8, N + 7);

    FittedFunction->SetNpx(150);

    for (int i = 0; i < n; ++i)
    {
        FittedFunction->SetParLimits(i, 0., 0.2);
        FittedFunction->SetParameter(i, 0.1);
        // FittedFunction->FixParameter(i, 5);
    }

    FittedFunction->SetParLimits(N, 400, 8000);
    FittedFunction->SetParameter(N, 500);

    FittedFunction->SetParLimits(N + 1, 1.0, 1.4);
    FittedFunction->SetParameter(N + 1, 1.2);

    FittedFunction->SetParLimits(N + 2, 1, 5000);
    FittedFunction->SetParameter(N + 2, 10);

    FittedFunction->SetParLimits(N + 3, 0, 10000);
    FittedFunction->SetParameter(N + 3, 10);

    FittedFunction->SetParLimits(N + 4, -0.4, 0.4);
    FittedFunction->SetParameter(N + 4, 0);

    FittedFunction->SetParLimits(N + 5, 0, 2);
    FittedFunction->SetParameter(N + 5, 0.5);

    FittedFunction->SetParLimits(N + 6, 0., 0.);
    FittedFunction->SetParameter(N + 6, 0.);

    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);
    H1D->Fit(FittedFunction, "MULTITHREAD R", "", -2, 2);
    H1D->Write();
    chi2 = FittedFunction->GetChisquare() / FittedFunction->GetNDF();

    return chi2;
}

void MyMinimizer()
{
    cout << "Minimization 1D" << endl;
    FunctionToMinimize1x();
    FunctionToMinimize1y();
}

/////////////////////////////////////////////////////////////////////

double MyFittedFunction2(double *x, double *par)
{
    double A[n][n];
    double sigma_x = par[n * n];
    double sigma_y = par[n * n + 1];
    double A_g = par[n * n + 2];
    double mu_gx = par[n * n + 3];
    double mu_gy = par[n * n + 4];
    double sigma_gx = par[n * n + 5];
    double sigma_gy = par[n * n + 6];
    double bkg = par[n * n + 7];
    double l = par[n * n + 8];

    double result_sum = 0.0;
    double grid_sum = 0.0;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            A[i][j] = par[i * n + j];
        }
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            double x_center = i * rho - rho * (n - 1) / 2 + 0.1;
            double y_center = j * rho - rho * (n - 1) / 2 + 0.15;

            double erf_x_pos = erf((x[0] - x_center + l / 2) / (sqrt(2) * sigma_x));
            double erf_x_neg = erf((x[0] - x_center - l / 2) / (sqrt(2) * sigma_x));
            double erf_y_pos = erf((x[1] - y_center + l / 2) / (sqrt(2) * sigma_y));
            double erf_y_neg = erf((x[1] - y_center - l / 2) / (sqrt(2) * sigma_y));

            result_sum += (erf_x_pos - erf_x_neg) * (erf_y_pos - erf_y_neg);
            grid_sum += A[i][j] * (erf_x_pos - erf_x_neg) * (erf_y_pos - erf_y_neg);
        }
    }

    // RAW EXPRESSION  OF 2D GAUSSIAN
    double gauss = A_g * exp(-0.5 * ((x[0] - mu_gx) * (x[0] - mu_gx) / (sigma_gx * sigma_gx) + (x[1] - mu_gy) * (x[1] - mu_gy) / (sigma_gy * sigma_gy)));

    result_sum *= gauss;

    return bkg + result_sum + grid_sum;
}

double FunctionToMinimize2()
{
    double chi2 = 0.0;
    int N = n * n;
    TF2 *FittedFunction = new TF2("FittedFunction", MyFittedFunction2, -2, 2, -2, 2, N + 9);
    FittedFunction->SetNpx(50);
    FittedFunction->SetNpy(50);

    // amplitude bkg in MCP
    for (int i = 0; i < N; ++i)
    {
        // FittedFunction->SetParLimits(i, 0., 500);
        // FittedFunction->SetParameter(i, 10);
        FittedFunction->FixParameter(i, 0);
        //
    }

    // //-1 -1
    // FittedFunction->FixParameter(0, 10);

    // // -1 1
    // FittedFunction->FixParameter(1, 20);

    // // 1 -1
    // FittedFunction->FixParameter(2, 10);

    // // 1 1
    // FittedFunction->FixParameter(3, 10);

    // FittedFunction->SetParLimits(0, 0., 900);
    // FittedFunction->SetParameter(0, 10);

    // sigma x
    FittedFunction->SetParLimits(N, 0., 0.3);
    FittedFunction->SetParameter(N, 0.15);
    // FittedFunction->FixParameter(N, 0.137);

    // sigma y
    FittedFunction->SetParLimits(N + 1, 0., 0.3);
    FittedFunction->SetParameter(N + 1, 0.15);
    // FittedFunction->FixParameter(N+1, 0.1);

    // amplitude gaus
    FittedFunction->SetParLimits(N + 2, 10, 5000);
    FittedFunction->SetParameter(N + 2, 100);
    // FittedFunction->FixParameter(N+2, 85);

    // mu gaus x
    FittedFunction->SetParLimits(N + 3, -2, 2);
    FittedFunction->SetParameter(N + 3, -0.2);
    // FittedFunction->FixParameter(N+3, -0.2);

    // mu gaus y
    FittedFunction->SetParLimits(N + 4, -2, 2);
    FittedFunction->SetParameter(N + 4, -0.2);
    // FittedFunction->FixParameter(N+4, -0.2);

    // sigma gaus x
    FittedFunction->SetParLimits(N + 5, 0.3, 2);
    FittedFunction->SetParameter(N + 5, 0.8);
    // FittedFunction->FixParameter(N+5, 0.8);

    // sigma gaus y
    FittedFunction->SetParLimits(N + 6, 0.2, 2);
    FittedFunction->SetParameter(N + 6, 0.6);
    // FittedFunction->FixParameter(N+6, 0.6);

    // bkg
    FittedFunction->SetParLimits(N + 7, 0., 100);
    FittedFunction->SetParameter(N + 7, 10);
    // FittedFunction->FixParameter(N+7, 10);

    // l
    FittedFunction->SetParLimits(N + 8, 0.8, 1.4);
    FittedFunction->SetParameter(N + 8, 1.2);
    // FittedFunction->FixParameter(N+8, 1.2);

    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);
    // H_precorrrected->Rebin2D(2, 2);
    H_precorrrected->GetXaxis()->SetRangeUser(-2, 2);
    H_precorrrected->GetYaxis()->SetRangeUser(-2, 2);

    H_precorrrected->SetBinContent(182, 178, 0);
    H_precorrrected->Fit(FittedFunction, "MULTITHREAD RN");

    cout << "chi2 = " << FittedFunction->GetChisquare() / FittedFunction->GetNDF() << endl;

    TCanvas *c = new TCanvas("c", "c", 800, 800);
    FittedFunction->Draw("SURF");
    H_precorrrected->Draw("SURF SAME");

    c->Write();
    FittedFunction->Write();

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    c1->Divide(2, 2);

    // suqare 1
    c1->cd(1);
    TH2D *Hz1 = (TH2D *)H_precorrrected->Clone("Hz1");
    Hz1->GetXaxis()->SetRangeUser(-2, 0);
    Hz1->GetYaxis()->SetRangeUser(-2, 0);
    Hz1->Draw("LEGO");
    TF2 *FittedFunction1 = (TF2 *)FittedFunction->Clone("FittedFunction1");
    for (int i = 0; i < N; ++i)
    {
        FittedFunction1->SetParameter(i, FittedFunction->GetParameter(i));
    }
    FittedFunction1->SetRange(-2, -2, 0, 0);
    FittedFunction1->Draw("SURF SAME");

    // suqare 2
    c1->cd(2);
    TH2D *Hz2 = (TH2D *)H_precorrrected->Clone("Hz2");
    Hz2->GetXaxis()->SetRangeUser(0, 2);
    Hz2->GetYaxis()->SetRangeUser(-2, 0);
    Hz2->Draw("LEGO");
    TF2 *FittedFunction2 = (TF2 *)FittedFunction->Clone("FittedFunction2");
    for (int i = 0; i < N; ++i)
    {
        FittedFunction2->SetParameter(i, FittedFunction->GetParameter(i));
    }
    FittedFunction2->SetRange(0, -2, 2, 0);
    FittedFunction2->Draw("SURF SAME");

    // //suqare 3
    c1->cd(3);
    TH2D *Hz3 = (TH2D *)H_precorrrected->Clone("Hz3");
    Hz3->GetXaxis()->SetRangeUser(-2, 0);
    Hz3->GetYaxis()->SetRangeUser(0, 2);
    Hz3->Draw("LEGO");
    TF2 *FittedFunction3 = (TF2 *)FittedFunction->Clone("FittedFunction3");
    for (int i = 0; i < N; ++i)
    {
        FittedFunction3->SetParameter(i, FittedFunction->GetParameter(i));
    }
    FittedFunction3->SetRange(-2, 0, 0, 2);
    FittedFunction3->Draw("SURF SAME");

    // //suqare 4
    c1->cd(4);
    TH2D *Hz4 = (TH2D *)H_precorrrected->Clone("Hz4");
    Hz4->GetXaxis()->SetRangeUser(0, 2);
    Hz4->GetYaxis()->SetRangeUser(0, 2);
    Hz4->Draw("LEGO");
    TF2 *FittedFunction4 = (TF2 *)FittedFunction->Clone("FittedFunction4");
    for (int i = 0; i < N; ++i)
    {
        FittedFunction4->SetParameter(i, FittedFunction->GetParameter(i));
    }
    FittedFunction4->SetRange(0, 0, 2, 2);
    FittedFunction4->Draw("SURF SAME");
    c1->Write();

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 800);
    H_precorrrected->Draw("COLZ");
    TF2 *xyg = new TF2("xygaus", "xygaus", -2, 2, -2, 2);
    xyg->SetParameter(0, FittedFunction->GetParameter(6));
    xyg->SetParameter(1, FittedFunction->GetParameter(7));
    xyg->SetParameter(2, FittedFunction->GetParameter(9));
    xyg->SetParameter(3, FittedFunction->GetParameter(8));
    xyg->SetParameter(4, FittedFunction->GetParameter(10));
    xyg->Draw("CONT1 SAME");
    c2->Write();

    chi2 = FittedFunction->GetChisquare() / FittedFunction->GetNDF();

    return chi2;
}

void MyMinimizer2D()
{
    cout << "Minimization 2D" << endl;
    FunctionToMinimize2();
}

/////////////////////////////////////////////////////////////////////:

double MyFittedFunction2D(double *x, double *par)
{
    double result = 0.0;
    double A[n][n];
    double sigma_x[n][n];
    double sigma_y[n][n];

    l = par[3 * n * n + 1];

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            A[i][j] = par[i * n + j];
            sigma_x[i][j] = par[n * n + i * n + j];
            sigma_y[i][j] = par[2 * n * n + i * n + j];
        }
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            double x_center = i * rho - rho * (n - 1) / 2;
            double y_center = j * rho - rho * (n - 1) / 2;

            double erf_x_pos = erf((x[0] - x_center + l / 2) / (sqrt(2) * sigma_x[i][j]));
            double erf_x_neg = erf((x[0] - x_center - l / 2) / (sqrt(2) * sigma_x[i][j]));
            double erf_y_pos = erf((x[1] - y_center + l / 2) / (sqrt(2) * sigma_y[i][j]));
            double erf_y_neg = erf((x[1] - y_center - l / 2) / (sqrt(2) * sigma_y[i][j]));

            result += A[i][j] * (erf_x_pos - erf_x_neg) * (erf_y_pos - erf_y_neg);
        }
    }

    return par[3 * n * n] + result;
}

double Fitting()
{
    double chi2 = 0.0;

    int N = 3 * n * n;
    TF2 *FittedFunction = new TF2("FittedFunction", MyFittedFunction2D, -10, 10, -10, 10, N + 2);
    FittedFunction->SetNpx(H_precorrrected->GetNbinsX());
    FittedFunction->SetNpy(H_precorrrected->GetNbinsY());
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            FittedFunction->SetParLimits(i * n + j, 1, 500);
            FittedFunction->SetParameter(i * n + j, 1);
            FittedFunction->SetParLimits(n * n + i * n + j, 0., 0.2);
            FittedFunction->SetParameter(n * n + i * n + j, 0.1);
            FittedFunction->SetParLimits(2 * n * n + i * n + j, 0, 0.2);
            FittedFunction->SetParameter(2 * n * n + i * n + j, 0.1);
        }
    }

    FittedFunction->SetParLimits(N, 0, 100);
    FittedFunction->SetParameter(N, 1);

    FittedFunction->SetParLimits(N + 1, 1.2, 1.2);
    FittedFunction->SetParameter(N + 1, 1);

    cout << "Fitting" << endl;
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
    H_precorrrected->Rebin2D(2, 2);
    H_precorrrected->Fit(FittedFunction, "MULTITHREAD R", "", -(n - 1) * rho, (n - 1) * rho);

    FittedFunction->Write();
    H_precorrrected->Write();
    chi2 = FittedFunction->GetChisquare() / FittedFunction->GetNDF();

    TCanvas *c = new TCanvas("c", "c", 800, 800);
    TH1D *H1D = H_precorrrected->ProjectionY("H1D", 80, 89);
    c->cd();
    H1D->Draw("HIST");
    TH2D *FH = (TH2D *)FittedFunction->GetHistogram();
    TH1D *FH1D = FH->ProjectionY("FH1D", 80, 89);
    FH1D->Draw("SAME");
    c->Write();

    // // ouvrir fichier txt
    // ofstream file("chi2.txt", ios::app);
    // file << chi2 << " " << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << endl;
    // file.close();

    cout << "Chi2 : " << chi2 << endl;

    delete FittedFunction;

    return chi2;
}

/////////////////////////////////////////////////////////////////////
//////////// NEW 2D FULL FITTING ////////////////////////////////////
/////////////////////////////////////////////////////////////////////

///////////////// FIITING DEFORMED IMAGE /////////////////////////////
double MyFullFittedFunction2D(double *x, double *par)
{
    double A = par[0];
    double sigma_x = par[1];
    double sigma_y = par[2];
    double A_g = par[3];
    double mu_gx = par[4];
    double mu_gy = par[5];
    double sigma_gx = par[6];
    double sigma_gy = par[7];
    double bkg = par[8];
    double l2 = par[9] / 2;
    double beta_x[n][n];
    double beta_y[n][n];
    // double a[4][4] = { {par[2 * n*n + 10], par[2 * n*n + 11], par[2 * n*n + 12], par[2 * n*n + 13]},
    //                    {par[2 * n*n + 14], par[2 * n*n + 15], par[2 * n*n + 16], par[2 * n*n + 17]},
    //                    {par[2 * n*n + 18], par[2 * n*n + 19], par[2 * n*n + 20], par[2 * n*n + 21]},
    //                    {par[2 * n*n + 22], par[2 * n*n + 23], par[2 * n*n + 24], par[2 * n*n + 25]}};

    double result_sum = 0.0;
    double grid_sum = 0.0;

    double r = sqrt(x[0] * x[0] + x[1] * x[1]);

    if (r > 0.4)
    {
        return 0;
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {

            if (!M8[i][j])
                continue;

            if (abs(sqrt(M8_center[i * n + j].first * M8_center[i * n + j].first + M8_center[i * n + j].second * M8_center[i * n + j].second) - r) > 2 * rho)
            {
                continue;
            }

            beta_x[i][j] = par[10 + i * n + j];
            beta_y[i][j] = par[10 + n * n + i * n + j];

            double erf_x_pos = erf((x[0] - (beta_x[i][j] - l2)) / (sqrt(2) * sigma_x));
            double erf_x_neg = erf((x[0] - (beta_x[i][j] + l2)) / (sqrt(2) * sigma_x));
            double erf_y_pos = erf((x[1] - (beta_y[i][j] - l2)) / (sqrt(2) * sigma_y));
            double erf_y_neg = erf((x[1] - (beta_y[i][j] + l2)) / (sqrt(2) * sigma_y));

            // result_sum += (erf_x_pos - erf_x_neg) * (erf_y_pos - erf_y_neg) ;
            grid_sum += A * (erf_x_pos - erf_x_neg) * (erf_y_pos - erf_y_neg);
        }
    }

    // RAW EXPRESSION  OF 2D GAUSSIAN
    double gauss = 0; // A_g * exp(-0.5 * ((x[0] - mu_gx) * (x[0] - mu_gx) / (sigma_gx * sigma_gx) + (x[1] - mu_gy) * (x[1] - mu_gy) / (sigma_gy * sigma_gy)));

    // result_sum *= gauss;

    double result = result_sum + grid_sum;

    // EDGE
    // if (n > 5 && r > 0.25)
    // {

    //     double edge_line;
    //     double R;
    //     if (x[1] > x[0] && x[1] > -x[0])
    //     {
    //         edge_line = a[0][0] + a[0][1] * x[0] + a[0][2] * pow(x[0], 2); // + a[0][3] * pow(x[0], 3);
    //         // R = sqrt(edge_line * edge_line + x[0] * x[0]);
    //         if (x[1] > edge_line)
    //         {return 0;}
    //     }
    //     else if (x[1] > x[0] && x[1] < -x[0])
    //     {
    //         edge_line = a[1][0] + a[1][1] * x[1] + a[1][2] * pow(x[1], 2); // + a[1][3] * pow(x[1], 3);
    //         // R = sqrt(x[1] * x[1] + edge_line * edge_line);
    //         if (x[0] < edge_line)
    //         {return 0;}
    //     }
    //     else if (x[1] < x[0] && x[1] < -x[0])
    //     {
    //         edge_line = a[2][0] + a[2][1] * x[0] + a[2][2] * pow(x[0], 2); // + a[2][3] * pow(x[0], 3);
    //         // R = sqrt(edge_line * edge_line + x[0] * x[0]);
    //         if (x[1] < edge_line)
    //         {return 0;}
    //     }
    //     else
    //     {
    //         edge_line = a[3][0] + a[3][1] * x[1] + a[3][2] * pow(x[1], 2); // + a[3][3] * pow(x[1], 3);
    //         if (x[0] > edge_line)
    //         {return 0;}
    //     }

    //     // double theta = atan(x[1] / x[0]);
    //     // double sigma_r = 1 / sqrt(2) * sqrt(pow(cos(theta) * sigma_x, 2) + pow(sin(theta) * sigma_y, 2));
    //     // double edge = 0.5 * erfc(((r - R) / (sqrt(2) * sigma_r)));

    //     // result *= edge;
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
    FittedFunction = new TF2("FittedFunction", MyFullFittedFunction2D, -0.4, 0.4, -0.4, 0.4, 2 * N + 10);
    FittedFunction->SetNpx(500);
    FittedFunction->SetNpy(500);

    for (int i = 0; i < H_precorrrected->GetNbinsX(); i++)
    {
        for (int j = 0; j < H_precorrrected->GetNbinsY(); j++)
        {
            if (H_precorrrected->GetBinContent(i, j) > MAXI)
            {
                H_precorrrected->SetBinContent(i, j, MAXI);
            }
        }
    }

    // Amplitude
    FittedFunction->SetParLimits(0, 0., 100);
    FittedFunction->SetParameter(0, 10);

    // sigma x
    FittedFunction->SetParLimits(1, 0., 0.02);
    FittedFunction->SetParameter(1, 0.005);
    // FittedFunction->FixParameter(1, 0.0);

    // sigma y
    FittedFunction->SetParLimits(2, 0., 0.02);
    FittedFunction->SetParameter(2, 0.005);
    // FittedFunction->FixParameter(2, 0.0);

    // amplitude gaus
    FittedFunction->SetParLimits(3, 10, 5000);
    FittedFunction->SetParameter(3, 100);
    FittedFunction->FixParameter(3, 0);

    // mu gaus x
    FittedFunction->SetParLimits(4, -2, 2);
    FittedFunction->SetParameter(4, -0.2);
    FittedFunction->FixParameter(4, 0);

    // mu gaus y
    FittedFunction->SetParLimits(5, -2, 2);
    FittedFunction->SetParameter(5, -0.2);
    FittedFunction->FixParameter(5, 0);

    // sigma gaus x
    FittedFunction->SetParLimits(6, 0.3, 2);
    FittedFunction->SetParameter(6, 0.8);
    FittedFunction->FixParameter(6, 0);

    // sigma gaus y
    FittedFunction->SetParLimits(7, 0.2, 2);
    FittedFunction->SetParameter(7, 0.6);
    FittedFunction->FixParameter(7, 0);

    // bkg
    // FittedFunction->SetParLimits(8, 0., 20);
    // FittedFunction->SetParameter(8, 2);
    FittedFunction->FixParameter(8, 0);

    // l
    FittedFunction->SetParLimits(9, 0.05, 0.1);
    FittedFunction->SetParameter(9, 0.08);
    // FittedFunction->FixParameter(9, 0.08);

    for (int i = 0; i < N; ++i)
    {
        // amplitude bkg in MCP
        // FittedFunction->SetParLimits(i, 0., 500);
        // FittedFunction->SetParameter(i, 10);

        if (!M8[i / n][i % n])
        {
            // FittedFunction->FixParameter(i, 0);
            // beta_x
            FittedFunction->FixParameter(10 + i, 0.);
            // beta_y
            FittedFunction->FixParameter(N + 10 + i, 0.);
        }
        else
        {
            // FittedFunction->FixParameter(i, 0);
            // double x_center = i / n * rho - rho * (n-1)/2-0.01;
            // double y_center = i % n * rho - rho * (n-1)/2+0.017;

            // beta_x
            FittedFunction->SetParLimits(10 + i, M8_center[i].first - 0.03, M8_center[i].first + 0.03);
            FittedFunction->SetParameter(10 + i, M8_center[i].first);

            // beta_y
            FittedFunction->SetParLimits(N + 10 + i, M8_center[i].second - 0.03, M8_center[i].second + 0.03);
            FittedFunction->SetParameter(N + 10 + i, M8_center[i].second);
        }
    }

    // // // edge a[line][coef]
    // if (n > 5)
    // {
    //     for (int line = 0; line < 4; line++)
    //     {
    //         // FittedFunction->SetParLimits(3 * N + 9 + 4*line, -0.1, 0.1);
    //         double guess = 0.0;
    //         if (line == 0 || line == 3) guess = 0.3;
    //         else guess = -0.3;
    //         FittedFunction->FixParameter(2 * N + 10 + 0 + 4*line, 0.5);

    //         // FittedFunction->SetParLimits(3 * N + 10 + 4*line, -0.1, 0.1);
    //         FittedFunction->FixParameter(2 * N + 10 + 1 + 4*line, 0);

    //         // FittedFunction->SetParLimits(3 * N + 11 + 4*line, -0.1, 0.1);
    //         FittedFunction->FixParameter(2 * N + 10 + 2 + 4*line, 0);

    //         // FittedFunction->SetParLimits(3 * N + 12 + 4*line, -0.1, 0.1);
    //         FittedFunction->FixParameter(2 * N + 10 + 3 + 4*line, 0);
    //     }
    // }
    // else
    // {
    //     FittedFunction->FixParameter(3*N+9, 0);
    //     FittedFunction->FixParameter(3*N+10, 0);
    //     FittedFunction->FixParameter(3*N+11, 0);
    //     FittedFunction->FixParameter(3*N+12, 0);
    // }

    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0);
    // H_precorrrected->Rebin2D(2, 2);
    H_precorrrected->GetXaxis()->SetRangeUser(-0.4, 0.4);
    H_precorrrected->GetYaxis()->SetRangeUser(-0.4, 0.4);

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
}

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
    if (FittedFunction != nullptr)
    {
        for (int i = 0; i < n * n; ++i)
        {
            if (M8[i / n][i % n])
            {
                G_coord_fitgrid->SetPoint(counter, FittedFunction->GetParameter(10 + i), FittedFunction->GetParameter(n * n + 10 + i));
                G_coord_fitgrid->SetPointError(counter, FittedFunction->GetParError(10 + i), FittedFunction->GetParError(n * n + 10 + i));
                counter++;
            }
        }
        TCanvas *c_coord_fitgrid = new TCanvas("c_coord_fitgrid", "c_coord_fitgrid", 800, 800);
        G_coord_fitgrid->SetMarkerStyle(20);
        G_coord_fitgrid->GetXaxis()->SetTitle("x [mm]");
        G_coord_fitgrid->GetYaxis()->SetTitle("y [mm]");
        G_coord_fitgrid->Draw("AP");
        c_coord_fitgrid->Write();
    }
    else
    {
        std::ifstream file("guess_center.txt");
        std::string line;
        int counter = 0;
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
            // cout << N << " " << x << " " << y << endl;
            counter++;
        }
        TCanvas *c_coord_fitgrid = new TCanvas("c_coord_fitgrid", "c_coord_fitgrid", 800, 800);
        G_coord_fitgrid->SetMarkerStyle(20);
        G_coord_fitgrid->GetXaxis()->SetTitle("x [mm]");
        G_coord_fitgrid->GetYaxis()->SetTitle("y [mm]");
        G_coord_fitgrid->Draw("AP");
        c_coord_fitgrid->Write();
    }

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
    TCanvas *c_reconstruction = new TCanvas("c_reconstruction", "c_reconstruction", 800, 800);
    H_reconstruction = new TH2D("H_reconstruction", "H_reconstruction", 80, -4, 4, 80, -4, 4);
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

///////////////// Verify grid /////////////////////////////////////////
double MyFinalFittedFunction2D(double *x, double *par)
{
    double A = par[0];
    double sigma_x = par[1];
    double sigma_y = par[2];
    double A_g = par[3];
    double mu_gx = par[4];
    double mu_gy = par[5];
    double sigma_gx = par[6];
    double sigma_gy = par[7];
    double bkg = par[8];
    double l2 = par[9] / 2;
    double beta_x[n][n];
    double beta_y[n][n];

    double result_sum = 0.0;
    double grid_sum = 0.0;

    double r = sqrt(x[0] * x[0] + x[1] * x[1]);

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {

            if (!M_final[i][j])
                continue;

            beta_x[i][j] = par[10 + i * n + j];
            beta_y[i][j] = par[10 + n * n + i * n + j];

            double erf_x_pos = erf((x[0] - (beta_x[i][j] - l2)) / (sqrt(2) * sigma_x));
            double erf_x_neg = erf((x[0] - (beta_x[i][j] + l2)) / (sqrt(2) * sigma_x));
            double erf_y_pos = erf((x[1] - (beta_y[i][j] - l2)) / (sqrt(2) * sigma_y));
            double erf_y_neg = erf((x[1] - (beta_y[i][j] + l2)) / (sqrt(2) * sigma_y));

            // result_sum += (erf_x_pos - erf_x_neg) * (erf_y_pos - erf_y_neg) ;
            grid_sum += A * (erf_x_pos - erf_x_neg) * (erf_y_pos - erf_y_neg);
        }
    }

    // RAW EXPRESSION  OF 2D GAUSSIAN
    double gauss = 0; // A_g * exp(-0.5 * ((x[0] - mu_gx) * (x[0] - mu_gx) / (sigma_gx * sigma_gx) + (x[1] - mu_gy) * (x[1] - mu_gy) / (sigma_gy * sigma_gy)));

    // result_sum *= gauss;

    double result = result_sum + grid_sum;

    return result;
}

void FinalFunctionToMinimize2D()
{
    double chi2 = 0.0;
    int N = n * n;
    rho = 2.;
    FinalFunction = new TF2("FinalFunction", MyFinalFittedFunction2D, x_final_min, x_final_max, x_final_min, x_final_max, 2 * N + 10);
    FinalFunction->SetNpx(75);
    FinalFunction->SetNpy(75);

    for (int i = 0; i < H_reconstruction->GetNbinsX(); i++)
    {
        for (int j = 0; j < H_reconstruction->GetNbinsY(); j++)
        {
            if (H_reconstruction->GetBinContent(i, j) > 100)
            {
                H_reconstruction->SetBinContent(i, j, 100);
            }
        }
    }

    // Amplitude
    FinalFunction->SetParLimits(0, 0., 100);
    FinalFunction->SetParameter(0, 30);

    // sigma x
    FinalFunction->SetParLimits(1, 0., 0.5);
    FinalFunction->SetParameter(1, 0.3);
    // FinalFunction->FixParameter(1, 0.0);

    // sigma y
    FinalFunction->SetParLimits(2, 0., 0.5);
    FinalFunction->SetParameter(2, 0.3);
    // FinalFunction->FixParameter(2, 0.0);

    // amplitude gaus
    FinalFunction->SetParLimits(3, 10, 5000);
    FinalFunction->SetParameter(3, 100);
    FinalFunction->FixParameter(3, 0);

    // mu gaus x
    FinalFunction->SetParLimits(4, -2, 2);
    FinalFunction->SetParameter(4, -0.2);
    FinalFunction->FixParameter(4, 0);

    // mu gaus y
    FinalFunction->SetParLimits(5, -2, 2);
    FinalFunction->SetParameter(5, -0.2);
    FinalFunction->FixParameter(5, 0);

    // sigma gaus x
    FinalFunction->SetParLimits(6, 0.3, 2);
    FinalFunction->SetParameter(6, 0.8);
    FinalFunction->FixParameter(6, 0);

    // sigma gaus y
    FinalFunction->SetParLimits(7, 0.2, 2);
    FinalFunction->SetParameter(7, 0.6);
    FinalFunction->FixParameter(7, 0);

    // bkg
    FinalFunction->SetParLimits(8, 0., 20);
    FinalFunction->SetParameter(8, 2);
    // FinalFunction->FixParameter(8, 0);

    // l
    FinalFunction->SetParLimits(9, 1., 1.4);
    FinalFunction->SetParameter(9, 1.2);
    // FinalFunction->FixParameter(9, 1.2);

    for (int i = 0; i < N; ++i)
    {
        // amplitude bkg in MCP
        // FinalFunction->SetParLimits(i, 0., 500);
        // FinalFunction->SetParameter(i, 10);

        if (!M_final[i / n][i % n])
        {
            // FinalFunction->FixParameter(i, 0);
            // beta_x
            FinalFunction->FixParameter(10 + i, 0.);
            // beta_y
            FinalFunction->FixParameter(N + 10 + i, 0.);
        }
        else
        {
            // FinalFunction->FixParameter(i, 0);
            double x_center = i / n * rho - rho * (n - 1) / 2;
            double y_center = i % n * rho - rho * (n - 1) / 2;

            // beta_x
            FinalFunction->SetParLimits(10 + i, x_center - 0.4, x_center + 0.4);
            FinalFunction->SetParameter(10 + i, x_center);

            // beta_y
            FinalFunction->SetParLimits(N + 10 + i, y_center - 0.4, y_center + 0.4);
            FinalFunction->SetParameter(N + 10 + i, y_center);
        }
    }

    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0);
    // H_precorrrected->Rebin2D(2, 2);
    H_reconstruction->Write();
    // H_reconstruction->GetXaxis()->SetRangeUser(-5, 5);
    // H_reconstruction->GetYaxis()->SetRangeUser(-5, 5);

    // H_reconstruction->SetBinContent(182, 178, 0);
    H_reconstruction->Fit(FinalFunction, "MULTITHREAD RN", "", x_final_min, x_final_max);

    cout << "chi2 = " << FinalFunction->GetChisquare() / FinalFunction->GetNDF() << endl;

    TCanvas *c = new TCanvas("c", "c", 800, 800);
    H_reconstruction->Draw("COLZ");
    FinalFunction->Draw("SAME");
    c->Write();

    /////// SAVING PARAMETER ////////
    SIGMA_X = make_pair(FinalFunction->GetParameter(1), FinalFunction->GetParError(1));
    SIGMA_Y = make_pair(FinalFunction->GetParameter(2), FinalFunction->GetParError(2));
    /////////////////////////////////
    // TH2D* hf = (TH2D*)FinalFunction->GetHistogram("hf");
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
    H_reconstruction->Draw("COLZ");
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (!M8[i][j])
                continue;
            double x_center = i * rho - rho * (n - 1) / 2;
            double y_center = j * rho - rho * (n - 1) / 2;

            TGraph *c_center = new TGraph(1);
            c_center->SetPoint(0, FinalFunction->GetParameter(10 + i * n + j), FinalFunction->GetParameter(N + 10 + i * n + j));
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
    H_reconstruction->Draw("COLZ");
    l = FinalFunction->GetParameter(9);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            // double x_center = i * rho - rho * (n-1)/2-0.01;
            // double y_center = j * rho - rho * (n-1)/2+0.017;

            TGraph *c_corner = new TGraph(4);
            c_corner->SetPoint(0, FinalFunction->GetParameter(10 + i * n + j) + l / 2, FinalFunction->GetParameter(N + 10 + i * n + j) + l / 2);
            c_corner->SetPoint(1, FinalFunction->GetParameter(10 + i * n + j) - l / 2, FinalFunction->GetParameter(N + 10 + i * n + j) + l / 2);
            c_corner->SetPoint(2, FinalFunction->GetParameter(10 + i * n + j) + l / 2, FinalFunction->GetParameter(N + 10 + i * n + j) - l / 2);
            c_corner->SetPoint(3, FinalFunction->GetParameter(10 + i * n + j) - l / 2, FinalFunction->GetParameter(N + 10 + i * n + j) - l / 2);

            c_corner->SetMarkerStyle(20);
            c_corner->SetMarkerSize(2);
            c_corner->Draw("P SAME");
        }
    }
    c_corner->Write();

    FinalFunction->Write();
}

///////////////// 2D beam measurement /////////////////////////////////////////

double MyGaussian(double *x, double *par)
{
    double A_g = par[0];
    double mu_gx = par[1];
    double sigma_gx = par[2];
    double mu_gy = par[3];
    double sigma_gy = par[4];

    // RAW EXPRESSION  OF 2D GAUSSIAN (R)
    double gauss = A_g * exp(-0.5 * ((x[0] - mu_gx) * (x[0] - mu_gx) / (sigma_gx * sigma_gx) + (x[1] - mu_gy) * (x[1] - mu_gy) / (sigma_gy * sigma_gy)));

    return gauss;
}

double MeasurementFittedFunction2D(double *x, double *par)
{

    // double A = par[0];
    double sigma_x = par[1];
    double sigma_y = par[2];
    double A_g = par[3];
    double mu_gx = par[4];
    double mu_gy = par[5];
    double sigma_gx = par[6];
    double sigma_gy = par[7];
    double bkg = par[8];
    double l2 = par[9] / 2;
    double beta_x[n][n];
    double beta_y[n][n];
    double A[n][n];
    double a[6] = {par[3*n*n+10], par[3*n*n+11], par[3*n*n+21], par[3*n*n+13], par[3*n*n+14], par[3*n*n+15]};

    double result_sum = 0.0;
    double grid_sum = 0.0;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {

            if (!M_measurement[i][j])
                continue;

            A[i][j] = par[2 * n * n + 10 + i * n + j];
            beta_x[i][j] = par[10 + i * n + j];
            beta_y[i][j] = par[10 + n * n + i * n + j];

            double erf_x_pos = erf((x[0] - (beta_x[i][j] - l2)) / (sqrt(2) * sigma_x));
            double erf_x_neg = erf((x[0] - (beta_x[i][j] + l2)) / (sqrt(2) * sigma_x));
            double erf_y_pos = erf((x[1] - (beta_y[i][j] - l2)) / (sqrt(2) * sigma_y));
            double erf_y_neg = erf((x[1] - (beta_y[i][j] + l2)) / (sqrt(2) * sigma_y));

            result_sum += (erf_x_pos - erf_x_neg) * (erf_y_pos - erf_y_neg);
            grid_sum += A[i][j] * (erf_x_pos - erf_x_neg) * (erf_y_pos - erf_y_neg);
        }
    }

    // RAW EXPRESSION  OF 2D GAUSSIAN
    double gauss = A_g * exp(-0.5 * ((x[0] - mu_gx) * (x[0] - mu_gx) / (sigma_gx * sigma_gx) + (x[1] - mu_gy) * (x[1] - mu_gy) / (sigma_gy * sigma_gy)));

    // RAW EXPRESSION  OF 2D GAUSSIAN (R)
    // double GaussianParameters[9] = {A_g, mu_gx, mu_gy, a[0], a[1], a[2], a[3], a[4], a[5]};
    
    result_sum *= gauss;

    double result = result_sum + grid_sum;

    return result;
}

void MeasurementFunctionToMinimize2D()
{

    double chi2 = 0.0;
    int N = n * n;
    MeasurementFunction = new TF2("MeasurementFunction", MeasurementFittedFunction2D, x_measurement_min, x_measurement_max, x_measurement_min, x_measurement_max, 3 * N + 10);
    MeasurementFunction->SetNpx(75);
    MeasurementFunction->SetNpy(75);

    // Amplitude
    MeasurementFunction->SetParLimits(0, 0, 1000);
    MeasurementFunction->SetParameter(0, 50);
    MeasurementFunction->FixParameter(0, 0);

    // sigma x
    MeasurementFunction->SetParLimits(1, 0.1, 0.4);
    // MeasurementFunction->FixParameter(1, FinalFunction->GetParameter(1));

    // sigma y
    MeasurementFunction->SetParLimits(2, 0.1, 0.4);
    // MeasurementFunction->FixParameter(2, FinalFunction->GetParameter(2));

    // amplitude gaus
    MeasurementFunction->SetParLimits(3, 10, 1000);
    MeasurementFunction->SetParameter(3, 100);
    //MeasurementFunction->FixParameter(3, 1000);

    // mu gaus x
    MeasurementFunction->SetParLimits(4, -1, 0.25);
    MeasurementFunction->SetParameter(4, -0.2);
    // MeasurementFunction->FixParameter(4,-0.2);

    // mu gaus y
    MeasurementFunction->SetParLimits(5, -1, 0.);
    MeasurementFunction->SetParameter(5, -0.1);
    // MeasurementFunction->FixParameter(5, -0.1);

    // sigma gaus x
    MeasurementFunction->SetParLimits(6, 0.3, 1.5);
    MeasurementFunction->SetParameter(6, 0.5);
    // MeasurementFunction->FixParameter(6, 1);

    // sigma gaus y
    MeasurementFunction->SetParLimits(7, 0.3, 1.5);
    MeasurementFunction->SetParameter(7, 0.5);
    // MeasurementFunction->FixParameter(7, 1);

    // bkg
    MeasurementFunction->SetParLimits(8, 0, 100);
    MeasurementFunction->SetParameter(8, 10);
    // MeasurementFunction->FixParameter(8, 0);

    // l
    MeasurementFunction->FixParameter(9, FinalFunction->GetParameter(9));

    for (int i = 0; i < N; ++i)
    {

        if (!M_measurement[i / n][i % n])
        {
            MeasurementFunction->FixParameter(10 + i, 0.);
            MeasurementFunction->FixParameter(N + 10 + i, 0.);
            MeasurementFunction->FixParameter(2*N + 10 + i, 0.);
        }
        else
        {
            // MeasurementFunction->FixParameter(i, 0);
            double x_center = i / n * rho - rho * (n - 1) / 2;
            double y_center = i % n * rho - rho * (n - 1) / 2;

            // beta_x
            MeasurementFunction->FixParameter(10 + i, x_center);

            // beta_y
            MeasurementFunction->FixParameter(N + 10 + i, y_center);

            // amp
            MeasurementFunction->SetParLimits(2 * N + 10 + i, 0, 200);
            // MeasurementFunction->FixParameter(2*N + 10 + i, 0.);
        }
    }    


    // MeasurementFunction->SetParLimits(3 * N + 10, 0, 20);
    // // MeasurementFunction->SetParameter(3 * N + 10, 0.3);

    // MeasurementFunction->SetParLimits(3 * N + 11, 0, 20);

    // MeasurementFunction->FixParameter(3 * N + 12, 0);

    // MeasurementFunction->FixParameter(3 * N + 13, 0);

    // MeasurementFunction->FixParameter(3 * N + 14, 0);

    // MeasurementFunction->FixParameter(3 * N + 15, 0);


    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0);
    H_measurement_reconstruted->Rebin2D(2, 2);

    H_measurement_reconstruted->GetXaxis()->SetRangeUser(x_measurement_min, x_measurement_max);
    H_measurement_reconstruted->GetYaxis()->SetRangeUser(x_measurement_min, x_measurement_max);

    H_measurement_reconstruted->SetBinContent(57, 51, 0);
    H_measurement_reconstruted->Fit(MeasurementFunction, "MULTITHREAD RN");
    MeasurementFunction->Write();

    cout << "chi2 = " << MeasurementFunction->GetChisquare() / MeasurementFunction->GetNDF() << endl;

    TCanvas *c1 = new TCanvas("MeassurementFitted_2D_View", "MeassurementFitted_2D_View", 800, 800);
    H_measurement_reconstruted->Draw("COLZ");
    MeasurementFunction->Draw("SAME");
    c1->Write();

    TCanvas *c = new TCanvas("MeassurementFitted_3D_View", "MeassurementFitted_3D_View", 800, 800);
    H_measurement_reconstruted->GetXaxis()->SetRangeUser(-2, 2);
    H_measurement_reconstruted->GetYaxis()->SetRangeUser(-2, 2);
    H_measurement_reconstruted->Draw("LEGO2");
    MeasurementFunction->SetNpx(75);
    MeasurementFunction->SetNpy(75);
    MeasurementFunction->Draw("SURF SAME");
    c->Write();

    

    /////// SAVING PARAMETER ////////
    SIGMA_BEAM_X = make_pair(MeasurementFunction->GetParameter(6), MeasurementFunction->GetParError(6));
    SIGMA_BEAM_Y = make_pair(MeasurementFunction->GetParameter(7), MeasurementFunction->GetParError(7));
    MEAN_BEAM_X = make_pair(MeasurementFunction->GetParameter(4), MeasurementFunction->GetParError(4));
    MEAN_BEAM_Y = make_pair(MeasurementFunction->GetParameter(5), MeasurementFunction->GetParError(5));
    AMPLITUDE_BEAM = make_pair(MeasurementFunction->GetParameter(3), MeasurementFunction->GetParError(3));
    /////////////////////////////////

    TF2 *fMyGaussian = new TF2("Beam_Profile", MyGaussian, -2, 2, -2, 2, 5);
    fMyGaussian->SetParameter(0, MeasurementFunction->GetParameter(3));
    fMyGaussian->SetParameter(1, MeasurementFunction->GetParameter(4));
    fMyGaussian->SetParameter(2, MeasurementFunction->GetParameter(6));
    fMyGaussian->SetParameter(3, MeasurementFunction->GetParameter(5));
    fMyGaussian->SetParameter(4, MeasurementFunction->GetParameter(7));
    fMyGaussian->Write();

}

void Measurement2D()
{
    H_measurement_reconstruted = new TH2D("H_measurement_reconstruted", "H_measurement_reconstruted", 100, x_measurement_min, x_measurement_max, 100, x_measurement_min, x_measurement_max);
    for (int i = 0; i < 60202; ++i)
    {
        double x, y;
        H_measurement->GetRandom2(x, y);
        double x_fit = f_X->Eval(x, y)-0.27;
        double y_fit = f_Y->Eval(y, x)+0.30;
        H_measurement_reconstruted->Fill(x_fit, y_fit);
    }

    H_measurement_reconstruted->SetBinContent(57, 51, 0);


    TCanvas *c_measurement_reconstruted = new TCanvas("Measurement_2D_View", "Measurement_2D_View", 800, 800);
    H_measurement_reconstruted->Draw("COLZ");
    c_measurement_reconstruted->Write();
    
}

/// 1D diag

double func1d(double *x, double *par)
{
    double A1 = par[0];
    double A2 = par[8];
    double mu1 = par[1];
    double mu2 = par[2];
    double sigmar = par[3];

    double A_g = par[4];
    double mu_g = par[5];
    double sigma_g = par[6];

    double grid1 = ( erf((x[0] - mu1 + l/2)/(sqrt(2) * sigmar)) - erf((x[0] - mu1 - l/2)/(sqrt(2) * sigmar)));
    double grid2 = (erf((x[0] - mu2 + l/2)/(sqrt(2) * sigmar)) - erf((x[0] - mu2 - l/2)/(sqrt(2) * sigmar)) );
    double gauss = A_g * exp(-0.5 * (x[0] - mu_g) * (x[0] - mu_g) / (sigma_g * sigma_g));

    return par[7]+ A1 * grid1 + A2 * grid2 + gauss*(grid1+grid2);
}

void testdiag()
{
    TH2D* H = new TH2D("H", "H", 80, -2, 2, 80, -2, 2);
    for (int i = 0; i < 60202; ++i)
    {
        double x, y;
        H_measurement->GetRandom2(x, y);
        double x_fit = f_X->Eval(x, y)-0.27;;
        double y_fit = f_Y->Eval(y, x)+0.30;;

        double x_corr = x_fit;
        double y_corr = y_fit;

        double r_corr = sqrt(x_corr * x_corr + y_corr * y_corr);
        double phi_corr = atan2(y_corr, x_corr);

        double x_diag = r_corr * cos(phi_corr - 45 * TMath::Pi() / 180);
        double y_diag = r_corr * sin(phi_corr - 45 * TMath::Pi() / 180);
        H->Fill(x_diag, y_diag);
    }

    H->Write();


    TH1D* H_diag = H->ProjectionX("H_diag", 37, 42);

    TF1* f1 = new TF1("f1", func1d, -3, 3, 9);
    f1->SetParLimits(0, 0, 300);
    f1->SetParLimits(1, -1.5, -0.8);
    f1->SetParameter(1, -1);
    f1->SetParLimits(2, 0.8, 1.4);
    f1->SetParameter(2, 1);
    f1->SetParLimits(3, 0., 1);
    f1->SetParameter(3, 0.3);
    f1->SetParLimits(4, 0, 1000000);
    f1->SetParameter(4, 1000);
    f1->SetParLimits(5, -1, 1.0);
    // f1->SetParameter(5, -0.5);
    f1->SetParLimits(6, 0.2, 3);
    f1->SetParameter(7, 10);
    f1->SetParLimits(8, 0, 500);

    TCanvas *c = new TCanvas("fit", "fit", 800, 800);
    H_diag->Fit(f1, "RN");
    H_diag->Draw("HIST");
    f1->Draw("SAME");
    c->Write();


    TH1D* H_diagy = H->ProjectionY("H_diagy", 37, 42);

    TF1* f1y = new TF1("f1", func1d, -3, 3, 9);
    f1y->SetParLimits(0, 0, 1000);
    f1y->SetParLimits(1, -1.5, -0.8);
    f1y->SetParameter(1, -1);
    f1y->SetParLimits(2, 0.8, 1.4);
    f1y->SetParameter(2, 1);
    f1y->SetParLimits(3, 0., 1);
    f1y->SetParameter(3, 0.3);
    f1y->SetParLimits(4, 0, 100000);
    f1y->SetParameter(4, 10000);
    f1y->SetParLimits(5, -1, 2.0);
    // f1->SetParameter(5, -0.5);
    f1y->SetParLimits(6, 0.2, 10);
    f1y->SetParameter(7, 10);
    f1y->SetParLimits(8, 0, 1000);

    TCanvas *cy = new TCanvas("fity", "fity", 800, 800);
    H_diagy->Fit(f1y, "RN");
    H_diagy->Draw("HIST");
    f1y->Draw("SAME");
    cy->Write();
}

/////// print results
void PrintResults()
{
    //// 2D
    cout << "########## 2D ##########" << endl;
    cout << "-----MCP-----" << endl;
    cout << "Resolution X = " << SIGMA_X.first << " +/- " << SIGMA_X.second << endl;
    cout << "Resolution Y = " << SIGMA_Y.first << " +/- " << SIGMA_Y.second << endl;
    cout << "-----BEAM-----" << endl;
    cout << "AMPLITUDE = " << AMPLITUDE_BEAM.first << " +/- " << AMPLITUDE_BEAM.second << endl;
    cout << "Mean X = " << MEAN_BEAM_X.first << " +/- " << MEAN_BEAM_X.second << endl;
    cout << "Mean Y = " << MEAN_BEAM_Y.first << " +/- " << MEAN_BEAM_Y.second << endl;
    cout << "Size X = " << SIGMA_BEAM_X.first << " +/- " << SIGMA_BEAM_X.second << endl;
    cout << "Size Y = " << SIGMA_BEAM_Y.first << " +/- " << SIGMA_BEAM_Y.second << endl;
    
    //// 1D
    cout << endl;
    cout << "########## 1D ##########" << endl;
    cout << "-----BEAM-----" << endl;
    // cout << "AMPLITUDE = " << AMPLITUDE_BEAM_1D.first << " +/- " << AMPLITUDE_BEAM_1D.second << endl;


}
#endif