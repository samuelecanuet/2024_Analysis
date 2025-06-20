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
#include <Minuit2/Minuit2Minimizer.h>
#include "TFitResult.h"
#include "TVirtualFitter.h"
#include <gsl/gsl_statistics.h>

#include "../../Grouper/include/Detectors.hh"

using namespace std;
using namespace ROOT::Math;
default_random_engine generator;

template <size_t ROWS, size_t COLS>
vector<vector<bool>> ConvertToVector(bool (&array)[ROWS][COLS]) {
    vector<vector<bool>> vec(ROWS, vector<bool>(COLS));
    for (size_t i = 0; i < ROWS; i++)
        for (size_t j = 0; j < COLS; j++)
            vec[i][j] = array[i][j];
    return vec;
}

/// FILE ///

TFile *ROOT_file;
TFile *FINAL_file;
TTreeReader *Reader;
TTree *Tree;
TTree *Cleaned_Tree;

string Calibration_Filename;
string Measurement_Filename;

TFile* fSaved;
bool WRITTING;
string year;

TH1D *H_RAW[5];
TH2D *H_RAW_2D;
TH2D *H;
TH2D *H_precorrrected;
TH2D *H_corrected;
TH2D *H_reconstruction;
TH2D *H_measurement;
TH2D *H_measurement_reconstructed;

int A = 1;
int B = 3;
int C = 4;
int D = 2;

double rho_real = 2;
double l_real;
double rho = 0.1;
double l = 0.08;
double limit_real = 7.5;
int n;
double limit = limit_real * rho / rho_real;
int MAXI = 30;

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

bool M8_2024[8][8] = {
    {false, false, false, false, false, false, false, false},
    {false, false, true, true, true, true, false, false},
    {false, true, true, true, true, true, true, false},
    {false, true, true, true, true, true, true, false},
    {false, true, true, true, true, true, true, false},
    {false, true, true, true, true, true, true, false},
    {false, false, true, true, true, true, false, false},
    {false, false, false, false, false, false, false, false}};

bool M7_2025[7][7] = {
    {false, true, true, true, true, true, false},
    {true, true, true, true, true, true, true},
    {true, true, true, true, true, true, true},
    {true, true, true, true, true, true, true},
    {true, true, true, true, true, true, true},
    {true, true, true, true, true, true, true},
    {false, true, true, true, true, true, false}};

    // {false, false, false, false, false, false, false},
    // {false, true, true, true, true, true, false},
    // {false, true, true, true, true, true, false},
    // {false, true, true, true, true, true, false},
    // {false, true, true, true, true, true, false},
    // {false, true, true, true, true, true, false},
    // {false, false, false, false, false, false, false}};

    // {false, false, false, false, false, false, false},
    // {false, false, false, false, false, false, false},
    // {false, false, true, true, true, false, false},
    // {false, false, true, true, true, false, false},
    // {false, false, true, true, true, false, false},
    // {false, false, false, false, false, false, false},
    // {false, false, false, false, false, false, false}};

vector<vector<bool>> M;

double x_final_max = 8;
double x_final_min = -8;
double y_final_max = 8;
double y_final_min = -8;

bool M_final_2024[8][8] = {
    {false, false, false, false, false, false, false, false},
    {false, false, false, false, false, false, false, false},
    {false, false, true, true, true, true, false, false},
    {false, false, true, true, true, true, false, false},
    {false, false, true, true, true, true, false, false},
    {false, false, true, true, true, true, false, false},
    {false, false, false, false, false, false, false, false},
    {false, false, false, false, false, false, false, false}};

bool M_final_2025[7][7] = {
    {false, false, false, false, false, false, false},
    {false, false, false, false, false, false, false},
    {false, false, true, true, true, false, false},
    {false, false, true, true, true, false, false},
    {false, false, true, true, true, false, false},
    {false, false, false, false, false, false, false},
    {false, false, false, false, false, false, false}};

vector<vector<bool>> M_final;

double x_measurement_min = -4.;
double x_measurement_max = 4.;
double y_measurement_min = -4.;
double y_measurement_max = 4.;

bool M_measurement_2024[8][8] = {
    {false, false, false, false, false, false, false, false},
    {false, false, false, false, false, false, false, false},
    {false, false, false, false, false, false, false, false},
    {false, false, false, true, true, false, false, false},
    {false, false, false, true, true, false, false, false},
    {false, false, false, false, false, false, false, false},
    {false, false, false, false, false, false, false, false},
    {false, false, false, false, false, false, false, false}};

bool M_measurement_2025[7][7] = {
    {false, false, false, false, false, false, false},
    {false, false, false, false, false, false, false},
    {false, false, true, true, true, false, false},
    {false, false, true, true, true, false, false},
    {false, false, true, true, true, false, false},
    {false, false, false, false, false, false, false},
    {false, false, false, false, false, false, false}};

vector<vector<bool>> M_measurement;

pair<double, double> M_center[64];
vector<pair<double, double>> M_corner[64];

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

void writefit()
{
    ofstream outfile(("fit_params_"+year+".txt").c_str());
    // outfile << xmax << " " << xmin << " " << ymax << " " << ymin << endl;
    // for (int i = 0; i < Nparam; i++)
    // {
    //     outfile << f_X->GetParameter(i) << " ";
    // }

    for (int i = 0; i <= NparX; i++)
    {
        for (int j = 0; j <= NparY; j++)
        {
            outfile << f_X->GetParameter(i * (NparX + 1) + j) << " ";
        }
    }
    outfile << endl;
    // for (int i = 0; i < Nparam; i++)
    // {
    //     outfile << f_Y->GetParameter(i) << " ";
    // }
    for (int i = 0; i <= NparX; i++)
    {
        for (int j = 0; j <= NparY; j++)
        {
            outfile << f_Y->GetParameter(i * (NparX + 1) + j) << " ";
        }
    }
    outfile << endl;
    outfile.close();
}

void readfit(string Year = "2024")
{
    ifstream infile(("fit_params_"+Year+".txt").c_str());
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

    double norm =1. / log((Va * Vb * Vc * Vd) / pow(Va + Vb + Vc + Vd, 4));
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


void InitYearConfiguration(string Year = "2024")
{   
    if (Year == "2024")
    {
        M = ConvertToVector(M8_2024);
        M_final = ConvertToVector(M_final_2024);
        M_measurement = ConvertToVector(M_measurement_2024);
        n = 8;
        l_real = 1.2;
        Calibration_Filename = "MCP_008_4T_0001.root";
        Measurement_Filename = "/RAW/MCP_010_4T.root";

    }
    else if (Year == "2025")
    {
        M = ConvertToVector(M7_2025);
        M_final = ConvertToVector(M_final_2025);
        M_measurement = ConvertToVector(M_measurement_2025);
        n = 7;
        l_real = 1.4;
        Calibration_Filename = "/mnt/hgfs/shared-2/2025_DATA/MCP_DATA/run_001_StableBeamScan_grouped.root";  
        // Measurement_Filename = "run_041_MCP_32Ar_Beam_4T_grouped.root";
        Measurement_Filename = "run_079_MCP_32Ar_Heinz14kV_grouped.root";

        // Calibration_Filename = "/mnt/hgfs/shared-2/2025_DATA/MCP_DATA/03_TEST/005_MCP_1p9kV_BeamScan.fast/005_MCP_1p9kV_BeamScan_0001.root";  
        // Measurement_Filename = "/mnt/hgfs/shared-2/2025_DATA/MCP_DATA/03_TEST/005_MCP_1p9kV_BeamScan.fast/005_MCP_1p9kV_BeamScan_0001.root";  
    }
    else
    {
        Error("InitYearConfiguration", "Year not found");
    }
}
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
    double l2x;
    double l2y;
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

            if (!M[i][j])
                continue;

            if (abs(sqrt(M_center[i * n + j].first * M_center[i * n + j].first + M_center[i * n + j].second * M_center[i * n + j].second) - r) > 2 * rho)
            {
                continue;
            }

            beta_x[i][j] = par[10 + i * n + j];
            beta_y[i][j] = par[10 + n * n + i * n + j];

            if (par[10 + 2 * n * n + i * n + j] != 0)
                l2x = par[10 + 2 * n * n + i * n + j] / 2;
            else
                l2x = l2;

            if (par[10 + 3 * n * n + i * n + j] != 0)
                l2y = par[10 + 3 * n * n + i * n + j] / 2;
            else
                l2y = l2;

            double erf_x_pos = erf((x[0] - (beta_x[i][j] - l2x)) / (sqrt(2) * sigma_x));
            double erf_x_neg = erf((x[0] - (beta_x[i][j] + l2x)) / (sqrt(2) * sigma_x));
            double erf_y_pos = erf((x[1] - (beta_y[i][j] - l2y)) / (sqrt(2) * sigma_y));
            double erf_y_neg = erf((x[1] - (beta_y[i][j] + l2y)) / (sqrt(2) * sigma_y));

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
    //     // double sigma_r =1. / sqrt(2) * sqrt(pow(cos(theta) * sigma_x, 2) + pow(sin(theta) * sigma_y, 2));
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
    Info("Fitting grid with the guess");
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
    // FittedFunction->FixParameter(0, 30);

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
    FittedFunction->SetParLimits(9, 0.04, 0.1);
    FittedFunction->SetParameter(9, 0.08);
    // FittedFunction->FixParameter(9, 0.08);

    for (int i = 0; i < N; ++i)
    {
        // amplitude bkg in MCP
        // FittedFunction->SetParLimits(i, 0., 500);
        // FittedFunction->SetParameter(i, 10);

        if (!M[i / n][i % n])
        {
            // FittedFunction->FixParameter(i, 0);
            // beta_x
            FittedFunction->FixParameter(10 + i, 0.);
            // beta_y
            FittedFunction->FixParameter(N + 10 + i, 0.);
        }
        else
        {
            double delta = 0.03;
            // FittedFunction->FixParameter(i, 0);
            // double x_center = i / n * rho - rho * (n-1)/2-0.01;
            // double y_center = i % n * rho - rho * (n-1)/2+0.017;

            // beta_x
            FittedFunction->SetParLimits(10 + i, M_center[i].first -delta, M_center[i].first +delta);
            FittedFunction->SetParameter(10 + i, M_center[i].first);
            // FittedFunction->FixParameter(10 + i, M_center[i].first);


            // beta_y
            FittedFunction->SetParLimits(N + 10 + i, M_center[i].second -delta, M_center[i].second +delta);
            FittedFunction->SetParameter(N + 10 + i, M_center[i].second);
            // FittedFunction->FixParameter(N + 10 + i, M_center[i].second);
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

    // H_precorrrected->SetBinContent(182, 178, 0);
    H_precorrrected->Fit(FittedFunction, "MULTITHREAD RN", "", -0.4, 0.4);

    cout << "chi2 = " << FittedFunction->GetChisquare() / FittedFunction->GetNDF() << endl;

    TCanvas *c = new TCanvas("c", "c", 800, 800);
    H_precorrrected->Draw("COLZ");
    FittedFunction->Draw("SAME");
    c->Write();

    //// writing center of grid cell in a txt file with cell number x an y
    ofstream file;
    file.open(("out_centers_"+year+".txt").c_str());
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
            if (!M[i][j])
                continue;
            double x_center = M_center[i * n + j].first;
            double y_center = M_center[i * n + j].second;

            // cout << "Cell " << i * n + j << " : " << x_center << " " << y_center << endl;

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

////////// TRYING WITH CORNERS /////////////////////////////////

// Compute signed distance to a line
double signed_distance(double x, double y, double a, double b) {
    return (a * x - y + b) / std::sqrt(a * a + 1);
}
double MyFullFittedFunction2D_CORNER(double *x, double *par)
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

    // if (r > 0.5)
    // {
    //     return 0;
    // }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {

            if (!M[i][j])
                continue;

            double Ax = par[10 + i * n + j];
            double Ay = par[10 + n * n + i * n + j];
            double Bx = par[10 + 2 * n * n + i * n + j];
            double By = par[10 + 3 * n * n + i * n + j];
            double Cx = par[10 + 4 * n * n + i * n + j];
            double Cy = par[10 + 5 * n * n + i * n + j];
            double Dx = par[10 + 6 * n * n + i * n + j];
            double Dy = par[10 + 7 * n * n + i * n + j];

            // if (abs(sqrt(Ax*Ax + Ay*Ay) - r) > 2 * rho)
            // {
            //     continue;
            // }

            double aAB = (Ay - By) / (Ax - Bx);
            double aBC = (Bx - Cx) / (By - Cy);
            double aCD = (Cy - Dy) / (Cx - Dx);
            double aDA = (Dx - Ax) / (Dy - Ay);

            double bAB = Ay - aAB * Ax;
            double bBC = Bx - aBC * By;
            double bCD = Cy - aCD * Cx;
            double bDA = Dx - aDA * Dy;

            double edgeAB = 0.5*(1+erf((x[1]-aAB*x[0]-bAB)/(sqrt(2)*sigma_y)));
            double edgeBC = 0.5*(1-erf((x[0]-aBC*x[1]-bBC)/(sqrt(2)*sigma_x)));
            double edgeCD = 0.5*(1-erf((x[1]-aCD*x[0]-bCD)/(sqrt(2)*sigma_y)));
            double edgeDA = 0.5*(1+erf((x[0]-aDA*x[1]-bDA)/(sqrt(2)*sigma_x)));

            grid_sum += A * edgeAB * edgeBC * edgeCD * edgeDA;
        }
    }

    double result = grid_sum;

    if (result > MAXI)
    {
        return MAXI;
    }

    return result;
}

double Function(double *x, double *par)
{
    double res = 0;
    for (int i = 0; i <= NparX; i++)
    {
        for (int j = 0; j <= NparY; j++)
        {
            res += par[i * (NparX + 1) + j] * pow(x[0], i) * pow(x[1], j);
        }
    }

    return res;

    // return 1 * (par[1] * x[0] + par[2] + pow(x[0], 2) + par[3] * pow(x[0], 3) + par[6] * x[1] * x[1] * x[1] * x[0] + par[0] * pow(x[0], 4) + par[7] * pow(x[0], 5));
}

void FullFunctionToMinimize2D_CORNER()
{
    Info("Fitting grid with the guess");
    double chi2 = 0.0;
    int N = n * n;
    FittedFunction = new TF2("FittedFunction", MyFullFittedFunction2D_CORNER, -0.4, 0.4, -0.4, 0.4, 8 * N + 10);
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
    FittedFunction->SetParLimits(9, 0.04, 0.1);
    // FittedFunction->SetParameter(9, 0.08);
    // FittedFunction->FixParameter(9, 0.08);

    for (int i = 0; i < N; ++i)
    {
        // amplitude bkg in MCP
        // FittedFunction->SetParLimits(i, 0., 500);
        // FittedFunction->SetParameter(i, 10);

        if (!M[i / n][i % n])
        {
            // FittedFunction->FixParameter(i, 0);
            // a_x
            FittedFunction->FixParameter(10 + i, 0.);
            // a_y
            FittedFunction->FixParameter(N + 10 + i, 0.);
            // b_x
            FittedFunction->FixParameter(2 * N + 10 + i, 0.);
            // b_y
            FittedFunction->FixParameter(3 * N + 10 + i, 0.);
            // c_x
            FittedFunction->FixParameter(4 * N + 10 + i, 0.);
            // c_y
            FittedFunction->FixParameter(5 * N + 10 + i, 0.);
            // d_x
            FittedFunction->FixParameter(6 * N + 10 + i, 0.);
            // d_y
            FittedFunction->FixParameter(7 * N + 10 + i, 0.);

        }
        else
        {
            double delta = 0.02;
            // a_x
            FittedFunction->SetParLimits(10 + i, M_corner[i][0].first -delta, M_corner[i][0].first +delta);
            FittedFunction->SetParameter(10 + i, M_corner[i][0].first);
            // FittedFunction->FixParameter(10 + i, M_corner[i][0].first);

            // a_y
            FittedFunction->SetParLimits(N + 10 + i, M_corner[i][0].second -delta, M_corner[i][0].second +delta);
            FittedFunction->SetParameter(N + 10 + i, M_corner[i][0].second);
            // FittedFunction->FixParameter(N + 10 + i, M_corner[i][0].second);

            // b_x
            FittedFunction->SetParLimits(2 * N + 10 + i, M_corner[i][1].first -delta, M_corner[i][1].first +delta);
            FittedFunction->SetParameter(2 * N + 10 + i, M_corner[i][1].first);
            // FittedFunction->FixParameter(2 * N + 10 + i, M_corner[i][1].first);

            // b_y
            FittedFunction->SetParLimits(3 * N + 10 + i, M_corner[i][1].second -delta, M_corner[i][1].second +delta);
            FittedFunction->SetParameter(3 * N + 10 + i, M_corner[i][1].second);
            // FittedFunction->FixParameter(3 * N + 10 + i, M_corner[i][1].second);

            // c_x
            FittedFunction->SetParLimits(4 * N + 10 + i, M_corner[i][2].first -delta, M_corner[i][2].first +delta);
            FittedFunction->SetParameter(4 * N + 10 + i, M_corner[i][2].first);
            // FittedFunction->FixParameter(4 * N + 10 + i, M_corner[i][2].first);

            // c_y
            FittedFunction->SetParLimits(5 * N + 10 + i, M_corner[i][2].second -delta, M_corner[i][2].second +delta);
            FittedFunction->SetParameter(5 * N + 10 + i, M_corner[i][2].second);
            // FittedFunction->FixParameter(5 * N + 10 + i, M_corner[i][2].second);

            // d_x
            FittedFunction->SetParLimits(6 * N + 10 + i, M_corner[i][3].first -delta, M_corner[i][3].first +delta);
            FittedFunction->SetParameter(6 * N + 10 + i, M_corner[i][3].first);
            // FittedFunction->FixParameter(6 * N + 10 + i, M_corner[i][3].first);

            // d_y
            FittedFunction->SetParLimits(7 * N + 10 + i, M_corner[i][3].second -delta, M_corner[i][3].second +delta);
            FittedFunction->SetParameter(7 * N + 10 + i, M_corner[i][3].second);
            // FittedFunction->FixParameter(7 * N + 10 + i, M_corner[i][3].second);
        }
    }

    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0);
    H_precorrrected->Rebin2D(2, 2);
    H_precorrrected->GetXaxis()->SetRangeUser(-0.4, 0.4);
    H_precorrrected->GetYaxis()->SetRangeUser(-0.4, 0.4);

    // H_precorrrected->SetBinContent(182, 178, 0);
    H_precorrrected->Fit(FittedFunction, "MULTITHREAD RN", "", -0.4, 0.4);

    cout << "chi2 = " << FittedFunction->GetChisquare() / FittedFunction->GetNDF() << endl;

    TCanvas *c = new TCanvas("c", "c", 800, 800);
    H_precorrrected->Draw("COLZ");
    FittedFunction->Draw("SAME");
    c->Write();

    //// writing center of grid cell in a txt file with cell number x an y
    ofstream file;
    file.open(("out_corner_"+year+".txt").c_str());
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            int I = i * n + j;
            file << i * n + j << " " << FittedFunction->GetParameter(10 + I) << " " << FittedFunction->GetParameter(N + 10 + I) << " " << FittedFunction->GetParameter(10 + 2 * N + I) << " " << FittedFunction->GetParameter(N + 10 + 2 * N + I) << " " << FittedFunction->GetParameter(10 + 4 * N + I) << " " << FittedFunction->GetParameter(N + 10 + 4 * N + I) << " " << FittedFunction->GetParameter(10 + 6 * N + I) << " " << FittedFunction->GetParameter(N + 10 + 6 * N + I) << endl;
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

    TCanvas *c_corner = new TCanvas("c_corner", "c_corner", 800, 800);
    H_precorrrected->Draw("COLZ");
    for (int i = 0; i < N; ++i)
    {
        if (M[i / n][i % n])
        {

            TGraph *c_corner = new TGraph(4);
            // a
            c_corner->SetPoint(0, FittedFunction->GetParameter(10 + i), FittedFunction->GetParameter(N + 10 + i));
            // b
            c_corner->SetPoint(1, FittedFunction->GetParameter(10 + 2 * N + i), FittedFunction->GetParameter(N + 10 + 2 * N + i));
            // c
            c_corner->SetPoint(2, FittedFunction->GetParameter(10 + 4 * N + i), FittedFunction->GetParameter(N + 10 + 4 * N + i));
            // d
            c_corner->SetPoint(3, FittedFunction->GetParameter(10 + 6 * N + i), FittedFunction->GetParameter(N + 10 + 6 * N + i));


            c_corner->SetMarkerStyle(20);
            c_corner->SetMarkerSize(2);
            c_corner->Draw("P SAME");
        }
    }
    c_corner->Write();

    FittedFunction->Write();
}

void FittingReconstruction_CORNER()
{
    int counter = 0;
    int N = n * n;
    TGraphErrors *G_coord_fitgrid = new TGraphErrors();
    if (FittedFunction != nullptr)
    {
        for (int i = 0; i < n * n; ++i)
        {
            if (M[i / n][i % n])
            {
                G_coord_fitgrid->SetPoint(counter, FittedFunction->GetParameter(10 + i), FittedFunction->GetParameter(N + 10 + i));
                G_coord_fitgrid->SetPointError(counter, FittedFunction->GetParError(10 + i), FittedFunction->GetParError(N + 10 + i));
                counter++;

                G_coord_fitgrid->SetPoint(counter, FittedFunction->GetParameter(2*N + 10 + i), FittedFunction->GetParameter(N + 10 + 2*N + i));
                G_coord_fitgrid->SetPointError(counter, FittedFunction->GetParError(2*N + 10 + i), FittedFunction->GetParError(N + 10 + 2*N + i));
                counter++;

                G_coord_fitgrid->SetPoint(counter, FittedFunction->GetParameter(4*N + 10 + i), FittedFunction->GetParameter(N + 10 + 4*N + i));
                G_coord_fitgrid->SetPointError(counter, FittedFunction->GetParError(4*N + 10 + i), FittedFunction->GetParError(N + 10 + 4*N + i));
                counter++;

                G_coord_fitgrid->SetPoint(counter, FittedFunction->GetParameter(6*N + 10 + i), FittedFunction->GetParameter(N + 10 + 6*N + i));
                G_coord_fitgrid->SetPointError(counter, FittedFunction->GetParError(6*N + 10 + i), FittedFunction->GetParError(N + 10 + 6*N + i));
                counter++;
            }
        }
        TCanvas *c_coord_fitgrid = new TCanvas("c_coord_fitgrid", "c_coord_fitgrid", 800, 800);
        G_coord_fitgrid->SetMarkerStyle(20);
        G_coord_fitgrid->GetXaxis()->SetTitle("x [mm]");
        G_coord_fitgrid->GetYaxis()->SetTitle("y [mm]");
        G_coord_fitgrid->Draw("AP");
        c_coord_fitgrid->cd();
        for (int i = 0; i < G_coord_fitgrid->GetN(); ++i)
        {
            double x = G_coord_fitgrid->GetX()[i];
            double y = G_coord_fitgrid->GetY()[i];

            TText *t = new TText(x, y, Form("%d", i));
            t->SetTextSize(0.02);
            t->Draw("SAME");
        }
        c_coord_fitgrid->Write();
    }
    else
    {
        std::ifstream file(("out_corner_"+year+".txt").c_str());
        std::string line;
        int counter = 0;
        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            double Ax, Ay, Bx, By, Cx, Cy, Dx, Dy;
            int i;
            if (!(iss >> i >> Ax >> Ay >> Bx >> By >> Cx >> Cy >> Dx >> Dy))
            {
                break;
            }
            i = i - 1;
            if (!M[i / n][i % n])
            {
                continue;
            }
            G_coord_fitgrid->SetPoint(counter, Ax, Ay);
            G_coord_fitgrid->SetPoint(counter + 1, Bx, By);
            G_coord_fitgrid->SetPoint(counter + 2, Cx, Cy);
            G_coord_fitgrid->SetPoint(counter + 3, Dx, Dy);
            // cout << N << " " << x << " " << y << endl;
            counter+=4;
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
    for (int i = 0; i < n * n; i++)
    {
        if (M[i / n][i % n])
        {
            double X_real = i % n * rho - rho * (n-1) / 2 - l_real / 2;
            double Y_real = i / n * rho - rho * (n-1) / 2 - l_real / 2;

            double x, y;
            G_coord_fitgrid->GetPoint(counter, x, y);
            G2D_X->SetPoint(counter, x, y, X_real);
            G2D_Y->SetPoint(counter, y, x, Y_real);

            G_coord_fitgrid->GetPoint(counter + 1, x, y);
            G2D_X->SetPoint(counter + 1, x, y, X_real + l_real);
            G2D_Y->SetPoint(counter + 1, y, x, Y_real);

            G_coord_fitgrid->GetPoint(counter + 2, x, y);
            G2D_X->SetPoint(counter + 2, x, y, X_real + l_real);
            G2D_Y->SetPoint(counter + 2, y, x, Y_real + l_real);

            G_coord_fitgrid->GetPoint(counter + 3, x, y);
            G2D_X->SetPoint(counter + 3, x, y, X_real);
            G2D_Y->SetPoint(counter + 3, y, x, Y_real + l_real);

            counter += 4;

            // cout << "Cell " << i << " : " << X_real << " " << Y_real << endl;
        }
    }
    
    // fitting X
    TCanvas *c_fit_X = new TCanvas("c_fit_X", "c_fit_X", 800, 800);
    G2D_X->SetMarkerStyle(20);
    f_X = new TF2("f_X", Function, G2D_X->GetXmin(), G2D_X->GetXmax(), G2D_X->GetYmin(), G2D_X->GetYmax(), (NparX + 1) * (NparY + 1));
    G2D_X->Fit("f_X", "MULTITHREAD");
    f_X->Draw("SURF2");
    G2D_X->Draw("AP SAME");
    c_fit_X->Write();

    // fitting Y
    TCanvas *c_fit_Y = new TCanvas("c_fit_Y", "c_fit_Y", 800, 800);
    G2D_Y->SetMarkerStyle(20);

    f_Y = new TF2("f_Y", Function, G2D_Y->GetXmin(), G2D_Y->GetXmax(), G2D_Y->GetYmin(), G2D_Y->GetYmax(), (NparX + 1) * (NparY + 1));
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

        c_residus_XY->cd();
        TText *t = new TText(x, y, Form("%d", i));
        t->SetTextSize(0.02);
        t->Draw();

    }

    G_real_XY->SetMarkerStyle(20);
    G_real_XY->Draw("AP");
    G_residus_XY->SetMarkerStyle(20);
    G_residus_XY->SetMarkerColor(kRed);
    G_residus_XY->Draw("P SAME");
    

    for (int i = 0; i < G2D_X->GetN(); ++i)
    {
        double x = G2D_X->GetZ()[i];
        double y = G2D_Y->GetZ()[i];

        TText *t = new TText(x, y, Form("%d", i));
        t->SetTextSize(0.02);
        t->Draw("SAME");
    }

    c_residus_XY->Write();

    // Reconstruction of the grid
    TCanvas *c_reconstruction = new TCanvas("c_reconstruction", "c_reconstruction", 800, 800);
    H_reconstruction = new TH2D("H_reconstruction", "H_reconstruction", 800, -8, 8, 800, -8, 8);

    if (year == "2024")
    {
        for (int i = 0; i < H_precorrrected->GetEntries(); ++i)
        {
            double x, y;
            H_precorrrected->GetRandom2(x, y);
            double x_fit = f_X->Eval(x, y);
            double y_fit = f_Y->Eval(y, x);
            H_reconstruction->Fill(x_fit, y_fit);
        }
    }
    else
    {
        TFile *f = MyTFile((Calibration_Filename).c_str(), "READ");
        TTree *tree = (TTree *)f->Get("treeMCP");
        TTreeReader *Reader = new TTreeReader(tree);
        TTreeReaderValue<double> *X0 = new TTreeReaderValue<double>(*Reader, "X");
        TTreeReaderValue<double> *Y0 = new TTreeReaderValue<double>(*Reader, "Y");

        int Entries = tree->GetEntries();
        clock_t start = clock(), Current;
        while(Reader->Next())
        {
            // ProgressBar()
            double x = **X0;
            double y = **Y0;
            double x_fit = f_X->Eval(x, y);
            double y_fit = f_Y->Eval(y, x);
            H_reconstruction->Fill(x_fit, y_fit);
        }
        f->Close();
        FINAL_file->cd();
    }
    H_reconstruction->Draw("COLZ");
    c_reconstruction->Write();

    writefit();
}

///////////////// Fitting Points /////////////////////////////////////
void LoadCenters(string Year = "2024")
{
    if (Year == "2024")
    {
        ifstream file("guess_center_2024.txt");
        string line;
        while (getline(file, line))
        {
            istringstream iss(line);
            int cell;
            double x, y;
            iss >> cell >> x >> y;
            M_center[cell] = make_pair(x, y);
        }
        file.close();
        Success("Centers loaded");   
    }
    else if (Year == "2025")
    {
        ifstream file("guess_center_2025.txt");
        string line;
        while (getline(file, line))
        {
            istringstream iss(line);
            int cell;
            double x, y;
            iss >> cell >> x >> y;
            M_center[cell] = make_pair(x, y);
        }
        file.close();
        Success("Centers loaded"); 

        ifstream file2("guess_corner_2025.txt");  
        if (!file2.is_open())
        {
            Error("LoadCenters", "File guess_corner_2025.txt not found");
            return;
        }
        string line2;
        while (getline(file2, line2))
        {
            istringstream iss(line2);
            int cell;
            double ax, ay, bx, by, cx, cy, dx, dy;
            iss >> cell >> ax >> ay >> bx >> by >> cx >> cy >> dx >> dy;
            cell = cell - 1;
            M_corner[cell] = {make_pair(ax, ay), make_pair(bx, by), make_pair(cx, cy), make_pair(dx, dy)};

            // cout << cell << " " << ax << " " << ay << " " << bx << " " << by << " " << cx << " " << cy << " " << dx << " " << dy << endl;
        }
        


        FINAL_file->cd();
        TCanvas *c_PointVerification = new TCanvas("c_PointVerification", "c_PointVerification", 800, 800);
        H_precorrrected->Draw("COLZ");
        TGraph *G_coord_fitgrid = new TGraph();
        for (int i = 0; i < n * n; ++i)
        {
            if (M[i / n][i % n])
            {
                G_coord_fitgrid->AddPoint(M_corner[i][0].first, M_corner[i][0].second);
                G_coord_fitgrid->AddPoint(M_corner[i][1].first, M_corner[i][1].second);
                G_coord_fitgrid->AddPoint(M_corner[i][2].first, M_corner[i][2].second);
                G_coord_fitgrid->AddPoint(M_corner[i][3].first, M_corner[i][3].second);
            }
        }
        G_coord_fitgrid->SetMarkerStyle(20);
        G_coord_fitgrid->SetMarkerSize(2);
        G_coord_fitgrid->SetMarkerColor(kRed);
        G_coord_fitgrid->Draw("P SAME");
        c_PointVerification->Write();

        Success("Corner loaded"); 
    }
    else
    {
        Error("LoadCenters", "Year not found");
    }
}



void FittingReconstruction()
{
    int counter = 0;
    TGraphErrors *G_coord_fitgrid = new TGraphErrors();
    if (FittedFunction != nullptr)
    {
        for (int i = 0; i < n * n; ++i)
        {
            if (M[i / n][i % n])
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
        std::ifstream file(("out_centers_"+year+".txt").c_str());
        std::string line;
        int counter = 0;
        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            double x, y;
            int i;
            if (!(iss >> i >> x >> y))
            {
                break;
            }
            if (!M[i / n][i % n])
            {
                continue;
            }
            G_coord_fitgrid->SetPoint(counter, x, y);
            // cout << i << " " << x << " " << y << endl;
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
        if (M[i / n][i % n])
        {
            double X_real = i / n * rho - rho * (n - 1) / 2;
            double Y_real = i % n * rho - rho * (n - 1) / 2;

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
    f_X = new TF2("f_X", Function, G2D_X->GetXmin(), G2D_X->GetXmax(), G2D_X->GetYmin(), G2D_X->GetYmax(), (NparX + 1) * (NparY + 1));
    G2D_X->Fit("f_X", "MULTITHREAD");
    f_X->Draw("SURF2");
    G2D_X->Draw("AP SAME");
    c_fit_X->Write();

    // fitting Y
    TCanvas *c_fit_Y = new TCanvas("c_fit_Y", "c_fit_Y", 800, 800);
    G2D_Y->SetMarkerStyle(20);

    f_Y = new TF2("f_Y", Function, G2D_Y->GetXmin(), G2D_Y->GetXmax(), G2D_Y->GetYmin(), G2D_Y->GetYmax(), (NparX + 1) * (NparY + 1));
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
    H_reconstruction = new TH2D("H_reconstruction", "H_reconstruction", 800, -8, 8, 800, -8, 8);

    // if (year == "2024")
    // {
        for (int i = 0; i < H_precorrrected->GetEntries(); ++i)
        {
            double x, y;
            H_precorrrected->GetRandom2(x, y);
            double x_fit = f_X->Eval(x, y);
            double y_fit = f_Y->Eval(y, x);
            H_reconstruction->Fill(x_fit, y_fit);
        }
    // }
    // else
    // {
    //     TFile *f = new TFile(Calibration_Filename.c_str(), "READ");
    //     TTree *tree = (TTree *)f->Get("treeMCP");
    //     TTreeReader *Reader = new TTreeReader(tree);
    //     TTreeReaderValue<double> *X0 = new TTreeReaderValue<double>(*Reader, "X0");
    //     TTreeReaderValue<double> *Y0 = new TTreeReaderValue<double>(*Reader, "Y0");

    //     while(Reader->Next())
    //     {
    //         double x = **X0;
    //         double y = **Y0;
    //         double x_fit = f_X->Eval(x, y);
    //         double y_fit = f_Y->Eval(y, x);
    //         H_reconstruction->Fill(x_fit, y_fit);
    //     }
    //     f->Close();
    //     FINAL_file->cd();
    // }
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
    // FinalFunction->SetParameter(9, 1.2);
    FinalFunction->FixParameter(9, l_real);

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
            if (!M[i][j])
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
    double gauss = A_g /(2*M_PI*sigma_gx*sigma_gy) * exp(-0.5 * ((x[0] - mu_gx) * (x[0] - mu_gx) / (2 * sigma_gx * sigma_gx) + (x[1] - mu_gy) * (x[1] - mu_gy) / (2 * sigma_gy * sigma_gy)));

    return gauss;
}

double MeasurementFittedFunction2D(double *x, double *par)
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
    double l2 = par[9]/2;
    double beta_x[n][n];
    double beta_y[n][n];
    // double A[n][n];
    double a[6] = {par[3*n*n+10], par[3*n*n+11], par[3*n*n+21], par[3*n*n+13], par[3*n*n+14], par[3*n*n+15]};

    double result_sum = 0.0;
    double grid_sum = 0.0;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {

            if (!M_measurement[i][j])
                continue;

            // A[i][j] = par[2 * n * n + 10 + i * n + j];
            beta_x[i][j] = par[10 + i * n + j];
            beta_y[i][j] = par[10 + n * n + i * n + j];

            double erf_x_pos = 0.5*(1+erf((x[0] - (beta_x[i][j] - l2)) / (sqrt(2) * sigma_x)));
            double erf_x_neg = 0.5*(1+erf((x[0] - (beta_x[i][j] + l2)) / (sqrt(2) * sigma_x)));
            double erf_y_pos = 0.5*(1+erf((x[1] - (beta_y[i][j] - l2)) / (sqrt(2) * sigma_y)));
            double erf_y_neg = 0.5*(1+erf((x[1] - (beta_y[i][j] + l2)) / (sqrt(2) * sigma_y)));

            result_sum += (erf_x_pos - erf_x_neg) * (erf_y_pos - erf_y_neg);
            grid_sum += (erf_x_pos - erf_x_neg) * (erf_y_pos - erf_y_neg);
        }
    }

    // RAW EXPRESSION  OF 2D GAUSSIAN
    double theta = A * TMath::Pi() / 180.;
    double Xr = (x[0]-mu_gx)*cos(theta) + (x[1]-mu_gy)*sin(theta);
    double Yr = -(x[0]-mu_gx)*sin(theta) + (x[1]-mu_gy)*cos(theta);

    // double rho = 0.0;
    double gaussr = A_g *1./(sigma_gx*sigma_gy) * exp(-0.5 * ((Xr) * (Xr) / (sigma_gx * sigma_gx) + (Yr) * (Yr) / (sigma_gy * sigma_gy)));
    // double gaussr_corr = bkg + A_g * exp(-0.5 *1./(1-rho*rho) * ((Xr) * (Xr) / (sigma_gx * sigma_gx) + (Yr) * (Yr) / (sigma_gy * sigma_gy)) + 2 * rho * Xr * Yr / (sigma_gx * sigma_gy));

    double gauss = A_g *1./(sigma_gx*sigma_gy) * exp(-0.5 * ((x[0] - mu_gx) * (x[0] - mu_gx) / (sigma_gx * sigma_gx) + (x[1] - mu_gy) * (x[1] - mu_gy) / (sigma_gy * sigma_gy)));

    // RAW EXPRESSION  OF 2D GAUSSIAN (R)
    // double GaussianParameters[9] = {A_g, mu_gx, mu_gy, a[0], a[1], a[2], a[3], a[4], a[5]};
    
    result_sum *= gaussr;

    double result = result_sum + bkg*grid_sum;

    return result;
}

double f_CORNER(double *x, double *par)
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
    double l2 = par[9]/2;
    // double beta_x[n][n];
    // double beta_y[n][n];
    // double A[n][n];

    double result_sum = 0.0;
    double grid_sum = 0.0;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {

            if (!M_measurement[i][j])
                continue;

            // A[i][j] = par[2 * n * n + 10 + i * n + j];
            // beta_x[i][j] = par[10 + i * n + j];
            // beta_y[i][j] = par[10 + n * n + i * n + j];

            // double erf_x_pos = 0.5*(1+erf((x[0] - (beta_x[i][j] - l2)) / (sqrt(2) * sigma_x)));
            // double erf_x_neg = 0.5*(1+erf((x[0] - (beta_x[i][j] + l2)) / (sqrt(2) * sigma_x)));
            // double erf_x = (x[0] >= (beta_x[i][j] - l2) && x[0] <= (beta_x[i][j] + l2)) ? 1.0 : 0.0;

            // double erf_y_pos = 0.5*(1+erf((x[1] - (beta_y[i][j] - l2)) / (sqrt(2) * sigma_y)));
            // double erf_y_neg = 0.5*(1+erf((x[1] - (beta_y[i][j] + l2)) / (sqrt(2) * sigma_y)));
            // double erf_y = (x[1] >= (beta_y[i][j] - l2) && x[1] <= (beta_y[i][j] + l2)) ? 1.0 : 0.0;

            // result_sum += (erf_x_pos - erf_x_neg) * (erf_y_pos - erf_y_neg);
            // grid_sum += (erf_x_pos - erf_x_neg) * (erf_y_pos - erf_y_neg);

            double Ax = par[10 + i * n + j];
            double Ay = par[10 + n * n + i * n + j];
            double Bx = par[10 + 2 * n * n + i * n + j];
            // double By = par[10 + 3 * n * n + i * n + j];
            // double Cx = par[10 + 4 * n * n + i * n + j];
            // double Cy = par[10 + 5 * n * n + i * n + j];
            // double Dx = par[10 + 6 * n * n + i * n + j];
            double Dy = par[10 + 7 * n * n + i * n + j];

            double erf_x = (x[0] >= (Ax) && x[0] <= (Bx)) ? 1.0 : 0.0;
            double erf_y = (x[1] >= (Ay) && x[1] <= (Dy)) ? 1.0 : 0.0;

            result_sum += erf_x * erf_y;
            grid_sum += erf_x * erf_y;
        }
    }

    // RAW EXPRESSION  OF 2D GAUSSIAN
    // double theta = A * TMath::Pi() / 180.;
    // double Xr = (x[0]-mu_gx)*cos(theta) + (x[1]-mu_gy)*sin(theta);
    // double Yr = -(x[0]-mu_gx)*sin(theta) + (x[1]-mu_gy)*cos(theta);

    // double rho = 0.0;
    // double gaussr = A_g / (2*M_PI*sigma_gx*sigma_gy) * exp(-0.5 * ((Xr) * (Xr) / (sigma_gx * sigma_gx) + (Yr) * (Yr) / (sigma_gy * sigma_gy)));
    // double gaussr_corr = bkg + A_g * exp(-0.5 *1./(1-rho*rho) * ((Xr) * (Xr) / (sigma_gx * sigma_gx) + (Yr) * (Yr) / (sigma_gy * sigma_gy)) + 2 * rho * Xr * Yr / (sigma_gx * sigma_gy));

    double par_gauss[5] = {A_g, mu_gx, sigma_gx, mu_gy, sigma_gy};
    double gauss = MyGaussian(x, par_gauss);
    // double gauss = A_g / (2*M_PI*sigma_gx*sigma_gy) * exp(-0.5 * ((x[0] - mu_gx) * (x[0] - mu_gx) / (sigma_gx * sigma_gx) + (x[1] - mu_gy) * (x[1] - mu_gy) / (sigma_gy * sigma_gy)));
    
    result_sum *= gauss;

    return result_sum + bkg*grid_sum;
}

double f(double *x, double *par)
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
    double l2 = par[9]/2;
    double beta_x[n][n];
    double beta_y[n][n];
    // double A[n][n];
    double a[6] = {par[3*n*n+10], par[3*n*n+11], par[3*n*n+21], par[3*n*n+13], par[3*n*n+14], par[3*n*n+15]};

    double result_sum = 0.0;
    double grid_sum = 0.0;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {

            if (!M_measurement[i][j])
                continue;

            // A[i][j] = par[2 * n * n + 10 + i * n + j];
            beta_x[i][j] = par[10 + i * n + j];
            beta_y[i][j] = par[10 + n * n + i * n + j];

            // double erf_x_pos = 0.5*(1+erf((x[0] - (beta_x[i][j] - l2)) / (sqrt(2) * sigma_x)));
            // double erf_x_neg = 0.5*(1+erf((x[0] - (beta_x[i][j] + l2)) / (sqrt(2) * sigma_x)));
            double erf_x = (x[0] >= (beta_x[i][j] - l2) && x[0] <= (beta_x[i][j] + l2)) ? 1.0 : 0.0;

            // double erf_y_pos = 0.5*(1+erf((x[1] - (beta_y[i][j] - l2)) / (sqrt(2) * sigma_y)));
            // double erf_y_neg = 0.5*(1+erf((x[1] - (beta_y[i][j] + l2)) / (sqrt(2) * sigma_y)));
            double erf_y = (x[1] >= (beta_y[i][j] - l2) && x[1] <= (beta_y[i][j] + l2)) ? 1.0 : 0.0;

            // result_sum += (erf_x_pos - erf_x_neg) * (erf_y_pos - erf_y_neg);
            // grid_sum += (erf_x_pos - erf_x_neg) * (erf_y_pos - erf_y_neg);

            result_sum += erf_x * erf_y;
            grid_sum += erf_x * erf_y;
        }
    }

    // RAW EXPRESSION  OF 2D GAUSSIAN
    // double theta = A * TMath::Pi() / 180.;
    // double Xr = (x[0]-mu_gx)*cos(theta) + (x[1]-mu_gy)*sin(theta);
    // double Yr = -(x[0]-mu_gx)*sin(theta) + (x[1]-mu_gy)*cos(theta);

    // double rho = 0.0;
    // double gaussr = A_g / (2*M_PI*sigma_gx*sigma_gy) * exp(-0.5 * ((Xr) * (Xr) / (sigma_gx * sigma_gx) + (Yr) * (Yr) / (sigma_gy * sigma_gy)));
    // double gaussr_corr = bkg + A_g * exp(-0.5 *1./(1-rho*rho) * ((Xr) * (Xr) / (sigma_gx * sigma_gx) + (Yr) * (Yr) / (sigma_gy * sigma_gy)) + 2 * rho * Xr * Yr / (sigma_gx * sigma_gy));

    double par_gauss[5] = {A_g, mu_gx, sigma_gx, mu_gy, sigma_gy};
    double gauss = MyGaussian(x, par_gauss);
    // double gauss = A_g / (2*M_PI*sigma_gx*sigma_gy) * exp(-0.5 * ((x[0] - mu_gx) * (x[0] - mu_gx) / (sigma_gx * sigma_gx) + (x[1] - mu_gy) * (x[1] - mu_gy) / (sigma_gy * sigma_gy)));
    
    result_sum *= gauss;

    return result_sum + bkg*grid_sum;
}

double MeasurementFittedFunction2D_Convoluted_CORNER(double *x, double *par)
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
    double l2 = par[9]/2;
    double beta_x[n][n];
    double beta_y[n][n];
    // double A[n][n];
    // double a[6] = {par[3*n*n+10], par[3*n*n+11], par[3*n*n+21], par[3*n*n+13], par[3*n*n+14], par[3*n*n+15]};

    // double result = f(x, par);
    /// CONVOLUTING result by a 2D GAUSSIAN
    double result_conv = 0.0;
    double par_g[5] = {1, 0, sigma_x, 0, sigma_y};

    // double loop for convolution gaussian 2d on value 
    int Step = 100;
    double fStep = (x_measurement_max - x_measurement_min) / Step;  
    for (int i = 0 ; i < Step; i++)
    {
        double tx = x_measurement_min + i * fStep;
        for (int j = 0 ; j < Step; j++)
        {
            double ty = x_measurement_min + j * fStep;
            double xx[2] = {tx, ty};
            double xxx[2] = {x[0]-tx, x[1]-ty};
            // result_conv += Gaussian->Eval(x[0]-tx, x[1]-ty) * f(xx, par) * fStep * fStep;
            result_conv += MyGaussian(xxx, par_g) * f_CORNER(xx, par) * fStep * fStep;
        }
    }

    return result_conv;
}


double MeasurementFittedFunction2D_Convoluted(double *x, double *par)
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
    double l2 = par[9]/2;
    double beta_x[n][n];
    double beta_y[n][n];
    // double A[n][n];
    // double a[6] = {par[3*n*n+10], par[3*n*n+11], par[3*n*n+21], par[3*n*n+13], par[3*n*n+14], par[3*n*n+15]};

    // double result = f(x, par);
    /// CONVOLUTING result by a 2D GAUSSIAN
    double result_conv = 0.0;
    double par_g[5] = {1, 0, sigma_x, 0, sigma_y};

    // double loop for convolution gaussian 2d on value 
    int Step = 100;
    double fStep = (x_measurement_max - x_measurement_min) / Step;  
    for (int i = 0 ; i < Step; i++)
    {
        double tx = x_measurement_min + i * fStep;
        for (int j = 0 ; j < Step; j++)
        {
            double ty = x_measurement_min + j * fStep;
            double xx[2] = {tx, ty};
            double xxx[2] = {x[0]-tx, x[1]-ty};
            // result_conv += Gaussian->Eval(x[0]-tx, x[1]-ty) * f(xx, par) * fStep * fStep;
            result_conv += MyGaussian(xxx, par_g) * f(xx, par) * fStep * fStep;
        }
    }

    return result_conv;
}

double MeasurementFunctionToMinimize2D_CORNER()
{

    ////////// SOLUTION OF FIT //////////////
    // ### ERROR AT THE END THE FILE ### // 
    double sigma_x = 1.70116e-01;
    double sigma_y = 1.28192e-01;
    double Amplitude_gauss = 2.62887e+02;
    double mu_gx = -7.91786e-02;
    double mu_gy = 6.31885e-02;
    double sigma_gx = 6.59898e-01;
    double sigma_gy = 7.73420e-01;
    double bkg = 3.40290e+00;
    /////////////////////////////////////////   

    double chi2 = 0.0;
    int N = n * n;

    // ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(2); // 2 for verbose
    ROOT::Math::MinimizerOptions::SetDefaultErrorDef(9.30*9.30);
    MeasurementFunction = new TF2("MeasurementFunction", MeasurementFittedFunction2D_Convoluted_CORNER, x_measurement_min, x_measurement_max, x_measurement_min, x_measurement_max, 8 * N + 10);
    MeasurementFunction->SetNpx(75);
    MeasurementFunction->SetNpy(75);

    // Amplitude
    // MeasurementFunction->SetParLimits(0, 0, 180);
    // MeasurementFunction->SetParameter(0, 50);
    MeasurementFunction->FixParameter(0, 0);

    // sigma x
    MeasurementFunction->SetParLimits(1, 0.1, 0.4);
    MeasurementFunction->SetParameter(1, 0.17);
    // MeasurementFunction->FixParameter(1, sigma_x);

    // sigma y
    MeasurementFunction->SetParLimits(2, 0.1, 0.4);
    MeasurementFunction->SetParameter(2, 0.12);
    // MeasurementFunction->FixParameter(2, sigma_y);

    // amplitude gaus
    MeasurementFunction->SetParLimits(3, 0, 1000);
    MeasurementFunction->SetParameter(3, 50);
    // MeasurementFunction->FixParameter(3, Amplitude_gauss);

    // mu gaus x
    MeasurementFunction->SetParLimits(4, -1, 1);
    MeasurementFunction->SetParameter(4, 0);
    // MeasurementFunction->FixParameter(4, mu_gx);

    // mu gaus y
    MeasurementFunction->SetParLimits(5, -1, 1);
    MeasurementFunction->SetParameter(5, 0);
    // MeasurementFunction->FixParameter(5, mu_gy);

    // sigma gaus x
    MeasurementFunction->SetParLimits(6, 0.1, 1.5);
    MeasurementFunction->SetParameter(6, 0.6);
    // MeasurementFunction->FixParameter(6, sigma_gx);

    // sigma gaus y
    MeasurementFunction->SetParLimits(7, 0.1, 1.5);
    MeasurementFunction->SetParameter(7, 0.6);
    // MeasurementFunction->FixParameter(7, sigma_gy);

    // bkg
    MeasurementFunction->SetParLimits(8, 0, 100);
    MeasurementFunction->SetParameter(8, 4);
    // MeasurementFunction->FixParameter(8, bkg);

    // l
    // MeasurementFunction->SetParLimits(9, 0.8, 1.4);
    MeasurementFunction->FixParameter(9, l_real);

    for (int i = 0; i < N; ++i)
    {

        if (!M_measurement[i / n][i % n])
        {
            // a_x
            MeasurementFunction->FixParameter(10 + i, 0.);
            // a_y
            MeasurementFunction->FixParameter(N + 10 + i, 0.);
            // b_x
            MeasurementFunction->FixParameter(2 * N + 10 + i, 0.);
            // b_y
            MeasurementFunction->FixParameter(3 * N + 10 + i, 0.);
            // c_x
            MeasurementFunction->FixParameter(4 * N + 10 + i, 0.);
            // c_y
            MeasurementFunction->FixParameter(5 * N + 10 + i, 0.);
            // d_x
            MeasurementFunction->FixParameter(6 * N + 10 + i, 0.);
            // d_y
            MeasurementFunction->FixParameter(7 * N + 10 + i, 0.);
        }
        else
        {
            // MeasurementFunction->FixParameter(i, 0);
            double x_center = i / n * rho_real - rho_real * (n - 1) / 2;
            double y_center = i % n * rho_real - rho_real * (n - 1) / 2;

            // a_x
            MeasurementFunction->FixParameter(10 + i, x_center - l_real / 2);

            // a_y
            MeasurementFunction->FixParameter(N + 10 + i, y_center - l_real / 2);

            // b_x
            MeasurementFunction->FixParameter(2 * N + 10 + i, x_center + l_real / 2);

            // b_y
            MeasurementFunction->FixParameter(3 * N + 10 + i, y_center - l_real / 2);

            // c_x
            MeasurementFunction->FixParameter(4 * N + 10 + i, x_center + l_real / 2);

            // c_y
            MeasurementFunction->FixParameter(5 * N + 10 + i, y_center + l_real / 2);

            // d_x
            MeasurementFunction->FixParameter(6 * N + 10 + i, x_center - l_real / 2);

            // d_y
            MeasurementFunction->FixParameter(7 * N + 10 + i, y_center + l_real / 2);
            
        }
    }    

    
    H_measurement_reconstructed->Rebin2D(2, 2);

    H_measurement_reconstructed->GetXaxis()->SetRangeUser(x_measurement_min, x_measurement_max);
    H_measurement_reconstructed->GetYaxis()->SetRangeUser(x_measurement_min, x_measurement_max);

    // H_measurement_reconstructed->SetBinContent(57, 51, 0);
    H_measurement_reconstructed->Fit(MeasurementFunction, "MULTITHREAD RNM");
    
    chi2 = H_measurement_reconstructed->Chisquare(MeasurementFunction) / MeasurementFunction->GetNDF();
    cout << "chi2 = " << chi2 << endl;

    

    if (WRITTING)
    {
        FINAL_file->cd();

        MeasurementFunction->Write();

        TCanvas *c1 = new TCanvas("MeassurementFitted_2D_View", "MeassurementFitted_2D_View", 800, 800);
        H_measurement_reconstructed->Draw("COLZ");
        MeasurementFunction->Draw("SAME");
        c1->Write();

        TCanvas *c = new TCanvas("MeassurementFitted_3D_View", "MeassurementFitted_3D_View", 800, 800);
        c->Divide(2, 1);
        c->cd(1);
        H_measurement_reconstructed->GetXaxis()->SetRangeUser(x_measurement_min, x_measurement_max);
        H_measurement_reconstructed->GetYaxis()->SetRangeUser(x_measurement_min, x_measurement_max);
        MeasurementFunction->Draw("SURF");
        c->cd(2);
        // MeasurementFunction->SetNpx(75);
        // MeasurementFunction->SetNpy(75);
        H_measurement_reconstructed->Draw("LEGO2");
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
    return chi2;
}

double MeasurementFunctionToMinimize2D()
{
    // double sigma_x = par[0];
    // double sigma_y = par[1];
    // double Amplitude_gauss = par[2];
    // double mu_gx = par[3];
    // double mu_gy = par[4];
    // double sigma_gx = par[5];
    // double sigma_gy = par[6];
    // double bkg = par[7];


    ////////// SOLUTION OF FIT //////////////
    // ### ERROR AT THE END THE FILE ### // 
    double sigma_x = 1.70116e-01;
    double sigma_y = 1.28192e-01;
    double Amplitude_gauss = 2.62887e+02;
    double mu_gx = -7.91786e-02;
    double mu_gy = 6.31885e-02;
    double sigma_gx = 6.59898e-01;
    double sigma_gy = 7.73420e-01;
    double bkg = 3.40290e+00;
    /////////////////////////////////////////   

    double chi2 = 0.0;
    int N = n * n;

    // ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(2); // 2 for verbose
    ROOT::Math::MinimizerOptions::SetDefaultErrorDef(9.30*9.30);
    MeasurementFunction = new TF2("MeasurementFunction", MeasurementFittedFunction2D_Convoluted, x_measurement_min, x_measurement_max, x_measurement_min, x_measurement_max, 3 * N + 10);
    MeasurementFunction->SetNpx(75);
    MeasurementFunction->SetNpy(75);

    // Amplitude
    // MeasurementFunction->SetParLimits(0, 0, 180);
    // MeasurementFunction->SetParameter(0, 50);
    MeasurementFunction->FixParameter(0, 0);

    // sigma x
    MeasurementFunction->SetParLimits(1, 0.1, 0.4);
    MeasurementFunction->SetParameter(1, 0.17);
    MeasurementFunction->FixParameter(1, sigma_x);

    // sigma y
    MeasurementFunction->SetParLimits(2, 0.1, 0.4);
    MeasurementFunction->SetParameter(2, 0.12);
    MeasurementFunction->FixParameter(2, sigma_y);

    // amplitude gaus
    MeasurementFunction->SetParLimits(3, 0, 1000);
    MeasurementFunction->SetParameter(3, 50);
    MeasurementFunction->FixParameter(3, Amplitude_gauss);

    // mu gaus x
    MeasurementFunction->SetParLimits(4, -1, 1);
    MeasurementFunction->SetParameter(4, 0);
    MeasurementFunction->FixParameter(4, mu_gx);

    // mu gaus y
    MeasurementFunction->SetParLimits(5, -1, 1);
    MeasurementFunction->SetParameter(5, 0);
    MeasurementFunction->FixParameter(5, mu_gy);

    // sigma gaus x
    MeasurementFunction->SetParLimits(6, 0.1, 1.5);
    MeasurementFunction->SetParameter(6, 0.6);
    MeasurementFunction->FixParameter(6, sigma_gx);

    // sigma gaus y
    MeasurementFunction->SetParLimits(7, 0.1, 1.5);
    MeasurementFunction->SetParameter(7, 0.6);
    MeasurementFunction->FixParameter(7, sigma_gy);

    // bkg
    MeasurementFunction->SetParLimits(8, 0, 100);
    MeasurementFunction->SetParameter(8, 4);
    MeasurementFunction->FixParameter(8, bkg);

    // l
    // MeasurementFunction->SetParLimits(9, 0.8, 1.4);
    MeasurementFunction->FixParameter(9, l_real);

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
            double x_center = i / n * rho_real - rho_real * (n - 1) / 2;
            double y_center = i % n * rho_real - rho_real * (n - 1) / 2;

            // beta_x
            MeasurementFunction->FixParameter(10 + i, x_center);

            // beta_y
            MeasurementFunction->FixParameter(N + 10 + i, y_center);

            // amp
            // MeasurementFunction->SetParLimits(2 * N + 10 + i, 0, 200);
            MeasurementFunction->FixParameter(2*N + 10 + i, 0.);
        }
    }    

    H_measurement_reconstructed->Rebin2D(2, 2);

    H_measurement_reconstructed->GetXaxis()->SetRangeUser(x_measurement_min, x_measurement_max);
    H_measurement_reconstructed->GetYaxis()->SetRangeUser(x_measurement_min, x_measurement_max);

    H_measurement_reconstructed->SetBinContent(57, 51, 0);
    H_measurement_reconstructed->Fit(MeasurementFunction, "MULTITHREAD RNM");
    
    chi2 = H_measurement_reconstructed->Chisquare(MeasurementFunction) / MeasurementFunction->GetNDF();
    cout << "chi2 = " << chi2 << endl;

    

    if (WRITTING)
    {
        FINAL_file->cd();

        MeasurementFunction->Write();

        TCanvas *c1 = new TCanvas("MeassurementFitted_2D_View", "MeassurementFitted_2D_View", 800, 800);
        H_measurement_reconstructed->Draw("COLZ");
        MeasurementFunction->Draw("SAME");
        c1->Write();

        TCanvas *c = new TCanvas("MeassurementFitted_3D_View", "MeassurementFitted_3D_View", 800, 800);
        c->Divide(2, 1);
        c->cd(1);
        H_measurement_reconstructed->GetXaxis()->SetRangeUser(x_measurement_min, x_measurement_max);
        H_measurement_reconstructed->GetYaxis()->SetRangeUser(x_measurement_min, x_measurement_max);
        MeasurementFunction->Draw("SURF");
        c->cd(2);
        // MeasurementFunction->SetNpx(75);
        // MeasurementFunction->SetNpy(75);
        H_measurement_reconstructed->Draw("LEGO2");
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
    return chi2;
}

void Measurement2D()
{
    FINAL_file->cd();
    if (YEAR == 2024)
    {
        // recreate data in histogram
        H_measurement_reconstructed = new TH2D("H_measurement_reconstructed", "H_measurement_reconstructed", 1000, x_measurement_min, x_measurement_max, 1000, x_measurement_min, x_measurement_max);
        for (int i = 0; i < H_measurement->GetEntries(); ++i)
        {
            double x, y;
            H_measurement->GetRandom2(x, y);
            double x_fit = f_X->Eval(x, y) - 0.27;
            double y_fit = f_Y->Eval(y, x) + 0.30;
            H_measurement_reconstructed->Fill(x_fit, y_fit);
        }

        // H_measurement_reconstructed->SetBinContent(57, 51, 0);
    }
    else
    {
        // recreate data in histogram from Tree
        H_measurement_reconstructed = new TH2D("H_measurement_reconstructed", "H_measurement_reconstructed", 300, x_measurement_min, x_measurement_max, 300, x_measurement_min, x_measurement_max);
        TFile *file_measurement = MyTFile((DIR_ROOT_DATA_MCP+Measurement_Filename).c_str(), "READ");

        TTree *tree = (TTree*)file_measurement->Get("treeMCP");
        TTreeReader *Reader = new TTreeReader(tree);
        TTreeReaderValue<double> *X_Tree = new TTreeReaderValue<double>(*Reader, "X");
        TTreeReaderValue<double> *Y_Tree = new TTreeReaderValue<double>(*Reader, "Y");

        while (Reader->Next())
        {
            double x = **X_Tree;
            double y = **Y_Tree;

            double x_fit = f_X->Eval(x, y);
            double y_fit = f_Y->Eval(y, x);

            H_measurement_reconstructed->Fill(x_fit, y_fit);
        }
        file_measurement->Close();
    }
    FINAL_file->cd();

    TCanvas *c_measurement_reconstructed = new TCanvas("Measurement_2D_View", "Measurement_2D_View", 800, 800);
    H_measurement_reconstructed->Draw("COLZ");
    c_measurement_reconstructed->Write();

    fSaved = new TFile("Calibration_Saved.root", "RECREATE");
    fSaved->cd();
    H_measurement_reconstructed->SetName("H_measurement_reconstructed");
    H_measurement_reconstructed->Write();
    fSaved->Close();

    FINAL_file->cd();
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




// 1 SIGMA 8 PARAMETERS

// EXTERNAL ERROR MATRIX.    NDIM= 203    NPAR=  8    ERR DEF=86.49
//   1.090e-03  1.717e-04  1.769e-02 -1.672e-05  2.030e-04  2.956e-04  3.294e-04 -7.549e-03 
//   1.717e-04  1.065e-03 -1.146e-02 -3.404e-04  3.019e-04  6.300e-04  1.234e-04 -1.506e-03 
//   1.769e-02 -1.146e-02  4.445e+01  9.188e-02 -3.663e-02  4.397e-01  4.347e-01 -1.731e+01 
//  -1.672e-05 -3.404e-04  9.188e-02  3.317e-03 -7.914e-04  2.142e-03  1.888e-03 -6.069e-02 
//   2.030e-04  3.019e-04 -3.663e-02 -7.914e-04  5.406e-03  2.339e-03 -2.235e-04 -1.686e-02 
//   2.956e-04  6.300e-04  4.397e-01  2.142e-03  2.339e-03  2.649e-02  1.675e-02 -5.411e-01 
//   3.294e-04  1.234e-04  4.347e-01  1.888e-03 -2.235e-04  1.675e-02  2.039e-02 -4.447e-01 
//  -7.549e-03 -1.506e-03 -1.731e+01 -6.069e-02 -1.686e-02 -5.411e-01 -4.447e-01  1.460e+01 
//  PARAMETER  CORRELATION COEFFICIENTS  
//        NO.  GLOBAL      2      3      4      5      6      7      8      9
//         2  0.20601   1.000  0.159  0.080 -0.009  0.084  0.055  0.070 -0.060
//         3  0.32994   0.159  1.000 -0.053 -0.181  0.126  0.119  0.026 -0.012
//         4  0.79627   0.080 -0.053  1.000  0.239 -0.075  0.405  0.457 -0.680
//         5  0.39200  -0.009 -0.181  0.239  1.000 -0.187  0.229  0.230 -0.276
//         6  0.39340   0.084  0.126 -0.075 -0.187  1.000  0.196 -0.021 -0.060
//         7  0.91705   0.055  0.119  0.405  0.229  0.196  1.000  0.721 -0.870
//         8  0.83256   0.070  0.026  0.457  0.230 -0.021  0.721  1.000 -0.815
//         9  0.95979  -0.060 -0.012 -0.680 -0.276 -0.060 -0.870 -0.815  1.000
//  FCN=5943.79 FROM HESSE     STATUS=OK             61 CALLS         594 TOTAL
//                      EDM=1.7357e-06    STRATEGY= 1      ERROR MATRIX ACCURATE 
//   EXT PARAMETER                                   STEP         FIRST   
//   NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
//    1  p0           0.00000e+00     fixed    
//    2  p1           1.70116e-01   3.26478e-02   2.07529e-04   1.67458e-02
//    3  p2           1.28192e-01   3.18885e-02   2.87029e-04  -3.31911e-02
//    4  p3           4.18400e+01   6.66552e+00   3.28492e-06  -5.75860e-02
//    5  p4          -7.91798e-02   5.75601e-02   4.33378e-05  -5.69671e-02
//    6  p5           6.31906e-02   7.34606e-02   5.52322e-05   2.39150e-02
//    7  p6           6.59905e-01   1.61221e-01   7.71534e-05   3.67518e-02
//    8  p7           7.73425e-01   1.41796e-01   9.22056e-05  -6.31119e-02
//    9  p8           3.40275e+00   3.79227e+00   9.64843e-06  -6.75542e-02

//// 2 SIGMA 8 PARAMETERS

//                             ERR DEF= 249.64
//  EXTERNAL ERROR MATRIX.    NDIM= 203    NPAR=  8    ERR DEF=249.64
//   3.147e-03  4.955e-04  5.105e-02 -4.826e-05  5.861e-04  8.533e-04  9.507e-04 -2.179e-02 
//   4.955e-04  3.075e-03 -3.308e-02 -9.825e-04  8.715e-04  1.818e-03  3.561e-04 -4.347e-03 
//   5.105e-02 -3.308e-02  1.283e+02  2.652e-01 -1.057e-01  1.269e+00  1.255e+00 -4.996e+01 
//  -4.826e-05 -9.825e-04  2.652e-01  9.574e-03 -2.284e-03  6.183e-03  5.449e-03 -1.752e-01 
//   5.861e-04  8.715e-04 -1.057e-01 -2.284e-03  1.560e-02  6.752e-03 -6.449e-04 -4.865e-02 
//   8.533e-04  1.818e-03  1.269e+00  6.183e-03  6.752e-03  7.645e-02  4.836e-02 -1.562e+00 
//   9.507e-04  3.561e-04  1.255e+00  5.449e-03 -6.449e-04  4.836e-02  5.885e-02 -1.283e+00 
//  -2.179e-02 -4.347e-03 -4.996e+01 -1.752e-01 -4.865e-02 -1.562e+00 -1.283e+00  4.213e+01 
//  PARAMETER  CORRELATION COEFFICIENTS  
//        NO.  GLOBAL      2      3      4      5      6      7      8      9
//         2  0.20601   1.000  0.159  0.080 -0.009  0.084  0.055  0.070 -0.060
//         3  0.32994   0.159  1.000 -0.053 -0.181  0.126  0.119  0.026 -0.012
//         4  0.79627   0.080 -0.053  1.000  0.239 -0.075  0.405  0.457 -0.680
//         5  0.39199  -0.009 -0.181  0.239  1.000 -0.187  0.229  0.230 -0.276
//         6  0.39340   0.084  0.126 -0.075 -0.187  1.000  0.196 -0.021 -0.060
//         7  0.91705   0.055  0.119  0.405  0.229  0.196  1.000  0.721 -0.870
//         8  0.83256   0.070  0.026  0.457  0.230 -0.021  0.721  1.000 -0.815
//         9  0.95979  -0.060 -0.012 -0.680 -0.276 -0.060 -0.870 -0.815  1.000
//  FCN=5943.79 FROM HESSE     STATUS=OK             61 CALLS         609 TOTAL
//                      EDM=1.74129e-06    STRATEGY= 1      ERROR MATRIX ACCURATE 
//   EXT PARAMETER                                   STEP         FIRST   
//   NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
//    1  p0           0.00000e+00     fixed    
//    2  p1           1.70116e-01   5.42888e-02   2.10318e-04   1.67659e-02
//    3  p2           1.28192e-01   5.18172e-02   2.90886e-04  -3.32327e-02
//    4  p3           4.18400e+01   1.13203e+01   3.32906e-06  -5.77193e-02
//    5  p4          -7.91798e-02   9.76877e-02   4.39202e-05  -5.68426e-02
//    6  p5           6.31906e-02   1.24591e-01   5.59744e-05   2.37812e-02
//    7  p6           6.59905e-01   2.69061e-01   7.81901e-05   3.69810e-02
//    8  p7           7.73425e-01   2.37748e-01   9.34446e-05  -6.31110e-02
//    9  p8           3.40275e+00   6.35295e+00   9.77808e-06  -6.74828e-02
