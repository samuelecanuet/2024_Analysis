#include "Detectors.hh"
#include "Math/Minimizer.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <map>

#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TProfile.h"
#include "TLine.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TRandom.h"
#include "TLatex.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TSpline.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TTree.h"
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
#include "TVirtualFitter.h"

TFile * MERGED_File;
    TFile * SIMULATED_File;
    TFile * TEST_File;

    TH1D* H_Exp;
    TH1D* H_Sim;
    TH1D* H_Sim_Conv;
    TH1D* H_Sim_Conv2;
    TF1 *fBC;
using namespace ROOT::Math;

double gaussBortelsCollaers(double *x, double *p)
{
    double A = p[0];
    double mu = p[1];
    double sigma = p[2]; 
    double eta = p[3];
    double tau1 = p[4];
    double tau2 = p[5];

    double first  = (1 - eta) / tau1    * exp((x[0]-mu)/tau1 + pow(sigma, 2)/(2*pow(tau1, 2))) * erfc(1/sqrt(2) * (x[0] - mu)/sigma + sigma/tau1);
    double second =      eta  / tau2    * exp((x[0]-mu)/tau2 + pow(sigma, 2)/(2*pow(tau2, 2))) * erfc(1/sqrt(2) * (x[0] - mu)/sigma + sigma/tau2);
    // double second = 0;
    // return erfc(x[0]-mu);
    return A/2 * (first+second);
}

double gauss(double *x, double *p)
{
    double A = p[0];
    double mu = p[1];
    double sigma = p[2];
    return A * exp(-pow(x[0]-mu, 2)/(2*pow(sigma, 2)));
}

double FunctionToMinimize(const double *par)
{
    H_Sim_Conv->Reset();
    fBC->SetParameters(par);
    for (int i = 0; i < H_Sim->GetNbinsX(); i++)
    {
        double value = H_Sim->GetBinCenter(i);
        if (value < 2000 || value > 3500) continue;
        double content = H_Sim->GetBinContent(i);
        fBC->SetParameter(1, value+par[1]);
        H_Sim_Conv->Add(fBC->GetHistogram(), content);
    }

    H_Exp->GetXaxis()->SetRangeUser(3275, 3375);
    H_Sim_Conv->GetXaxis()->SetRangeUser(3275, 3375);

    double chi2 = H_Exp->Chi2Test(H_Sim_Conv, "UW CHI2/NDF");
    H_Sim_Conv->Scale(H_Exp->Integral() / H_Sim_Conv->Integral());

    cout << "Chi2: " << chi2 << endl;
    cout << "A: " << par[0] << " mu: " << par[1] << " sigma: " << par[2] << " eta: " << par[3] << " tau1: " << par[4] << " tau2: " << par[5] << endl;


    return chi2;
}

void Chi2Mini()
{

    Minimizer *minimizer = Factory::CreateMinimizer("Minuit2", "Migrad");
    ROOT::Math::Functor functor(&FunctionToMinimize, 6);
    minimizer->SetFunction(functor);
    // minimizer->SetStrategy(2);
    minimizer->SetLimitedVariable(0, "A", 2e6, 1000, 9e5, 5e6);
    minimizer->SetLimitedVariable(1, "mu", 0, 1, -15, 15);
    minimizer->SetLimitedVariable(2, "sigma", 0.1, 0.01, 0., 2);
    minimizer->SetLimitedVariable(3, "eta", 0.008, 0.001, 0, 0.1);
    minimizer->SetLimitedVariable(4, "tau1", 4, 1, 1, 10);
    minimizer->SetLimitedVariable(5, "tau2", 60, 1, 10, 200);
    minimizer->SetPrecision(0.001);
    minimizer->SetTolerance(0.001);
    minimizer->SetMaxFunctionCalls(1000000);
    minimizer->SetMaxIterations(1000000);

    // minimizer->Minimize();

    // const double *par = minimizer->X();
    // cout << "A: " << par[0] << " mu: " << par[1] << " sigma: " << par[2] << " eta: " << par[3] << " tau1: " << par[4] << " tau2: " << par[5] << endl;
    // FunctionToMinimize(new double[6]{par[0], par[1], par[2], par[3], par[4], par[5]});
    FunctionToMinimize(new double[6]{900124, 4, 2, 0.008, 4, 200});

    
}
