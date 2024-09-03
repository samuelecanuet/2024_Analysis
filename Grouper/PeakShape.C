#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <string>
using namespace std;
double f(double *x, double *p)
{
    double A = p[0];
    double mu = p[1]/2;
    double sigma = p[2];
    double sigma_gate = p[3];

    double res = A * sigma * ( erf ( (4*mu + sigma_gate - 2*x[0])/(2*sqrt(2)*sigma) ) - erf((4*mu - sigma_gate - 2*x[0])/(2*sqrt(2)*sigma)  ) );
    return res;
}

double ff (double *x, double *par)
{
    double A = par[4];
    double mu = par[1];
    double sigma = par[5];

    double par1[4]  = {par[0], par[1], par[2], par[3]};

    return 1 * f(x, par1);
} 

void PeakShape()
{
    TFile* file = new TFile("../../../../../../mnt/hgfs/shared-2/2024_DATA/SIMULATED_DATA/32Ar__CS0_CSP0_CV1_CVP1.root", "READ");
    TH1D* h = (TH1D*)file->Get("D5.1_single");
    h->GetXaxis()->SetRangeUser(3300, 3380);




    TF1 *function = new TF1("test", ff, 0, 10000, 6);
    function->SetNpx(10000);
    function->SetParameters(1000, 3350, 5, 50);
    function->SetParLimits(0, 0, 10000);
    function->SetParLimits(1, 3320, 3380);
    function->SetParLimits(2, 1, 1000);
    function->SetParLimits(3, 1, 1000);
    function->SetParLimits(4, 0, 1000);
    function->SetParLimits(5, 0, 1000);

    h->Fit("test", "R", "", 3300, 3380);
       

    // Draw the convoluted function
    TCanvas *c1 = new TCanvas("c1", "Convolution", 800, 600);
    h->Draw("HIST");
    function->Draw("SAME");
}

