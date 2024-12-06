#ifndef __FULLFUNCTION_HH__
#define __FULLFUNCTION_HH__


#include <iostream>
#include <cmath>
#include <vector>
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"   
#include "TVirtualFFT.h"
#include "TF1Convolution.h" 
#include "TCanvas.h"
#include "TGraph.h"

#include "Constants.hh"
#include "Peak.hh"
#include "Addition.hh"  
#include "Convolution.hh"
#include "Spectrum.hh"


// RESOLUTION //
double fGauss(double *x, double *par)
{
    double E = x[0];
    double E_p = 0;
    double A = par[0];
    double Sigma = par[1];

    return A * exp(-0.5 * pow((E - E_p) / Sigma, 2)) / (Sigma * sqrt(2 * M_PI));
}

double fGaussBortelsCollaers(double *x, double *p)
{
    double A = p[0];
    double mu = p[1];
    double sigma = p[2]; 
    double eta = p[3];
    double tau1 = p[4];
    double tau2 = p[5];

    double first  = (1 - eta) / tau1    * exp((x[0]-mu)/tau1 + pow(sigma, 2)/(2*pow(tau1, 2))) * erfc(1/sqrt(2) * (x[0] - mu)/sigma + sigma/tau1);
    double second =      eta  / tau2    * exp((x[0]-mu)/tau2 + pow(sigma, 2)/(2*pow(tau2, 2))) * erfc(1/sqrt(2) * (x[0] - mu)/sigma + sigma/tau2);

    return A/2 * (first+second);
}

#endif