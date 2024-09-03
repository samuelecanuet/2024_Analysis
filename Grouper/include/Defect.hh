


#ifndef DEFECT_HH
#define DEFECT_HH

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

#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TMatrixD.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
using namespace std;
using namespace ROOT::Math;
default_random_engine generator;

/// FOLDERS ///
string DIR_FAST_DATA = "../../../../../run/media/local1/Disque_Dur/2024_DATA/DETECTOR_DATA/DATA/";
string DIR_ROOT_DATA = "../../../../../mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/ROOT/";
string DIR_ROOT_DATA_GROUPED = "../../../../../mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/GROUPED/";
string DIR_ROOT_DATA_CLEANED = "../../../../../mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/CLEANED/";
string DIR_ROOT_DATA_MATCHED = "../../../../../mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/MATCHED/";
string DIR_ROOT_DATA_MERGED  = "../../../../../mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/MERGED/";
string DIR_ROOT_DATA_SIMULATED = "../../../../../mnt/hgfs/shared-2/2024_DATA/SIMULATED_DATA/";
///////////////



map<string, TFile*> GROUPED_File;
map<string, TFile*> SIMULATED_File;
TFile *DEFECT_File;

// map<string, pair<double, double>> WindowExp = {
//     {"4alpha", {0, 100000}},
//     {"32Ar", {0, 100000}},
// };
map<string, pair<double, double>> WindowExp = {
    {"32Ar", {0, 100000}},
    {"148Gd", {43000, 45000}},
    {"239Pu", {70000, 74000}},
    {"241Am", {74000, 78000}},
    {"244Cm", {79000, 82000}},
};  

//D1.1
// map<string, vector<int>> Peaks_Exp = {
//     {"18N", {11800, 16750}},
//     {"148Gd", {44500}},
//     {"239Pu", {71250, 72000, 72000}},
//     {"241Am", {75250, 76000, 76750}},
//     {"244Cm", {80400, 81100}},
// };

// map<string, vector<int>> Peaks_Sim = {
//     {"18N", {670, 1058}},
//     {"148Gd", {3163}},
//     {"239Pu", {5091, 5142, 5129}},
//     {"241Am", {5374, 5429, 5472}},
//     {"244Cm", {5749, 5792}},
// };

//D5.1
map<string, vector<int>> Peaks_Exp = {
    {"18N", {14800, 19500}},
    {"148Gd", {44750}},
    {"239Pu", {71850, 72500, 72550}},
    {"241Am", {75750, 76500, 77250}},
    {"244Cm", {81000, 81600}},
};

map<string, vector<int>> Peaks_Sim = {
    {"18N", {1050, 1350}},
    {"148Gd", {3163}},
    {"239Pu", {5095, 5130, 5150}},
    {"241Am", {5374, 5429, 5472}},
    {"244Cm", {5749, 5792}},
};

double gauss(double *x, double *p)
{
    double A = p[0];
    double mu = p[1];
    double sigma = p[2];
    return A * exp(-pow(x[0]-mu, 2)/(2*pow(sigma, 2)));
}

double doublegauss(double *x, double *p)
{
    double p1[3] = {p[0], p[1], p[2]};
    double p2[3] = {p[3], p[4], p[2]};
    return gauss(x, p1) + gauss(x, p2);
}

double triplegauss(double *x, double *p)
{
    double p1[3] = {p[0], p[1], p[2]};
    double p2[3] = {p[3], p[4], p[2]};
    double p3[3] = {p[5], p[6], p[2]};
    return gauss(x, p1) + gauss(x, p2) + gauss(x, p3);
}

double gaussBortelsCollaersRight(double *x, double *p)
{
    double A = p[0];
    double mu = p[1];
    double sigma = p[2]; 
    double eta = p[3];
    double tau1 = p[4];
    double tau2 = p[5];

    double first  = (1 - eta) / tau1    * exp((x[0]-mu)/tau1 + pow(sigma, 2)/(2*pow(tau1, 2))) * 0.5*(1+erf(1/sqrt(2) * (x[0] - mu)/sigma + sigma/tau1));
    double second =      eta  / tau2    * exp((x[0]-mu)/tau2 + pow(sigma, 2)/(2*pow(tau2, 2))) * 0.5*(1+erf(1/sqrt(2) * (x[0] - mu)/sigma + sigma/tau2));
    // double second = 0;

    return A/2 * (first+second);
}

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

    return A/2 * (first+second);
}

double doublegaussBortelsCollaers(double *x, double *p)
{
    double p1[6] = {p[0], p[1], p[2], p[3], p[4], p[5]};
    double p2[6] = {p[6], p[7], p[2], p[3], p[4], p[5]};
    return gaussBortelsCollaers(x, p1) + gaussBortelsCollaers(x, p2);
}

double triplegaussBortelsCollaers(double *x, double *p)
{
    double p1[6] = {p[0], p[1], p[2], p[3], p[4], p[5]};
    double p2[6] = {p[6], p[7], p[2], p[3], p[4], p[5]};
    double p3[6] = {p[8], p[9], p[2], p[3], p[4], p[5]};
    return gaussBortelsCollaers(x, p1) + gaussBortelsCollaers(x, p2) + gaussBortelsCollaers(x, p3);
}

double fullgaussBortelsCollaers(double *x, double *p)
{
    // nine fucntions
    double p1[6] = {p[0], p[1], p[2], p[3], p[4], p[5]};
    double p2[6] = {p[6], p[7], p[2], p[3], p[4], p[5]};
    double p3[6] = {p[8], p[9], p[2], p[3], p[4], p[5]};
    double p4[6] = {p[10], p[11], p[2], p[3], p[4], p[5]};
    double p5[6] = {p[12], p[13], p[2], p[3], p[4], p[5]};
    double p6[6] = {p[14], p[15], p[2], p[3], p[4], p[5]};
    double p7[6] = {p[16], p[17], p[2], p[3], p[4], p[5]};
    double p8[6] = {p[18], p[19], p[2], p[3], p[4], p[5]};
    double p9[6] = {p[20], p[21], p[2], p[3], p[4], p[5]};

    return gaussBortelsCollaers(x, p1) + gaussBortelsCollaers(x, p2) + gaussBortelsCollaers(x, p3) + gaussBortelsCollaers(x, p4) + gaussBortelsCollaers(x, p5) + gaussBortelsCollaers(x, p6) + gaussBortelsCollaers(x, p7) + gaussBortelsCollaers(x, p8) + gaussBortelsCollaers(x, p9);
}


#endif