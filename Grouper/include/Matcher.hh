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

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
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

#include "../../../lib/SignalDict/Signal.h"
#include "Detectors.hh"

using namespace std;
using namespace ROOT::Math;
default_random_engine generator;

/// FILE ///
string GROUPED_filename;
string GROUPED_basefilename;
TFile *GROUPED_File;
string REFERENCE_filename;
TFile *REFERENCE_File;
TFile *MATCHED_File;

int Run;
int REFERENCE_Run = 77;

TTree *GROUPED_Tree_Detectors[SIGNAL_MAX];
TTreeReader *Reader;
TTreeReaderArray<Signal> *Silicon;
double Channel;

//HISTOGRAMS
TH1D *H_Run[SIGNAL_MAX];
TH1D *H_Run_Ref[SIGNAL_MAX];
TH1D *H[SIGNAL_MAX];
TGraph *G[SIGNAL_MAX];

int current_detector;
double CHI2;
double parameter[2];

void InitHistograms()
{
    for (int i = 0; i < detectorNum; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            H_Run[i] = new TH1D(("H_Run_" + detectorName[i]).c_str(), ("H_Run_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
            H_Run[i]->GetXaxis()->SetTitle("Channel");
            H_Run[i]->GetYaxis()->SetTitle("Counts");
            H_Run[i]->GetXaxis()->CenterTitle();
            H_Run[i]->GetYaxis()->CenterTitle();

            H[i] = new TH1D(("H_" + detectorName[i]).c_str(), ("H_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
            H[i]->GetXaxis()->SetTitle("Reference");
            H[i]->GetYaxis()->SetTitle("Run");
            H[i]->GetXaxis()->CenterTitle();
            H[i]->GetYaxis()->CenterTitle();

            G[i] = new TGraph();
            
        }
    }
}

double FunctionToMinimize(const double *par)
{

    // cout << setprecision(5) << "PAR : " << par[0] << " " << par[1] << endl;
    
    Reader->Restart();
    H_Run[current_detector]->Reset();
    TTreeReaderValue<double> ChannelDet(*Reader, "Channel");
    while (Reader->Next())
    {
        H_Run[current_detector]->Fill(par[0] + *ChannelDet * par[1] );
    }

    double chi2 = 0;
    chi2 = H_Run[current_detector]->Chi2Test(H_Run_Ref[current_detector], "CHI2/NDF");

    if (CHI2 > chi2)
    {
        CHI2 = chi2;
        parameter[0] = par[0];
        parameter[1] = par[1];
    }

    // cout << "Chi2 : " << chi2 << endl;
    // cout << endl;
    
    return chi2;
}

void CHI2Minimization(int i)
{
    current_detector = i;

    Minimizer *minimizer = Factory::CreateMinimizer("Minuit2", "Migrad");
    ROOT::Math::Functor functor(&FunctionToMinimize, 2);
    minimizer->SetFunction(functor);
    minimizer->SetFixedVariable(0, "Offset", 0);
    minimizer->SetFixedVariable(1, "Proportionnality", 1);
    minimizer->SetPrecision(0.00000001);
    minimizer->SetTolerance(0.00000001);
    minimizer->SetMaxFunctionCalls(1000000);
    minimizer->SetMaxIterations(1000000);

    minimizer->Minimize();
    const double *par = minimizer->X();
    
    // const double par[2] = {0, 1.0};
    // cout << par[0] << " " << par[1] << endl;
    FunctionToMinimize(par);

    TCanvas *c = new TCanvas((detectorName[i]).c_str(), (detectorName[i]).c_str(), 800, 800);
    if (GetDetector(i) <= 4)
    {
        H_Run[i]->GetXaxis()->SetRangeUser(45000, 48000);
        H_Run_Ref[i]->GetXaxis()->SetRangeUser(45000, 48000);
    }
    else
    {
        H_Run[i]->GetXaxis()->SetRangeUser(45000, 48000);
        H_Run_Ref[i]->GetXaxis()->SetRangeUser(45000, 48000);
    }
    H_Run[i]->Scale(H_Run_Ref[i]->Integral() / H_Run[i]->Integral());
    H_Run[i]->Draw("HIST");
    H_Run_Ref[i]->SetLineColor(kRed);
    H_Run_Ref[i]->Draw("SAME");
    c->Write();

    H_Run_Ref[i]->Write();
    H_Run[i]->Write();
}

#endif