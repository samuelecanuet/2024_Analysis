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

#include "../../../lib/SignalDict/Signal.h"
#include "Detectors.hh"

using namespace std;
using namespace ROOT::Math;
default_random_engine generator;


TFile *file;
TFile *FINAL_File;

TH1D* H_Channel[SIGNAL_MAX];
TGraphErrors *G_Linearity[SIGNAL_MAX];
TGraphErrors *G_Global =  new TGraphErrors();


void InitTGraphs()
{
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (!IsDetectorSiliStrip(i))
            continue;
        G_Linearity[i] = new TGraphErrors();
    }
}

void InitHistogramsForRun(int run)
{
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (!IsDetectorSiliStrip(i))
            continue;
        H_Channel[i] = new TH1D(Form("H_Channel_run%d_Det%d", run, i), Form("H_Channel_run%d_Det%d", run, i), 20000, 0, 100);
        H_Channel[i]->GetXaxis()->SetTitle("Channel");
        H_Channel[i]->GetYaxis()->SetTitle("Counts");
    }
}