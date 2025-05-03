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


TH1D* Channel[SIGNAL_MAX];
TGraphErrors *G_Resolution = new TGraphErrors();
double Mean = 0;
bool FLAG_PEAK_FINDER = false;
vector<int> DetectorExlude;

void BuildHistograms()
{
    Start("Build Histograms");
    // Initialize the histograms
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            Channel[i] = new TH1D(detectorName[i].c_str(), detectorName[i].c_str(), eSiliN, eSiliMin, eSiliMax);
            Channel[i]->GetXaxis()->SetTitle("Channel");
            Channel[i]->GetYaxis()->SetTitle("Counts");

        }
    }

    // Reader TreeGroup from FAST2ROOT file
    if (YEAR != 2021)
    {
        TTree *tree = (TTree *)file->Get("Tree_Group");
        TTreeReader *Reader = new TTreeReader(tree);
        TTreeReaderArray<Signal> *signals = new TTreeReaderArray<Signal>(*Reader, "Signal");

        int Entries = tree->GetEntries();
        while (Reader->Next())
        {
            for (int i = 0; i < signals->GetSize(); i++)
            {
                Signal signal = signals->At(i);
                if (IsDetectorSiliStrip(signal.Label))
                {
                    Channel[signal.Label]->Fill(signal.Channel);
                }
            }
        }
    }
    else
    {
        TTree *tree = (TTree *)file->Get("Tree");
        TTreeReader *Reader = new TTreeReader(tree);
        TTreeReaderValue<Signal> *signal = new TTreeReaderValue<Signal>(*Reader, "Signal");

        int Entries = tree->GetEntries();
        while (Reader->Next())
        {
            if (IsDetectorSiliStrip((**signal).Label))
            {
                Channel[(**signal).Label]->Fill((**signal).Channel);
            }
        }
    }
}

void FittingHistograms()
{
    Start("Fitting Histograms");
    ofstream out(("Config_Files/" + to_string(YEAR) + "/Electronic_Resolution_" + to_string(YEAR) + ".txt").c_str());

    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        
        if (IsDetectorSiliStrip(i))
        {
            if (FLAG_PEAK_FINDER) // using peak finder for 1 peak
            {
                TSpectrum *s = new TSpectrum();
                int nfound = s->Search(Channel[i], 2, "", 0.1);
                double *xpeaks = s->GetPositionX();
                Mean = xpeaks[0];
            }
            
            // fit gauss
            Channel[i]->GetXaxis()->SetRangeUser(Mean - 4000, Mean + 4000);
            TFitResultPtr r = Channel[i]->Fit("gaus", "S RQ", "", Mean - 4000, Mean + 4000);

            if (r != 0)
            {
                out << i << " " << 0 << " " << 0 << endl;
                Warning("Fit failed for " + detectorName[i]);
                continue;
            }

            out << i << " " << Channel[i]->GetFunction("gaus")->GetParameter(2) << " " << Channel[i]->GetFunction("gaus")->GetParError(2) << endl;

            G_Resolution->AddPoint(i, Channel[i]->GetFunction("gaus")->GetParameter(2));
            G_Resolution->SetPointError(G_Resolution->GetN()-1, 0, Channel[i]->GetFunction("gaus")->GetParError(2));

            Info(detectorName[i] + " : " + to_string(Channel[i]->GetFunction("gaus")->GetParameter(2)) + " CH", 1);
        }
    }

    out.close();
}


void WriteHistograms()
{
    Start("Write Histograms");
    FINAL_File->cd();

    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            Channel[i]->Write();
        }
    }
    G_Resolution->Write();

    FINAL_File->Close();
}   