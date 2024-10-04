#include "Detectors.hh"
#include "../../../lib/SignalDict/Signal.h"

TH1D* H_SiPM[SIGNAL_MAX];
TF1* f_matching[SIGNAL_MAX];
TF1* f_fitted[SIGNAL_MAX];   
TGraph *g_calib[SIGNAL_MAX];
int det;

double guess_mean[10] = {250e3, 0, 0, 0, 0, 0, 0, 0, 0, 0};


void FitHistogram(TH1D *hist)
{
    f_fitted[det] = new TF1("fit", "gaus(0)+gaus(3)", eHighMin, eHighMax);
    f_fitted[det]->SetParameter(0, 5000);
    f_fitted[det]->SetParameter(1, 250e3);
    f_fitted[det]->SetParLimits(1, 200e3, 400e3);
    f_fitted[det]->SetParameter(2, 150e3);
    f_fitted[det]->SetParLimits(2, 50e3, 250e3);
    f_fitted[det]->SetParameter(3, 5000);
    f_fitted[det]->SetParameter(4, 600e3);
    f_fitted[det]->SetParLimits(4, 400e3, 900e3);
    f_fitted[det]->SetParameter(5, 100e3);
    f_fitted[det]->SetParLimits(5, 50e3, 250e3);
    hist->Fit(f_fitted[det], "QR", "", 200e3, 1200e3);
}