#include "Detectors.hh"


TFile *SIMULATED_File;
TFile *MERGED_File;

TH1D* H_Time_Sim;
TH1D* H_Time_Exp;

TGraph* G_N = new TGraph();
TH1D* H_SiPM_LeftLeft;
TH1D* H_SiPM_Left;
TH1D* H_SiPM_Right;
TH1D* H_SiPM_Center;
TH1D* H_S;
TH1D* H_Time_MC = new TH1D("H_Time_MC", "H_Time_MC", 300, -300, 300);
TH1D* H_Fake = new TH1D("H_Fake", "H_Fake", 300, -300, 300);

TH1D* H_Ar = new TH1D("H_Ar", "H_Ar", 300, -300, 300);

default_random_engine generator;

double start_gate;
double end_gate;

double time_resolution;
////////////////////////////

int counter_real;
int counter_fake;

TF1* fScintillator = new TF1("fScintillator", "exp(-x*0.6/2)", 0, 300);
TH1D* HScintillator;

void Fill(double time, bool Ar, bool Beta)
{
    normal_distribution<double> distribution(0, time_resolution);
    time += distribution(generator);
    time += H_Time_Sim->GetRandom();
    if (time > 240)
    {
        return;
    }
    H_Time_MC->Fill(time);


    if ((Beta && Ar) || (Beta && !Ar))
    {
        H_Fake->Fill(time);
        if (time > start_gate && time < end_gate)
            counter_fake++;
    }
    else
    {
        if (time > start_gate && time < end_gate)
            counter_real++;
    }

    if (Ar && !Beta)
    {
        H_Ar->Fill(time);
    }
}
