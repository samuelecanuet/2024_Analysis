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


default_random_engine generator;

double start_gate;
double end_gate;
////////////////////////////

int counter_real;
int counter_fake;

void Fill(double time, bool Ar, bool Beta)
{
    H_Time_MC->Fill(time);

    if (Beta && !Ar)
    {
        H_Fake->Fill(time);
        counter_fake++;
    }
    else
    {
        counter_real++;
    }
}

double FunctionToMinimize(const double *par)
{

  double guess_N = par[0];
  int mul = par[1];
  double mean_height_left_left = par[2];
  double mean_height_left = par[3];

  double alpha = (mean_height_left-mean_height_left_left);

    // compute f(t)
    TH1D* H_Fortuitous = (TH1D*)H_SiPM_Center->Clone("H_Fortuitous_");
    H_Fortuitous->Reset();
    TH1D* H_fx = (TH1D*)H_SiPM_Center->Clone("H_fx_");
    TH1D* H_Sx = (TH1D*)H_SiPM_Center->Clone("H_Sx_");
    H_fx->Reset();
    H_Sx->Reset();
    // cout << "alpha = " << alpha << endl;
    double Fx = 0;
    for (int bin = 0; bin < H_SiPM_Center->GetNbinsX(); bin++)
    {
        if (H_S->GetBinCenter(bin) < -20 || H_S->GetBinCenter(bin) > 70)
        {
            H_Fortuitous->SetBinContent(bin, mean_height_left_left);
        }
        else
        {
            double fx = (H_S->GetBinContent(bin) - alpha * (1 - Fx / guess_N)) / (1 - alpha / guess_N);
            H_Fortuitous->SetBinContent(bin, alpha * (1 - Fx / guess_N) + mean_height_left_left);
            H_fx->SetBinContent(bin, fx);

            Fx += fx;
        }
    }

    cout << "Guess N = " << guess_N << "     Integral = " << H_fx->Integral() << "    Diff = "  << abs(guess_N-H_fx->Integral()) << endl;
    


    H_Sx->Add(H_fx, 1);
    // cout << "Integral = " << H_fx->Integral() << endl;
    H_Sx->Add(H_Fortuitous, 1);
    // cout << "Integral = " << H_Sx->Integral() << endl;

    G_N->AddPoint(guess_N, guess_N-H_fx->Integral());

    


    return abs(guess_N-H_fx->Integral());
}
