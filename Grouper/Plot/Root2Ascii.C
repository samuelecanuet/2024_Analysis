#include <iostream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <ostream>
using namespace std;

int Root2Ascii()
{
    TFile *file = new TFile("/mnt/hgfs/shared-2/2025_DATA/DETECTOR_DATA/CALIBRATED/Calibrated_2025.root", "READ");
    TCanvas *c = (TCanvas *)file->Get("32Ar_SUM_Down");
    TH1D *h = (TH1D *)c->GetPrimitive("32Ar_Exp_All_Down");

    h->Rebin(10);

    // in ascii
    ofstream outfile("32Ar_Exp_All_Down.txt");
    for (int i = 1; i <= h->GetNbinsX(); i++)
    {
        outfile << h->GetBinCenter(i) << "\t" << h->GetBinContent(i) << endl;
    }
    outfile.close();
    return 0;
}