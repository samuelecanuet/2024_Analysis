#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "myRooFit.h"
using namespace RooFit;


double ErrorDef(int nPar, double sigma)
{
    // run a terminal command: python3 ConfidenceLevel.py -n 1000000 -p 3.526740380330011 -e 0.001 -s 0.01
    string command = "python3 ../ConfidenceLevel.py --dir ../ --sigma " + to_string(sigma) + " --npar " + to_string(nPar);
    system(command.c_str());

    // get the print() python output written in the terminal
    double result = 0.0;
    ifstream infile("ConfidenceLevel_Dim.txt");
    if (infile.is_open()) {
        infile >> result;
        infile.close();
    } else {
        std::cerr << "Error opening file." << std::endl;
    }

    return pow(result, 2);
}

double GetNFittedPar(RooAbsReal *nll, RooDataSet *data = nullptr)
{
    // get number of non fixed paramaters
    int nPar = 0;
    RooArgSet *params = nll->getParameters(data);
    RooFIter iter = params->fwdIterator();
    RooAbsArg *arg;
    while ((arg = iter.next())) {
        if (!arg->isConstant()) {
            nPar++;
        }
    }
    return nPar;
}

void myRooFit()
{
    // Contruct obs
    RooRealVar E("E", "E", 3250, 3450);

    // contruct BreitWigner
    RooRealVar mu_l("mu_l", "mean BreitWigner", 0);
    RooRealVar sigma_l("sigma_l", "sigma BreitWigner", 1);
    RooBreitWigner BreitWigner("BreitWigner", "BreitWigner", E, mu_l, sigma_l);

    // contrsuct gauss
    RooRealVar mu_g("mu_g", "mean Gauss", 0);
    RooRealVar sigma_g("sigma_g", "sigma Gauss", 10, 0, 50);
    RooGaussian gauss("gauss", "gauss", E, mu_g, sigma_g);

    E.setBins(1000, "cache");

    RooFFTConvPdf convPdf("convPdf", "convolution of Gauss and BreitWigner", E, BreitWigner, gauss);

    RooDataSet *data = convPdf.generate(E, 1000000);

    // --- Fit the data with the convolution PDF ---
    RooRealVar N("PeakN", "Peak number", 1);
    RooRealVar E_mean("E_mean", "Energy mean", 3352);
    RooRealVar HalfLife("HalfLife", "Half-life", 20e-3);
    RooRealVar Qb("Qb", "Qbeta", sqrt(pow(5046 + 1022 - 511, 2) + pow(511, 2)));
    RooRealVar A("A", "Mass number", 32);
    RooRealVar Z("Z", "Atomic number", 18);
    RooRealVar a("a", "Fermi parameter", 1);
    RooRealVar Phi_min("Phi_min", "Minimum angle", 0);
    RooRealVar Phi_max("Phi_max", "Maximum angle", TMath::Pi());

    PeakShape peakShape("peak", "Peak",
        E, N, E_mean, Qb, A, Z, a, Phi_min, Phi_max);

    Peak peak("peak", "Peak",
        E, N, E_mean, HalfLife, Qb, A, Z, a, Phi_min, Phi_max);

    RooFFTConvPdf conv = RooFFTConvPdf("conv", "Convolution", E, peakShape, );

    new TCanvas ("myRooFit", "myRooFit", 600, 600);
    gPad->SetLeftMargin(0.15);
    RooPlot *frame = E.frame(Title("Convolution of Gauss and BreitWigner"));
    conv.plotOn(frame);
    // data->plotOn(frame);
    frame->Draw();
    
    




}