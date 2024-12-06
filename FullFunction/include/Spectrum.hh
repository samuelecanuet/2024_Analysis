#ifndef __SPECTRUM_HH__
#define __SPECTRUM_HH__

#include "Peak.hh"
#include "Constants.hh"
#include "Addition.hh"

class Spectrum
{
public:
    Spectrum(int, int, string);
    ~Spectrum();

    string PATH = "../../CRADLE/RadiationData/";
    string Parent_filename;
    string Daughter_filename;

    int fA;
    int fZ;

    // map<int, Peak*> Peaks;
    vector<Peak *> Peaks;
    vector<TF1 *> fPeaks;
    TF1 *fSpectrum;

    void GenerateSpectrum();
    void ReadFile();
    TF1 *GetFunction();
    vector<TF1 *> GetTF1Peaks();
    vector<Peak *> GetPeaks();

    double SecondTokeV(double t)
    {
        if (t == 0)
        {
            return 0;
        }   
        return hbar / t * log(2) / keV;
    }

    double f(double *x, double *par)
    {
        return SpectrumAdd->Evaluate(x, par);
    }   

    Addition *SpectrumAdd;
};

#endif
