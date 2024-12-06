#include "Spectrum.hh"

Spectrum::Spectrum(int Z, int A, string decay)
{
    fA = A;
    fZ = Z;
    Parent_filename = PATH+"z"+to_string(fZ)+".a"+to_string(fA);
    if (decay == "proton")
    {
        Daughter_filename = PATH+"z"+to_string(fZ-1)+".a"+to_string(fA);
    }
    else
    {
        cout << "ERROR : ONLY PROTON DECAY IMPLEMENTED" << endl;
    }
}

Spectrum::~Spectrum()
{
}

TF1* Spectrum::GetFunction()
{
    return fSpectrum;
}

vector<TF1*> Spectrum::GetTF1Peaks()
{
    return fPeaks;
}

vector<Peak*> Spectrum::GetPeaks()
{
    return Peaks;
}

void Spectrum::GenerateSpectrum()
{
    ReadFile();
    // vector<TF1* > fPeakss = {fPeaks[0], fPeaks[1], fPeaks[2]};
    SpectrumAdd = new Addition(fPeaks);
    fSpectrum = new TF1("fSpectrum", [this](double *x, double *p) {
        return this->f(x, p);
    }, SpectrumAdd->GetXmin(), SpectrumAdd->GetXmax(), SpectrumAdd->GetNpar());
}


void Spectrum::ReadFile()
{
    ifstream file(Daughter_filename);
    string line;

    while (getline(file, line))
    {
        istringstream iss(line);

        string p;
        string flag;
        double LevelEnergy;
        double LevelHalfLife;
        string proton;
        double Daughter_energy;
        double Intensity;
        double Qp;

        if (line[0] == '#')
        {
            continue;
        }
        else if (!line.compare(0, 1, "P"))
        {
            iss >> p >> LevelEnergy >> flag >> LevelHalfLife;
            continue;
        }
        else
        {
            iss >> proton >> Daughter_energy >> flag >> Intensity >> Qp;

            double Energy = Qp*(fA-1)/fA;

            ifstream file1(Parent_filename);
            string line1;
            while (getline(file1, line1))
            {
                istringstream iss1(line1);
                string p1;
                double LevelEnergy1;
                double LevelHalfLife1;
                string proton1;
                double Daughter_energy1;
                string flag1;
                double Intensity1;
                double Q_beta;

                if (line1[0] == '#')
                {
                    continue;
                }
                else if (!line1.compare(0, 1, "P"))
                {
                    continue;
                }
                else
                {
                    iss1 >> proton1 >> Daughter_energy1 >> flag1 >> Intensity1 >> Q_beta;
                    if (abs(Daughter_energy1 - LevelEnergy) < 1)
                    {

                        if (Energy > 0)
                        {   
                            double a;
                            if (abs(Q_beta - 6088) < 1)
                            {
                                a = 1;
                            }
                            else
                            {
                                a = -1./3.;
                            }
                            Peak* p = new Peak(Energy, Intensity, SecondTokeV(LevelHalfLife), sqrt((Q_beta-511) * (Q_beta-511) + 511 * 511), fA, fZ, a, M_PI, 0, Energy-250, Energy+250, 0.1);
                            Peaks.push_back(p);
                            fPeaks.push_back(p->GetFunction());
                        }
                        break;
                    }
                }
            }
        }
    }
}
