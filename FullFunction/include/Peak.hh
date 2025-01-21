#ifndef __PEAK_HH__
#define __PEAK_HH__

#include "TF1.h"
#include "TCanvas.h"
#include "Convolution.hh"
#include "Constants.hh"

class Peak
{
public:
    Peak(double Energy, double BranchingRatio, double HalfLife, double Q_beta, int A, int Z, double a, double Phi_max, double Phi_min, double E_min, double E_max, double step);
    ~Peak();

    TF1 *GetFunction();
    TF1 *GetLifeTimeFunction();
    TF1 *GetKinematicFunction();

    // TF1 *GetBetaSpectrum();

    Convolution *ConvolutedNucKin;
    double fE_min_func;
    double fE_max_func;

    double GetK()
    {
        return K(fEnergy, fA * m);
    }

    vector<double> Parameters;
    bool FinalParameters = false;
    void FixAllParameters(TF1 *, int);

    void ExternalFixAllParameters(TF1 *, int, vector<int>);
    void ExternalSetParameter(TF1 *, int, int);
    void ExternalSetParameter(TF1 *, int, int, double);
    void ExternalSetParLimits(TF1 *, int, int, double, double);

    double GetBR();

    int GetNPeak();

    // fucntion to display in cout a peak : cout << Peak << endl;
    friend ostream &operator<<(ostream &os, const Peak &p)
    {
        os << "Peak : " << p.fEnergy << " MeV" << endl;
        os << "Branching Ratio : " << p.fBranchingRatio << endl;
        os << "Half Life : " << p.fHalfLife << " s" << endl;
        os << "Q_beta : " << p.fQ_beta << " MeV" << endl;
        os << "A : " << p.fA << endl;
        os << "Z : " << p.fZ << endl;
        os << "a : " << p.fa << endl;
        os << "Phi_max : " << p.fPhi_max << endl;
        os << "Phi_min : " << p.fPhi_min << endl;
        os << "E_min : " << p.fE_min << endl;
        os << "E_max : " << p.fE_max << endl;
        os << "Step : " << p.fStep << endl;
        return os;
    }

private:
    double fEnergy;
    double fBranchingRatio;
    double fHalfLife;
    double fQ_beta;
    int fA;
    int fZ;
    double fa;
    double fPhi_max;
    double fPhi_min;
    double fE_min;
    double fE_max;

    TF1 *fFunction;
    TF1 *TF1NuclearBroadering;
    TF1 *TF1KinematicFunction;

    int NParKinematic;
    int NParNuclear;
    int NPar;
    double fStep;

    double f(double *x, double *par)
    {
        return ConvolutedNucKin->Evaluate(x, par);
    }

    /// COMMON FUNCTION ///
    double EnergyToMomentum(double E, double Mass)
    {
        return sqrt(E * E - Mass * Mass);
    }

    double MomentumToEnergy(double p, double Mass)
    {
        return sqrt(p * p + Mass * Mass);
    }

    double EnergyToVelocity(double E, double Mass)
    {
        return E / Mass;
    }

    /// KINEMATICS ///

    double K(double Ep, double M)
    {
        return sqrt(2 * Ep / M * (m_p * (M - m_p)) / (M * M));
    }

    double Fermi(int Z, double E_b, double M)
    {
        double v = EnergyToVelocity(E_b, M);
        double nu = Z * alpha / v * c;
        return 2 * M_PI * nu / (1 - exp(-2 * M_PI * nu));
    }

    double fWmin(double W0, double k, double t)
    {
        if (t / k < -(W0 - m_e))
            return (m_e * m_e + (W0 + t / k) * (W0 + t / k)) / (2 * (W0 + t / k));
        else if (abs(t / k) <= W0 - m_e)
            return m_e;
        else
            return (m_e * m_e + (W0 - t / k) * (W0 - t / k)) / (2 * (W0 - t / k));
    }

    double WeightKinematicShift(double t, double k, double W0, double W, double a, double M, double Z, double Phi_min, double Phi_max)
    {
        double p_b = EnergyToMomentum(W, m_e);
        if (p_b == 0)
            return 0;

        double c1 = min(cos(Phi_max), max(cos(Phi_min), (W - W0 - t / k) / p_b));
        double c2 = max(cos(Phi_min), min(cos(Phi_max), (W0 - W - t / k) / p_b));

        double I1 = (c2 - c1) * W * (W0 - W);
        double I2a = -a * (pow(c2, 2) - pow(c1, 2)) / 2.0 * p_b * t / k;
        double I2b = -a * (pow(c2, 3) - pow(c1, 3)) / 3.0 * pow(p_b, 2);

        return Fermi(Z, W, M) * p_b / (2 * k) * (I1 + I2a + I2b) * 1e-28;
    }

    double fTotalWeightKinematicShift(double *t, double *par)
    {
        double A = par[0];
        double E_p = par[1];
        double W0 = par[2];
        double k = par[3];
        double a = par[4];
        double A = par[5];
        double Z = par[6];
        double Phi_min = par[7];
        double Phi_max = par[8];

        double t_ = t[0] - E_p;
        double M = A * m;

        if (abs(t_ / k) > sqrt(W0 * W0 - m_e * m_e))
        {
            return 0;
        }

        double wt = 0;
        double Wmin = fWmin(W0, k, t_);
        double step = (W0 - Wmin) / 1000.0;

        for (double Wb = Wmin; Wb <= W0; Wb += step)
        {
            wt += A*WeightKinematicShift(t_, k, W0, Wb, a, M, Z, Phi_min, Phi_max) * step;
        }

        return wt;
    }

    /// NUCLEAR BROADENING ///
    double fNuclearBroadering(double *x, double *par)
    {
        double E = x[0];
        double A = 1;
        double Gamma = par[1];

        // condition if there is no half-life
        if (Gamma == 0)
        {
            if (E == 0)
            {
                return A;
            }
            else
            {
                return 0;
            }
        }

        return A * Gamma / (2 * M_PI) / (E * E + Gamma * Gamma / 4);
    }
};

#endif