#ifndef __PEAK_HH__
#define __PEAK_HH__

#include "TF1.h"
#include "Convolution.hh"
#include "Constants.hh"
#include "MyTF1.hh"

class Peak
{
    public:
        Peak(double, double, double, double, int , int , double , double, double, double, double);
        ~Peak();

        MyTF1 *GetFunction();
        TF1 *GetLifeTimeFunction();
        TF1 *GetKinematicFunction(); 

        double GetHandlerValue(int i); 
        void SetParameter(int i, double value);

        // TF1 *GetBetaSpectrum();


        // void FixAllParameters();     
    
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

        MyTF1 *fFunction;
        TF1 *fLifeTimeFunction;
        MyTF1 *fKinematicFunction;

        // double EnergyToMomentum(double E, double Mass);
        // double MomentumToEnergy(double p, double Mass);
        // double EnergyToVelocity(double E, double Mass);
        // double K(double E, double M);
        // double Fermi(int Z, double E_b, double M);
        // double WeightKinematicShift(double t, double k, double W0, double W, double a, double M, double Z, double Phi_min, double Phi_max);
        // // double fTotalWeightKinematicShift(double *t, double *par);
        // double fWmin(double W0, double k, double t);
        // // double fNuclearBroadering(double *x, double *par);
};

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
    double I2a = - a * (pow(c2, 2) - pow(c1, 2)) / 2.0 * p_b * t / k;
    double I2b = - a * (pow(c2, 3) - pow(c1, 3)) / 3.0 * pow(p_b, 2);

    return Fermi(Z, W, M) * p_b / (2 * k) * (I1 + I2a + I2b) * 1e-28;
}

double fTotalWeightKinematicShift(double *t, double *par)
{   
    double Amp = par[0];
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
        return 0;

    double wt = 0;
    double Wmin = fWmin(W0, k, t_);
    double step = (W0 - Wmin) / 1000.0;

    for (double Wb = Wmin; Wb <= W0; Wb += step)
    {
        wt += Amp*WeightKinematicShift(t_, k, W0, Wb, a, M, Z, Phi_min, Phi_max) * step;
    }

    // KinematicShift_HANDLER_VALUE[(int)t[0]] = wt;
    
    return wt;
}

/// NUCLEAR BROADENING ///
double fNuclearBroadering(double *x, double *par)
{
    double E = x[0];
    double A = par[0];
    double Gamma = par[1];

    return A * Gamma / (2 * M_PI) / (E * E + Gamma * Gamma / 4);
}

#endif