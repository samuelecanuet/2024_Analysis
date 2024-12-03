#include "Peak.hh"

Peak::Peak(double Energy, double BranchingRatio, double HalfLife, double Q_beta, int A, int Z, double a, double Phi_max, double Phi_min, double E_min, double E_max)
{
    fEnergy = Energy;
    fBranchingRatio = BranchingRatio;
    fHalfLife = HalfLife;
    fQ_beta = Q_beta;
    fA = A;
    fZ = Z;
    fa = a;
    fPhi_max = Phi_max;
    fPhi_min = Phi_min;
    fE_min = E_min;
    fE_max = E_max;

    MyTF1 *KinematicShift = new MyTF1("KinematicShift", fTotalWeightKinematicShift, E_min, E_max, 9);
    KinematicShift->SetParameters(1, fEnergy, fQ_beta, K(fEnergy, fA * m), fa, fA, fZ, fPhi_min, fPhi_max);

    TF1 *NuclearBroadering = new TF1("NuclearBroadering", fNuclearBroadering, E_min, E_max, 2);
    NuclearBroadering->SetParameters(1, fHalfLife);

    Convolution ConvolutedNucKin(KinematicShift, NuclearBroadering, E_min, E_max);
    fFunction = new MyTF1("fFunction", [&ConvolutedNucKin](double *x, double *p) {
        return ConvolutedNucKin.Evaluate(x, p);
    }, E_min, E_max, 11);
    fFunction->SetParameters(1, fEnergy, fQ_beta, K(fEnergy, fA * m), fa, fA, fZ, fPhi_min, fPhi_max, 0, fHalfLife);
}

Peak::~Peak()
{
}

void Peak::SetParameter(int i, double value)
{
    fFunction->SetParameter(i, value);
}

MyTF1* Peak::GetFunction()
{
    return fFunction;
}

TF1* Peak::GetLifeTimeFunction()
{
    return fLifeTimeFunction;
}

TF1* Peak::GetKinematicFunction()
{
    return fKinematicFunction;
}

