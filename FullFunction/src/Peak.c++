#include "Peak.hh"

Peak::Peak(double Energy, double BranchingRatio, double HalfLife, double Q_beta, int A, int Z, double a, double Phi_max, double Phi_min, double E_min, double E_max, double step)
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
    fStep = step;

    NParKinematic = 9;
    NParNuclear = 2;
    NPar = NParKinematic + NParNuclear;

    Parameters = {1, fEnergy, fQ_beta, K(fEnergy, fA * m), fa, (double)fA, (double)fZ, fPhi_min, fPhi_max, 1, fHalfLife};

    fE_min_func = -K(fEnergy, fA * m) * sqrt(fQ_beta * fQ_beta - m_e * m_e) + fEnergy;
    fE_max_func = K(fEnergy, fA * m) * sqrt(fQ_beta * fQ_beta - m_e * m_e) + fEnergy;

    TF1KinematicFunction = new TF1("KinematicShift", [this](double *x, double *par)
                                   { return this->fTotalWeightKinematicShift(x, par); }, fE_min, fE_max, NParKinematic);
    TF1KinematicFunction->SetParameter(0, 1);
    TF1KinematicFunction->FixParameter(1, fEnergy);
    TF1KinematicFunction->FixParameter(2, fQ_beta);
    TF1KinematicFunction->FixParameter(3, K(fEnergy, fA * m));
    TF1KinematicFunction->FixParameter(4, fa);
    TF1KinematicFunction->FixParameter(5, fA);
    TF1KinematicFunction->FixParameter(6, fZ);
    TF1KinematicFunction->FixParameter(7, fPhi_min);
    TF1KinematicFunction->FixParameter(8, fPhi_max);

    TF1NuclearBroadering = new TF1("NuclearBroadering", [this](double *x, double *par)
                                   { return this->fNuclearBroadering(x, par); }, -300, 300, NParNuclear);
    TF1NuclearBroadering->SetParameters(1, fHalfLife);

    ConvolutedNucKin = new Convolution(TF1KinematicFunction, TF1NuclearBroadering, fE_min, fE_max, fStep);
    fFunction = new TF1("fFunction", [this](double *x, double *par)
                        { return par[9]*this->f(x, par); }, fE_min, fE_max, NPar);
    fFunction->SetParameter(0, 1);
    fFunction->FixParameter(1, fEnergy);
    fFunction->FixParameter(2, fQ_beta);
    fFunction->FixParameter(3, K(fEnergy, fA * m));
    fFunction->FixParameter(4, fa);
    fFunction->FixParameter(5, fA);
    fFunction->FixParameter(6, fZ);
    fFunction->FixParameter(7, fPhi_min);
    fFunction->FixParameter(8, fPhi_max);
    fFunction->SetParameter(9, 1e5);
    fFunction->SetParLimits(9, 1e2, 1e7);
    fFunction->FixParameter(10, HalfLife);
    // fFunction->SetParLimits(10, 0, 10);
}

Peak::~Peak()
{
}

TF1 *Peak::GetFunction()
{
    return fFunction;
}

TF1 *Peak::GetLifeTimeFunction()
{
    return TF1NuclearBroadering;
}

TF1 *Peak::GetKinematicFunction()
{
    return TF1KinematicFunction;
}

void Peak::FixAllParameters(TF1 *f, int peak_number)
{
    FinalParameters = true;
    for (int ipar = 0; ipar < NPar; ipar++)
    {
        cout << "Fixing parameter " << ipar << " to " << f->GetParameter(peak_number * NPar + ipar) << endl;
        fFunction->FixParameter(ipar, f->GetParameter(peak_number * NPar + ipar));
    }

    for (int ipar = 0; ipar < NParKinematic; ipar++)
    {
        TF1KinematicFunction->FixParameter(ipar, f->GetParameter(peak_number * NPar + ipar));
    }

    for (int ipar = 0; ipar < NParNuclear; ipar++)
    {
        TF1NuclearBroadering->FixParameter(ipar, f->GetParameter(peak_number * NPar + NParKinematic + ipar));
    }
}

void Peak::ExternalFixAllParameters(TF1 *f, int peak_number, vector<int> except)
{
    for (int ipar = 0; ipar < NPar; ipar++)
    {
        if (find(except.begin(), except.end(), ipar) == except.end())
        {
            f->FixParameter(NPar * peak_number + ipar, Parameters[ipar]);
        }
    }
}

void Peak::ExternalSetParameter(TF1 *f, int peak_number, int ipar)
{
    f->SetParameter(NPar * peak_number + ipar, Parameters[ipar]);
}

void Peak::ExternalSetParameter(TF1 *f, int peak_number, int ipar, double value)
{
    f->SetParameter(NPar * peak_number + ipar, value);
}

void Peak::ExternalSetParLimits(TF1 *f, int peak_number, int ipar, double low, double high)
{
    f->SetParLimits(NPar * peak_number + ipar, low, high);
}

int Peak::GetNPeak()
{
    return 1;
}

double Peak::GetBR()
{
    return -1;
}