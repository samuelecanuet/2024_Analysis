#include "Convolution.hh"

Convolution::Convolution(TF1 *f1, TF1 *f2, double xmin, double xmax, double step)
    : fFunction1(f1), fFunction2(f2), fXmin(xmin), fXmax(xmax), fStep(step)
{
    fNSteps = (int)((fXmax - fXmin) / fStep);
    fFunction1->SetNpx(fNSteps);
    fFunction2->SetNpx(fNSteps);
}

Convolution::~Convolution()
{
}

double Convolution::Evaluate(double *x, double *params)
{
    bool f1Changed = false;
    for (int i = 1; i < fFunction1->GetNpar(); ++i)
    {
        if (fFunction1->GetParameter(i) != params[i])
        {
            f1Changed = true;
            f1Saved = false;
            fFunction1->SetParameter(i, params[i]);
        }
    }
    bool f2Changed = false;
    for (int i = 0; i < fFunction2->GetNpar(); ++i)
    {
        if (fFunction2->GetParameter(i) != params[fFunction1->GetNpar() + i])
        {
            f2Changed = true;
            f2Saved = false;
            fFunction2->SetParameter(i, params[fFunction1->GetNpar() + i]);
        }
    }

    if (f1Changed || !f1Saved)
    {
        fFunction1->Save(fXmin, fXmax, 0, 0, 0, 0);
        f1Saved = true;
    }

    if (f2Changed || !f2Saved)
    {
        fFunction2->Save(fXmin-fXmax, fXmax-fXmin, 0, 0, 0, 0);
        f2Saved = true;
    }

    double result = 0.0;
    for (int i = 0; i < fNSteps; ++i)
    {
        double t = fXmin + i * fStep;
        double y1 = fFunction1->GetSave(&t);
        double xx = x[0] - t;
        double y2 = fFunction2->GetSave(&xx);
        result += y1 * y2 * fStep;
    }
        
    return result;
}

TF1 *Convolution::GetFunction1()
{
    return fFunction1;
}

TF1 *Convolution::GetFunction2()
{
    return fFunction2;
}

