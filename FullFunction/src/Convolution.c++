#include "Convolution.hh"


Convolution::Convolution(MyTF1 *f1, TF1 *f2, double xmin, double xmax, int nsteps)
    : fFunction1(f1), fFunction2(f2), fXmin(xmin), fXmax(xmax), fNSteps(nsteps) 
{
}

Convolution::~Convolution()
{
}

double Convolution::Evaluate(double *x, double *params) 
    {
        bool changed = false;
        for (int i = 0; i < fFunction1->GetNpar(); ++i) 
        {
            if (params[i] != fFunction1->GetParameter(i)) 
            {
                changed = true;
                fFunction1->SetParameter(i, params[i]);
            }            
        }
        for (int i = 0; i < fFunction2->GetNpar(); ++i) 
        {
            fFunction2->SetParameter(i, params[i + fFunction1->GetNpar()]);
        }

        double result = 0.0;
        double step = (fXmax - fXmin) / fNSteps;
        for (int i = 0; i < fNSteps; ++i) 
        {
            double t = fXmin + i * step;
            double y1;               
            if (!changed)
                y1 = fFunction1->GetHandlerValue((int)t);
            else
                y1 = fFunction1->Eval(t);
            double y2 = fFunction2->Eval(x[0]-t);
            result +=  y1 * y2 * step;
        }
        return result;
    }

TF1* Convolution::GetFunction1()
{
    return fFunction1;
}

TF1* Convolution::GetFunction2()
{
    return fFunction2;
}
