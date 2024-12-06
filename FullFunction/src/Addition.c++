#include "Addition.hh"

Addition::Addition(vector<TF1 *> functions) : fFunctions(functions)
{
    fSaved.resize(fFunctions.size(), false);

    cout << &fFunctions << endl;
    
}

Addition::~Addition()
{
}

double Addition::Evaluate(double *x, double *params)
{
    double res = 0;
    for (int i = 0; i < fFunctions.size(); ++i)
    {
        TF1 *Function = fFunctions[i];
        bool changed = false;

        for (int ipar = 0; ipar < Function->GetNpar(); ipar++)
        {
            if (Function->GetParameter(ipar) != params[i * Function->GetNpar() + ipar])
            {
                changed = true;
                fSaved[i] = false;
                Function->SetParameter(ipar, params[i * Function->GetNpar() + ipar]);
            }
        }

        if (changed || !fSaved[i])
        {
            if (x[0] > Function->GetXmin() && x[0] < Function->GetXmax())
            {
                Function->Save(Function->GetXmin(), Function->GetXmax(), 0, 0, 0, 0);
                fSaved[i] = true;
            }
            else
            {
                continue;   
            }
        }

        res += Function->GetSave(x);
    }
    return res;
}

TF1 *Addition::GetFunction(int i)
{
    return fFunctions[i];
}

int Addition::GetNPeak()
{
    return fFunctions.size();
}

int Addition::GetXmin()
{
    double xmin = fFunctions[0]->GetXmin();
    for(int i = 0; i < fFunctions.size(); ++i)
    {
        if(fFunctions[i]->GetXmin() < xmin)
        {
            xmin = fFunctions[i]->GetXmin();
        }
    }
    return xmin;
}

int Addition::GetXmax()
{
    double xmax = fFunctions[0]->GetXmax();
    for(int i = 0; i < fFunctions.size(); ++i)
    {
        if(fFunctions[i]->GetXmax() > xmax)
        {
            xmax = fFunctions[i]->GetXmax();
        }
    }
    return xmax;
}

int Addition::GetNpar()
{
    int npar = 0;
    for(int i = 0; i < fFunctions.size(); ++i)
    {
        npar += fFunctions[i]->GetNpar();
    }
    return npar;
}