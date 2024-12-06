#ifndef __ADDITION_HH__
#define __ADDITION_HH__


#include "TF1.h"
#include "Constants.hh"

class Addition
{
public:
    Addition(vector<TF1*>);
    ~Addition();

    vector<TF1*> fFunctions;
    vector<bool> fSaved; 

    int fNSteps;
        
    double Evaluate(double *x, double *params);
    int GetNPeak();
    int GetXmin();
    int GetXmax();
    int GetNpar();

    TF1* GetFunction(int);

};

#endif
