#ifndef __CONVOLUTION_HH__
#define __CONVOLUTION_HH__


#include "TF1.h"
#include "Constants.hh"



class Convolution
{
public:
    Convolution(TF1 *f1, TF1 *f2, double xmin, double xmax, double step = 0.1);
    ~Convolution();

    TF1* fFunction1;
    TF1* fFunction2;

    bool f1Saved = false;
    bool f2Saved = false;

    double fXmin;
    double fXmax;
    int fNSteps;
    double fStep;
        
    double Evaluate(double *x, double *params);

    TF1* GetFunction1();
    TF1* GetFunction2();


};

#endif
