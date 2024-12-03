#ifndef __CONVOLUTION_HH__
#define __CONVOLUTION_HH__


#include "TF1.h"
#include "MyTF1.hh"


class Convolution
{
public:
    Convolution(MyTF1 *f1, TF1 *f2, double xmin, double xmax, int nsteps = 1000);
    ~Convolution();

    MyTF1* fFunction1;
    TF1* fFunction2;

    double fXmin;
    double fXmax;
    int fNSteps;
        
    double Evaluate(double *x, double *params);

    TF1* GetFunction1();
    TF1* GetFunction2();

};

#endif
