#ifndef __MYTF1_HH__
#define __MYTF1_HH__

#include "TF1.h"
#include <functional>
#include <vector>
#include <string>

class MyTF1 : public TF1
{
public:
    MyTF1(const char *name, double (*fcn)(double *, double *),
          double xmin, double xmax, int npar = 0);

    MyTF1(const char *name, std::function<double(double*, double*)> func,
          double xmin, double xmax, int npar = 0);
   
    
    ~MyTF1();
    

    double KinematicShift_HANDLER_VALUE[10000] = {-1};

    virtual Double_t Eval(Double_t x)
    {
        Double_t result = TF1::Eval(x);
        KinematicShift_HANDLER_VALUE[(int)x] = result;
        return result;
    }

    double GetHandlerValue(int i)
    {
        if (KinematicShift_HANDLER_VALUE[i] == -1)
            KinematicShift_HANDLER_VALUE[i] = Eval(i);
        return KinematicShift_HANDLER_VALUE[i];
    }
};

#endif