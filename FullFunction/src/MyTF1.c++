#include "MyTF1.hh"

MyTF1::MyTF1(const char *name, double (*fcn)(double *, double *),
             double xmin, double xmax, int npar)
    : TF1(name, fcn, xmin, xmax, npar)
{
}

MyTF1::MyTF1(const char *name, std::function<double(double *, double *)> func,
             double xmin, double xmax, int npar)
    : TF1(name, func, xmin, xmax, npar)
{
}

MyTF1::~MyTF1()
{
}