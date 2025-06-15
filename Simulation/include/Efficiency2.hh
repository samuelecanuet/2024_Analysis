#include "../../Grouper/include/Detectors.hh"

int strip_ref = 1;
int detector_ref = 1;

double Chi2(const TGraphErrors* g1, const TGraphErrors* g2) {
    if (!g1 || !g2) {
        Error("One of the graphs is null.");
        return -1;
    }

    int n1 = g1->GetN();
    int n2 = g2->GetN();

    if (n1 != n2) {
        Error("Graphs have different number of points.");
        return -1;
    }

    double chi2 = 0.0;
    for (int i = 0; i < n1; ++i) {
        double x1, y1, x2, y2;
        g1->GetPoint(i, x1, y1);
        g2->GetPoint(i, x2, y2);

        if (std::abs(x1 - x2) > 1e-6) {
            Error("x-values at point " + std::to_string(i) + " differ: " + std::to_string(x1) + " vs " + std::to_string(x2));   
            continue;
        }

        if (x1 == strip_ref)
            continue;

        double err1 = g1->GetErrorY(i);
        double err2 = g2->GetErrorY(i);
        double err_tot2 = err1 * err1 + err2 * err2;

        if (err_tot2 <= 0) {
            Error("Total error squared is zero or negative at point " + std::to_string(i));
            continue;
        }

        double delta = y1 - y2;
        double contrib = (delta * delta) / err_tot2;
        chi2 += contrib;
    }

    return chi2;
}


struct Position
{
    double x;
    double y;
    double z;
    double theta;

    bool operator<(const Position& other) const {
        if (x*x + y*y + z*z <= other.x*other.x + other.y*other.y + other.z*other.z)
            return true;
        else
            return false;
    }

    bool operator==(const Position& other) const {
        return (x == other.x && y == other.y && z == other.z && theta == other.theta);
    }

    bool operator!=(const Position& other) const {
        return !(*this == other);
    }

    friend ostream& operator<<(ostream& os, const Position& pos) {
        os << "x: " << pos.x << ", y: " << pos.y << ", z: " << pos.z << ", theta: " << pos.theta;
        return os;
    }
};

bool Accepting_Position(Position pos)
{
    // if (pos.y > 0.0)
    //     return false;
    // if (abs(pos.y - pos.z) > 1.)
    //     return false;
    // if (pos.y != 0.5 || pos.z != 0.0)
    //     return false;
    // if (abs(pos.theta) > 1.0)
    //     return false;
    // if (pos.z != 0.0)
    //     return false;
    // if (pos.y == 0.0)
    //     return false;
    return true;
}


#include <TFile.h>
#include <TTree.h>
#include <TF3.h>
#include <Fit/Fitter.h>
#include <Math/Functor.h>
#include <iostream>
#include <vector>

double model(double *x, double *p) {
    double y = x[0];
    double z = x[1];
    double theta = x[2];

    return p[0] + p[1]*y + p[2]*z + p[3]*theta ;
    // +
    //        p[4]*y*z + p[5]*y*theta + p[6]*z*theta +
    //        p[7]*y*y + p[8]*z*z + p[9]*theta*theta;
}

void FitChi2Surface(TTree* tree) {
    double y, z, theta, chi2;
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("z", &z);
    tree->SetBranchAddress("theta", &theta);
    tree->SetBranchAddress("chi2", &chi2);

    std::vector<std::array<double, 3>> coords;
    std::vector<double> values;

    Long64_t N = tree->GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);
        
         // Skip out-of-bounds values
        if (y < -2.0 || y > 0.0 || z < 0.0 || z > 2.0 || theta < -1.0 || theta > 1.0) {
            continue;
        }

        coords.push_back({y, z, theta});
        values.push_back(chi2);
    }

    
    const int npar = 4;
    ROOT::Fit::Fitter fitter;

    // Allocate space for parameters first
    double par[npar] = {0.0};
    fitter.Config().SetParamsSettings(npar, par);

    for (int i = 0; i < npar; ++i)
        fitter.Config().ParSettings(i).SetValue(0.0);  // initial guesses

    auto chi2_func = [&](const double *par) {
        double sum = 0.0;
        for (size_t i = 0; i < coords.size(); ++i) {
            double x[3] = {coords[i][0], coords[i][1], coords[i][2]};
            double res = values[i] - model(x, (double*)par);
            sum += res * res;
        }
        return sum;
    };

    ROOT::Math::Functor f(chi2_func, npar);
    fitter.SetFCN(f);

    Start("Fitting Chi2 Surface");
    bool ok = fitter.FitFCN();
    if (!ok) {
        std::cerr << "Fit failed!" << std::endl;
        return;
    }

    const ROOT::Fit::FitResult &result = fitter.Result();
    result.Print(std::cout);

    TF3 *fit_func = new TF3("fit_func", "[0] + [1]*x + [2]*y + [3]*z"
    //  + "
    //                           "[4]*x*y + [5]*x*z + [6]*y*z + "
    //                           "[7]*x*x + [8]*y*y + [9]*z*z"
                              ,-2.0, 0.0, 0.0, 2.0, -1.0, 1.0);
    for (int i = 0; i < npar; ++i) {
        fit_func->SetParameter(i, result.Parameter(i));
    }

    // Give minimum of the fucntion with y, z, theta
    double yy, zz, thetaa;
    double chi2min = fit_func->GetMinimumXYZ(yy, zz, thetaa);

    cout << "Minimum Chi2: " << chi2min << " at (y, z, theta) = (" 
         << yy << ", " << zz << ", " << thetaa << ")" << endl;
}

