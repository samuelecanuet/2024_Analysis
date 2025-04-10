// // #include <iostream>
// // #include <vector>
// // #include <cmath>
// // #include <complex>
// // #include <gsl/gsl_sf_coulomb.h>
// // #include "TMath.h"
// // #include "TH1.h"
// // #include "TCanvas.h"
// // #include "TGraph.h"
// // #include "TF1.h"
// // #include "TFile.h"
// // #include "RMATRIX.hh"

// // class RMatrix {
// // public:
// //     RMatrix(double Z1, double Z2, double reducedMass, double radius);
// //     double CalculateCrossSection(double energy, const std::vector<double>& parameters);
// //     double Penetrability(double energy);
// //     double PhaseShift(double energy);

// // private:
// //     double Z1, Z2;           // Charges of the particles
// //     double reducedMass;       // Reduced mass of the system
// //     double radius;            // Nuclear interaction radius
// // };

// // RMatrix::RMatrix(double Z1, double Z2, double reducedMass, double radius)
// //     : Z1(Z1), Z2(Z2), reducedMass(reducedMass), radius(radius) {}

// // double RMatrix::CalculateCrossSection(double energy, const std::vector<double>& parameters) {
// //     size_t nResonances = parameters.size() / 3; // Each resonance has 3 parameters (energy, width, branching ratio)
// //     std::complex<double> totalAmplitude(0.0, 0.0);

// //     for (size_t i = 0; i < nResonances; ++i) {
// //         double resonanceEnergy = parameters[i * 3];        // Resonance energy
// //         double gammaSquared = parameters[i * 3 + 1];       // Width
// //         double branchingRatio = parameters[i * 3 + 2];     // Branching ratio

// //         // Penetrability factor
// //         double penetrationFactor = Penetrability(energy);

// //         // Width, scaled by branching ratio
// //         double width = 2 * gammaSquared * penetrationFactor * branchingRatio;

// //         // Phase shift
// //         double delta = PhaseShift(energy);

// //         // Resonance amplitude (complex number: real part and imaginary part)
// //         std::complex<double> resonanceAmplitude(width, (energy - resonanceEnergy));

// //         // Apply phase shift as a complex exponential
// //         std::complex<double> phaseFactor = std::exp(std::complex<double>(0.0, -2 * delta));

// //         // Interference occurs when we sum the complex amplitudes
// //         totalAmplitude += resonanceAmplitude * phaseFactor;
// //     }

// //     // The cross-section is proportional to the square of the magnitude of the total amplitude
// //     double crossSection = std::norm(totalAmplitude); // norm() returns |totalAmplitude|^2
// //     return crossSection;
// // }

// // double RMatrix::Penetrability(double energy) {
// //     // Sommerfeld parameter: Z1 * Z2 * e^2 / (ħc)
// //     double eta = Z1 * Z2 * 1.4399764 / energy;

// //     double rho = sqrt(2 * reducedMass * energy) * radius;
// //     int L = 0;  // For s-wave (L=0)

// //     gsl_sf_result F_L, G_L, dF_L, dG_L;
// //     double a, b;

// //     // Call the GSL function correctly and check the status
// //     int status = gsl_sf_coulomb_wave_FG_e(eta, rho, L, 1, &F_L, &G_L, &dF_L, &dG_L, &a, &b);

// //     if (status != 0) {
// //         std::cerr << "Error: Coulomb wave function calculation failed!" << std::endl;
// //         return 0;
// //     }

// //     // Penetrability is the square of the modulus of the wave functions
// //     double penetrability = pow(F_L.val, 2) + pow(G_L.val, 2);
// //     return penetrability;
// // }


// // double RMatrix::PhaseShift(double energy) {
// //     // Sommerfeld parameter: Z1 * Z2 * e^2 / (ħc)
// //     double eta = Z1 * Z2 * 1.4399764 / energy;
    
// //     // Coulomb wave parameters
// //     double rho = sqrt(2 * reducedMass * energy) * radius;
// //     int L = 0;  // For s-wave (L=0)

// //     gsl_sf_result F_L, G_L, dF_L, dG_L;
// //     double a, b;

// //     // // Calculate Coulomb wave functions F_L and G_L
// //     // int status = gsl_sf_coulomb_wave_FG_e(eta, rho, L, L, &F_L, &G_L, &dF_L, &dG_L, &a, &b);

// //     // if (status != 0) {
// //     //     std::cerr << "Error: Coulomb wave function calculation failed!" << std::endl;
// //     //     return 0;
// //     // }

// //     // // Phase shift is the arctangent of the ratio of the regular and irregular Coulomb functions
// //     // double phaseShift = atan2(F_L.val, G_L.val);

// //     return 1;
// // }

// // // Global RMatrix object for use in the fitting function
// // RMatrix rmatrix(18, 1, 1.007276, 5.0); // Example for Z1=18 (Ar), Z2=1 (proton), reduced mass in u, radius in fm

// // // This function is the callable for ROOT's TF1
// // Double_t FitFunction(Double_t *x, Double_t *par) {
// //     // The energy at which we want to calculate the cross-section
// //     double energy = x[0];

// //     // Convert the ROOT parameters to a std::vector
// //     std::vector<double> parameters(par, par + 3 * (sizeof(par) / sizeof(par[0])));

// //     // Calculate the cross-section using the RMatrix class
// //     double crossSection = rmatrix.CalculateCrossSection(energy, parameters);

// //     return crossSection;
// // }

// // // Function to fit experimental data
// // void FitExperimentalHistogram(const char* filename) {
// //     // Open the ROOT file with the experimental histogram
// //     // TFile *file = TFile::Open(filename);
// //     // if (!file || file->IsZombie()) {
// //     //     std::cerr << "Error: Could not open file " << filename << std::endl;
// //     //     return;
// //     // }

// //     // // Get the experimental histogram
// //     // TH1 *hist = dynamic_cast<TH1*>(file->Get("histogram_name"));
// //     // if (!hist) {
// //     //     std::cerr << "Error: Could not find histogram in file " << filename << std::endl;
// //     //     return;
// //     // }
// //     TFile *f = new TFile("data.root", "RECREATE");
// //     // Create the TF1 for fitting
// //     int nResonances = 3; // Example: 3 resonances
// //     TF1 *fitFunc = new TF1("fitFunc", FitFunction, 0, 5000, 1);

// //     // Set initial guess parameters: resonance energies, widths, and branching ratios
// //     fitFunc->SetParameters(2e3, 0.1, 0.8); // Example initial parameters

// //     // // Perform the fit
// //     fitFunc->SetNpx(1000); // Number of points for the fit
// //     fitFunc->Write();


// //     // Draw the result
// //     TCanvas *c1 = new TCanvas("c1", "Fit Result", 800, 600);
// //     // hist->Draw();
// //     fitFunc->Draw("same");

// //     // Save the fit result
// //     c1->SaveAs("fit_result.png");

// //     // Close the file
// //     f->Close();
// // }

// #include <gsl/gsl_sf_coulomb.h>
// #include <gsl/gsl_sf.h>
// #include <cmath>
// #include <iostream>
// #include <vector>
// using namespace std;


// double Penetrability(double eta, double rho, double L) {
//     gsl_sf_result F;
//     gsl_sf_result G;
//     gsl_sf_result Fp;
//     gsl_sf_result Gp;
//     double exp_F;
//     double exp_G;
//     int status = gsl_sf_coulomb_wave_FG_e(eta, rho, L, 0, &F, &Fp, &G, &Gp, &exp_F, &exp_G);

//     double Fval = F.val;
//     double Gval = G.val;

//     if (status == GSL_EOVRFLW)
//     {
//         Fval = exp(exp_F);
//         Gval = exp(exp_G);
//     }

//     return rho / (pow(Fval, 2) + pow(Gval, 2));
// }


// double HardSpherePhaseShift(double eta, double rho, double L)
// {
//     gsl_sf_result F;
//     gsl_sf_result G;
//     gsl_sf_result Fp;
//     gsl_sf_result Gp;
//     double exp_F;
//     double exp_G;
//     int status = gsl_sf_coulomb_wave_FG_e(eta, rho, L, 0, &F, &Fp, &G, &Gp, &exp_F, &exp_G);

//     double Fval = F.val;
//     double Gval = G.val;

//     if (status == GSL_EOVRFLW)
//     {
//         Fval = exp(exp_F);
//         Gval = exp(exp_G);
//     }

//     return -atan(Fval / Gval);
// }

// double ShiftFactor(double eta, double rho, double L) 
// {
//     gsl_sf_result F;
//     gsl_sf_result G;
//     gsl_sf_result Fp;
//     gsl_sf_result Gp;
//     double exp_F;
//     double exp_G;
//     int status = gsl_sf_coulomb_wave_FG_e(eta, rho, L, 0, &F, &Fp, &G, &Gp, &exp_F, &exp_G);

//     double Fval = F.val;
//     double Gval = G.val;
//     double Fpval = Fp.val;
//     double Gpval = Gp.val;

//     if (status == GSL_EOVRFLW)
//     {
//         Fval = exp(exp_F);
//         Gval = exp(exp_G);
//     }

//     double P = Penetrability(eta, rho, L);

//     return P * (Fval * Fp.val - Gval * Gp.val);
    
// }

// double PhaseShift(double E, double eta, double rho, double L, double gamma_nl, double E_nl) 
// {
//     double hard_sphere = HardSpherePhaseShift(eta, rho, L);
//     double penetrability = Penetrability(eta, rho, L);
//     double shift_factor = ShiftFactor(eta, rho, L); 

//     return hard_sphere - atan( ( pow(gamma_nl, 2) * penetrability) / (E_nl - pow(gamma_nl, 2) * shift_factor - E) );
// }

// double R(double E, vector<double> gamma_nl, vector<double> E_nl)
// {
//     double res = 0;
//     for (int i = 0; i < gamma_nl.size(); i++)
//     {
//         res += pow(gamma_nl[i], 2) / (E_nl[i] - E);
//     }
//     return res;
// }

// double SquaredTerm(double E, double eta, double rho, double L, vector<double> gamma_nl, vector<double> E_nl, vector<double> A_nl, double B)
// {
//     double res = 0;
    
//     double nominator = 0;
//     for (int i = 0; i < gamma_nl.size(); i++)
//     {
//         nominator += A_nl[i] / (E_nl[i] - E);
//     }

//     double denominator = 1 - ( ShiftFactor(eta, rho, L) - B + 


//     res = pow (abs(nominator/denominator), 2);

//     return res;
// }

// double Integrated_Fermi_Function(double E)
// {
//     return 1;
// }

// double Weight_Alpha(double E, double eta, double rho, vector<double> L_vec, vector<double> gamma_nl, vector<vector<double>> E_nl, vector<vector<double>> A_nl, vector<double> B_vec)
// {
//     double res = 0;
    
//     for (int i = 0; i < L_vec.size(); i++)
//     {
//         res += Penetrability(eta, rho, L_vec[i]) * SquaredTerm(E, eta, rho, L_vec[i], gamma_nl, E_nl[i], A_nl[i], B_vec[i]);
//     }

//     return Integrated_Fermi_Function(E) * res;
// }

// int main() 
// {
//     double Z1 = 18;
//     double Z2 = 1;
//     double e = 1.60217662e-19;
//     double hbar = 1.0545718e-34;
//     double a = 5.6; //fm


//     double Energy = 2e3; //keV


//     double eta = Z1 * Z2 * e * e / Energy / hbar;





//     return 0;
// }

