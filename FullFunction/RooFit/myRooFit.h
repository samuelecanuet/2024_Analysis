#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooFit.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"

#include "../../Grouper/include/Utilities.hh"

const double m = 931000.5;              // keV/c^2
const double m_e = 510.9989461;        // Electron mass in keV/c^2
const double m_p = 938272.0813;         // Proton mass in keV/c^2
const double c = 299792458;             // Speed of light in m/s
const double alpha = 1 / 137.035999084; // Fine-structure constant

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// CLASS PEAK SHAPE //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
class PeakShape : public RooAbsPdf
{
public:
    PeakShape(const char *name, const char *title,
              RooRealVar &E, RooRealVar &PeakN, RooRealVar &E_mean,
              RooRealVar &Qb, RooRealVar &A,
              RooRealVar &Z, RooRealVar &a, RooRealVar &Phi_min,
              RooRealVar &Phi_max, RooRealVar &Halflife)
        : RooAbsPdf(name, title),
          E_("E", "Observable", this, E),
          PeakN_("PeakN", "Peak number", this, PeakN),
          E_mean_("E_mean", "Peak mean", this, E_mean),
          Qb_("Qb", "Qbeta", this, Qb),
          A_("A", "Mass number", this, A),
          Z_("Z", "Atomic number", this, Z),
          a_("a", "Fermi parameter", this, a),
          Phi_min_("Phi_min", "Minimum angle", this, Phi_max),
          Phi_max_("Phi_max", "Maximum angle", this, Phi_min),
          HalfLife_("HalfLife", "Half-life", this, Halflife) // Default value for Half-life

    {
        cout << "PeakShape created with parameters:" << endl;
        cout << "E: " << E_mean_ << " keV" << endl;
        cout << "PeakN: " << PeakN_ << endl;
        cout << "Qb: " << Qb_ << " keV" << endl;
        cout << "A: " << A_ << endl;
        cout << "Z: " << Z_ << endl;
        cout << "a: " << a_ << endl;
        cout << "Phi_min: " << Phi_min_ << " rad" << endl;
        cout << "Phi_max: " << Phi_max_ << " rad" << endl;
        cout << "HalfLife: " << HalfLife_ << " keV" << endl;
        cout << "K: " << K_ << endl;
        cout << "----------------------------------------" << endl;

        BW = new RooBreitWigner("BreitWigner", "Breit-Wigner PDF", E, E_mean, Halflife);
        

    }

    PeakShape(const PeakShape &other, const char *name = nullptr)
        : RooAbsPdf(other, name),
          E_("E", this, other.E_), PeakN_("PeakN", this, other.PeakN_),
          E_mean_("E_mean", this, other.E_mean_),
          Qb_("Qb", this, other.Qb_), A_("A", this, other.A_), Z_("Z", this, other.Z_),
          a_("a", this, other.a_), Phi_min_("Phi_min", this, other.Phi_min_),
          Phi_max_("Phi_max", this, other.Phi_max_), HalfLife_("HalfLife", this, other.HalfLife_)
    {
    }

    TObject *clone(const char *newname) const override
    {
        return new PeakShape(*this, newname);
    }

    double evaluate() const override
    {
        return fTotalWeightKinematicShift(E_, 1, E_mean_, Qb_, K(E_mean_, A_*m), a_, A_, Z_, Phi_min_, Phi_max_);
    }

    double Test(double x, double a) const
    {
        return a * x;
    }

    /// COMMON FUNCTION ///
    double EnergyToMomentum(double E, double Mass) const
    {
        return sqrt(E * E - Mass * Mass);
    }

    double MomentumToEnergy(double p, double Mass) const
    {
        return sqrt(p * p + Mass * Mass);
    }

    double EnergyToVelocity(double E, double Mass) const
    {
        return E / Mass;
    }

    /// KINEMATICS ///

    double K(double Ep, double M) const 
    {
        return sqrt(2 * Ep / M * (m_p * (M - m_p)) / (M * M));
    }

    double Fermi(int Z, double E_b, double M) const
    {
        double v = EnergyToVelocity(E_b, M);
        double nu = Z * alpha / v * c;
        return 2 * M_PI * nu / (1 - exp(-2 * M_PI * nu));
    }

    double fWmin(double W0, double k, double t) const
    {
        if (t / k < -(W0 - m_e))
            return (m_e * m_e + (W0 + t / k) * (W0 + t / k)) / (2 * (W0 + t / k));
        else if (abs(t / k) <= W0 - m_e)
            return m_e;
        else
            return (m_e * m_e + (W0 - t / k) * (W0 - t / k)) / (2 * (W0 - t / k));
    }

    double WeightKinematicShift(double t, double k, double W0, double W, double a, double M, double Z, double Phi_min, double Phi_max) const
    {
        double p_b = EnergyToMomentum(W, m_e);
        if (p_b == 0)
            return 0;

        double c1 = min(cos(Phi_max), max(cos(Phi_min), (W - W0 - t / k) / p_b));
        double c2 = max(cos(Phi_min), min(cos(Phi_max), (W0 - W - t / k) / p_b));

        double I1 = (c2 - c1) * W * (W0 - W);
        double I2a = -a * (pow(c2, 2) - pow(c1, 2)) / 2.0 * p_b * t / k;
        double I2b = -a * (pow(c2, 3) - pow(c1, 3)) / 3.0 * pow(p_b, 2);

        return Fermi(Z, W, M) * p_b / (2 * k) * (I1 + I2a + I2b) * 1e-28;
    }

    double fTotalWeightKinematicShift(double t, double AA, double E_p, double W0, double k, double a, double A, double Z, double Phi_min, double Phi_max) const
    {
        // double AA = par[0];
        // double E_p = par[1];
        // double W0 = par[2];
        // double k = par[3];
        // double a = par[4];
        // double A = par[5];
        // double Z = par[6];
        // double Phi_min = par[7];
        // double Phi_max = par[8];

        double t_ = t - E_p;
        double M = A * m;

        // cout << abs(t_ / k) << " " << sqrt(W0 * W0 - m_e * m_e) << endl;

        if (abs(t_ / k) > sqrt(W0 * W0 - m_e * m_e))
        {
            return 0;
        }

        // cout << "ok" << endl;

        double wt = 0;
        double Wmin = fWmin(W0, k, t_);
        double step = (W0 - Wmin) / 1000.0;

        // cout << "ok" << endl;

        for (double Wb = Wmin; Wb <= W0; Wb += step)
        {
            wt += AA * WeightKinematicShift(t_, k, W0, Wb, a, M, Z, Phi_min, Phi_max) * step;
        }

        // cout << "Total weight for E = " << t << ": " << wt << endl;

        return wt;
    }

    RooBreitWigner *GetBreitWigner()
    {
        return BW;
    }

private:
    RooRealProxy E_;        // VARIABLE
    RooRealProxy PeakN_;    // Peak Number
    RooRealProxy E_mean_;   // Mean energy of the peak
    RooRealProxy Qb_;       // Q_beta value
    RooRealProxy A_;        // Mass number
    RooRealProxy Z_;        // Atomic number
    RooRealProxy a_;        // Fermi parameter
    RooRealProxy Phi_min_;  // Minimum angle
    RooRealProxy Phi_max_;  // Maximum angle
    RooRealProxy HalfLife_; // Nuclear level LifeTime

    RooRealVar K_; // Kinematic factor

    RooBreitWigner *BW;

};

/////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// CLASS PEAK x BREITWIGNER //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
class Peak :
{
public:
    Peak(const char *name, const char *title,
         RooRealVar &E, RooRealVar &PeakN, RooRealVar &E_mean,
         RooRealVar &HalfLife, RooRealVar &Qb, RooRealVar &A,
         RooRealVar &Z, RooRealVar &a, RooRealVar &Phi_min,
         RooRealVar &Phi_max)
        : RooAbsPdf(name, title),
          E_("E", "Observable", this, E),
          PeakN_("PeakN", "Peak number", this, PeakN),
          E_mean_("E_mean", "Peak mean", this, E_mean),
          HalfLife_("HalfLife", "Half-life", this, HalfLife),
          Qb_("Qb", "Qbeta", this, Qb),
          A_("A", "Mass number", this, A),
          Z_("Z", "Atomic number", this, Z),
          a_("a", "Fermi parameter", this, a),
          Phi_min_("Phi_min", "Minimum angle", this, Phi_min),
          Phi_max_("Phi_max", "Maximum angle", this, Phi_max)

    {

        cout << "Peak created with parameters:" << endl;
        cout << "E: " << E_mean_ << " keV" << endl;
        cout << "PeakN: " << PeakN_ << endl;
        cout << "HalfLife: " << HalfLife_ << " keV" << endl;
        cout << "Qb: " << Qb_ << " keV" << endl;
        cout << "A: " << A_ << endl;
        cout << "Z: " << Z_ << endl;
        cout << "a: " << a_ << endl;
        cout << "Phi_min: " << Phi_min_ << " rad" << endl;
        cout << "Phi_max: " << Phi_max_ << " rad" << endl;
        cout << "----------------------------------------" << endl;
        // Build the components and convolution
        peakShape_ = new PeakShape(Form("%s_peakShape", name), Form("%s Peak Shape", title),
                                   (RooRealVar &)E, (RooRealVar &)PeakN, (RooRealVar &)E_mean,
                                   (RooRealVar &)Qb, (RooRealVar &)A, (RooRealVar &)Z,
                                   (RooRealVar &)a, (RooRealVar &)Phi_min, (RooRealVar &)Phi_max);

        zeroMean = new RooRealVar("zeroMean", "Zero mean", 0);
        breitWigner_ = new RooBreitWigner(Form("%s_breitWigner", name), Form("%s Breit-Wigner", title),
                                         (RooRealVar &)E, *zeroMean, (RooRealVar &)HalfLife);

        conv_ = new RooFFTConvPdf(Form("%s_conv", name), Form("%s convolution", title),
                                  (RooRealVar &)E, *peakShape_, *breitWigner_);
    }

    Peak(const Peak &other, const char *name = nullptr)
        : RooAbsPdf(other, name),
          E_("E", this, other.E_), PeakN_("PeakN", this, other.PeakN_),
          E_mean_("E_mean", this, other.E_mean_), HalfLife_("HalfLife", this, other.HalfLife_),
          Qb_("Qb", this, other.Qb_), A_("A", this, other.A_), Z_("Z", this, other.Z_),
          a_("a", this, other.a_), Phi_min_("Phi_min", this, other.Phi_min_),
          Phi_max_("Phi_max", this, other.Phi_max_)
    {
        // Must rebuild sub-PDFs with correct references
        peakShape_ = new PeakShape(*other.peakShape_);
        breitWigner_ = new RooBreitWigner(*other.breitWigner_);
        conv_ = new RooFFTConvPdf(*other.conv_);
    }

    virtual ~Peak()
    {
        delete peakShape_;
        delete breitWigner_;
        delete conv_;
    }

    TObject *clone(const char *newname) const override
    {
        return new Peak(*this, newname);
    }

    RooAbsPdf* GetPDF() const
    {
        if (HalfLife_ == 0) {
            return peakShape_;
        }
        return conv_.get();
    }

    double GetBR() const
    {
        // Assuming all peaks have a branching ratio of 1.0
        return n1.0;
    }

private:
    RooRealProxy E_;       // VARIABLE
    RooRealProxy PeakN_;   // Peak Number
    RooRealProxy E_mean_;  // Mean energy of the peak
    RooRealProxy HalfLife_; // Nuclear level LifeTime
    RooRealProxy Qb_;      // Q_beta value
    RooRealProxy A_;       // Mass number
    RooRealProxy Z_;       // Atomic number
    RooRealProxy a_;       // Fermi parameter
    RooRealProxy Phi_min_; // Minimum angle
    RooRealProxy Phi_max_; // Maximum angle

    PeakShape *peakShape_;
    RooBreitWigner *breitWigner_;
    RooFFTConvPdf *conv_;

    RooRealVar *zeroMean;
};
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

class Nucleus:
{
public:
    Nucleus(string name, RooRealVar &E)
        : name_(name), E_(E)
    {
        cout << "Nucleus created with name: " << name_ << endl;

        int PDG_code = NametoCode_map(name);
        if (PDG_code == 0) {
            cout << "Error: Invalid nucleus name: " << name_ << endl;
            exit(0);
        }
        A_ = RooRealVar("A", "Mass number", (PDG_code - 100000000) / 1000);
        Z_ = RooRealVar("Z", "Atomic number", (PDG_code - 100000000) % 1000);

        if (name == "32Ar")
            IAS = 14;
        else if (name = "32Cl")
            IAS = 21;
        
    }

    void LoadData()
    {
        cout << "Loading data for nucleus: " << name_ << endl;
        
        ifstream file("../../Simulation/CRADLE_FILE_GENERATOR/" + name_ + "proton.txt");
        if (!file.is_open()) {
            cout << "Error: Could not open file for nucleus: " << name_ << endl;
            exit(0);
        }
        string line;
        while (getline(file, line)) {
            istringstream iss(line);
            double code, energy, energy_err, intensity, intensity_err;
            iss >> code >> energy >> energy_err >> intensity >> intensity_err;
            if (iss.fail()) {
                cout << "Error: Invalid data format in file for nucleus: " << name_ << endl;
                exit(0);
            }

            double a = code == IAS ? 1.0 : -1./3.; // Fermi parameter, set to 1 for IAS, 0 otherwise

            Peak *peak = new Peak(Form("%s_peak_%f", name_.c_str(), code),
                                   Form("%s Peak %f", name_.c_str(), code),
                                   E_, RooRealVar("PeakN", "Peak number", code),
                                   RooRealVar("E_mean", "Energy mean", energy, energy-5*energy_err, energy+5*energy_err),
                                   RooRealVar("HalfLife", "Half-life", 20e-3), // Default value for Half-life
                                   RooRealVar("Qb", "Qbeta", sqrt(pow(5046 + 1022 - 511, 2) + pow(511, 2))),
                                   A_, Z_,
                                   RooRealVar("a", "Fermi parameter", a),
                                   RooRealVar("Phi_min", "Minimum angle", 0),
                                   RooRealVar("Phi_max", "Maximum angle", TMath::Pi()));
            AddPeak(code, peak);
        }
        file.close();
    }

    void AddPeak(double code, Peak *peak)
    {
        if (peaks_.find(code) != peaks_.end()) {
            cout << "Peak with code " << code << " already exists in nucleus: " << name_ << endl;
            exit(0);
        }
        peaks_[code] = peak;
        cout << "Peak with code " << code << " added to nucleus: " << name_ << endl;
    }


    RooAbsPdf* GetPDF()
    {
        // Summing in a new pDF all the peaks
        RooArgList pdfList;
        RooArgList BRList;
        for (const auto &pair : peaks_) {
            double code = pair.first;
            Peak *peak = pair.second;
            pdfList.add(*peak->GetPDF());  
            BRList.add(RooRealVar(Form("BR_%f", code), Form("Branching Ratio for code %f", code), *peak->GetBR())); // Assuming all peaks have a branching ratio of 1.0 
        }

        if (pdfList.getSize() == 0) {
            cout << "Error: No peaks available in nucleus: " << name_ << endl;
            exit(0);
        }

        if (pdfList.getSize() == 1) {
            Total_ = (RooAbsPdf*) pdfList.first();
        } else {
            Total_ = new RooAddPdf(Form("%s_Total", name_.c_str()), Form("%s Total PDF", name_.c_str()), pdfList, BRList);
        }

        return Total_.get();
    }



private:
    string name_;   // Name of the nucleus
    RooRealVar &E_; // Energy observable
    RooRealVar &A_; // Mass number
    RooRealVar &Z_; // Atomic number
    int IAS = 0; // Isobaric Analog State (IAS) code


    map<double, PeakShape*> peaks_; // Map of peaks by code

    std::unique_ptr<RooFFTConvPdf> Total_;
};
