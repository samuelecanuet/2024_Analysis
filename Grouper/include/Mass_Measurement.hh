#include "Detectors.hh"
#include "Analyser_New.hh"
// #include "../../Simulation/include/Systematic.hh"

pair<double, double> ComputeEshift(TFile *file, int det)
{
    TCanvas *c = nullptr;
    TH1D *h_single = nullptr, *h_coinc = nullptr, *h_nocoinc = nullptr;
    if (string(file->GetName()).find("analysed_new") != std::string::npos)
    {
        c = (TCanvas *)file->Get(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_Coinc").c_str());
        h_single = (TH1D *)c->GetPrimitive(("Silicon_Detector_Energy_Deposit_SINGLE_" + detectorName[det] + "_All").c_str());
        h_coinc = (TH1D *)c->GetPrimitive(("Silicon_Detector_Energy_Deposit_COINC_" + detectorName[det] + "_All").c_str());
        h_nocoinc = (TH1D *)c->GetPrimitive(("Silicon_Detector_Energy_Deposit_NOCOINC_" + detectorName[det] + "_All").c_str());
    }
    else if (string(file->GetName()).find("result") != std::string::npos)
    {
        c = (TCanvas *)file->Get(("H_Single_" + detectorName[det]).c_str());
        h_single = (TH1D *)c->GetPrimitive(("H_Single_" + detectorName[det]).c_str());
        h_coinc = (TH1D *)c->GetPrimitive(("H_Coinc_" + detectorName[det]).c_str());
        h_nocoinc = (TH1D *)c->GetPrimitive(("H_NoCoinc_" + detectorName[det]).c_str());
    }
    else
    {
        Error("Unknown file format for file: " + string(file->GetName()));
    }

    if (h_single == nullptr || h_coinc == nullptr || h_nocoinc == nullptr)
    {
        Error("Histograms not found for detector: " + detectorName[det]);
    }

    return ComputeEshift(det, h_single, h_coinc, h_nocoinc);
}