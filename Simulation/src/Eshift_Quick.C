

#include "../../Grouper/include/Detectors.hh"
#include "../../Grouper/include/Analyser.hh"

int Eshift_Quick()
{
    FLAG2024 = true;
    NUCLEUS = "32Ar";
    IAS = 14;
    InitDetectors("../../Grouper/Config_Files/sample.pid");
    TFile *f = MyTFile((DIR_DATA_HDD + "../SIMULATED_DATA/05-16/32Ar_IAS_1mm_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
    
    InitWindows("../../Grouper/");


    TGraphErrors *G_Eshift = new TGraphErrors();
    TGraphErrors *G_Single = new TGraphErrors();
    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliStrip(det))
        {
            TCanvas *c = (TCanvas *)f->Get(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_Coinc").c_str());            
            TH1D *h_single = (TH1D *)c->GetPrimitive(("Silicon_Detector_Energy_Deposit_SINGLE_" + detectorName[det] + "_All").c_str());
            TH1D *h_coinc = (TH1D *)c->GetPrimitive(("Silicon_Detector_Energy_Deposit_COINC_" + detectorName[det] + "_All").c_str());
            TH1D *h_nocoinc = (TH1D *)c->GetPrimitive(("Silicon_Detector_Energy_Deposit_NOCOINC_" + detectorName[det] + "_All").c_str());

            pair<double, double> e = ComputeEshift(det, h_single, h_coinc, h_nocoinc);

            double eshift = e.first;
            double eshift_error = e.second;

            G_Eshift->AddPoint(det, eshift);
            G_Eshift->SetPointError(G_Eshift->GetN() - 1, 0, eshift_error);

            G_Single->AddPoint(det, h_single->Integral());
            G_Single->SetPointError(G_Single->GetN() - 1, 0, sqrt(h_single->Integral()));
        }
    }
    
    TFile *output = MyTFile(("Eshift_" + to_string(YEAR) + ".root").c_str(), "RECREATE");
    TCanvas *c = new TCanvas("Eshift", "Eshift", 800, 600); 
    c->cd();
    G_Eshift->SetTitle("Eshift");
    G_Eshift->GetXaxis()->SetTitle("Strip");
    G_Eshift->GetYaxis()->SetTitle("Eshift [keV]");
    G_Eshift->SetMarkerStyle(20);
    G_Eshift->SetMarkerSize(1);
    G_Eshift->Draw("AP");
    c->Write();

    TCanvas *c1 = new TCanvas("Single", "Single", 800, 600);
    c1->cd();
    G_Single->SetTitle("Single");
    G_Single->GetXaxis()->SetTitle("Strip");
    G_Single->GetYaxis()->SetTitle("Counts");
    G_Single->SetMarkerStyle(20);
    G_Single->SetMarkerSize(1); 
    G_Single->Draw("AP");
    c1->Write();

    output->Close();
    f->Close();
    
    return 0;

}