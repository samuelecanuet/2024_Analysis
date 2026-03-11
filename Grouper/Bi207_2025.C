// plotting BI207 runs in 2025 with different col 5mm, 10mm and CF

#include "include/Detectors.hh"
#include "include/Utilities.hh"

int Bi207_2025()
{

    FLAG2025 = true;
    InitDetectors("Config_Files/sample.pid");
    
    // 116 - 10mm
    TFile *file_116 = MyTFile((DIR_ROOT_DATA_GROUPED + "run_116_data_207Bi_4T_grouped.root").c_str(), "READ");
    double time116 = Diff_Date(GetTime(file_116).first, GetTime(file_116).second, FLAG2021);

    // 117 - 5mm
    TFile *file_117 = MyTFile((DIR_ROOT_DATA_GROUPED + "run_117_data_207Bi_4T_small_collimator_grouped.root").c_str(), "READ");
    double time117 = Diff_Date(GetTime(file_117).first, GetTime(file_117).second, FLAG2021);

    // 118 - CF
    TFile *file_118 = MyTFile((DIR_ROOT_DATA_GROUPED + "run_118_data_207Bi_4T_FC_grouped.root").c_str(), "READ");
    double time118 = Diff_Date(GetTime(file_118).first, GetTime(file_118).second, FLAG2021);                                                                        


    TFile *fout = MyTFile((DIR_ROOT_DATA_GROUPED + "Comparison_207Bi_4T.root").c_str(), "RECREATE");

    map<string, TH1D*> files;
    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorBetaHigh(det))
        {
            TCanvas *chigh = new TCanvas(("c_" + detectorName[det] + "_High").c_str(), ("c_" + detectorName[det] + "_High").c_str(), 1920, 1080);
        
            cout << det << "   " << GetDetectorChannel(det) << endl;

            TH1D *chigh116 = (TH1D*)file_116->Get(Form("SiPM_High/SiPMHigh_Channel/BetaHi%d/Channel_Cleaned_BetaHi%d", GetDetectorChannel(det), GetDetectorChannel(det)));
            TH1D *chigh117 = (TH1D*)file_117->Get(Form("SiPM_High/SiPMHigh_Channel/BetaHi%d/Channel_Cleaned_BetaHi%d", GetDetectorChannel(det), GetDetectorChannel(det)));
            chigh117->Scale(time116 / time117);
            TH1D *chigh118 = (TH1D*)file_118->Get(Form("SiPM_High/SiPMHigh_Channel/BetaHi%d/Channel_Cleaned_BetaHi%d", GetDetectorChannel(det), GetDetectorChannel(det)));
            chigh118->Scale(time116 / time118);
            chigh116->SetLineColor(kRed);
            chigh117->SetLineColor(kBlue);
            chigh118->SetLineColor(kGreen);
            chigh116->Draw("HIST");
            chigh117->Draw("SAME");
            chigh118->Draw("SAME");
            chigh->BuildLegend();
            chigh->Write();

            TCanvas *clow = new TCanvas(("c_" + detectorName[det] + "_Low").c_str(), ("c_" + detectorName[det] + "_Low").c_str(), 1920, 1080);
            TH1D *clow116 = (TH1D*)file_116->Get(Form("SiPMLow/SiPMLow_Channel/BetaLo%d/Channel_Cleaned_BetaLo%d", GetDetectorChannel(det), GetDetectorChannel(det)));
            TH1D *clow117 = (TH1D*)file_117->Get(Form("SiPMLow/SiPMLow_Channel/BetaLo%d/Channel_Cleaned_BetaLo%d", GetDetectorChannel(det), GetDetectorChannel(det)));
            clow117->Scale(time116 / time117);
            TH1D *clow118 = (TH1D*)file_118->Get(Form("SiPMLow/SiPMLow_Channel/BetaLo%d/Channel_Cleaned_BetaLo%d", GetDetectorChannel(det), GetDetectorChannel(det)));
            clow118->Scale(time116 / time118);
            clow116->SetLineColor(kRed);
            clow117->SetLineColor(kBlue);
            clow118->SetLineColor(kGreen);
            clow116->Draw("HIST");
            clow117->Draw("SAME");
            clow118->Draw("SAME");
            clow->BuildLegend();
            clow->Write();
        }
    }

    fout->Close();


    return 0;
}