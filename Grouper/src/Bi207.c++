#include "Bi207.hh"

int main()
{
    InitDetectors("Config_Files/sample.pid");
    TFile *ROOT_File = new TFile((DIR_ROOT_DATA+"run_137_multifast_207Bi.root").c_str(), "READ");
    TTreeReader *Reader = new TTreeReader((TTree *)ROOT_File->Get("Tree_Group"));
    TTreeReaderArray<Signal> signals(*Reader, "Signal");

    TFile *ffile = new TFile((DIR_ROOT_DATA_CALIBRATED + "SiPM_Calibrated.root").c_str(), "READ");
    for (int i = 1; i <= 9; i++)
    {
        f_matching[i] = (TF1*)ffile->Get(("SiPM_"+to_string(i) + "/MatchingSiPM_32Ar_SiPM" + to_string(i)).c_str());
    }

    TFile *Outfile = new TFile((DIR_ROOT_DATA_GROUPED + "Bi207.root").c_str(), "RECREATE");

    // Init histograms
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBeta(i))
        {
            H_SiPM[i] = new TH1D(Form("H_SiPM_%d", i), Form("H_SiPM_%d", i), eHighN, eHighMin, eHighMax);
        }
    }

    // read data
    clock_t start = clock(), Current;
    int entries = Reader->GetEntries();
    cout << "Selecting Groups : " << endl;
    while (Reader->Next() && Reader->GetCurrentEntry() < 10000000)
    {
        ProgressBar(Reader->GetCurrentEntry(), entries, start, Current, "Selecting Groups : ");
        for (Signal &signal : signals)
        {
            if (IsDetectorBeta(signal.Label))
            {
                
                H_SiPM[signal.Label]->Fill(f_matching[GetDetectorChannel(signal.Label)]->Eval(signal.Channel));
            }
        }
    }

    //fit histograms
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBetaHigh(i))
        {
            det = GetDetectorChannel(i);
            FitHistogram(H_SiPM[i]);
            g_calib[i] = new TGraph();
            g_calib[i]->SetPoint(0, f_fitted[det]->GetParameter(1), 569.698-100);
            g_calib[i]->SetPoint(1, f_fitted[det]->GetParameter(4), 1063.656-100);
            TF1 *f_calib = new TF1("f_calib", "[0]*x", 0, 1200e3);
            g_calib[i]->Fit(f_calib, "Q");
            g_calib[i]->Write();
            g_calib[i]->GetFunction("f_calib")->Write();
            cout << g_calib[i]->GetFunction("f_calib")->Eval(40e3) << endl;
            cout << g_calib[i]->GetFunction("f_calib")->GetParameter(0) << endl;

        }
    }

    // write histograms
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBeta(i))
        {
            TCanvas *c = new TCanvas(Form("c_SiPM_%d", i), Form("c_SiPM_%d", i), 800, 600);
            H_SiPM[i]->GetXaxis()->SetRangeUser(0, 1200e3);
            H_SiPM[i]->GetYaxis()->SetRangeUser(0, 8000);
            H_SiPM[i]->Draw("HIST");
            f_fitted[GetDetectorChannel(i)]->SetLineColor(kRed);
            f_fitted[GetDetectorChannel(i)]->Draw("SAME");
            c->Write();
        }
    }

    Outfile->Close();


}