#include "Bi207.hh"

int main()
{
    InitDetectors("Config_Files/sample.pid");
    TFile *ROOT_File = new TFile((DIR_ROOT_DATA+"run_137_multifast_207Bi.root").c_str(), "READ");
    TTreeReader *Reader = new TTreeReader((TTree *)ROOT_File->Get("Tree_Group"));
    TTreeReaderArray<Signal> signals(*Reader, "Signal");

    TFile *ffile = new TFile((DIR_ROOT_DATA_CALIBRATED + "SiPM_Calibrated_114_QDC2.root").c_str(), "READ");
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
            
            H_SiPM_time[i] = new TH2D(Form("H_SiPM_time_%d", i), Form("H_SiPM_time_%d", i), 500, 0, 1000, eHighN/100, eHighMin, eHighMax);
        }
        if (IsDetectorBetaHigh(i))
        {
            H2D_SiPM_HighLow[GetDetectorChannel(i)] = new TGraph();
        }
    }

    // read data
    clock_t start = clock(), Current;
    int entries = Reader->GetEntries();
    cout << "Selecting Groups : " << endl;
    while (Reader->Next() && Reader->GetCurrentEntry() < 10000000)
    {
        ProgressBar(Reader->GetCurrentEntry(), entries, start, Current, "Selecting Groups : ");

        Signal SiPM_High[10] = {Signal()};
        Signal SiPM_Low[10] = {Signal()};
        for (Signal &signal : signals)
        {

            if (IsDetectorBeta(signal.Label))
            {
                H_SiPM[signal.Label]->Fill(f_matching[GetDetectorChannel(signal.Label)]->Eval(signal.Channel));
            }

            if (IsDetectorBetaHigh(signal.Label))
            {
                SiPM_High[GetDetectorChannel(signal.Label)] = signal;
            }
            if (IsDetectorBetaLow(signal.Label))
            {
                SiPM_Low[GetDetectorChannel(signal.Label)] = signal;
            }
        }

        for (int det = 1; det <= 9; det++)
        {
            if (SiPM_High[det].isValid && SiPM_Low[det].isValid)
            {
                H2D_SiPM_HighLow[det]->AddPoint(SiPM_Low[det].Channel, SiPM_High[det].Channel);
            }
        }

        // time
        double time_ref;
        int counter = 0;
        for (Signal &signal : signals)
        {

            if (IsDetectorBetaHigh(signal.Label))
            {
                if (counter == 0)
                {
                    time_ref = signal.Time;
                    counter++;
                    continue;
                }
                H_SiPM_time[signal.Label]->Fill(signal.Time - time_ref, signal.Channel);
            }
        }
    }

    // fit histograms
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

            TCanvas *c = new TCanvas(Form("c_SiPM_%d", i), Form("c_SiPM_%d", i), 800, 600);
            H2D_SiPM_HighLow[det]->Draw("AP");
            c->Write();

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

            H_SiPM_time[i]->Write();
        }
    }

    Outfile->Close();


}