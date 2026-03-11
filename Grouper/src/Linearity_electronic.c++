#include "Linearity_electronic.hh"
#include "Analyser.hh"
int main()
{

    /////// DEFINING A TF1 CALLED Linearity_detname ///////
    // reading run 017 with pulser test 
    // detecting peaks, finding with gaussian, calibrating on a "new" observable near channel
    // function wrote in Linearity_electronic.root
    // -- as to be used before fitting in Calibration

    FLAG2025 = true;
    InitDetectors("Config_Files/sample.pid");

    // OUTPUT file
    TFile *fout = MyTFile(Form("%sLinearity_electronic.root", DIR_ROOT_DATA_CALIBRATED.c_str()), "RECREATE");

    int run = 17;
    InitHistogramsForRun(run);
    InitTGraphs();

    // Reading
    TFile *f = MyTFile(Form("%s/run_0%d_208Po_test_pulser.root", DIR_ROOT_DATA.c_str(), run), "READ");
    TTreeReader *Reader = new TTreeReader("Tree", f);
    TTreeReaderValue<Signal> *signal = new TTreeReaderValue<Signal>(*Reader, "Signal");

    Start("Reading data");
    while (Reader->Next())
    {
        int Label = (**signal).Label;

        if (IsDetectorSiliStrip(Label))
        {
            H_Channel[Label]->Fill((**signal).Channel/1000);
        }
    }

    // Fitting  
    Start("Fitting data");
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (!IsDetectorSiliStrip(i))
            continue;
        Info("Fining peaks...", 1);
        TSpectrum s;
        int nfound = s.Search(H_Channel[i], 10, "", 0.01);
        double *xpeaks = s.GetPositionX();
        sort(xpeaks, xpeaks + nfound);

        Info(Form("Fitting detector %d", i), 1);
        int x = 10;
        for (int j = 0; j < nfound; j++)
        {
            if (xpeaks[j] < 8)
                continue;
            Info(Form("Fitting peak at %.2f CH", xpeaks[j]), 2);
            H_Channel[i]->GetXaxis()->SetRangeUser(xpeaks[j]-1, xpeaks[j]+1);
            if (H_Channel[i]->Integral() == 0)
            {   
                continue;
            }
            TF1 *f1 = new TF1("f1", "gaus", xpeaks[j]-1, xpeaks[j]+1);
            f1->SetParameter(0, H_Channel[i]->GetMaximum());
            f1->SetParameter(1, H_Channel[i]->GetMean());
            f1->SetParameter(2, 35/1000);
            TFitResultPtr r = H_Channel[i]->Fit(f1, "RQS");

            if (r->Status() != 0)
            {
                Warning(Form("Fit failed for detector %d ", i), 1);
                continue;
            }

            G_Linearity[i]->AddPoint(x, f1->GetParameter(1));
            G_Linearity[i]->SetPointError(G_Linearity[i]->GetN()-1, 0, f1->GetParError(1));
            x+=10;
        }
        
        if (!G_Linearity[i]->GetN())
        {
            Warning(Form("No points for detector %d, skipping...", i), 1);
            continue;
        }
        fout->cd();
        TCanvas *c1 = new TCanvas(Form("c1_Det%d", i), Form("c1_Det%d", i), 800, 600);
        c1->Divide(1,3);
        c1->cd(1);
        G_Linearity[i]->SetTitle(Form("Linearity_Det%d", i));
        G_Linearity[i]->GetXaxis()->SetTitle("Pulser Value (V)");
        G_Linearity[i]->GetYaxis()->SetTitle("Energy");
        G_Linearity[i]->SetMarkerStyle(20);
        G_Linearity[i]->SetMarkerSize(1);
        G_Linearity[i]->Draw("AP");
        string func_name = "pol4";
        G_Linearity[i]->Fit(func_name.c_str(), "Q");
        c1->cd(3);
        H_Channel[i]->GetXaxis()->SetRangeUser(-1111, -1111);
        H_Channel[i]->Draw();

        double chi2_red = G_Linearity[i]->GetFunction(func_name.c_str())->GetChisquare() / G_Linearity[i]->GetFunction(func_name.c_str())->GetNDF();
        G_Global->AddPoint(i, chi2_red);   

        c1->cd(2); // residuals of tgraph fit
        TGraphErrors *g_res = new TGraphErrors();
        for (int p = 0; p < G_Linearity[i]->GetN(); p++)
        {
            double x, y;
            G_Linearity[i]->GetPoint(p, x, y);
            double y_fit = G_Linearity[i]->GetFunction(func_name.c_str())->Eval(x);
            g_res->SetPoint(p, x, y - y_fit);
            g_res->SetPointError(p, 0, G_Linearity[i]->GetErrorY(p));
        }
        g_res->SetTitle("Residuals");
        g_res->GetXaxis()->SetTitle("Pulser Value (V)");
        g_res->GetYaxis()->SetTitle("Residuals (Channels)");
        g_res->SetMarkerStyle(20);
        g_res->SetMarkerSize(1);
        g_res->Draw("AP");
        TLine *l = new TLine(G_Linearity[i]->GetXaxis()->GetXmin(), 0, G_Linearity[i]->GetXaxis()->GetXmax(), 0);
        l->SetLineColor(kRed);
        l->Draw("same");

        
        c1->cd(1);
        TLatex t;
        t.SetNDC();
        t.SetTextSize(0.04);
        t.DrawLatex(0.15, 0.85, Form("#chi^{2}_{#nu} = %.2f", chi2_red));
        c1->Write();

        TF1 *fout_func = (TF1*) G_Linearity[i]->GetFunction(func_name.c_str())->Clone(Form("Linearity_%s", detectorName[i].c_str()));
        fout_func->SetTitle(Form("Linearity_%s", detectorName[i].c_str()));
        fout_func->Write();
    }

    TCanvas *c2 = new TCanvas("Global", "Global", 800, 600);
    G_Global->SetTitle("Global_Chi2");
    G_Global->GetXaxis()->SetTitle("Detector Number");
    G_Global->GetYaxis()->SetTitle("Reduced Chi2");
    G_Global->SetMarkerStyle(20);
    G_Global->SetMarkerSize(1);
    G_Global->Draw("AP");
    c2->Write();

    fout->Close();   


    return 0;
}