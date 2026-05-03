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

    InitCalib();

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

    //ref
    int det_ref = 12;
    bool ref_found = false;

    Start("Fitting data");
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (i == 0) i = det_ref;
        if (!IsDetectorSiliStrip(i))
            continue;
        Info("Finding peaks...", 1);
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
            f1->SetParameter(2, 35./1000);
            TFitResultPtr r = H_Channel[i]->Fit(f1, "RQS");

            if (r->Status() != 0)
            {
                Warning(Form("Fit failed for detector %d ", i), 1);
                continue;
            }

            G_Linearity[i]->AddPoint(f1->GetParameter(1), x);
            G_Linearity[i]->SetPointError(G_Linearity[i]->GetN()-1, f1->GetParError(1), 0);
            
            
            xpeaks[j] /= 10;
            int h = round(xpeaks[j]);
            if (G_Global_FWHM[h] == nullptr)
                G_Global_FWHM[h] = new TGraphErrors();
            G_Global_FWHM[h]->AddPoint(i, f1->GetParameter(2)*2.355*Calibration[i]->GetParameter(1));
            G_Global_FWHM[h]->SetPointError(G_Global_FWHM[h]->GetN()-1, 0, f1->GetParError(2)*2.355*Calibration[i]->GetParameter(1));
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
        TF1 *fit = new TF1(func_name.c_str(), func_name.c_str(), 0, 100);
        // fit->FixParameter(0, 0);
        G_Linearity[i]->Fit(func_name.c_str(), "QW");
        c1->cd(3);
        H_Channel[i]->GetXaxis()->SetRangeUser(-1111, -1111);
        H_Channel[i]->Draw();

        double chi2_red = G_Linearity[i]->GetFunction(func_name.c_str())->GetChisquare() / G_Linearity[i]->GetFunction(func_name.c_str())->GetNDF();
        G_Global->AddPoint(i, chi2_red);   

        c1->cd(2); // residuals of tgraph fit
        TGraphErrors *g_res = new TGraphErrors();
        double maxi = 0;
        double maxi_err = 0;
        for (int p = 0; p < G_Linearity[i]->GetN(); p++)
        {
            double x, y;
            G_Linearity[i]->GetPoint(p, x, y);
            double y_fit = G_Linearity[i]->GetFunction(func_name.c_str())->Eval(x);
            g_res->SetPoint(p, x, y - y_fit);
            g_res->SetPointError(p, 0, G_Linearity[i]->GetErrorY(p));
            if (abs(y - y_fit) > maxi)
            {
                maxi = abs(y - y_fit)*Calibration[i]->GetParameter(1); // in keV
                maxi_err = G_Linearity[i]->GetErrorY(p)*Calibration[i]->GetParameter(1);
            }
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

        G_Global_Delta->AddPoint(i, maxi);
        G_Global_Delta->SetPointError(G_Global_Delta->GetN()-1, 0, maxi_err);

        
        c1->cd(1);
        // TLatex t;
        // t.SetNDC();
        // t.SetTextSize(0.04);
        // t.DrawLatex(0.15, 0.85, Form("#chi^{2}_{#nu} = %.2f", chi2_red));
        c1->Write();

       

        
        // matching 
        G_Electronic_Matching[i] = new TGraphErrors();
        for (int p = 0; p < G_Linearity[det_ref]->GetN(); p++)
        {
            if (p >= G_Linearity[i]->GetN()-1 || p > 7)
                break;
            double x_ref, y_ref;
            G_Linearity[det_ref]->GetPoint(p, x_ref, y_ref);
            double x_ref_err = G_Linearity[det_ref]->GetErrorX(p);
            double x, y;
            G_Linearity[i]->GetPoint(p, x, y);
            G_Electronic_Matching[i]->AddPoint(x, x_ref);
            double x_err = G_Linearity[i]->GetErrorX(p);
            G_Electronic_Matching[i]->SetPointError(G_Electronic_Matching[i]->GetN()-1, x_err, x_ref_err);
        }
        //

        TCanvas *c_match = new TCanvas(Form("Matching_Det%d", i), Form("Matching_Det%d", i), 800, 600);
        c_match->cd();
        G_Electronic_Matching[i]->Fit("pol1", "Q");
        TPad *pad1 = new TPad("graph", "graph", 0, 0.3, 1, 1);
        pad1->Draw();
        pad1->cd();
        G_Electronic_Matching[i]->SetTitle(Form("Electronic Matching for %s", detectorName[i].c_str()));
        G_Electronic_Matching[i]->GetXaxis()->SetTitle("Pulser Value (V)");
        G_Electronic_Matching[i]->GetYaxis()->SetTitle("Pulser Value (V) (Reference)");
        G_Electronic_Matching[i]->SetMarkerStyle(20);
        G_Electronic_Matching[i]->SetMarkerSize(1);
        G_Electronic_Matching[i]->Draw("AP");
        c_match->cd();
        TPad *pad2 = new TPad("residuals", "residuals", 0, 0, 1, 0.3);
        pad2->Draw();
        pad2->cd();
        TGraphErrors *g_match_res = GetResiduals(G_Electronic_Matching[i], G_Electronic_Matching[i]->GetFunction("pol1"));
        g_match_res->SetTitle("Matching Residuals");
        g_match_res->GetXaxis()->SetTitle("Pulser Value (V)");
        g_match_res->GetYaxis()->SetTitle("Residuals (V)");
        g_match_res->SetMarkerStyle(20);
        g_match_res->SetMarkerSize(1);
        g_match_res->Draw("AP");
        c_match->Write();    

        TF1 *fout_func = (TF1*) G_Linearity[i]->GetFunction(func_name.c_str())->Clone(Form("Linearity_%s", detectorName[i].c_str()));
        // TF1 *fout_func = (TF1*) G_Electronic_Matching[i]->GetFunction("pol1")->Clone(Form("Linearity_%s", detectorName[i].c_str()));
        fout_func->SetTitle(Form("Linearity_%s", detectorName[i].c_str()));
        fout_func->Write();    

        if (i == det_ref && !ref_found) 
        {
            i = 10;
            ref_found = true;
        }
    }

    TCanvas *c2 = new TCanvas("Global", "Global", 800, 600);
    G_Global->SetTitle("Global_Chi2");
    G_Global->GetXaxis()->SetTitle("Detector Number");
    G_Global->GetYaxis()->SetTitle("Reduced Chi2");
    G_Global->SetMarkerStyle(20);
    G_Global->SetMarkerSize(1);
    G_Global->Draw("AP");
    c2->Write();

    TCanvas *c3 = new TCanvas("Delta", "Delta", 800, 600);
    G_Global_Delta->SetTitle("Global_Delta");
    G_Global_Delta->GetXaxis()->SetTitle("Detector Number");
    G_Global_Delta->GetYaxis()->SetTitle("Maximum Deviation (keV)");
    G_Global_Delta->SetMarkerStyle(20);
    G_Global_Delta->SetMarkerSize(1);
    G_Global_Delta->Draw("AP");
    c3->Write();

    TCanvas *c4 = new TCanvas("Global_FWHM", "Global_FWHM", 800, 600);
    TMultiGraph *mg = new TMultiGraph();
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    for (int j = 0; j < 10; j++)
    {
        if (G_Global_FWHM[j] == nullptr)
            continue;
        G_Global_FWHM[j]->SetTitle(Form("Global_FWHM_Peak%d", j));
        G_Global_FWHM[j]->GetXaxis()->SetTitle("Detector Number");
        G_Global_FWHM[j]->GetYaxis()->SetTitle("FWHM (keV)");
        G_Global_FWHM[j]->SetMarkerStyle(20);
        G_Global_FWHM[j]->SetMarkerSize(1);
        G_Global_FWHM[j]->SetMarkerColor(j+1);
        mg->Add(G_Global_FWHM[j]);
        legend->AddEntry(G_Global_FWHM[j], Form("~%.0f keV", Calibration[51]->Eval(j*10)), "p");

    }
    mg->GetXaxis()->SetTitle("Detector");
    mg->GetYaxis()->SetTitle("FWHM (keV)");
    mg->SetTitle("Electronic FWHM for each detector");
    mg->Draw("AP");
    legend->Draw("SAME");
    c4->Write();

    TCanvas *c5 = new TCanvas("Global_MeanSigma_Resolution", "Global_MeanSigma_Resolution", 800, 600);
    TGraphErrors *g_res = new TGraphErrors();
    for (int i = 0; i < G_Global_FWHM[3]->GetN(); i++)
    {
        int det = G_Global_FWHM[3]->GetX()[i];
        vector<double> values;
        for (int j = 0; j < 10; j++)
        {
            if (G_Global_FWHM[j] == nullptr)
                continue;
            
            double x, y;
            for (int p = 0; p < G_Global_FWHM[j]->GetN(); p++)
            {
                if (G_Global_FWHM[j]->GetX()[p] == det)
                {
                    G_Global_FWHM[j]->GetPoint(p, x, y);
                    values.push_back(y/2.355);
                    break;
                }
            }
        }
        double mean = accumulate(values.begin(), values.end(), 0.0) / values.size();
        double sigma = sqrt(accumulate(values.begin(), values.end(), 0.0, [mean](double acc, double val) { return acc + pow(val - mean, 2); }) / values.size());
        g_res->AddPoint(det, mean);
        g_res->SetPointError(g_res->GetN()-1, 0, sigma);
    }
    g_res->SetTitle("G_Electronic_Resolution_Sigma");
    g_res->SetName("G_Electronic_Resolution_Sigma");
    g_res->GetXaxis()->SetTitle("Detector");
    g_res->GetYaxis()->SetTitle("#sigma_{e} (keV)");
    g_res->SetMarkerStyle(20);
    g_res->SetMarkerSize(1);
    g_res->Draw("AP");
    c5->Write();

    fout->Close();   


    return 0;
}