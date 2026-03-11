#include "include/Detectors.hh"
#include "include/Utilities.hh"

// working with 2024 on run 115
// - checking 2024 all runs
// - checkiong 2025
// - find check for 1 or 2 lines
// - excluding points at edges
// (- fitting at the same time both line)

double startt = 10e3;
double endd = 170e3;
TGraph *DATA;
TH2D *H_SiPM_HighLow;
int Fitting2Lines()
{

    FLAG2025 = true;
    InitDetectors("Config_Files/sample.pid");
    InitRuns();

    double guess_a, guess_b;

    TFile *f = MyTFile((DIR_ROOT_DATA_MATCHED + "SiPM_Matching_test.root").c_str(), "READ");

    // TCanvas *cRaw = new TCanvas ("cRaw", "cRaw", 800, 600);
    TCanvas *c = (TCanvas *)f->Get("105/Graph_HighLowSiPM_Fitting_105");
    // lopping in TCnVAS 
    for (auto key : *c->GetListOfPrimitives())
    {
        TPad *pad = (TPad *)key;
        string name = pad->GetName();
        cout << name << endl;
        if (name.find("105_1") != string::npos)
        {
            for (auto key2 : *pad->GetListOfPrimitives())
            {
                if (string(key2->ClassName()) == "TF1")
                {
                    TF1 * fkey2 = (TF1 *)key2;
                    guess_a = fkey2->GetParameter(0);
                    guess_b = fkey2->GetParameter(1);
                }

                // if type  == TGRaphErrors
                if (string(key2->ClassName()) != "TGraph")
                    continue;
                TGraph * graph = (TGraph *)key2;
                // cRaw->cd();
                // graph->Draw("AP");
                DATA = (TGraph*)graph->Clone(); 
            }
            break;
        }
    }

    // fill hist 
    H_SiPM_HighLow = new TH2D("H_SiPM_HighLow", "H_SiPM_HighLow", 10000, startt, endd, 10000, startt*10, endd*10);
    for (int i = 0; i < DATA->GetN(); i++)
    {
        double x, y;
        DATA->GetPoint(i, x, y);
        H_SiPM_HighLow->Fill(x, y);
    }
    H_SiPM_HighLow->Rebin2D(4);

    TGraphErrors *G1 = new TGraphErrors();
    TGraphErrors *G2 = new TGraphErrors();
    // fitting slices or TSpectrum
    int bin_width = 25;
    for (int binx = 1; binx <= H_SiPM_HighLow->GetNbinsX(); binx+= bin_width)
    {
        TH1D *proj = H_SiPM_HighLow->ProjectionY(Form("proj_%d", binx), binx, binx+bin_width -1);
        if (proj->GetEntries() < 50)
            continue;

        // proj->Rebin(2);
        TF1 *gaus = new TF1("double_gauss", "[0]/sqrt(2*3.14)/[2] * exp(-0.5*((x - [1])/[2])**2) + [3]/sqrt(2*3.14)/[5] * exp(-0.5*((x -( [1]+[4] ))/[5])**2)", startt*10, endd*10);
        // find position mean of the 2 peaks with TSpectrum
        TSpectrum *s = new TSpectrum(2);
        int nfound = s->Search(proj, 20, "", 0.1);
        if (nfound != 2)
            continue;
        double *peaks = s->GetPositionX();
        if (peaks[1] < peaks[0])
            swap(peaks[0], peaks[1]);
        
        gaus->SetParLimits(0, 0, 10e4);
        gaus->SetParameter(1, peaks[0]);
        gaus->SetParLimits(1, peaks[0]-20e3, peaks[0]+20e3);
        gaus->SetParameter(2, 10e3);
        gaus->SetParLimits(2, 0, 10e4);
        gaus->SetParLimits(3, 0, 10e4);
        gaus->SetParameter(4, abs(peaks[0]-peaks[1]));
        gaus->SetParLimits(4, 10e3, 3*abs(peaks[0]-peaks[1]));
        gaus->SetParameter(5, 10e3);
        gaus->SetParLimits(5, 0, 10e4);

        TFitResultPtr r = proj->Fit(gaus, "SQNR");
        if (r != 0)
            continue;

        // rejecting if error overlapping points 
        if ((gaus->GetParError(1)) > abs(gaus->GetParameter(4)))
            continue;
        
        if ((gaus->GetParError(4)) > abs(gaus->GetParameter(4)))
            continue;

        if (abs(10e3 - gaus->GetParameter(4)) < 10)
            continue;
        
        G1->AddPoint(H_SiPM_HighLow->GetXaxis()->GetBinCenter(binx), gaus->GetParameter(1) + gaus->GetParameter(4));
        G1->SetPointError(G1->GetN() -1, 0, sqrt(pow(gaus->GetParError(1),2) + pow(gaus->GetParError(4),2)) );

        G2->AddPoint(H_SiPM_HighLow->GetXaxis()->GetBinCenter(binx), gaus->GetParameter(1));
        G2->SetPointError(G2->GetN() -1, 0, gaus->GetParError(1)); 
    }


    // delta condition to rejecting points
    TCanvas *cFitting = new TCanvas ("cFitting", "cFitting", 800, 600);
    TH1D *h1 = new TH1D("h1", "h1", 100, 0, 100e3);
    for (int i = 0; i < G1->GetN(); i++)
    {
        double x1, y1, x2, y2;
        G1->GetPoint(i, x1, y1);
        G2->GetPoint(i, x2, y2);
        if (y1 - y2 <= 0)
            continue;
        h1->Fill(y1 - y2);
    }

    TF1 *gaus = new TF1("gaus", "[0]/sqrt(2*3.14)/[2] * exp(-0.5*((x - [1])/[2])**2)+[3]", 0, h1->GetXaxis()->GetXmax());
    gaus->SetParameter(0, h1->Integral()* h1->GetBinWidth(1));
    gaus->SetParameter(1, h1->GetBinCenter(h1->GetMaximumBin()));
    gaus->SetParLimits(1, h1->GetBinCenter(h1->GetMaximumBin()) - 10e3, h1->GetBinCenter(h1->GetMaximumBin()) + 10e3);
    gaus->SetParameter(2, h1->GetRMS());
    gaus->SetParLimits(2, 1, 10e4);
    gaus->SetParLimits(3, 0, h1->GetMaximum());
    TFitResultPtr r = h1->Fit(gaus, "SR");
    h1->Draw("HIST");
    cFitting->Draw();

    double delta;
    double sigma_delta;
    if (r != 0)
    {
        delta = h1->GetBinCenter(h1->GetMaximumBin());
        sigma_delta = h1->GetRMS();
        Warning("Gaussian fit failed, using histogram peak and RMS for delta and sigma_delta", 1);
    }
    else
    {
        delta = gaus->GetParameter(1);
        sigma_delta = 3*gaus->GetParameter(2);
    }

    cout << "Delta: " << delta << " Sigma Delta: " << sigma_delta << endl;


    // cleaning points - rejection
    TGraphErrors *G1_clean = new TGraphErrors();
    TGraphErrors *G1_rejected = new TGraphErrors();

    TGraphErrors *G2_clean = new TGraphErrors();
    TGraphErrors *G2_rejected = new TGraphErrors();

    for (int i = 0; i < G1->GetN(); i++)
    {
        double x1, y1, x2, y2;
        G1->GetPoint(i, x1, y1);
        G2->GetPoint(i, x2, y2);
        double y1_err = G1->GetErrorY(i);
        double y2_err = G2->GetErrorY(i);
        if (abs((y1 - y2) - delta) < sigma_delta)
        {
            G1_clean->AddPoint(x1, y1);
            G1_clean->SetPointError(G1_clean->GetN() -1, 0, y1_err);
            G2_clean->AddPoint(x2, y2);
            G2_clean->SetPointError(G2_clean->GetN() -1, 0, y2_err);
        }
        else
        {
            G1_rejected->AddPoint(x1, y1);
            G2_rejected->AddPoint(x2, y2);
        }
    }
    TCanvas *cFitClean = new TCanvas ("cFitClean", "cFitClean", 800, 600);
    TPad *pad = new TPad("pad", "pad", 0, 0.3, 1, 1);
    pad->Draw();
    pad->cd();
    G1_clean->SetMarkerColor(kRed);
    G1_clean->SetMarkerStyle(20);
    TF1 *f1 = new TF1("pol1", "[1]*x + [0] + [2]*x*x", startt, endd-30e3);
    // f1->FixParameter(0, 10.8037);
    // f1->FixParameter(1, 0);
    G1_clean->Fit(f1, "QR", "", startt, endd-30e3);
    cout << "Fit 1: " << f1->GetParameter(0) << " +/- " << f1->GetParError(0) << ", " << f1->GetParameter(1) << " +/- " << f1->GetParError(1) << endl;
    f1->SetLineColor(kBlack);
    G2_clean->SetMarkerColor(kBlue);
    G2_clean->SetMarkerStyle(20);
    TF1 *f2 = new TF1("pol1", "[0]*x + [1] + [2]*x*x", startt, endd-30e3);
    f2->FixParameter(0, f1->GetParameter(0));
    G2_clean->Fit(f2, "QR", "", startt, endd-30e3);
    cout << "Fit 2: " << f2->GetParameter(0) << " +/- " << f2->GetParError(0) << ", " << f2->GetParameter(1) << " +/- " << f2->GetParError(1) << endl;
    f2->SetLineColor(kBlack);
    G1_rejected->SetMarkerColor(kGray);
    G1_rejected->SetMarkerStyle(20);
    G2_rejected->SetMarkerColor(kGray);
    G2_rejected->SetMarkerStyle(20);
    TMultiGraph *mg_clean = new TMultiGraph();
    mg_clean->Add(G1_clean, "AP");
    mg_clean->Add(G2_clean, "AP");
    mg_clean->Add(G1_rejected, "AP");
    mg_clean->Add(G2_rejected, "AP");
    mg_clean->Draw("AP");

    cFitClean->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
    pad2->Draw();
    pad2->cd();
    TGraphErrors *G1_clean_residuals = new TGraphErrors();
    TGraphErrors *G2_clean_residuals = new TGraphErrors();
        for (int i = 0; i < G1_clean->GetN(); i++)
    {
        double x1, y1, x2, y2;
        G1_clean->GetPoint(i, x1, y1);
        G2_clean->GetPoint(i, x2, y2);
        double resid1 = y1 - f1->Eval(x1);
        double resid2 = y2 - f2->Eval(x2);
        G1_clean_residuals->AddPoint(x1, resid1/f1->Eval(x1)*100);
        G1_clean_residuals->SetPointError(G1_clean_residuals->GetN() -1, 0, G1_clean->GetErrorY(i)/f1->Eval(x1)*100);
        G2_clean_residuals->AddPoint(x2, resid2/f2->Eval(x2)*100);
        G2_clean_residuals->SetPointError(G2_clean_residuals->GetN() -1, 0, G2_clean->GetErrorY(i)/f2->Eval(x2)*100);
    }
    G1_clean_residuals->SetMarkerColor(kRed);
    G1_clean_residuals->SetMarkerStyle(20);
    G1_clean_residuals->GetXaxis()->SetTitle("Channel");
    G1_clean_residuals->GetYaxis()->SetTitle("Residuals (\%)");
    G2_clean_residuals->SetMarkerColor(kBlue);
    G2_clean_residuals->SetMarkerStyle(20);
    G2_clean_residuals->GetXaxis()->SetTitle("Channel");
    G2_clean_residuals->GetYaxis()->SetTitle("Residuals (\%)");
    TMultiGraph *mg_clean_residuals = new TMultiGraph();
    mg_clean_residuals->Add(G1_clean_residuals, "AP");
    mg_clean_residuals->Add(G2_clean_residuals, "AP");
    mg_clean_residuals->Draw("AP");
    TLine *line_diag = new TLine(startt, 0, endd, 0);
    line_diag->SetLineColor(kGray);
    line_diag->SetLineStyle(2);
    line_diag->Draw("SAME");
    cFitClean->Draw();


    // plotting results
    TCanvas *cResults = new TCanvas ("cResults", "cResults", 800, 600);
    H_SiPM_HighLow->Draw("COLZ");
    f1->SetLineColor(kRed);
    f1->Draw("SAME");
    TLatex *lat = new TLatex();
    lat->SetNDC();
    lat->SetTextSize(0.04);
    lat->DrawLatex(0.15, 0.85, Form("%.6f x + %.2f", f1->GetParameter(0), f1->GetParError(0)));
    lat->Draw("SAME");
    f2->SetLineColor(kRed);
    f2->Draw("SAME");
    TLatex *lat2 = new TLatex();
    lat2->SetNDC();
    lat2->DrawLatex(0.15, 0.80, Form("%.6f x + %.2f", f2->GetParameter(0), f2->GetParError(0)));
    lat2->Draw("SAME");
    cResults->Draw();

    // plotting projected histogram along the fitetd fucntion f1 on the H2D
    TCanvas *cProjected = new TCanvas ("cProjected", "cProjected", 800, 600);
    TH2D *H_rot_f1 = new TH2D("H_rot_f1", "H_rot_f1", 2000, -100e3, 100e3, 200, startt, endd*10);
    TF1 *invertedf1 = InvertFunction(f1);
    for (int i = 0; i<DATA->GetN(); i++)
    {
        double x, y;
        DATA->GetPoint(i, x, y);
        double x_fit = invertedf1->Eval(y);
        H_rot_f1->Fill(x-x_fit, y);
    }
    H_rot_f1->Draw("COLZ");
    cProjected->Draw();


    
    double a = f1->GetParameter(0);
    double b = f1->GetParameter(1);
    


    double best_max = 0;
    double best_a = a;
    TH1D *H_rot_f1_final = new TH1D("H_rot_f1_final", "H_rot_f1_final", 200, -20e3, 20e3);
    // loop 
    double delta_a = 0.1;
    double step_a = delta_a / 20;
    for (double a_try = a - delta_a; a_try <= a + delta_a; a_try += step_a)
    {
        TF1 *invertedf1_try = new TF1("bilinf1_try", "[0]*x + [1]", 0, 100e7);
        invertedf1_try->SetParameter(0, 1/a_try);
        invertedf1_try->SetParameter(1, -b/a_try);

        TH1D *H_rot_f1_try = new TH1D("H_rot_f1_try", "H_rot_f1_try", 200, -1e3, 20e3);
        for (int i = 0; i<DATA->GetN(); i++)
        {
            double x, y;
            DATA->GetPoint(i, x, y);
            double x_fit = invertedf1_try->Eval(y);
            H_rot_f1_try->Fill(x - x_fit);
        }

        int max_bin = H_rot_f1_try->GetMaximumBin();
        double sum = H_rot_f1_try->GetBinContent(max_bin) + H_rot_f1_try->GetBinContent(max_bin -1) + H_rot_f1_try->GetBinContent(max_bin +1);

        if (best_max == 0 || sum > best_max)
        {
            best_max = sum;
            a = a_try;
            H_rot_f1_final = (TH1D *)H_rot_f1_try->Clone("H_rot_f1_final");
        }
    }

    TCanvas *cFinal = new TCanvas ("cFinal", "cFinal", 800, 600);
    H_rot_f1_final->Draw("HIST");
    TLatex *lat_final = new TLatex();
    lat_final->SetNDC();
    lat_final->SetTextSize(0.04);
    lat_final->DrawLatex(0.15, 0.85, Form("Best a: %.6f", a));
    lat_final->Draw("SAME");

    cFinal->Draw();


    
    

    return 0;
}