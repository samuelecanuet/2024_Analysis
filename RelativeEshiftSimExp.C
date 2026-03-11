#include "Grouper/include/Detectors.hh"
#include "Grouper/include/Utilities.hh"



int RelativeEshiftSimExp()
{
    FLAG2025 = true;
    InitDetectors("Grouper/Config_Files/sample.pid");




    

    TFile *Exp_FILE = MyTFile((DIR_ROOT_DATA_ANALYSED + "32Ar_" + to_string(YEAR) + "_analysed.root").c_str(), "READ");   
    TCanvas *cEshiftCompare_Corrected = (TCanvas *)Exp_FILE->Get("Eshift/Eshift_M3/Eshift_3");
    TGraphErrors* G_Exp_Eshift = nullptr;
    for (auto key : *cEshiftCompare_Corrected->GetListOfPrimitives())
    {
        if (string(key->ClassName()).find("TGraphErrors") != string::npos)
        {
            G_Exp_Eshift = (TGraphErrors *)key;
        }
    }
    if (G_Exp_Eshift == nullptr) Error("Eshift graph not found in experimental data");
    Info("Eshift graph loaded from experimental data", 1);

    double mean_exp_strip[SIGNAL_MAX] = {0};
    double mean_exp_strip_err[SIGNAL_MAX] = {0};

    double mean_exp = 0;
    double err_mean_exp = 0;
    for (int i = 0; i < G_Exp_Eshift->GetN(); i++)
    {

        double x = G_Exp_Eshift->GetX()[i];
        mean_exp += G_Exp_Eshift->GetY()[i];
        err_mean_exp += pow(G_Exp_Eshift->GetErrorY(i), 2);

        mean_exp_strip[GetDetectorChannel(x)] += G_Exp_Eshift->GetY()[i];
        mean_exp_strip_err[GetDetectorChannel(x)] += pow(G_Exp_Eshift->GetErrorY(i), 2);
    }
    mean_exp /= G_Exp_Eshift->GetN();
    err_mean_exp = sqrt(err_mean_exp) / G_Exp_Eshift->GetN();

    for (int strip = 0; strip < 20; strip++)
    {
        mean_exp_strip[strip] /= 8;
        mean_exp_strip_err[strip] = sqrt(mean_exp_strip_err[strip]) / 4;
    }


    //// SIM
    // a=1
    TFile *f1 = MyTFile("/run/media/local1/DATANEX/Samuel-G4/32Ar_IAS_a1.0_b0.0_result.root", "READ");
    TGraphErrors* G_Sim_Eshift1 = nullptr;
    TCanvas *c_Sim = (TCanvas *)f1->Get("Eshift");
    if (c_Sim == nullptr)
    {
        Warning("Eshift canvas not found in simulated file");
    }
    for (auto key : *c_Sim->GetListOfPrimitives())
    {
        if (string(key->ClassName()).find("TGraphErrors") != string::npos)
        {
            G_Sim_Eshift1 = (TGraphErrors *)key;
        }
    }

    double mean_sim1_strip[SIGNAL_MAX] = {0};
    double mean_sim1_strip_err[SIGNAL_MAX] = {0};
    double mean_sim1 = 0;
    double err_mean_sim1 = 0;
    for (int i = 0; i < G_Sim_Eshift1->GetN(); i++)
    {
        mean_sim1 += G_Sim_Eshift1->GetY()[i];
        err_mean_sim1 += pow(G_Sim_Eshift1->GetErrorY(i), 2);

        double x = G_Sim_Eshift1->GetX()[i];
        mean_sim1_strip[GetDetectorChannel(x)] += G_Sim_Eshift1->GetY()[i];
        mean_sim1_strip_err[GetDetectorChannel(x)] += pow(G_Sim_Eshift1->GetErrorY(i), 2);
    }
    mean_sim1 /= G_Sim_Eshift1->GetN();
    err_mean_sim1 = sqrt(err_mean_sim1) / G_Sim_Eshift1->GetN();

    for (int strip = 0; strip < 20; strip++)
    {
        mean_sim1_strip[strip] /= 8;
        mean_sim1_strip_err[strip] = sqrt(mean_sim1_strip_err[strip]) / 4;
    }

    TGraphErrors* G_Sim_Eshift1_Relative = new TGraphErrors();
    for (int i = 0; i < G_Sim_Eshift1->GetN(); i++)
    {
        double x = G_Sim_Eshift1->GetX()[i];
        double y = G_Sim_Eshift1->GetY()[i];
        double err_y = G_Sim_Eshift1->GetErrorY(i);
        // double exp_y = G_Exp_Eshift->GetY()[i];
        // double err_exp_y = G_Exp_Eshift->GetErrorY(i);  
        double relative_y = (y - (mean_sim1_strip[GetDetectorChannel(x)] - mean_exp_strip[GetDetectorChannel(x)]));
        double relative_err_y = sqrt( pow(err_y / mean_exp_strip[GetDetectorChannel(x)], 2) + pow(mean_sim1_strip_err[GetDetectorChannel(x)] * y / (mean_exp_strip[GetDetectorChannel(x)] * mean_exp_strip[GetDetectorChannel(x)]), 2) + pow(mean_exp_strip_err[GetDetectorChannel(x)] * y / (mean_exp_strip[GetDetectorChannel(x)] * mean_exp_strip[GetDetectorChannel(x)]), 2) + pow(err_mean_sim1 / mean_exp_strip[GetDetectorChannel(x)], 2) + pow(err_mean_exp * y / (mean_exp_strip[GetDetectorChannel(x)] * mean_exp_strip[GetDetectorChannel(x)]), 2) ) ;
        G_Sim_Eshift1_Relative->AddPoint(x, relative_y);
        G_Sim_Eshift1_Relative->SetPointError(i, 0, relative_err_y);
    }



    ///

    TFile *f07 = MyTFile("/run/media/local1/DATANEX/Samuel-G4/32Ar_IAS_a0.7_b0.0_result.root", "READ");
    TGraphErrors* G_Sim_Eshift07 = nullptr;
    c_Sim = (TCanvas *)f07->Get("Eshift");
    if (c_Sim == nullptr)
    {
        Error("Eshift canvas not found in simulated file");
    }
    for (auto key : *c_Sim->GetListOfPrimitives())
    {
        if (string(key->ClassName()).find("TGraphErrors") != string::npos)
        {
            G_Sim_Eshift07 = (TGraphErrors *)key;
        }
    }


    double mean_sim07_strip[SIGNAL_MAX] = {0};
    double mean_sim07_strip_err[SIGNAL_MAX] = {0};

    double mean_sim07 = 0;
    double err_mean_sim07 = 0;
    for (int i = 0; i < G_Sim_Eshift07->GetN(); i++)
    {
        mean_sim07 += G_Sim_Eshift07->GetY()[i];
        err_mean_sim07 += pow(G_Sim_Eshift07->GetErrorY(i), 2);

        double x = G_Sim_Eshift07->GetX()[i];
        mean_sim07_strip[GetDetectorChannel(x)] += G_Sim_Eshift07->GetY()[i];
        mean_sim07_strip_err[GetDetectorChannel(x)] += pow(G_Sim_Eshift07->GetErrorY(i), 2);
    }
    mean_sim07 /= G_Sim_Eshift07->GetN();
    err_mean_sim07 = sqrt(err_mean_sim07) / G_Sim_Eshift07->GetN();

    for (int strip = 0; strip < 20; strip++)
    {
        mean_sim07_strip[strip] /= 8;
        mean_sim07_strip_err[strip] = sqrt(mean_sim07_strip_err[strip]) / 4;
    }


    TGraphErrors *G_Sim_Eshift07_Relative = new TGraphErrors();
    for (int i = 0; i < G_Sim_Eshift07->GetN(); i++)
    {
        double x = G_Sim_Eshift07->GetX()[i];
        double y = G_Sim_Eshift07->GetY()[i];
        double err_y = G_Sim_Eshift07->GetErrorY(i);
        // double exp_y = G_Exp_Eshift->GetY()[i];
        // double err_exp_y = G_Exp_Eshift->GetErrorY(i);  
        double relative_y = (y - (mean_sim07_strip[GetDetectorChannel(x)] - mean_exp_strip[GetDetectorChannel(x)]));
        double relative_err_y = sqrt( pow(err_y / mean_exp_strip[GetDetectorChannel(x)], 2) + pow(mean_sim07_strip_err[GetDetectorChannel(x)] * y / (mean_exp_strip[GetDetectorChannel(x)] * mean_exp_strip[GetDetectorChannel(x)]), 2) + pow(mean_exp_strip_err[GetDetectorChannel(x)] * y / (mean_exp_strip[GetDetectorChannel(x)] * mean_exp_strip[GetDetectorChannel(x)]), 2) + pow(err_mean_sim07 / mean_exp_strip[GetDetectorChannel(x)], 2) + pow(err_mean_exp * y / (mean_exp_strip[GetDetectorChannel(x)] * mean_exp_strip[GetDetectorChannel(x)]), 2) ) ;
        G_Sim_Eshift07_Relative->AddPoint(x, relative_y);
        G_Sim_Eshift07_Relative->SetPointError(i, 0, relative_err_y);
    }

    TCanvas *cEshiftCompare = new TCanvas("cEshiftCompare", "cEshiftCompare", 1200, 800);
    G_Sim_Eshift1_Relative->SetName("a=1.0");
    G_Sim_Eshift1_Relative->SetMarkerStyle(20);
    G_Sim_Eshift1_Relative->SetMarkerColor(kRed);
    G_Sim_Eshift1_Relative->SetLineColor(kRed);
    G_Sim_Eshift1_Relative->SetMarkerSize(1);
    G_Sim_Eshift1_Relative->Draw("AP");
    G_Sim_Eshift07_Relative->SetName("a=0.7");
    G_Sim_Eshift07_Relative->SetMarkerStyle(20);
    G_Sim_Eshift07_Relative->SetMarkerColor(kBlue);
    G_Sim_Eshift07_Relative->SetLineColor(kBlue);
    G_Sim_Eshift07_Relative->SetMarkerSize(1);
    G_Sim_Eshift07_Relative->Draw("P SAME");
    cEshiftCompare->BuildLegend(0.6, 0.7, 0.9, 0.9);
    cEshiftCompare->Draw();


    TCanvas *cResidual = new TCanvas("cResidual", "cResidual", 1200, 800);
    TGraphErrors *G_Residual = new TGraphErrors();
    for (int i = 0; i < G_Sim_Eshift1_Relative->GetN(); i++)
    {
        double x = G_Sim_Eshift1_Relative->GetX()[i];
        double y_sim1 = G_Sim_Eshift1_Relative->GetY()[i];
        double err_y_sim1 = G_Sim_Eshift1_Relative->GetErrorY(i);
        double y_sim07 = G_Sim_Eshift07_Relative->GetY()[i];
        double err_y_sim07 = G_Sim_Eshift07_Relative->GetErrorY(i);
        double residual = (y_sim1 - y_sim07);
        double err_residual = sqrt( pow(err_y_sim1, 2) + pow(err_y_sim07, 2) );
        G_Residual->AddPoint(x, residual);
        G_Residual->SetPointError(i, 0, err_residual);
    }
    G_Residual->SetMarkerStyle(20);
    G_Residual->SetMarkerColor(kGreen+2);
    G_Residual->SetLineColor(kGreen+2);
    G_Residual->SetMarkerSize(1);
    G_Residual->Draw("AP");
    TLine *line = new TLine(G_Residual->GetXaxis()->GetXmin(), 0, G_Residual->GetXaxis()->GetXmax(), 0);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->Draw("SAME");
    // computing chi2 of residuals 
    double chi2_residual = 0;
    for (int i = 0; i < G_Residual->GetN(); i++)
    {
        double y = G_Residual->GetY()[i];
        double err_y = G_Residual->GetErrorY(i);
        if (err_y > 0)
            chi2_residual += pow(y / err_y, 2);
    }
    double ndf_residual = G_Residual->GetN();
    double p_value_residual = TMath::Prob(chi2_residual, ndf_residual);
    TLatex *latex = new TLatex(0.2, 0.8, Form("Chi2: %.2f, p-value: %.4f", chi2_residual/ndf_residual, p_value_residual));
    latex->SetNDC();
    latex->SetTextColor(kGreen+2);
    latex->Draw("SAME");

    cResidual->Draw();

    
    return 0;    
}