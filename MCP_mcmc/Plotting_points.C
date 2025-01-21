
double Function(double *x, double *par)
{
    return 1 * (par[1] * x[0] + par[2] + pow(x[0], 2) + par[3] * pow(x[0], 3) + par[4] * x[1] * x[0] + par[5] * x[1] * x[1] * x[0] + par[6] * x[1] * x[1] * x[1] * x[0] + par[0] * pow(x[0], 4) + par[7] * pow(x[0], 5));
}

void Plotting_points()
{

    TFile *f = new TFile("res.root", "RECREATE");

    double rho = 2;
    TGraph *gr = new TGraph();

    // open txt file
    std::ifstream file("out_centers.txt");
    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        double x, y;
        int N;
        if (!(iss >> N >> x >> y))
        {
            break;
        }
        gr->SetPoint(N, x, y);
    }

    // TCanvas *c = new TCanvas("","",800,600);
    // gr->SetMarkerStyle(20);
    // gr->Draw("AP");
    // c->Write();

    TGraph2D *grx = new TGraph2D();
    TGraph2D *gry = new TGraph2D();

    TGraph *gx = new TGraph();
    TGraph *gy = new TGraph();
    int counter = -1;
    for (int i = 0; i < gr->GetN(); i++)
    {

        double x_center = i % 8 * rho - rho * (8 - 1) / 2;
        double y_center = i / 8 * rho - rho * (8 - 1) / 2;

        if (gr->GetY()[i] == 0 && gr->GetX()[i] == 0)
        {
            continue;
        }

        i+=1;
        if (i != 18 && i != 47 && i != 26 && i != 34 && i != 42 && i != 51 && i != 52 && i != 53 && i != 54 && i != 47 && i != 39 && i != 31 && i != 23 && i != 11 && i != 12 && i != 13 && i != 14 && i != 19 && i != 20 && i != 21 && i != 22 && i != 23 && i != 27 && i != 28 && i != 29 && i != 30 && i != 31 && i != 35 && i != 36 && i != 37 && i != 38 && i != 39 && i != 43 && i != 44 && i != 45 && i != 46)
        {
            continue;
        }
        i-=1;
        counter++;

        double x = gr->GetX()[i];
        double y = gr->GetY()[i];
        double r = sqrt(pow(gr->GetX()[i], 2) + pow(gr->GetY()[i], 2));
        cout << "r: " << r << "  " << i << endl;
        double phi = atan2(gr->GetY()[i], gr->GetX()[i]);
        double R = sqrt(pow(x_center, 2) + pow(y_center, 2));
        double Phi = atan2(y_center, x_center);

        gx->SetPoint(counter, x, x_center);
        grx->SetPoint(counter, x, y, x_center);
        gy->SetPoint(counter, y, y_center);
        gry->SetPoint(counter, y, x, y_center);
    }

    TCanvas *c0 = new TCanvas("", "", 800, 600);
    gx->SetMarkerStyle(20);
    gx->GetXaxis()->SetTitle("x");
    gx->GetYaxis()->SetTitle("X (mm)");
    gx->Draw("AP");
    // gx->Fit(Function, "R");
    // f0->Draw("SAME");
    c0->Write();

    // TCanvas *c1y = new TCanvas("","",800,600);
    // gy->SetMarkerStyle(20);
    // gy->GetXaxis()->SetTitle("y");
    // gy->GetYaxis()->SetTitle("X (mm)");
    // counter = 0;
    // for (int i = 0; i < grx->GetN(); i++)
    // {
    //     gry->SetPoint(i, grx->GetX()[i], grx->GetY()[i], grx->GetZ()[i]-f0->Eval(grx->GetX()[i]));
    //     gy -> SetPoint(i, grx->GetX()[i], grx->GetZ()[i]-f0->Eval(grx->GetX()[i]));
    // }
    // gry->SetMarkerStyle(20);
    // gry->Draw("AP");
    // c1y->Write();

    TCanvas *c1 = new TCanvas("grx", "", 800, 600);
    grx->SetMarkerStyle(20);
    grx->GetXaxis()->SetTitle("x");
    grx->GetYaxis()->SetTitle("y");
    grx->GetZaxis()->SetTitle("X (mm)");

    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
    // fitting by a 2D polynomial function of degree of 5
    Double_t *x_values = grx->GetX();
    Int_t n_points = grx->GetN();
    Double_t min_x = TMath::MinElement(n_points, x_values);
    Double_t max_x = TMath::MaxElement(n_points, x_values);

    Double_t *y_values = grx->GetY();
    Double_t min_y = TMath::MinElement(n_points, y_values);
    Double_t max_y = TMath::MaxElement(n_points, y_values);

    TF2 *f1 = new TF2("f1", Function, min_x, max_x, min_y, max_y, 8);
    f1->SetParameters(0.58, 10.6, -1.4, 13);
    /// f1->SetParameters(13, -1.4, 10.6, 0.58);
    TFitResultPtr rx = grx->Fit(f1, "MULTITHREAD S");
    f1->Draw("SURF2");
    grx->Draw("AP SAME");
    c1->Write();

    // compute residuals and plot in a hist2d
    TCanvas *c2 = new TCanvas("", "", 800, 600);
    TH2D *h = new TH2D("h", "h", 50, min_x, max_x, 50, min_y, max_y);
    for (int i = 0; i < grx->GetN(); i++)
    {
        double x = grx->GetX()[i];
        double y = grx->GetY()[i];
        double z = grx->GetZ()[i];
        double res = abs(z - f1->Eval(x, y));

        h->SetBinContent(h->FindBin(x, y), res);
    }

    h->Draw("LEGO2");
    c2->Write();

    TCanvas *c3 = new TCanvas("", "", 800, 600);
    gry->SetMarkerStyle(20);
    gry->GetXaxis()->SetTitle("x");
    gry->GetYaxis()->SetTitle("y");
    gry->GetZaxis()->SetTitle("Y (mm)");

    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
    // fitting by a 2D polynomial function of degree of 5
    Double_t *x_values_y = gry->GetX();
    Int_t n_points_y = gry->GetN();
    Double_t min_x_y = TMath::MinElement(n_points_y, x_values_y);
    Double_t max_x_y = TMath::MaxElement(n_points_y, x_values_y);

    Double_t *y_values_y = gry->GetY();
    Double_t min_y_y = TMath::MinElement(n_points_y, y_values_y);
    Double_t max_y_y = TMath::MaxElement(n_points_y, y_values_y);

    TF2 *f1_y = new TF2("f1_y", Function, min_x_y, max_x_y, min_y_y, max_y_y, 8);
    f1_y->SetParameters(0.58, 10.6, -1.4, 13);
    /// f1->SetParameters(13, -1.4, 10.6, 0.58);
    TFitResultPtr ry = gry->Fit(f1_y, "MULTITHREAD S");
    f1_y->Draw("SURF2");
    gry->Draw("AP SAME");
    c3->Write();

    // compute residuals and plot in a hist2d
    TCanvas *c4 = new TCanvas("", "", 800, 600);
    TH2D *h_y = new TH2D("h_y", "h_y", 50, min_x_y, max_x_y, 50, min_y_y, max_y_y);
    for (int i = 0; i < gry->GetN(); i++)
    {
        double x = gry->GetX()[i];
        double y = gry->GetY()[i];
        double z = gry->GetZ()[i];
        double res = abs(z - f1_y->Eval(x, y));

        h_y->SetBinContent(h_y->FindBin(x, y), res);
    }

    h_y->Draw("LEGO2");
    c4->Write();

    TCanvas *c5 = new TCanvas("", "", 800, 600);
    auto *gr_real = new TGraph();
    auto *gr_fit = new TGraph();
    for (int i = 0; i < grx->GetN(); i++)
    {
        double x = grx->GetX()[i];
        double y = gry->GetX()[i];
        // double r = x;
        // double phi = y;
        double X = grx->GetZ()[i];
        double Y = gry->GetZ()[i];
        // double R = X;
        // double Phi = Y;

        gr_real->SetPoint(i, X, Y);
        gr_fit->SetPoint(i, f1->Eval(x, y), f1_y->Eval(y, x));
    }
    gr_real->SetMarkerStyle(20);
    gr_real->SetMarkerColor(kRed);
    gr_fit->SetMarkerStyle(20);
    gr_fit->SetMarkerColor(kBlue);
    gr_real->Draw("AP");
    gr_fit->Draw("P SAME");
    c5->Write();

    TCanvas *c6 = new TCanvas("", "", 800, 600);
    TH2D *h_r = new TH2D("h_r", "h_r", 50, -8, 8, 50, -8, 8);
    for (int i = 0; i < gr_real->GetN(); i++)
    {
        double x_real = gr_real->GetX()[i];
        double y_real = gr_real->GetY()[i];

        double x_fit = gr_fit->GetX()[i];
        double y_fit = gr_fit->GetY()[i];

        double diff_r = sqrt(pow(x_real - x_fit, 2) + pow(y_real - y_fit, 2));
        h_r->SetBinContent(h_r->FindBin(x_real, y_real), diff_r);
    }

    h_r->Draw("LEGO2");

    c6->Write();

}