#include "include/Detectors.hh"

void Test()
{
    TFile *f = new TFile("test.root", "RECREATE");

    double b12 = -10e3;
    double b31 = -5e3;  
    

    double a12 = 1.;
    double a31 = 1.;
    double a23 =1./(a31*a12);

    double b23 = - (a31*b12+b31 ) / (a31*a12) ;


    double alpha1 = 30;
    double alpha2 = 20;
    double alpha3 = 10;

    double beta1 = 10e3;
    double beta2 = 20e3;
    double beta3 = 30e3;

    double gamma1 = 0.001;
    double gamma2 = 0.002;
    double gamma3 = 0.003;

    default_random_engine generator;

    TH2D* H12 = new TH2D("H12", "H12", 1000, 0, 100e5, 1000, 0, 100e5);
    TH2D* H13 = new TH2D("H31", "H31", 1000, 0, 100e5, 1000, 0, 100e5);
    TH2D* H23 = new TH2D("H23", "H23", 1000, 0, 100e5, 1000, 0, 100e5);

    for (int i = -20e5; i < 120e5; i+= 1e4)
    {
        ProgressBar(i, 100e5, 0, 0, "Filling histogram");
        double x = i;
        double y = a12*x + b12;

        double sigma1 = gamma1*x+alpha1*sqrt(x)+beta1;
        double sigma2 = gamma2*x+alpha2*sqrt(y)+beta2;

        normal_distribution<double> distribution1(x, sigma1);
        normal_distribution<double> distribution2(y, sigma2);
        for (int j = 0; j < 10e3; j++)
        {
            double xx = distribution1(generator);
            double yy = distribution2(generator);

            H12->Fill(xx, yy);
        }
    }

    for (int i = i = -20e5; i < 120e5; i+= 1e4)
    {
        ProgressBar(i, 100e5, 0, 0, "Filling histogram");
        double x = i;
        double y = a31*x + b31;

        double sigma3 = gamma3*x+alpha3*sqrt(x)+beta3;
        double sigma1 = gamma1*x+alpha1*sqrt(y)+beta1;

        normal_distribution<double> distribution1(x, sigma3);
        normal_distribution<double> distribution2(y, sigma1);
        for (int j = 0; j < 10e3; j++)
        {
            double xx = distribution1(generator);
            double yy = distribution2(generator);

            H13->Fill(xx, yy);
        }
    }

    for (int i = i = -20e5; i < 120e5; i+= 1e4)
    {
        ProgressBar(i, 100e5, 0, 0, "Filling histogram");
        double x = i;
        double y = a23*x + b23;

        double sigma2 = gamma2*x+alpha2*sqrt(x)+beta2;
        double sigma3 = gamma3*x+alpha3*sqrt(y)+beta3;

        normal_distribution<double> distribution1(x, sigma2);
        normal_distribution<double> distribution2(y, sigma3);
        for (int j = 0; j < 10e3; j++)
        {
            double xx = distribution1(generator);
            double yy = distribution2(generator);

            H23->Fill(xx, yy);
        }
    }

    TCanvas* c = new TCanvas("c12", "c12", 1920, 1080);
    H12->Draw("COLZ");
    c->Write();

    //fit slices Y
    int counter = 0;
    TGraph* g12 = new TGraph();
    for (int bin_y = 0; bin_y < H12->GetNbinsY(); bin_y++)
    {
        TH1D* slice = H12->ProjectionX("slice", bin_y, bin_y);
        if (slice->GetEntries() == 0)
            continue;
        TF1* f = new TF1("f", "gaus");
        double y = H12->GetYaxis()->GetBinCenter(bin_y);
        double x = (y-b12)/a12;
        f->FixParameter(1, x);
        slice->Fit(f);
        delete slice;
        g12->SetPoint(counter, x, f->GetParameter(2));
        counter++;
        delete f;
    }
    
    TCanvas *c12 = new TCanvas("c12_res", "c12_res", 1920, 1080);
    g12->Draw("AP");
    c12->Write();

    TCanvas* c1 = new TCanvas("c31", "c31", 1920, 1080);
    H13->Draw("COLZ");
    c1->Write();

    //fit slices Y
    counter=0;
    TGraph* g31 = new TGraph();
    for (int bin_y = 0; bin_y < H13->GetNbinsY(); bin_y++)
    {
        TH1D* slice = H13->ProjectionX("slice", bin_y, bin_y);
        if (slice->GetEntries() == 0)
            continue;
        TF1* f = new TF1("f", "gaus");
        double y = H13->GetYaxis()->GetBinCenter(bin_y);
        double x = (y-b31)/a31;
        f->FixParameter(1, x);
        slice->Fit(f);
        delete slice;
        g31->SetPoint(counter, x, f->GetParameter(2));
        counter++;
        delete f;
    }

    TCanvas *c31 = new TCanvas("c31_res", "c31_res", 1920, 1080);
    g31->Draw("AP");
    c31->Write();

    TCanvas* c2 = new TCanvas("c23", "c23", 1920, 1080);
    H23->Draw("COLZ");
    c2->Write();  

    //fit slices Y
    counter=0;
    TGraph* g23 = new TGraph();
    for (int bin_y = 0; bin_y < H23->GetNbinsY(); bin_y++)
    {
        TH1D* slice = H23->ProjectionX("slice", bin_y, bin_y);
        if (slice->GetEntries() == 0)
            continue;
        TF1* f = new TF1("f", "gaus");
        double y = H23->GetYaxis()->GetBinCenter(bin_y);
        double x = (y-b23)/a23;
        f->FixParameter(1, x);
        slice->Fit(f);
        delete slice;
        g23->SetPoint(counter, x, f->GetParameter(2));
        counter++;
        delete f;

    }

    TCanvas *c23 = new TCanvas("c23_res", "c23_res", 1920, 1080);
    g23->Draw("AP");
    c23->Write();

    f->Close();


    


}