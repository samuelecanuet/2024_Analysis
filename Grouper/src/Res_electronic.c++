#include "Res_electronic.hh"



int main()
{
    file = new TFile("../../../../../../../../mnt/hgfs/shared-2/010_SiPMs_SiDet_Pulser_Po208_-15C_4T.root", "READ");
    TFile *newfile = new TFile("res.root", "RECREATE");
    InitDetectors("Config_Files/sample.pid");

    std::ofstream outfile("Config_Files/res_electronic.txt");

    TCanvas *c = new TCanvas("c", "c", 800, 600);
    TGraphErrors *g = new TGraphErrors();
    int counter = 0;

    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            int max = 62000;
            channel[i] = (TH1D *)file->Get((detectorName[i] + "_SPECTRO_GRAPH").c_str());
            channel[i]->GetXaxis()->SetRangeUser(45000, max);
            TF1 *gaus = new TF1("gaus", "gaus", 45000, max);
            gaus->SetNpx(100000);

            // gaus->SetParLimits(0, 100, 10000);
            // gaus->SetParameter(0, 1800);
            // gaus->SetParLimits(1, 45000, max);
            // // gaus->SetParameter(1, channel[i]->GetMean());
            gaus->SetParLimits(2, 10, 300);
            // gaus->SetParameter(2, 50);

            channel[i]->Fit(gaus, "RN", "", 45000, max);

            if (gaus->GetParError(2) > 1)
            {
                channel[i]->Fit(gaus, "RNL", "", 45000, max);
            }

            TCanvas *c = new TCanvas(detectorName[i].c_str(), detectorName[i].c_str(), 800, 600);
            channel[i]->Draw();
            gaus->Draw("same");
            c->Write();


            g->SetPoint(counter, i, gaus->GetParameter(2));
            g->SetPointError(counter, 0, gaus->GetParError(2));

            counter++;

            outfile << i << " " << gaus->GetParameter(2) << " " << gaus->GetParError(2) << endl;

            delete c;
            delete gaus;
        }
    }

    c->cd();
    g->GetYaxis()->SetRangeUser(0, 400);
    g->Draw("*AP");
    c->Write();

    newfile->Close();
    outfile.close();

    return 0;
}