#include "InterStrip.hh"

int main()
{
    InitDetectors("../Grouper/Config_Files/sample.pid");
    InitWindows();
    TFile *final_f = MyTFile("InterStrip.root", "RECREATE");
    string path = DIR_DATA_HDD + "../SIMULATED_DATA/03-17/";
    // FILES.push_back(make_pair(250, MyTFile(path + "32Ar_inter025_a1_b0_analysed.root", "READ")));
    //FILES.push_back(make_pair(200.150, MyTFile(path + "32Ar_inter02000150_a1_b0_analysed.root", "READ")));
    // FILES.push_back(make_pair(150, MyTFile(path + "32Ar_inter015_a1_b0_analysed.root", "READ")));
    // FILES.push_back(make_pair(50, MyTFile(path + "32Ar_inter005_a1_b0_analysed3.root", "READ")));
    //FILES.push_back(make_pair(25, MyTFile(path + "32Ar_inter0025_a1_b0_analysed.root", "READ")));
    // FILES.push_back(make_pair(0, MyTFile(path + "32Ar_inter000_a1_b0_analysed.root", "READ")));
    FILES.push_back(make_pair(0, MyTFile(DIR_ROOT_DATA_SIMULATED + "03-09/32Ar_ION_a1_b0_analysed.root", "READ")));
    TFile *exp_f = MyTFile(DIR_ROOT_DATA_CALIBRATED + "Calibrated (copy).root", "READ");

    InitCalib();
    InitElectronicResolution();

    int counter = 0;
    for (auto &pair : FILES)
    {
        double interstrip = pair.first;
        TFile *file = pair.second;

        Info("File: " + to_string(interstrip));

        for (int det = 11; det < 86; det++)
        {
            if (IsDetectorSiliStrip(det))
            {
                Info("Detector: " + detectorName[det], 1);

                if (counter == 0)
                {
                    c[det] = new TCanvas(detectorName[det].c_str(), detectorName[det].c_str(), 800, 600);
                    leg[det] = new TLegend(0.7, 0.7, 0.9, 0.9);

                    H[10][det] = (TH1D *)exp_f->Get((detectorName[det] + "/32Ar/H_Exp_32Ar_" + detectorName[det]).c_str());
                    H[10][det]->SetLineColor(kBlack);
                    H[10][det]->Rebin(10);
                    H[10][det]->Draw("HIST");
                }

                H[counter][det] = (TH1D *)file->Get(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_All").c_str());
                H[counter][det]->SetLineColor(counter+2);
                H[counter][det]->GetXaxis()->SetTitle("Energy [keV]");
                H[counter][det]->GetYaxis()->SetTitle("Counts");
                H[counter][det]->GetXaxis()->CenterTitle();
                H[counter][det]->GetYaxis()->CenterTitle();

                leg[det]->AddEntry(H[counter][det], (to_string(interstrip) + " um").c_str(), "l");

                c[det]->cd();

                TH1D* copy = (TH1D*)H[counter][det]->Clone();
                H[counter][det]->Reset();
                double prob =1./10;
                for (int i = 1; i <= copy->GetEntries(); i++)
                {
                    double value = copy->GetRandom();
                    double offset = Calibration_Function[det]->GetParameter(1) * Detector_Resolution[det] + Calibration_Function[det]->GetParameter(2) * pow(Detector_Resolution[det], 2);
                    double res = gRandom->Gaus(value, sqrt(3.6 * 1e-3 * 0.134) * sqrt(value) + offset);

                    // cout << gRandom->Uniform() << endl;

                    // // random from 0 to 1
                    if (gRandom->Uniform() < 0.01 && det > 50)
                    {
                        res += abs(gRandom->Gaus(0, 25));
                    }

                    H[counter][det]->Fill(res);
                }

                H[10][det]->GetXaxis()->SetRangeUser(WindowsMap["32Ar"][14][det].first, WindowsMap["32Ar"][14][det].second);
                H[counter][det]->GetXaxis()->SetRangeUser(WindowsMap["32Ar"][14][det].first, WindowsMap["32Ar"][14][det].second);
                H[counter][det]->Scale(H[10][det]->Integral() / H[counter][det]->Integral());
                H[counter][det]->Rebin(10);

                
                H[counter][det]->Draw("HIST SAME");
                
                if (counter == FILES.size() - 1)
                {
                    final_f->cd();
                    leg[det]->Draw("SAME");
                    c[det]->Write();
                }
            }
        }

        counter++;
    }

    final_f->Close();
}