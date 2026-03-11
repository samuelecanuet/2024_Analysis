

#include "../../Grouper/include/Utilities.hh"

int RateCorrectionForPlot()
{
    Info("Rate Correction For Plot");

    TFile *f = MyTFile("../run_034_MCP_32Ar_BeamScan_2T_calibrated_secondstep.root", "READ");

    TTree *tree = (TTree *)f->Get("treeMCP");
    TTreeReader reader(tree);
    TTreeReaderValue<double> X_Tree(reader, "X");
    TTreeReaderValue<double> Y_Tree(reader, "Y");
    TTreeReaderValue<double> Time_Tree(reader, "Time");

    TH1D *H = new TH1D("H", "H", 10000, 0, 1000); // in s
    TH3D *H2 = new TH3D("H", "H", 100, 0, 1000, 80, -8, 8, 80, -8, 8); // in s

    TH2D *H2D = new TH2D("H2D", "H2D", 800, -8, 8, 800, -8, 8);


    while (reader.Next())
    {
        double x = *X_Tree;
        double y = *Y_Tree;
        double time = *Time_Tree / 1e9; // convert to seconds

        H->Fill(time);
        // H2D->Fill(x, y);
        H2->Fill(time, x, y);
    }

    // H->Draw("HIST");

    reader.Restart();
    while (reader.Next())
    {
        double x = *X_Tree;
        double y = *Y_Tree;
        double time = *Time_Tree / 1e9; // convert to seconds

        // double rate = H->GetBinContent(H->FindBin(time));
        double rate = H2->GetBinContent(H2->FindBin(time, x, y));
        if (rate > 0)
        {
            double weight = 1.0 / rate;
            H2D->Fill(x, y, weight);
        }
    }

    TFile *fout = MyTFile("RateCorrectedPlot.root", "RECREATE");
    fout->cd();

    H2D->SetName("H_reconstruction_interpolated_second_measurement");
    H2D->Write();



    return 0;
}