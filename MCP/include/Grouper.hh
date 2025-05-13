#include "../../Grouper/include/Detectors.hh"

string year;
string ROOT_filename;
string ROOT_basefilename;

TH1D *H_MCP_CORNERS_TIME[5];
TH1D *H_Channel[6];
TH2D *H_QDC1_QDC2;
TH2D *H_RAW_2D;
TH2D *H_LOG_2D;

TFile *OUTPUT_File;
double X_Tree, Y_Tree;

TTree *OutTree;

void InitHistograms()
{
    
    Info("Init Histograms");
    OUTPUT_File->cd();

    OutTree = new TTree("treeMCP", "treeMCP");
    OutTree->Branch("X", &X_Tree, "X/D");
    OutTree->Branch("Y", &Y_Tree, "Y/D");

    for (int i = 1; i <= 5; i++)
    {
        H_MCP_CORNERS_TIME[i] = new TH1D(("H_MCP_CORNERS_TIME_" + to_string(i)).c_str(), ("H_MCP_CORNERS_TIME_" + to_string(i)).c_str(), 500, -500, 500);
        H_MCP_CORNERS_TIME[i]->GetXaxis()->SetTitle("Time [ns]");
        H_MCP_CORNERS_TIME[i]->GetYaxis()->SetTitle("Counts");
        H_MCP_CORNERS_TIME[i]->GetXaxis()->CenterTitle();
        H_MCP_CORNERS_TIME[i]->GetYaxis()->CenterTitle();

        H_Channel[i] = new TH1D(("H_Channel_" + to_string(i)).c_str(), ("H_Channel_" + to_string(i)).c_str(), 10000, 0, 100000);
        H_Channel[i]->GetXaxis()->SetTitle("Channel");
        H_Channel[i]->GetYaxis()->SetTitle("Counts");
        H_Channel[i]->GetXaxis()->CenterTitle();
        H_Channel[i]->GetYaxis()->CenterTitle();
    }

    H_QDC1_QDC2 = new TH2D("H_QDC1_QDC2", "H_QDC1_QDC2", 10000, 0, 100000, 10000, 0, 100000);

    H_RAW_2D = new TH2D("H_RAW_2D", "H_RAW_2D", 1000, -1, 1, 1000, -1, 1);
    H_RAW_2D->GetXaxis()->SetTitle("X");
    H_RAW_2D->GetYaxis()->SetTitle("Y");
    H_RAW_2D->GetXaxis()->CenterTitle();
    H_RAW_2D->GetYaxis()->CenterTitle();

    H_LOG_2D = new TH2D("h_Image", "h_Image", 1000, -1, 1, 1000, -1, 1);
    H_LOG_2D->GetXaxis()->SetTitle("X");
    H_LOG_2D->GetYaxis()->SetTitle("Y");
    H_LOG_2D->GetXaxis()->CenterTitle();
    H_LOG_2D->GetYaxis()->CenterTitle();

    Info("Init Histograms Done");
}

void WriteHistograms()
{
    Info("Write Histograms");

    OUTPUT_File->cd();
    TCanvas *c = new TCanvas("TIME_CORNERS", "TIME_CORNERS", 800, 800);
    TLegend *leg = new TLegend(0.1, 0.7, 0.3, 0.9);
    for (int i = 1; i < 5; i++)
    {
        H_MCP_CORNERS_TIME[i]->Write();
        H_Channel[i]->Write();
        c->cd();
        H_MCP_CORNERS_TIME[i]->SetLineColor(i);
        H_MCP_CORNERS_TIME[i]->Draw("SAME");
        leg->AddEntry(H_MCP_CORNERS_TIME[i], ("Corner " + to_string(i)).c_str(), "l");
    }
    leg->Draw("SAME");
    c->Write();

    OUTPUT_File->cd();
    H_Channel[5]->Write();
    H_RAW_2D->Write();
    H_LOG_2D->Write();
    H_QDC1_QDC2->Write();

    OutTree->Write();

    OUTPUT_File->Close();
    Info("Write Histograms Done");
}

pair<double, double> ComputeXY(double A, double B, double C, double D)
{
    double norm = A + B + C + D;
    double X = (A + B - C - D) / (A + B + C + D);
    double Y = (A - B + C - D) / (A + B + C + D);
    return make_pair(X, Y);
}

pair<double, double> ComputeXY_log(double A, double B, double C, double D)
{
    double norm = A + B + C + D;
    double X = -(log(A/norm) + log(B/norm) - log(C/norm) - log(D/norm)) / (log(A/norm) + log(B/norm) + log(C/norm) + log(D/norm));
    double Y = -(log(A/norm) - log(B/norm) + log(C/norm) - log(D/norm)) / (log(A/norm) + log(B/norm) + log(C/norm) + log(D/norm));
    return make_pair(X, Y);
}

void ReadTree_MCP(TTreeReaderArray<Signal> &signals)
{
    double MCP_Time = 0;
    int counter_corner=0;
    vector<Signal> signals_corner(5);
    for (int i = 0; i < signals.GetSize(); i++)
    {
        Signal signal = signals[i];
        
        // IS MCP
        if (signal.Label == 5)
        {
            MCP_Time = signal.Time;
            signals_corner[0] = signal;
        }

        if (signal.Label >= 1 && signal.Label <= 4)
        {
            if (signals_corner[signal.Label].isValid)
                continue;

            counter_corner ++;
            signals_corner[signal.Label] = signal;
        }
    }

    if (abs(signals_corner[0].Channel-signals_corner[0].Pileup) > 3000)
    {
        return;
    }

    if (counter_corner == 4)
    {
        for (int i = 1; i < 5; i++)
        {

            if ((signals_corner[i].Time-MCP_Time) > -100 || (signals_corner[i].Time-MCP_Time) < -180)
            {
                return;
            }

            H_MCP_CORNERS_TIME[i]->Fill((signals_corner[i].Time-MCP_Time));
            H_Channel[signals_corner[i].Label]->Fill(signals_corner[i].Channel);

            
        }
        H_Channel[5]->Fill(signals_corner[0].Channel);
        H_QDC1_QDC2->Fill(signals_corner[0].Channel, signals_corner[0].Pileup);

        

        if (signals_corner[0].Channel <= 0 || signals_corner[1].Channel <= 0 || signals_corner[2].Channel <= 0 || signals_corner[3].Channel <= 0 || signals_corner[4].Channel <= 0)
        {
            return;
        }
        // cout << signals_corner[0].Channel << " " << signals_corner[1].Channel << " " << signals_corner[2].Channel << " " << signals_corner[3].Channel << " " << signals_corner[4].Channel << endl;

        pair<double, double> XY = ComputeXY(signals_corner[2].Channel, signals_corner[1].Channel, signals_corner[3].Channel, signals_corner[4].Channel);
        pair<double, double> XY_log = ComputeXY_log(signals_corner[2].Channel, signals_corner[1].Channel, signals_corner[3].Channel, signals_corner[4].Channel);

        double X = XY.first;
        double Y = XY.second;
        double X_log = XY_log.first;
        double Y_log = XY_log.second;

        H_RAW_2D->Fill(X, Y);
        H_LOG_2D->Fill(X_log, Y_log);

        X_Tree = X;
        Y_Tree = Y;
        OutTree->Fill();
    }

}