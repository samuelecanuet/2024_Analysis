#include "../../Grouper/include/Detectors.hh"

string year;
string ROOT_filename;
string ROOT_basefilename;

TH1D *H_MCP_CORNERS_TIME[5];
TH1D* H_MCP_CORNERS_TIME_FULL[5];
TH1D *H_Channel[6];
TH2D *H_QDC1_QDC2;
TH2D *H_RAW_2D;
TH2D *H_RAW_2D_uncleaned;
TH2D *H_LOG_2D;

TFile *OUTPUT_File;
double X_Tree, Y_Tree;

TTree * Signals_Tree;
double C1, C2, C3, C4;
double X_s, Y_s;

TTree *OutTree;

void InitHistograms()
{
    
    Info("Init Histograms");
    OUTPUT_File->cd();

    OutTree = new TTree("treeMCP", "treeMCP");
    OutTree->Branch("X", &X_Tree, "X/D");
    OutTree->Branch("Y", &Y_Tree, "Y/D");

    Signals_Tree = new TTree("Signals_Tree", "Signals_Tree");
    Signals_Tree->Branch("X", &X_s, "X/D");
    Signals_Tree->Branch("Y", &Y_s, "Y/D");
    Signals_Tree->Branch("C1", &C1, "C1/D");
    Signals_Tree->Branch("C2", &C2, "C2/D");
    Signals_Tree->Branch("C3", &C3, "C3/D");
    Signals_Tree->Branch("C4", &C4, "C4/D");

    for (int i = 1; i <= 5; i++)
    {
        H_MCP_CORNERS_TIME[i] = new TH1D(("H_MCP_CORNERS_TIME_" + to_string(i)).c_str(), ("H_MCP_CORNERS_TIME_" + to_string(i)).c_str(), 500, -500, 500);
        H_MCP_CORNERS_TIME[i]->GetXaxis()->SetTitle("Time [ns]");
        H_MCP_CORNERS_TIME[i]->GetYaxis()->SetTitle("Counts / 2 ns");
        H_MCP_CORNERS_TIME[i]->GetXaxis()->CenterTitle();
        H_MCP_CORNERS_TIME[i]->GetYaxis()->CenterTitle();

        H_MCP_CORNERS_TIME_FULL[i] = new TH1D(("H_MCP_CORNERS_TIME_FULL_" + to_string(i)).c_str(), ("H_MCP_CORNERS_TIME_FULL_" + to_string(i)).c_str(), 1000, -1000, 1000);
        H_MCP_CORNERS_TIME_FULL[i]->GetXaxis()->SetTitle("Time [ns]");
        H_MCP_CORNERS_TIME_FULL[i]->GetYaxis()->SetTitle("Counts / 2 ns");
        H_MCP_CORNERS_TIME_FULL[i]->GetXaxis()->CenterTitle();
        H_MCP_CORNERS_TIME_FULL[i]->GetYaxis()->CenterTitle();

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

    H_RAW_2D_uncleaned = new TH2D("H_RAW_2D_uncleaned", "H_RAW_2D_uncleaned", 1000, -1, 1, 1000, -1, 1);
    H_RAW_2D_uncleaned->GetXaxis()->SetTitle("X");
    H_RAW_2D_uncleaned->GetYaxis()->SetTitle("Y");
    H_RAW_2D_uncleaned->GetXaxis()->CenterTitle();
    H_RAW_2D_uncleaned->GetYaxis()->CenterTitle();

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

    TCanvas *c_full = new TCanvas("TIME_CORNERS_FULL", "TIME_CORNERS_FULL", 800, 800);
    TLegend *leg_full = new TLegend(0.1, 0.7, 0.3, 0.9);
    for (int i = 1; i < 5; i++)
    {
        H_MCP_CORNERS_TIME_FULL[i]->Write();
        c_full->cd();
        H_MCP_CORNERS_TIME_FULL[i]->SetLineColor(i);
        H_MCP_CORNERS_TIME_FULL[i]->Draw("SAME");
        leg_full->AddEntry(H_MCP_CORNERS_TIME_FULL[i], ("Corner " + to_string(i)).c_str(), "l");
    }
    leg_full->Draw("SAME");
    c_full->Write();

    OUTPUT_File->cd();
    H_Channel[5]->Write();
    H_RAW_2D->Write();
    H_RAW_2D_uncleaned->Write();
    H_LOG_2D->Write();
    H_QDC1_QDC2->Write();

    OutTree->Write();

    OUTPUT_File->Close();
    Info("Write Histograms Done");
}

pair<double, double> ComputeXY(double A, double B, double C, double D)
{
    // double norm = A + B + C + D;
    // double X = (-A + B + C - D) / (A + B + C + D);
    // double Y = (-A - B + C + D) / (A + B + C + D);
    // return make_pair(X, Y);
    double norm = A + B + C + D;
    double X = (A + B - C - D) / (A + B + C + D);
    double Y = (A - B + C - D) / (A + B + C + D);
    return make_pair(X, Y);
}

pair<double, double> ComputeXY_log(double A, double B, double C, double D)
{
    // double norm = A + B + C + D;
    // double X = -(-log(A/norm) + log(B/norm) + log(C/norm) - log(D/norm)) / (log(A/norm) + log(B/norm) + log(C/norm) + log(D/norm));
    // double Y = -(-log(A/norm) - log(B/norm) + log(C/norm) + log(D/norm)) / (log(A/norm) + log(B/norm) + log(C/norm) + log(D/norm));
    // return make_pair(X, Y);
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

        if (signals_corner[0].Channel <= 0 || signals_corner[1].Channel <= 0 || signals_corner[2].Channel <= 0 || signals_corner[3].Channel <= 0 || signals_corner[4].Channel <= 0)
        {
            return;
        }

        // full time window
        for (int i = 1; i < 5; i++)
        {
            H_MCP_CORNERS_TIME_FULL[i]->Fill((signals_corner[i].Time-MCP_Time));
        }
        pair<double, double> XY_uncleaned = ComputeXY(signals_corner[2].Channel, signals_corner[1].Channel, signals_corner[3].Channel, signals_corner[4].Channel);
        H_RAW_2D_uncleaned->Fill(XY_uncleaned.first, XY_uncleaned.second);

        // Fill histograms for MCP corners
        for (int i = 1; i < 5; i++)
        {

            if ((signals_corner[i].Time-MCP_Time) > -100 || (signals_corner[i].Time-MCP_Time) < -190)
            {
                return;
            }

            H_MCP_CORNERS_TIME[i]->Fill((signals_corner[i].Time-MCP_Time));
            H_Channel[signals_corner[i].Label]->Fill(signals_corner[i].Channel);

            
        }
        H_Channel[5]->Fill(signals_corner[0].Channel);
        H_QDC1_QDC2->Fill(signals_corner[0].Channel, signals_corner[0].Pileup);

        

        
        // cout << signals_corner[0].Channel << " " << signals_corner[1].Channel << " " << signals_corner[2].Channel << " " << signals_corner[3].Channel << " " << signals_corner[4].Channel << endl;

        pair<double, double> XY = ComputeXY(signals_corner[2].Channel, signals_corner[1].Channel, signals_corner[3].Channel, signals_corner[4].Channel);
        pair<double, double> XY_log = ComputeXY_log(signals_corner[2].Channel, signals_corner[1].Channel, signals_corner[3].Channel, signals_corner[4].Channel);

        double X = XY.first;
        double Y = XY.second;
        double X_log = XY_log.first;
        double Y_log = XY_log.second;

        H_RAW_2D->Fill(X, Y);
        H_LOG_2D->Fill(X_log, Y_log);

        X_Tree = X_log;
        Y_Tree = Y_log;
        OutTree->Fill();

        X_s = X;
        Y_s = Y;
        C1 = signals_corner[2].Channel;
        C2 = signals_corner[1].Channel;
        C3 = signals_corner[3].Channel;
        C4 = signals_corner[4].Channel;
        Signals_Tree->Fill();
    }
}

struct Corner
{
    double X_real;
    double Y_real;
    double C1_real;
    double C2_real;
    double C3_real;
    double C4_real;
    double X_measured;
    double Y_measured;
    int C_counter = 0;
    double C1_measured;
    double C2_measured;
    double C3_measured;
    double C4_measured;
};

vector<Corner> M_corner[64];
double rho_real = 2.0;
double l_real = 1.4;
void InitSignalsValueforeachCorner()
{
    int n = 7;
    ifstream file2("guess_corner_2025_dis.txt");
    if (!file2.is_open())
    {
        Error("LoadCenters", "File guess_corner_2025.txt not found");
        return;
    }
    string line2;
    while (getline(file2, line2))
    {
        istringstream iss(line2);
        int cell;
        double ax, ay, bx, by, cx, cy, dx, dy;
        iss >> cell >> ax >> ay >> bx >> by >> cx >> cy >> dx >> dy;
        cell = cell - 1;

        double X_real = cell % n * rho_real - rho_real * (n-1) / 2 - l_real / 2;
        double Y_real = cell / n * rho_real - rho_real * (n-1) / 2 - l_real / 2;

        M_corner[cell] = vector<Corner>(4);

        // A
        M_corner[cell][0].X_measured = ax;
        M_corner[cell][0].Y_measured = ay;
        M_corner[cell][0].X_real = X_real;
        M_corner[cell][0].Y_real = Y_real;
        M_corner[cell][0].C1_real = 0.25 + (X_real + Y_real) / 4;
        M_corner[cell][0].C2_real = 0.25 + (-X_real + Y_real) / 4;
        M_corner[cell][0].C3_real = 0.25 + (X_real - Y_real) / 4;
        M_corner[cell][0].C4_real = 0.25 + (-X_real - Y_real) / 4;

        // B
        X_real += l_real / 2;
        M_corner[cell][1].X_measured = bx;
        M_corner[cell][1].Y_measured = by;
        M_corner[cell][1].X_real = X_real;
        M_corner[cell][1].Y_real = Y_real;
        M_corner[cell][1].C1_real = 0.25 + (X_real + Y_real) / 4;
        M_corner[cell][1].C2_real = 0.25 + (-X_real + Y_real) / 4;
        M_corner[cell][1].C3_real = 0.25 + (X_real - Y_real) / 4;
        M_corner[cell][1].C4_real = 0.25 + (-X_real - Y_real) / 4;

        // C
        Y_real += l_real / 2;
        M_corner[cell][2].X_measured = cx;
        M_corner[cell][2].Y_measured = cy;
        M_corner[cell][2].X_real = X_real;
        M_corner[cell][2].Y_real = Y_real;
        M_corner[cell][2].C1_real = 0.25 + (X_real + Y_real) / 4;   
        M_corner[cell][2].C2_real = 0.25 + (-X_real + Y_real) / 4;
        M_corner[cell][2].C3_real = 0.25 + (X_real - Y_real) / 4;
        M_corner[cell][2].C4_real = 0.25 + (-X_real - Y_real) / 4;

        // D
        X_real -= l_real / 2;
        M_corner[cell][3].X_measured = dx;
        M_corner[cell][3].Y_measured = dy;
        M_corner[cell][3].X_real = X_real;
        M_corner[cell][3].Y_real = Y_real;
        M_corner[cell][3].C1_real = 0.25 + (X_real + Y_real) / 4;
        M_corner[cell][3].C2_real = 0.25 + (-X_real + Y_real) / 4;
        M_corner[cell][3].C3_real = 0.25 + (X_real - Y_real) / 4;
        M_corner[cell][3].C4_real = 0.25 + (-X_real - Y_real) / 4;

        cout << cell << " " << ax << " " << ay << " " << bx << " " << by << " " << cx << " " << cy << " " << dx << " " << dy << endl;
    }

    TTreeReader *Reader = new TTreeReader(Signals_Tree);
    TTreeReaderValue<double> X_s(*Reader, "X");
    TTreeReaderValue<double> Y_s(*Reader, "Y");
    TTreeReaderValue<double> C1_s(*Reader, "C1");
    TTreeReaderValue<double> C2_s(*Reader, "C2");
    TTreeReaderValue<double> C3_s(*Reader, "C3");
    TTreeReaderValue<double> C4_s(*Reader, "C4");
    while (Reader->Next())
    {
        for (auto &hole : M_corner)
        {
            for (auto &corner : hole)
            {

                if (abs(corner.X_measured - *X_s) < 0.05 && abs(corner.Y_measured - *Y_s) < 0.05)
                {

                    // cout << "Found corner: " << corner.X_measured << " " << corner.Y_measured << " " << *X_s << " " << *Y_s << endl;
                    double sum = *C1_s + *C2_s + *C3_s + *C4_s;
                    corner.C1_measured += *C1_s / sum;
                    corner.C2_measured += *C2_s / sum;
                    corner.C3_measured += *C3_s / sum;
                    corner.C4_measured += *C4_s / sum;
                    corner.C_counter++;
                }
            }
        }
    }
}

double FunctionToMinimize()
{
    TGraph *G_C1 = new TGraph();
    TGraph *G_C2 = new TGraph();
    TGraph *G_C3 = new TGraph();
    TGraph *G_C4 = new TGraph();


    double chi2 = 0;
    for (int i = 0; i < 64; i++)
    {
        if (M_corner[i].size() < 4)
            continue;

        for (auto &corner : M_corner[i])
        {

            double X_real = corner.X_real;
            double Y_real = corner.Y_real;
            double C1_real = corner.C1_real;
            double C2_real = corner.C2_real;
            double C3_real = corner.C3_real;
            double C4_real = corner.C4_real;

            double X_measured = corner.X_measured;
            double Y_measured = corner.Y_measured;
            double C1_measured = corner.C1_measured / corner.C_counter;
            double C2_measured = corner.C2_measured / corner.C_counter;
            double C3_measured = corner.C3_measured / corner.C_counter;
            double C4_measured = corner.C4_measured / corner.C_counter;

            G_C1->AddPoint(C1_measured, C1_real);
            G_C2->AddPoint(C2_measured, C2_real);
            G_C3->AddPoint(C3_measured, C3_real);
            G_C4->AddPoint(C4_measured, C4_real);

            // Apply the parameters to the measured values
            // double C1_fit = par[0] + par[1] * C1_measured + par[2] * pow(C1_measured, 2) + par[3] * pow(C1_measured, 3);
            // double C2_fit = par[4] + par[5] * C2_measured + par[6] * pow(C2_measured, 2) + par[7] * pow(C2_measured, 3);
            // double C3_fit = par[8] + par[9] * C3_measured + par[10] * pow(C3_measured, 2) + par[11] * pow(C3_measured, 3);
            // double C4_fit = par[12] + par[13] * C4_measured + par[14] * pow(C4_measured, 2) + par[15] * pow(C4_measured, 3);

            // double sum = C1_fit + C2_fit + C3_fit + C4_fit;

            // double C1_fit_norm = C1_fit / sum;
            // double C2_fit_norm = C2_fit / sum;
            // double C3_fit_norm = C3_fit / sum;
            // double C4_fit_norm = C4_fit / sum;

            // // Compute the chi2
            // chi2 += pow((X_real - (C1_fit_norm + C2_fit_norm - C3_fit_norm - C4_fit_norm)), 2);
            // chi2 += pow((Y_real - (C1_fit_norm - C2_fit_norm + C3_fit_norm - C4_fit_norm)), 2);
        }
    }

    TCanvas *c = new TCanvas("C_Corners", "C_Corners", 800, 800);
    c->Divide(2, 2);
    c->cd(1);
    G_C1->SetTitle("C1");
    G_C1->GetXaxis()->SetTitle("C1 measured");
    G_C1->GetYaxis()->SetTitle("C1 real");
    G_C1->GetXaxis()->CenterTitle();
    G_C1->GetYaxis()->CenterTitle();
    G_C1->Draw("AP");
    c->cd(2);
    G_C2->SetTitle("C2");
    G_C2->GetXaxis()->SetTitle("C2 measured");
    G_C2->GetYaxis()->SetTitle("C2 real");
    G_C2->GetXaxis()->CenterTitle();
    G_C2->GetYaxis()->CenterTitle();
    G_C2->Draw("AP");
    c->cd(3);
    G_C3->SetTitle("C3");
    G_C3->GetXaxis()->SetTitle("C3 measured");
    G_C3->GetYaxis()->SetTitle("C3 real");
    G_C3->GetXaxis()->CenterTitle();
    G_C3->GetYaxis()->CenterTitle();
    G_C3->Draw("AP");
    c->cd(4);
    G_C4->SetTitle("C4");
    G_C4->GetXaxis()->SetTitle("C4 measured");
    G_C4->GetYaxis()->SetTitle("C4 real");
    G_C4->GetXaxis()->CenterTitle();
    G_C4->GetYaxis()->CenterTitle();
    G_C4->Draw("AP");
    c->Write();

    return chi2;
}

void FitTry()
{
    InitSignalsValueforeachCorner();


    FunctionToMinimize();

    // int Npar = 4;

    // Minimizer *minimizer = Factory::CreateMinimizer("Minuit2", "Migrad");
    // ROOT::Math::Functor functor(&FunctionToMinimize, 4*Npar);
    // minimizer->SetFunction(functor);
    // for (int i = 1; i <= 4; i++)
    // {
    //     for (int j = 0; j < Npar; j++)
    //     {
    //         minimizer->SetVariable(i*Npar+j, "par" + to_string(i*Npar+j), 0.0, 0.1);
    //     }
    // }
    // minimizer->SetMaxFunctionCalls(1000000);
    // minimizer->SetMaxIterations(1000000);

    // minimizer->Minimize();
    // const double *par = minimizer->X();

    // double chi2 = FunctionToMinimize(par);

    // cout << "Chi2: " << chi2 << endl;


}