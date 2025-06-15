#include "Analyser.hh"

int main(int argc, char **argv)
{
    
    NUCLEUS = "32Ar";
    IAS[NUCLEUS] = 14;
    FLAG2025 = true;
    VERBOSE = 0;


    InitDetectors("Config_Files/sample.pid");
    
    
    ///////// READING FILE //////////
    TTreeReader *Reader;
    TTreeReaderArray<Signal> *Silicon;
    TTreeReaderValue<vector<vector<pair<Signal, Signal>>>> *SiPM_Groups;

    MERGED_File = MyTFile((DIR_ROOT_DATA_MERGED + NUCLEUS + "_merged.root").c_str(), "READ");
   
    
    ///////// NEW FILE //////////
    ANALYSED_File = MyTFile((DIR_ROOT_DATA_ANALYSED + NUCLEUS + "_analysed.root").c_str(), "RECREATE");
    TTree *ANALYSED_Tree = new TTree("ANALYSED_Tree", "ANALYSED_Tree");

    double TH[10] = {0, 300, 300, 300, 400, 200, 400, 400, 400, 300};

    InitCalib();
    InitWindows();
    InitHistograms();
    

    ReaderData();
    RandomCorrection();

    
    for (int mul = 1; mul <= BETA_SIZE; mul++)
    {
        cRatioCoinc_NoCoinc[mul] = new TCanvas(("RatioCoinc_NoCoinc_" + to_string(mul)).c_str(), ("RatioCoinc_NoCoinc_" + to_string(mul)).c_str(), 1920, 1080);
        G_RatioCoinc_NoCoinc[mul] = new TGraph();
        cRatioCoinc_NoCoinc_Corrected[mul] = new TCanvas(("RatioCoinc_NoCoinc_Corrected_" + to_string(mul)).c_str(), ("RatioCoinc_NoCoinc_Corrected_" + to_string(mul)).c_str(), 1920, 1080);
        G_RatioCoinc_NoCoinc_Corrected[mul] = new TGraph();
    }

    /// DRAWING
    ANALYSED_File->cd();
    double all_coinc = 0;
    double all_nocoinc = 0;
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            if (H_Single[i]->GetEntries() == 0)
                continue;
            // cout << endl;
            if (VERBOSE == 1) Info(detectorName[i]);
            for (int mul = 1; mul <= BETA_SIZE; mul++)
            {
                // Info("Multiplicity: " + to_string(mul), 1);
                dir_RandomCorrection_Strip[i]->cd();
                TCanvas *c = new TCanvas(("RAW_" + detectorName[i] + "_" + to_string(mul)).c_str(), ("RAW_" + detectorName[i] + "_" + to_string(mul)).c_str(), 1920, 1080);
                H_Single[i]->SetLineColor(kBlack);
                H_Single[i]->Draw("HIST");
                H_Coinc[i][mul]->SetLineColor(kRed);
                H_Coinc[i][mul]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][IAS[NUCLEUS]][i].first, WindowsMap[NUCLEUS][IAS[NUCLEUS]][i].second);
                H_Coinc[i][mul]->Draw("HIST SAME");
                H_NoCoinc[i][mul]->SetLineColor(kBlue);
                H_NoCoinc[i][mul]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][IAS[NUCLEUS]][i].first, WindowsMap[NUCLEUS][IAS[NUCLEUS]][i].second);
                H_NoCoinc[i][mul]->Draw("HIST SAME");
                H_Random[i][mul]->SetLineColor(kGreen);
                H_Random[i][mul]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][IAS[NUCLEUS]][i].first, WindowsMap[NUCLEUS][IAS[NUCLEUS]][i].second);
                H_Random[i][mul]->Draw("HIST SAME");
                c->Write();

                TCanvas *c_Corrected = new TCanvas(("Corrected_" + detectorName[i] + "_" + to_string(mul)).c_str(), ("Corrected_" + detectorName[i] + "_" + to_string(mul)).c_str(), 1920, 1080);
                H_Single[i]->SetLineColor(kBlack);
                H_Single[i]->Draw("HIST");
                H_Coinc_Corrected[i][mul]->SetLineColor(kRed);
                H_Coinc_Corrected[i][mul]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][IAS[NUCLEUS]][i].first, WindowsMap[NUCLEUS][IAS[NUCLEUS]][i].second);
                H_Coinc_Corrected[i][mul]->Draw("HIST SAME");
                H_NoCoinc_Corrected[i][mul]->SetLineColor(kBlue);
                H_NoCoinc_Corrected[i][mul]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][IAS[NUCLEUS]][i].first, WindowsMap[NUCLEUS][IAS[NUCLEUS]][i].second);
                H_NoCoinc_Corrected[i][mul]->Draw("HIST SAME");
                c_Corrected->Write();

                TCanvas *c_Time = new TCanvas(("Time_" + detectorName[i] + "_" + to_string(mul)).c_str(), ("Time_" + detectorName[i] + "_" + to_string(mul)).c_str(), 1920, 1080);
                H_SiPM_Time[i][mul]->Draw("HIST");
                H_SiPM_Time_Coinc[i][mul]->SetLineColor(kRed);
                H_SiPM_Time_Coinc[i][mul]->Draw("HIST SAME");
                c_Time->Write();

                // printing //
                
                // Info("## RAW ##", 2);
                int nocoinc = H_NoCoinc[i][mul]->Integral();
                int coinc = H_Coinc[i][mul]->Integral();
                int single = H_Single[i]->Integral();
                // if (mul == 9)
                    // Info("Single: " + to_string(single), 3);
            
                // Info("No Coinc: " + to_string(nocoinc), 3);
                // Info("Coinc: " + to_string(coinc), 3);
                // Info("Coinc/NoCoinc: " + to_string(coinc / (double)nocoinc), 3);
                G_RatioCoinc_NoCoinc[mul]->AddPoint(i, coinc / (double)nocoinc);
                // Info("## CORRECTED ##", 2);
                // Info("## Random: " + to_string(NRandom[i][mul]), 3);
                nocoinc = H_NoCoinc_Corrected[i][mul]->Integral();
                coinc = H_Coinc_Corrected[i][mul]->Integral();
                // Info("No Coinc: " + to_string(nocoinc), 3);
                // Info("Coinc: " + to_string(coinc), 3);
                // Info("Coinc/NoCoinc: " + to_string(coinc / (double)nocoinc), 3);
                G_RatioCoinc_NoCoinc_Corrected[mul]->AddPoint(i, coinc / (double)nocoinc);
                // Info("## SHIFT CORRECTED ##", 2);
                // Info("No Coinc: " + to_string(H_NoCoinc_Corrected[i][mul]->GetMean()), 3);
                // Info("Coinc: " + to_string(H_Coinc_Corrected[i][mul]->GetMean()), 3);
                // Info("E: " + to_string(0.5 * abs(H_NoCoinc_Corrected[i][mul]->GetMean() - H_Coinc_Corrected[i][mul]->GetMean())), 3);

                dir_RandomCorrection_Strip_Write[i]->cd();
                H_NoCoinc[i][mul]->Write();
                H_Coinc[i][mul]->Write();
                H_NoCoinc_Corrected[i][mul]->Write();
                H_Coinc_Corrected[i][mul]->Write();
                H_SiPM_Time[i][mul]->Write();
                H_SiPM_Time_Coinc[i][mul]->Write();

                // all_coinc+= H_Coinc_Corrected[i][mul]->Integral();
                // all_nocoinc += H_NoCoinc_Corrected[i][mul]->Integral();
            }

            // Compare coinc with multiplicity
            dir_RandomCorrection_Strip[i]->cd();
            TCanvas *c = new TCanvas(("H_Coinc_AtLeastMulitplicity_" + detectorName[i]).c_str(), ("H_Coinc_AtLeastMulitplicity_" + detectorName[i]).c_str(), 1920, 1080);
            TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);
            H_Single[i]->SetLineColor(kBlack);
            H_Single[i]->Draw("HIST");
            for (int mul = 1; mul <= BETA_SIZE; mul++)
            {
                H_Coinc[i][mul]->SetLineColor(mul);
                H_Coinc[i][mul]->Draw("HIST SAME");
                legend->AddEntry(H_Coinc[i][mul], ("Multiplicity " + to_string(mul)).c_str(), "l");
            }
            legend->Draw("SAME");
            c->Write();

            TCanvas *c1 = new TCanvas(("H_Coinc_EqualMulitplicity_" + detectorName[i]).c_str(), ("H_Coinc_EqualMulitplicity_" + detectorName[i]).c_str(), 1920, 1080);
            TLegend *legend1 = new TLegend(0.1, 0.7, 0.48, 0.9);
            H_Single[i]->SetLineColor(kBlack);
            H_Single[i]->Draw("HIST");
            for (int mul = 1; mul <= BETA_SIZE; mul++)
            {
                H_Coinc_Mulitplicity[i][mul]->SetLineColor(mul);
                H_Coinc_Mulitplicity[i][mul]->Draw("HIST SAME");
                legend1->AddEntry(H_Coinc_Mulitplicity[i][mul], ("Multiplicity " + to_string(mul)).c_str(), "l");
            }
            legend1->Draw("SAME");
            c1->Write();

            H_Time_Channel[i]->Write();
        }
    }

    cout << "All coinc: " << all_coinc << endl;
    cout << "All nocoinc: " << all_nocoinc << endl;
    cout << "ratio: " << all_coinc / all_nocoinc << endl;

    ANALYSED_File->cd();
    // Compare NRandom with multiplicity
    for (int mul = 1; mul <= BETA_SIZE; mul++)
    {
        TCanvas *c = new TCanvas(("H_NRandom_" + to_string(mul)).c_str(), ("H_NRandom_" + to_string(mul)).c_str(), 1920, 1080);
        G_NRandom[mul]->SetMarkerStyle(20);
        G_NRandom[mul]->SetMarkerSize(1);
        G_NRandom[mul]->GetXaxis()->SetTitle("Strip");
        G_NRandom[mul]->GetYaxis()->SetTitle("Random coincidences");
        G_NRandom[mul]->Draw("AP");
        c->Write();
    }

    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBeta(i))
        {
            TCanvas *c = new TCanvas(("H_SiPM_Channel_nocoinc" + detectorName[i]).c_str(), ("H_SiPM_Channel_nocoinc" + detectorName[i]).c_str(), 1920, 1080);

            TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);
            c->cd();
            for (int m = 1; m <= MAX_MULTIPLICTY; m++)
            {
                H_SiPM_Channel_M_nocoinc[i][m]->SetLineColor(m);
                if (m == 1)
                    H_SiPM_Channel_M_nocoinc[i][m]->Draw("HIST");
                else
                    H_SiPM_Channel_M_nocoinc[i][m]->Draw("HIST SAME");
                legend->AddEntry(H_SiPM_Channel_M_nocoinc[i][m], ("Multiplicity " + to_string(m)).c_str(), "l");
            }
            legend->Draw("SAME");
            c->Write();

            TCanvas *c1 = new TCanvas(("H_SiPM_Channel_coinc" + detectorName[i]).c_str(), ("H_SiPM_Channel_coinc" + detectorName[i]).c_str(), 1920, 1080);
            TLegend *legend1 = new TLegend(0.1, 0.7, 0.48, 0.9);
            c1->cd();
            for (int m = 1; m <= MAX_MULTIPLICTY; m++)
            {
                H_SiPM_Channel_M_coinc[i][m]->SetLineColor(m);
                if (m == 1)
                    H_SiPM_Channel_M_coinc[i][m]->Draw("HIST");
                else
                    H_SiPM_Channel_M_coinc[i][m]->Draw("HIST SAME");
                legend1->AddEntry(H_SiPM_Channel_M_coinc[i][m], ("Multiplicity " + to_string(m)).c_str(), "l");
            }
            legend1->Draw("SAME");
            c1->Write();
        }
    }

    // for each multiplciity
    for (int m = 1; m <= MAX_MULTIPLICTY; m++)
    {
        TCanvas *cH = new TCanvas(("H_SiPMHigh_Channel_coinc_M" + to_string(m)).c_str(), ("H_SiPMHigh_Channel_coinc_M" + to_string(m)).c_str(), 1920, 1080);
        TCanvas *cL = new TCanvas(("H_SiPMLow_Channel_coinc_M" + to_string(m)).c_str(), ("H_SiPMLow_Channel_coinc_M" + to_string(m)).c_str(), 1920, 1080);
        TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorBetaHigh(i))
            {
                cH->cd();
                H_SiPM_Channel_M_coinc[i][m]->SetLineColor(GetDetectorChannel(i));
                if (i == 0)
                    H_SiPM_Channel_M_coinc[i][m]->Draw("HIST");
                else
                    H_SiPM_Channel_M_coinc[i][m]->Draw("HIST SAME");

                legend->AddEntry(H_SiPM_Channel_M_coinc[i][m], ("SiPM" + to_string(GetDetectorChannel(i))).c_str(), "l");
            }

            if (IsDetectorBetaLow(i))
            {
                cL->cd();
                H_SiPM_Channel_M_coinc[i][m]->SetLineColor(GetDetectorChannel(i));
                if (i == 0)
                    H_SiPM_Channel_M_coinc[i][m]->Draw("HIST");
                else
                    H_SiPM_Channel_M_coinc[i][m]->Draw("HIST SAME");
            }
        }
        cH->cd();
        legend->Draw("SAME");
        cH->Write();
        cL->cd();
        legend->Draw("SAME");
        cL->Write();
    }

    int counter = 0;

    for (int mul = 1; mul <= BETA_SIZE; mul++)
    {
        cRatioCoinc_NoCoinc[mul]->cd();
        G_RatioCoinc_NoCoinc[mul]->SetMarkerStyle(20);
        G_RatioCoinc_NoCoinc[mul]->SetMarkerSize(1);
        G_RatioCoinc_NoCoinc[mul]->GetXaxis()->SetTitle("Strip");
        G_RatioCoinc_NoCoinc[mul]->GetYaxis()->SetTitle("Ratio Coinc/NoCoinc");
        G_RatioCoinc_NoCoinc[mul]->Draw("AP");
        cRatioCoinc_NoCoinc[mul]->Write();

        cRatioCoinc_NoCoinc_Corrected[mul]->cd();
        G_RatioCoinc_NoCoinc_Corrected[mul]->SetMarkerStyle(20);
        G_RatioCoinc_NoCoinc_Corrected[mul]->SetMarkerSize(1);
        G_RatioCoinc_NoCoinc_Corrected[mul]->GetXaxis()->SetTitle("Strip");
        G_RatioCoinc_NoCoinc_Corrected[mul]->GetYaxis()->SetTitle("Ratio Coinc/NoCoinc Corrected");
        G_RatioCoinc_NoCoinc_Corrected[mul]->Draw("AP");
        cRatioCoinc_NoCoinc_Corrected[mul]->Write();
    }

    // COUNTS
    TGraphErrors *G_Single_strip = new TGraphErrors();
    TGraphErrors *G_Coinc_strip = new TGraphErrors();
    TGraphErrors *G_Single_det = new TGraphErrors();
    for (int det = 0; det <= SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliStrip(det))
        {
            if (VERBOSE == 1) Info("Detector: " + detectorName[det], 2);

            //SINGLE
            double single = H_Single[det]->Integral();
            double single_error = sqrt(single);
            G_Single_strip->AddPoint(det, single);
            G_Single_strip->SetPointError(G_Single_strip->GetN() - 1, 0, single_error);

            if (GetDetectorChannel(det) == 1)
            {
                G_Single_det->SetPoint(GetDetector(det)-1, GetDetector(det), single);
            }
            else
            {
                double s = G_Single_det->GetPointY(GetDetector(det)-1);
                G_Single_det->SetPoint(GetDetector(det)-1, GetDetector(det), s + single);
                G_Single_det->SetPointError(GetDetector(det)-1, 0, sqrt(s + single));
            }
            
            //COINC
            double coinc = H_Coinc[det][3]->Integral();
            double coinc_error = sqrt(coinc);
            G_Coinc_strip->AddPoint(det, coinc);
            G_Coinc_strip->SetPointError(G_Coinc_strip->GetN() - 1, 0, coinc_error);

            // if (GetDetectorChannel(det) == 1)
            // {
            //     G_Coinc_det->SetPoint(GetDetector(det)-1, GetDetector(det), coinc);
            // }
            // else
            // {
            //     double s = G_Coinc_det->GetPointY(GetDetector(det)-1);
            //     G_Coinc_det->SetPoint(GetDetector(det)-1, GetDetector(det), s + coinc);
            //     G_Coinc_det->SetPointError(GetDetector(det)-1, 0, sqrt(s + coinc));
            // }
        }
    }

    TCanvas *c = new TCanvas("H_Single_Strip", "H_Single_Strip", 1920, 1080);
    G_Single_strip->SetMarkerStyle(20);
    G_Single_strip->SetMarkerSize(1);
    G_Single_strip->GetXaxis()->SetTitle("Strip");
    G_Single_strip->GetYaxis()->SetTitle("Counts");
    G_Single_strip->Draw("AP");
    c->Write();

    TCanvas *c1 = new TCanvas("H_Coinc_Strip", "H_Coinc_Strip", 1920, 1080);
    G_Coinc_strip->SetMarkerStyle(20);
    G_Coinc_strip->SetMarkerSize(1);
    G_Coinc_strip->GetXaxis()->SetTitle("Strip");
    G_Coinc_strip->GetYaxis()->SetTitle("Counts");
    G_Coinc_strip->Draw("AP");
    c1->Write();

    TCanvas *c2 = new TCanvas("H_Single_Det", "H_Single_Det", 1920, 1080);
    G_Single_det->SetMarkerStyle(20);
    G_Single_det->SetMarkerSize(1);
    G_Single_det->GetXaxis()->SetTitle("Detector");
    G_Single_det->GetYaxis()->SetTitle("Counts");
    G_Single_det->Draw("AP");
    c2->Write();

    ///EHISFT Calculation
    TDirectory *dir_Eshift = ANALYSED_File->mkdir("Eshift");
    TDirectory *dir_Eshift_M[BETA_SIZE+1] = {nullptr};
    Start("Compute Eshift");
    TGraphErrors *G_Eshift[BETA_SIZE+1] = {nullptr};
    TGraphErrors *G_Eshift_det[2][SILI_SIZE][BETA_SIZE+1] = {nullptr};
    int dir;
    for (int det = 0 ; det <= SIGNAL_MAX ; det++)
    {
        if (IsDetectorSiliStrip(det))
        {
            Info("Detector: " + detectorName[det], 2);
            for (int mul = 1; mul <= BETA_SIZE; mul++)
            {
                pair<double, double> E = ComputeEshift(det, H_Single[det], H_Coinc_Corrected[det][mul], H_NoCoinc_Corrected[det][mul]);
                double Eshift = E.first;
                double Eshift_error = E.second;
                
                // printing
                if (mul == 3) Info("Eshift: " + to_string(Eshift) + " +/- " + to_string(Eshift_error) + " keV");

                // global
                // Info("Global");
                if (G_Eshift[mul] == nullptr)  G_Eshift[mul] = new TGraphErrors();
                G_Eshift[mul]->AddPoint(det, Eshift);
                G_Eshift[mul]->SetPointError(G_Eshift[mul]->GetN()-1, 0, Eshift_error);

                // per hemisphere and strip
                // Info("Hemi");
                if (GetDetector(det) < 5 ) dir = 0;
                else dir = 1;
                // cout << "dir: " << dir << " Strip: " << GetDetectorChannel(det) << " mul: " << mul << endl;
                if (G_Eshift_det[dir][GetDetectorChannel(det)][mul] == nullptr) G_Eshift_det[dir][GetDetectorChannel(det)][mul] = new TGraphErrors();
                G_Eshift_det[dir][GetDetectorChannel(det)][mul]->AddPoint(GetDetector(det), Eshift);
                G_Eshift_det[dir][GetDetectorChannel(det)][mul]->SetPointError(G_Eshift_det[dir][GetDetectorChannel(det)][mul]->GetN()-1, 0, Eshift_error);
            }
        }
    }


    TMultiGraph *MG_Eshift[SILI_SIZE][BETA_SIZE+1] = {nullptr};
    TLine *line_sigma_plus[2][SILI_SIZE][BETA_SIZE+1] = {nullptr};
    TLine *line_sigma_minus[2][SILI_SIZE][BETA_SIZE+1] = {nullptr};
    TLine *line[2][SILI_SIZE][BETA_SIZE+1] = {nullptr};

    for (int mul = 1; mul <= BETA_SIZE; mul++)
    {
        // Info("Multiplicity: " + to_string(mul), 1);
        dir_Eshift_M[mul] = dir_Eshift->mkdir(("Eshift_M" + to_string(mul)).c_str());
        dir_Eshift_M[mul]->cd();
        TCanvas *c = new TCanvas(("Eshift_" + to_string(mul)).c_str(), ("Eshift_M" + to_string(mul)).c_str(), 1920, 1080);
        G_Eshift[mul]->SetMarkerStyle(20);
        G_Eshift[mul]->SetMarkerSize(1);
        G_Eshift[mul]->GetXaxis()->SetTitle("Strip");
        G_Eshift[mul]->GetYaxis()->SetTitle("Eshift");
        G_Eshift[mul]->Draw("AP");
        c->Write();

        TCanvas *c1[SILI_SIZE];
        for (int strip = 1; strip < SILI_SIZE; strip++)
        {
            // Info("Strip: " + to_string(strip), 2);
            c1[strip] = new TCanvas(("Eshift_M" + to_string(mul) + "_Strip" + to_string(strip)).c_str(), ("Eshift_" + to_string(mul) + "_Strip" + to_string(strip)).c_str(), 1920, 1080);
            MG_Eshift[strip][mul] = new TMultiGraph();
            for (int dir = 0; dir <= 1; dir++)
            {
                // Info("Dir: " + to_string(dir), 3);
                c1[strip]->cd();
                TF1 *f = new TF1("pol0", "pol0", 0, 10);
                G_Eshift_det[dir][strip][mul]->Fit(f, "QN");
                G_Eshift_det[dir][strip][mul]->SetMarkerStyle(20);
                double from = (dir == 0 ? 0.8 : 4.5);
                double to = (dir == 0 ? 4.5 : 8.2);
                
                double Eshift_err = f->GetParError(0);
                // Info("Lines", 4);
                line_sigma_plus[dir][strip][mul] = new TLine(from, f->GetParameter(0) + Eshift_err, to, f->GetParameter(0) + Eshift_err);
                line_sigma_plus[dir][strip][mul]->SetLineColor(kRed);
                line_sigma_plus[dir][strip][mul]->SetLineStyle(2);
                line_sigma_minus[dir][strip][mul] = new TLine(from, f->GetParameter(0) - Eshift_err, to, f->GetParameter(0) - Eshift_err);
                line_sigma_minus[dir][strip][mul]->SetLineColor(kRed);
                line_sigma_minus[dir][strip][mul]->SetLineStyle(2);
                line[dir][strip][mul] = new TLine(from, f->GetParameter(0), to, f->GetParameter(0));
                line[dir][strip][mul]->SetLineColor(kRed);

                MG_Eshift[strip][mul]->Add(G_Eshift_det[dir][strip][mul]);
            }

            // Info("Draw", 3);
            c1[strip]->cd();
            MG_Eshift[strip][mul]->GetYaxis()->SetRangeUser(0, 5);
            MG_Eshift[strip][mul]->GetXaxis()->SetRangeUser(0.8, 8.2);
            MG_Eshift[strip][mul]->Draw("AP");
            line_sigma_minus[0][strip][mul]->Draw("SAME");
            line_sigma_plus[0][strip][mul]->Draw("SAME");
            line[0][strip][mul]->Draw("SAME");
            line_sigma_minus[1][strip][mul]->Draw("SAME");
            line_sigma_plus[1][strip][mul]->Draw("SAME");
            line[1][strip][mul]->Draw("SAME");
            c1[strip]->Write();
        }
    }


    ANALYSED_File->Close();

    return 0;
}