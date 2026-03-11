#include "Detectors.hh"

int VERBOSE = 0;
string NUCLEUS;

TFile *ANALYSED_File;
TFile *SIM_File;

TTreeReader *Reader;    
TTreeReaderValue<Signal> *Silicon;
TTreeReaderArray<Signal> *SiPM;

TF1* F_Resolution_SiPM[SIGNAL_MAX];

TH1D *H_Single[SIGNAL_MAX];
TH1D *H_Coinc[SIGNAL_MAX];
TH1D *H_NoCoinc[SIGNAL_MAX];
TH1D* H_SiPM[SIGNAL_MAX];
TH1D* H_SiPM_Conv[SIGNAL_MAX][MAX_MULTIPLICTY];

TGraphErrors *G_Eshift = new TGraphErrors();
TGraphErrors *G_Ratio = new TGraphErrors();




void InitResolution_SiPM()
{
    Info("Init Resolution SiPM start");

    TFile *CALIBRATED_File = MyTFile((DIR_ROOT_DATA_CALIBRATED + "SiPM_Calibrated_" + to_string(YEAR) + ".root").c_str(), "READ");

    for (int sipm = 1; sipm <= 9; sipm++)
    {
        F_Resolution_SiPM[sipm] = (TF1 *)CALIBRATED_File->Get(("SiPM_" + to_string(sipm) + "/F_Resolution_SiPM" + to_string(sipm)).c_str());
        if (F_Resolution_SiPM[sipm] == NULL)
        {
            Error("No SiPM resolution found for SiPM: " + to_string(sipm));
        }
    }

    Info("Init Resolution SiPM end");
}

void InitHistograms()
{
    Info("Init Histograms start");

    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliStrip(det))
        {
            H_Single[det] = new TH1D(("H_Single_" + detectorName[det]).c_str(), ("H_Single_" + detectorName[det]).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
            H_Single[det]->GetXaxis()->SetTitle("Energy [keV]");
            H_Single[det]->GetYaxis()->SetTitle("Counts");
            H_Single[det]->GetXaxis()->CenterTitle();
            H_Single[det]->GetYaxis()->CenterTitle();

            H_Coinc[det] = new TH1D(("H_Coinc_" + detectorName[det]).c_str(), ("H_Coinc_" + detectorName[det]).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
            H_Coinc[det]->GetXaxis()->SetTitle("Energy [keV]");
            H_Coinc[det]->GetYaxis()->SetTitle("Counts");
            H_Coinc[det]->GetXaxis()->CenterTitle();
            H_Coinc[det]->GetYaxis()->CenterTitle();

            H_NoCoinc[det] = new TH1D(("H_NoCoinc_" + detectorName[det]).c_str(), ("H_NoCoinc_" + detectorName[det]).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
            H_NoCoinc[det]->GetXaxis()->SetTitle("Energy [keV]");
            H_NoCoinc[det]->GetYaxis()->SetTitle("Counts");
            H_NoCoinc[det]->GetXaxis()->CenterTitle();
            H_NoCoinc[det]->GetYaxis()->CenterTitle();
        }


        if (IsDetectorBetaHigh(det))
        {
            H_SiPM[det] = new TH1D(("H_SiPM_" + detectorName[det]).c_str(), ("H_SiPM_" + detectorName[det]).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
            H_SiPM[det]->GetXaxis()->SetTitle("Energy [keV]");
            H_SiPM[det]->GetYaxis()->SetTitle("Counts");
            H_SiPM[det]->GetXaxis()->CenterTitle();
            H_SiPM[det]->GetYaxis()->CenterTitle();

            for (int mul = 1; mul <= BETA_SIZE; mul++)
            {
                H_SiPM_Conv[GetDetectorChannel(det)][mul] = new TH1D(("H_SiPM_Conv_" + detectorName[det] + "_M" + to_string(mul)).c_str(), ("H_SiPM_Conv_" + detectorName[det] + "_M" + to_string(mul)).c_str(), eSiliN_cal/10, eSiliMin_cal, eSiliMax_cal);
                H_SiPM_Conv[GetDetectorChannel(det)][mul]->GetXaxis()->SetTitle("Energy [keV]");
                H_SiPM_Conv[GetDetectorChannel(det)][mul]->GetYaxis()->SetTitle("Counts");
                H_SiPM_Conv[GetDetectorChannel(det)][mul]->GetXaxis()->CenterTitle();
                H_SiPM_Conv[GetDetectorChannel(det)][mul]->GetYaxis()->CenterTitle();
            }
        }
    }

    Info("Init Histograms end");
}

pair<double, double> ComputeEshift(int det, TH1D *H_Single, TH1D *H_Coinc, TH1D *H_NoCoinc)
{
    H_Single->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][IAS[NUCLEUS]][det].first, WindowsMap[NUCLEUS][IAS[NUCLEUS]][det].second);
    H_Coinc->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][IAS[NUCLEUS]][det].first, WindowsMap[NUCLEUS][IAS[NUCLEUS]][det].second);
    H_NoCoinc->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][IAS[NUCLEUS]][det].first, WindowsMap[NUCLEUS][IAS[NUCLEUS]][det].second);

    double Ns = H_Single->Integral();
    double Nc = H_Coinc->Integral();
    double Nnc = H_NoCoinc->Integral();

    // cout << "Ns: " << Ns << endl;
    // cout << "Nc: " << Nc << endl;
    // cout << "Nnc: " << Nnc << endl;

    double Es = H_Single->GetMean();
    double Ec = H_Coinc->GetMean();
    double Enc = H_NoCoinc->GetMean();

    // cout << "Es: " << Es << endl;
    // cout << "Ec: " << Ec << endl;
    // cout << "Enc: " << Enc << endl;

    double SEs = H_Single->GetMeanError();
    double SEc = H_Coinc->GetMeanError();
    double SEnc = H_NoCoinc->GetMeanError();

    // cout << "SEs: " << SEs << endl;
    // cout << "SEc: " << SEc << endl;
    // cout << "SEnc: " << SEnc << endl;


    // ESHIFT
    double Eshift = abs(Es - Ec);

    // ESHIFT EROOR // ## APPROXIMATED

    double Eshift_error_1A = 0.5 * pow(SEc, 2);
    double Eshift_error_1B = 1./Ns * pow(Eshift, 2);
    double Eshift_error_1 = sqrt(Eshift_error_1A + Eshift_error_1B);



    double Eshift_error_2A = pow( 1. - Nc / Ns , 2) * (pow(SEc, 2) + pow(SEnc, 2));
    double Eshift_error_2B = Nc*Nnc/pow(Ns, 3) * pow(Eshift, 2);
    double Eshift_error_2 = sqrt(Eshift_error_2A + Eshift_error_2B);


    // cout << "Eshift: " << Eshift << endl;

    // cout << "Eshift_error_1: " << Eshift_error_1 << endl;
    // cout << "Eshift_error_2: " << Eshift_error_2 << endl;

    // cout << "Eshift_error_1A: " << sqrt(Eshift_error_1A) << endl;
    // cout << "Eshift_error_2A: " << sqrt(Eshift_error_2A) << endl;

    // cout << "Eshift_error_1B: " << sqrt(Eshift_error_1B) << endl;
    // cout << "Eshift_error_2B: " << sqrt(Eshift_error_2B) << endl;
        
    
    return make_pair(Eshift, Eshift_error_2);
}

void WriteHistograms()
{

    Info("Writing Histograms");
    ANALYSED_File->cd();
    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliStrip(det))
        {
            TCanvas *c = new TCanvas(("H_Single_" + detectorName[det]).c_str(), ("H_Single_" + detectorName[det]).c_str(), 1920, 1080);
            H_Single[det]->SetLineColor(kBlack);
            H_Single[det]->Draw("HIST");
            H_NoCoinc[det]->SetLineColor(kBlue);
            H_NoCoinc[det]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][IAS[NUCLEUS]][det].first, WindowsMap[NUCLEUS][IAS[NUCLEUS]][det].second);
            H_NoCoinc[det]->Draw("HIST SAME");
            H_Coinc[det]->SetLineColor(kRed);
            H_Coinc[det]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][IAS[NUCLEUS]][det].first, WindowsMap[NUCLEUS][IAS[NUCLEUS]][det].second);
            H_Coinc[det]->Draw("HIST SAME");
            c->Write();

            pair<double, double> Eshift_result = ComputeEshift(det, H_Single[det], H_Coinc[det], H_NoCoinc[det]);
            double Eshift = Eshift_result.first;
            double Eshift_error = Eshift_result.second;

            G_Eshift->AddPoint(det, Eshift);
            G_Eshift->SetPointError(G_Eshift->GetN() - 1, 0, Eshift_error);

            
            double Nc = H_Coinc[det]->Integral();
            double Nnc = H_NoCoinc[det]->Integral();
            double ratio = Nc / (Nc + Nnc);
            double ratio_err = sqrt( pow ( sqrt(Nc) * (Nnc)/( (Nc + Nnc)*(Nc + Nnc) ),2) + pow( sqrt(Nnc) * Nc / ( (Nc + Nnc)*(Nc + Nnc) ),2) );

            G_Ratio->AddPoint(det, ratio);
            G_Ratio->SetPointError(G_Ratio->GetN() - 1, 0, ratio_err);


        }
    }

    TCanvas *cEshift = new TCanvas("Eshift", "Eshift", 1920, 1080);
    G_Eshift->SetMarkerStyle(20);
    G_Eshift->SetMarkerSize(1);
    G_Eshift->GetXaxis()->SetTitle("Strip");
    G_Eshift->GetYaxis()->SetTitle("Eshift [keV]");
    G_Eshift->Draw("AP");
    cEshift->Write();

    TCanvas *cRatio = new TCanvas("Ratio_Coinc_Single", "Ratio_Coinc_Single", 1920, 1080);
    G_Ratio->SetMarkerStyle(20);
    G_Ratio->SetMarkerSize(1);
    G_Ratio->GetXaxis()->SetTitle("Strip");
    G_Ratio->GetYaxis()->SetTitle("Ratio Coinc/Single");
    G_Ratio->Draw("AP");
    cRatio->Write();

    for (int mul = 1; mul <= BETA_SIZE; mul++)
    {
        TCanvas *cSiPM_Mul = new TCanvas(("SiPM_Energy_Calibrated_Multiplicity_" + to_string(mul)).c_str(), ("SiPM_Energy_Calibrated_Multiplicity_" + to_string(mul)).c_str(), 1920, 1080);
        for (int sipm = 1; sipm <= 9; sipm++)
        {
            H_SiPM_Conv[sipm][mul]->SetLineColor(sipm);
            if (sipm == 1)
                H_SiPM_Conv[sipm][mul]->Draw("HIST");
            else
                H_SiPM_Conv[sipm][mul]->Draw("HIST SAME");
        }
        cSiPM_Mul->Write();
    }
    // TCanvas *cBeta = new TCanvas("SiPM_Energy_Calibrated", "SiPM_Energy_Calibrated", 1920, 1080);
    // H_SiPM[104]->SetLineColor(kBlack);
    // H_SiPM[104]->Draw("HIST");
    // H_SiPM_Conv[104]->SetLineColor(kRed);
    // H_SiPM_Conv[104]->Draw("HIST SAME");
    // cBeta->Write();
}



void ReaderData()
{

    Start("Reader Data");

    Reader = new TTreeReader((TTree*)SIM_File->Get("Tree"));
    Silicon = new TTreeReaderValue<Signal>(*Reader, "Silicon");
    SiPM = new TTreeReaderArray<Signal>(*Reader, "SiPM");

    clock_t start_time = clock(), Current;

    Long64_t nEntries = Reader->GetEntries();
    while (Reader->Next())
    {
        ProgressBar(Reader->GetCurrentEntry(), nEntries, start_time, Current, "Reading Tree");

        int Strip_Label = (**Silicon).Label;
        double energy = (**Silicon).Channel;// + gRandom->Gaus(0., 10./2.35);
        double Time_Silicon = (**Silicon).Time;

        H_Single[Strip_Label]->Fill(energy);


        bool coinc = false;
        double SiPM_Energy = 0.;
        double SiPM_Energy_Conv[10];
        int Multiplicity = 0;
        for (size_t i = 0; i < (*SiPM).GetSize(); i++)
        {
            double SiPM_Time = (*SiPM)[i].Time;
            SiPM_Energy = (*SiPM)[i].Channel;
            
            
            // multiplicity condition
            Multiplicity = 0;
            for (int sipm = 1; sipm <= 9; sipm++)
            {
                if (YEAR == 2025 && GetDetectorChannel((*SiPM)[i].Label) == 6) // excluding SiPM6
                    continue;
                SiPM_Energy_Conv[sipm] = gRandom->Gaus(SiPM_Energy, F_Resolution_SiPM[sipm]->Eval(SiPM_Energy));
                if (SiPM_Energy_Conv[sipm] > 100.) // threshold
                    Multiplicity++;
            }

            if (abs((*SiPM)[i].Time - Time_Silicon) < 100. && Multiplicity >=3)
            {
                coinc = true;
                break;
            }
        }


        if (coinc)
        {
            H_Coinc[Strip_Label]->Fill(energy);
            H_SiPM[104]->Fill(SiPM_Energy);
            for (int sipm = 1; sipm <= 9; sipm++)
            {
                for (int mul = 1; mul <= BETA_SIZE; mul++)
                {
                    if (Multiplicity >= mul)
                    {
                        H_SiPM_Conv[sipm][mul]->Fill(SiPM_Energy_Conv[sipm]);
                    }
                }
            }
        }
        else
        {
            H_NoCoinc[Strip_Label]->Fill(energy);
        }
    }

    Info("Reader Data");

}

