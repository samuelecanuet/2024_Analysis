#include "Detectors.hh"

bool FULL = false;
int Run;

TFile *GROUPED_File;
TTree *CLEANED_Tree;

vector<Signal> CLEANED_Tree_Silicon;
vector<vector<pair<Signal, Signal>>> CLEANED_Tree_SiPMGroup;
// TTreeReaderArray<Signal> *signals;

int TOTAL_COUNTER[BETA_SIZE+1] = {0};
pair<int, int> Rear_IAS[9] = {make_pair(0, 0), make_pair(40000, 44500), make_pair(36000, 39500), make_pair(36000, 39500), make_pair(33000, 36500), make_pair(36000, 39000), make_pair(33000, 36500), make_pair(35500, 38500), make_pair(35000, 38500)};


// HIST
TH1D* H_TOTAL_MULTIPLICITY;
TH1D* H_Channel[SIGNAL_MAX];

TH1D* H_Channel_Coinc_All[SIGNAL_MAX][BETA_SIZE+1];
TH1D* H_Channel_NoCoinc_All[SIGNAL_MAX][BETA_SIZE+1];
TH1D* H_RearGroup_Time_All[SIGNAL_MAX];
TH1D* H_RearSiPM_Time_All[SIGNAL_MAX][SIGNAL_MAX];

TH1D* H_Channel_Coinc_Restricted[SIGNAL_MAX][BETA_SIZE+1];
TH1D* H_Channel_NoCoinc_Restricted[SIGNAL_MAX][BETA_SIZE+1];
TH1D* H_RearGroup_Time_Restricted[SIGNAL_MAX];
TH1D* H_RearSiPM_Time_Restricted[SIGNAL_MAX][SIGNAL_MAX];

TH1D* H_Channel_Coinc_May[SIGNAL_MAX];
TH1D* H_Channel_NoCoinc_May[SIGNAL_MAX];

TDirectory* dir_detector[SIGNAL_MAX];

void InitHistograms()
{
    Info("InitHistograms");
    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliBack(det))
        {
            // Info(detectorName[det], 1);
            dir_detector[det] = GROUPED_File->mkdir(detectorName[det].c_str());
            H_Channel[det] = new TH1D(("H_Channel_" + detectorName[det]).c_str(), ("H_Channel_" + detectorName[det]).c_str(), eSiliN, eSiliMin, eSiliMax);
            
            H_RearGroup_Time_All[det] = new TH1D(("H_RearGroup_Time_" + detectorName[det]).c_str(), ("H_RearGroup_Time_" + detectorName[det]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
            H_RearGroup_Time_Restricted[det] = new TH1D(("H_RearGroup_Time_Restricted_" + detectorName[det]).c_str(), ("H_RearGroup_Time_Restricted_" + detectorName[det]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
            
            H_Channel_Coinc_May[det] = new TH1D(("H_Channel_Coinc_May_" + detectorName[det]).c_str(), ("H_Channel_Coinc_May_" + detectorName[det]).c_str(), eSiliN, eSiliMin, eSiliMax);
            H_Channel_NoCoinc_May[det] = new TH1D(("H_Channel_NoCoinc_May_" + detectorName[det]).c_str(), ("H_Channel_NoCoinc_May_" + detectorName[det]).c_str(), eSiliN, eSiliMin, eSiliMax);
            
            for (int mul = 0; mul <= BETA_SIZE; mul++)
            {
                // Info("Multiplicity: " + to_string(mul), 2);
                H_Channel_Coinc_All[det][mul] = new TH1D(("H_Channel_Coinc_" + detectorName[det] + "_" + to_string(mul)).c_str(), ("H_Channel_Coinc_" + detectorName[det] + "_" + to_string(mul)).c_str(), eSiliN, eSiliMin, eSiliMax);
                H_Channel_NoCoinc_All[det][mul] = new TH1D(("H_Channel_NoCoinc_" + detectorName[det] + "_" + to_string(mul)).c_str(), ("H_Channel_NoCoinc_" + detectorName[det] + "_" + to_string(mul)).c_str(), eSiliN, eSiliMin, eSiliMax);
            
                H_Channel_Coinc_Restricted[det][mul] = new TH1D(("H_Channel_Coinc_Restricted_" + detectorName[det] + "_" + to_string(mul)).c_str(), ("H_Channel_Coinc_Restricted_" + detectorName[det] + "_" + to_string(mul)).c_str(), eSiliN, eSiliMin, eSiliMax);
                H_Channel_NoCoinc_Restricted[det][mul] = new TH1D(("H_Channel_NoCoinc_Restricted_" + detectorName[det] + "_" + to_string(mul)).c_str(), ("H_Channel_NoCoinc_Restricted_" + detectorName[det] + "_" + to_string(mul)).c_str(), eSiliN, eSiliMin, eSiliMax);
            
            }

            for (int det2 = 0; det2 < SIGNAL_MAX; det2++)
            {
                if (IsDetectorBetaHigh(det2))
                {
                    H_RearSiPM_Time_All[det][det2] = new TH1D(("H_RearSiPM_Time_" + detectorName[det] + "_" + detectorName[det2]).c_str(), ("H_RearSiPM_Time_" + detectorName[det] + "_" + detectorName[det2]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
                    H_RearSiPM_Time_Restricted[det][det2] = new TH1D(("H_RearSiPM_Time_Restricted_" + detectorName[det] + "_" + detectorName[det2]).c_str(), ("H_RearSiPM_Time_Restricted_" + detectorName[det] + "_" + detectorName[det2]).c_str(), winGroupN_Beta, winGroupMin_Beta, winGroupMax_Beta);
                }
            }
        }
    }

    H_TOTAL_MULTIPLICITY = new TH1D("H_TOTAL_MULTIPLICITY", "H_TOTAL_MULTIPLICITY", 20, 0, 20);
    Info("InitHistograms Done");
}

void WriteHistograms()
{
    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliBack(det))
        {
            dir_detector[det]->cd();
            H_Channel[det]->Write();
            H_RearGroup_Time_All[det]->Write();

            TCanvas *cMay = new TCanvas(("H_Channel_May_" + detectorName[det]).c_str(), ("H_Channel_May_" + detectorName[det]).c_str(), 1920, 1080);
            H_Channel[det]->SetLineColor(kBlack);
            H_Channel[det]->Draw("HIST");
            H_Channel_Coinc_May[det]->SetLineColor(kRed);
            H_Channel_Coinc_May[det]->Draw("HIST SAME");
            H_Channel_NoCoinc_May[det]->SetLineColor(kBlue);
            H_Channel_NoCoinc_May[det]->Draw("HIST SAME");
            cMay->Write();
            

            /// FULL GROUP TIME
            TCanvas *c_RearGroupTime_All = new TCanvas(("H_RearGroup_Time_" + detectorName[det]).c_str(), ("H_RearGroup_Time_" + detectorName[det]).c_str(), 1920, 1080);
            TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);
            H_RearGroup_Time_All[det]->SetLineColor(kBlack);
            H_RearGroup_Time_All[det]->Draw("HIST");
            for (int det2 = 0; det2 < SIGNAL_MAX; det2++)
            {
                if (IsDetectorBetaHigh(det2))
                {
                    H_RearSiPM_Time_All[det][det2]->SetLineColor(GetDetectorChannel(det2)+1);
                    H_RearSiPM_Time_All[det][det2]->Draw("HIST SAME");
                    legend->AddEntry(H_RearSiPM_Time_All[det][det2], detectorName[det2].c_str(), "l");
                }
            }
            legend->Draw("SAME");
            c_RearGroupTime_All->Write();

            for (int mul = 0; mul <= BETA_SIZE; mul++)
            {
                H_Channel_Coinc_All[det][mul]->Write();
                H_Channel_NoCoinc_All[det][mul]->Write();

                TCanvas *c = new TCanvas(("H_Channel_" + detectorName[det] + "_" + to_string(mul)).c_str(), ("H_Channel_" + detectorName[det] + "_" + to_string(mul)).c_str(), 1920, 1080);
                H_Channel[det]->SetLineColor(kBlack);
                H_Channel[det]->Draw("HIST");
                H_Channel_Coinc_All[det][mul]->SetLineColor(kRed);
                H_Channel_Coinc_All[det][mul]->Draw("HIST SAME");
                H_Channel_NoCoinc_All[det][mul]->SetLineColor(kBlue);
                H_Channel_NoCoinc_All[det][mul]->Draw("HIST SAME");
                c->Write();
            }
            ///

            /// TIME RESTRICTED
            TCanvas *c_RearGroupTime_Restricted = new TCanvas(("H_RearGroup_Time_Restricted_" + detectorName[det]).c_str(), ("H_RearGroup_Time_Restricted_" + detectorName[det]).c_str(), 1920, 1080);
            TLegend *legend_Restricted = new TLegend(0.1, 0.7, 0.48, 0.9);
            H_RearGroup_Time_Restricted[det]->SetLineColor(kBlack);
            H_RearGroup_Time_Restricted[det]->Draw("HIST");
            for (int det2 = 0; det2 < SIGNAL_MAX; det2++)
            {
                if (IsDetectorBetaHigh(det2))
                {
                    H_RearSiPM_Time_Restricted[det][det2]->SetLineColor(GetDetectorChannel(det2)+1);
                    H_RearSiPM_Time_Restricted[det][det2]->Draw("HIST SAME");
                    legend_Restricted->AddEntry(H_RearSiPM_Time_Restricted[det][det2], detectorName[det2].c_str(), "l");
                }
            }
            legend_Restricted->Draw("SAME");
            c_RearGroupTime_Restricted->Write();

            for (int mul = 0; mul <= BETA_SIZE; mul++)
            {
                H_Channel_Coinc_Restricted[det][mul]->Write();
                H_Channel_NoCoinc_Restricted[det][mul]->Write();

                TCanvas *c = new TCanvas(("H_Channel_Restricted_" + detectorName[det] + "_" + to_string(mul)).c_str(), ("H_Channel_Restricted_" + detectorName[det] + "_" + to_string(mul)).c_str(), 1920, 1080);
                H_Channel[det]->SetLineColor(kBlack);
                H_Channel[det]->Draw("HIST");
                H_Channel_Coinc_Restricted[det][mul]->SetLineColor(kRed);
                H_Channel_Coinc_Restricted[det][mul]->Draw("HIST SAME");
                H_Channel_NoCoinc_Restricted[det][mul]->SetLineColor(kBlue);
                H_Channel_NoCoinc_Restricted[det][mul]->Draw("HIST SAME");
                c->Write();
            }
            ///

        }
    }

    GROUPED_File->cd();
    H_TOTAL_MULTIPLICITY->Write();
}

int CleaningGroups(TTreeReaderArray<Signal> &signals, int verbose)
{

    vector<int> Rear_Position;
    vector<int> Strip_Position;

    vector<Signal> Silicon;
    vector<Signal> SiPM_High;
    vector<Signal> SiPM_Low;

    for (int index = 0; index < signals.GetSize(); index++)
    {
        int current_label = signals[index].Label;
        // cout << current_label << endl;
        if (IsDetectorSiliBack(current_label))
        {
            signals[index].Channel = signals[index].Channel - 0.5 + gRandom->Rndm();
            signals[index].Time = signals[index].Time + (-0.5 + gRandom->Rndm()) * 8;
            Rear_Position.push_back(index);
        }
        else if (IsDetectorSiliStrip(current_label))
        {
            signals[index].Channel = signals[index].Channel - 0.5 + gRandom->Rndm();
            signals[index].Time = signals[index].Time + (-0.5 + gRandom->Rndm()) * 8;
            Strip_Position.push_back(index);
        }
        else if (IsDetectorBetaHigh(current_label))
        {
            signals[index].Channel = signals[index].Channel - 0.5 + gRandom->Rndm();
            signals[index].Time = signals[index].Time + (-0.5 + gRandom->Rndm()) * 2;
            SiPM_High.push_back(signals[index]);
        }
        else if (IsDetectorBetaLow(current_label))
        {
            signals[index].Channel = signals[index].Channel - 0.5 + gRandom->Rndm();
            signals[index].Time = signals[index].Time + (-0.5 + gRandom->Rndm()) * 2;
            SiPM_Low.push_back(signals[index]);
        }
    }

    if (Rear_Position.size() != 1)
    {
        return 0;
    }

    Signal Rear = signals[Rear_Position[0]];

    bool IAS = false;
    if (Rear.Channel > Rear_IAS[GetDetector(Rear.Label)].first && Rear.Channel < Rear_IAS[GetDetector(Rear.Label)].second)
    {
        IAS = true;
    }

    if (!IAS)
    {
        return 0;
    }

    /////////////////////////////
    // counting total multiplicity in faster group
    bool SiPM_Triggered[BETA_SIZE+1] = {false};
    for (int i_h = 0; i_h < SiPM_High.size(); i_h++)
    {
        SiPM_Triggered[GetDetectorChannel(SiPM_High[i_h].Label)] = true;
    }


    int TOTAL_MULT = 0;
    for (int i = 0; i <= BETA_SIZE; i++)
    {
        if (SiPM_Triggered[i])
        {
            TOTAL_MULT++;
        }
    }

    H_TOTAL_MULTIPLICITY->Fill(TOTAL_MULT);
    
    TOTAL_COUNTER[TOTAL_MULT]++;
    /////////////////////////////

    /////////////////////////////
    // TIME WITHOUT CORRECTION ALL GROUP TIME //
    double mean_time = 0;
    if (SiPM_High.size() > 0)
    {
        for (int i_h = 0; i_h < SiPM_High.size(); i_h++)
        {
            mean_time += (Rear.Time - SiPM_High[i_h].Time);

            H_RearSiPM_Time_All[Rear.Label][SiPM_High[i_h].Label]->Fill(-(Rear.Time - SiPM_High[i_h].Time));
        }

        mean_time = mean_time / SiPM_High.size();

        H_RearGroup_Time_All[Rear.Label]->Fill(-mean_time);
    }

    // FILLING REAR CHANNEL HIST
    H_Channel[Rear.Label]->Fill(Rear.Channel);
    for (int mul = 0; mul <= BETA_SIZE; mul++)
    {
        if (mul <= TOTAL_MULT)
        {
            H_Channel_Coinc_All[Rear.Label][mul]->Fill(Rear.Channel);
        }
        else
        {
            H_Channel_NoCoinc_All[Rear.Label][mul]->Fill(Rear.Channel);
        }
    }
    ////////////////////////////


    /////////////////////////////
    // TIME WITHOUT CORRECTION WITH TIME RESTRICTION //
    vector<Signal> SiPM_High_Restricted;
    for (int i_h = 0; i_h < SiPM_High.size(); i_h++)
    {
        if (-(Rear.Time - SiPM_High[i_h].Time) > 100 && -(Rear.Time - SiPM_High[i_h].Time) < 200)
        {
            SiPM_High_Restricted.push_back(SiPM_High[i_h]);
        }
    }

    mean_time = 0;
    if (SiPM_High_Restricted.size() > 0)
    {
        for (int i_h = 0; i_h < SiPM_High_Restricted.size(); i_h++)
        {
            mean_time += (Rear.Time - SiPM_High_Restricted[i_h].Time);

            H_RearSiPM_Time_Restricted[Rear.Label][SiPM_High_Restricted[i_h].Label]->Fill(-(Rear.Time - SiPM_High_Restricted[i_h].Time));
        }

        mean_time = mean_time / SiPM_High_Restricted.size();

        H_RearGroup_Time_Restricted[Rear.Label]->Fill(-mean_time);
    }

    for (int mul = 0; mul <= BETA_SIZE; mul++)
    {
        if (mul <= SiPM_High_Restricted.size())
        {
            H_Channel_Coinc_Restricted[Rear.Label][mul]->Fill(Rear.Channel);
        }
        else
        {
            H_Channel_NoCoinc_Restricted[Rear.Label][mul]->Fill(Rear.Channel);
        }
    }

    /////////////////////////////


    /////////////////////////////
    // MAY 2024
    if (SiPM_High.size() == 9)
    {
        H_Channel_Coinc_May[Rear.Label]->Fill(Rear.Channel);
    }
    else
    {
        H_Channel_NoCoinc_May[Rear.Label]->Fill(Rear.Channel);
    }
    /////////////////////////////
    
}


