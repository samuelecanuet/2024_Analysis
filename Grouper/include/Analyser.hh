#include "Detectors.hh"

int VERBOSE = 0;

TFile *ANALYSED_File;
TFile *MERGED_File;
TFile *CALIBRATED_File;

TTreeReader *Reader;
TTreeReaderArray<Signal> *Silicon;
TTreeReaderValue<vector<vector<pair<Signal, Signal>>>> *SiPM_Groups;
TTreeReaderValue<Signal> *HRS;

string NUCLEUS = "32Ar";

TF1* Calibration[SIGNAL_MAX];
TH1D* H_NoCoinc[SIGNAL_MAX][MAX_MULTIPLICTY];
TH1D* H_Coinc[SIGNAL_MAX][MAX_MULTIPLICTY];
TH1D* H_Coinc_Mulitplicity[SIGNAL_MAX][MAX_MULTIPLICTY];
TH1D *H_Coinc_Mulitplicity_Corrected[SIGNAL_MAX][MAX_MULTIPLICTY];
TH1D *H_Single[SIGNAL_MAX];
TH1D *H_NoCoinc_Corrected[SIGNAL_MAX][MAX_MULTIPLICTY];
TH1D *H_Coinc_Corrected[SIGNAL_MAX][MAX_MULTIPLICTY];
TH1D *H_SiPM_Time[SIGNAL_MAX][MAX_MULTIPLICTY];
TH1D *H_SiPM_Time_Coinc[SIGNAL_MAX][MAX_MULTIPLICTY];
TH1D *H_Random[SIGNAL_MAX][MAX_MULTIPLICTY];
TGraph *G_NRandom[MAX_MULTIPLICTY];
TCanvas *cRatioCoinc_NoCoinc[BETA_SIZE];
TGraph *G_RatioCoinc_NoCoinc[BETA_SIZE];
TCanvas *cRatioCoinc_NoCoinc_Corrected[BETA_SIZE];
TGraph *G_RatioCoinc_NoCoinc_Corrected[BETA_SIZE];
TDirectory *dir_RandomCorrection;
TDirectory *dir_RandomCorrection_Strip[SIGNAL_MAX];
TDirectory *dir_RandomCorrection_Strip_Write[SIGNAL_MAX];
TDirectory *dir_RandomCorrection_Strip_Time[SIGNAL_MAX];
double NRandom[SIGNAL_MAX][MAX_MULTIPLICTY];
double HRandom[SIGNAL_MAX][MAX_MULTIPLICTY];

TH1D* H_SiPM_Channel_M_nocoinc[SIGNAL_MAX][MAX_MULTIPLICTY+1];
TH1D* H_SiPM_Channel_M_coinc[SIGNAL_MAX][MAX_MULTIPLICTY+1];

TH2D* H_Time_Channel[SIGNAL_MAX];

TH2D* H_Coinc_SiPMChannel[SIGNAL_MAX][MAX_MULTIPLICTY+1];

TH2D* H_Channel_Time_Coinc[SIGNAL_MAX][MAX_MULTIPLICTY+1];

/// params
double start_gate = -20;
double end_gate = 40;
double start_gate_fake = -300;
double end_gate_fake = -100;

int Run=0;
string Run_string;

void InitHistograms()
{
    Info("Init Histograms start");
    dir_RandomCorrection = ANALYSED_File->mkdir("RandomCorrection");
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            dir_RandomCorrection_Strip[i] = dir_RandomCorrection->mkdir(detectorName[i].c_str());
            dir_RandomCorrection_Strip_Write[i] = dir_RandomCorrection_Strip[i]->mkdir("Write");


            H_Single[i] = new TH1D(("H_Single_" + detectorName[i]).c_str(), ("H_Single_" + detectorName[i]).c_str(), 50000, 0, 10000);
            H_Single[i]->GetXaxis()->SetTitle("Energy [keV]");
            H_Single[i]->GetYaxis()->SetTitle("Counts");
            H_Single[i]->GetXaxis()->CenterTitle();
            H_Single[i]->GetYaxis()->CenterTitle();

            H_Time_Channel[i] = new TH2D(("H_Time_Channel_" + detectorName[i]).c_str(), ("H_Time_Channel_" + detectorName[i]).c_str(), 300, -300, 300, (WindowsMap[NUCLEUS][IAS[NUCLEUS]][i].second-WindowsMap[NUCLEUS][IAS[NUCLEUS]][i].first)*10, WindowsMap[NUCLEUS][IAS[NUCLEUS]][i].first, WindowsMap[NUCLEUS][IAS[NUCLEUS]][i].second);
            H_Time_Channel[i]->GetXaxis()->SetTitle("Time [ns]");
            H_Time_Channel[i]->GetYaxis()->SetTitle("Energy [keV]");
            H_Time_Channel[i]->GetXaxis()->CenterTitle();
            H_Time_Channel[i]->GetYaxis()->CenterTitle();


            for (int mul = 1; mul <= BETA_SIZE; mul++)
            {
                H_NoCoinc[i][mul] = new TH1D(("H_NoCoinc_" + detectorName[i] + "_" + to_string(mul)).c_str(), ("H_NoCoinc_" + detectorName[i] + "_" + to_string(mul)).c_str(), 50000, 0, 10000);
                H_NoCoinc[i][mul]->GetXaxis()->SetTitle("Energy [keV]");
                H_NoCoinc[i][mul]->GetYaxis()->SetTitle("Counts");
                H_NoCoinc[i][mul]->GetXaxis()->CenterTitle();
                H_NoCoinc[i][mul]->GetYaxis()->CenterTitle();

                H_Coinc[i][mul] = new TH1D(("H_Coinc_" + detectorName[i] + "_" + to_string(mul)).c_str(), ("H_Coinc_" + detectorName[i] + "_" + to_string(mul)).c_str(), 50000, 0, 10000);
                H_Coinc[i][mul]->GetXaxis()->SetTitle("Energy [keV]");
                H_Coinc[i][mul]->GetYaxis()->SetTitle("Counts");
                H_Coinc[i][mul]->GetXaxis()->CenterTitle();
                H_Coinc[i][mul]->GetYaxis()->CenterTitle();

                H_Coinc_Mulitplicity[i][mul] = new TH1D(("H_Coinc_Mulitplicity_" + detectorName[i] + "_" + to_string(mul)).c_str(), ("H_Coinc_Mulitplicity_" + detectorName[i] + "_" + to_string(mul)).c_str(), 50000, 0, 10000);
                H_Coinc_Mulitplicity[i][mul]->GetXaxis()->SetTitle("Energy [keV]");
                H_Coinc_Mulitplicity[i][mul]->GetYaxis()->SetTitle("Counts");
                H_Coinc_Mulitplicity[i][mul]->GetXaxis()->CenterTitle();
                H_Coinc_Mulitplicity[i][mul]->GetYaxis()->CenterTitle();    

                H_Random[i][mul] = new TH1D(("H_Random_" + detectorName[i] + "_" + to_string(mul)).c_str(), ("H_Random_" + detectorName[i] + "_" + to_string(mul)).c_str(), 50000, 0, 10000);
                H_Random[i][mul]->GetXaxis()->SetTitle("Energy [keV]");
                H_Random[i][mul]->GetYaxis()->SetTitle("Counts");
                H_Random[i][mul]->GetXaxis()->CenterTitle();
                H_Random[i][mul]->GetYaxis()->CenterTitle();         

                H_SiPM_Time[i][mul] = new TH1D(("H_SiPM_Time_" + detectorName[i] + "_" + to_string(mul)).c_str(), ("H_SiPM_Time_" + detectorName[i] + "_" + to_string(mul)).c_str(), 300, -300, 300);
                H_SiPM_Time[i][mul]->GetXaxis()->SetTitle("Time [ns]");
                H_SiPM_Time[i][mul]->GetYaxis()->SetTitle("Counts");
                H_SiPM_Time[i][mul]->GetXaxis()->CenterTitle();
                H_SiPM_Time[i][mul]->GetYaxis()->CenterTitle();

                H_SiPM_Time_Coinc[i][mul] = new TH1D(("H_SiPM_Time_Coinc_" + detectorName[i] + "_" + to_string(mul)).c_str(), ("H_SiPM_Time_Coinc_" + detectorName[i] + "_" + to_string(mul)).c_str(), 300, -300, 300);
                H_SiPM_Time_Coinc[i][mul]->GetXaxis()->SetTitle("Time [ns]");
                H_SiPM_Time_Coinc[i][mul]->GetYaxis()->SetTitle("Counts");
                H_SiPM_Time_Coinc[i][mul]->GetXaxis()->CenterTitle();
                H_SiPM_Time_Coinc[i][mul]->GetYaxis()->CenterTitle();

                

            }
        }

        if (IsDetectorBetaHigh(i))
        {
            for (int m = 1; m <= MAX_MULTIPLICTY; m++)
            {
                H_SiPM_Channel_M_nocoinc[i][m] = new TH1D(("H_SiPM_Channel_nocoinc_" + detectorName[i] + "_M" + to_string(m)).c_str(), ("H_SiPM_Channel_nocoinc_" + detectorName[i] + "_M" + to_string(m)).c_str(), eHighN, eHighMin, eHighMax);
                H_SiPM_Channel_M_nocoinc[i][m]->GetXaxis()->SetTitle("Channel");
                H_SiPM_Channel_M_nocoinc[i][m]->GetYaxis()->SetTitle("Counts");
                H_SiPM_Channel_M_nocoinc[i][m]->GetXaxis()->CenterTitle();
                H_SiPM_Channel_M_nocoinc[i][m]->GetYaxis()->CenterTitle();

                H_SiPM_Channel_M_coinc[i][m] = new TH1D(("H_SiPM_Channel_coinc_" + detectorName[i] + "_M" + to_string(m)).c_str(), ("H_SiPM_Channel_coinc_" + detectorName[i] + "_M" + to_string(m)).c_str(), eHighN, eHighMin, eHighMax);
                H_SiPM_Channel_M_coinc[i][m]->GetXaxis()->SetTitle("Channel");
                H_SiPM_Channel_M_coinc[i][m]->GetYaxis()->SetTitle("Counts");
                H_SiPM_Channel_M_coinc[i][m]->GetXaxis()->CenterTitle();
                H_SiPM_Channel_M_coinc[i][m]->GetYaxis()->CenterTitle();
            }
        }

        if (IsDetectorBetaLow(i))
        {
            for (int m = 1; m <= MAX_MULTIPLICTY; m++)
            {
                H_SiPM_Channel_M_nocoinc[i][m] = new TH1D(("H_SiPM_Channel_nocoinc_" + detectorName[i] + "_M" + to_string(m)).c_str(), ("H_SiPM_Channel_nocoinc_" + detectorName[i] + "_M" + to_string(m)).c_str(), eLowN, eLowMin, eLowMax);
                H_SiPM_Channel_M_nocoinc[i][m]->GetXaxis()->SetTitle("Channel");
                H_SiPM_Channel_M_nocoinc[i][m]->GetYaxis()->SetTitle("Counts");
                H_SiPM_Channel_M_nocoinc[i][m]->GetXaxis()->CenterTitle();
                H_SiPM_Channel_M_nocoinc[i][m]->GetYaxis()->CenterTitle();

                H_SiPM_Channel_M_coinc[i][m] = new TH1D(("H_SiPM_Channel_coinc_" + detectorName[i] + "_M" + to_string(m)).c_str(), ("H_SiPM_Channel_coinc_" + detectorName[i] + "_M" + to_string(m)).c_str(), eLowN, eLowMin, eLowMax);
                H_SiPM_Channel_M_coinc[i][m]->GetXaxis()->SetTitle("Channel");
                H_SiPM_Channel_M_coinc[i][m]->GetYaxis()->SetTitle("Counts");
                H_SiPM_Channel_M_coinc[i][m]->GetXaxis()->CenterTitle();
                H_SiPM_Channel_M_coinc[i][m]->GetYaxis()->CenterTitle();
            }
        }
    }

    for (int mul = 1; mul <= BETA_SIZE; mul++)
    {
        G_NRandom[mul] = new TGraph();
    }

    Info("Init Histograms end");
}

void InitCalib()
{
    CALIBRATED_File = MyTFile((DIR_ROOT_DATA_CALIBRATED + "Calibrated_" + to_string(YEAR) + ".root").c_str(), "READ");
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            Calibration[i] = (TF1 *)CALIBRATED_File->Get(("Calibration_" + detectorName[i]).c_str());

            if (Calibration[i] == NULL)
            {
                Calibration[i] = new TF1(("Calibration_" + detectorName[i]).c_str(), "x", 0, 10000);
                Warning("No calibration found for " + detectorName[i]);
            }
        }
    }

    Success("Calibration loaded");
}

void ReaderData()
{
    Info("Reading data");
    Reader = new TTreeReader((TTree *)MERGED_File->Get("MERGED_Tree"));
    Silicon = new TTreeReaderArray<Signal>(*Reader, "MERGED_Tree_Silicon");
    SiPM_Groups = new TTreeReaderValue<vector<vector<pair<Signal, Signal>>>>(*Reader, "MERGED_Tree_SiPMGroup");
    if (FLAG2025) HRS = new TTreeReaderValue<Signal>(*Reader, "MERGED_Tree_HRS");

    clock_t start = clock(), Current;
    int Entries = Reader->GetEntries();
    int limit = Entries;
    while (Reader->Next() && Reader->GetCurrentEntry() < limit)
    {
        ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "Reading Tree");
        if (FLAG2025)
        {
            if (IsHRSProton((**HRS).Label))
                continue;
        }

        int Strip_Label = (*Silicon)[1].Label;
        double energy = Calibration[Strip_Label]->Eval((*Silicon)[1].Channel / 1000);

        // cout << "Energy: " << energy << "    " << Strip_Label << endl;
       

        if (WindowsMap[NUCLEUS][IAS[NUCLEUS]][Strip_Label].first > energy || energy > WindowsMap[NUCLEUS][IAS[NUCLEUS]][Strip_Label].second) // only for ias
        {
            continue;
        }
        H_Single[Strip_Label]->Fill(energy);
        double NEAREST = 1e9;
        int NEAREST_GROUP_INDEX = -1;
        vector<double> nearest_group_time;
        vector<double> mean_group_time = vector<double>(10, 0);
        if ((**SiPM_Groups).size() > 0)
        {
            // cout << "Size: " << (**SiPM_Groups).size() << endl;
            // Lopping on subgroups
            for (int i_group = 0; i_group < (**SiPM_Groups).size(); i_group++)
            {
                nearest_group_time.push_back(1e9);
                int counter_mean_group_time = 0;
                // if ((**SiPM_Groups)[i_group].size() > 9)
                // {
                //     for (int i_pair = 0; i_pair < (**SiPM_Groups)[i_group].size(); i_pair++)
                //     {
                //         cout << (**SiPM_Groups)[i_group][i_pair].first << endl;
                //     }
                // }
                // Looping on pair of the subgroup
                for (int i_pair = 0; i_pair < (**SiPM_Groups)[i_group].size(); i_pair++)
                {

                    // Taking valid element of pair with High priority
                    if ((**SiPM_Groups)[i_group][i_pair].first.isValid && !isnan((**SiPM_Groups)[i_group][i_pair].first.Time))
                    {
                        if (abs(nearest_group_time[i_group]) > abs((**SiPM_Groups)[i_group][i_pair].first.Time))
                        {
                            nearest_group_time[i_group] = (**SiPM_Groups)[i_group][i_pair].first.Time;
                        }

                        mean_group_time[i_group] += (**SiPM_Groups)[i_group][i_pair].first.Time;
                        counter_mean_group_time++;
                    }
                    else if ((**SiPM_Groups)[i_group][i_pair].second.isValid && !isnan((**SiPM_Groups)[i_group][i_pair].second.Time))
                    {
                        if (abs(nearest_group_time[i_group]) > abs((**SiPM_Groups)[i_group][i_pair].second.Time))
                        {
                            nearest_group_time[i_group] = (**SiPM_Groups)[i_group][i_pair].second.Time;
                        }
                        mean_group_time[i_group] += (**SiPM_Groups)[i_group][i_pair].first.Time;
                        counter_mean_group_time++;
                    }
                }
                mean_group_time[i_group] /= counter_mean_group_time;
                // cout << "Mean Time: " << mean_group_time[i_group] << endl;
            }

            // Saving the nearest group in time and index group
            for (int i_group = 0; i_group < nearest_group_time.size(); i_group++)
            {
                if (abs(NEAREST) > abs(nearest_group_time[i_group]))
                {
                    NEAREST_GROUP_INDEX = i_group;
                    NEAREST = nearest_group_time[i_group];
                }

                for (int mul = 1; mul <= BETA_SIZE; mul++)
                {
                    if (mul <= (**SiPM_Groups)[i_group].size())
                    {
                        H_SiPM_Time[Strip_Label][mul]->Fill(nearest_group_time[i_group]);
                        if (mul == 9)
                        {
                            H_Time_Channel[Strip_Label]->Fill(mean_group_time[i_group], energy);
                        }
                    }
                }
            }
        }

        // cout << "Nearest: " << NEAREST << "    " << NEAREST_GROUP_INDEX << endl;

        int Multiplicity;
        if (NEAREST_GROUP_INDEX == -1)
            Multiplicity = 0;
        else
        {
            Multiplicity = (**SiPM_Groups)[NEAREST_GROUP_INDEX].size();
        }

        // Filling Hist with repect to Multiplicity
        for (int mul = 1; mul <= BETA_SIZE; mul++)
        {
            if (mul <= Multiplicity                              // Mulitplicity sorting condition
                && (NEAREST > start_gate && NEAREST < end_gate)) // Time gate condition
            {

                H_Coinc[Strip_Label][mul]->Fill(energy);
                H_SiPM_Time_Coinc[Strip_Label][mul]->Fill(NEAREST);
                
                if (mul == Multiplicity) // for equal Multiplicity
                    H_Coinc_Mulitplicity[Strip_Label][mul]->Fill(energy);

                for (int i = 0; i < Multiplicity; i++)
                {
                    if (!(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.isValid)
                        continue;
                    H_SiPM_Channel_M_coinc[(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.Label][mul]->Fill((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.Channel);
                }
                for (int i = 0; i < Multiplicity; i++)
                {
                    if (!(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.isValid)
                        continue;
                    H_SiPM_Channel_M_coinc[(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.Label][mul]->Fill((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.Channel);
                }
            }
            else
            {
                H_NoCoinc[Strip_Label][mul]->Fill(energy);

                // for (int i = 0; i < Multiplicity; i++)
                // {
                //     if (!(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.isValid)
                //         continue;
                //     H_SiPM_Channel_M_nocoinc[(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.Label][mul]->Fill((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.Channel);
                // }
                // for (int i = 0; i < Multiplicity; i++)
                // {
                //     if (!(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.isValid)
                //         continue;
                //     H_SiPM_Channel_M_nocoinc[(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.Label][mul]->Fill((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.Channel);
                // }
            }
        }

        for (int i_group = 0; i_group < (**SiPM_Groups).size(); i_group++)
        {
            for (int mul = 1; mul <= BETA_SIZE; mul++)
            {
                if (mul <= (**SiPM_Groups)[i_group].size()                                                                         // Mulitplicity sorting condition
                    && ((**SiPM_Groups)[i_group][0].first.Time < start_gate || (**SiPM_Groups)[i_group][0].first.Time > end_gate)) // Time gate condition
                {
                    for (int i = 0; i < (**SiPM_Groups)[i_group].size(); i++)
                    {
                        if ((**SiPM_Groups)[i_group][i].first.isValid)
                            H_SiPM_Channel_M_nocoinc[(**SiPM_Groups)[i_group][i].first.Label][mul]->Fill((**SiPM_Groups)[i_group][i].first.Channel);
                        if ((**SiPM_Groups)[i_group][i].second.isValid)
                            H_SiPM_Channel_M_nocoinc[(**SiPM_Groups)[i_group][i].second.Label][mul]->Fill((**SiPM_Groups)[i_group][i].second.Channel);
                    }
                }
            }
        }
    }
}

void RandomCorrection()
{
    Info("Random Coincidence Correction");
    ///////// COMPUTE COINCIDENCE CORRECTION //////////
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            for (int mul = 1; mul <= BETA_SIZE; mul++)
            {
                H_SiPM_Time[i][mul]->GetXaxis()->SetRangeUser(start_gate_fake, end_gate_fake);
                HRandom[i][mul] = H_SiPM_Time[i][mul]->Integral() / (abs(start_gate_fake - end_gate_fake));
                H_SiPM_Time[i][mul]->GetXaxis()->SetRangeUser(-1111, -1111);
                NRandom[i][mul] = HRandom[i][mul] * abs(end_gate - start_gate);
            }
        }
    }

    ///////// APPLY CORRECTION OF RANDOM COINCIDENCES //////////

    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            for (int mul = 1; mul <= BETA_SIZE; mul++)
            {
                H_NoCoinc_Corrected[i][mul] = (TH1D *)H_NoCoinc[i][mul]->Clone(("H_NoCoinc_Corrected" + detectorName[i] + "_" + to_string(mul)).c_str());
                H_NoCoinc_Corrected[i][mul]->SetTitle(("H_NoCoinc_Corrected" + detectorName[i] + "_" + to_string(mul)).c_str());
                H_Coinc_Corrected[i][mul] = (TH1D *)H_Coinc[i][mul]->Clone(("H_Coinc_Corrected" + detectorName[i] + "_" + to_string(mul)).c_str());
                H_Coinc_Corrected[i][mul]->SetTitle(("H_Coinc_Corrected" + detectorName[i] + "_" + to_string(mul)).c_str());

                // compute
                H_NoCoinc[i][mul]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][IAS[NUCLEUS]][i].first, WindowsMap[NUCLEUS][IAS[NUCLEUS]][i].second);
                double integral = H_NoCoinc[i][mul]->Integral();
                H_NoCoinc[i][mul]->Scale(1. / integral);
                H_NoCoinc[i][mul]->GetXaxis()->SetRangeUser(-1111, -1111);

                // apply
                H_Coinc_Corrected[i][mul]->Add(H_NoCoinc[i][mul], -NRandom[i][mul]);
                H_NoCoinc_Corrected[i][mul]->Add(H_NoCoinc[i][mul], NRandom[i][mul]);

                H_Random[i][mul] = (TH1D *)H_NoCoinc[i][mul]->Clone(("H_Random" + detectorName[i] + "_" + to_string(mul)).c_str());
                H_Random[i][mul]->Scale(NRandom[i][mul]);
                H_NoCoinc[i][mul]->Scale(integral);

                // for display
                G_NRandom[mul]->AddPoint(i, NRandom[i][mul]);
            }
        }
    }
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
    double Eshift_error_2B = Nc/Ns * pow(Eshift, 2) / Ns * ( 1 - Nc/Ns );
    double Eshift_error_2 = sqrt(Eshift_error_2A + Eshift_error_2B);


    cout << "Eshift: " << Eshift << endl;

    cout << "Eshift_error_1: " << Eshift_error_1 << endl;
    cout << "Eshift_error_2: " << Eshift_error_2 << endl;

    cout << "Eshift_error_1A: " << Eshift_error_1A << endl;
    cout << "Eshift_error_2A: " << Eshift_error_2A << endl;

    cout << "Eshift_error_1B: " << Eshift_error_1B << endl;
    cout << "Eshift_error_2B: " << Eshift_error_2B << endl;
        
    
    return make_pair(Eshift, Eshift_error_2);
}