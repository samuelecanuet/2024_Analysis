#include "Detectors.hh"

int VERBOSE = 0;

TFile *ANALYSED_File;
TFile *MERGED_File;
TFile *CALIBRATED_File;
TFile *CALIBRATED_SiPM_File;

TTreeReader *Reader;
TTreeReaderArray<Signal> *Silicon;
TTreeReaderValue<vector<vector<pair<Signal, Signal>>>> *SiPM_Groups;
TTreeReaderValue<Signal> *HRS;

string NUCLEUS = "32Ar";
double PEAK;

int MULTIPLICITY; // at least
double THRESHOLD; //keV

// Loading calibrations
TF1* Calibration[SIGNAL_MAX];
TF1* SiPM_Calibration[SIGNAL_MAX];

// -- raw -- //
// Proton spectrum
TH1D *H_Single[SIGNAL_MAX];
TH1D* H_NoCoinc[SIGNAL_MAX];
TH1D* H_Coinc[SIGNAL_MAX];
TH1D* H_Coinc_Mequal[SIGNAL_MAX];

// Time
TH1D *H_Time[SIGNAL_MAX];
TH2D* H_Time_Energy[SIGNAL_MAX];

// Beta spectrum
TH1D* H_SiPM_Energy_M_coinc[SIGNAL_MAX];
TH1D* H_SiPM_Energy_M_nocoinc[SIGNAL_MAX];
TH1D* H_SiPM_Energy_M_coinc_cal[SIGNAL_MAX];
TH1D* H_SiPM_Energy_M_nocoinc_cal[SIGNAL_MAX];

// Random
TGraphErrors *G_NRandom;
TH1D* H_Random[SIGNAL_MAX];
double NRandom[SIGNAL_MAX];

// -- corrected -- //
TH1D *H_Coinc_Corrected[SIGNAL_MAX];
TH1D *H_NoCoinc_Corrected[SIGNAL_MAX];


// Directories
TDirectory *dir_RandomCorrection;
TDirectory *dir_RandomCorrection_Strip[SIGNAL_MAX];
TDirectory *dir_RandomCorrection_Strip_Write[SIGNAL_MAX];
TDirectory *dir_RandomCorrection_Strip_Time[SIGNAL_MAX];

/// params
double start_gate = -30;
double end_gate = 50;
double start_gate_fake = -300;
double end_gate_fake = -100;

int Run=0;
string Run_string;

/////////////////////////////////
TF1* SpreadAcceptance[SIGNAL_MAX];
double MeanAcceptance[SIGNAL_MAX];
void LoadFitParameters()
{
  string add = "_2024";
  if (FLAG2021) add = "_2021";
  if (FLAG2025) add = "_2025";
  TFile *file = MyTFile((DIR_ROOT_DATA_GROUPED + "Grouping_FitParameters"+add+".root").c_str(), "READ");
  for (int i = 0; i < detectorNum; i++)
  {
    if (IsDetectorSiliStrip(i))
    {
      //cut
      SpreadAcceptance[i] = (TF1 *)file->Get(("SpreadAcceptance_" + detectorName[i]).c_str());
      MeanAcceptance[i] = ((TGraph *)file->Get("MeanAcceptance"))->GetY()[i];
      if (SpreadAcceptance[i] == nullptr)
      {
        Error(("SpreadAcceptance_" + detectorName[i]).c_str());
      }
      if (MeanAcceptance[i] ==  0.)
      {
        Warning(("MeanAcceptance_" + detectorName[i]).c_str());
      }
    }
  }
  file->Close();
  Success("Fit parameters loaded");
}

bool IsAccepted(double RearChannel, double Strip_Channel, double sigma, int det)
{
    double diff = (RearChannel - Strip_Channel) / Strip_Channel;
    double spread = SpreadAcceptance[det]->Eval(Strip_Channel); 
    double mean = MeanAcceptance[det];
    if (diff > mean + sigma*spread || diff < mean - sigma*spread)
    {
        return true;
    }
    else
    {
        return false;
    }
}


////////////////////////////////

void InitCalib(int sigma_calib = 0)
{
    CALIBRATED_File = MyTFile((DIR_ROOT_DATA_CALIBRATED + "Calibrated_" + to_string(YEAR) + ".root").c_str(), "READ");
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            Calibration[i] = (TF1 *)CALIBRATED_File->Get(("Calibration_" + detectorName[i]).c_str());
            

            if (Calibration[i] == NULL)
            {
                Error("No calibration found for " + detectorName[i]);
            }

            // error progagtion
            Calibration[i]->SetParameter(1, Calibration[i]->GetParameter(1) + sigma_calib * Calibration[i]->GetParError(1));
            Calibration[i]->SetParameter(0, Calibration[i]->GetParameter(0) + sigma_calib * Calibration[i]->GetParError(0));            
        }
    }

    Success("Silicon Calibration loaded");


    CALIBRATED_SiPM_File = MyTFile((DIR_ROOT_DATA_CALIBRATED + "SiPM_Calibrated_" + to_string(YEAR) + ".root").c_str(), "READ");
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorBetaHigh(i))
        {
            SiPM_Calibration[GetDetectorChannel(i)] = (TF1 *)CALIBRATED_SiPM_File->Get(("SiPM_" + to_string(GetDetectorChannel(i)) + "/F_Calibration_SiPM" + to_string(GetDetectorChannel(i))).c_str());

            if (SiPM_Calibration[GetDetectorChannel(i)] == NULL)
            {
                Error("No SiPM calibration found for " + detectorName[i]);
            }
        }
    }

    Success("SiPM Calibration loaded");
}


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

            // Proton Spectrum
            // Single
            H_Single[i] = new TH1D(("H_Single_" + detectorName[i]).c_str(), ("H_Single_" + detectorName[i]).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
            H_Single[i]->GetXaxis()->SetTitle("Energy (keV)");
            H_Single[i]->GetYaxis()->SetTitle("Counts");
            H_Single[i]->GetXaxis()->CenterTitle();
            H_Single[i]->GetYaxis()->CenterTitle();

            // Coinc at least
            H_Coinc[i] = new TH1D(("H_Coinc_" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str(), ("H_Coinc_" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
            H_Coinc[i]->GetXaxis()->SetTitle("Energy (keV)");
            H_Coinc[i]->GetYaxis()->SetTitle("Counts");
            H_Coinc[i]->GetXaxis()->CenterTitle();
            H_Coinc[i]->GetYaxis()->CenterTitle();

            // No coinc at least
            H_NoCoinc[i] = new TH1D(("H_NoCoinc_" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str(), ("H_NoCoinc_" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
            H_NoCoinc[i]->GetXaxis()->SetTitle("Energy (keV)");
            H_NoCoinc[i]->GetYaxis()->SetTitle("Counts");
            H_NoCoinc[i]->GetXaxis()->CenterTitle();
            H_NoCoinc[i]->GetYaxis()->CenterTitle();

            // Coinc equal
            H_Coinc_Mequal[i] = new TH1D(("H_Coinc_Mequal_" + detectorName[i] + "_M=" + to_string(MULTIPLICITY)).c_str(), ("H_Coinc_Mequal_" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
            H_Coinc_Mequal[i]->GetXaxis()->SetTitle("Energy (keV)");
            H_Coinc_Mequal[i]->GetYaxis()->SetTitle("Counts");
            H_Coinc_Mequal[i]->GetXaxis()->CenterTitle();
            H_Coinc_Mequal[i]->GetYaxis()->CenterTitle();

            H_Time[i] = new TH1D(("H_Time_" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str(), ("H_Time_" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str(), maxi*10, mini, maxi);
            H_Time[i]->GetXaxis()->SetTitle("Time (ns)");
            H_Time[i]->GetYaxis()->SetTitle("Counts");
            H_Time[i]->GetXaxis()->CenterTitle();
            H_Time[i]->GetYaxis()->CenterTitle();

            H_Time_Energy[i] = new TH2D(("H_Time_Energy_" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str(), ("H_Time_Energy_" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str(), maxi, mini, maxi, (WindowsMap[NUCLEUS][PEAK][i].second-WindowsMap[NUCLEUS][PEAK][i].first)*10, WindowsMap[NUCLEUS][PEAK][i].first, WindowsMap[NUCLEUS][PEAK][i].second);
            H_Time_Energy[i]->GetXaxis()->SetTitle("Time (ns)");
            H_Time_Energy[i]->GetYaxis()->SetTitle("Energy (keV)");
            H_Time_Energy[i]->GetXaxis()->CenterTitle();
            H_Time_Energy[i]->GetYaxis()->CenterTitle();

            H_Random[i] = new TH1D(("H_Random_" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str(), ("H_Random_" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str(), eSiliN_cal, eSiliMin_cal, eSiliMax_cal);
            H_Random[i]->GetXaxis()->SetTitle("Energy (keV)");
            H_Random[i]->GetYaxis()->SetTitle("Counts");
            H_Random[i]->GetXaxis()->CenterTitle();
            H_Random[i]->GetYaxis()->CenterTitle();         
        }

        if (IsDetectorBetaHigh(i))
        {

                H_SiPM_Energy_M_nocoinc[i] = new TH1D(("H_SiPM_Energy_nocoinc_" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str(), ("H_SiPM_Energy_nocoinc_" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str(), eHighN, eHighMin, eHighMax);
                H_SiPM_Energy_M_nocoinc[i]->GetXaxis()->SetTitle("Channel");
                H_SiPM_Energy_M_nocoinc[i]->GetYaxis()->SetTitle("Counts");
                H_SiPM_Energy_M_nocoinc[i]->GetXaxis()->CenterTitle();
                H_SiPM_Energy_M_nocoinc[i]->GetYaxis()->CenterTitle();

                H_SiPM_Energy_M_coinc[i] = new TH1D(("H_SiPM_Energy_coinc_" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str(), ("H_SiPM_Energy_coinc_" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str(), eHighN, eHighMin, eHighMax);
                H_SiPM_Energy_M_coinc[i]->GetXaxis()->SetTitle("Channel");
                H_SiPM_Energy_M_coinc[i]->GetYaxis()->SetTitle("Counts");
                H_SiPM_Energy_M_coinc[i]->GetXaxis()->CenterTitle();
                H_SiPM_Energy_M_coinc[i]->GetYaxis()->CenterTitle();

                //calibrated 
                H_SiPM_Energy_M_nocoinc_cal[i] = new TH1D(("H_SiPM_Energy_nocoinc_" + detectorName[i] + "_M" + to_string(MULTIPLICITY) + "_cal").c_str(), ("H_SiPM_Energy_nocoinc_" + detectorName[i] + "_M" + to_string(MULTIPLICITY) + "_cal").c_str(), eSiliN_cal/10, eSiliMin_cal, eSiliMax_cal);
                H_SiPM_Energy_M_nocoinc_cal[i]->GetXaxis()->SetTitle("Energy (keV)");
                H_SiPM_Energy_M_nocoinc_cal[i]->GetYaxis()->SetTitle("Counts");
                H_SiPM_Energy_M_nocoinc_cal[i]->GetXaxis()->CenterTitle();
                H_SiPM_Energy_M_nocoinc_cal[i]->GetYaxis()->CenterTitle();

                H_SiPM_Energy_M_coinc_cal[i] = new TH1D(("H_SiPM_Energy_coinc_" + detectorName[i] + "_M" + to_string(MULTIPLICITY) + "_cal").c_str(), ("H_SiPM_Energy_coinc_" + detectorName[i] + "_M" + to_string(MULTIPLICITY) + "_cal").c_str(), eSiliN_cal/10, eSiliMin_cal, eSiliMax_cal);
                H_SiPM_Energy_M_coinc_cal[i]->GetXaxis()->SetTitle("Energy (keV)");
                H_SiPM_Energy_M_coinc_cal[i]->GetYaxis()->SetTitle("Counts");
                H_SiPM_Energy_M_coinc_cal[i]->GetXaxis()->CenterTitle();
                H_SiPM_Energy_M_coinc_cal[i]->GetYaxis()->CenterTitle();
            
        }

        if (IsDetectorBetaLow(i))
        {

                H_SiPM_Energy_M_nocoinc[i] = new TH1D(("H_SiPM_Energy_nocoinc_" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str(), ("H_SiPM_Energy_nocoinc_" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str(), eLowN, eLowMin, eLowMax);
                H_SiPM_Energy_M_nocoinc[i]->GetXaxis()->SetTitle("Channel");
                H_SiPM_Energy_M_nocoinc[i]->GetYaxis()->SetTitle("Counts");
                H_SiPM_Energy_M_nocoinc[i]->GetXaxis()->CenterTitle();
                H_SiPM_Energy_M_nocoinc[i]->GetYaxis()->CenterTitle();

                H_SiPM_Energy_M_coinc[i] = new TH1D(("H_SiPM_Energy_coinc_" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str(), ("H_SiPM_Energy_coinc_" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str(), eLowN, eLowMin, eLowMax);
                H_SiPM_Energy_M_coinc[i]->GetXaxis()->SetTitle("Channel");
                H_SiPM_Energy_M_coinc[i]->GetYaxis()->SetTitle("Counts");
                H_SiPM_Energy_M_coinc[i]->GetXaxis()->CenterTitle();
                H_SiPM_Energy_M_coinc[i]->GetYaxis()->CenterTitle();

                //calibrated
                H_SiPM_Energy_M_nocoinc_cal[i] = new TH1D(("H_SiPM_Energy_nocoinc_" + detectorName[i] + "_M" + to_string(MULTIPLICITY) + "_cal").c_str(), ("H_SiPM_Energy_nocoinc_" + detectorName[i] + "_M" + to_string(MULTIPLICITY) + "_cal").c_str(), eSiliN_cal/10, eSiliMin_cal, eSiliMax_cal);
                H_SiPM_Energy_M_nocoinc_cal[i]->GetXaxis()->SetTitle("Energy (keV)");
                H_SiPM_Energy_M_nocoinc_cal[i]->GetYaxis()->SetTitle("Counts");
                H_SiPM_Energy_M_nocoinc_cal[i]->GetXaxis()->CenterTitle();
                H_SiPM_Energy_M_nocoinc_cal[i]->GetYaxis()->CenterTitle();  

                H_SiPM_Energy_M_coinc_cal[i] = new TH1D(("H_SiPM_Energy_coinc_" + detectorName[i] + "_M" + to_string(MULTIPLICITY) + "_cal").c_str(), ("H_SiPM_Energy_coinc_" + detectorName[i] + "_M" + to_string(MULTIPLICITY) + "_cal").c_str(), eSiliN_cal/10, eSiliMin_cal, eSiliMax_cal);
                H_SiPM_Energy_M_coinc_cal[i]->GetXaxis()->SetTitle("Energy (keV)");
                H_SiPM_Energy_M_coinc_cal[i]->GetYaxis()->SetTitle("Counts");
                H_SiPM_Energy_M_coinc_cal[i]->GetXaxis()->CenterTitle();
                H_SiPM_Energy_M_coinc_cal[i]->GetYaxis()->CenterTitle();
            
        }
    }

    G_NRandom = new TGraphErrors();

    Info("Init Histograms end");
}

void ReaderData()
{
    Info("Reading data");
    Reader = new TTreeReader((TTree *)MERGED_File->Get("MERGED_Tree"));
    Silicon = new TTreeReaderArray<Signal>(*Reader, "MERGED_Tree_Silicon");
    SiPM_Groups = new TTreeReaderValue<vector<vector<pair<Signal, Signal>>>>(*Reader, "MERGED_Tree_SiPMGroup");
    HRS = new TTreeReaderValue<Signal>(*Reader, "MERGED_Tree_HRS");

    clock_t start = clock(), Current;
    int Entries = Reader->GetEntries();
    int limit = Entries;
    while (Reader->Next() && Reader->GetCurrentEntry() < limit)
    {
        ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "Reading Tree");

        if (IsHRSProton((**HRS).Label))
            continue;

        int Strip_Label = (*Silicon)[1].Label;
        double Strip_Energy = Calibration[Strip_Label]->Eval((*Silicon)[1].Channel / 1000.);

        // if (IsAccepted((*Silicon)[0].Channel, (*Silicon)[1].Channel, 1, Strip_Label))
        // {
        //     continue;
        // }

        // Proton selection
        if (WindowsMap[NUCLEUS][PEAK][Strip_Label].first > Strip_Energy || Strip_Energy > WindowsMap[NUCLEUS][PEAK][Strip_Label].second)
            continue;
        
        H_Single[Strip_Label]->Fill(Strip_Energy);

        // ------------------------------------------------ //
        // IF 2025 REMOVING SIPM6
        if (YEAR == 2025 && (**SiPM_Groups).size() > 0)
        {
            for (int i_group = 0; i_group < (**SiPM_Groups).size(); i_group++)
            {
                // where is SiPM6 in the group
                vector<int> index_delete = {};
                for (int i = 0; i < (**SiPM_Groups)[i_group].size(); i++)
                {
                    if (GetDetectorChannel((**SiPM_Groups)[i_group][i].first.Label) == 6 || GetDetectorChannel((**SiPM_Groups)[i_group][i].second.Label) == 6)
                    {
                        index_delete.push_back(i);
                    }
                }

                // removing SiPM6 from the group
                for (int i = index_delete.size() - 1; i >= 0; i--)
                {
                    (**SiPM_Groups)[i_group].erase((**SiPM_Groups)[i_group].begin() + index_delete[i]);
                }
            }           
        }
        // ------------------------------------------------ //

        // ------------------------------------------------ //
        // GETTING MULTIPLCITY FOR EACH SUBGROUP
        for (int i_group = (**SiPM_Groups).size() - 1; i_group >= 0; i_group--)
        {
            // calculating the multiplicity of the subgroup
            int Multiplicity = 0;
            for (int i = 0; i < (**SiPM_Groups)[i_group].size(); i++)
            {
                double e_sipm = 0;
                if ((**SiPM_Groups)[i_group][i].first.isValid)
                {
                    e_sipm = SiPM_Calibration[GetDetectorChannel((**SiPM_Groups)[i_group][i].first.Label)]->Eval((**SiPM_Groups)[i_group][i].first.Channel / 1000.);
                }
                // else if ((**SiPM_Groups)[i_group][i].second.isValid)
                // {
                //     e_sipm = SiPM_Calibration[GetDetectorChannel((**SiPM_Groups)[i_group][i].second.Label)]->Eval((**SiPM_Groups)[i_group][i].second.Channel / 1000.);
                // }

                if (e_sipm > THRESHOLD) // 100keV
                {
                    Multiplicity++;
                }
            }   

            // removing the subgroup if multiplicity < M3
            if (Multiplicity < MULTIPLICITY)
            {
                (**SiPM_Groups).erase((**SiPM_Groups).begin() + i_group);
            }
        }

        // ------------------------------------------------ //
        // GETTING THE NEAREST SUBGROUP IF MORE THAN 1
        double Nearest_Mean_Time_Diff = -1111;
        int NEAREST_GROUP_INDEX = -1;
        if ((**SiPM_Groups).size() > 0)
        {

            // MEAN TIME SUBGROUP
            vector<double> Mean_Time_Diff_SubGroup = vector<double>((**SiPM_Groups).size(), -1111);
            for (int i_group = 0; i_group < (**SiPM_Groups).size(); i_group++)
            {
                double Mean_Time_SiPM = 0;
                int N_SiPM = 0;
                for (int i = 0; i < (**SiPM_Groups)[i_group].size(); i++)
                {
                    if ((**SiPM_Groups)[i_group][i].first.isValid && !isnan((**SiPM_Groups)[i_group][i].first.Time))
                    {
                        Mean_Time_SiPM += (**SiPM_Groups)[i_group][i].first.Time;
                        N_SiPM++;
                    }
                    // else if ((**SiPM_Groups)[i_group][i].second.isValid && !isnan((**SiPM_Groups)[i_group][i].second.Time))
                    // {
                    //     Mean_Time_SiPM += (**SiPM_Groups)[i_group][i].second.Time;
                    //     N_SiPM++;
                    // }
                }
                if (N_SiPM > 0)
                {
                    Mean_Time_SiPM /= N_SiPM;
                    Mean_Time_Diff_SubGroup[i_group] = Mean_Time_SiPM;
                }
            }

            // find nearest from 0
            // NEAREST_GROUP_INDEX = distance(Mean_Time_Diff_SubGroup.begin(), min_element(Mean_Time_Diff_SubGroup.begin(), Mean_Time_Diff_SubGroup.end(), [](double a, double b) { return abs(a) < abs(b); }));
            // own min loop
            double min_time = 1111;
            for (int i = 0; i < Mean_Time_Diff_SubGroup.size(); i++)
            {
                if (abs(Mean_Time_Diff_SubGroup[i]) < abs(min_time))
                {
                    min_time = Mean_Time_Diff_SubGroup[i];
                    NEAREST_GROUP_INDEX = i;
                }
            }
            Nearest_Mean_Time_Diff = Mean_Time_Diff_SubGroup[NEAREST_GROUP_INDEX];

            // NEAREST TIME SUBGROUP
            // vector<double> Nearest_Time_Diff_SubGroup = vector<double>((**SiPM_Groups).size(), -1111);
            // for (int i_group = 0; i_group < (**SiPM_Groups).size(); i_group++)
            // {
            //     double Nearest_Time_Diff_SiPM = 1111;
            //     for (int i = 0; i < (**SiPM_Groups)[i_group].size(); i++)
            //     {
            //         if ((**SiPM_Groups)[i_group][i].first.isValid && !isnan((**SiPM_Groups)[i_group][i].first.Time))
            //         {
            //             if (abs((**SiPM_Groups)[i_group][i].first.Time) < abs(Nearest_Time_Diff_SiPM))
            //             {
            //                 Nearest_Time_Diff_SiPM = (**SiPM_Groups)[i_group][i].first.Time;
            //             }
            //         }
            //         else if ((**SiPM_Groups)[i_group][i].second.isValid && !isnan((**SiPM_Groups)[i_group][i].second.Time))
            //         {
            //             if (abs((**SiPM_Groups)[i_group][i].second.Time) < abs(Nearest_Time_Diff_SiPM))
            //             {
            //                 Nearest_Time_Diff_SiPM = (**SiPM_Groups)[i_group][i].second.Time;
            //             }
            //         }
            //     }
            //     Nearest_Time_Diff_SubGroup[i_group] = Nearest_Time_Diff_SiPM;
            // }
            // // find nearest from 0
            // NEAREST_GROUP_INDEX = distance(Nearest_Time_Diff_SubGroup.begin(), min_element(Nearest_Time_Diff_SubGroup.begin(), Nearest_Time_Diff_SubGroup.end(), [](double a, double b) { return abs(a) < abs(b); }));
            // Nearest_Mean_Time_Diff = Nearest_Time_Diff_SubGroup[NEAREST_GROUP_INDEX];
        }
        // ------------------------------------------------ //


        // ------------------------------------------------ //
        // ADD THE CONDITION ON THE TIME GATE //
        if ((Nearest_Mean_Time_Diff > start_gate && Nearest_Mean_Time_Diff < end_gate) && Nearest_Mean_Time_Diff != -1111) 
        {
            int Multiplicity = 0;
            
            // Filling at least multiplicity
            H_Coinc[Strip_Label]->Fill(Strip_Energy);           
            H_Time_Energy[Strip_Label]->Fill(Nearest_Mean_Time_Diff, Strip_Energy);

            // Beta Energy
            if ((**SiPM_Groups).size() != 0)
            {
                for (int i = 0; i < (**SiPM_Groups)[NEAREST_GROUP_INDEX].size(); i++)
                {
                    double e_sipm = 0;
                    if ((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.isValid)
                    {
                        e_sipm = SiPM_Calibration[GetDetectorChannel((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.Label)]->Eval((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.Channel / 1000.);
                        H_SiPM_Energy_M_coinc[(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.Label]->Fill((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.Channel);
                        H_SiPM_Energy_M_coinc_cal[(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.Label]->Fill(e_sipm);
                    }
                    else if ((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.isValid)
                    {
                        e_sipm = SiPM_Calibration[GetDetectorChannel((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.Label)]->Eval((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.Channel / 1000.);
                        H_SiPM_Energy_M_coinc[(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.Label]->Fill((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.Channel);
                        H_SiPM_Energy_M_coinc_cal[(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.Label]->Fill(e_sipm);
                    }

                    if (e_sipm > THRESHOLD) // 100keV
                    {
                        Multiplicity++;
                    }
                }
            }

            // Filling equal multiplicity
            if (Multiplicity == MULTIPLICITY)
                H_Coinc_Mequal[Strip_Label]->Fill(Strip_Energy);
        }
        else 
        {
            H_NoCoinc[Strip_Label]->Fill(Strip_Energy);
            // Beta Energy
            if ((**SiPM_Groups).size() != 0)
            {
                for (int i = 0; i < (**SiPM_Groups)[NEAREST_GROUP_INDEX].size(); i++)
                {
                    double e_sipm = 0;
                    if ((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.isValid)
                    {
                        e_sipm = SiPM_Calibration[GetDetectorChannel((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.Label)]->Eval((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.Channel / 1000.);

                        H_SiPM_Energy_M_nocoinc[(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.Label]->Fill((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.Channel);
                        H_SiPM_Energy_M_nocoinc_cal[(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].first.Label]->Fill(e_sipm);
                    }
                    else if ((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.isValid)
                    {
                        e_sipm = SiPM_Calibration[GetDetectorChannel((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.Label)]->Eval((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.Channel / 1000.);

                        H_SiPM_Energy_M_nocoinc[(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.Label]->Fill((**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.Channel);
                        H_SiPM_Energy_M_nocoinc_cal[(**SiPM_Groups)[NEAREST_GROUP_INDEX][i].second.Label]->Fill(e_sipm);
                    }                
                }
            }

            // Filling Random
            // if (Nearest_Mean_Time_Diff > start_gate_fake && Nearest_Mean_Time_Diff < end_gate_fake)
            // {
            //     H_Random[Strip_Label]->Fill(Strip_Energy); // not used
            // }
        }
        H_Time[Strip_Label]->Fill(Nearest_Mean_Time_Diff);
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
            H_Time[i]->GetXaxis()->SetRangeUser(start_gate_fake, end_gate_fake);
            double scaling = H_Time[i]->Integral() / (abs(start_gate_fake - end_gate_fake));
            H_Time[i]->GetXaxis()->SetRangeUser(-1111, -1111);
            NRandom[i] = scaling * abs(end_gate - start_gate);
        }
    }
    ///////// APPLY CORRECTION OF RANDOM COINCIDENCES //////////
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
                H_NoCoinc_Corrected[i] = (TH1D *)H_NoCoinc[i]->Clone(("H_NoCoinc_Corrected" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str());
                H_NoCoinc_Corrected[i]->SetTitle(("H_NoCoinc_Corrected" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str());
                H_Coinc_Corrected[i] = (TH1D *)H_Coinc[i]->Clone(("H_Coinc_Corrected" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str());
                H_Coinc_Corrected[i]->SetTitle(("H_Coinc_Corrected" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str());

                // compute
                H_NoCoinc[i]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][PEAK][i].first, WindowsMap[NUCLEUS][PEAK][i].second);
                double integral = H_NoCoinc[i]->Integral();
                H_NoCoinc[i]->Scale(1. / integral);
                H_NoCoinc[i]->GetXaxis()->SetRangeUser(-1111, -1111);

                // apply
                H_Coinc_Corrected[i]->Add(H_NoCoinc[i], -NRandom[i]);
                H_NoCoinc_Corrected[i]->Add(H_NoCoinc[i], NRandom[i]);

                H_Random[i] = (TH1D *)H_NoCoinc[i]->Clone(("H_Random" + detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str());
                H_Random[i]->Scale(NRandom[i]);
                H_NoCoinc[i]->Scale(integral);

                // for display
                G_NRandom->AddPoint(i, NRandom[i]);
            
        }
    }
}

pair<double, double> ComputeEshift(int det, TH1D *H_Single, TH1D *H_Coinc, TH1D *H_NoCoinc)
{
    H_Single->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][PEAK][det].first, WindowsMap[NUCLEUS][PEAK][det].second);
    H_Coinc->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][PEAK][det].first, WindowsMap[NUCLEUS][PEAK][det].second);
    H_NoCoinc->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][PEAK][det].first, WindowsMap[NUCLEUS][PEAK][det].second);

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
        
    
    return make_pair(Eshift, Eshift_error_1);
}

void WriteHistograms()
{
    Info("Writing Histograms");

    ANALYSED_File->cd();
    TDirectory *dir_Eshift = ANALYSED_File->mkdir("Eshift");
    TDirectory *dir_Eshift_m = dir_Eshift->mkdir(Form("Eshift_M%d", MULTIPLICITY));
    TGraphErrors *G_Eshift = new TGraphErrors();
    TGraphErrors *G_Ratio_Coinc_Single = new TGraphErrors();
    TGraphErrors *G_Ratio_Coinc_Single_CORRECTED = new TGraphErrors();
    //high 
    TCanvas *c_SiPMHigh_Energy_Coinc = new TCanvas(("SiPMHigh_Energy_M" + to_string(MULTIPLICITY) + "_Coinc").c_str(), ("SiPMHigh_Energy_M" + to_string(MULTIPLICITY) + "_Coinc").c_str(), 1920, 1080);
    TCanvas *c_SiPMHigh_Channel_Coinc = new TCanvas(("SiPMHigh_Channel_M" + to_string(MULTIPLICITY) + "_Coinc").c_str(), ("SiPMHigh_Channel_M" + to_string(MULTIPLICITY) + "_Coinc").c_str(), 1920, 1080);
    TCanvas *c_SiPMHigh_Energy_Nocoinc = new TCanvas(("SiPMHigh_Energy_M" + to_string(MULTIPLICITY) + "_NoCoinc").c_str(), ("SiPMHigh_Energy_M" + to_string(MULTIPLICITY) + "_NoCoinc").c_str(), 1920, 1080);
    TCanvas *c_SiPMHigh_Channel_Nocoinc = new TCanvas(("SiPMHigh_Channel_M" + to_string(MULTIPLICITY) + "_NoCoinc").c_str(), ("SiPMHigh_Channel_M" + to_string(MULTIPLICITY) + "_NoCoinc").c_str(), 1920, 1080);
    //low
    TCanvas *c_SiPMLow_Energy_Coinc = new TCanvas(("SiPMLow_Energy_M" + to_string(MULTIPLICITY) + "_Coinc").c_str(), ("SiPMLow_Energy_M" + to_string(MULTIPLICITY) + "_Coinc").c_str(), 1920, 1080);
    TCanvas *c_SiPMLow_Channel_Coinc = new TCanvas(("SiPMLow_Channel_M" + to_string(MULTIPLICITY) + "_Coinc").c_str(), ("SiPMLow_Channel_M" + to_string(MULTIPLICITY) + "_Coinc").c_str(), 1920, 1080);
    TCanvas *c_SiPMLow_Energy_Nocoinc = new TCanvas(("SiPMLow_Energy_M" + to_string(MULTIPLICITY) + "_NoCoinc").c_str(), ("SiPMLow_Energy_M" + to_string(MULTIPLICITY) + "_NoCoinc").c_str(), 1920, 1080);
    TCanvas *c_SiPMLow_Channel_Nocoinc = new TCanvas(("SiPMLow_Channel_M" + to_string(MULTIPLICITY) + "_NoCoinc").c_str(), ("SiPMLow_Channel_M" + to_string(MULTIPLICITY) + "_NoCoinc").c_str(), 1920, 1080);

    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            if (H_Single[i]->GetEntries() == 0)
                continue;
            // cout << endl;
            if (VERBOSE == 1) Info(detectorName[i]);

            // Info("Multiplicity: " + to_string(mul), 1);
            dir_RandomCorrection_Strip[i]->cd();
            TCanvas *c = new TCanvas((detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str(), (detectorName[i] + "_M" + to_string(MULTIPLICITY)).c_str(), 1920, 1080);
            c->Divide(3, 1);
            c->cd(1);
            H_Single[i]->SetLineColor(kBlack);
            H_Single[i]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][PEAK][i].first, WindowsMap[NUCLEUS][PEAK][i].second);
            H_Single[i]->Draw("HIST");
            H_Coinc[i]->SetLineColor(kRed);
            H_Coinc[i]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][PEAK][i].first, WindowsMap[NUCLEUS][PEAK][i].second);
            H_Coinc[i]->Draw("HIST SAME");
            H_NoCoinc[i]->SetLineColor(kBlue);
            H_NoCoinc[i]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][PEAK][i].first, WindowsMap[NUCLEUS][PEAK][i].second);
            H_NoCoinc[i]->Draw("HIST SAME");
            H_Random[i]->SetLineColor(kGreen);
            H_Random[i]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][PEAK][i].first, WindowsMap[NUCLEUS][PEAK][i].second);
            H_Random[i]->Draw("HIST SAME");

            c->cd(2);
            H_Single[i]->SetLineColor(kBlack);
            H_Single[i]->Draw("HIST");
            H_Coinc_Corrected[i]->SetLineColor(kRed);
            H_Coinc_Corrected[i]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][PEAK][i].first, WindowsMap[NUCLEUS][PEAK][i].second);
            H_Coinc_Corrected[i]->Draw("HIST SAME");
            H_NoCoinc_Corrected[i]->SetLineColor(kBlue);
            H_NoCoinc_Corrected[i]->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][PEAK][i].first, WindowsMap[NUCLEUS][PEAK][i].second);
            H_NoCoinc_Corrected[i]->Draw("HIST SAME");

            c->cd(3);
            H_Time[i]->Draw("HIST");
            TLine *line_gate_start = new TLine(start_gate, 0, start_gate, H_Time[i]->GetMaximum());
            line_gate_start->SetLineColor(kRed);
            line_gate_start->SetLineStyle(2);
            line_gate_start->Draw("SAME");
            TLine *line_gate_end = new TLine(end_gate, 0, end_gate, H_Time[i]->GetMaximum());
            line_gate_end->SetLineColor(kRed);
            line_gate_end->SetLineStyle(2);
            line_gate_end->Draw("SAME");
            TLine *line_gate_fake_start = new TLine(start_gate_fake, 0, start_gate_fake, H_Time[i]->GetMaximum());
            line_gate_fake_start->SetLineColor(kGreen);
            line_gate_fake_start->SetLineStyle(2);
            line_gate_fake_start->Draw("SAME");
            TLine *line_gate_fake_end = new TLine(end_gate_fake, 0, end_gate_fake, H_Time[i]->GetMaximum());
            line_gate_fake_end->SetLineColor(kGreen);
            line_gate_fake_end->SetLineStyle(2);
            line_gate_fake_end->Draw("SAME");
            c->Write(); 
            
            
            // --- Ratio Coinc / Single --- //
            // RAW
            double ratio = H_Coinc[i]->Integral() / H_Single[i]->Integral();
            double ratio_error = sqrt( pow(sqrt(H_Coinc[i]->Integral()) / H_Single[i]->Integral(), 2) + pow(H_Coinc[i]->Integral() * sqrt(H_Single[i]->Integral()) / pow(H_Single[i]->Integral(), 2), 2) );
            G_Ratio_Coinc_Single->AddPoint(i, ratio);
            G_Ratio_Coinc_Single->SetPointError(G_Ratio_Coinc_Single->GetN()-1, 0, ratio_error);

            Info("Detector: " + detectorName[i] + ", Ratio Coinc/Single: " + to_string(ratio) + " +/- " + to_string(ratio_error), 1);
            // CORRECTED
            double ratio_corrected = H_Coinc_Corrected[i]->Integral() / H_Single[i]->Integral();
            double ratio_corrected_error = sqrt( pow(sqrt(H_Coinc_Corrected[i]->Integral()) / H_Single[i]->Integral(), 2) + pow(H_Coinc_Corrected[i]->Integral() * sqrt(H_Single[i]->Integral()) / pow(H_Single[i]->Integral(), 2), 2) );
            // Info("Detector: " + detectorName[i] + ", Ratio Coinc_Corrected/Single: " + to_string(ratio_corrected) + " +/- " + to_string(ratio_corrected_error), 1);
            G_Ratio_Coinc_Single_CORRECTED->AddPoint(i, ratio_corrected);
            G_Ratio_Coinc_Single_CORRECTED->SetPointError(G_Ratio_Coinc_Single_CORRECTED->GetN()-1, 0, ratio_corrected_error);

            dir_RandomCorrection_Strip_Write[i]->cd();
            H_Single[i]->Write();
            H_Coinc[i]->Write();
            H_NoCoinc[i]->Write();
            H_Random[i]->Write();
            H_Coinc_Corrected[i]->Write();
            H_NoCoinc_Corrected[i]->Write();

            ///Eshift Calculation
            pair<double, double> Eshift_result = ComputeEshift(i, H_Single[i], H_Coinc_Corrected[i], H_NoCoinc_Corrected[i]);
            double Eshift = Eshift_result.first;
            double Eshift_error = Eshift_result.second;
            Info("Eshift: " + to_string(Eshift) + " +/- " + to_string(Eshift_error) + " keV");
            G_Eshift->AddPoint(i, Eshift);
            G_Eshift->SetPointError(G_Eshift->GetN()-1, 0, Eshift_error);
        }

        if (IsDetectorBetaHigh(i))
        {
            c_SiPMHigh_Energy_Coinc->cd();
            H_SiPM_Energy_M_coinc_cal[i]->Draw("HIST SAME");
            c_SiPMHigh_Channel_Coinc->cd();
            H_SiPM_Energy_M_coinc[i]->Draw("HIST SAME");
            c_SiPMHigh_Energy_Nocoinc->cd();
            H_SiPM_Energy_M_nocoinc_cal[i]->Draw("HIST SAME");
            c_SiPMHigh_Channel_Nocoinc->cd();
            H_SiPM_Energy_M_nocoinc[i]->Draw("HIST SAME");  
        }
         
        if (IsDetectorBetaLow(i))
        {
            c_SiPMLow_Energy_Coinc->cd();
            H_SiPM_Energy_M_coinc_cal[i]->Draw("HIST SAME");
            c_SiPMLow_Channel_Coinc->cd();
            H_SiPM_Energy_M_coinc[i]->Draw("HIST SAME");
            c_SiPMLow_Energy_Nocoinc->cd();
            H_SiPM_Energy_M_nocoinc_cal[i]->Draw("HIST SAME");
            c_SiPMLow_Channel_Nocoinc->cd();
            H_SiPM_Energy_M_nocoinc[i]->Draw("HIST SAME");
        }
    }

    dir_Eshift_m->cd();
    TCanvas *c_Eshift = new TCanvas(Form("Eshift_M%d", MULTIPLICITY), Form("Eshift_M%d", MULTIPLICITY), 1920, 1080);
    G_Eshift->GetXaxis()->SetTitle("Detector");
    G_Eshift->GetYaxis()->SetTitle("Eshift (keV)");
    G_Eshift->GetXaxis()->SetLimits(-1, SIGNAL_MAX);
    G_Eshift->SetMarkerStyle(20);
    G_Eshift->Draw("AP");
    c_Eshift->Write();

    ANALYSED_File->cd();
    // ratio
    TCanvas *c_Ratio_Coinc_Single = new TCanvas(Form("Ratio_Coinc_Single_M%d", MULTIPLICITY), Form("Ratio_Coinc_Single_M%d", MULTIPLICITY), 1920, 1080);
    TMultiGraph *mg_Ratio_Coinc_Single = new TMultiGraph();
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(G_Ratio_Coinc_Single, "Raw", "P");
    legend->AddEntry(G_Ratio_Coinc_Single_CORRECTED, "Corrected", "P");
    G_Ratio_Coinc_Single->SetMarkerStyle(20);
    G_Ratio_Coinc_Single->SetMarkerColor(kRed);
    G_Ratio_Coinc_Single->SetLineColor(kRed);
    G_Ratio_Coinc_Single_CORRECTED->SetMarkerStyle(20);
    G_Ratio_Coinc_Single_CORRECTED->SetMarkerColor(kBlue);
    G_Ratio_Coinc_Single_CORRECTED->SetLineColor(kBlue);
    mg_Ratio_Coinc_Single->Add(G_Ratio_Coinc_Single, "AP");
    mg_Ratio_Coinc_Single->Add(G_Ratio_Coinc_Single_CORRECTED, "AP");
    mg_Ratio_Coinc_Single->GetXaxis()->SetTitle("Detector");
    mg_Ratio_Coinc_Single->GetYaxis()->SetTitle("Ratio Coinc/Single");
    mg_Ratio_Coinc_Single->Draw("A");
    legend->Draw("SAME");
    c_Ratio_Coinc_Single->Write();

    // Beta
    c_SiPMHigh_Energy_Coinc->Write();
    c_SiPMHigh_Channel_Coinc->Write();
    c_SiPMHigh_Energy_Nocoinc->Write();
    c_SiPMHigh_Channel_Nocoinc->Write();
    c_SiPMLow_Energy_Coinc->Write();
    c_SiPMLow_Channel_Coinc->Write();
    c_SiPMLow_Energy_Nocoinc->Write();
    c_SiPMLow_Channel_Nocoinc->Write();    





    Info("Writing Histograms done");
}