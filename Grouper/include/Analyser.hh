#include "Detectors.hh"

TFile *ANALYSED_File;
TFile *MERGED_File;
TFile *CALIBRATED_File;


map<string, pair<double, double>[100][SIGNAL_MAX]> WindowsMap;

TF1* Calibration[SIGNAL_MAX];
TH1D* H_NoCoinc[SIGNAL_MAX][MAX_MULTIPLICTY];
TH1D* H_Coinc[SIGNAL_MAX][MAX_MULTIPLICTY];
TH1D* H_Coinc_Mulitplicity[SIGNAL_MAX][MAX_MULTIPLICTY];
TH1D* H_Coinc_Mulitplicity_Corrected[SIGNAL_MAX][MAX_MULTIPLICTY];
TH1D* H_Single[SIGNAL_MAX];
TH1D* H_NoCoinc_Corrected[SIGNAL_MAX][MAX_MULTIPLICTY];
TH1D* H_Coinc_Corrected[SIGNAL_MAX][MAX_MULTIPLICTY];
TH1D* H_SiPM_Time[SIGNAL_MAX][MAX_MULTIPLICTY];
TH1D* H_SiPM_Time_Coinc[SIGNAL_MAX][MAX_MULTIPLICTY];
TGraph* G_NFake[MAX_MULTIPLICTY];
TDirectory *dir_FakeCorrection;
TDirectory *dir_FakeCorrection_Strip[SIGNAL_MAX];   
TDirectory *dir_FakeCorrection_Strip_Write[SIGNAL_MAX];
TDirectory *dir_FakeCorrection_Strip_Time[SIGNAL_MAX];
double NFake[SIGNAL_MAX][MAX_MULTIPLICTY];
double HFake[SIGNAL_MAX][MAX_MULTIPLICTY];


/// params
double start_gate = -20;
double end_gate = 50;
double start_gate_fake = -300;
double end_gate_fake = -240;

void InitHistograms()
{
    dir_FakeCorrection = ANALYSED_File->mkdir("FakeCorrection");
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            dir_FakeCorrection_Strip[i] = dir_FakeCorrection->mkdir(detectorName[i].c_str());
            dir_FakeCorrection_Strip_Write[i] = dir_FakeCorrection_Strip[i]->mkdir("Write");


            H_Single[i] = new TH1D(("H_Single_" + detectorName[i]).c_str(), ("H_Single_" + detectorName[i]).c_str(), 50000, 0, 10000);
            H_Single[i]->GetXaxis()->SetTitle("Energy [keV]");
            H_Single[i]->GetYaxis()->SetTitle("Counts");
            H_Single[i]->GetXaxis()->CenterTitle();
            H_Single[i]->GetYaxis()->CenterTitle();

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
    }

    for (int mul = 1; mul <= BETA_SIZE; mul++)
    {
        G_NFake[mul] = new TGraph();
    }

    
}

void InitWindows()
{
    string direction[2] = {"Up", "Down"};
    for (auto dir : direction)
    {
        for (int strip = 1; strip <= 5; strip++)
        {
            ifstream file("Config_Files/Detector_Window/" + dir + "_" + to_string(strip) + ".txt");
            if (!file.is_open())
            {
                Error("Impossible to open " + dir + "_" + to_string(strip) + ".txt");
            }

            string line;
            double energy_low;
            double energy_high;
            int number;
            string nuclei;
            while (getline(file, line))
            {
                energy_high = -1;
                energy_low = -1;

                if (line.empty())
                {
                    continue;
                }

                if (line.find("#") != string::npos)
                {
                    nuclei = line.substr(1);
                    continue;
                }
                stringstream ss(line);
                ss >> number >> energy_low >> energy_high;

                // init silicon detector window
                for (int i : Dir2Det(dir, strip))
                {
                    WindowsMap[nuclei][number][i] = make_pair(energy_low, energy_high);
                }
            }
        }
    }
}


void InitCalib()
{
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            Calibration[i] = (TF1*)CALIBRATED_File->Get(("Calibration_" + detectorName[i]).c_str());

            if (Calibration[i] == NULL)
            {
                Calibration[i] = new TF1(("Calibration_" + detectorName[i]).c_str(), "x", 0, 10000);
            }
        }
    }
}
