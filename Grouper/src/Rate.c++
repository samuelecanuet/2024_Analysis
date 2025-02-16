#include "Rate.hh"

int main(int argc, char *argv[])
{
    int year;
    if (argc == 1)
    {
        Error("No year given");
    }
    if (argc == 2)
    {
        if (string(argv[1]) == "2021")
        {
            FLAG2021 = true;
            year = 2021;
        }
        else if (string(argv[1]) == "2024")
        {
            FLAG2021 = false;
            year = 2024;
        }
        else
        {
            Error("Wrong year given");
        }
    }

    if (FLAG2021)
        DIR_ROOT_DATA_RATE += "2021/";
    
    TFile *FINAL_file = MyTFile(DIR_ROOT_DATA_RATE + "Rate.root", "RECREATE");
    Info("Year: " + to_string(year));   
    InitDetectors("Config_Files/sample.pid", year);
    InitHistograms();
    Init(); 
    
    int counter_run = -1;
    int counter_time[SIGNAL_MAX] = {0};
    for (string run : Runs[year])
    {
        counter_run ++;
        Info("Processing run " + run);

        //// ****SCALER PLOT**** ////
        TFile *file;
        if (FLAG2021)
            DIR_ROOT_DATA += "2021/";

        string filename = SearchFiles(DIR_ROOT_DATA, run);
        if (filename == "")
        {
            continue;
        }
        file = MyTFile(DIR_ROOT_DATA + filename, "READ");

        // getting of the run
        pair<string, string> time = GetTime(file);

        if (time.second == " ")
        {
            string runafter = Runs[year][counter_run+1];
            string filenameafter = SearchFiles(DIR_ROOT_DATA, runafter);
            TFile *fileafter = new TFile((DIR_ROOT_DATA + filenameafter).c_str(), "READ");
            pair<string, string> timeafter = GetTime(fileafter);
            time.second = timeafter.first;
            fileafter->Close();           
        }

        double Time_Run = Convert_DatetoTime(time.second, time.first);
        double Time_Run_Reference = Convert_DatetoTimeSec(time.first);     

        // Reference time for the Year
        if (run == Runs[year][0])
        {
            Time_Year_Reference = Convert_DatetoTimeSec(time.first);
        }

        Info("Run time: " + to_string(Time_Run) + " s", 1);
        Info("Run time start: " + time.first, 1);
        Info("Run time end: " + time.second, 1);


        // Silicon Scalers
        Info("Processing Silicon Scalers/Rate");
        int counter=0;
        int TOTAL_Silicon_Rate = 0;
        int TOTAL_Silicon_Scaler = 0;
        int TOTAL_SiPM_Rate = 0;
        int TOTAL_SiPM_Scaler = 0;
        for (int det = 0; det < SIGNAL_MAX; det++)
        {
            
            if (IsDetectorSiliBack(det))
            {
                if (FLAG2021)
                    detectorName[det] = "Si" + to_string(GetDetector(det)) + "_R";

                TH1D* H_TOTAL_SCALER = (TH1D*)file->Get((detectorName[det] + "/" + detectorName[det] + "_Scaler").c_str());   
                double scaler = H_TOTAL_SCALER->GetBinContent(1);
 
                G_Scaler[det]->SetPoint( counter_run, (Time_Run/2 + Time_Run_Reference)/3600, scaler);
                G_Scaler[det]->SetPointError( counter_run, Time_Run/2/3600, sqrt(scaler));

                TOTAL_Silicon_Scaler += scaler;

                TH1D* H_TOTAL_RATE = (TH1D*)file->Get((detectorName[det] + "/" + detectorName[det] + "_DiffScalerTrigTime").c_str());
                for (int bin = 0; bin < H_TOTAL_RATE->GetNbinsX(); bin++)
                {
                    G_Rate[det]->SetPoint( counter_time[det], H_TOTAL_RATE->GetBinCenter(bin), H_TOTAL_RATE->GetBinContent(bin));
                    
                    counter_time[det]++;
                }

            }
        }

        G_Scaler_Silicons->SetPoint( counter_run, (Time_Run/2 + Time_Run_Reference)/3600, TOTAL_Silicon_Scaler);
        G_Scaler_Silicons->SetPointError( counter_run, Time_Run/2/3600, sqrt(TOTAL_Silicon_Scaler));
        
        file->Close();
    }

    FINAL_file->cd();
    WriteHistograms();
}