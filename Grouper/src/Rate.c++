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
    InitDetectors("Config_Files/sample.pid");
    InitHistograms();
    Init();

    // Reference time for the Year

    //// ****SCALER PLOT**** ////
    TFile *file;
    if (FLAG2021)
        DIR_ROOT_DATA += "2021/";
    if (FLAG2021)
            DIR_ROOT_DATA_ANALYSED += "2021/";
    string filename;
    for (string run : Runs[year])
    {
        filename = SearchFiles(DIR_ROOT_DATA, run, "QW");
        if (filename != "")
        {
           break;
        }
    }
    file = MyTFile(DIR_ROOT_DATA + filename, "READ", "Q");
    // getting of the run
    pair<string, string> time = GetTime(file);
    Time_Year_Reference = Convert_DatetoTimeSec(time.first, FLAG2021);
    Info("Time reference of the year: " + to_string(Time_Year_Reference) + " s");

    int counter_run = -1;
    int counter_time[SIGNAL_MAX] = {0};
    for (string run : Runs[year])
    {
        
        Info("Processing run " + run);

        //// ****SCALER PLOT**** ////
        TFile *file;

        string filename = SearchFiles(DIR_ROOT_DATA, run);
        if (filename == "")
        {
            continue;
        }
        counter_run ++;
        file = MyTFile(DIR_ROOT_DATA + filename, "READ");

        // getting of the run
        pair<string, string> time = GetTime(file);

        if (time.second == " ")
        {
            Warning("No end time found for run " + run);
            string runafter = Runs[year][counter_run+1];
            string filenameafter = SearchFiles(DIR_ROOT_DATA, runafter);
            TFile *fileafter = new TFile((DIR_ROOT_DATA + filenameafter).c_str(), "READ");
            pair<string, string> timeafter = GetTime(fileafter);
            time.second = timeafter.first;
            fileafter->Close();           
        }

        double Time_Run = Diff_Date(time.second, time.first, FLAG2021);
        double Time_Run_Reference = Convert_DatetoTimeSec(time.first, FLAG2021);     

        

        Info("Run time: " + to_string(Time_Run) + " s", 1);
        Info("Run time reference: " + to_string(Time_Run_Reference) + " s", 1);
        Info("Run time start: " + time.first, 1);
        Info("Run time end: " + time.second, 1);


        // Silicon Scalers
        Info("Processing Silicon Scalers/Rate");
        int counter=0;
        int TOTAL_Silicon_Rate = 0;
        int TOTAL_Silicon_Scaler = 0;
        for (int det = 0; det < SIGNAL_MAX; det++)
        {
            if (IsDetectorSiliBack(det))
            {
                if (FLAG2021)
                    detectorName[det] = "Si" + to_string(GetDetector(det)) + "_R";

                TH1D* H_TOTAL_SCALER = (TH1D*)file->Get((detectorName[det] + "/" + detectorName[det] + "_Scaler").c_str());   
                double scaler = H_TOTAL_SCALER->GetBinContent(1);
 
                G_Scaler[det]->SetPoint( counter_run, (Time_Run/2 + Time_Run_Reference - Time_Year_Reference)/3600, scaler);
                G_Scaler[det]->SetPointError( counter_run, Time_Run/2/3600, sqrt(scaler));

                TOTAL_Silicon_Scaler += scaler;

                TH1D* H_TOTAL_RATE = (TH1D*)file->Get((detectorName[det] + "/" + detectorName[det] + "_DiffScalerTrigTime").c_str());
                for (int bin = 0; bin < H_TOTAL_RATE->GetNbinsX(); bin++)
                {
                    G_Rate[det]->SetPoint( counter_time[det], (Time_Run_Reference+H_TOTAL_RATE->GetBinCenter(bin)-Time_Year_Reference)/3600, H_TOTAL_RATE->GetBinContent(bin));

                    counter_time[det]++;
                }
            }
        }

        G_Scaler_Silicons->SetPoint(counter_run, (Time_Run / 2 + Time_Run_Reference - Time_Year_Reference) / 3600, TOTAL_Silicon_Scaler);
        G_Scaler_Silicons->SetPointError(counter_run, Time_Run /2. / 3600, sqrt(TOTAL_Silicon_Scaler));

        G_MeanRate_Silicons->SetPoint(counter_run, (Time_Run / 2 + Time_Run_Reference - Time_Year_Reference) / 3600, TOTAL_Silicon_Scaler/Time_Run);
        G_MeanRate_Silicons->SetPointError(counter_run, Time_Run /2. / 3600, sqrt(TOTAL_Silicon_Scaler)/Time_Run);

        // adding all the G_rate[det] to the total
        for (int det = 0; det < SIGNAL_MAX; det++)
        {
            if (IsDetectorSiliBack(det))
            {
                for (int i = 0; i < counter_time[det]; i++)
                {
                    double x, y;
                    G_Rate[det]->GetPoint(i, x, y);
                    G_Rate_Silicons->SetPoint(i, x, G_Rate_Silicons->GetPointY(i) + y);
                }
            }
        }

        // SiPM Scalers
        counter=0;
        int TOTAL_SiPM_Rate = 0;
        int TOTAL_SiPM_Scaler = 0;
        Info("Processing SiPM Scalers/Rate");
        for (int det = 0; det < SIGNAL_MAX; det++)
        {
            if (IsDetectorBetaHigh(det))
            {
                TH1D* H_TOTAL_SCALER = (TH1D*)file->Get((detectorName[det] + "/" + detectorName[det] + "_Scaler").c_str());   
                double scaler = H_TOTAL_SCALER->GetBinContent(1);
 
                G_Scaler[det]->SetPoint( counter_run, (Time_Run/2 + Time_Run_Reference - Time_Year_Reference)/3600, scaler);
                G_Scaler[det]->SetPointError( counter_run, Time_Run/2/3600, sqrt(scaler));

                TOTAL_SiPM_Scaler += scaler;

                TH1D* H_TOTAL_RATE = (TH1D*)file->Get((detectorName[det] + "/" + detectorName[det] + "_DiffScalerCalcTime").c_str());
                for (int bin = 0; bin < H_TOTAL_RATE->GetNbinsX(); bin++)
                {
                    G_Rate[det]->SetPoint( counter_time[det], (Time_Run_Reference+H_TOTAL_RATE->GetBinCenter(bin)-Time_Year_Reference)/3600, H_TOTAL_RATE->GetBinContent(bin));

                    counter_time[det]++;
                }
            }
        }

        G_Scaler_SiPMs->SetPoint(counter_run, (Time_Run / 2 + Time_Run_Reference - Time_Year_Reference) / 3600, TOTAL_SiPM_Scaler);
        G_Scaler_SiPMs->SetPointError(counter_run, Time_Run /2. / 3600, sqrt(TOTAL_SiPM_Scaler));

        G_MeanRate_SiPMs->SetPoint(counter_run, (Time_Run / 2 + Time_Run_Reference - Time_Year_Reference) / 3600, TOTAL_SiPM_Scaler/Time_Run);
        G_MeanRate_SiPMs->SetPointError(counter_run, Time_Run /2. / 3600, sqrt(TOTAL_SiPM_Scaler)/Time_Run);

        // adding all the G_rate[det] to the total
        for (int det = 0; det < SIGNAL_MAX; det++)
        {
            if (IsDetectorBetaHigh(det))
            {
                for (int i = 0; i < counter_time[det]; i++)
                {
                    double x, y;
                    G_Rate[det]->GetPoint(i, x, y);
                    G_Rate_SiPMs->SetPoint(i, x, G_Rate_SiPMs->GetPointY(i) + y);
                }
            }
        }

        file->Close();


        // Proportion
        Info("Processing Proportion");
        
        filename = SearchFiles(DIR_ROOT_DATA_ANALYSED, run);
        file = MyTFile(DIR_ROOT_DATA_ANALYSED + filename, "READ");

        // getting TGraph drawn in a TCanvas in file
        TCanvas *c = (TCanvas*)file->Get("RatioCoinc_NoCoinc_3");
        //loop on primitive of the canvas to get the TGraph
        for (int i = 0; i < c->GetListOfPrimitives()->GetSize(); i++)
        {
            if (c->GetListOfPrimitives()->At(i)->InheritsFrom("TGraph"))
            {
                TGraph *g = (TGraph*)c->GetListOfPrimitives()->At(i);
                g->Fit("pol0", "RQ");
                double p0;
                if (!FLAG2021)
                    p0 = g->GetFunction("pol0")->GetParameter(0);
                else
                {
                    // get numbe rof the poinbt with x = 41
                    int point = 0;
                    for (int i = 0; i < g->GetN(); i++)
                    {
                        double x, y;
                        g->GetPoint(i, x, y);
                        if (x == 45)
                        {
                            point = i;
                            break;
                        }
                    }
                    p0 = g->GetY()[point];
                }
                G_MeanProportion->SetPoint(counter_run, (Time_Run / 2 + Time_Run_Reference - Time_Year_Reference) / 3600, p0);
                G_MeanProportion->SetPointError(counter_run, Time_Run /2. / 3600, 0);
            }
        }

        
    }

    FINAL_file->cd();
    WriteHistograms();
    FINAL_file->Close();
}