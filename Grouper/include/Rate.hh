
#include "Detectors.hh"

map<int , vector<string>> Runs;
double Time_Year_Reference;


TGraphErrors *G_Scaler[SIGNAL_MAX];
TGraphErrors *G_Rate[SIGNAL_MAX];
TGraphErrors *G_Scaler_Silicons;


void Init()
{
    Runs[2021] = {"Run1", "Run2", "Run3"};
    Runs[2024] = {
                           "057", "058", "059",
                         "061", "062", 
                         "064", "065", "066", "067", "068", "069",
                         "070", "071", "072", "074", "077",
                        
                        
                        
                        "112", "113", "114", "115", "116", "118"
                        };
}   

void InitHistograms()
{
    Info("Init Histograms start");
    G_Scaler_Silicons = new TGraphErrors();
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliBack(i))
        {
            G_Scaler[i] = new TGraphErrors();
            G_Rate[i] = new TGraphErrors();
        }
    }
    Info("Init Histograms end");
}

void WriteHistograms()
{
    Info("Write Histograms start");
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliBack(i))
        {
            G_Scaler[i]->Write(("G_Scaler_" + detectorName[i]).c_str());
            G_Rate[i]->Write(("G_Rate_" + detectorName[i]).c_str());
        }
    }
    G_Scaler_Silicons->Write("G_Scaler_Silicons");
    Info("Write Histograms end");
}