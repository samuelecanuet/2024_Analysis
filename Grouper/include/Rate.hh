
#include "Detectors.hh"

map<int , vector<string>> Runs;
double Time_Year_Reference;

//all
TGraphErrors *G_Scaler[SIGNAL_MAX];
TGraphErrors *G_Rate[SIGNAL_MAX];

//silicon detectors
TGraphErrors *G_Scaler_Silicons;
TGraphErrors *G_Rate_Silicons;
TGraphErrors *G_MeanRate_Silicons;

//sipm
TGraphErrors *G_Scaler_SiPMs;
TGraphErrors *G_Rate_SiPMs;
TGraphErrors *G_MeanRate_SiPMs;

//physics
TGraphErrors *G_IAS[SIGNAL_MAX];
TGraphErrors *G_IAS_Silicons;

//proportion
TGraphErrors *G_Proportion[SIGNAL_MAX];
TGraphErrors *G_MeanProportion;



void Init()
{
    Runs[2021] = {"001", "002", "003", "004", "005", "006", "007", "008", "009", 
                    "010", "011", "012", "013", "014", "015", "016", "017", "018", "019",
                    "020", "021", "022", "023", "024", "025", "026", "027", "028", "029",
                    "030", "031", "032", "033", "034", "035", "036", "037", "038", "039",
                    "040", "041", "042", "043", "044", "045", "046", "047", "048", "049",
                    "050", "051", "052", "053", "054", "055", "056", "057", "058", "059",
                    "060", "061", "062", "063", "064", "065", "066", "067", "068", "069",
    };

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
    G_Rate_Silicons = new TGraphErrors();
    G_Scaler_SiPMs = new TGraphErrors();
    G_Rate_SiPMs = new TGraphErrors();
    G_MeanRate_Silicons = new TGraphErrors();
    G_MeanRate_SiPMs = new TGraphErrors();
    G_MeanProportion = new TGraphErrors();
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            G_Proportion[i] = new TGraphErrors();
        }
        if (IsDetectorSiliBack(i))
        {
            G_Scaler[i] = new TGraphErrors();
            G_Rate[i] = new TGraphErrors();
        }

        if (IsDetectorBetaHigh  (i))
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
    TCanvas *c_Silicons_Rate = new TCanvas("c_Silicons_Rate", "c_Silicons_Rate", 1920, 1080);
    TCanvas *c_Silicons_Scaler = new TCanvas("c_Silicons_Scaler", "c_Silicons_Scaler", 1920, 1080);
    TCanvas *c_SiPMs_Rate = new TCanvas("c_SiPMs_Rate", "c_SiPMs_Rate", 1920, 1080);
    TCanvas *c_SiPMs_Scaler = new TCanvas("c_SiPMs_Scaler", "c_SiPMs_Scaler", 1920, 1080);
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliBack(i))
        {
            G_Scaler[i]->Write(("G_Scaler_" + detectorName[i]).c_str());
            G_Rate[i]->Write(("G_Rate_" + detectorName[i]).c_str());
        }
        if (IsDetectorBetaHigh(i))
        {
            G_Scaler[i]->Write(("G_Scaler_" + detectorName[i]).c_str());
            G_Rate[i]->Write(("G_Rate_" + detectorName[i]).c_str());
        }
    }
    G_Scaler_Silicons->Write("G_Scaler_Silicons");
    G_Rate_Silicons->Write("G_Rate_Silicons");
    G_Scaler_SiPMs->Write("G_Scaler_SiPMs");
    G_Rate_SiPMs->Write("G_Rate_SiPMs");
    G_MeanRate_Silicons->Write("G_MeanRate_Silicons");
    G_MeanRate_SiPMs->Write("G_MeanRate_SiPMs");
    G_MeanProportion->Write("G_MeanProportion");
    Info("Write Histograms end");
}