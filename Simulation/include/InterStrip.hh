#include "../../Grouper/include/Detectors.hh"

vector<pair<double, TFile *>> FILES;

TH1D* H[20][SIGNAL_MAX];
TCanvas *c[SIGNAL_MAX];
TLegend *leg[SIGNAL_MAX];   

TF1 *Calibration_Function[SIGNAL_MAX];
void InitCalib()
{
    TFile *CALIBRATED_File = MyTFile(DIR_ROOT_DATA_CALIBRATED + "Calibrated (copy).root", "READ");
    for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorSiliStrip(i))
            {
                Calibration_Function[i] = (TF1 *)CALIBRATED_File->Get(("Calibration_" + detectorName[i]).c_str());

                if (Calibration_Function[i] == NULL)
                {
                    Calibration_Function[i] = new TF1(("Calibration_" + detectorName[i]).c_str(), "x", 0, 10000);
                    Warning("No calibration found for " + detectorName[i]);
                }
            }
        }
}

double Detector_Resolution[SIGNAL_MAX];

void InitElectronicResolution()
{
    ifstream file(("../Grouper/Config_Files/" + to_string(YEAR) + "/Electronic_Resolution_" + to_string(YEAR) + ".txt").c_str());
    if (!file.is_open())
    {
        Error("Impossible to open Electronic_Resolution_" + to_string(YEAR) + ".txt");
    }
    Info("Reading Electronic_Resolution_" + to_string(YEAR) + ".txt");

    string line;
    int code;
    double resolution;
    double err_resolution;
    int counter = 0;
    while (getline(file, line))
    {
        counter++;
        stringstream ss(line);
        ss >> code >> resolution >> err_resolution;

        Detector_Resolution[code] = resolution / 1000.;
    }

    file.close();
}