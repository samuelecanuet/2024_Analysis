#include "../../Grouper/include/Detectors.hh"

vector<pair<double, TFile *>> FILES;

TH1D* H[20][SIGNAL_MAX];
TCanvas *c[SIGNAL_MAX];
TLegend *leg[SIGNAL_MAX];   
map<string, pair<double, double>[100][SIGNAL_MAX]> WindowsMap;

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

void InitWindows()
{
    string direction[2] = {"Up", "Down"};
    for (auto dir : direction)
    {
        for (int strip = 1; strip <= 5; strip++)
        {
            ifstream file("../Grouper/Config_Files/Detector_Window/" + dir + "_" + to_string(strip) + ".txt");
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

                for (int i : Dir2Det(dir, strip))
                {
                    WindowsMap[nuclei][number][i] = make_pair(energy_low, energy_high);
                    // cout << "Nuclei : " << nuclei << " Number : " << number << " Detector : " << detectorName[i] << " Energy Low : " << energy_low << " Energy High : " << energy_high << endl;
                }
            }
        }
    }
}