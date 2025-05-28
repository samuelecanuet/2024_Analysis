#include "../../Grouper/include/Detectors.hh"


int Verbosee = 0;

map<string, pair<double, double>[100][SIGNAL_MAX]> WindowsMap;

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

double Chi2(const TGraphErrors* g1, const TGraphErrors* g2) {
    if (!g1 || !g2) {
        Error("One of the graphs is null.");
        return -1;
    }

    int n1 = g1->GetN();
    int n2 = g2->GetN();

    if (n1 != n2) {
        Error("Graphs have different number of points.");
        return -1;
    }

    double chi2 = 0.0;
    for (int i = 0; i < n1; ++i) {
        double x1, y1, x2, y2;
        g1->GetPoint(i, x1, y1);
        g2->GetPoint(i, x2, y2);

        if (std::abs(x1 - x2) > 1e-6) {
            Error("x-values at point " + std::to_string(i) + " differ: " + std::to_string(x1) + " vs " + std::to_string(x2));   
            continue;
        }

        double err1 = g1->GetErrorY(i);
        double err2 = g2->GetErrorY(i);
        double err_tot2 = err1 * err1 + err2 * err2;

        if (err_tot2 <= 0) {
            Error("Total error squared is zero or negative at point " + std::to_string(i));
            continue;
        }

        double delta = y1 - y2;
        double contrib = (delta * delta) / err_tot2;
        chi2 += contrib;
    }

    return chi2;
}

