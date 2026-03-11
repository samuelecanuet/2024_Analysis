#include "../../Grouper/include/Detectors.hh"
#include "../../Grouper/include/Utilities.hh"
#include "Math/Minimizer.h"

// TFile
TFile *Exp_FILE;
TFile *Exp_FILE_grouped;
TFile *OUTPUT_FILE;

TDirectory *DIR_SIM;

// TF1
TF1 *Silicon_Calibration[SIGNAL_MAX];

// TGraphErrors
TGraphErrors *G_Exp_ProtonCounts = new TGraphErrors();
TGraphErrors *G_Exp_Eshift;
TGraphErrors *G_Exp_CoincSingle;

// Maps
map<string, TGraphErrors *> G_Sim_Eshift;
map<string, TGraphErrors *> G_Sim_Eshift_Relative;
map<string, TGraphErrors *> G_Sim_CoincSingle;
map<string, TGraphErrors *> G_Sim_ProtonCounts;
map<string, TGraphErrors *> G_Sim_ProtonCounts_Relative;

map<string, double> Chi2_Eshift;
map<string, double> Chi2_CoincSingle;
map<string, double> Chi2_ProtonCounts;
map<string, double> Chi2_Total;

TH3D *H3_Eshift;
TH3D *H3_CoincSingle;
TH3D *H3_ProtonCounts;
TH3D *H3_Total;

TGraph2DErrors *G2_Eshift_YZ = new TGraph2DErrors();
TGraph2D *G2_CoincSingle_YZ = new TGraph2D();
TGraph2DErrors *G2_ProtonCounts_YZ = new TGraph2DErrors();
TGraph2D *G2_Total_YZ = new TGraph2D();

TGraph2DErrors *G2_Eshift_XY = new TGraph2DErrors();
TGraph2D *G2_CoincSingle_XY = new TGraph2D();
TGraph2DErrors *G2_ProtonCounts_XY = new TGraph2DErrors();
TGraph2D *G2_Total_XY = new TGraph2D();

TGraph2DErrors *G2_MeanDiff_CoincSingle_XY = new TGraph2DErrors();

// Variables
string NUCLEUS;
string TYPE;

double GetMin(TH2D *H)
{
    // find the minimum bin content

    double mini = -1111;
    for (int i = 1; i <= H->GetNbinsX(); i++)
    {
        for (int j = 1; j <= H->GetNbinsY(); j++)
        {
            double content = H->GetBinContent(i, j);
            if ((mini == -1111 || content < mini) && content != 0)
                mini = content;
        }
    }

    return mini;
}

int GetColor(double value, double min, double max)
{
    if (value < min)
        value = min;
    if (value > max)
        value = max;
    double ratio = (value - min) / (max - min);
    int r = static_cast<int>(255. * ratio);
    int g = static_cast<int>(255. * (1 - ratio));
    int b = 0;
    return TColor::GetColor(r, g, b);
}

double Chi2(const TGraphErrors *g1, const TGraphErrors *g2, string option = "")
{
    if (!g1 || !g2)
    {
        Error("One of the graphs is null.");
        return -1;
    }

    int n1 = g1->GetN();
    int n2 = g2->GetN();

    if (n1 != n2)
    {
        Error("Graphs have different number of points.");
        return -1;
    }

    double chi2 = 0.0;
    for (int i = 0; i < n1; ++i)
    {
        double x1, y1, x2, y2;
        g1->GetPoint(i, x1, y1);
        g2->GetPoint(i, x2, y2);

        if (std::abs(x1 - x2) > 1e-6)
        {
            Error("x-values at point " + std::to_string(i) + " differ: " + std::to_string(x1) + " vs " + std::to_string(x2));
            continue;
        }

        if (option == "Up" && x1 < 50)
            continue;
        if (option == "Down" && x1 > 50)
            continue;

        double err1 = g1->GetErrorY(i);
        double err2 = g2->GetErrorY(i);
        double err_tot2 = err1 * err1 + err2 * err2;

        if (err_tot2 <= 0)
        {
            Error("Total error squared is zero or negative at point " + std::to_string(i));
            continue;
        }

        double delta = y1 - y2;
        double contrib = (delta * delta) / err_tot2;
        chi2 += contrib;
    }

    // if nan warning
    if (std::isnan(chi2))
    {
        Warning("Chi2 calculation resulted in NaN.");
        return 1e6;
    }

    return chi2;
}

double PValue(const TGraphErrors *g1, const TGraphErrors *g2, string option = "")
{
    double chi2 = Chi2(g1, g2, option);
    if (chi2 < 0)
    {
        return -1; // Error case
    }

    int n = g1->GetN();
    int dof = n - 1; // degrees of freedom

    // Calculate p-value using the chi-squared distribution
    double p_value = TMath::Prob(chi2, dof);
    return p_value;
}

double Sum(const TGraphErrors *g, string option = "")
{
    if (!g)
    {
        Error("Graph is null.");
        return -1;
    }

    double sum = 0.0;
    for (int i = 0; i < g->GetN(); ++i)
    {
        if (option == "Up" && g->GetX()[i] < 50)
            continue;
        if (option == "Down" && g->GetX()[i] > 50)
            continue;
        double x, y;
        g->GetPoint(i, x, y);
        sum += std::abs(y);
    }
    return sum;
}

pair<double, double> MeanDistance(TGraphErrors *g)
{
    if (!g)
    {
        Error("Graph is null.");
        return {-1, -1};
    }

    double mean = 0.0;
    double mean_err = 0.0;
    double weight_sum = 0.0;

    for (int i = 0; i < g->GetN(); ++i)
    {
        double x, y;
        double y_err = g->GetErrorY(i);
        g->GetPoint(i, x, y);   
        
        double weight = 1./(y_err*y_err);
        mean += abs(weight*y);
        weight_sum += weight;
    }

    mean /= weight_sum;
    mean_err = 1./sqrt(weight_sum);

    return {mean, mean_err};
}

pair<double, double> MeanDiff(TGraphErrors *g1, TGraphErrors *g2)
{
    // Calculate the mean diff between two graphs, and the mean error on the diff
    if (!g1 || !g2)
    {
        Error("One of the graphs is null.");
        return {-1, -1};
    }

    int n1 = g1->GetN();
    int n2 = g2->GetN();
    if (n1 != n2)
    {
        Error("Graphs have different number of points.");
        return {-1, -1};
    }

    double mean_g1 = 0.0;
    double mean_g2 = 0.0;
    double mean_err_g1 = 0.0;
    double mean_err_g2 = 0.0;

    for (int i = 0; i < n1; ++i)
    {
        double x1, y1, x2, y2;
        g1->GetPoint(i, x1, y1);
        g2->GetPoint(i, x2, y2);

        if (std::abs(x1 - x2) > 1e-6)
        {
            Error("x-values at point " + std::to_string(i) + " differ: " + std::to_string(x1) + " vs " + std::to_string(x2));
            continue;
        }

        mean_g1 += y1;
        mean_g2 += y2;
        mean_err_g1 += g1->GetErrorY(i) * g1->GetErrorY(i);
        mean_err_g2 += g2->GetErrorY(i) * g2->GetErrorY(i);
    }

    mean_g1 /= n1;
    mean_g2 /= n2;
    mean_err_g1 = sqrt(mean_err_g1) / n1;
    mean_err_g2 = sqrt(mean_err_g2) / n2;

    double mean_diff = mean_g1 - mean_g2;
    double mean_diff_err = sqrt(mean_err_g1 * mean_err_g1 + mean_err_g2 * mean_err_g2);

    return {mean_diff, mean_diff_err};
}

void LoadCalibrations()
{
    TFile *CALIBRATED_File = MyTFile(DIR_ROOT_DATA_CALIBRATED + "Calibrated_" + to_string(YEAR) + ".root", "READ");
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            Silicon_Calibration[i] = (TF1 *)CALIBRATED_File->Get(("Calibration_" + detectorName[i]).c_str());

            if (Silicon_Calibration[i] == NULL)
            {
                Silicon_Calibration[i] = new TF1(("Calibration_" + detectorName[i]).c_str(), "x", 0, 10000);
                Warning("No calibration found for " + detectorName[i]);
            }

            Silicon_Calibration[i] = InvertFunction((TF1 *)CALIBRATED_File->Get(("Calibration_" + detectorName[i]).c_str()));
        }
    }
}

vector<string> FindAllSimulatedFiles(const string &directory, const string &year)
{
    vector<string> filenames;
    DIR *dir = opendir(directory.c_str());
    if (dir == nullptr)
    {
        Error("Unable to open directory: " + directory);
    }

    // analysed files
    struct dirent *entry;
    while ((entry = readdir(dir)) != nullptr)
    {
        string file_name = entry->d_name;
        if (file_name == "." || file_name == "..")
            continue;

        if (file_name.find("result.root") != string::npos && file_name.find(year) != string::npos)
        {
            TFile *file = new TFile((directory + file_name).c_str(), "READ");
            if (file != nullptr)
            {
                filenames.push_back(directory + file_name);
                file->Close();
            }
        }
    }
    closedir(dir);

    return filenames;
}

double FindMinimalStep(const vector<double> &values)
{
    if (values.size() < 2)
        return 0;

    vector<double> sorted_values = values;
    sort(sorted_values.begin(), sorted_values.end());
    sorted_values.erase(unique(sorted_values.begin(), sorted_values.end()), sorted_values.end());

    double min_step = sorted_values[1] - sorted_values[0];
    for (size_t i = 1; i < sorted_values.size(); i++)
    {
        double step = sorted_values[i] - sorted_values[i - 1];
        if (step < min_step)
            min_step = step;
    }
    return min_step;
}

string AnalyseSim(string filename, double &x, double &y, double &z)
{
    TFile *Sim = MyTFile(filename.c_str(), "READ", "Q");

    // extracting parametrs from filename with regex
    std::regex rgx;
    if (TYPE == "Catcher")
        rgx = std::regex("catcherx([-+]?[0-9]*\\.?[0-9]+)_catchery([-+]?[0-9]*\\.?[0-9]+)_catcherz([-+]?[0-9]*\\.?[0-9]+)");
    else if (TYPE == "Detectors")
        rgx = std::regex("catcherz([-+]?[0-9]*\\.?[0-9]+)_Det_x([-+]?[0-9]*\\.?[0-9]+)_y([-+]?[0-9]*\\.?[0-9]+)");
    else
    {
        Error("Unknown TYPE: " + TYPE);
        return "";
    }
    std::smatch match;
    if (std::regex_search(filename, match, rgx) && match.size() == 4)
    {
        string x_pos = TYPE == "Catcher" ? to_string(stod(match[1])) : to_string(stod(match[2]));
        string y_pos = TYPE == "Catcher" ? to_string(stod(match[2])) : to_string(stod(match[3]));
        string z_pos = TYPE == "Catcher" ? to_string(stod(match[3])) : to_string(stod(match[1]));
        if (stod(x_pos) == 0.25 && stod(y_pos) == 4.0 && stod(z_pos) == 1.15)
            return "";

        x = TYPE == "Catcher" ? stod(match[1]) : stod(match[2]);
        y = TYPE == "Catcher" ? stod(match[2]) : stod(match[3]);
        z = TYPE == "Catcher" ? stod(match[3]) : stod(match[1]);
    }

    string par_string = "x=" + to_string(x) + " y=" + to_string(y) + " z=" + to_string(z);

    // - Eshift
    TCanvas *c_Sim = (TCanvas *)Sim->Get("Eshift");
    if (c_Sim == nullptr)
    {
        Warning("Eshift canvas not found in simulated file: " + filename);
        return "";
    }
    for (auto key : *c_Sim->GetListOfPrimitives())
    {
        if (string(key->ClassName()).find("TGraphErrors") != string::npos)
        {
            G_Sim_Eshift[par_string] = (TGraphErrors *)key;
        }
    }
    // Info("Eshift graph loaded", 1);

    // - Coinc/Single
    TCanvas *cCoincSingle = (TCanvas *)Sim->Get("Ratio_Coinc_Single");
    if (cCoincSingle == nullptr)
    {
        Warning("Coinc/Single canvas not found in simulated file: " + filename);
        return "";
    }
    for (auto key : *cCoincSingle->GetListOfPrimitives())
    {
        if (string(key->ClassName()).find("TGraphErrors") != string::npos)
        {
            G_Sim_CoincSingle[par_string] = (TGraphErrors *)key;
        }
    }
    // Info("Coinc/Single graph loaded", 1);

    // - Proton counts
    G_Sim_ProtonCounts[par_string] = new TGraphErrors();

    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliStrip(det))
        {
            TCanvas *c_Sim = (TCanvas *)Sim->Get(("H_Single_" + detectorName[det]).c_str());
            if (c_Sim == nullptr)
            {
                Warning("H_Single canvas not found for " + detectorName[det] + " in simulated file: " + filename);
                continue;
            }
            TH1D *H_Exp = nullptr;
            for (auto key : *c_Sim->GetListOfPrimitives())
            {
                if (string(key->GetName()).find("H_Single") != string::npos)
                {
                    H_Exp = (TH1D *)key;
                }
            }
            if (H_Exp == nullptr)
            {
                Warning("H_Single histogram not found in simulated file: " + filename);
                continue;
            }
            H_Exp->GetXaxis()->SetRangeUser(WindowsMap[NUCLEUS][IAS[NUCLEUS]][det].first, WindowsMap[NUCLEUS][IAS[NUCLEUS]][det].second);
            double integral = H_Exp->Integral();
            G_Sim_ProtonCounts[par_string]->AddPoint(det, integral);
            G_Sim_ProtonCounts[par_string]->SetPointError(G_Sim_ProtonCounts[par_string]->GetN() - 1, 0, sqrt(integral));
            delete H_Exp;
            delete c_Sim;
        }
    }

    // normalizing to the experimental number of counts
    double exp_total_counts = 0;
    for (int i = 0; i < G_Exp_ProtonCounts->GetN(); i++)
    {
        exp_total_counts += G_Exp_ProtonCounts->GetY()[i];
    }
    double sim_total_counts = 0;
    for (int i = 0; i < G_Sim_ProtonCounts[par_string]->GetN(); i++)
    {
        sim_total_counts += G_Sim_ProtonCounts[par_string]->GetY()[i];
    }
    if (sim_total_counts > 0)
    {
        double scale_factor = exp_total_counts / sim_total_counts;
        for (int i = 0; i < G_Sim_ProtonCounts[par_string]->GetN(); i++)
        {
            double y = G_Sim_ProtonCounts[par_string]->GetY()[i] * scale_factor;
            double err = sqrt(pow(exp_total_counts / sim_total_counts * sqrt(G_Sim_ProtonCounts[par_string]->GetY()[i]), 2) +
                              pow(G_Sim_ProtonCounts[par_string]->GetY()[i] / sim_total_counts * sqrt(exp_total_counts), 2) +
                              pow(G_Sim_ProtonCounts[par_string]->GetY()[i] * exp_total_counts / (sim_total_counts * sim_total_counts) * sqrt(sim_total_counts), 2));
            G_Sim_ProtonCounts[par_string]->SetPoint(i, G_Sim_ProtonCounts[par_string]->GetX()[i], y);
            G_Sim_ProtonCounts[par_string]->SetPointError(i, 0, err);
        }
    }
    // Info("Proton counts graph loaded", 1);

    return par_string;
}

void PlottingSimulation(string par_string)
{
    OUTPUT_FILE->cd();

    if (DIR_SIM == nullptr)
        DIR_SIM = OUTPUT_FILE->mkdir("Simulation");
    DIR_SIM->cd();

    TCanvas *C_Sim = new TCanvas(("C_Sim_" + par_string).c_str(), ("C_Sim_" + par_string).c_str(), 1200, 800);
    C_Sim->Divide(1, 3);

    // Proton counts
    C_Sim->cd(1);
    G_Sim_ProtonCounts[par_string]->SetTitle("Proton Counts;Detector;Counts");
    G_Sim_ProtonCounts[par_string]->SetMarkerStyle(20);
    G_Sim_ProtonCounts[par_string]->SetMarkerColor(kRed);
    G_Sim_ProtonCounts[par_string]->SetMarkerSize(2);
    G_Sim_ProtonCounts[par_string]->Draw("AP");

    // Coinc/Single
    C_Sim->cd(2);
    G_Sim_CoincSingle[par_string]->SetTitle("Coinc/Single Ratio;Code;Ratio");
    G_Sim_CoincSingle[par_string]->SetMarkerStyle(20);
    G_Sim_CoincSingle[par_string]->SetMarkerColor(kRed);
    G_Sim_CoincSingle[par_string]->SetMarkerSize(2);
    G_Sim_CoincSingle[par_string]->Draw("AP");

    // Eshift
    C_Sim->cd(3);
    G_Sim_Eshift[par_string]->SetTitle("Eshift;Code;Eshift (keV)");
    G_Sim_Eshift[par_string]->SetMarkerStyle(20);
    G_Sim_Eshift[par_string]->SetMarkerColor(kRed);
    G_Sim_Eshift[par_string]->SetLineColor(kRed);
    G_Sim_Eshift[par_string]->SetMarkerSize(2);
    G_Sim_Eshift[par_string]->Draw("AP");

    C_Sim->Write();
    OUTPUT_FILE->cd();
}

void PlottingExperimental()
{
    OUTPUT_FILE->cd();

    TCanvas *C_Exp = new TCanvas("C_Exp", "C_Exp", 1200, 800);
    C_Exp->Divide(1, 3);

    // Proton counts
    C_Exp->cd(1);
    G_Exp_ProtonCounts->SetTitle("Proton Counts;Detector;Counts");
    G_Exp_ProtonCounts->SetMarkerStyle(20);
    G_Exp_ProtonCounts->SetMarkerColor(kBlack);
    G_Exp_ProtonCounts->SetMarkerSize(2);
    G_Exp_ProtonCounts->Draw("AP");

    // Coinc/Single
    C_Exp->cd(2);
    G_Exp_CoincSingle->SetTitle("Coinc/Single Ratio;Code;Ratio");
    G_Exp_CoincSingle->SetMarkerStyle(20);
    G_Exp_CoincSingle->SetMarkerColor(kBlack);
    G_Exp_CoincSingle->SetMarkerSize(2);
    G_Exp_CoincSingle->Draw("AP");

    // Eshift
    C_Exp->cd(3);
    G_Exp_Eshift->SetTitle("Eshift;Code;Eshift (keV)");
    G_Exp_Eshift->SetMarkerStyle(20);
    G_Exp_Eshift->SetMarkerColor(kBlack);
    G_Exp_Eshift->SetLineColor(kBlack);
    G_Exp_Eshift->SetMarkerSize(2);
    G_Exp_Eshift->Draw("AP");

    C_Exp->Write();
}

void ComputeRelativeEshift(string par_string)
{
    double mean_sim_strip[SIGNAL_MAX] = {0};
    double mean_sim_strip_err[SIGNAL_MAX] = {0};
    double mean_exp_strip[SIGNAL_MAX] = {0};
    double mean_exp_strip_err[SIGNAL_MAX] = {0};
    for (int i = 0; i < G_Exp_Eshift->GetN(); i++)
    {

        double x = G_Exp_Eshift->GetX()[i];

        int strip_hemi = GetDetector(x) < 5 ? GetDetectorChannel(x) : GetDetectorChannel(x) + 10;

        mean_exp_strip[strip_hemi] += G_Exp_Eshift->GetY()[i];
        mean_exp_strip_err[strip_hemi] += pow(G_Exp_Eshift->GetErrorY(i), 2);
        mean_sim_strip[strip_hemi] += G_Sim_Eshift[par_string]->GetY()[i];
        mean_sim_strip_err[strip_hemi] += pow(G_Sim_Eshift[par_string]->GetErrorY(i), 2);
    }
    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        int strip_hemi = GetDetector(det) < 5 ? GetDetectorChannel(det) : GetDetectorChannel(det) + 10;
        mean_exp_strip[strip_hemi] /= (double)SILI_NUM / 2;
        mean_exp_strip_err[strip_hemi] = sqrt(mean_exp_strip_err[strip_hemi]) / (double)SILI_NUM / 2;
        mean_sim_strip[strip_hemi] /= (double)SILI_NUM / 2;
        mean_sim_strip_err[strip_hemi] = sqrt(mean_sim_strip_err[strip_hemi]) / (double)SILI_NUM / 2;
    }

    G_Sim_Eshift_Relative[par_string] = new TGraphErrors();
    for (int i = 0; i < G_Exp_Eshift->GetN(); i++)
    {
        double x = G_Exp_Eshift->GetX()[i];
        int strip_hemi = GetDetector(x) < 5 ? GetDetectorChannel(x) : GetDetectorChannel(x) + 10;

        double sim_y = G_Sim_Eshift[par_string]->GetY()[i];
        double sim_err = G_Sim_Eshift[par_string]->GetErrorY(i);
        double exp_y = G_Exp_Eshift->GetY()[i];
        double exp_err = G_Exp_Eshift->GetErrorY(i);
        double relative_y = (sim_y - exp_y - (mean_sim_strip[strip_hemi] - mean_exp_strip[strip_hemi]))/exp_y;
        double relative_err = sqrt(pow(sim_err/exp_y, 2) + pow(exp_err/exp_y, 2) + pow(mean_sim_strip_err[strip_hemi]/exp_y, 2) + pow(mean_exp_strip_err[strip_hemi]/exp_y, 2) + pow((sim_y - exp_y - (mean_sim_strip[strip_hemi] - mean_exp_strip[strip_hemi]))*exp_err/pow(exp_y, 2), 2));
        G_Sim_Eshift_Relative[par_string]->AddPoint(x, relative_y);
        G_Sim_Eshift_Relative[par_string]->SetPointError(G_Sim_Eshift_Relative[par_string]->GetN() - 1, 0, relative_err);
    }
}

void ComputingRelativeProtonCounts(string par_string)
{
    G_Sim_ProtonCounts_Relative[par_string] = new TGraphErrors();
    for (int i = 0; i < G_Exp_ProtonCounts->GetN(); i++)
    {
        double exp_y = G_Exp_ProtonCounts->GetY()[i];
        double sim_y = G_Sim_ProtonCounts[par_string]->GetY()[i];
        double exp_err = G_Exp_ProtonCounts->GetErrorY(i);
        double sim_err = G_Sim_ProtonCounts[par_string]->GetErrorY(i);
        double relative_y = (sim_y - exp_y) / exp_y;
        double relative_err = sqrt(pow(sim_err / exp_y, 2) + pow(exp_err * sim_y / (exp_y * exp_y), 2));
        G_Sim_ProtonCounts_Relative[par_string]->AddPoint(G_Exp_ProtonCounts->GetX()[i], relative_y);
        G_Sim_ProtonCounts_Relative[par_string]->SetPointError(G_Sim_ProtonCounts_Relative[par_string]->GetN() - 1, 0, relative_err);
    }
}

void PlottingExpvsSim(string par_string, string type = "")
{

    OUTPUT_FILE->cd();
    if (DIR_SIM == nullptr)
        DIR_SIM = OUTPUT_FILE->mkdir("Simulation");
    DIR_SIM->cd();

    if (type == "")
    {
        // contructing the canvas
        DIR_SIM->cd();
        TCanvas *C_Compare = new TCanvas(("C_Compare_" + par_string).c_str(), ("C_Compare_" + par_string).c_str(), 1200, 800);
        C_Compare->Divide(1, 3);

        // Proton counts
        C_Compare->cd(1);
        G_Sim_ProtonCounts_Relative[par_string]->SetTitle("Relative Difference in Proton Counts;Detector;(Sim - Exp) / Exp");
        G_Sim_ProtonCounts_Relative[par_string]->SetMarkerStyle(20);
        G_Sim_ProtonCounts_Relative[par_string]->SetMarkerColor(kRed);
        G_Sim_ProtonCounts_Relative[par_string]->SetMarkerSize(1);
        G_Sim_ProtonCounts_Relative[par_string]->Draw("AP");
        TLine *lineproton = new TLine(G_Exp_ProtonCounts->GetXaxis()->GetXmin(), 0, G_Exp_ProtonCounts->GetXaxis()->GetXmax(), 0);
        lineproton->SetLineColor(kBlack);
        lineproton->SetLineStyle(2);
        lineproton->Draw("SAME");
        TLatex *latex1 = new TLatex(0.2, 0.8, Form("Chi2: %.2f", Chi2_ProtonCounts[par_string]));
        latex1->SetNDC();
        latex1->SetTextColor(kRed);
        latex1->Draw("SAME");

        // Coinc/Single
        C_Compare->cd(2);
        G_Exp_CoincSingle->SetTitle("Coinc/Single Ratio;Code;Ratio");
        G_Exp_CoincSingle->SetMarkerStyle(20);
        G_Exp_CoincSingle->SetMarkerColor(kBlack);
        G_Exp_CoincSingle->SetMarkerSize(1);
        G_Exp_CoincSingle->Draw("AP");
        G_Sim_CoincSingle[par_string]->SetMarkerStyle(20);
        G_Sim_CoincSingle[par_string]->SetMarkerColor(kRed);
        G_Sim_CoincSingle[par_string]->SetMarkerSize(1);
        G_Sim_CoincSingle[par_string]->Draw("P SAME");
        TLatex *latex2 = new TLatex(0.2, 0.8, Form("Chi2: %.2f", Chi2_CoincSingle[par_string]));
        latex2->SetNDC();
        latex2->SetTextColor(kRed);
        latex2->Draw("SAME");
        // set range
        double max_y = 0;
        double min_y = 1.;
        for (int i = 0; i < G_Exp_CoincSingle->GetN(); i++)
        {
            if (G_Exp_CoincSingle->GetY()[i] + G_Exp_CoincSingle->GetErrorY(i) > max_y)
                max_y = G_Exp_CoincSingle->GetY()[i] + G_Exp_CoincSingle->GetErrorY(i);
            if (G_Exp_CoincSingle->GetY()[i] - G_Exp_CoincSingle->GetErrorY(i) < min_y)
                min_y = G_Exp_CoincSingle->GetY()[i] - G_Exp_CoincSingle->GetErrorY(i);
        }
        for (int i = 0; i < G_Sim_CoincSingle[par_string]->GetN(); i++)
        {
            if (G_Sim_CoincSingle[par_string]->GetY()[i] + G_Sim_CoincSingle[par_string]->GetErrorY(i) > max_y)
                max_y = G_Sim_CoincSingle[par_string]->GetY()[i] + G_Sim_CoincSingle[par_string]->GetErrorY(i);
            if (G_Sim_CoincSingle[par_string]->GetY()[i] - G_Sim_CoincSingle[par_string]->GetErrorY(i) < min_y)
                min_y = G_Sim_CoincSingle[par_string]->GetY()[i] - G_Sim_CoincSingle[par_string]->GetErrorY(i);
        }
        G_Exp_CoincSingle->GetYaxis()->SetRangeUser(min_y * 0.9, max_y * 1.1);

        // Eshift
        C_Compare->cd(3);
        G_Sim_Eshift_Relative[par_string]->SetTitle("Relative Eshift;Code;Eshift_{sim} - (Eshift_{sim,strip} - Eshift_{exp,strip}) (keV)");
        G_Sim_Eshift_Relative[par_string]->SetMarkerStyle(20);
        G_Sim_Eshift_Relative[par_string]->SetMarkerColor(kRed);
        G_Sim_Eshift_Relative[par_string]->SetLineColor(kRed);
        G_Sim_Eshift_Relative[par_string]->SetMarkerSize(1);
        G_Sim_Eshift_Relative[par_string]->Draw("AP");
        TLine *line = new TLine(G_Sim_Eshift_Relative[par_string]->GetXaxis()->GetXmin(), 0, G_Sim_Eshift_Relative[par_string]->GetXaxis()->GetXmax(), 0);
        line->SetLineColor(kBlack);
        line->SetLineStyle(2);
        line->Draw("SAME");
        TLatex *latex3 = new TLatex(0.2, 0.8, Form("Chi2: %.2f", Chi2_Eshift[par_string]));
        latex3->SetNDC();
        latex3->SetTextColor(kRed);
        latex3->Draw("SAME");

        C_Compare->Write();
    }
    else
    {
        // rewritting the canvas
        TDirectory *dir = OUTPUT_FILE->mkdir(type.c_str());
        dir->cd();
        TCanvas *C_Compare = (TCanvas *)OUTPUT_FILE->Get(("Simulation/C_Compare_" + par_string).c_str());
        if (C_Compare == nullptr)
        {
            Warning("Canvas not found for " + par_string + " and type " + type);
            return;
        }
        C_Compare->Write();
    }
}

void PlottingScans()
{
    OUTPUT_FILE->cd();

    // flat on x
    vector<double> x_values, y_values;
    for (int i = 0; i < G2_CoincSingle_YZ->GetN(); i++)
    {
        double x = G2_CoincSingle_YZ->GetX()[i];
        double y = G2_CoincSingle_YZ->GetY()[i];
        x_values.push_back(x);
        y_values.push_back(y);
    }
    sort(x_values.begin(), x_values.end());
    sort(y_values.begin(), y_values.end());
    auto last_x = unique(x_values.begin(), x_values.end());
    auto last_y = unique(y_values.begin(), y_values.end());
    x_values.erase(last_x, x_values.end());
    y_values.erase(last_y, y_values.end());
    double min_step_x = FindMinimalStep(x_values);
    double min_step_y = FindMinimalStep(y_values);

    double x_maxi = abs(x_values.back()) > abs(x_values.front()) ? abs(x_values.back()) : abs(x_values.front());
    double x_min = -x_maxi - min_step_x / 2;
    double x_max = x_maxi + min_step_x / 2;

    double y_maxi = abs(y_values.back()) > abs(y_values.front()) ? abs(y_values.back()) : abs(y_values.front());
    double y_min = -y_maxi - min_step_y / 2;
    double y_max = y_maxi + min_step_y / 2;

    int n_bins_x = (int)((x_max - x_min) / min_step_x) * 10;
    int n_bins_y = (int)((y_max - y_min) / min_step_y) * 10;

    if (TYPE == "Detectors")
    {

        TCanvas *C_Scan_YZ = new TCanvas("C_Scan_YZ", "C_Scan_YZ", 1200, 800);
        C_Scan_YZ->Divide(4, 1);
        // Proton Counts
        TH2D *H_ProtonCounts_YZ = new TH2D("H_ProtonCounts_YZ", "Proton Counts;Y (mm);Z (mm)", n_bins_x, x_min, x_max, n_bins_y, y_min, y_max);
        // Coinc/Single
        TH2D *H_CoincSingle_YZ = new TH2D("H_CoincSingle_YZ", "Coinc/Single Ratio;Y (mm);Z (mm)", n_bins_x, x_min, x_max, n_bins_y, y_min, y_max);
        // Eshift
        TH2D *H_Eshift_YZ = new TH2D("H_Eshift_YZ", "Eshift;Y (mm);Z (mm)", n_bins_x, x_min, x_max, n_bins_y, y_min, y_max);
        // TOTAL
        TH2D *H_Total_YZ = new TH2D("H_Total_YZ", "Total;Y (mm);Z (mm)", n_bins_x, x_min, x_max, n_bins_y, y_min, y_max);

        // TGraph to display real points
        TGraph *G_YZ = new TGraph();
        G_YZ->SetMarkerStyle(20);
        G_YZ->SetMarkerSize(1);
        for (int i = 0; i < G2_ProtonCounts_YZ->GetN(); i++)
        {
            double x = G2_ProtonCounts_YZ->GetX()[i];
            double y = G2_ProtonCounts_YZ->GetY()[i];
            G_YZ->AddPoint(x, y);
        }
        for (int i = 1; i <= H_ProtonCounts_YZ->GetNbinsX(); i++)
        {
            for (int j = 1; j <= H_ProtonCounts_YZ->GetNbinsY(); j++)
            {
                if (G2_ProtonCounts_YZ->GetN() < 3)
                    continue;
                double chi2;

                // Proton Counts
                chi2 = G2_ProtonCounts_YZ->Interpolate(H_ProtonCounts_YZ->GetXaxis()->GetBinCenter(i), H_ProtonCounts_YZ->GetYaxis()->GetBinCenter(j));
                H_ProtonCounts_YZ->SetBinContent(i, j, chi2);

                // Coinc/Single
                chi2 = G2_CoincSingle_YZ->Interpolate(H_CoincSingle_YZ->GetXaxis()->GetBinCenter(i), H_CoincSingle_YZ->GetYaxis()->GetBinCenter(j));
                H_CoincSingle_YZ->SetBinContent(i, j, chi2);

                // Eshift
                chi2 = G2_Eshift_YZ->Interpolate(H_Eshift_YZ->GetXaxis()->GetBinCenter(i), H_Eshift_YZ->GetYaxis()->GetBinCenter(j));
                H_Eshift_YZ->SetBinContent(i, j, chi2);

                // TOTAL
                chi2 = G2_Total_YZ->Interpolate(H_Total_YZ->GetXaxis()->GetBinCenter(i), H_Total_YZ->GetYaxis()->GetBinCenter(j));
                H_Total_YZ->SetBinContent(i, j, chi2);
            }
        }
        C_Scan_YZ->cd(1);
        H_ProtonCounts_YZ->SetMinimum(GetMin(H_ProtonCounts_YZ));
        H_ProtonCounts_YZ->Draw("COLZ");
        G_YZ->Draw("P SAME");

        C_Scan_YZ->cd(2);
        H_CoincSingle_YZ->SetMinimum(GetMin(H_CoincSingle_YZ));
        H_CoincSingle_YZ->Draw("COLZ");
        G_YZ->Draw("P SAME");

        C_Scan_YZ->cd(3);
        H_Eshift_YZ->SetMinimum(GetMin(H_Eshift_YZ));
        H_Eshift_YZ->Draw("COLZ");
        G_YZ->Draw("P SAME");

        C_Scan_YZ->cd(4);
        H_Total_YZ->SetMinimum(GetMin(H_Total_YZ));
        H_Total_YZ->Draw("COLZ");
        G_YZ->Draw("P SAME");

        C_Scan_YZ->Write();
    }


    // flat on z
    x_values = vector<double>();
    y_values = vector<double>();
    for (int i = 0; i < G2_CoincSingle_XY->GetN(); i++)
    {
        double x = G2_CoincSingle_XY->GetX()[i];
        double y = G2_CoincSingle_XY->GetY()[i];
        x_values.push_back(x);
        y_values.push_back(y);
    }
    sort(x_values.begin(), x_values.end());
    sort(y_values.begin(), y_values.end());
    last_x = unique(x_values.begin(), x_values.end());
    last_y = unique(y_values.begin(), y_values.end());
    x_values.erase(last_x, x_values.end());
    y_values.erase(last_y, y_values.end());
    min_step_x = FindMinimalStep(x_values);
    min_step_y = FindMinimalStep(y_values);

    x_maxi = abs(x_values.back()) > abs(x_values.front()) ? abs(x_values.back()) : abs(x_values.front());
    x_min = -x_maxi - min_step_x / 2;
    x_max = x_maxi + min_step_x / 2;

    y_maxi = abs(y_values.back()) > abs(y_values.front()) ? abs(y_values.back()) : abs(y_values.front());
    y_min = -y_maxi - min_step_y / 2;
    y_max = y_maxi + min_step_y / 2;

    n_bins_x = (int)((x_max - x_min) / min_step_x) * 10;
    n_bins_y = (int)((y_max - y_min) / min_step_y) * 10;

    TCanvas *C_Scan_XY = new TCanvas("C_Scan_XY", "C_Scan_XY", 1200, 800);
    C_Scan_XY->Divide(4, 1);
    // Proton Counts
    TH2D *H_ProtonCounts_XY = new TH2D("H_ProtonCounts_XY", "Proton Counts;X (mm);Y (mm)", n_bins_x, x_min, x_max, n_bins_y, y_min, y_max);
    // Coinc/Single
    TH2D *H_CoincSingle_XY = new TH2D("H_CoincSingle_XY", "Coinc/Single Ratio;X (mm);Y (mm)", n_bins_x, x_min, x_max, n_bins_y, y_min, y_max);
    // Eshift
    TH2D *H_Eshift_XY = new TH2D("H_Eshift_XY", "Eshift;X (mm);Y (mm)", n_bins_x, x_min, x_max, n_bins_y, y_min, y_max);
    // TOTAL
    TH2D *H_Total_XY = new TH2D("H_Total_XY", "Total;X (mm);Y (mm)", n_bins_x, x_min, x_max, n_bins_y, y_min, y_max);

    // TGraph to display real points
    TGraph *G_XY = new TGraph();
    G_XY->SetMarkerStyle(20);
    G_XY->SetMarkerSize(1);
    for (int i = 0; i < G2_ProtonCounts_XY->GetN(); i++)
    {
        double x = G2_ProtonCounts_XY->GetX()[i];
        double y = G2_ProtonCounts_XY->GetY()[i];
        G_XY->AddPoint(x, y);
    }

    for (int i = 1; i <= H_ProtonCounts_XY->GetNbinsX(); i++)
    {
        for (int j = 1; j <= H_ProtonCounts_XY->GetNbinsY(); j++)
        {
            double chi2;

            // Proton Counts
            chi2 = G2_ProtonCounts_XY->Interpolate(H_ProtonCounts_XY->GetXaxis()->GetBinCenter(i), H_ProtonCounts_XY->GetYaxis()->GetBinCenter(j));
            H_ProtonCounts_XY->SetBinContent(i, j, chi2);

            // Coinc/Single
            chi2 = G2_CoincSingle_XY->Interpolate(H_CoincSingle_XY->GetXaxis()->GetBinCenter(i), H_CoincSingle_XY->GetYaxis()->GetBinCenter(j));
            if (TYPE == "Catcher")
                H_CoincSingle_XY->SetBinContent(i, j, G2_MeanDiff_CoincSingle_XY->Interpolate(H_ProtonCounts_XY->GetXaxis()->GetBinCenter(i), H_ProtonCounts_XY->GetYaxis()->GetBinCenter(j)));
            else
                H_CoincSingle_XY->SetBinContent(i, j, chi2);

            // Eshift
            chi2 = G2_Eshift_XY->Interpolate(H_Eshift_XY->GetXaxis()->GetBinCenter(i), H_Eshift_XY->GetYaxis()->GetBinCenter(j));
            H_Eshift_XY->SetBinContent(i, j, chi2);

            // TOTAL
            chi2 = G2_Total_XY->Interpolate(H_Total_XY->GetXaxis()->GetBinCenter(i), H_Total_XY->GetYaxis()->GetBinCenter(j));
            H_Total_XY->SetBinContent(i, j, chi2);
        }
    }
    C_Scan_XY->cd(1);
    H_ProtonCounts_XY->SetMinimum(GetMin(H_ProtonCounts_XY));
    H_ProtonCounts_XY->Draw("COLZ");
    G_XY->Draw("P SAME");

    C_Scan_XY->cd(2);
    H_CoincSingle_XY->SetMinimum(GetMin(H_CoincSingle_XY));
    H_CoincSingle_XY->Draw("COLZ");
    G_XY->Draw("P SAME");

    C_Scan_XY->cd(3);
    H_Eshift_XY->SetMinimum(GetMin(H_Eshift_XY));
    H_Eshift_XY->Draw("COLZ");
    G_XY->Draw("P SAME");

    C_Scan_XY->cd(4);
    H_Total_XY->SetMinimum(GetMin(H_Total_XY));
    H_Total_XY->Draw("COLZ");
    G_XY->Draw("P SAME");

    C_Scan_XY->Write();

    if (TYPE == "Catcher")
    {

        //////////// FITTING ON CHI2 MAPS ////////////
        Info("Fitting on Chi2", 2);
        // 1D (x = 0)
        TCanvas *C_Scan_Y = new TCanvas("C_Scan_Y", "C_Scan_Y", 1200, 800);
        TGraph *G_Y = new TGraph();
        for (int i = 0; i < G2_CoincSingle_XY->GetN(); i++)
        {
            double x = G2_CoincSingle_XY->GetX()[i];
            if (x != 0)
                continue;
            double y = G2_CoincSingle_XY->GetY()[i];
            double chi2 = G2_CoincSingle_XY->GetZ()[i];
            G_Y->AddPoint(y, chi2);
        }
        G_Y->SetTitle("Coinc/Single Ratio at x=0;Y (mm);Ratio");
        G_Y->SetMarkerStyle(20);
        G_Y->SetMarkerColor(kBlack);
        G_Y->SetMarkerSize(2);

        TF1 *f = new TF1("asympol2", AsymetricPol2, 1., 3.0, 4);
        f->SetParameters(2.1, 0, 0, 0);
        TFitResultPtr r = G_Y->Fit(f, "RS");
        G_Y->Draw("AP");
        C_Scan_Y->Write();

        if (r->IsValid())
        {
            Info("1D Fit on Coinc/Single Ratio Results:");
            Warning("Z is fixed!");
            Info("χ2 / NDF = " + to_string(r->Chi2() / r->Ndf()));
            Info(Form("Y0 = %.2f +/- %.2f mm", f->GetParameter(0), f->GetParError(0)));
        }
         else
        {
            Warning("1D Fit on Coinc/Single Ratio failed!");
        }

        // 2D fit on a circle
        double x_min_fit = -1;
        double x_max_fit = 1;
        double y_min_fit = 1.5;
        double y_max_fit = 2.5;
        TGraph2D *G2_CoincSingle_XY_Fit = new TGraph2D();
        for (int i = 0; i < G2_CoincSingle_XY->GetN(); i++)
        {
            double x = G2_CoincSingle_XY->GetX()[i];
            double y = G2_CoincSingle_XY->GetY()[i];
            double z = G2_CoincSingle_XY->GetZ()[i];
            if (x < x_min_fit || x > x_max_fit || y < y_min_fit || y > y_max_fit)
                continue;
            G2_CoincSingle_XY_Fit->AddPoint(x, y, z);
        }


        TF2 * f2 = new TF2("asympol2_2D", AsymetricPol2_Circle2D, x_min_fit, x_max_fit, y_min_fit, y_max_fit, 6);
        f2->SetParLimits(0, 0., 100000);
        f2->SetParLimits(1, 0., 100000);
        f2->SetParLimits(2, 0., 100000);
        f2->SetParameter(2, 55);
        // X0
        f2->SetParLimits(3, -1, 1);
        f2->SetParameter(3, 0);
        // Y0
        f2->SetParLimits(4, -3, 2);
        f2->SetParameter(4, -0.5);
        // R
        f2->SetParLimits(5, 0.1, 4);
        f2->SetParameter(5, 2);
        TCanvas *C_Scan_2D = new TCanvas("C_Scan_2D", "C_Scan_2D", 1200, 800);
        r = G2_CoincSingle_XY_Fit->Fit(f2, "RS", "");
        G2_CoincSingle_XY_Fit->SetMarkerStyle(20);
        G2_CoincSingle_XY_Fit->SetMarkerColor(kBlack);
        G2_CoincSingle_XY_Fit->SetMarkerSize(2);
        G2_CoincSingle_XY_Fit->Draw("P");  
        f2->Draw("SURF SAME");
        C_Scan_2D->Write();

        if (r->IsValid())
        {
            Info("2D Fit Results:");
            Warning("Z is fixed!");
            Info("χ2 / NDF = " + to_string(r->Chi2() / r->Ndf()));
            Info(Form("Y0 = %.2f +/- %.2f mm", sqrt(pow(f2->GetParameter(5), 2) + pow(f2->GetParameter(3), 2)) + f2->GetParameter(4), sqrt(pow(f2->GetParError(5), 2) + pow(f2->GetParError(3), 2) + pow(f2->GetParError(4), 2))));
        }

        Info("Fitting on Diff Ratio", 2);
        // 1D (x = 0)
        TCanvas *C_Scan_Y_Diff = new TCanvas("C_Scan_Y_Diff", "C_Scan_Y_Diff", 1200, 800);
        TGraphErrors *G_Y_Diff = new TGraphErrors();
        for (int i = 0; i < G2_MeanDiff_CoincSingle_XY->GetN(); i++)
        {
            double x = G2_MeanDiff_CoincSingle_XY->GetX()[i];
            if (x != 0)
                continue;
            double y = G2_MeanDiff_CoincSingle_XY->GetY()[i];
            double value = G2_MeanDiff_CoincSingle_XY->GetZ()[i];
            double err = G2_MeanDiff_CoincSingle_XY->GetErrorZ(i);
            G_Y_Diff->AddPoint(y, value);
            G_Y_Diff->SetPointError(G_Y_Diff->GetN() - 1, 0, err);
        }

        G_Y_Diff->SetTitle("Mean Diff Coinc/Single Ratio at x=0;Y (mm);Mean Diff Ratio");
        G_Y_Diff->SetMarkerStyle(20);
        G_Y_Diff->SetMarkerColor(kBlack);
        G_Y_Diff->SetMarkerSize(2);
        TF1 *f_diff = new TF1("asympol2_diff", AsymetricPol2, 1.25, 2.5, 4);
        f_diff->SetParameters(2.1, 0, 0, 0);
        r = G_Y_Diff->Fit(f_diff, "RS");
        G_Y_Diff->Draw("AP");
        C_Scan_Y_Diff->Write();

        if (r->IsValid())
        {
            Info("1D Fit on Diff Ratio Results:");
            Warning("Z is fixed!");
            Info("χ2 / NDF = " + to_string(r->Chi2() / r->Ndf()));
            Info(Form("Y0 = %.2f +/- %.2f mm", f_diff->GetParameter(0), f_diff->GetParError(0)));
        }
         else
        {
            Warning("1D Fit on Diff Ratio did not converge!");
        }
        // 2D fit on a circle
        TCanvas *C_Scan_2D_Diff = new TCanvas("C_Scan_2D_Diff", "C_Scan_2D_Diff", 1200, 800);
        TGraph2DErrors *G2_CoincSingle_XY_Diff_Fit = new TGraph2DErrors(); 

        x_min_fit = -3;
        x_max_fit = 3;
        y_min_fit = -3;
        y_max_fit = 3;
        for (int i = 0; i < G2_MeanDiff_CoincSingle_XY->GetN(); i++)
        {
            double x = G2_MeanDiff_CoincSingle_XY->GetX()[i];
            double y = G2_MeanDiff_CoincSingle_XY->GetY()[i];
            double z = G2_MeanDiff_CoincSingle_XY->GetZ()[i];
            double err = G2_MeanDiff_CoincSingle_XY->GetErrorZ(i);
            if (x < x_min_fit || x > x_max_fit || y < y_min_fit || y > y_max_fit)
                continue;
            G2_CoincSingle_XY_Diff_Fit->AddPoint(x, y, z);
            G2_CoincSingle_XY_Diff_Fit->SetPointError(G2_CoincSingle_XY_Diff_Fit->GetN() - 1, 0, 0, err);
        }


        TF2 * f2_diff = new TF2("asympol2_2D_diff", AsymetricPol2_Circle2D, x_min_fit, x_max_fit, y_min_fit, y_max_fit, 6);
        f2_diff->SetParLimits(0, 0., 1);
        f2_diff->SetParLimits(1, 0., 1);
        f2_diff->SetParLimits(2, 0., 1);
        // f2_diff->SetParameter(2, 55);
        // X0
        f2_diff->SetParLimits(3, -0.25, 0.25);
        f2_diff->SetParameter(3, 0);
        // Y0
        f2_diff->SetParLimits(4, -2, 1);
        f2_diff->SetParameter(4, -0.5);
        // R
        f2_diff->SetParLimits(5, 1.5, 3);
        f2_diff->SetParameter(5, 2);
        r = G2_CoincSingle_XY_Diff_Fit->Fit(f2_diff, "RS", "");
        G2_CoincSingle_XY_Diff_Fit->SetMarkerStyle(20);
        G2_CoincSingle_XY_Diff_Fit->SetMarkerColor(kBlack);
        G2_CoincSingle_XY_Diff_Fit->SetMarkerSize(2);
        G2_CoincSingle_XY_Diff_Fit->Draw("P");
        G2_CoincSingle_XY_Diff_Fit->GetXaxis()->SetTitle("X (mm)");
        G2_CoincSingle_XY_Diff_Fit->GetYaxis()->SetTitle("Y (mm)");
        G2_CoincSingle_XY_Diff_Fit->GetZaxis()->SetTitle("Mean Diff Ratio");
        // draw axis

        f2_diff->Draw("SURF SAME");
        C_Scan_2D_Diff->Write();
        if (r->IsValid())
        {
            Info("2D Fit on Diff Ratio Results:");
            Warning("Z is fixed!");
            Info("χ2 / NDF = " + to_string(r->Chi2() / r->Ndf()));
            double R = f2_diff->GetParameter(5);
            double X0 = f2_diff->GetParameter(3);
            double Y0 = f2_diff->GetParameter(4);
            double err_R = f2_diff->GetParError(5);
            double err_X0 = f2_diff->GetParError(3);
            double err_Y0 = f2_diff->GetParError(4);
            double value = sqrt(pow(R, 2) - pow(X0, 2)) + Y0;
            double err = sqrt(  1./(pow(R, 2) + pow(X0, 2)) * (pow(R * err_R, 2) + pow(X0 * err_X0, 2)) + pow(err_Y0, 2) );
            Info(Form("Y0 = %.2f +/- %.2f mm", value, err));
        }
        else
        {
            Warning("2D Fit on Diff Ratio did not converge!");
        }        
        

    }

    if (TYPE == "Detectors")
    {

        /// PROTON COUNTS ///
        Info(" --- Fitting Proton Counts --- ");
        Info("Fitting XY", 2);
        // 2D fit in a global minimum 2D asym XY 
        double x_min_fit = -2;
        double x_max_fit = 2;
        double y_min_fit = 0.;
        double y_max_fit = 3.;
        TGraph2DErrors *G2_ProtonCounts_XY_Fit = new TGraph2DErrors();
        for (int i = 0; i < G2_ProtonCounts_XY->GetN(); i++)
        {
            double x = G2_ProtonCounts_XY->GetX()[i];
            double y = G2_ProtonCounts_XY->GetY()[i];
            double z = G2_ProtonCounts_XY->GetZ()[i];
            double err = G2_ProtonCounts_XY->GetErrorZ(i);
            if (x < x_min_fit || x > x_max_fit || y < y_min_fit || y > y_max_fit)
                continue;
            G2_ProtonCounts_XY_Fit->AddPoint(x, y, z);
            G2_ProtonCounts_XY_Fit->SetPointError(G2_ProtonCounts_XY_Fit->GetN() - 1, 0, 0, err);
        }

        
        TCanvas *C_Scan_XY = new TCanvas("C_Scan_XY_Fit", "C_Scan_XY_Fit", 1200, 800);
        TF2 * f2 = new TF2("asympol2_2D", AsymetricPol2_2D, x_min_fit, x_max_fit, y_min_fit, y_max_fit, 7);
        f2->SetParLimits(0, -1, 1);
        f2->SetParameter(0, 0.);
        f2->SetParLimits(1, 0.5, 3);
        f2->SetParameter(1, 2.);
        f2->SetParLimits(2, 0.0, 10);
        f2->SetParLimits(3, 0.0, 10);
        f2->SetParLimits(4, 0.0, 10);
        f2->SetParLimits(5, 0.0, 10);
        f2->SetParLimits(6, 0.0, 10);
        // set call limit
        TFitResultPtr r = G2_ProtonCounts_XY_Fit->Fit(f2, "RS", "");
        G2_ProtonCounts_XY_Fit->SetMarkerStyle(20);
        G2_ProtonCounts_XY_Fit->SetMarkerColor(kBlack);
        G2_ProtonCounts_XY_Fit->SetMarkerSize(2);
        G2_ProtonCounts_XY_Fit->Draw("P");
        f2->Draw("SURF SAME");
        C_Scan_XY->Write();

        if (r->IsValid())
        {
            Info("2D Fit Results:");
            Warning("Z is fixed!");
            Info("χ2 / NDF = " + to_string(r->Chi2() / r->Ndf()));
            Info(Form("X0 = %.2f +/- %.2f mm", f2->GetParameter(0), f2->GetParError(0)));
            Info(Form("Y0 = %.2f +/- %.2f mm", f2->GetParameter(1), f2->GetParError(1)));   
        } 
        else
        {
            Warning("2D Fit on XY did not converge!");
        }

        y_min_fit = 1.3;
        y_max_fit = 2.5;

        Info("Fitting YZ", 2);  
        // 2D fit in a global minimum 2D asym YZ
        double z_min_fit = 0;
        double z_max_fit = 0.8 ;
        TGraph2DErrors *G2_ProtonCounts_YZ_Fit = new TGraph2DErrors();
        for (int i = 0; i < G2_ProtonCounts_YZ->GetN(); i++)
        {
            double y = G2_ProtonCounts_YZ->GetX()[i];
            double z = G2_ProtonCounts_YZ->GetY()[i];
            double chi2 = G2_ProtonCounts_YZ->GetZ()[i];
            double err = G2_ProtonCounts_YZ->GetErrorZ(i);
            if (y < y_min_fit || y > y_max_fit || z < z_min_fit || z > z_max_fit)
                continue;
            G2_ProtonCounts_YZ_Fit->AddPoint(y, z, chi2);
            G2_ProtonCounts_YZ_Fit->SetPointError(G2_ProtonCounts_YZ_Fit->GetN() - 1, 0, 0, err);
        }
        TF2 * f2_yz = new TF2("asympol2_2D_YZ", AsymetricPol2_2D, y_min_fit, y_max_fit, z_min_fit, z_max_fit, 7);
        f2_yz->SetParLimits(0, 1.6, 2.5);
        f2_yz->SetParameter(0, 2.);
        f2_yz->SetParLimits(1, 0.2, 0.7);
        f2_yz->SetParameter(1, 0.);
        f2_yz->SetParLimits(2, 0.0, 10000);
        f2_yz->SetParameter(2, 600);
        f2_yz->SetParLimits(3, 0.0, 10000);
        f2_yz->SetParameter(3, 600);
        f2_yz->SetParLimits(4, 0.0, 10000);
        f2_yz->SetParameter(4, 600);
        f2_yz->SetParLimits(5, 0.0, 10000);
        f2_yz->SetParameter(5, 600);
        f2_yz->SetParLimits(6, 0., 10000);
        TCanvas *C_Scan_YZ = new TCanvas("C_Scan_YZ_Fit", "C_Scan_YZ_Fit", 1200, 800);
        r = G2_ProtonCounts_YZ_Fit->Fit(f2_yz, "RS", "");
        G2_ProtonCounts_YZ_Fit->SetMarkerStyle(20);
        G2_ProtonCounts_YZ_Fit->SetMarkerColor(kBlack);
        G2_ProtonCounts_YZ_Fit->SetMarkerSize(2);
        G2_ProtonCounts_YZ_Fit->Draw("P");
        f2_yz->Draw("SURF SAME");   
        C_Scan_YZ->Write();

        if (r->IsValid())
        {
            Info("2D Fit Results:");
            Warning("X is fixed!");
            Info("χ2 / NDF = " + to_string(r->Chi2() / r->Ndf()));
            Info(Form("Y0 = %.2f +/- %.2f mm", f2_yz->GetParameter(0), f2_yz->GetParError(0)));
            Info(Form("Z0 = %.2f +/- %.2f mm", f2_yz->GetParameter(1), f2_yz->GetParError(1)));   
        } 
        else
        {
            Warning("2D Fit on YZ did not converge!");
        }
                   

        Info("Fitting XYZ", 2);    
        // 3D Fit on a global minimum 3D asym
        vector<double> z_values;
        for (auto &entry : Chi2_ProtonCounts)
        {
            string par_string = entry.first;
            double chi2 = entry.second;
            double x, y, z;
            std::regex rgx("x=([-+]?[0-9]*\\.?[0-9]+) y=([-+]?[0-9]*\\.?[0-9]+) z=([-+]?[0-9]*\\.?[0-9]+)");
            std::smatch match;
            if (std::regex_search(par_string, match, rgx) && match.size() == 4)
            {
                x = std::stod(match[1].str());
                y = std::stod(match[2].str());
                z = std::stod(match[3].str());
                z_values.push_back(z);
            }
        }

        double z_maxi = abs(z_max_fit) > abs(z_min_fit) ? abs(z_max_fit) : abs(z_min_fit);
        double z_min = -z_maxi - min_step_y / 2;
        double z_max = z_maxi + min_step_y / 2;

        double min_step_z = FindMinimalStep(z_values);
        int n_bins_z = (int)((z_max - z_min) / min_step_z);
        
        TH3D *H_ProtonCounts_XYZ = new TH3D("H_ProtonCounts_XYZ", "Proton Counts;X (mm);Y (mm);Z", n_bins_x, x_min, x_max, n_bins_y, y_min, y_max, n_bins_z, z_min, z_max);
        for (auto &entry : Chi2_ProtonCounts)
        {
            string par_string = entry.first;
            double chi2 = entry.second;
            double x, y, z;
            std::regex rgx("x=([-+]?[0-9]*\\.?[0-9]+) y=([-+]?[0-9]*\\.?[0-9]+) z=([-+]?[0-9]*\\.?[0-9]+)");
            std::smatch match;
            if (std::regex_search(par_string, match, rgx) && match.size() == 4)
            {
                x = std::stod(match[1].str());
                y = std::stod(match[2].str());
                z = std::stod(match[3].str());
                // H_ProtonCounts_XYZ->Fill(x, y, z, chi2);
                pair<double, double > value = MeanDistance(G_Sim_ProtonCounts_Relative[par_string]);
                H_ProtonCounts_XYZ->Fill(x, y, z, value.first);
                H_ProtonCounts_XYZ->SetBinError(H_ProtonCounts_XYZ->GetXaxis()->FindBin(x), H_ProtonCounts_XYZ->GetYaxis()->FindBin(y), H_ProtonCounts_XYZ->GetZaxis()->FindBin(z), value.second);
            }
        }

        
        TF3 *f3 = new TF3("asympol2_3D", AsymetricPol2_3D, x_min_fit, x_max_fit, y_min_fit, y_max_fit, z_min_fit, z_max_fit, 10);
        f3->SetParLimits(0, -1, 1);
        f3->SetParameter(0, 0.);
        f3->SetParLimits(1, 1.5, 3);
        f3->SetParameter(1, 2.);
        f3->SetParLimits(2, -1, 1.);
        f3->SetParameter(2, 0.3);
        f3->SetParLimits(3, 0.0, 10);
        f3->SetParLimits(4, 0.0, 10);
        f3->SetParLimits(5, 0.0, 10);
        f3->SetParLimits(6, 0.0, 10);
        f3->SetParLimits(7, 0.0, 10);
        f3->SetParLimits(8, 0.0, 10);
        f3->SetParLimits(9, 0., 10000);
        ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 10e7);

        r = H_ProtonCounts_XYZ->Fit(f3, "RS", "");

        if (r->IsValid())
        {
            Info("3D Fit Results:");
            Info("χ2 / NDF = " + to_string(r->Chi2() / r->Ndf()));
            Info(Form("X0 = %.2f +/- %.2f mm", f3->GetParameter(0), f3->GetParError(0)));
            Info(Form("Y0 = %.2f +/- %.2f mm", f3->GetParameter(1), f3->GetParError(1)));
            Info(Form("Z0 = %.2f +/- %.2f mm", f3->GetParameter(2), f3->GetParError(2)));
        }
        else
        {
            Warning("3D Fit on XYZ did not converge!");
        }

        x_min_fit = -1.5;
        x_max_fit = 2;
        y_min_fit = 0.5;
        y_max_fit = 2.5;
        z_min_fit = 0.5;
        z_max_fit = 2;


        /// ESHIFT ///
        Info(" ----- Fitting ESHIFT ----- ");
            // 2D fit in a global minimum 2D asym XY
        TGraph2DErrors *G2_Eshift_XY_Fit = new TGraph2DErrors();
        for (int i = 0; i < G2_Eshift_XY->GetN(); i++)
        {
            double x = G2_Eshift_XY->GetX()[i];
            double y = G2_Eshift_XY->GetY()[i];
            double z = G2_Eshift_XY->GetZ()[i];
            double z_err = G2_Eshift_XY->GetErrorZ(i);
            if (x < x_min_fit || x > x_max_fit || y < y_min_fit || y > y_max_fit)
                continue;
            G2_Eshift_XY_Fit->AddPoint(x, y, z);
            G2_Eshift_XY_Fit->SetPointError(G2_Eshift_XY_Fit->GetN() - 1, 0, 0, z_err);
        }

        f2 = new TF2("asympol2_2D_Eshift", AsymetricPol2_2D, x_min_fit, x_max_fit, y_min_fit, y_max_fit, 7);
        f2->SetParLimits(0, -1, 1);
        f2->SetParameter(0, 0.);
        f2->SetParLimits(1, 1.5, 2.5);
        f2->SetParameter(1, 2.);
        f2->SetParLimits(2, 0.0, 1000);
        f2->SetParLimits(3, 0.0, 1000);
        f2->SetParLimits(4, 0.0, 1000);
        f2->SetParLimits(5, 0.0, 1000);
        f2->SetParLimits(6, 0., 10000);
        TCanvas *C_Scan_XY_Eshift = new TCanvas("C_Scan_XY_Eshift_Fit", "C_Scan_XY_Eshift_Fit", 1200, 800);
        r = G2_Eshift_XY_Fit->Fit(f2, "RS", "");
        G2_Eshift_XY_Fit->SetMarkerStyle(20);
        G2_Eshift_XY_Fit->SetMarkerColor(kBlack);
        G2_Eshift_XY_Fit->SetMarkerSize(2);
        G2_Eshift_XY_Fit->Draw("P");
        f2->Draw("SURF SAME");
        C_Scan_XY_Eshift->Write();

        if (r->IsValid())
        {
            Info("2D Fit Results:");
            Warning("Z is fixed!");
            Info("χ2 / NDF = " + to_string(r->Chi2() / r->Ndf()));
            Info(Form("X0 = %.2f +/- %.2f mm", f2->GetParameter(0), f2->GetParError(0)));
            Info(Form("Y0 = %.2f +/- %.2f mm", f2->GetParameter(1), f2->GetParError(1)));   
        } 
        else
        {
            Warning("2D Fit on XY did not converge!");
        }

        // 2D fit in a global minimum 2D asym YZ
        TGraph2DErrors *G2_Eshift_YZ_Fit = new TGraph2DErrors();
        for (int i = 0; i < G2_Eshift_YZ->GetN(); i++)
        {
            double y = G2_Eshift_YZ->GetX()[i];
            double z = G2_Eshift_YZ->GetY()[i];
            double chi2 = G2_Eshift_YZ->GetZ()[i];
            double err = G2_Eshift_YZ->GetErrorZ(i);
            if (y < y_min_fit || y > y_max_fit || z < z_min_fit || z > z_max_fit)
                continue;
            G2_Eshift_YZ_Fit->AddPoint(y, z, chi2);
            G2_Eshift_YZ_Fit->SetPointError(G2_Eshift_YZ_Fit->GetN() - 1, 0, 0, err);
        }
        TF2 * f2_yz_e = new TF2("asympol2_2D_YZ_Eshift", AsymetricPol2_2D, y_min_fit, y_max_fit, z_min_fit, z_max_fit, 7);
        f2_yz_e->SetParLimits(0, 1.5, 2.5);
        f2_yz_e->SetParameter(0, 2.);
        f2_yz_e->SetParLimits(1, 0.2, 2.0);
        f2_yz_e->SetParameter(1, 0.);
        f2_yz_e->SetParLimits(2, 0.0, 10000);
        f2_yz_e->SetParameter(2, 600);
        f2_yz_e->SetParLimits(3, 0.0, 10000);
        f2_yz_e->SetParameter(3, 600);
        f2_yz_e->SetParLimits(4, 0.0, 10000);
        f2_yz_e->SetParameter(4, 600);
        f2_yz_e->SetParLimits(5, 0.0, 10000);
        f2_yz_e->SetParameter(5, 600);
        f2_yz_e->SetParLimits(6, 0., 10000);
        TCanvas *C_Scan_YZ_Eshift = new TCanvas("C_Scan_YZ_Eshift_Fit", "C_Scan_YZ_Eshift_Fit", 1200, 800);
        r = G2_Eshift_YZ_Fit->Fit(f2_yz_e, "RS", "");
        G2_Eshift_YZ_Fit->SetMarkerStyle(20);
        G2_Eshift_YZ_Fit->SetMarkerColor(kBlack);
        G2_Eshift_YZ_Fit->SetMarkerSize(2);
        G2_Eshift_YZ_Fit->Draw("P");
        f2_yz_e->Draw("SURF SAME");
        C_Scan_YZ_Eshift->Write();
        if (r->IsValid())
        {
            Info("2D Fit Results:");
            Warning("X is fixed!");
            Info("χ2 / NDF = " + to_string(r->Chi2() / r->Ndf()));
            Info(Form("Y0 = %.2f +/- %.2f mm", f2_yz_e->GetParameter(0), f2_yz_e->GetParError(0)));
            Info(Form("Z0 = %.2f +/- %.2f mm", f2_yz_e->GetParameter(1), f2_yz_e->GetParError(1)));   
        } 
        else
        {
            Warning("2D Fit on YZ did not converge!");
        }

        // 3D Fit on a global minimum 3D asym
        TH3D *H_Eshift_XYZ = new TH3D("H_Eshift_XYZ", "Eshift;X (mm);Y (mm);Z", n_bins_x, x_min, x_max, n_bins_y, y_min, y_max, n_bins_z, z_min, z_max);
        for (auto &entry : Chi2_Eshift)
        {
            string par_string = entry.first;
            double chi2 = entry.second;
            double x, y, z;
            std::regex rgx("x=([-+]?[0-9]*\\.?[0-9]+) y=([-+]?[0-9]*\\.?[0-9]+) z=([-+]?[0-9]*\\.?[0-9]+)");
            std::smatch match;
            if (std::regex_search(par_string, match, rgx) && match.size() == 4)
            {
                x = std::stod(match[1].str());
                y = std::stod(match[2].str());
                z = std::stod(match[3].str());
                H_Eshift_XYZ->Fill(x, y, z, chi2);
            }
        }

        f3 = new TF3("asympol2_3D_Eshift", AsymetricPol2_3D, x_min_fit, x_max_fit, y_min_fit, y_max_fit, z_min_fit, z_max_fit, 10);
        f3->SetParLimits(0, -1, 1);
        f3->SetParameter(0, 0.);
        f3->SetParLimits(1, 0.5, 3);
        f3->SetParameter(1, 2.);
        f3->SetParLimits(2, -1, 2.);
        f3->SetParameter(2, 0.);
        f3->SetParLimits(3, 0.0, 1000000);
        f3->SetParLimits(4, 0.0, 1000000);
        f3->SetParLimits(5, 0.0, 1000000);
        f3->SetParLimits(6, 0.0, 1000000);
        f3->SetParLimits(7, 0.0, 1000000);
        f3->SetParLimits(8, 0.0, 1000000);
        f3->SetParLimits(9, 0., 10000);
        r = H_Eshift_XYZ->Fit(f3, "RSL", "");

        if (r->IsValid())
        {
            Info("3D Fit Results:");
            Info("χ2 / NDF = " + to_string(r->Chi2() / r->Ndf()));
            Info(Form("X0 = %.2f +/- %.2f mm", f3->GetParameter(0), f3->GetParError(0)));
            Info(Form("Y0 = %.2f +/- %.2f mm", f3->GetParameter(1), f3->GetParError(1)));
            Info(Form("Z0 = %.2f +/- %.2f mm", f3->GetParameter(2), f3->GetParError(2)));
        }
        else
        {
            Warning("3D Fit on XYZ did not converge!");
        }


    }

}

void Plotting3D()
{

    vector<pair<string, double>> Chi2_ProtonCounts_vec(Chi2_ProtonCounts.begin(), Chi2_ProtonCounts.end());
    sort(Chi2_ProtonCounts_vec.begin(), Chi2_ProtonCounts_vec.end(), [](const pair<string, double> &a, const pair<string, double> &b)
         { return a.second < b.second; });

    // find ranges or x, y, z
    double x_min = 1e9, x_max = -1e9;
    double y_min = 1e9, y_max = -1e9;
    double z_min = 1e9, z_max = -1e9;
    for (auto &entry : Chi2_ProtonCounts_vec)
    {
        string par_string = entry.first;
        double chi2 = entry.second;
        double x, y, z;
        std::regex rgx("x=([-+]?[0-9]*\\.?[0-9]+) y=([-+]?[0-9]*\\.?[0-9]+) z=([-+]?[0-9]*\\.?[0-9]+)");
        std::smatch match;
        if (std::regex_search(par_string, match, rgx) && match.size() == 4)
        {
            x = std::stod(match[1].str());
            y = std::stod(match[2].str());
            z = std::stod(match[3].str());
        }
        if (x < x_min)
            x_min = x;
        if (x > x_max)
            x_max = x;
        if (y < y_min)
            y_min = y;
        if (y > y_max)
            y_max = y;
        if (z < z_min)
            z_min = z;
        if (z > z_max)
            z_max = z;
    }

    TCanvas *C_3D = new TCanvas("C_3D", "C_3D", 1200, 800);
    TGraph2D *dump = new TGraph2D();
    dump->SetPoint(0, x_min, y_min, z_min);
    dump->SetPoint(1, x_max, y_max, z_max);
    dump->GetXaxis()->SetTitle("X (mm)");
    dump->GetYaxis()->SetTitle("Y (mm)");
    dump->GetZaxis()->SetTitle("Z (mm)");
    dump->Draw("P");

    map<string, TGraph2D *> Chi2_ProtonCounts_3D;
    for (auto &entry : Chi2_ProtonCounts)
    {
        string par_string = entry.first;

        double chi2 = entry.second;

        double x, y, z;

        std::regex rgx("x=([-+]?[0-9]*\\.?[0-9]+) y=([-+]?[0-9]*\\.?[0-9]+) z=([-+]?[0-9]*\\.?[0-9]+)");
        std::smatch match;
        if (std::regex_search(par_string, match, rgx) && match.size() == 4)
        {
            x = std::stod(match[1].str());
            y = std::stod(match[2].str());
            z = std::stod(match[3].str());
        }

        Chi2_ProtonCounts_3D[entry.first] = new TGraph2D();
        Chi2_ProtonCounts_3D[entry.first]->SetPoint(0, x, y, z);
        Chi2_ProtonCounts_3D[entry.first]->SetMarkerStyle(20);
        Chi2_ProtonCounts_3D[entry.first]->SetMarkerColor(GetColor(chi2, Chi2_ProtonCounts_vec.front().second, Chi2_ProtonCounts_vec.back().second));
        Chi2_ProtonCounts_3D[entry.first]->Draw("P SAME");
    }

    C_3D->Write();
}

void AddorReplace(TGraph2D *G, double X, double Y, double CHI2)
{
    for (int i = 0; i < G->GetN(); i++)
    {
        double x, y, chi2;
        G->GetPoint(i, x, y, chi2);
        if (x == X && y == Y && chi2 > CHI2)
        {
            G->SetPoint(i, X, Y, CHI2);
            return;
        }
    }
    G->AddPoint(X, Y, CHI2);
}

void PlottingStackedScan()
{
    if (TYPE != "Detectors")
        return;
    // same purpose as PlottingScan(); but with all the data (projection xy and yz but with the best value not plotted)
    vector<double> x_values, y_values, z_values;
    TGraph2D *G2_ESHIFT_XY_BEST = new TGraph2D();
    TGraph2D *G2_ESHIFT_YZ_BEST = new TGraph2D();
    for (auto &entry : Chi2_Eshift)
    {
        string par_string = entry.first;

        double chi2 = entry.second;

        double x, y, z;

        std::regex rgx("x=([-+]?[0-9]*\\.?[0-9]+) y=([-+]?[0-9]*\\.?[0-9]+) z=([-+]?[0-9]*\\.?[0-9]+)");
        std::smatch match;
        if (std::regex_search(par_string, match, rgx) && match.size() == 4)
        {
            x = std::stod(match[1].str());
            y = std::stod(match[2].str());
            z = std::stod(match[3].str());
        }


        // Replace if better
        AddorReplace(G2_ESHIFT_XY_BEST, x, y, chi2);
        AddorReplace(G2_ESHIFT_YZ_BEST, y, z, chi2);
        x_values.push_back(x);
        y_values.push_back(y);
        z_values.push_back(z);
    }
    sort(x_values.begin(), x_values.end());
    sort(y_values.begin(), y_values.end());
    sort(z_values.begin(), z_values.end());
    auto last_x = unique(x_values.begin(), x_values.end());
    auto last_y = unique(y_values.begin(), y_values.end());
    auto last_z = unique(z_values.begin(), z_values.end());
    x_values.erase(last_x, x_values.end());
    y_values.erase(last_y, y_values.end());
    z_values.erase(last_z, z_values.end());
    double min_step_x = FindMinimalStep(x_values);
    double min_step_y = FindMinimalStep(y_values);
    double min_step_z = FindMinimalStep(z_values);

    double x_maxi = abs(x_values.back()) > abs(x_values.front()) ? abs(x_values.back()) : abs(x_values.front());
    double x_min = -x_maxi - min_step_x / 2;
    double x_max = x_maxi + min_step_x / 2;

    double y_maxi = abs(y_values.back()) > abs(y_values.front()) ? abs(y_values.back()) : abs(y_values.front());
    double y_min = -y_maxi - min_step_y / 2;
    double y_max = y_maxi + min_step_y / 2;

    double z_maxi = abs(z_values.back()) > abs(z_values.front()) ? abs(z_values.back()) : abs(z_values.front());
    double z_min = -z_maxi - min_step_z / 2;
    double z_max = z_maxi + min_step_z / 2;

    int n_bins_x = (int)((x_max - x_min) / min_step_x) * 10;
    int n_bins_y = (int)((y_max - y_min) / min_step_y) * 10;
    int n_bins_z = (int)((z_max - z_min) / min_step_z) * 10;

    TCanvas *C_Scan = new TCanvas("C_Scan_BEST", "C_Scan_BEST", 1200, 800);
    C_Scan->Divide(2, 1);
    // Eshift
    TH2D *H_Eshift_XY = new TH2D("H_Eshift_XY_BEST", "Eshift;X (mm);Y (mm)", n_bins_x, x_min, x_max, n_bins_y, y_min, y_max);
    TH2D *H_Eshift_YZ = new TH2D("H_Eshift_YZ_BEST", "Eshift;Y (mm);Z (mm)", n_bins_y, y_min, y_max, n_bins_z, z_min, z_max);

    for (int i = 1; i <= H_Eshift_XY->GetNbinsX(); i++)
    {
        for (int j = 1; j <= H_Eshift_XY->GetNbinsY(); j++)
        {
            double chi2;

            // Proton Counts
            chi2 = G2_ESHIFT_XY_BEST->Interpolate(H_Eshift_XY->GetXaxis()->GetBinCenter(i), H_Eshift_XY->GetYaxis()->GetBinCenter(j));
            H_Eshift_XY->SetBinContent(i, j, chi2);
        }
    }

    for (int i = 1; i <= H_Eshift_YZ->GetNbinsX(); i++)
    {
        for (int j = 1; j <= H_Eshift_YZ->GetNbinsY(); j++)
        {
            double chi2;

            // Proton Counts
            chi2 = G2_ESHIFT_YZ_BEST->Interpolate(H_Eshift_YZ->GetXaxis()->GetBinCenter(i), H_Eshift_YZ->GetYaxis()->GetBinCenter(j));
            H_Eshift_YZ->SetBinContent(i, j, chi2);
        }
    }

    C_Scan->cd(1);
    H_Eshift_XY->SetMinimum(GetMin(H_Eshift_XY));
    H_Eshift_XY->Draw("COLZ");

    C_Scan->cd(2);
    H_Eshift_YZ->SetMinimum(GetMin(H_Eshift_YZ));
    H_Eshift_YZ->Draw("COLZ");

    C_Scan->Write();
}
