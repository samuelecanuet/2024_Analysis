#ifndef SIPM_DARK_HH
#define SIPM_DARK_HH

#include "../../../lib/SignalDict/Signal.h"
#include "Detectors.hh"

vector<string> Runs = {
    "092", "093", "094", "095", "096", "097", "098", "099", "100",
    "101", "102", "103", "104", "105", "106", "107", "108", "109"
    }; 
// vector<string> Runs = {
//     "092", "093", "097", "098", "099", "108", "109"
//     }; 

map<string, double[10]> Thresholds;
map<string, TH1D *[10]> H_High;
map<string, vector<bool>> Taking;

double osc_function(double *x, double *par)
{
    double th_level = par[0];
    double th_sigma = par[1];
    
    double offset = par[2];
    double delta_x = par[3];

    double A[10] = {par[6], par[7], par[8], par[9], par[10], par[11], par[12], par[13], par[14], par[15]};
    double result = 0;

    for (int i = 0; i < 10; i++)
    {
        double sigma = par[4] + i * par[5];
        double AA = par[6] * exp(i* par[7]);
        result += AA/(sqrt(2*M_PI)*sigma) * exp(-pow(delta_x*i + offset - x[0], 2) / (2 * pow(sigma,  2)));
        // result += A[i] * exp(-pow(delta_x*i + offset - x[0], 2) / (2 * pow(sigma,  2)));
    }

    double th_prob = 0.5*(1+erf((x[0] - th_level)/(sqrt(2) * th_sigma)));       

    return th_prob * result;
}
TF1 *osc = new TF1("osc_function", osc_function, 0, 200000, 8);

void InitHistograms()
{
    Taking["092"] = {false, true,   true,   true,   true,   true,   true,   true,   true,   true};
    Taking["093"] = {false, true,   false,  true,   false,  false,  true,   false,  false,  false};
    Taking["094"] = {false, false,  true,       false,  true,   false,   true,   false,   true,   false};
    Taking["095"] = {false, false,  false,      false,  false,   false,   false,   false,   false,   true};
    Taking["096"] = {false, false,  false,   false,  false,   false,   false,   false,   false,   true};
    Taking["097"] = {false, true,   true,   true,  true,   true,   true,   true,   true,   true};
    Taking["098"] = {false, true,   true,   true,  true,   true,   true,   true,   true,   true};
    Taking["099"] = {false, true,   false,  false,  false,   false,   false,   false,   false,   false};
    Taking["100"] = {false, false,  true,   false,  false,   false,   false,   false,   false,   false};
    Taking["101"] = {false, false,  false,   true,  false,   false,   false,   false,   false,   false};
    Taking["102"] = {false, false,  false,   false,  true,   false,   false,   false,   false,   false};
    Taking["103"] = {false, false,  false,   false,  false,   true,   false,   false,   false,   false};
    Taking["104"] = {false, false,  false,   false,  false,   false,   true,   false,   false,   false};
    Taking["105"] = {false, false,  false,   false,  false,   false,   false,   true,   false,   false};
    Taking["106"] = {false, false,  false,   false,  false,   false,   false,   false,   true,   false};
    Taking["107"] = {false, false,  false,   false,  false,   false,   false,   false,   false,   true};
    Taking["108"] = {false, true,   true,   true,  true,   true,   true,   true,   true,   true};
    Taking["109"] = {false, true,   true,   true,  true,   true,   true,   true,   true,   true};


    for (string run : Runs)
    {
        for (int i = 1; i < 10; i++)
        {
            H_High[run][i] = new TH1D(("H_High_" + to_string(i) + "_" + run).c_str(), ("H_High_" + to_string(i) + "_" + run).c_str(), 10000, 0, 200000);
            H_High[run][i]->GetXaxis()->SetTitle("Channel");
            H_High[run][i]->GetYaxis()->SetTitle("Counts");
            H_High[run][i]->GetXaxis()->CenterTitle();
            H_High[run][i]->GetYaxis()->CenterTitle();
        }
    }
}

void WriteInFile()
{
    TFile *file = new TFile((DIR_ROOT_DATA_CALIBRATED + "SIPM_Dark_tmp.root").c_str(), "RECREATE");
    for (string run : Runs)
    {
        for (int i = 1; i < 10; i++)
        {
            H_High[run][i]->Write();
        }
    }
    file->Close();
}

void ReadInFile()
{
    TFile *file = new TFile((DIR_ROOT_DATA_CALIBRATED + "SIPM_Dark_tmp.root").c_str(), "READ");
    for (string run : Runs)
    {
        for (int i = 1; i < 10; i++)
        {
            H_High[run][i] = (TH1D *)file->Get(("H_High_" + to_string(i) + "_" + run).c_str());
        }
    }
    // file->Close();
}
void WriteHistograms()
{
    TFile *file = new TFile((DIR_ROOT_DATA_CALIBRATED + "SIPM_Dark.root").c_str(), "RECREATE");
    cout << "Writing histograms to " << DIR_ROOT_DATA_CALIBRATED + "SIPM_Dark.root" << endl;

    TGraph *G_Thresholdi = new TGraph();
    G_Thresholdi->SetTitle("Threshold");
    G_Thresholdi->GetYaxis()->SetTitle("CHANNEL");
    G_Thresholdi->GetXaxis()->SetTitle("SiPM");
    G_Thresholdi->GetXaxis()->CenterTitle();
    G_Thresholdi->GetYaxis()->CenterTitle();

    TGraph *G_ThresholdmV = new TGraph();
    G_ThresholdmV->SetTitle("Threshold");
    G_ThresholdmV->GetYaxis()->SetTitle("mV");
    G_ThresholdmV->GetXaxis()->SetTitle("SiPM");
    G_ThresholdmV->GetXaxis()->CenterTitle();
    G_ThresholdmV->GetYaxis()->CenterTitle();

    TGraph *G_Resolution = new TGraph();
    G_Resolution->SetTitle("Resolution");
    G_Resolution->GetYaxis()->SetTitle("Channel");
    G_Resolution->GetXaxis()->SetTitle("SiPM");
    G_Resolution->GetXaxis()->CenterTitle();
    G_Resolution->GetYaxis()->CenterTitle();

    for (int i = 1; i < 10; i++)
    {
        Info("Processing SiPM " + to_string(i));
        vector<pair<int, TH1D*>> entries;
        TCanvas *cSiPMi = new TCanvas(("cSiPM" + to_string(i)).c_str(), ("cSiPM" + to_string(i)).c_str(), 800, 600);
        TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);
        for (string run : Runs)
        {
            if (!Taking[run][i])
                continue;
            // H_High[run][i]->SetLineColor(stoi(run) - 91);
            // H_High[run][i]->Draw("SAME HIST");
            entries.emplace_back(Thresholds[run][i], H_High[run][i]);
        }
        sort(entries.begin(), entries.end());
        for (const auto &entry : entries)
        {
            entry.second->SetLineColor(static_cast<int>(round(entry.first)));
            entry.second->Draw("SAME HIST");
            legend->AddEntry(entry.second, (to_string(static_cast<int>(round(entry.first))) + " mV").c_str(), "l");
        }
        legend->Draw("SAME");
        cSiPMi->Write();

        /// fitting lowest threshold 
        //th level
        osc->SetParLimits(0, 0, 1e4);
        osc->SetParameter(0, 4e3);
        //th sigma
        osc->SetParLimits(1, 0, 10e3);
        osc->SetParameter(1, 6e3);

        //offset
        osc->SetParLimits(2, -1e4, 1e4);
        osc->SetParameter(2, 0);
        //delta x
        osc->SetParLimits(3, 10e3, 15e3);
        osc->SetParameter(3, 12e3);

        //sigma offset
        osc->SetParLimits(4, 0, 1e4);
        osc->SetParameter(4, 3.7e3);
        //sigma i
        osc->SetParLimits(5, 0, 1e4);
        osc->SetParameter(5, 0);

        // A
        // osc->SetParLimits(6, 0, 5e6);
        // osc->SetParameter(6, 20e3);
        // osc->SetParLimits(7, 0, 50e4);
        // osc->SetParameter(7, 6e3);

        // osc->SetParLimits(6, 0, 1e15);
        osc->SetParameter(6, 4.17e9);
        osc->SetParLimits(7, -5, 0);
        osc->SetParameter(7, -1.87);

        // osc->SetParLimits(8, 0, 5e4);
        // osc->SetParameter(8, 0.9e3);
        // // osc->FixParameter(8, 0.9e3);
        // osc->SetParLimits(9, 0, 50e4);
        // osc->SetParameter(9, 0.1e3);
        // // osc->FixParameter(9, 0.1e3);
        // osc->SetParLimits(10, 0, 50e4);
        // osc->SetParameter(10, 0.1e3);
        // // osc->FixParameter(10, 0.1e3);
        // osc->SetParLimits(11, 0, 50e4);
        // osc->SetParameter(11, 0.1e3);
        // // osc->FixParameter(11, 0.1e3);
        // osc->SetParLimits(12, 0, 50e4);
        // osc->SetParameter(12, 0.1e3);
        // // osc->FixParameter(12, 0.1e3);
        // osc->SetParLimits(13, 0, 50e4);
        // osc->SetParameter(13, 0.1e3);
        // // osc->FixParameter(13, 0.1e3);
        // osc->SetParLimits(14, 0, 50e2);
        // osc->SetParameter(14, 0.1e1);
        // // osc->FixParameter(14, 0.1e1);
        // osc->SetParLimits(15, 0, 50e2);
        // osc->SetParameter(15, 0.1e1);
        // // osc->FixParameter(15, 0.1e1);

                
        osc->SetNpx(10000);
        entries[0].second->Fit(osc, "MULTITHREAD EQ");
        G_Thresholdi->SetPoint(i-1, i, osc->GetParameter(0));
        G_ThresholdmV->SetPoint(i-1, i, entries[0].first);
        G_Resolution->SetPoint(i-1, i, osc->GetParameter(1));


        TCanvas *cSiPMiFit = new TCanvas(("cSiPM" + to_string(i) + "Fit").c_str(), ("cSiPM" + to_string(i) + "Fit").c_str(), 800, 600);
        entries[0].second->Draw();
        // osc->Draw("SAME");
        for (int n = 0 ; n < 10 ; n++)
        {
            TF1 *gaus = new TF1(("gaus" + to_string(n)).c_str(), "gaus", 0, 200000);
            double AA = osc->GetParameter(6) *exp(n * osc->GetParameter(7));
            // double A = osc->GetParameter(6+n);  
            double sigma = osc->GetParameter(4) + n * osc->GetParameter(5);
            gaus->SetParameters(AA/(sqrt(2*M_PI)*sigma), osc->GetParameter(2) + n * osc->GetParameter(3), sigma);
            gaus->SetLineColor(kBlue);
            gaus->Draw("SAME");
        }
        //add on canvas threshold plot

        cSiPMiFit->Write();

        double lower_threshold = osc->GetParameter(0);
        double lower_threshold_mV = entries[0].first;    
        //fitting other threshold with the previous function
        TCanvas *cSiPMiFit_all = new TCanvas(("cSiPM" + to_string(i) + "Fit_all").c_str(), ("cSiPM" + to_string(i) + "Fit_all").c_str(), 800, 600);
        int number = sqrt(entries.size())+1;
        cSiPMiFit_all->Divide(number, number);
        TGraphErrors *G_Threshold = new TGraphErrors();
        G_Threshold->SetTitle("Threshold");
        G_Threshold->GetYaxis()->SetTitle("mV");
        G_Threshold->GetXaxis()->SetTitle("Channel");

        cSiPMiFit_all->cd(1);
        cSiPMiFit_all->SetLogy();
        entries[0].second->Draw();
        G_Threshold->SetPoint(0, lower_threshold/1000, lower_threshold_mV);
        int counter = 2;
        for (const auto &entry : entries)
        {
            if (entry.first == entries[0].first)
                continue;
            cSiPMiFit_all->cd(counter);
            cSiPMiFit_all->SetLogy();
            osc->SetParLimits(0, 0, 2e5);
            osc->SetParameter(0, lower_threshold * entry.first / lower_threshold_mV);
            osc->FixParameter(1, osc->GetParameter(1));

            osc->SetParameter(2, osc->GetParameter(2));
            osc->FixParameter(3, osc->GetParameter(3));
            osc->FixParameter(4, osc->GetParameter(4));
            osc->FixParameter(5, osc->GetParameter(5));
           
            osc->SetParameter(6, osc->GetParameter(6) / entries[0].second->Integral() * entry.second->Integral());
            osc->SetParameter(7, osc->GetParameter(7));
            entry.second->Fit(osc, "MULTITHREAD EQ");
            entry.second->SetTitle((to_string(static_cast<int>(round(entry.first))) + " mV").c_str());
            entry.second->Draw();

            G_Threshold->SetPoint(counter-1, osc->GetParameter(0)/1000, entry.first);
            G_Threshold->SetPointError(counter-1, osc->GetParError(0)/1000 * osc->GetChisquare() / osc->GetNDF(), 0);

            counter++;

            
        }
        cSiPMiFit_all->cd(number*number);
        G_Threshold->SetMarkerStyle(20);
        G_Threshold->SetMarkerSize(1);
        G_Threshold->Fit("pol1", "N");
        G_Threshold->GetYaxis()->SetRangeUser(0, -1111);
        G_Threshold->Draw("AP");
        cSiPMiFit_all->Write();
    }

    string run = "092"; 
    TCanvas *cSiPMi = new TCanvas("cSiPMs", "cSiPM", 800, 600);
    TLegend *legend= new TLegend(0.7, 0.7, 0.9, 0.9);
    double integral;
    for (int i = 1; i < 10; i++)
    {        
        H_High[run][i]->SetLineColor(i);
        if (i == 1)
            integral = H_High[run][i]->Integral();
        H_High[run][i]->Scale(integral / H_High[run][i]->Integral());
        H_High[run][i]->Draw("SAME HIST");
        legend->AddEntry(H_High[run][i], (detectorName[100 + i] + " (" + to_string(static_cast<int>(round(Thresholds[run][i]))) + " mV)").c_str(), "l");
    }
    legend->Draw("SAME");
    cSiPMi->Write();

    TCanvas *cThreshold = new TCanvas("cThreshold", "cThreshold", 800, 600);
    G_Thresholdi->SetMarkerStyle(20);
    G_Thresholdi->SetMarkerSize(1);
    G_Thresholdi->Draw("AP");
    cThreshold->Write();

    TCanvas *cThresholdmV = new TCanvas("cThresholdmV", "cThresholdmV", 800, 600);
    G_ThresholdmV->SetMarkerStyle(20);
    G_ThresholdmV->SetMarkerSize(1);
    G_ThresholdmV->Draw("SAME AP");
    cThresholdmV->Write();

    TCanvas *cResolution = new TCanvas("cResolution", "cResolution", 800, 600);
    G_Resolution->SetMarkerStyle(20);   
    G_Resolution->SetMarkerSize(1);
    G_Resolution->Draw("SAME AP");
    cResolution->Write();
    file->Close();
}




void ReadThreshold(string run)
{
    // replac e.root by .setup in ROOTFILE name
    string filename = DIR_DATA_HDD + "run_" + run + "_lossless_NoBeam.fast/run_" + run + "_lossless_NoBeam.root";
    filename.replace(filename.find(".root"), 5, ".setup");

    ifstream file(filename);
    if (!file.is_open())
    {
        Error("Could not open file " + filename);
    }

    double threshold;

    while (!file.eof())
    {
        string line;
        getline(file, line);
        

        string name;
        string equal;
        double value;
        istringstream iss(line);

        iss >> name >> equal >> value;

        if (name == "Trigger_Threshold_Mode_Threshold")
        {
            threshold = value;
        }

        if (name == "DSP_Label" && IsDetectorBetaHigh((int)value))
        {
            Thresholds[run][GetDetectorChannel((int)value)] = threshold;
        }
    }
}


#endif