#include "TFile.h"
#include "../../Grouper/include/Detectors.hh"
#include "TCanvas.h"
#include "TH1D.h"
#include <dirent.h>
#include <sys/stat.h>
#include <cstring>
#include <string>
#include <vector>
#include <iostream> 


void Search_files(const std::string &dir_search, vector<std::string> filename_base, std::vector<std::string> &filenames)
{
    // Info("Searching for files in directory: " + dir_search + " with base name: " + filename_base, 1);
    DIR *dir = opendir(dir_search.c_str());
    if (!dir)
    {
        Error("Error opening directory: " + dir_search);
        return;
    }

    struct dirent *entry;
    while ((entry = readdir(dir)) != nullptr)
    {
        std::string name = entry->d_name;

        // cout << "Checking file: " << name << endl;

        // Skip '.' and '..'
        // if (name == "." || name == "..") continue;

        std::string full_path = dir_search + "/" + name;

        // cout << "Full path: " << full_path << endl;

        struct stat file_stat;
        if (stat(full_path.c_str(), &file_stat) == 0 && S_ISREG(file_stat.st_mode))
        {
            bool found = false;
            for (const auto &filename_base : filename_base)
            {
                if (name.find(filename_base) != std::string::npos && name.find(".root") != std::string::npos)
                {
                    found = true;
                }
                else 
                {
                    found = false;
                    break;
                }
            }

            if (found)
            {
                // Info("Found file: " + full_path, 1);
                filenames.push_back(full_path);
            }
        }
    }

    closedir(dir);
}


int CompareSilicon()
{
    FLAG2025 = true;
    InitDetectors("../../Grouper/Config_Files/sample.pid");    
    vector<std::string> filenames;
    Search_files("/run/media/local1/DATANEX/Samuel-G4/Systematics/Cuts/StepMax/", {"_analysed_new.root"}, filenames);

    TFile *fout = new TFile("/run/media/local1/DATANEX/Samuel-G4/Systematics/Cuts/StepMax/Comparison_StepMax.root", "RECREATE");

    TCanvas *c[SIGNAL_MAX];
    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (!IsDetectorSiliStrip(det))
            continue;
        c[det] = new TCanvas(("c_" + detectorName[det]).c_str(), ("Comparison " + detectorName[det]).c_str(), 800, 600);
    }

    TCanvas *cdown = new TCanvas("c_down", "Comparison down", 800, 600);
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    

    int counter = 0;

    for (auto &filename : filenames)
    {
        if (filename != filenames[filenames.size() - 1] && filename != filenames[0])
            continue;
        cout << "Processing file: " << filename << endl;
        counter++;
        TFile *file = TFile::Open(filename.c_str(), "READ");
        if (!file || file->IsZombie())
        {
            Error("Error opening file: " + filename);
            continue;
        }

        TH1D* hdown = nullptr;

        double cuts = stod(filename.substr(filename.find("_StepMax") + 8, filename.find("_analysed") - (filename.find("_StepMax") + 8)));
        if (cuts == 0)
            cuts = 1e-4;

        for (int det = 0; det < SIGNAL_MAX; det++)
        {
            if (!IsDetectorSiliStrip(det))
                continue;
            // TCanvas *cfile = (TCanvas*)file->Get(("H_Single_" + detectorName[det]).c_str());
            // if (!cfile)
            // {
            //     Error("Canvas 'c_" + detectorName[det] + "' not found in file: " + filename);
            //     continue;
            // }

            // get hist 
            // TH1D *hfile = (TH1D*)cfile->GetPrimitive(("H_Single_" + detectorName[det]).c_str());
            TH1D *hfile = (TH1D*)file->Get(("p/Silicon_Detector/D." + to_string(GetDetector(det)) + "/Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_p").c_str());
            if (!hfile)
            {
                Error("Histogram 'H_Single_" + detectorName[det] + "' not found in canvas: ");
                continue;
            }
            hfile->Rebin(2);

            c[det]->cd();
            hfile->SetLineColor(counter);

            if (counter == 1)
                hfile->Draw("HIST");
            else
                hfile->Draw("SAME");
            

            if (det > 50)
            {
                if (hdown == nullptr)
                {
                    // hdown = (TH1D*)cfile->GetPrimitive(("H_Single_" + detectorName[det]).c_str());
                    hdown = (TH1D*)file->Get(("p/Silicon_Detector/D." + to_string(GetDetector(det)) + "/Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_p").c_str());
                    if (!hdown)
                    {
                        Error("Histogram 'H_Down' not found in canvas");
                        continue;
                    }
                }
                else
                {
                    cout << "Adding histogram from file: " << filename << " to hdown" << endl;
                    hdown->Add((TH1D*)file->Get(("p/Silicon_Detector/D." + to_string(GetDetector(det)) + "/Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_p").c_str()));
                }
            }
        }

        cdown->cd();
        hdown->SetLineColor(counter);
        hdown->Draw("SAME");
        legend->AddEntry(hdown, ("Cuts = " + to_string(cuts)).c_str(), "l");
        // file->Close();

        cout << "Finished processing file: " << filename << endl;

    }

    fout->cd();
    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (!IsDetectorSiliStrip(det))
            continue;
        c[det]->Write(("Comparison_" + detectorName[det]).c_str());
    }
    

    // legend->Draw("SAME");
    // cdown->Draw();

    fout->Close();
    return 1; 
}