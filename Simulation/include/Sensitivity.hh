#include "../../Grouper/include/Detectors.hh"
#include "../../Grouper/include/Analyser.hh"
#include "TKey.h"

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

pair<double, double> ComputeEshift(TFile *file, int det)
{
    TCanvas *c = nullptr;
    TH1D *h_single = nullptr, *h_coinc = nullptr, *h_nocoinc = nullptr;
    if (string(file->GetName()).find("analysed_new") != std::string::npos)
    {
        c = (TCanvas *)file->Get(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_Coinc").c_str());
        h_single = (TH1D *)c->GetPrimitive(("Silicon_Detector_Energy_Deposit_SINGLE_" + detectorName[det] + "_All").c_str());
        h_coinc = (TH1D *)c->GetPrimitive(("Silicon_Detector_Energy_Deposit_COINC_" + detectorName[det] + "_All").c_str());
        h_nocoinc = (TH1D *)c->GetPrimitive(("Silicon_Detector_Energy_Deposit_NOCOINC_" + detectorName[det] + "_All").c_str());
    }
    else if (string(file->GetName()).find("result") != std::string::npos)
    {
        c = (TCanvas *)file->Get(("H_Single_" + detectorName[det]).c_str());
        h_single = (TH1D *)c->GetPrimitive(("H_Single_" + detectorName[det]).c_str());
        h_coinc = (TH1D *)c->GetPrimitive(("H_Coinc_" + detectorName[det]).c_str());
        h_nocoinc = (TH1D *)c->GetPrimitive(("H_NoCoinc_" + detectorName[det]).c_str());
    }
    else
    {
        Error("Unknown file format for file: " + string(file->GetName()));
    }

    if (h_single == nullptr || h_coinc == nullptr || h_nocoinc == nullptr)
    {
        Error("Histograms not found for detector: " + detectorName[det]);
    }

    return ComputeEshift(det, h_single, h_coinc, h_nocoinc);
}