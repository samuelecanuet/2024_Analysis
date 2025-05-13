#include "Merger.hh"
#include <ctime>

int main(int argc, char *argv[])
{
    FLAG2021 = true;
    InitDetectors("Config_Files/sample.pid");
    InitRuns();

    MATCHED_File = MyTFile((DIR_ROOT_DATA_MATCHED + "matched.root").c_str(), "READ");

    //////////////////// MATCHING REAR STRIP ////////////////////
    // TFile *File_Fit = MyTFile((DIR_ROOT_DATA_MATCHED + "RearStrip_Fit.root").c_str(), "RECREATE");
    // File_Fit->cd();

    ////////////////////////////////////////////////////////////
    // using correction of run gaindrifting and correction from strip to strip and merge all the runs in a single file for each nucleus
    //////////////////// MERGING ////////////////////
    Info("Merging");
    for (auto &pairr : Map_RunFiles)
    {
        string NUCLEUS = pairr.first;
        if (NUCLEUS.find("32Ar") != string::npos && NUCLEUS.find("33Ar") != string::npos)
            continue;
        pair<string, vector<string>> NUCLEUS_Run = make_pair(NUCLEUS, Map_RunFiles[NUCLEUS]);
        Start("Merging " + NUCLEUS_Run.first);
        MERGED_File = MyTFile((DIR_ROOT_DATA_MERGED + NUCLEUS_Run.first + "_merged.root").c_str(), "RECREATE");
        InitGraph();
        TDirectory *dir = MERGED_File->mkdir("Strips");
        TTree *MERGED_Tree = new TTree("MERGED_Tree", "MERGED_Tree");
        MERGED_Tree->Branch("MERGED_Tree_Silicon", &MERGED_Tree_Silicon);
        MERGED_Tree->Branch("MERGED_Tree_SiPMGroup", &MERGED_Tree_SiPMGroup);
        MERGED_Tree->Branch("MERGED_Tree_HRS", &MERGED_Tree_HRS);

        TTree *MERGED_Tree_Detector[SIGNAL_MAX];

        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorSiliStrip(i))
            {
                MERGED_Tree_Detector[i] = new TTree(("MERGED_Tree_" + detectorName[i]).c_str(), ("MERGED_Tree_" + detectorName[i]).c_str());
                MERGED_Tree_Detector[i]->Branch("Channel", &Channel);
            }
        }

        for (int i = 0; i < NUCLEUS_Run.second.size(); i++)
        {
            string Run = NUCLEUS_Run.second[i];
            Info("Current Run : " + Run);
            GROUPED_filename = SearchFiles(DIR_ROOT_DATA_GROUPED, Run);
            GROUPED_File = MyTFile((DIR_ROOT_DATA_GROUPED + GROUPED_filename).c_str(), "READ");
            if (GROUPED_File == NULL)
            {
                continue;
            }

            TTree *Tree = (TTree *)GROUPED_File->Get("CLEANED_Tree");
            if (Tree == NULL)
            {
                Warning("No Tree found in " + GROUPED_filename);
                continue;
            }

            LoadMatchingFunction(atoi(Run.c_str()));
            TTreeReader *Reader = new TTreeReader(Tree);
            Silicon = new TTreeReaderArray<Signal>(*Reader, "CLEANED_Tree_Silicon");
            SiPM_Groups = new TTreeReaderValue<vector<vector<pair<Signal, Signal>>>>(*Reader, "CLEANED_Tree_SiPMGroup");
            if (YEAR == 2025)
                HRS = new TTreeReaderValue<Signal>(*Reader, "CLEANED_Tree_HRS");

            // MAKING THE SEPARATION BETWEEN PRE AND POST GAIN CHANGE IN SIPMS
            // if (atoi(Run.c_str()) == 112 && YEAR == 2024)
            // {
            //     MERGED_Tree_Silicon = {Signal(), Signal()};
            //     MERGED_Tree_SiPMGroup = {{make_pair(Signal(), Signal())}};
            //     MERGED_Tree->Fill();
            //     MERGED_Tree_Silicon.clear();
            //     MERGED_Tree_SiPMGroup.clear();
            // }
            //////////////////////////////

            clock_t start = clock(), Current;
            int Entries = Reader->GetEntries();
            while (Reader->Next())
            {
                ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "");

                MERGED_Tree_HRS = Signal();
                MERGED_Tree_SiPMGroup = vector<vector<pair<Signal, Signal>>>();
                MERGED_Tree_Silicon = vector<Signal>();
                if (YEAR == 2025)
                {
                    if ((**HRS).isValid)
                    {
                        MERGED_Tree_HRS = **HRS;
                        MERGED_Tree->Fill();
                        continue;
                    }
                }

                MERGED_Tree_Silicon.push_back((*Silicon)[0]);

                // run matching correction
                (*Silicon)[1].Channel = Matching_function[(*Silicon)[1].Label]->Eval((*Silicon)[1].Channel);

                MERGED_Tree_Silicon.push_back((*Silicon)[1]);
                Channel = (*Silicon)[1].Channel;

                // !!! ADDING SIPM MATCHING !!! //
                MERGED_Tree_SiPMGroup = **SiPM_Groups;

                MERGED_Tree->Fill();
                MERGED_Tree_Detector[(*Silicon)[1].Label]->Fill();

                MERGED_Tree_Silicon.clear();
                MERGED_Tree_SiPMGroup.clear();
                Channel = 0;
            }
        }

        Info("Merging done");

        MERGED_File->cd();
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorSiliStrip(i))
            {
                // Info("Writing Tree " + detectorName[i], 1);
                MERGED_File->cd();
                MERGED_Tree_Detector[i]->Write();
            }
        }

        MERGED_File->cd();
        MERGED_Tree->Write();
        MERGED_File->Close();
    }

    return 0;
}