#include "Merger.hh"
#include <ctime>

int main(int argc, char *argv[])
{

    InitDetectors("Config_Files/sample.pid");
    Init();
    MATCHED_File = new TFile((DIR_ROOT_DATA_MATCHED + "matched.root").c_str(), "READ");

    for (string NUCLEUS : NUCLEI)
    {
        pair<string, vector<string>> NUCLEUS_Run = make_pair(NUCLEUS, Map_RunFiles[NUCLEUS]);
        Start("Merging " + NUCLEUS_Run.first);
        MERGED_File = new TFile((DIR_ROOT_DATA_MERGED + NUCLEUS_Run.first + "_merged.root").c_str(), "RECREATE");
        TDirectory *dir = MERGED_File->mkdir("Strips");

        TTree *MERGED_Tree = new TTree("MERGED_Tree", "MERGED_Tree");
        MERGED_Tree->Branch("MERGED_Tree_Silicon", &MERGED_Tree_Silicon);
        MERGED_Tree->Branch("MERGED_Tree_SiPMHigh", &MERGED_Tree_SiPMHigh);
        MERGED_Tree->Branch("MERGED_Tree_SiPMLow", &MERGED_Tree_SiPMLow);

        TTree* MERGED_Tree_Detector[SIGNAL_MAX];

        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorSiliStrip(i))
            {
                H[i] = new TH1D(detectorName[i].c_str(), detectorName[i].c_str(), 10000, 0, 100000);
                MERGED_Tree_Detector[i] = new TTree(("MERGED_Tree_"+detectorName[i]).c_str(), ("MERGED_Tree_"+detectorName[i]).c_str());
                MERGED_Tree_Detector[i]->Branch("Channel", &Channel);
            }
        }

        for (int i = 0; i < NUCLEUS_Run.second.size(); i++)
        {
            string Run = NUCLEUS_Run.second[i];
            Info("Current Run : " + Run);
            GROUPED_filename = SearchFiles(DIR_ROOT_DATA_GROUPED, Run);
            GROUPED_File = new TFile((DIR_ROOT_DATA_GROUPED + GROUPED_filename).c_str(), "READ");
        
            TTree *Tree = (TTree *)GROUPED_File->Get("CLEANED_Tree");
            if (Tree == NULL)
            {
                Warning("No Tree found in " + GROUPED_filename);
                continue;
            }

            LoadMatchingFunction(atoi(Run.c_str()));
            
            TTreeReader *Reader = new TTreeReader(Tree);
            Silicon = new TTreeReaderArray<Signal>(*Reader, "CLEANED_Tree_Silicon");
            SiPM_High = new TTreeReaderArray<Signal>(*Reader, "CLEANED_Tree_SiPMHigh");
            SiPM_Low = new TTreeReaderArray<Signal>(*Reader, "CLEANED_Tree_SiPMLow");

            clock_t start = clock(), Current;
            int Entries = Reader->GetEntries();
            while (Reader->Next())
            {
                ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "");

                MERGED_Tree_Silicon.push_back((*Silicon)[0]);


                (*Silicon)[1].Channel = Matching_function[(*Silicon)[1].Label]->Eval((*Silicon)[1].Channel);
                MERGED_Tree_Silicon.push_back((*Silicon)[1]);
                Channel = (*Silicon)[1].Channel;

                for (Signal s : *SiPM_High)
                    MERGED_Tree_SiPMHigh.push_back(s);
                for (Signal s : *SiPM_Low)
                    MERGED_Tree_SiPMLow.push_back(s);


                MERGED_Tree->Fill();
                MERGED_Tree_Detector[(*Silicon)[1].Label]->Fill();
                
                MERGED_Tree_Silicon.clear();
                MERGED_Tree_SiPMHigh.clear();
                MERGED_Tree_SiPMLow.clear();
                Channel = 0;

                H[(*Silicon)[1].Label]->Fill((*Silicon)[1].Channel);
            }
        }

        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorSiliStrip(i))
            {
                dir_det[i] = dir->mkdir(detectorName[i].c_str());
                dir_det[i]->cd();
                H[i]->Write();
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