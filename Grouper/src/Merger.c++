#include "Merger.hh"
#include <ctime>

int main(int argc, char *argv[])
{

    InitDetectors("Config_Files/sample.pid");
    Init();
    
    MATCHED_File = new TFile((DIR_ROOT_DATA_MATCHED + "matched.root").c_str(), "READ");
    
    //////////////////// MATCHING REAR STRIP ////////////////////
    TFile *File_Fit = new TFile((DIR_ROOT_DATA_MATCHED + "RearStrip_Fit.root").c_str(), "RECREATE");
    File_Fit->cd();
    
    Info("Matching Rear Strip");
    // read all run and all nuclei to fill a grph and then fit with a linear fucntion the gain atching betweens strip and Rear
    // for (string NUCLEUS : NUCLEI)
    // {
    //     pair<string, vector<string>> NUCLEUS_Run = make_pair(NUCLEUS, Map_RunFiles[NUCLEUS]);
    //     for (int i = 0; i < NUCLEUS_Run.second.size(); i++)
    //     {
    //         string Run = NUCLEUS_Run.second[i];
    //         Info("Current Run : " + Run);
    //         GROUPED_filename = SearchFiles(DIR_ROOT_DATA_GROUPED, Run);
    //         GROUPED_File = new TFile((DIR_ROOT_DATA_GROUPED + GROUPED_filename).c_str(), "READ");
        
    //         TTree *Tree = (TTree *)GROUPED_File->Get("CLEANED_Tree");
    //         if (Tree == NULL)
    //         {
    //             Warning("No Tree found in " + GROUPED_filename);
    //             continue;
    //         }

    //         LoadMatchingFunction(atoi(Run.c_str()));
            
    //         TTreeReader *Reader = new TTreeReader(Tree);
    //         Silicon = new TTreeReaderArray<Signal>(*Reader, "CLEANED_Tree_Silicon");
    //         SiPM_High = new TTreeReaderArray<Signal>(*Reader, "CLEANED_Tree_SiPMHigh");
    //         SiPM_Low = new TTreeReaderArray<Signal>(*Reader, "CLEANED_Tree_SiPMLow");

    //         clock_t start = clock(), Current;
    //         int Entries = Reader->GetEntries();
    //         while (Reader->Next())
    //         {
    //             ProgressBar(Reader->GetCurrentEntry(), Entries, start, Current, "");

    //             (*Silicon)[1].Channel = Matching_function[(*Silicon)[1].Label]->Eval((*Silicon)[1].Channel);
    //             G_RearStrip[(*Silicon)[1].Label]->SetPoint(counter_graph[(*Silicon)[1].Label], (*Silicon)[0].Channel, (*Silicon)[1].Channel);
    //             counter_graph[(*Silicon)[1].Label]++;
    //         }
    //     }
    // }

    // fit the graph with linear fucntion, save the fucntion/graph in a root file and compute in MatchingRearStrip the gain match fucntion to correct strip 2,3,4,5 on 1
    // File_Fit->cd();
    // for (int i = 0; i < SIGNAL_MAX; i++)
    // {
    //     if (IsDetectorSiliStrip(i))
    //     {
    //         TCanvas *c = new TCanvas(("RearStrip_" + detectorName[i]).c_str(), ("RearStrip_" + detectorName[i]).c_str(), 1920, 1080);
    //         TF1 *f = new TF1(("pol1_" + detectorName[i]).c_str(), "[0]*x", 0, 200000);
    //         G_RearStrip[i]->Fit(f, "QR");
    //         G_RearStrip[i]->Draw("AP");
    //         c->Write();     


    //         // compute fucntion for matching rear strip for all strip except 1 (took as reference)
    //         if (GetDetectorChannel(i) != 1)
    //         {
    //             // formula : 
    //             // E_ref =   a1 * x
    //             // E     =   a2 * x
    //             // E_corrected = a1/a2 * x

    //             double a1 = G_RearStrip[GetDetector(i)*10+1]->GetFunction(("pol1_D" + to_string(GetDetector(i)) + ".1").c_str())->GetParameter(0);
    //             double a2 = G_RearStrip[i]->GetFunction(("pol1_" + detectorName[i]).c_str())->GetParameter(0);

    //             double a_new = a1/a2;

    //             MatchingRearStrip[i] = new TF1(("MatchingRearStrip_" + detectorName[i]).c_str(), "[0]*x", 0, 200000);
    //             MatchingRearStrip[i]->SetParameter(0, a_new);
                
    //         }   
    //         else
    //         {
    //             MatchingRearStrip[i] = new TF1(("MatchingRearStrip_" + detectorName[i]).c_str(), "[0]*x", 0, 200000);
    //             MatchingRearStrip[i]->SetParameter(0, 1);
    //         }   

    //         // // compute fucntion for matching rear strip for all strip except 1 (took as reference)
    //         // if (GetDetectorChannel(i) != 1)
    //         // {
    //         //     // formula : 
    //         //     // E_ref =   a1 * x + b1
    //         //     // E     =   a2 * x + b2
    //         //     // E_corrected = a1/a2 * x - a1/a2 * b2 + b1

    //         //     double a1 = G_RearStrip[GetDetector(i)*10+1]->GetFunction("pol1")->GetParameter(1);
    //         //     double b1 = G_RearStrip[GetDetector(i)*10+1]->GetFunction("pol1")->GetParameter(0);
    //         //     double a2 = G_RearStrip[i]->GetFunction("pol1")->GetParameter(1);
    //         //     double b2 = G_RearStrip[i]->GetFunction("pol1")->GetParameter(0);

    //         //     double a_new = a1/a2;
    //         //     double b_new = -a1/a2 * b2 + b1;

    //         //     MatchingRearStrip[i] = new TF1(("MatchingRearStrip_" + detectorName[i]).c_str(), "pol1", 0, 200000);
    //         //     MatchingRearStrip[i]->SetParameter(0, b_new);
    //         //     MatchingRearStrip[i]->SetParameter(1, a_new);
                
    //         // }   
    //         // else
    //         // {
    //         //     MatchingRearStrip[i] = new TF1(("MatchingRearStrip_" + detectorName[i]).c_str(), "pol1", 0, 200000);
    //         //     MatchingRearStrip[i]->SetParameters(0, 1);
    //         // }   
    //     }

    //     if (IsDetectorSiliBack(i))
    //     {
    //         for (int strip = 1; strip <= 5; strip++)
    //         {
    //             delete G_RearStrip[GetDetector(i)*10+strip];
    //         }
    //     }
    // }


    ////////////////////////////////////////////////////////////
    // using correction of run gaindrifting and correction from strip to strip and merge all the runs in a single file for each nucleus
    //////////////////// MERGING ////////////////////
    Info("Merging");
    for (string NUCLEUS : NUCLEI)
    {
        pair<string, vector<string>> NUCLEUS_Run = make_pair(NUCLEUS, Map_RunFiles[NUCLEUS]);
        Start("Merging " + NUCLEUS_Run.first);
        MERGED_File = new TFile((DIR_ROOT_DATA_MERGED + NUCLEUS_Run.first + "_merged.root").c_str(), "RECREATE");
        InitGraph();
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
                
                // run matching correction
                // (*Silicon)[1].Channel = Matching_function[(*Silicon)[1].Label]->Eval((*Silicon)[1].Channel);
                H_RearStrip[(*Silicon)[0].Label]->Fill((*Silicon)[0].Channel, (*Silicon)[1].Channel);
                
                // rear-strip matching correction
                // (*Silicon)[1].Channel = MatchingRearStrip[(*Silicon)[1].Label]->Eval((*Silicon)[1].Channel);
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

                // H_Strip_Matched[(*Silicon)[1].Label]->Fill(MatchingRearStrip[(*Silicon)[1].Label]->Eval((*Silicon)[1].Channel));
                // H_RearStrip_Matched[(*Silicon)[0].Label]->Fill((*Silicon)[0].Channel, MatchingRearStrip[(*Silicon)[1].Label]->Eval((*Silicon)[1].Channel));
            }

            /// adding histograms for fake coincidences
            if (Run == "057")
            {
                for (int mul = 1; mul <= 9; mul++)
                {
                    H_SiPM_Full[mul] = (TH1D *)GROUPED_File->Get(("Coincidences/RearSiPM_Time_New_Nearest_M" + to_string(mul)).c_str());
                }
            }
            else
            {
                for (int mul = 1; mul <= 9; mul++)
                {
                    TH1D * H = (TH1D *)GROUPED_File->Get(("Coincidences/RearSiPM_Time_New_Nearest_M" + to_string(mul)).c_str())->Clone();
                    H_SiPM_Full[mul]->Add((TH1D*)H->Clone());
                }
            }

            H_Sum->Add((TH1D*)GROUPED_File->Get("Coincidences/RearSiPM_Time_New_Nearest_M9")->Clone());
        }
        
        MERGED_File->cd();
        H_Sum->Write();
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorSiliStrip(i))
            {
                MERGED_File->cd();
                MERGED_Tree_Detector[i]->Write();
            }
        }

       

       if (NUCLEUS == "32Ar")
       {
           ComputeFakeCoincidences();

           //writting rear and strip
           File_Fit->cd();
           for (int i = 0; i < SIGNAL_MAX; i++)
           {
               if (IsDetectorSiliBack(i))
               {
                   H_RearStrip[i]->Write();
                   H_RearStrip_Matched[i]->Write();

                   TCanvas *c = new TCanvas(("D" + to_string(GetDetector(i))).c_str(), ("D" + to_string(GetDetector(i))).c_str(), 1920, 1080);
                   TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);
                   for (int strip = 1; strip <= 5; strip++)
                   {
                       H_Strip_Matched[GetDetector(i) * 10 + strip]->SetLineColor(strip);
                       H_Strip_Matched[GetDetector(i) * 10 + strip]->Draw("HIST SAME");
                       legend->AddEntry(H_Strip_Matched[GetDetector(i) * 10 + strip], ("Strip " + to_string(strip)).c_str(), "l");
                   }
                   legend->Draw("SAME");
                   c->Write();
               }
           }
           File_Fit->Close();
       }
    

       MERGED_File->cd();
       MERGED_Tree->Write();
       MERGED_File->Close();
    }  
    

    // File_Fit->cd();
    // for (int i = 0; i < SIGNAL_MAX; i++)
    // {
    //     if (IsDetectorSiliBack(i))
    //     {
    //         H_RearStrip[i]->Write();
    //         H_RearStrip_Matched[i]->Write();

    //         TCanvas *c = new TCanvas(("D" + to_string(GetDetector(i))).c_str(), ("D" + to_string(GetDetector(i))).c_str(), 1920, 1080);
    //         TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);
    //         for (int strip = 1; strip <= 5; strip ++)
    //         {
    //             H_Strip_Matched[GetDetector(i)*10+strip]->SetLineColor(strip);
    //             H_Strip_Matched[GetDetector(i)*10+strip]->Draw("HIST SAME");
    //             legend->AddEntry(H_Strip_Matched[GetDetector(i)*10+strip], ("Strip " + to_string(strip)).c_str(), "l");
    //         }
    //         legend->Draw("SAME");
    //         c->Write();
    //     }
    // }
    // File_Fit->Close();

    return 0;
}