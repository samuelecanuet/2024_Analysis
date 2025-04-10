#include "Grouper.hh"
#include <ctime>

int main(int argc, char *argv[])
{
    InitDetectors("Config_Files/sample.pid");
    string Run_string;

    if (argc < 2)
    {
        Error("No Run Number Given");
    }
    else
    {
        // option
        if (argc == 3)
        {
            if (string(argv[2]) == "full")
            {
                FULL = true;
            }
        }


        // Run number
        Run = atoi(argv[1]);
        if (Run < 10)
            Run_string = "00" + to_string(Run);
        else if (Run < 100)
            Run_string = "0" + to_string(Run);
        else
            Run_string = to_string(Run);
    }

    Info("Current Run : " + Run_string);
    ///////////////////////////////////  INPUT ///////////////////////////////////
    // DIR_ROOT_DATA = "../../../../../run/media/local1/T7/Samuel/Regrouped/ROOT/";
    ROOT_filename = SearchFiles(DIR_ROOT_DATA, Run_string);
    ROOT_basefilename = ROOT_filename.substr(0, ROOT_filename.find(".root"));

    TFile *ROOT_File = MyTFile((DIR_ROOT_DATA + ROOT_filename).c_str(), "READ");
    TTree *Tree = (TTree *)ROOT_File->Get("Tree_Group");
    TTreeReader *Reader = new TTreeReader(Tree);
    TTreeReaderArray<Signal> signals(*Reader, "Signal");
    // TTreeReaderValue<Signal> signals(*Reader, "Signal");

    ///////////////////////////////////  Grouped ///////////////////////////////////
    if (FULL)
        GROUPED_File = MyTFile((DIR_ROOT_DATA_GROUPED + ROOT_basefilename + "_grouped_full.root").c_str(), "RECREATE");
    else
        GROUPED_File = MyTFile((DIR_ROOT_DATA_GROUPED + ROOT_basefilename + "_grouped.root").c_str(), "RECREATE");

    // GROUPED_File = MyTFile("../../../../../../../mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/run_054_lossless_grouped.root", "RECREATE");
    GROUPED_File->cd();
    WriteTime(ROOT_File, GROUPED_File);
    ///////////////////////////////////  INITIALISATION ///////////////////////////////////
    
    InitHistograms_Grouped();
    InitCalibration();

    clock_t start = clock(), Current;
    if (FULL) ///////// SAVING FIT PARAMETER IF FULL MODE ELSE LOADING THEM
    {
    
        /////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////  SELECTING TREE /////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////
        GROUPED_Tree = new TTree("GROUPED_Tree", "GROUPED_Tree");
        GROUPED_Tree->Branch("GROUPED_Tree_Silicon", &GROUPED_Tree_Silicon);
        GROUPED_Tree->Branch("GROUPED_Tree_SiPMHigh", &GROUPED_Tree_SiPMHigh);
        GROUPED_Tree->Branch("GROUPED_Tree_SiPMLow", &GROUPED_Tree_SiPMLow);
        ULong64_t TotalEntries = Reader->GetEntries();
        while (Reader->Next())
        {
            ProgressBar(Reader->GetCurrentEntry(), TotalEntries, start, Current, "Selecting Groups : ");

            // for (int i = 0; i < signals.GetSize(); i++)
            //     Verbose(signals[i], VERBOSE, 2);

            SearchForCoincidence(signals);
        }

        WriteHistograms_Grouped();
        /////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////  CUTTING GROUPS /////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////
        start = clock(), Current;

        CUTTED_Tree = new TTree("CUTTED_Tree", "CUTTED_Tree");
        CUTTED_Tree->Branch("CUTTED_Tree_Silicon", &CUTTED_Tree_Silicon);
        CUTTED_Tree->Branch("CUTTED_Tree_SiPMHigh", &CUTTED_Tree_SiPMHigh);
        CUTTED_Tree->Branch("CUTTED_Tree_SiPMLow", &CUTTED_Tree_SiPMLow);

        Reader = new TTreeReader(GROUPED_Tree);
        Silicon = new TTreeReaderArray<Signal>(*Reader, "GROUPED_Tree_Silicon");
        SiPM_High = new TTreeReaderArray<Signal>(*Reader, "GROUPED_Tree_SiPMHigh");
        SiPM_Low = new TTreeReaderArray<Signal>(*Reader, "GROUPED_Tree_SiPMLow");

        if (Run != 114 || Run != 37) LoadFitParameters();
        Reader->Restart();
        TotalEntries = Reader->GetEntries();
        while (Reader->Next())
        {
            ProgressBar(Reader->GetCurrentEntry(), TotalEntries, start, Current, "Cutting Groups : ");

            // for (int i = 0; i < signals.GetSize(); i++)
            //     Verbose(signals[i], VERBOSE, 2);

            CuttingGroups();
        }

        WriteHistograms_Cutted();

        /////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////// SILICON WALK GROUPS /////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////
        start = clock(), Current;

        Reader = new TTreeReader(CUTTED_Tree);
        Silicon = new TTreeReaderArray<Signal>(*Reader, "CUTTED_Tree_Silicon");
        SiPM_High = new TTreeReaderArray<Signal>(*Reader, "CUTTED_Tree_SiPMHigh");
        SiPM_Low = new TTreeReaderArray<Signal>(*Reader, "CUTTED_Tree_SiPMLow");

        if (Run != 114 || Run != 37) LoadFitParameters();
        Reader->Restart();
        TotalEntries = Reader->GetEntries();
        while (Reader->Next())
        {
            ProgressBar(Reader->GetCurrentEntry(), TotalEntries, start, Current, "Removing Silicon walk : ");

            // for (int i = 0; i < signals.GetSize(); i++)
            //     Verbose(signals[i], VERBOSE, 2);

            SiliconWalkCorrection();
        }

        WriteHistograms_SiliconWalk();
        /////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////// SiPM WALK GROUPS ////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////
        start = clock(), Current;
        if (Run != 114 || Run != 37) LoadFitParameters();
        Reader->Restart();
        while (Reader->Next())
        {
            ProgressBar(Reader->GetCurrentEntry(), TotalEntries, start, Current, "Removing SiPM walk : ");

            // for (int i = 0; i < signals.GetSize(); i++)
            //     Verbose(signals[i], VERBOSE, 2);

            SiPMWalkCorrection();
        }

        WriteHistograms_SiPMWalk();
        /////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////  FINAL CLEANING /////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        

        if (Run == 114 || Run == 37) SaveFitParameters();
    }
    else
    {
        LoadFitParameters();
    }
    ///////////////////////////////////////////////////////////
    start = clock(), Current;
    CLEANED_Tree = new TTree("CLEANED_Tree", "CLEANED_Tree");
    CLEANED_Tree->Branch("CLEANED_Tree_Silicon", &CLEANED_Tree_Silicon);
    CLEANED_Tree->Branch("CLEANED_Tree_SiPMGroup", &CLEANED_Tree_SiPMGroup);

    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            CLEANED_Tree_detector[i] = new TTree(("CLEANED_Tree_" + detectorName[i]).c_str(), ("CLEANED_Tree_" + detectorName[i]).c_str());
            CLEANED_Tree_detector[i]->Branch("Channel", &Tree_Channel_detector);
        }
    }

    Reader = new TTreeReader(Tree);
    signals = TTreeReaderArray<Signal>(*Reader, "Signal");
    Reader->Restart();
    int TotalEntries = Reader->GetEntries();
    while (Reader->Next())
    {
        ProgressBar(Reader->GetCurrentEntry(), TotalEntries, start, Current, "Final Cleaning : ");

        // for (int i = 0; i < signals.GetSize(); i++)
        //     cout << signals[i] << endl;

        CleaningGroups(signals, 0);
    }

    // double past = -300;
    // int future = 350;

    // double start_gate = -50;
    // double stop_gate = 50;  

    // TH1D* H_Time_SIPMGROUPING_M[SIGNAL_MAX][BETA_SIZE+1];
    // TH1D* H_Time_SIPMGROUPING_M_equal[SIGNAL_MAX][BETA_SIZE+1];
    // TH1D* H_Channel_SIPMGROUPING_M[SIGNAL_MAX][BETA_SIZE+1];
    // TH1D* H_Channel_SIPMGROUPING_Mequal[SIGNAL_MAX][BETA_SIZE+1];

    // for (int i = 0; i < SIGNAL_MAX; i++)
    // {
    //     if (IsDetectorBetaHigh(i))
    //     {
    //         for (int m = 1; m <= BETA_SIZE; m++)
    //         {
    //             H_Time_SIPMGROUPING_M[i][m] = new TH1D(("H_Time_SIPMGROUPING_" + detectorName[i] + "_M>=" + to_string(m)).c_str(), ("H_Time_SIPMGROUPING_" + detectorName[i] + "_M>=" + to_string(m)).c_str(), 500, -500, 500);
    //             H_Time_SIPMGROUPING_M_equal[i][m] = new TH1D(("H_Time_SIPMGROUPING_" + detectorName[i] + "_M=" + to_string(m)).c_str(), ("H_Time_SIPMGROUPING_" + detectorName[i] + "_M=" + to_string(m)).c_str(), 500, -500, 500);

    //             H_Channel_SIPMGROUPING_M[i][m] = new TH1D(("H_Channel_SIPMGROUPING_" + detectorName[i] + "_M>=" + to_string(m)).c_str(), ("H_Channel_SIPMGROUPING_" + detectorName[i] + "_M>=" + to_string(m)).c_str(), eHighN, eHighMin, eHighMax);
    //             H_Channel_SIPMGROUPING_Mequal[i][m] = new TH1D(("H_Channel_SIPMGROUPING_" + detectorName[i] + "_M=" + to_string(m)).c_str(), ("H_Channel_SIPMGROUPING_" + detectorName[i] + "_M=" + to_string(m)).c_str(), eHighN, eHighMin, eHighMax);
    //         }
    //     }
    // }

    // while (Reader->Next() && Reader->GetCurrentEntry() < 10000000)
    // {
    //     ProgressBar(Reader->GetCurrentEntry(), TotalEntries, start, Current, "Final Cleaning : ");
        
    //     // cout << (*signals) << endl;
    //     if ((*signals).Label == 103)
    //     {
    //         // cout << "NEW TRIGGER" << endl;
    //         vector<Signal> signals_SiPM; 
    //         signals_SiPM.push_back(*signals);   
    //         double Trig_time = (*signals).Time;
    //         int Trig_Entry = Reader->GetCurrentEntry(); 

    //         // going in the past
    //         while ((*signals).Time - Trig_time > past && Reader->GetCurrentEntry() > 0)
    //         {
    //             Reader->SetEntry(Reader->GetCurrentEntry() - 1);
    //             if (!IsDetectorBetaHigh((*signals).Label))
    //                 continue;
    //             // cout << "pst : " << setprecision(15) << (*signals).Time << " - " << Trig_time << " = " << (*signals).Time - Trig_time << endl;
    //             signals_SiPM.push_back(*signals);
    //         }

    //         Reader->SetEntry(Trig_Entry);

    //         // going in the future
    //         while ((*signals).Time - Trig_time < future && Reader->GetCurrentEntry() < TotalEntries)
    //         {
    //             Reader->Next();
    //             if (!IsDetectorBetaHigh((*signals).Label))
    //                 continue;
    //             // cout << setprecision(15) << (*signals).Time << " - " << Trig_time << " = " << (*signals).Time - Trig_time << endl;
    //             signals_SiPM.push_back(*signals);
    //         }

    //         // cout << "Future ok" << endl;
        
    //         int Multiplicity = 0;
    //         for (auto signal : signals_SiPM)
    //         {

    //             if(signal.Time - Trig_time > future || signal.Time - Trig_time < past)
    //                 continue;
    //             // cout << "        " << signal << endl;

    //             // cout << setprecision(15) << signal.Time << " - " << Trig_time << " = " << signal.Time - Trig_time << endl;

    //             if (signal.Time - Trig_time > start_gate && signal.Time - Trig_time < stop_gate)
    //             {
    //                 Multiplicity++;
    //             }
    //         }

    //         for (auto signal : signals_SiPM)
    //         {
    //             if(signal.Time - Trig_time > future || signal.Time - Trig_time < past)
    //                 continue;
    //             for (int m = 1; m <= Multiplicity; m++)
    //             {
    //                 H_Channel_SIPMGROUPING_M[signal.Label][m]->Fill(signal.Channel);
    //                 H_Time_SIPMGROUPING_M[signal.Label][m]->Fill(signal.Time - Trig_time);
    //                 // cout << "    m    " << m << endl;
    //             }
    //             H_Channel_SIPMGROUPING_Mequal[signal.Label][Multiplicity]->Fill(signal.Channel);
    //             H_Time_SIPMGROUPING_M_equal[signal.Label][Multiplicity]->Fill(signal.Time - Trig_time);
    //         }
    //     }
    // }

    // TCanvas *c_M[BETA_SIZE+1];
    // TLegend *leg_M[BETA_SIZE+1];
    // TCanvas *c_M_equal[BETA_SIZE+1];
    // TLegend *leg_M_equal[BETA_SIZE+1];
    // TCanvas *c_det[SIGNAL_MAX];
    // TLegend *leg_det[SIGNAL_MAX];
    // TCanvas *c_det_equal[SIGNAL_MAX];
    // TLegend *leg_det_equal[SIGNAL_MAX];
    // TCanvas *c_Time[SIGNAL_MAX];
    // TLegend *leg_Time[SIGNAL_MAX];

    // TDirectory *dir_det[SIGNAL_MAX];
    // TDirectory *dir_M[BETA_SIZE+1];


    // for (int det = 0; det < SIGNAL_MAX; det++)
    // {
    //     if (IsDetectorBetaHigh(det))
    //     {
    //         dir_det[det] = GROUPED_File->mkdir(detectorName[det].c_str());  
    //         c_det[det] = new TCanvas(("c_det_" + detectorName[det]).c_str(), ("c_det_" + detectorName[det]).c_str(), 800, 600);
    //         leg_det[det] = new TLegend(0.7, 0.7, 0.9, 0.9);
    //         c_det_equal[det] = new TCanvas(("c_det_equal_" + detectorName[det]).c_str(), ("c_det_equal_" + detectorName[det]).c_str(), 800, 600);
    //         leg_det_equal[det] = new TLegend(0.7, 0.7, 0.9, 0.9);
    //     }
    // }
    // for (int m = 1; m <= BETA_SIZE; m++)
    // {
    //     dir_M[m] = GROUPED_File->mkdir(("M" + to_string(m)).c_str());
    //     c_M[m] = new TCanvas(("c_M_" + to_string(m)).c_str(), ("c_M_" + to_string(m)).c_str(), 800, 600);
    //     leg_M[m] = new TLegend(0.7, 0.7, 0.9, 0.9);
    //     c_M_equal[m] = new TCanvas(("c_M_equal_" + to_string(m)).c_str(), ("c_M_equal_" + to_string(m)).c_str(), 800, 600);
    //     leg_M_equal[m] = new TLegend(0.7, 0.7, 0.9, 0.9);
    // }

    // for (int det = 0; det < SIGNAL_MAX; det++)
    // {
    //     if (IsDetectorBetaHigh(det))
    //     {
    //         for (int m = 1; m <= BETA_SIZE; m++)
    //         {
    //             H_Channel_SIPMGROUPING_M[det][m]->Write();
    //             H_Channel_SIPMGROUPING_Mequal[det][m]->Write();
    //             H_Time_SIPMGROUPING_M[det][m]->Write();
    //             H_Time_SIPMGROUPING_M_equal[det][m]->Write();

    //             c_det[det]->cd();
    //             H_Channel_SIPMGROUPING_M[det][m]->SetLineColor(m*10+1);
    //             leg_det[det]->AddEntry(H_Channel_SIPMGROUPING_M[det][m], ("M>=" + to_string(m)).c_str());
    //             H_Channel_SIPMGROUPING_M[det][m]->Draw("HIST SAME");
    //         }
    //         for (int m = BETA_SIZE; m >= 1; m--)
    //         {
    //             c_det_equal[det]->cd();
    //             H_Channel_SIPMGROUPING_Mequal[det][m]->SetLineColor(m*10+1);
    //             leg_det_equal[det]->AddEntry(H_Channel_SIPMGROUPING_Mequal[det][m], ("M=" + to_string(m)).c_str());
    //             H_Channel_SIPMGROUPING_Mequal[det][m]->Draw("HIST SAME");
    //         }
    //     }
    // }
                

    // for (int det = 0; det < SIGNAL_MAX; det++)
    // {
    //     if (IsDetectorBetaHigh(det))
    //     {
    //         dir_det[det]->cd(); 
    //         c_det[det]->cd();
    //         leg_det[det]->Draw("SAME");
    //         c_det[det]->SetTitle(detectorName[det].c_str());
    //         c_det[det]->Write();
    //         c_det_equal[det]->SetTitle(detectorName[det].c_str());
    //         c_det_equal[det]->cd();
    //         leg_det_equal[det]->Draw("SAME");
    //         c_det_equal[det]->Write();
    //     }
    // }

    // for (int det = 0; det < SIGNAL_MAX; det++)
    // {
    //     if (IsDetectorBetaHigh(det))
    //     {
    //         for (int m = 1; m <= BETA_SIZE; m++)
    //         {

    //             c_M[m]->cd();
    //             H_Channel_SIPMGROUPING_M[det][m]->SetLineColor(GetDetectorChannel(det) * 10+1);
    //             leg_M[m]->AddEntry(H_Channel_SIPMGROUPING_M[det][m], detectorName[det].c_str());
    //             H_Channel_SIPMGROUPING_M[det][m]->Draw("HIST SAME");
    //         }

    //         for (int m = BETA_SIZE; m >= 1; m--)
    //         {
    //             c_M_equal[m]->cd();
    //             H_Channel_SIPMGROUPING_Mequal[det][m]->SetLineColor(GetDetectorChannel(det) * 10+1);
    //             leg_M_equal[m]->AddEntry(H_Channel_SIPMGROUPING_Mequal[det][m], detectorName[det].c_str());
    //             H_Channel_SIPMGROUPING_Mequal[det][m]->Draw("HIST SAME");
    //         }
    //     }
    // }

    // for (int m = 1; m <= BETA_SIZE; m++)
    // {
    //     dir_M[m]->cd(); 
    //     c_M[m]->cd();
    //     c_M[m]->SetTitle(("M>=" + to_string(m)).c_str());
    //     leg_M[m]->Draw("SAME");
    //     c_M[m]->Write();
    //     c_M_equal[m]->cd();
    //     c_M_equal[m]->SetTitle(("M=" + to_string(m)).c_str());
    //     leg_M_equal[m]->Draw("SAME");
    //     c_M_equal[m]->Write();
    // }

    WriteHistograms_Cleaned();

    //////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Counting IAS losses /////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    WriteIASLosses();

    delete GROUPED_Tree;
    delete CUTTED_Tree;
    WriteTree_Grouped();

    GROUPED_File->Close();
    ROOT_File->Close();
    Success("Grouped File Created");
}
