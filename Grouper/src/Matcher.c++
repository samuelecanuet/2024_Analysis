#include "Matcher.hh"


int main(int argc, char *argv[])
{
    string Run_string;

    if (argc == 1)
    {
        Error("No Run Number Given");
    }
    else
    {
        Run = atoi(argv[1]);
        if (Run < 100)
            Run_string = "0" + to_string(Run);
        else
            Run_string = to_string(Run);
    }

    Info("Current Run : " + Run_string);
    InitDetectors("Config_Files/sample.pid");
    ///////////////////////////////////  INPUT ///////////////////////////////////
    GROUPED_filename = SearchFiles(DIR_ROOT_DATA_GROUPED, Run_string);
    GROUPED_basefilename = GROUPED_filename.substr(0, GROUPED_filename.find(".root"));
    GROUPED_File = new TFile((DIR_ROOT_DATA_GROUPED+GROUPED_filename).c_str(), "READ");
    
    REFERENCE_filename = SearchFiles(DIR_ROOT_DATA_GROUPED, "0" + to_string(REFERENCE_Run));
    REFERENCE_File = new TFile((DIR_ROOT_DATA_GROUPED+REFERENCE_filename).c_str(), "READ");

    ///////////////////////////////////  OUTPUT ///////////////////////////////////
    MATCHED_File = new TFile((DIR_ROOT_DATA_MATCHED+GROUPED_basefilename+"_matched.root").c_str(), "RECREATE");
    MATCHED_File->cd();
    WriteTime(GROUPED_File, MATCHED_File);

    ///////////////////////////////////////////////////////////////////////////////
    TTree *Tree = (TTree *)GROUPED_File->Get("Fits/CLEANED_Tree");
    Reader = new TTreeReader(Tree);
    Silicon = new TTreeReaderArray<Signal>(*Reader, "CLEANED_Tree_Silicon");

    clock_t start = clock(), Current;
    int TotalEntries = Reader->GetEntries();

    ///////////////////////////////////  REFERENCE  //////////////////////////////
    for (int i = 0 ; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
            H_Run_Ref[i] = (TH1D*)REFERENCE_File->Get(("Strip/Strip_Channel/" + detectorName[i] + "/Channel_Cleaned_" + detectorName[i]).c_str());
        
    }
    ///////////////////////////////////////////////////////////////////////////////

    //Create exp tree for each detector
    for (int i = 0 ; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            GROUPED_Tree_Detectors[i] = new TTree(("GROUPED_Tree_" + detectorName[i]).c_str(), ("GROUPED_Tree_" + detectorName[i]).c_str());
            GROUPED_Tree_Detectors[i]->Branch("Channel", &Channel);
        }
    }

    //Filling exp tree for each detector
    while (Reader->Next())
    {
        ProgressBar(Reader->GetCurrentEntry(), TotalEntries, start, Current, "Reading Tree : ");
        Channel = (*Silicon)[1].Channel;
        GROUPED_Tree_Detectors[(*Silicon)[1].Label]->Fill();
    }

    //Init histograms
    InitHistograms();

    for (int i = 0 ; i < SIGNAL_MAX; i++)
    {
        if (IsDetectorSiliStrip(i))
        {
            Info(detectorName[i]);
            Reader = new TTreeReader(GROUPED_Tree_Detectors[i]);
            CHI2Minimization(i);
        }
    }
    
    MATCHED_File->Close();
    GROUPED_File->Close();

    return 0;
}
