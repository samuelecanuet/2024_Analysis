#include "SIPM_Calibration.hh"


int main()
{
    // INITIALIZATION
    InitDetectors("Config_Files/sample.pid");


    ///////////////////////////////////  FILES //////////////////////////////////
    GROUPED_File["32Ar"] = new TFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
    GROUPED_File["207Bi"] = new TFile((DIR_ROOT_DATA + "run_137_multifast_207Bi.root").c_str(), "READ");
    GROUPED_File["90Sr"] = new TFile((DIR_ROOT_DATA + "run_133_data_90Sr.root").c_str(), "READ");
    SIMULATED_File["32Ar"] = new TFile((DIR_ROOT_DATA_SIMULATED + "/30-09/32Ar_CS0_CSP0_CV1_CVP1_analysed.root").c_str(), "READ");
    SIMULATED_File["207Bi"] = new TFile((DIR_ROOT_DATA_SIMULATED + "../207Bi_a1_b0_analysed.root").c_str(), "READ");
    SIMULATED_File["90Sr"] = new TFile((DIR_ROOT_DATA_SIMULATED + "../90Sr_a1_b0_analysed.root").c_str(), "READ");

    ///////////////////////////////////  OUTPUT ///////////////////////////////////
    CALIBRATED_File = new TFile((DIR_ROOT_DATA_CALIBRATED + "SiPM_Calibrated.root").c_str(), "RECREATE");
    CALIBRATED_File->cd();

    NUCLEUS = "32Ar";
    f_tree = new TFile((DIR_ROOT_DATA_CALIBRATED + "Matching_SiPM_trees.root").c_str(), "READ");

    InitWindows();
    InitSiliconCalibration();
    InitHistograms();

    NUCLEUS = "32Ar";
    // Reader for grouped
    Tree = (TTree *)GROUPED_File[NUCLEUS]->Get("MERGED_Tree");
    Reader = new TTreeReader(Tree);
    Silicon = new TTreeReaderArray<Signal>(*Reader, "MERGED_Tree_Silicon");
    SiPM_High = new TTreeReaderArray<Signal>(*Reader, "MERGED_Tree_SiPMHigh");
    SiPM_Low = new TTreeReaderArray<Signal>(*Reader, "MERGED_Tree_SiPMLow");
    // Reader for simulated
    for (string NUCLEUS : Nucleis)
    {
        if (NUCLEUS == "32Ar")
            continue;
        Tree_SIMULATED = (TTree *)SIMULATED_File[NUCLEUS]->Get("PlasticIAS");
        Reader_SIMULATED = new TTreeReader(Tree_SIMULATED);
        Silicon_code = new TTreeReaderValue<int>(*Reader_SIMULATED, "Code");
        Silicon_energy = new TTreeReaderValue<double>(*Reader_SIMULATED, "Energy");
        SiPM_energy = new TTreeReaderValue<double>(*Reader_SIMULATED, "SiPM");

        clock_t start = clock(), Current;
        int Entries = Reader_SIMULATED->GetEntries();
        while (Reader_SIMULATED->Next())
        {
            ProgressBar(Reader_SIMULATED->GetCurrentEntry(), Entries, start, Current, "Reading Tree");
            int sili_code = **Silicon_code;
            double sili_e = **Silicon_energy;
            double SiPM_e = **SiPM_energy;

            if (NUCLEUS == "207Bi" || NUCLEUS == "90Sr")
            {
                H_Sim[NUCLEUS][0]->Fill(SiPM_e);
            }
            else
            {
                for (int peak_number = 0; peak_number < 50; peak_number++)
                {
                    if (WindowsMap[NUCLEUS][peak_number][sili_code].first == -1 || !WindowsMap[NUCLEUS][peak_number][sili_code].first)
                        continue;

                    if (sili_e > WindowsMap[NUCLEUS][peak_number][sili_code].first & sili_e < WindowsMap[NUCLEUS][peak_number][sili_code].second)
                    {
                        H_Sim[NUCLEUS][peak_number]->Fill(SiPM_e);
                    }
                }
            }
        }
    }

    NUCLEUS = "32Ar";
    ReadTree();    


    for (int det  = 1; det <= 9; det++)
    {
        current_detector = det;
        CalibrationSiPM();
    }
   
    
    WriteHistograms();
    CALIBRATED_File->Close();
    SIMULATED_File["32Ar"]->Close();
    GROUPED_File["32Ar"]->Close();
    return 0;
}