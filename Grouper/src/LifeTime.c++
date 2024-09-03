#include "LifeTime.hh"
#include <fftw3.h>

int main()
{
    /// FOR 32Ar
    InitDetectors("Config_Files/sample.pid");
    
    GROUPED_File = new TFile((DIR_ROOT_DATA_GROUPED + "run_114_multifast_32Ar_grouped.root").c_str(), "READ");
    TFile *fTime = new TFile("Time.root", "RECREATE");
    InitPeakWindow();
    TTree *tree = (TTree *)GROUPED_File->Get("CLEANED_Tree");
    TTreeReader *Reader = new TTreeReader(tree);
    TTreeReaderArray<Signal> *Silicon = new TTreeReaderArray<Signal>(*Reader, "CLEANED_Tree_Silicon");

    H_Time_All = new TH1D("Time_All", "Time_All", 12000, 0, 12000);
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (!IsDetectorSiliStrip(i))
        {
            continue;
        }
        H_Exp[i] = new TH1D(Form("H_Exp_%d", i), Form("H_Exp_%d", i), 10000, 0, 100000);
        H_Time[i] = new TH1D(Form("H_Time_%d", i), Form("H_Time_%d", i), 10000, 0, 30e3);
        H_Time[i]->GetXaxis()->SetTitle("Time [s]");
        H_Time[i]->GetYaxis()->SetTitle("Counts");
        H_Time[i]->GetXaxis()->CenterTitle();
        H_Time[i]->GetYaxis()->CenterTitle();
    }   

    vector<double> vec_Time;
    double last_time = 0;
    double add_time = 0;
    while (Reader->Next())
    {
        ProgressBar(Reader->GetCurrentEntry(), Reader->GetEntries(true), 0, 0, "Tree : ");
        double Time = (*Silicon)[1].Time;
        double Channel = (*Silicon)[1].Channel;
        int Label = (*Silicon)[1].Label;

        if (Channel > peaks_window_F[Label].first && Channel < peaks_window_F[Label].second)
        {
            if (Time < last_time)
            {
                cout << "Time : " << Time << " Last Time : " << last_time << " Add Time : " << add_time << endl;
            }
            H_Exp[Label]->Fill(Channel);
            H_Time[Label]->Fill(Time*1e-9);
            H_Time_All->Fill(Time*1e-9);
            vec_Time.push_back(Time*1e-9);
            last_time = Time;
        }
    }

    H_Time_All->Write();
    

    // make a fourier transform of the time
    int N = 100000;     // Number of points
    double dt = 1e-9; // Time step

    double *in = (double *)fftw_malloc(sizeof(double) * vec_Time.size());
    fftw_complex *out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    // Fill the input array with your data
    for (int i = 0; i < vec_Time.size(); i++)
    {
        in[i] = vec_Time[i];
    }

    // Create a plan
    fftw_plan plan = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);

    // Execute the plan
    fftw_execute(plan);

    // Now, 'out' contains the Fourier Transform of 'in'

    // plot the result with root
    TCanvas *c = new TCanvas("c", "c", 800, 600);
    TH1D *h = new TH1D("h", "h", N, 0, 10);
    for (int i = 0; i < N; i++)
    {
        h->SetBinContent(i, out[i][1]);
    }
    h->Draw();
    c->Write();

    TCanvas *c1= new TCanvas("c1", "c1", 800, 600);
    TH1D *h1 = new TH1D("h1", "h1", N, 0, 10);
    for (int i = 0; i < N; i++)
    {
        h1->SetBinContent(i, out[i][0]);
    }
    h1->Draw();
    c1->Write();

    for (int i = 0; i < SIGNAL_MAX; i++)
    {
        if (!IsDetectorSiliStrip(i))
        {
            continue;
        }
        H_Exp[i]->Write();
        H_Time[i]->Write();
    }

    // Don't forget to free the memory and destroy the plan
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);



    fTime->Close();
}