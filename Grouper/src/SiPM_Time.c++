#include "SiPM_Time.hh"

int main()
{
    
    SIMULATED_File = new TFile((DIR_ROOT_DATA_SIMULATED + "../Time_test_analysed.root").c_str(), "READ");
    MERGED_File = new TFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
    H_Time_Sim = (TH1D*) SIMULATED_File->Get("e+/PlasticScintillator/PlasticScintillator_Time_e+");
    

    
    

    // TCanvas *c = (TCanvas*)MERGED_File->Get("C_Fake_9");
    // for (int i = 0; i < c->GetListOfPrimitives()->GetSize(); i++)
    // {   
    //     string s = c->GetListOfPrimitives()->At(i)->GetName();
    //     if (s.find("RearSiPM_Time_New_Nearest_M9") != string::npos)
    //     {
    //         H_Time_Exp = (TH1D*)c->GetListOfPrimitives()->At(i);
    //         break;
    //     }
    // }

    // cout << H_Time_Exp->GetEntries() << endl;
    
    ///// INPUT PARAMETERS /////
    double rate_DC = 1e6;
    double rate_Beta = 102219*5; //1372.43/2*530;
    double rate_Ar =  5.54393e+06*0.9;
    double total_rate = rate_DC + rate_Beta + rate_Ar;
    
    time_resolution = 2.5;

    double time_start = 400;
    double time_end = 230;
    double time_integration = 230;


    start_gate = -240;
    end_gate = 130;
    ////////////////////////////

    counter_real = 0;
    counter_fake = 0;

    fScintillator->SetNpx(10000);
    HScintillator = (TH1D*)fScintillator->GetHistogram();


    for (int event = 0; event < 6.0e6; event++)
    {
        //DETERMINE BETA TIME
        // uniform distribution shooting
        double t_beta = (double)rand() / RAND_MAX;
        t_beta = t_beta * (time_start+time_end) - time_start;

        //DETERMINE BETA Ar TIME
        // Gaussian distribution shootting
        double t_Ar = 0;//H_Time_Sim->GetRandom();

        //
        bool Ar = false;
        bool Beta = false;

        int i = 0;

        double RAr = (double)rand() / RAND_MAX * total_rate;
        if (RAr < rate_Ar)
        {
            Ar = true;
        }
        double RBeta = (double)rand() / RAND_MAX * total_rate;
        if (RBeta < rate_Beta)
        {
            Beta = true;
        }
        i++;

        if (Ar && Beta)
        {
            if (abs(t_beta - t_Ar) < time_integration)
            {
                if (t_beta < t_Ar)
                {
                    Fill(t_beta, Ar, Beta);
                }
                else
                {
                    Fill(t_Ar, Ar, Beta);
                }
            }
            else
            {
                if (abs(t_beta) < abs(t_Ar))
                {
                    Fill(t_beta, Ar, Beta);
                }
                else
                {
                    Fill(t_Ar, Ar, Beta);
                }
            }
        }
        else
        {
            if (Ar)
            {
                Fill(t_Ar, Ar, Beta);
            }
            else if (Beta)
            {
                Fill(t_beta, Ar, Beta);
            }
        }
    }

    ////////////

    H_SiPM_LeftLeft = (TH1D*)H_Time_MC->Clone("H_SiPM_LeftLeft_");
    for (int bin = 0; bin < H_SiPM_LeftLeft->GetNbinsX(); bin++)
    {
      if (H_SiPM_LeftLeft->GetBinCenter(bin) > -240 || H_SiPM_LeftLeft->GetBinCenter(bin) < -300)
      {
        H_SiPM_LeftLeft->SetBinContent(bin, 0);
      }
    }

    double integral_left_left = H_SiPM_LeftLeft->Integral();
    double mean_height_left_left = integral_left_left/30;

    H_SiPM_Center = (TH1D*)H_Time_MC->Clone("H_SiPM_Center_");
    for (int bin = 0; bin < H_SiPM_Center->GetNbinsX(); bin++)
    {
      if (H_SiPM_Center->GetBinCenter(bin) > 230 || H_SiPM_Center->GetBinCenter(bin) < -240)
      {
        H_SiPM_Center->SetBinContent(bin, 0);
      }
    }

    TH1D* H_Fortuitous = (TH1D*)H_SiPM_Center->Clone("_Fortuitous_");
    H_Fortuitous->Reset();
    for (int bin = 0; bin < H_SiPM_Center->GetNbinsX(); bin++)
    {
      if (H_SiPM_Center->GetBinCenter(bin) < -240 && H_SiPM_Center->GetBinCenter(bin) > 130)
      {
        H_Fortuitous->SetBinContent(bin, mean_height_left_left);
      }
      if (H_SiPM_Center->GetBinCenter(bin) > 130 || H_SiPM_Center->GetBinCenter(bin) < -240)
      {
        H_Fortuitous->SetBinContent(bin, H_Time_MC->GetBinContent(bin));
      }
    }
    cout << "ok" << endl;

    TH1D* H_Real = (TH1D*)H_SiPM_Center->Clone("_Real_");
    for (int bin = 0; bin < H_SiPM_Center->GetNbinsX(); bin++)
    {
      H_Real->SetBinContent(bin, H_SiPM_Center->GetBinContent(bin) - mean_height_left_left );
    }

    cout << "Real : " << counter_real << "    Fake : " << counter_fake << endl;
    cout << "Ratio : " << (double)counter_fake/counter_real*100 << " %" << endl;

    TFile *output = new TFile("Time_MC.root", "RECREATE");

    TCanvas *C_FakeCoincidence = new TCanvas("C_Fake_", "C_Fake_", 800, 400);
    // H_Time_MC->SetFillColor(kBlack);
    // H_Time_MC->SetFillStyle(3244);
    H_Time_MC->GetXaxis()->SetTitle("Time [ns]");
    H_Time_MC->GetYaxis()->SetTitle("Counts");
    H_Time_MC->GetXaxis()->CenterTitle();
    H_Time_MC->GetYaxis()->CenterTitle();
    H_Time_MC->Draw("HIST");
    // H_Fortuitous->SetFillColor(kRed);
    // H_Fortuitous->SetFillStyle(3244);
    H_Fortuitous->Draw("SAME");

    TLatex *text_low = new TLatex();
    text_low->SetNDC();
    text_low->SetTextSize(0.03);
    H_Real->GetXaxis()->SetRangeUser(-240, 130);
    double pourcentage = mean_height_left_left/2*370 * 100 / H_Real->Integral();
    text_low->DrawLatex(0.15, 0.8, ("Fake coincidences: " + to_string(pourcentage) + " %").c_str());
    text_low->Draw("SAME");
    C_FakeCoincidence->Write(); 

    // TCanvas *C_Windows = new TCanvas("C_Windows_", "C_Windows_", 800, 400);
    // H_SiPM_Center->Draw("HIST");
    // H_SiPM_Left->SetLineColor(kRed);
    // H_SiPM_Left->Draw("SAME");
    // H_SiPM_Right->SetLineColor(kRed);
    // H_SiPM_Right->Draw("SAME");
    // C_Windows->Write();

    H_Ar->Write();
    H_Fake->Write();

    output->Close();






    return 0;
}