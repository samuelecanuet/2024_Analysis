#include "SiPM_Time.hh"

int main()
{
    
    SIMULATED_File = new TFile((DIR_ROOT_DATA_SIMULATED + "../Time_test_analysed.root").c_str(), "READ");
    cout << "ok" << endl;
    MERGED_File = new TFile((DIR_ROOT_DATA_MERGED + "32Ar_merged.root").c_str(), "READ");
    cout << "ok" << endl;
    H_Time_Sim = (TH1D*) SIMULATED_File->Get("e+/PlasticScintillator/PlasticScintillator_Time_e+");
    cout << H_Time_Sim->GetEntries() << endl;

    
    

    TCanvas *c = (TCanvas*)MERGED_File->Get("C_FakeCoincidence_9");
    for (int i = 0; i < c->GetListOfPrimitives()->GetSize(); i++)
    {   
        string s = c->GetListOfPrimitives()->At(i)->GetName();
        if (s.find("RearSiPM_Time_New_Nearest_M9") != string::npos)
        {
            H_Time_Exp = (TH1D*)c->GetListOfPrimitives()->At(i);
            break;
        }
    }

    cout << H_Time_Exp->GetEntries() << endl;
    
    ///// INPUT PARAMETERS /////
    double rate_DC = 0;
    double rate_Beta = 363693.95; //1372.43/2*530;
    double rate_Ar =  5718660*0.94;
    double total_rate = rate_DC + rate_Beta + rate_Ar;
    
    double time_resolution = 4;

    double time_start = 300;
    double time_end = 230;
    double time_integration = 230;


    start_gate = -20;
    end_gate = 70;
    ////////////////////////////

    counter_real = 0;
    counter_fake = 0;


    for (int event = 0; event < 5.7e8*1373/26441; event++)
    {
        //DETERMINE BETA TIME
        // uniform distribution shooting
        double t_beta = (double)rand() / RAND_MAX;
        t_beta = t_beta * (time_start+time_end) - time_start;

        //DETERMINE BETA Ar TIME
        // Gaussian distribution shootting
        normal_distribution<double> distribution(0, time_resolution);
        double t_Ar = distribution(generator) + H_Time_Sim->GetRandom();

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

    cout << H_Time_MC->GetEntries() << endl;
    H_SiPM_LeftLeft = (TH1D*)H_Time_MC->Clone("_SiPM_LeftLeft_");
    for (int bin = 0; bin < H_SiPM_LeftLeft->GetNbinsX(); bin++)
    {
      if (H_SiPM_LeftLeft->GetBinCenter(bin) > -250 || H_SiPM_LeftLeft->GetBinCenter(bin) < -300)
      {
        H_SiPM_LeftLeft->SetBinContent(bin, 0);
      }
    }

    H_SiPM_Left = (TH1D*)H_Time_MC->Clone("_SiPM_Left_");
    for (int bin = 0; bin < H_SiPM_Left->GetNbinsX(); bin++)
    {
      if (H_SiPM_Left->GetBinCenter(bin) > -20 || H_SiPM_Left->GetBinCenter(bin) < -200)
      {
        H_SiPM_Left->SetBinContent(bin, 0);
      }
    }

    H_SiPM_Right = (TH1D*)H_Time_MC->Clone("_SiPM_Right_");
    for (int bin = 0; bin < H_SiPM_Right->GetNbinsX(); bin++)
    {
      if (H_SiPM_Right->GetBinCenter(bin) < 70 || H_SiPM_Right->GetBinCenter(bin) > 200)
      {
        H_SiPM_Right->SetBinContent(bin, 0);
      }
    }

    H_SiPM_Center = (TH1D*)H_Time_MC->Clone("_SiPM_Center_");
    for (int bin = 0; bin < H_SiPM_Center->GetNbinsX(); bin++)
    {
      if (H_SiPM_Center->GetBinCenter(bin) > 70 || H_SiPM_Center->GetBinCenter(bin) < -20)
      {
        H_SiPM_Center->SetBinContent(bin, 0);
      }
    }

    cout << "ok1" << endl;

    double integral_left_left = H_SiPM_LeftLeft->Integral();
    double mean_height_left_left = integral_left_left/25;

    double integral_left = H_SiPM_Left->Integral();
    double mean_height_left = integral_left/90;

    double integral_right = H_SiPM_Right->Integral();
    double mean_height_right = integral_right/65;

    double integral_center = H_SiPM_Center->Integral();

    cout << mean_height_left_left << "    " << mean_height_left << "    " << mean_height_right << "    " << integral_center << endl;

    H_S = (TH1D*)H_SiPM_Center->Clone("_S_");
    for (int bin = 0; bin < H_S->GetNbinsX(); bin++)
    {
      H_S->SetBinContent(bin, H_SiPM_Center->GetBinContent(bin) - mean_height_left_left);
    }

    cout << "ok2" << endl;


    double guess_N = H_S->Integral();

    cout << "Guess N = " << guess_N << endl;

    Minimizer *minimizer = Factory::CreateMinimizer("Minuit2", "Migrad");
    ROOT::Math::Functor functor(&FunctionToMinimize, 4);
    minimizer->SetFunction(functor);
    minimizer->SetLimitedVariable(0, "Normalization", guess_N, 1, guess_N*0.9, guess_N);
    minimizer->SetFixedVariable(1, "Multiplicity", -1);
    minimizer->SetFixedVariable(2, "mean_height_left_left", mean_height_left_left);
    minimizer->SetFixedVariable(3, "mean_height_left", mean_height_left);
    // minimizer->SetPrecision(0.001);
    // minimizer->SetTolerance(0.001);
    minimizer->SetMaxFunctionCalls(1000000);
    minimizer->SetMaxIterations(1000000);

    minimizer->Minimize();
    const double *par = minimizer->X();

    cout << "ok3" << endl;

    G_N->Fit("pol1", "E");

    double a = G_N->GetFunction("pol1")->GetParameter(1);
    double b = G_N->GetFunction("pol1")->GetParameter(0);
    double err_a = G_N->GetFunction("pol1")->GetParError(1);
    double err_b = G_N->GetFunction("pol1")->GetParError(0);

    guess_N = -b/a;

    // compute error on guess_N
    double error_guess_N = sqrt(pow(err_b/a, 2) + pow(b*err_a/pow(a, 2), 2));

    cout << "Guess N = " << guess_N << " +/- " << error_guess_N << endl;

    guess_N = par[0];
    // compute alpha
    double alpha = (mean_height_left-mean_height_left_left);
    // compute f(t)
    TH1D* H_Fortuitous = (TH1D*)H_SiPM_Center->Clone("_Fortuitous_");
    H_Fortuitous->Reset();
    TH1D* H_fx = (TH1D*)H_SiPM_Center->Clone("_fx_");
    TH1D* H_Sx = (TH1D*)H_Time_MC->Clone("_Sx_");
    H_fx->Reset();
    H_Sx->Reset();
    cout << "alpha = " << alpha << endl;
    double Fx = 0;
    for (int bin = 0; bin < H_SiPM_Center->GetNbinsX(); bin++)
    {
        if (H_S->GetBinCenter(bin) < -20 || H_S->GetBinCenter(bin) > 70)
        {
            H_Fortuitous->SetBinContent(bin, mean_height_left_left);
        }
        else
        {
            double fx = (H_S->GetBinContent(bin) - alpha * (1 - Fx / guess_N)) / (1 - alpha / guess_N);
            H_Fortuitous->SetBinContent(bin, alpha * (1 - Fx / guess_N) + mean_height_left_left);
            H_fx->SetBinContent(bin, fx);

            Fx += fx;
        }
    }

    H_Sx->Add(H_fx, 1);
    H_Sx->Add(H_Fortuitous, 1);
    /////////////

    cout << "Real : " << counter_real << "    Fake : " << counter_fake << endl;
    cout << "Ratio : " << (double)counter_fake/counter_real*100 << " %" << endl;

    TFile *output = new TFile("Time_MC.root", "RECREATE");

    TCanvas *C_Fake = new TCanvas("C_FakeMC", "C_Fake", 800, 400);
    H_Time_MC->Draw("HIST");
    H_Fake->SetLineColor(kRed);
    H_Fake->Draw("HIST SAME");
    C_Fake->Write();




    int mul = 0 ;
    TCanvas *C_Fortuitous = new TCanvas(("C_FakeCoincidence_" + to_string(mul)).c_str(), ("C_FakeCoincidence_" + to_string(mul)).c_str(), 800, 400);
    H_Time_MC->Draw("HIST");
    H_Fortuitous->SetLineColor(kRed);
    H_Fortuitous->Draw("SAME");
    H_Sx->SetLineColor(kGreen);
    H_Sx->Draw("SAME");
    C_Fortuitous->Write();

    TCanvas *C_FakeCoincidence = new TCanvas(("C_Fake_" + to_string(mul)).c_str(), ("C_Fake_" + to_string(mul)).c_str(), 800, 400);
    H_Time_MC->Draw("HIST");
    H_SiPM_Center->SetFillColor(kBlack);
    H_SiPM_Center->SetFillStyle(3244);
    H_SiPM_Center->Draw("SAME");
    H_Fortuitous->GetXaxis()->SetRangeUser(-20, 70);
    H_Fortuitous->SetFillStyle(3244);
    H_Fortuitous->SetFillColor(kRed);
    H_Fortuitous->Draw("SAME");

    TLatex *text_low = new TLatex();
    text_low->SetNDC();
    text_low->SetTextSize(0.03);
    H_fx->GetXaxis()->SetRangeUser(-20, 70);
    H_Fortuitous->GetXaxis()->SetRangeUser(-20, 70);
    double pourcentage = H_Fortuitous->Integral() * 100 / (H_fx->Integral());
    text_low->DrawLatex(0.15, 0.8, ("Fake coincidences: " + to_string(pourcentage) + " %").c_str());
    text_low->Draw("SAME");
    C_FakeCoincidence->Write(); 

    TCanvas *C_Windows = new TCanvas(("C_Windows_" + to_string(mul)).c_str(), ("C_Windows_" + to_string(mul)).c_str(), 800, 400);
    H_SiPM_Center->Draw("HIST");
    H_SiPM_Left->SetLineColor(kRed);
    H_SiPM_Left->Draw("SAME");
    H_SiPM_Right->SetLineColor(kRed);
    H_SiPM_Right->Draw("SAME");
    C_Windows->Write();

    TCanvas *C_N = new TCanvas(("C_N_" + to_string(mul)).c_str(), ("C_N_" + to_string(mul)).c_str(), 800, 400);
    G_N->GetXaxis()->SetTitle("N");
    G_N->GetYaxis()->SetTitle("N - #\int f(t) dt");
    G_N->SetMarkerStyle(20);
    G_N->SetMarkerSize(1);
    G_N->Draw("AP");
    C_N->Write();


    output->Close();






    return 0;
}